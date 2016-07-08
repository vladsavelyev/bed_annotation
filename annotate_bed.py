#!/usr/bin/env python
from optparse import OptionParser, SUPPRESS_HELP
from collections import defaultdict
from os.path import isfile
from pybedtools import BedTool

import GeneAnnotation
from Utils import reference_data
from Utils.logger import warn, debug
from Utils.utils import OrderedDefaultDict
from Utils.bed_utils import verify_bed, bedtools_version, SortableByChrom, cut
from Utils.file_utils import verify_file, file_transaction, open_gzipsafe, which, intermediate_fname
from Utils.logger import critical, info


def main():
    options = [
        (['-o', '--output-file'], dict(
            dest='output_file',
            metavar='FILE',
            help='Output file',
        )),
        (['--reference', '--features'], dict(
            dest='features',
        )),
        (['--reuse'], dict(
            dest='reuse_intermediate',
            help='reuse intermediate non-empty files in the work dir from previous run',
            action='store_true',
        )),
        (['--work-dir'], dict(dest='work_dir', metavar='DIR', help=SUPPRESS_HELP)),
        (['--log-dir'], dict(dest='log_dir', metavar='DIR', help=SUPPRESS_HELP)),
        (['-g', '--genome'], dict(dest='genome')),
        (['--canonical'], dict(dest='only_canonical')),
    ]

    parser = OptionParser(description='Annotating BED file based on reference features annotations.')
    for args, kwargs in options:
        parser.add_option(*args, **kwargs)
    opts, args = parser.parse_args()

    if len(args) < 1:
        parser.exit('Usage: ' + __file__ + ' Input_BED_file -g hg19 -o Annotated_BED_file [--canonical]')
    input_bed_fpath = verify_bed(args[0], is_critical=True, description='Input BED file for ' + __file__)

    features_fpath = None
    if opts.features:
        features_fpath = verify_bed(opts.features, is_critical=True)
    elif opts.genome:
        if opts.only_canonical:
            features_fpath = verify_file(GeneAnnotation.get_all_features_canonical(opts.genome))
        else:
            features_fpath = verify_file(GeneAnnotation.get_all_features(opts.genome))
        if not features_fpath:
            critical('Genome ' + opts.genome + ' is not supported. Supported: ' + ', '.join(GeneAnnotation.SUPPORTED_GENOMES))
    else:
        critical('Error: neither --features nor --genome is specified. Features are required for annotation.')

    output_fpath = annotate(input_bed_fpath, features_fpath, opts.output_file, genome=opts.genome)

    info('Done, saved to ' + output_fpath)


def bed_chrom_order(bed_fpath):
    chroms = []
    chroms_set = set()
    with open(bed_fpath) as f:
        for l in f:
            l = l.strip()
            if l and not l.startswith('#'):
                chrom = l.split('\t')[0]
                if chrom not in chroms_set:
                    chroms_set.add(chrom)
                    chroms.append(chrom)
    chr_order = {c: i for i, c in enumerate(chroms)}
    return chr_order


def annotate(input_bed_fpath, features_bed_fpath, output_fpath, reuse=False, genome=None):
    if reuse and isfile(output_fpath) and verify_file(output_fpath):
        debug(output_fpath + ' exists, reusing.')
        return output_fpath

    if genome:
        fai_fpath = reference_data.get_fai(genome)
        chr_order = reference_data.get_chrom_order(genome)
    else:
        fai_fpath = None
        chr_order = bed_chrom_order(input_bed_fpath)

    bedtools = which('bedtools')
    bedtools_v = bedtools_version(bedtools)
    bed = BedTool(input_bed_fpath).cut([0, 1, 2])
    info()

    annotated = []
    off_targets = []
    if bedtools_v > 25 or not features_bed_fpath.endswith('.gz'):
        ref_bed = BedTool(features_bed_fpath)
    else:
        ref_bed = BedTool(open_gzipsafe(features_bed_fpath)).saveas()

    for feature in ['CDS', 'Exon', 'Transcript', 'Gene']:
        if bed is not None:
            info('Extracting RefSeq features')
            _ref_bed = ref_bed.filter(lambda x: x[6] == feature)

            info('Annotating based on ' + feature + '...')
            new_annotated, off_targets = _annotate(bed, _ref_bed, chr_order, fai_fpath)
            if not annotated:
                annotated = new_annotated
                for a in annotated:
                    a.feature = feature
            else:
                annotated.extend(new_annotated)

            if off_targets:
                bed = BedTool([(r.chrom, r.start, r.end) for r in off_targets])

                # off_target_fpath = _save_regions(off_targets, join(work_dirpath, 'off_target_1.bed'))
                # log('Saved off target1 to ' + str(off_target_fpath))
                info()
            else:
                break

    annotated.extend(off_targets)

    info('Saving annotated regions...')
    with file_transaction(None, output_fpath) as tx:
        with open(tx, 'w') as out:
            for region in sorted(annotated, key=lambda r: r.get_key()):
                out.write(str(region))
    return output_fpath


class Region(SortableByChrom):
    def __init__(self, chrom, start, end, ref_chrom_order, gene_symbol=None, exon=None, strand=None, feature=None, biotype=None):
        SortableByChrom.__init__(self, chrom, ref_chrom_order)
        self.chrom = chrom
        self.start = start
        self.end = end
        self.symbol = gene_symbol
        self.exon = exon
        self.strand = strand
        self.feature = feature
        self.biotype = biotype
        self.total_merged = 0

    def __str__(self):
        fs = [self.chrom,
              '{}'.format(self.start),
              '{}'.format(self.end),
              self.symbol or '.',
              self.exon or '.',
              self.strand or '.',
              self.biotype or '.']
        return '\t'.join(fs) + '\n'

    def get_key(self):
        return SortableByChrom.get_key(self), self.start, self.end, self.symbol


def _merge_fields(consensus_field, other_field):
    if not consensus_field:
        consensus_field = other_field
    else:
        consensus_field = ','.join(set(consensus_field.split(',')) | set(other_field.split(',')))
    return consensus_field


def _resolve_ambiguities(annotated_by_loc_by_gene, chrom_order):
    annotated = []
    for (chrom, start, end), overlaps_by_gene in annotated_by_loc_by_gene.iteritems():
        for g_name, overlaps in overlaps_by_gene.iteritems():
            consensus = Region(chrom, start, end, ref_chrom_order=chrom_order.get(chrom), gene_symbol=g_name, exon='', strand='', feature='', biotype='')
            for r, overlap_size in overlaps:
                if consensus.strand:
                    # RefSeq has exons from different strands with the same gene name (e.g. CTAGE4 for hg19),
                    # Such pair of exons may overlap with a single region, so taking strand from the first one
                    if consensus.strand != r.strand:
                        warn('Warning: different strands between consensus and next region (gene: ' + g_name + ')')
                    # assert consensus.strand == r.strand, 'Consensus strand is ' + \
                    #     consensus.strand + ', region strand is ' + r.strand
                else:
                    consensus.strand = r.strand
                consensus.exon = _merge_fields(consensus.exon, r.exon)
                consensus.feature = _merge_fields(consensus.feature, r.feature)
                consensus.biotype = _merge_fields(consensus.biotype, r.biotype)
                consensus.total_merged += 1

            annotated.append(consensus)

    return annotated


def _annotate(bed, ref_bed, chr_order, fai_fpath=None):
    # if genome:
        # genome_fpath = cut(fai_fpath, 2, output_fpath=intermediate_fname(work_dir, fai_fpath, 'cut2'))
        # intersection = bed.intersect(ref_bed, sorted=True, wao=True, g='<(cut -f1,2 ' + fai_fpath + ')')
        # intersection = bed.intersect(ref_bed, sorted=True, wao=True, genome=genome.split('-')[0])
    # else:

    intersection = bed.intersect(ref_bed, wao=True)
    # intersection = bed.intersect(ref_bed, wao=True, sorted=True, g=fai_fpath)

    total_annotated = 0
    total_uniq_annotated = 0

    met = set()

    annotated_by_loc_by_gene = OrderedDefaultDict(lambda: defaultdict(list))
    off_targets = list()

    for fs in intersection:
        fs = [f for f in fs]
        if len(fs) < 12:
            critical('Cannot parse the reference BED file - unexpected number of lines '
                     '(' + str(len(fs)) + ') in ' + str(fs))

        a_chr, a_start, a_end, e_chr, e_start, e_end, e_gene, e_exon, e_strand, \
            e_feature, e_biotype, e_transcript = fs[:12]
        overlap_size = int(fs[-1])

        assert e_chr == '.' or a_chr == e_chr, str((a_chr + ', ' + e_chr))

        if e_chr == '.':
            off_targets.append(Region(a_chr, int(a_start), int(a_end), ref_chrom_order=chr_order.get(a_chr)))
        else:
            total_annotated += 1
            if (a_chr, a_start, a_end) not in met:
                total_uniq_annotated += 1

            annotated_by_loc_by_gene[(a_chr, int(a_start), int(a_end))][e_gene].append((
                Region(chrom=e_chr, start=int(e_start), end=int(e_end), ref_chrom_order=chr_order.get(a_chr),
                       gene_symbol=e_gene, exon=e_exon, strand=e_strand, feature=e_feature, biotype=e_biotype),
                       overlap_size))

        met.add((a_chr, a_start, a_end))

    info('  Total annotated regions: ' + str(total_annotated))
    info('  Total uniq annotated regions: ' + str(total_uniq_annotated))
    info('  Total off target regions: ' + str(len(off_targets)))
    info('Resolving ambiguities...')
    annotated = _resolve_ambiguities(annotated_by_loc_by_gene, chr_order)

    return annotated, off_targets


def _save_regions(regions, fpath):
    with open(fpath, 'w') as off_target_f:
        for r in regions:
            off_target_f.write(str(r))

    return fpath


def _split_reference_by_priority(cnf, features_bed_fpath):
    features = ['CDS', 'Exon', 'Transcript', 'Gene']
    info('Splitting the reference file into ' + ', '.join(features))
    features_and_beds = []
    for f in features:
        features_and_beds.append((f, BedTool(features_bed_fpath).filter(lambda x: x[6] == f)))
    return features_and_beds


if __name__ == '__main__':
    main()
