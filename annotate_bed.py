#!/usr/bin/env python
from itertools import izip, count
from optparse import OptionParser, SUPPRESS_HELP
from collections import defaultdict
from os.path import isfile
from pybedtools import BedTool

import GeneAnnotation as ga
from Utils import reference_data
from Utils.logger import warn, debug
from Utils.utils import OrderedDefaultDict
from Utils.bed_utils import verify_bed, bedtools_version, SortableByChrom, cut, count_bed_cols
from Utils.file_utils import verify_file, file_transaction, open_gzipsafe, which, intermediate_fname
from Utils.logger import critical, info
from Utils import logger


def main():
    options = [
        (['-o', '--output-file'], dict(
            dest='output_file',
            metavar='FILE',
            help='Output file',
        )),
        (['--output-features'], dict(
            dest='output_features',
            action='store_true',
            default=False,
            help='Also output featues that used to annotate',
        )),
        (['--reference', '--features'], dict(
            dest='features',
        )),
        (['--reuse'], dict(
            dest='reuse_intermediate',
            action='store_true',
            help='reuse intermediate non-empty files in the work dir from previous run',
        )),
        (['-g', '--genome'], dict(
            dest='genome',
            help='Genome build. Accepted values: ' + ', '.join(ga.SUPPORTED_GENOMES),
        )),
        (['--canonical'], dict(
            dest='only_canonical',
            action='store_true',
            default=False,
            help='Use only features from canonical transcripts to annotate',
        )),
        (['--short'], dict(
            dest='short',
<<<<<<< HEAD
            action='store_true',
            default=True,
            help='Add only "Gene" column',
        )),
        (['--extended'], dict(
            dest='short',
            action='store_false',
            default=True,
            help='Add additional columns: transcript, GC, overlap size...',
        )),
        (['--high-confidence'], dict(
            dest='high_confidence',
            action='store_true',
            default=False,
            help='Annotate with only high confidence regions (TSL is 1 or NA, HUGO gene annotated, total overlap size > 50%)',
        )),
        (['--seq2c'], dict(
            dest='seq2c',
            action='store_true',
            default=False,
            help='Equals to 3 flags set: --canonical, --short, --high-confidence',
        )),
=======
            action='store_true',
            default=True,
            help='Add only "Gene" column',
        )),
        (['--extended'], dict(
            dest='short',
            action='store_false',
            default=True,
            help='Add additional columns: transcript, GC, overlap size...',
        )),
        (['--high-confidence'], dict(
            dest='high_confidence',
            action='store_true',
            default=False,
            help='Annotate with only high confidence regions (TSL is 1 or NA, HUGO gene annotated, total overlap size > 50%)',
        )),
        (['--seq2c'], dict(
            dest='seq2c',
            action='store_true',
            default=False,
            help='Equals to 3 flags set: --canonical, --short, --high-confidence',
        )),
>>>>>>> 83af7a28d38c14e9c1287475e1cbdc985eb836ce
        (['--debug'], dict(
            dest='debug',
            action='store_true',
            default=False,
            help=SUPPRESS_HELP,
         )),
        (['--not-collapse-exons'], dict(
            dest='collapse_exons',
            action='store_false',
            default=True,
            help=SUPPRESS_HELP,
         )),
        (['--work-dir'], dict(dest='work_dir', metavar='DIR', help=SUPPRESS_HELP)),
        (['--log-dir'], dict(dest='log_dir', metavar='DIR', help=SUPPRESS_HELP)),
    ]

    # TODO:
    # merge BED if not

    parser = OptionParser(description='Annotating BED file based on reference features annotations.')
    for args, kwargs in options:
        parser.add_option(*args, **kwargs)
    opts, args = parser.parse_args()
    logger.is_debug = opts.debug

    if opts.short and opts.output_features:
        critical('--short and --output-features can\'t be set both')

    if opts.seq2c:
        opts.short = opts.high_confidence = opts.only_canonical = True
        opts.output_features = False

    if len(args) < 1:
        parser.exit('Usage: ' + __file__ + ' Input_BED_file -g hg19 -o Annotated_BED_file [--canonical]')
    input_bed_fpath = verify_bed(args[0], is_critical=True, description='Input BED file for ' + __file__)

    features_bed = None
    if opts.features:
        debug('Getting features from ' + opts.features)
        features_bed = BedTool(verify_bed(opts.features, is_critical=True))
    elif opts.genome:
        debug('Getting features from storage')
        features_bed = ga.get_all_features(opts.genome)
        if not features_bed:
            critical('Genome ' + opts.genome + ' is not supported. Supported: ' + ', '.join(ga.SUPPORTED_GENOMES))
    else:
        critical('Error: neither --features nor --genome is specified. Features are required for annotation.')

    output_fpath = annotate(
        input_bed_fpath, features_bed, opts.output_file, genome=opts.genome,
        only_canonical=opts.only_canonical, short=opts.short, high_confidence=opts.high_confidence,
        collapse_exons=opts.collapse_exons, output_features=opts.output_features)

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


<<<<<<< HEAD
def get_sort_key(chr_order):
    return lambda fs: (chr_order[fs[ga.BedCols.CHROM]],
            int(fs[ga.BedCols.START]),
            int(fs[ga.BedCols.END]),
            0 if fs[ga.BedCols.FEATURE] == 'transcript' else 1,
            fs[ga.BedCols.ENSEMBL_ID],
            fs[ga.BedCols.GENE])


=======
>>>>>>> 83af7a28d38c14e9c1287475e1cbdc985eb836ce
def annotate(input_bed_fpath, features_bed, output_fpath, reuse=False, genome=None,
             only_canonical=False, short=False, high_confidence=False,
             collapse_exons=True, output_features=False):

    if reuse and isfile(output_fpath) and verify_file(output_fpath):
        debug(output_fpath + ' exists, reusing.')
        return output_fpath

    if genome:
        fai_fpath = reference_data.get_fai(genome)
        chr_order = reference_data.get_chrom_order(genome)
    else:
        fai_fpath = None
        chr_order = bed_chrom_order(input_bed_fpath)

    bed = BedTool(input_bed_fpath).cut([0, 1, 2])
    info()

<<<<<<< HEAD
    info('Extracting features')
    _ref_bed = features_bed.filter(lambda x:
        x[ga.BedCols.FEATURE] in ['exon', 'CDS', 'stop_codon', 'transcript']
        and (not high_confidence or x[ga.BedCols.TSL] in ['1', 'NA'] and x[ga.BedCols.HUGO]))

    info('Annotating...')
    annotated, off_targets = _annotate(bed, _ref_bed, chr_order, fai_fpath,
        high_confidence, short, collapse_exons, output_features)
    annotated.extend(off_targets)

    # annotated.sort(key=sort_key)
=======
    annotated = []
    off_targets = []

    filters = [
        lambda x: x[ga.BedCols.FEATURE] in ['exon', 'CDS', 'stop_codon', 'transcript']
                  and (not high_confidence or x[ga.BedCols.TSL] in ['1', 'NA']
                  and x[ga.BedCols.HUGO]),
    ]

    for filter_query in filters:
        if bed is not None:
            info('Extracting features')
            _ref_bed = features_bed.filter(filter_query)

            info('Annotating...')
            new_annotated, off_targets = _annotate(bed, _ref_bed, chr_order, fai_fpath,
                   high_confidence, short, collapse_exons)
            if not annotated:
                annotated = new_annotated
            else:
                annotated.extend(new_annotated)

            if off_targets:
                info()
            else:
                break
>>>>>>> 83af7a28d38c14e9c1287475e1cbdc985eb836ce

    header = [ga.BedCols.names[i] for i in ga.BedCols.cols]

    def sort_key(fs):
        return (chr_order[fs[ga.BedCols.CHROM]],
                int(fs[ga.BedCols.START]),
                int(fs[ga.BedCols.END]),
                0 if fs[ga.BedCols.FEATURE] == 'transcript' else 1,
                fs[ga.BedCols.ENSEMBL_ID],
                fs[ga.BedCols.GENE],
                )
    annotated.sort(key=sort_key)

    header = [ga.BedCols.names[i] for i in ga.BedCols.cols]

    info('Saving annotated regions...')
    with file_transaction(None, output_fpath) as tx:
        with open(tx, 'w') as out:
            out.write('\t'.join(header) + '\n')
<<<<<<< HEAD
            for fields in annotated:
=======
            for fields in sorted(annotated, key=sort_key):
>>>>>>> 83af7a28d38c14e9c1287475e1cbdc985eb836ce
                out.write('\t'.join(map(_format_field, fields)) + '\n')
    return output_fpath


class Region(SortableByChrom):
    def __init__(self, chrom, start, end, ref_chrom_order, gene_symbol=None, exon=None,
                 strand=None, other_fields=None):
        SortableByChrom.__init__(self, chrom, ref_chrom_order)
        self.chrom = chrom
        self.start = start
        self.end = end
        self.symbol = gene_symbol
        self.exon = exon
        self.strand = strand
        self.other_fields = other_fields or []
        self.total_merged = 0

    def __str__(self):
        fs = [
            self.chrom,
            self.start,
            self.end,
            self.symbol,
            self.exon,
            self.strand
        ] + self.other_fields
        fs = map(_format_field, fs)
        return '\t'.join(fs) + '\n'

    def get_key(self):
        return SortableByChrom.get_key(self), self.start, self.end, self.symbol

    def header(self):
        pass


def _format_field(value):
    if isinstance(value, list) or isinstance(value, set):
        return ', '.join(_format_field(v) for v in value)
    elif isinstance(value, float):
        return '{:.1f}'.format(value)
    else:
        return str(value or '.')
<<<<<<< HEAD
=======

>>>>>>> 83af7a28d38c14e9c1287475e1cbdc985eb836ce

def _resolve_ambiguities(annotated_by_loc_by_gene, chrom_order, collapse_exons=True):
    def _collapse(regions):
        cons = []
        for i in ga.BedCols.cols:
            rs = [r[i] for r in regions]
            if len(set(rs)) == 1:
                cons.append(rs[0])
            else:
                cons.append(rs)
        return cons

<<<<<<< HEAD
def _resolve_ambiguities(annotated_by_loc_by_tx_by_gene, chrom_order, collapse_exons=True,
                         high_confidence=False, short=False, output_features=False):
    def _collapse(regions):
        cons = []
        for i in ga.BedCols.cols:
            rs = [r[i] for r in regions]
            if len(set(rs)) == 1:
                cons.append(rs[0])
            else:
                cons.append(rs)
        return cons

    annotated = []
    for (chrom, start, end), overlaps_by_tx_by_gene in annotated_by_loc_by_tx_by_gene.iteritems():
        features = dict()
        for gname, overlaps_by_tx in overlaps_by_tx_by_gene.iteritems():
            tx_overlaps = [os for os in overlaps_by_tx for o in os if o[ga.BedCols.FEATURE] == 'transcript']
            tsl1_overlaps = [o for o in tx_overlaps if o[ga.BedCols.TSL] in ['1', 'NA']]
            if tsl1_overlaps:
                tx_overlaps = tsl1_overlaps
            hugo_overlaps = [o for o in tx_overlaps if o[ga.BedCols.HUGO]]
            if hugo_overlaps:
                tx_overlaps = hugo_overlaps
            best_overlap = max(sorted(tx_overlaps, key=lambda fs: int(fs[ga.BedCols.END]) - int(fs[ga.BedCols.START])))
            tx_id = best_overlap[ga.BedCols.ENSEMBL_ID]

            for overlaps in overlaps_by_tx[tx_id]:
=======
    annotated = []
    for (chrom, start, end), overlaps_by_tx in annotated_by_loc_by_gene.iteritems():
        for tx_name, overlaps in overlaps_by_tx.iteritems():
            if collapse_exons:
                # consensus = _collapse([_r for _r, _, _ in overlaps if _r[ga.BedCols.FEATURE] != 'transcript'])
>>>>>>> 83af7a28d38c14e9c1287475e1cbdc985eb836ce
                consensus = [None for _ in ga.BedCols.cols]
                consensus[:3] = chrom, start, end
                consensus[ga.BedCols.FEATURE] = 'capture'
                consensus[ga.BedCols.EXON_OVERLAPS_BASES] = 0
                consensus[ga.BedCols.EXON_OVERLAPS_PERCENTAGE] = 0
                consensus[ga.BedCols.EXON] = set()

<<<<<<< HEAD
                for fields, overlap_size in overlaps:
                    c_overlap_bp = overlap_size
                    c_overlap_pct = 100.0 * c_overlap_bp / (int(end) - int(start))
                    if high_confidence and c_overlap_pct < 50.0:
                        continue

                    if output_features:
                        f_start = int(fields[1])
                        f_end = int(fields[2])
                        feature = features.get((f_start, f_end))
                        if feature is None:
                            feature = [None for _ in ga.BedCols.cols]
                            feature[:len(fields)] = fields
                            # feature[ga.BedCols.TX_OVERLAP_BASES] = 0
                            # feature[ga.BedCols.TX_OVERLAP_PERCENTAGE] = 0
                            features[(f_start, f_end)] = feature
                        # feature[ga.BedCols.TX_OVERLAP_BASES] += c_overlap_bp
                        # feature[ga.BedCols.TX_OVERLAP_PERCENTAGE] += 100.0 * c_overlap_bp / (int(f_end) - int(f_start))
                        # TODO: don't forget to merge BED if not

=======
                for fields, overlap_bp, overlap_percentage in overlaps:
>>>>>>> 83af7a28d38c14e9c1287475e1cbdc985eb836ce
                    if fields[ga.BedCols.FEATURE] == 'transcript':
                        consensus[ga.BedCols.GENE] = fields[ga.BedCols.GENE]
                        consensus[ga.BedCols.STRAND] = fields[ga.BedCols.STRAND]
                        consensus[ga.BedCols.BIOTYPE] = fields[ga.BedCols.BIOTYPE]
                        consensus[ga.BedCols.ENSEMBL_ID] = fields[ga.BedCols.ENSEMBL_ID]
                        consensus[ga.BedCols.TSL] = fields[ga.BedCols.TSL]
                        consensus[ga.BedCols.HUGO] = fields[ga.BedCols.HUGO]
<<<<<<< HEAD
                        consensus[ga.BedCols.TX_OVERLAP_BASES] = c_overlap_bp
                        consensus[ga.BedCols.TX_OVERLAP_PERCENTAGE] = c_overlap_pct
                    elif fields[ga.BedCols.FEATURE] == 'exon':
                        consensus[ga.BedCols.EXON_OVERLAPS_BASES] += c_overlap_bp
                        consensus[ga.BedCols.EXON_OVERLAPS_PERCENTAGE] += c_overlap_pct
=======
                        consensus[ga.BedCols.TX_OVERLAP_BASES] = overlap_bp
                        consensus[ga.BedCols.TX_OVERLAP_PERCENTAGE] = overlap_percentage
                    elif fields[ga.BedCols.FEATURE] == 'exon':
                        consensus[ga.BedCols.EXON_OVERLAPS_BASES] += overlap_bp
                        consensus[ga.BedCols.EXON_OVERLAPS_PERCENTAGE] += overlap_percentage
>>>>>>> 83af7a28d38c14e9c1287475e1cbdc985eb836ce
                        consensus[ga.BedCols.EXON].add(int(fields[ga.BedCols.EXON]))
                consensus[ga.BedCols.EXON] = sorted(list(consensus[ga.BedCols.EXON]))

                annotated.append(consensus)

<<<<<<< HEAD
        if output_features:
            annotated.extend(sorted(features.values(), key=get_sort_key(chrom_order)))
=======
            else:
                for fields, overlap_bp, overlap_percentage in overlaps:
                    fields[:3] = chrom, start, end
                    fields[ga.BedCols.TX_OVERLAP_BASES] = overlap_bp
                    fields[ga.BedCols.TX_OVERLAP_PERCENTAGE] = overlap_percentage
                    annotated.append(fields)
>>>>>>> 83af7a28d38c14e9c1287475e1cbdc985eb836ce

    return annotated


<<<<<<< HEAD
def _annotate(bed, ref_bed, chr_order, fai_fpath=None, high_confidence=False, short=False,
              collapse_exons=True, output_features=False):
=======
def _annotate(bed, ref_bed, chr_order, fai_fpath=None, high_confidence=False, short=False, collapse_exons=True):
>>>>>>> 83af7a28d38c14e9c1287475e1cbdc985eb836ce
    # if genome:
        # genome_fpath = cut(fai_fpath, 2, output_fpath=intermediate_fname(work_dir, fai_fpath, 'cut2'))
        # intersection = bed.intersect(ref_bed, sorted=True, wao=True, g='<(cut -f1,2 ' + fai_fpath + ')')
        # intersection = bed.intersect(ref_bed, sorted=True, wao=True, genome=genome.split('-')[0])
    # else:

    if fai_fpath and count_bed_cols(fai_fpath) == 2:
        intersection = bed.intersect(ref_bed, wao=True, sorted=True, g=fai_fpath)
    else:
        intersection = bed.intersect(ref_bed, wao=True)

    total_annotated = 0
    total_uniq_annotated = 0

    met = set()

    annotated_by_loc_by_tx_by_gene = OrderedDefaultDict(lambda: OrderedDefaultDict(lambda: defaultdict(list)))
    off_targets = list()

    for intersection_fs in intersection:
        if len(intersection_fs) < 3 + len(ga.BedCols.cols) + 1:
            critical('Cannot parse the reference BED file - unexpected number of lines '
                     '(' + str(len(intersection_fs)) + ') in ' + str(intersection_fs))

        a_chr, a_start, a_end, e_chr = intersection_fs[:4]
        overlap_size = int(intersection_fs[-1])
        assert e_chr == '.' or a_chr == e_chr, str((a_chr + ', ' + e_chr))

        if e_chr == '.':
            fs = [None for _ in ga.BedCols.cols]
            fs[:3] = [a_chr, a_start, a_end]
            off_targets.append(fs)
        else:
            fs = [None for _ in ga.BedCols.cols]
            fs[:3] = [a_chr, a_start, a_end]
            fs[3:len(intersection_fs[6:-1])] = intersection_fs[6:-1]
            total_annotated += 1
            if (a_chr, a_start, a_end) not in met:
                total_uniq_annotated += 1
                met.add((a_chr, a_start, a_end))

<<<<<<< HEAD
            annotated_by_loc_by_tx_by_gene[(a_chr, int(a_start), int(a_end))][
                intersection_fs[3:ga.BedCols.GENE]
            ][
                intersection_fs[3+ga.BedCols.ENSEMBL_ID]
            ].append((intersection_fs[3:-1], overlap_size))
=======
            annotated_by_loc_by_gene[(a_chr, int(a_start), int(a_end))][fs[ga.BedCols.ENSEMBL_ID]].append(
                (fs, overlap_size, 100.0 * overlap_size / (int(a_end) - int(a_start))),
            )
>>>>>>> 83af7a28d38c14e9c1287475e1cbdc985eb836ce

    info('  Total annotated regions: ' + str(total_annotated))
    info('  Total unique annotated regions: ' + str(total_uniq_annotated))
    info('  Total off target regions: ' + str(len(off_targets)))
    info('Resolving ambiguities...')
<<<<<<< HEAD
    annotated = _resolve_ambiguities(annotated_by_loc_by_tx_by_gene, chr_order, collapse_exons,
         high_confidence, short, output_features)
=======
    annotated = _resolve_ambiguities(annotated_by_loc_by_gene, chr_order, collapse_exons)
>>>>>>> 83af7a28d38c14e9c1287475e1cbdc985eb836ce

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
