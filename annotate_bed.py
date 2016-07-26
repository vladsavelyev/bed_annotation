#!/usr/bin/env python
import os
from itertools import izip, count
from optparse import OptionParser, SUPPRESS_HELP
from collections import defaultdict, OrderedDict
from os.path import isfile, join, basename

from datetime import datetime
from pybedtools import BedTool

import GeneAnnotation as ga
from Utils import reference_data
from Utils.logger import warn, debug
from Utils.utils import OrderedDefaultDict
from Utils.bed_utils import verify_bed, bedtools_version, SortableByChrom, cut, count_bed_cols
from Utils.file_utils import verify_file, file_transaction, open_gzipsafe, which, intermediate_fname, adjust_path, \
    safe_mkdir
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
            action='store_true',
            default=False,
            help='Add only "Gene" column',
        )),
        (['--extended'], dict(
            dest='extended',
            action='store_true',
            default=False,
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

    parser = OptionParser(description='Annotating BED file based on reference features annotations.')
    for args, kwargs in options:
        parser.add_option(*args, **kwargs)
    opts, args = parser.parse_args()
    logger.is_debug = opts.debug

    if opts.short:
        if opts.extended:        critical('--short and --extended can\'t be set both')
        if opts.output_features: critical('--short and --output-features can\'t be set both')
    elif opts.output_features or opts.extended:
        opts.extended = True
        opts.short = False
    if opts.seq2c:
        opts.short = opts.high_confidence = opts.only_canonical = True
        opts.output_features = opts.extended = False

    if len(args) < 1:
        parser.exit('Usage: ' + __file__ + ' Input_BED_file -g hg19 -o Annotated_BED_file [--canonical]')
    input_bed_fpath = verify_bed(args[0], is_critical=True, description='Input BED file for ' + __file__)

    output_fpath = adjust_path(opts.output_file)

    work_dir = None
    prev_output_fpath = None
    if opts.debug:
        if isfile(output_fpath):
            prev_output_fpath = output_fpath + '_' + datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
            os.rename(output_fpath, prev_output_fpath)
        if opts.work_dir:
            work_dir = safe_mkdir(join(adjust_path(opts.work_dir), basename(input_bed_fpath)))

    output_fpath = annotate(
        input_bed_fpath, output_fpath, opts.genome, work_dir=work_dir, is_debug=opts.debug,
        only_canonical=opts.only_canonical, short=opts.short, extended=opts.extended, high_confidence=opts.high_confidence,
        collapse_exons=opts.collapse_exons, output_features=opts.output_features)

    if opts.debug:
        if prev_output_fpath:
            os.system('diff ' + prev_output_fpath + ' ' + output_fpath)
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


def get_sort_key(chr_order):
    return lambda fs: (chr_order[fs[ga.BedCols.CHROM]],
            int(fs[ga.BedCols.START]),
            int(fs[ga.BedCols.END]),
            0 if fs[ga.BedCols.FEATURE] == 'transcript' else 1,
            fs[ga.BedCols.ENSEMBL_ID],
            fs[ga.BedCols.GENE])


def annotate(input_bed_fpath, output_fpath, genome=None,
             work_dir=None, is_debug=False, reuse=False,
             only_canonical=False, short=False, extended=False, high_confidence=False,
             collapse_exons=True, output_features=False):

    if reuse and isfile(output_fpath) and verify_file(output_fpath):
        debug(output_fpath + ' exists, reusing.')
        return output_fpath

    debug('Getting features from storage')
    features_bed = ga.get_all_features(genome)
    if features_bed is None:
        critical('Genome ' + genome + ' is not supported. Supported: ' + ', '.join(ga.SUPPORTED_GENOMES))

    if genome:
        fai_fpath = reference_data.get_fai(genome)
        chr_order = reference_data.get_chrom_order(genome)
    else:
        fai_fpath = None
        chr_order = bed_chrom_order(input_bed_fpath)

    bed = BedTool(input_bed_fpath).cut([0, 1, 2])

    # features_bed = features_bed.saveas()
    # cols = features_bed.field_count()
    # if cols < 12:
    #     features_bed = features_bed.each(lambda f: f + ['.']*(12-cols))
    if high_confidence:
        features_bed = features_bed.filter(ga.high_confidence_filter)
    if work_dir and is_debug:
        features_bed = features_bed.saveas(join(work_dir, 'features_tmp.bed'))
        print features_bed.fn
    # unique_tx_by_gene = find_best_tx_by_gene(features_bed)

    info('Extracting features')
    features_bed = features_bed.filter(lambda x:
        x[ga.BedCols.FEATURE] in ['exon', 'CDS', 'stop_codon', 'transcript'])
        # x[ga.BedCols.ENSEMBL_ID] == unique_tx_by_gene[x[ga.BedCols.GENE]])

    info('Annotating...')
    annotated = _annotate(bed.saveas(join(work_dir, 'bed.bed')), features_bed.saveas(join(work_dir, 'features.bed')), chr_order, fai_fpath,
        high_confidence, collapse_exons, output_features, is_debug, work_dir)

    header = [ga.BedCols.names[i] for i in ga.BedCols.cols]

    info('Saving annotated regions...')
    with file_transaction(None, output_fpath) as tx:
        with open(tx, 'w') as out:
            if short:
                header = header[:4]
            if not extended:
                header = header[:6]
            out.write('\t'.join(header) + '\n')
            for fields in annotated:
                if short:
                    fields = fields[:4]
                if not extended:
                    fields = fields[:6]
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
        return ', '.join(_format_field(v) for v in value) or '.'
    elif isinstance(value, float):
        return '{:.1f}'.format(value)
    else:
        return str(value or '.')


# def find_best_tx_by_gene(features_bed):
#     all_transcripts = features_bed.filter(lambda x: x[ga.BedCols.FEATURE] in ['transcript'])
#
#     tx_by_gene = defaultdict(list)
#     for tx_region in all_transcripts:
#         tx_by_gene[tx_region.name].append(tx_region)
#     unique_tx_by_gene = dict()
#     for g, txs in tx_by_gene.iteritems():
#         tsl1_txs = [tx for tx in txs if tx[ga.BedCols.TSL] in ['1', 'NA']]
#         if tsl1_txs:
#             txs = tsl1_txs[:]
#         hugo_txs = [tx for tx in txs if tx[ga.BedCols.HUGO]]
#         if hugo_txs:
#             txs = hugo_txs[:]
#         best_tx = max(txs, key=lambda tx_: int(tx_[ga.BedCols.END]) - int(tx_[ga.BedCols.START]))
#         unique_tx_by_gene[g] = best_tx[ga.BedCols.ENSEMBL_ID]
#     return unique_tx_by_gene


def tx_sort_key(x):
    tsl_key = {'1': 0, '2': 2, '3': 3, '4': 4, '5': 5}.get(x[ga.BedCols.TSL], 1)
    hugo_key = 0 if x[ga.BedCols.HUGO] not in ['.', '', None] else 1
    length_key = -(int(x[ga.BedCols.END]) - int(x[ga.BedCols.START]))
    return tsl_key, hugo_key, length_key


# def select_best_tx(overlaps_by_tx):
#     tx_overlaps = [(o, size) for os in overlaps_by_tx.values() for (o, size) in os if o[ga.BedCols.FEATURE] == 'transcript']
#     tsl1_overlaps = [(o, size) for (o, size) in tx_overlaps if o[ga.BedCols.TSL] in ['1', 'NA']]
#     if tsl1_overlaps:
#         tx_overlaps = tsl1_overlaps
#     hugo_overlaps = [(o, size) for (o, size) in tx_overlaps if o[ga.BedCols.HUGO]]
#     if hugo_overlaps:
#         tx_overlaps = hugo_overlaps
#     (best_overlap, size) = max(sorted(tx_overlaps, key=lambda (o, size): int(o[ga.BedCols.END]) - int(o[ga.BedCols.START])))
#     tx_id = best_overlap[ga.BedCols.ENSEMBL_ID]
#     return tx_id


def _resolve_ambiguities(annotated_by_tx_by_gene_by_loc, chrom_order, collapse_exons=True,
                         high_confidence=False, output_features=False):
    # sort transcripts by "quality" and then select which tx id to report in case of several overlaps
    # (will be useful further when interating exons)
    # annotated_by_tx_by_gene_by_loc = OrderedDict([
    #     (loc, OrderedDict([
    #         (gname, sorted(annotated_by_tx.itervalues(), key=tx_sort_key))
    #             for gname, annotated_by_tx in annotated_by_tx_by_gene.iteritems()
    #     ])) for loc, annotated_by_tx_by_gene in annotated_by_tx_by_gene_by_loc.iteritems()
    # ])

    annotated = []
    for (chrom, start, end), overlaps_by_tx_by_gene in annotated_by_tx_by_gene_by_loc.iteritems():
        features = dict()
        for gname, overlaps_by_tx in overlaps_by_tx_by_gene.iteritems():
            consensus = [None for _ in ga.BedCols.cols]
            consensus[:3] = chrom, start, end
            consensus[ga.BedCols.FEATURE] = 'capture'
            # consensus[ga.BedCols.EXON_OVERLAPS_BASES] = 0
            consensus[ga.BedCols.EXON_OVERLAPS_PERCENTAGE] = 0
            consensus[ga.BedCols.EXON] = set()

            if not gname:
                annotated.append(consensus)
                continue

            # if high_confidence:
            #     confident_overlaps = (x for xx in overlaps_by_tx.itervalues()
            #                           for x, overlap_bp in xx if (100.0 * overlap_bp / (int(end) - int(start))) >= 50.0)
            # else:
            #     confident_overlaps = (x for xx in overlaps_by_tx.itervalues() for x, _ in xx)

            # Choosing confident transcript overlaps
            all_tx = (x for xx in overlaps_by_tx.itervalues() for x, overlap_bp in xx
                       if x[ga.BedCols.FEATURE] == 'transcript' and
                       (not high_confidence or (100.0 * overlap_bp / (int(end) - int(start))) >= 50.0))
            tx_sorted_list = [x[ga.BedCols.ENSEMBL_ID] for x in sorted(all_tx, key=tx_sort_key)]
            if not tx_sorted_list:
                continue
            tx_id = tx_sorted_list[0]
            overlaps = overlaps_by_tx[tx_id]

            # # sort transcripts by "quality" and then select which tx id to report in case of several overlaps
            # # (will be useful further when interating exons)
            # annotated_by_tx_by_gene_by_loc = OrderedDict([
            #     (loc, OrderedDict([
            #         (gname, sorted(annotated_by_tx.itervalues(), key=tx_sort_key))
            #             for gname, annotated_by_tx in annotated_by_tx_by_gene.iteritems()
            #     ])) for loc, annotated_by_tx_by_gene in annotated_by_tx_by_gene_by_loc.iteritems()
            # ])

            for fields, overlap_size in overlaps:
                c_overlap_bp = overlap_size
                c_overlap_pct = 100.0 * c_overlap_bp / (int(end) - int(start))

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

                if fields[ga.BedCols.FEATURE] == 'transcript':
                    consensus[ga.BedCols.GENE] = gname
                    consensus[ga.BedCols.STRAND] = fields[ga.BedCols.STRAND]
                    consensus[ga.BedCols.BIOTYPE] = fields[ga.BedCols.BIOTYPE]
                    consensus[ga.BedCols.ENSEMBL_ID] = fields[ga.BedCols.ENSEMBL_ID]
                    consensus[ga.BedCols.TSL] = fields[ga.BedCols.TSL]
                    consensus[ga.BedCols.HUGO] = fields[ga.BedCols.HUGO]
                    # consensus[ga.BedCols.TX_OVERLAP_BASES] = c_overlap_bp
                    consensus[ga.BedCols.TX_OVERLAP_PERCENTAGE] = c_overlap_pct
                elif fields[ga.BedCols.FEATURE] == 'exon':
                    # consensus[ga.BedCols.EXON_OVERLAPS_BASES] += c_overlap_bp
                    consensus[ga.BedCols.EXON_OVERLAPS_PERCENTAGE] += c_overlap_pct
                    consensus[ga.BedCols.EXON].add(int(fields[ga.BedCols.EXON]))
            consensus[ga.BedCols.EXON] = sorted(list(consensus[ga.BedCols.EXON]))

            annotated.append(consensus)

        features = sorted(features.values(), key=get_sort_key(chrom_order))
        if output_features:
            annotated.extend(features)
    return annotated


def _annotate(bed, ref_bed, chr_order, fai_fpath=None, high_confidence=False,
              collapse_exons=True, output_features=False, is_debug=False, work_dir=None):
    # if genome:
        # genome_fpath = cut(fai_fpath, 2, output_fpath=intermediate_fname(work_dir, fai_fpath, 'cut2'))
        # intersection = bed.intersect(ref_bed, sorted=True, wao=True, g='<(cut -f1,2 ' + fai_fpath + ')')
        # intersection = bed.intersect(ref_bed, sorted=True, wao=True, genome=genome.split('-')[0])
    # else:

    intersection = None
    intersection_fpath = None
    if work_dir and debug:
        intersection_fpath = join(work_dir, 'intersection.bed')
        if isfile(intersection_fpath):
            info('Loading from ' + intersection_fpath)
            intersection = BedTool(intersection_fpath)
    if not intersection:
        if fai_fpath and count_bed_cols(fai_fpath) == 2:
            intersection = bed.intersect(ref_bed, wao=True, sorted=True, g=fai_fpath)
        else:
            intersection = bed.intersect(ref_bed, wao=True)
    if work_dir and debug and not isfile(intersection_fpath):
        intersection.saveas(intersection_fpath)
        info('Saved to ' + intersection_fpath)

    total_annotated = 0
    total_uniq_annotated = 0
    total_off_target = 0

    met = set()

    annotated_by_tx_by_gene_by_loc = OrderedDefaultDict(lambda: OrderedDefaultDict(lambda: defaultdict(list)))
    # off_targets = list()

    for intersection_fs in intersection:
        if len(intersection_fs) < 3 + len(ga.BedCols.cols) + 1:
            critical('Cannot parse the reference BED file - unexpected number of lines '
                     '(' + str(len(intersection_fs)) + ') in ' + str(intersection_fs))

        a_chr, a_start, a_end, e_chr = intersection_fs[:4]
        overlap_size = int(intersection_fs[-1])
        assert e_chr == '.' or a_chr == e_chr, str((a_chr + ', ' + e_chr))

        fs = [None for _ in ga.BedCols.cols]
        fs[:3] = [a_chr, a_start, a_end]
        reg = (a_chr, int(a_start), int(a_end))

        if e_chr == '.':
            total_off_target += 1
            # off_targets.append(fs)
            annotated_by_tx_by_gene_by_loc[reg][''][''].append(([], 0))
        else:
            fs[3:len(intersection_fs[6:-1])] = intersection_fs[6:-1]
            total_annotated += 1
            if (a_chr, a_start, a_end) not in met:
                total_uniq_annotated += 1
                met.add((a_chr, a_start, a_end))

            gene = intersection_fs[3 + ga.BedCols.GENE] if not high_confidence else intersection_fs[3 + ga.BedCols.HUGO]
            tx = intersection_fs[3 + ga.BedCols.ENSEMBL_ID]
            annotated_by_tx_by_gene_by_loc[reg][gene][tx].append((intersection_fs[3:-1], overlap_size))

    info('  Total annotated regions: ' + str(total_annotated))
    info('  Total unique annotated regions: ' + str(total_uniq_annotated))
    info('  Total off target regions: ' + str(total_off_target))
    info('Resolving ambiguities...')
    annotated = _resolve_ambiguities(annotated_by_tx_by_gene_by_loc, chr_order, collapse_exons, high_confidence, output_features)

    return annotated


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
