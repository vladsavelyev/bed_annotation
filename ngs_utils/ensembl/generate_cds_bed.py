#!/usr/bin/env python

import ngs_utils.ensembl as ebl
import os
import shutil
from optparse import OptionParser, SUPPRESS_HELP
from os.path import isfile, join, basename, dirname, pardir
from ngs_utils import logger
from ngs_utils.file_utils import file_transaction, adjust_path, safe_mkdir, verify_file


''' Generates coding_regions BED file

    Example usage: 
    python {__file__} -g GRCh37 --canonical | grep -v ^MT | grep -v ^GL | sort -k1,1V -k2,2n | bedtools merge -i - > coding_regions.canonical.clean.sort.merged.bed
'''

def main():
    options = [
        (['-g', '--genome'], dict(
            dest='genome',
            help='Genome build. Accepted values: ' + ', '.join(ebl.SUPPORTED_GENOMES),
        )),
        (['-c', '--canonical'], dict(
            dest='canonical',
            action='store_true',
            help='Use canonical only',
        )),
    ]
    parser = OptionParser()
    for args, kwargs in options:
        parser.add_option(*args, **kwargs)
    opts, args = parser.parse_args()

    if not opts.genome:
        logger.critical('Error: please, specify genome build name with -g (e.g. `-g hg19`)')
    genome = opts.genome

    logger.debug('Getting features from storage')
    features_bed = ebl.get_all_features(genome)
    if features_bed is None:
        logger.critical('Genome ' + genome + ' is not supported. Supported: ' + ', '.join(ebl.SUPPORTED_GENOMES))

    logger.warn('Extracting features from Ensembl GTF')
    features_bed = features_bed.filter(lambda x: x[ebl.BedCols.FEATURE] == 'CDS')
    if opts.canonical:
        features_bed = features_bed.filter(ebl.get_only_canonical_filter(genome))

    logger.warn('Saving CDS regions...')
    output_fpath = adjust_path(join(dirname(__file__), pardir, genome, 'bed', 'CDS-canonical.bed'))
    with file_transaction(None, output_fpath) as tx:
        features_bed.cut(range(6)).saveas(tx)
    logger.warn('Done, saved to ' + output_fpath)


if __name__ == '__main__':
    main()
