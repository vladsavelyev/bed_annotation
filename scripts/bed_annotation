#!/usr/bin/env python

import click as click
import os
import shutil
from os.path import join, basename
from tempfile import mkdtemp

import ensembl as ebl
from ensembl.bed_annotation import annotate
from ngs_utils.bed_utils import verify_bed, clean_bed
from ngs_utils import logger
from ngs_utils.file_utils import adjust_path, safe_mkdir, verify_file
from ngs_utils.logger import info
from ngs_utils.logger import debug


@click.command()
@click.argument('input_bed', type=click.Path(exists=True))
@click.option('-o', '--output-file', type=click.Path(), help='Output file path')
@click.option('--output-features', is_flag=True, help='Also output featues that was used to annotate')
@click.option('-g', 'genome', default='GRCh37', type=click.Choice(ebl.SUPPORTED_GENOMES), help='Genome build')
@click.option('--canonical', '--only-canonical', is_flag=True, help='Use only features from canonical transcripts to annotate')
@click.option('--short', is_flag=True, help='Add only "Gene" column (outputa 4-col BED file instead of 6-col)')
@click.option('-e', '--extended', is_flag=True, help='Output additional columns: transcript, GC, overlap size...')
@click.option('--high-confidence', is_flag=True, help='Annotate with only high confidence regions (TSL is 1 or NA, with HUGO symbol, total overlap size > 50%)')
@click.option('-a', '--ambiguities', '--ambiguities-method', type=click.Choice(['best_one', 'best_ask', 'best_all', 'all_ask', 'all']), default='best_all',
              help='How to resolve ambuguios overlaps with reliable transcripts for a single region')
@click.option('--coding-only', is_flag=True, help='Use only protein coding genes to annotate')
@click.option('--collapse-exons', is_flag=True)
@click.option('--work-dir', default=None, type=click.Path())
@click.option('-d', '--debug', '--is-debug', is_flag=True)
def main(input_bed, output_file, output_features=False, genome=None,
         only_canonical=False, short=False, extended=False, high_confidence=False,
         ambiguities_method=False, coding_only=False, collapse_exons=False, work_dir=False, is_debug=False):
    """ Annotating BED file based on reference features annotations.
    """
    logger.init(is_debug_=is_debug)

    if not genome:
        raise click.BadParameter('Error: please, specify genome build name with -g (e.g. `-g hg19`)', param='genome')

    if short:
        if extended:        raise click.BadParameter('--short and --extended can\'t be set both', param='extended')
        if output_features: raise click.BadParameter('--short and --output-features can\'t be set both', param='output_features')
    elif output_features or extended:
        extended = True
        short    = False

    if not verify_file(input_bed):
        click.BadParameter(f'Usage: {__file__} Input_BED_file -g hg19 -o Annotated_BED_file [options]', param='input_bed')
    input_bed = verify_file(input_bed, is_critical=True, description=f'Input BED file for {__file__}')

    if work_dir:
        work_dir = join(adjust_path(work_dir), os.path.splitext(basename(input_bed))[0])
        safe_mkdir(work_dir)
        info(f'Created work directory {work_dir}')
    else:
        work_dir = mkdtemp('bed_annotate')
        debug('Created temporary work directory {work_dir}')

    input_bed = clean_bed(input_bed, work_dir)
    input_bed = verify_bed(input_bed, is_critical=True, description=f'Input BED file for {__file__} after cleaning')

    output_file = adjust_path(output_file)

    output_file = annotate(
        input_bed, output_file, work_dir, genome=genome,
        only_canonical=only_canonical, short=short, extended=extended,
        high_confidence=high_confidence, collapse_exons=collapse_exons,
        output_features=output_features,
        ambiguities_method=ambiguities_method, coding_only=coding_only,
        is_debug=is_debug)

    if not work_dir:
        debug(f'Removing work directory {work_dir}')
        shutil.rmtree(work_dir)

    info(f'Done, saved to {output_file}')


if __name__ == '__main__':
    main()



# TODO: prefer consecutive annotations
