from os.path import dirname, join, abspath, isfile, basename, splitext

import sys

from pybedtools import BedTool

SUPPORTED_GENOMES = ['hg19', 'hg19-noalt', 'hg38', 'hg38-noalt']
REFSEQ_DIR = 'RefSeq'


def check_genome(genome):
    if genome not in SUPPORTED_GENOMES:
        sys.stdout.write('Genome ' + genome + ' is not supported. Supported genomes: ' + ', '.join(SUPPORTED_GENOMES) + '\n')
        sys.exit(1)

def get_all_features(genome, gzip=True):
    return _get_refseq('all_features.{genome}.bed' + ('.gz' if gzip else ''), genome)

# ncRNA and protein coding CDS, Exons, Gene and Transcript features - only for canonical tracnscripts
# - used to annotate BED files
# - used to report in regional TargQC reports
def get_all_features_canonical(genome, gzip=True):
    return _get_refseq('all_features.{genome}.canon.bed' + ('.gz' if gzip else ''), genome)

# CDS for canonical transcripts, used as TargQC target when no BED available
def get_cds(genome):
    return _get_refseq('CDS.{genome}.bed', genome)

# CDS for canonical transcripts, one transcripts only for a gene. Used for Seq2C CNV calling
def get_seq2c_cds(genome):
    return _get_refseq('Seq2C_CDS.{genome}.bed', genome)

# TODO
# ?split features into smaller parts:
# hg19
#   CDS.canon.bed
#   CDS.bed
#   Exons.canon.bed
#   Exons.bed
#   Transcripts.canon.bed
#   Transcripts.bed
#   Genes.bed
#   CDS_miRNA.all_features.canon.bed

def get_refseq_knowngene(genome):
    return _get_refseq('RefSeq_knownGene.{genome}.txt', genome.split('-')[0])  # no -alt

def get_refseq_gene(genome):
    return _get_refseq('refGene.{genome}.txt.gz', genome.split('-')[0])  # no -alt

def get_hgnc_gene_synonyms():
    return _get_refseq('HGNC_gene_synonyms.txt')


def _get(relative_path, genome=None):
    """
    :param relative_path: relative path of the file inside the repository
    :param genome: genome name. Can contain chromosome name after comma, like hg19-chr20,
                   in case of BED, the returning BedTool will be with added filter.
    :return: BED fpath
    """
    chrom = None
    if genome:
        if '-chr' in genome:
            genome, chrom = genome.split('-')
        check_genome(genome)
        relative_path = relative_path.format(genome=genome)

    path = abspath(join(dirname(__file__), relative_path))

    if path.endswith('.bed') or path.endswith('.bed.gz'):
        if chrom:
            bed = BedTool(path)
            bed = bed.filter(lambda r: r.chrom == chrom)
            path = bed.saveas().fn
    return path

def get_refseq_dirpath():
    return abspath(join(dirname(__file__), REFSEQ_DIR))

def _get_refseq(fname, genome=None):
    return _get(join(REFSEQ_DIR, fname), genome)