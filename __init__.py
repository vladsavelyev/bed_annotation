from os.path import dirname, join, abspath

from Utils.logger import critical


SUPPORTED_GENOMES = ['hg19', 'hg19-noalt', 'hg38', 'hg38-noalt', 'hg19-chr20']
REFSEQ_DIR = 'RefSeq'


def check_genome(genome):
    if genome not in SUPPORTED_GENOMES:
        critical('Genome ' + genome + ' is not supported. Supported genomes: ' + ', '.join(SUPPORTED_GENOMES))

def get_all_features(genome):
    return _get_refseq('all_features.{genome}.bed', genome)

# ncRNA and protein coding CDS, Exons, Gene and Transcript features - only for canonical tracnscripts
# - used to annotate BED files
# - used to report in regional TargQC reports
def get_all_features_canonical(genome):
    return _get_refseq('all_features.{genome}.canon.bed', genome)

# CDS for canonical transcripts, used as TargQC target when no BED available
def get_cds(genome):
    return _get_refseq('CDS.{genome}.bed', genome)

# CDS for canonical transcripts, one transcripts only for a gene. Used for Seq2C CNV calling
def get_seq2c_cds(genome):
    return _get_refseq('Seq2C_CDS.{genome}.bed', genome)

    # features:                 /ngs/reference_data/genomes/Hsapiens/hg19/bed/Exons/RefSeq/
    # bed_annotation_features:  /ngs/reference_data/genomes/Hsapiens/hg19/bed/Exons/RefSeq/RefSeq.all_features.hg19.bed
    # cds:                      /ngs/reference_data/genomes/Hsapiens/hg19/bed/Exons/RefSeq/RefSeq_CDS.hg19.bed

# def get_features(genome, canonical=True, features=None, biotypes=None):
#     features = features or ['CDS', 'Exon', 'Transcript', 'Gene']
#
#     if biotypes:
#         biotypes = 'CDS'
#     fname = 'RefSeq' + ''
#     return _get(genome, join('RefSeq', 'RefSeq_CDS_miRNA.all_features.{genome}.canon.bed'))


# TODO
# split features into smaller parts:
#
# hg19
#   CDS.canon.bed
#   CDS.bed
#   Exons.canon.bed
#   Exons.bed
#   Transcripts.canon.bed
#   Transcripts.bed
#   Genes.bed
#   CDS_miRNA.all_features.canon.bed

def get_hgnc_gene_synonyms():
    return _get_refseq('HGNC_gene_synonyms.txt')

def get_refseq_knowngene(genome):
    return _get_refseq('RefSeq_knownGene.{genome}.txt', genome.split('-')[0])


def _get(relative_path, genome=None):
    if genome:
        check_genome(genome)
        relative_path = relative_path.format(genome=genome)
    return abspath(join(dirname(__file__), relative_path))

def get_refseq_dirpath():
    return abspath(join(dirname(__file__), REFSEQ_DIR))

def _get_refseq(fname, genome=None):
    return _get(join(REFSEQ_DIR, fname), genome)