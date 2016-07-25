from os.path import dirname, join, abspath, isfile, basename, splitext
import sys
from pybedtools import BedTool

from Utils.bed_utils import bedtools_version
from Utils.file_utils import which, open_gzipsafe


class BedCols:
    CHROM, \
    START, \
    END, \
    GENE, \
    EXON, \
    STRAND, \
    FEATURE, \
    BIOTYPE, \
    ENSEMBL_ID, \
    TSL, \
    HUGO, \
    TX_OVERLAP_BASES, \
    TX_OVERLAP_PERCENTAGE, \
    EXON_OVERLAPS_BASES, \
    EXON_OVERLAPS_PERCENTAGE \
        = cols = range(15)

    names = {
        CHROM: '#Chrom',
        START: 'Start',
        END: 'End',
        GENE: 'Gene',
        EXON: 'Exon',
        STRAND: 'Strand',
        FEATURE: 'Feature',
        BIOTYPE: 'Biotype',
        ENSEMBL_ID: 'Ensembl_ID',
        TSL: 'TSL',
        HUGO: 'HUGO',
        TX_OVERLAP_BASES: 'Tx overlap bp',
        TX_OVERLAP_PERCENTAGE: 'Tx overlap %',
        EXON_OVERLAPS_BASES: 'Exon overlaps bp',
        EXON_OVERLAPS_PERCENTAGE: 'Exon overlaps %',
    }

SUPPORTED_GENOMES = ['hg19', 'hg19-noalt', 'hg38', 'hg38-noalt']

def check_genome(genome):
    if genome not in SUPPORTED_GENOMES:
        sys.stdout.write('Genome ' + genome + ' is not supported. Supported genomes: ' + ', '.join(SUPPORTED_GENOMES) + '\n')
        sys.exit(1)


REFSEQ_DIR = 'RefSeq'

def refseq_knowngene_fpath(genome):
    return _get_refseq('RefSeq_knownGene.txt', genome.split('-')[0])  # no -alt

def get_refseq_gene(genome):
    return _get_refseq('refGene.txt.gz', genome.split('-')[0])  # no -alt

def _get_refseq(fname, genome=None):
    return _get(join(REFSEQ_DIR, genome, fname), genome)

def get_refseq_dirpath():
    return abspath(join(dirname(__file__), REFSEQ_DIR))


ENSEMBL_DIR = 'Ensembl'

def _get_ensembl(fname, genome=None):
    if genome:
        return _get(join(ENSEMBL_DIR, genome.split('-')[0], fname), genome)
    else:
        return _get(join(ENSEMBL_DIR, fname))

'''
This repository is made for storing genomic features coordinates and annotations.

Annotation in TargQC:
  - prepocess: generate BED file from ref-transcripts.gtf.db:
    - annotate with RefSeq IDs and TSL

  - annotate by priority:
    - use all RefSeq exons (TODO: check if need to separate processed_transcript)
    - use all RefSeq transcripts
    - use other exons
    - use other transcripts

Annotation in bcbio:
  - generate BED file from ref-transcripts.gtf.db, containing all features
    - annotate with RefSeq IDs and TSL

  - annotate by priority:



- BED files annotation.
  Stores known features:
   - Feature: Exon, CDS, Transcript, Gene
   - Biotype: protein-coding, RNA
   - TSL
   - Transcript ID in Ensembl and RefSeq

  Annotation priority:
   - known protein_coding CDS|stop_codon|ncRNA_exon
   - known protein_coding UTR
   - known protein_coding transcript (annotate as "intron")
   - predicted protein_coding CDS|stop_codon|ncRNA_exon
   - predicted protein_coding UTR
   - predicted protein_coding transcript (annotate as "intron")

   -

  After annotation, write all features used for annotation and make a file for TargQC coverage reports.

- CDS BED file for Seq2C, if capture panel is not known.
  Ensembl-based, contains
   - Known canonical CDS
'''

def ensembl_gtf_fpath(genome):
    return _get_ensembl(join('gtf', 'ref-transcripts.gtf'), genome.split('-')[0])  # no -alt

def get_all_features(genome):
    return _get_ensembl('ensembl.bed', genome)

def biomart_fpath(genome):
    return _get_ensembl('biomart.tsv')

def get_canonical_cds(genome):
    """
    CDS. Used:
    - as TargQC target when no BED available
    - for Seq2C CNV calling when no capture BED available
    """
    bed = get_all_features(genome)
    return bed.filter(lambda r: r.fields[BedCols.FEATURE] in ['CDS', 'stop_codon'])

def get_merged_cds(genome):
    """
    Returns all CDS merged, used for TargQC general reports CDS coverage statistics for WGS.
    """
    bed = get_all_features(genome)
    return bed.filter(lambda r: r.fields[BedCols.FEATURE] in ['CDS', 'stop_codon']
                                and r.fields[BedTool.canonical]).merge()


def _get(relative_path, genome=None):
    """
    :param relative_path: relative path of the file inside the repository
    :param genome: genome name. Can contain chromosome name after comma, like hg19-chr20,
                   in case of BED, the returning BedTool will be with added filter.
    :return: BedTools object if it's a BED file, or filepath
    """
    chrom = None
    if genome:
        if '-chr' in genome:
            genome, chrom = genome.split('-')
        check_genome(genome)
        relative_path = relative_path.format(genome=genome)

    path = abspath(join(dirname(__file__), relative_path))
    if not isfile(path) and isfile(path + '.gz'):
        path += '.gz'

    if path.endswith('.bed') or path.endswith('.bed.gz'):
        if path.endswith('.bed.gz'):
            bedtools = which('bedtools')
            bedtools_v = bedtools_version(bedtools)
            if bedtools_v > 25:
                bed = BedTool(path)
            else:
                bed = BedTool(open_gzipsafe(path)).saveas()
        else:
            bed = BedTool(path)

        if chrom:
            bed = bed.filter(lambda r: r.chrom == chrom)
        return bed
    else:
        return path

def get_hgnc_gene_synonyms():
    return _get('HGNC_gene_synonyms.txt')



