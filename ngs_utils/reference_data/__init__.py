from os.path import dirname, join, abspath, splitext, isfile

from ngs_utils.file_utils import verify_file, adjust_path, verify_dir
from ngs_utils.key_genes_utils import get_genes_from_file
from ngs_utils.logger import critical, debug


SUPPORTED_GENOMES = ['hg19', 'hg19-noalt', 'hg38', 'hg38-noalt', 'hg19-chr21', 'GRCh37', 'mm10']


def check_genome(genome):
    if genome not in SUPPORTED_GENOMES:
        critical('Genome ' + str(genome) + ' is not supported. Supported genomes: ' + ', '.join(SUPPORTED_GENOMES))

def _get(relative_path, genome=None, is_critical=False):
    if genome:
        check_genome(genome)
    else:
        genome = ''
    relative_path = relative_path.format(genome=genome)

    path = abspath(join(dirname(__file__), relative_path))
    if is_critical:
        return verify_file(path, is_critical=True)
    return path


######################
######## FAI #########
def get_fai(genome, is_critical=False):
    return _get(join('fai', '{genome}.fa.fai'), genome, is_critical=is_critical)

def get_chrom_lengths(genome=None, fai_fpath=None):
    assert genome or fai_fpath, f'One of genome or fai_fpath should be not None: genome={genome}, fai_fpath={fai_fpath}'

    if not fai_fpath:
        check_genome(genome)
        fai_fpath = get_fai(genome)
    else:
        fai_fpath = verify_file(fai_fpath, is_critical=True)
        if not fai_fpath.endswith('.fai') and not fai_fpath.endswith('.fa'):
            critical('Error: .fai or .fa is accepted.')

    chr_lengths = []

    if fai_fpath.endswith('.fa'):
        debug('Reading genome sequence (.fa) to get chromosome lengths')
        with open(fai_fpath, 'r') as handle:
            from Bio import SeqIO
            reference_records = SeqIO.parse(handle, 'fasta')
            for record in reference_records:
                chrom = record.id
                chr_lengths.append((chrom, len(record.seq)))

    else:
        debug('Reading genome index file (.fai) to get chromosome lengths')
        with open(fai_fpath, 'r') as handle:
            for line in handle:
                line = line.strip()
                if line:
                    chrom, length = line.split()[0], int(line.split()[1])
                    chr_lengths.append((chrom, length))

    return chr_lengths

def get_chrom_order(genome=None, fai_fpath=None):
    chr_lengths = get_chrom_lengths(genome, fai_fpath)
    chr_order = {c: i for i, (c, l) in enumerate(chr_lengths)}
    return chr_order

def ucsc_to_ensembl(genome, is_critical=False):
    """ mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -N -e "select * from ucscToEnsembl;" hg19 > hg19.ucscToEnsembl.tsv
        mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -N -e "select * from ucscToEnsembl;" hg38 > hg38.ucscToEnsembl.tsv
    """
    return _get(join('fai', '{genome}.ucscToEnsembl.tsv'), genome, is_critical=is_critical)


#############################
######## Gene lists #########

def get_signatures_probabilities(is_critical=False):
    return _get('signatures_probabilities.txt', is_critical=is_critical)

def get_key_genes():
    return _get('key_genes/umccr_cancer_genes.latest.tsv')

def get_key_genes_txt():
    return _get('key_genes/umccr_cancer_genes.latest.genes')

def get_key_tsgenes_txt():
    return _get('key_genes/umccr_cancer_genes.tsgenes.latest.genes')

def get_key_genes_set(fpath=get_key_genes()):
    genes = set()
    with open(fpath) as f:
        for i, l in enumerate(f):
            if i != 0:
                genes.add(l.strip().split('\t')[0])
    return genes

def get_key_genes_bed(genome, is_critical=False, coding_only=False):
    return _get(f'key_genes/umccr_cancer_genes.{genome}.{"transcript" if not coding_only else "coding"}.bed',
                is_critical=is_critical)

def get_predispose_genes_txt():
    return _get('key_genes/sources/predispose_genes.txt')

def get_predispose_genes_bed(genome, is_critical=False, coding_only=False):
    return _get(f'key_genes/predispose_genes.{genome}.{"transcript" if not coding_only else "coding"}.bed',
                is_critical=is_critical)


#############################
###### Known fusions ########

"""
Known fusions are curated by [Hartwig Medical Foundation](https://nc.hartwigmedicalfoundation.nl/index.php/s/a8lgLsUrZI5gndd?path=%2FHMF-Pipeline-Resources)

[Identification of gene fusions](https://www.biorxiv.org/content/biorxiv/early/2018/09/20/415133.full.pdf):

    For each structural variant, every combination of annotated overlapping transcripts from each breakend
    was tested to see if it could potentially form an intronic inframe fusion. A list of 411 curated known fusion
    pairs was sourced by taking the union of known fusions from the following external databases:
        - Cosmic curated fusions (v83)
        - OncoKb (download = 01-mar-2018)
        - CGI (update: 17-jan-2018)
        - CIViC (download = 01-mar-2018)

    We then also created a list of promiscuous fusion partners, defined as any gene which appears on the same side in more 
    than 3 of the curated fusion pairs OR is marked as promiscuous in either OncoKb, CGI or CIVIC.
    
    For each promiscuous partner we also curated a list of essential domains that must be preserved to form
    a viable fusion partner.
    
    Finally, we report an intronic inframe fusion if it matches an exact fusion from the curated list 
    OR is intergenic and matches 5’ promiscuous OR matches 3’ promiscuous gene.

Finally, we report an intronic inframe fusion if the following conditions are met:
    - Matches an exact fusion from the curated list 
      OR is intergenic and matches 5’ promiscuous 
      OR matches 3’ promiscuous gene 
    - Curated domains are preserved
    - Does not involve the 3’UTR region of either gene
    - For intragenic fusions, must start and end in coding regions of the gene.

To download HMF fusions:

```
wget https://nc.hartwigmedicalfoundation.nl/index.php/s/a8lgLsUrZI5gndd/download?path=%2FHMF-Pipeline-Resources&files=KnownFusions.zip -O hmf.zip
unzip hmf.zip
ls
$ knownFusionPairs.csv
$ knownPromiscuousFive.csv
$ knownPromiscuousThree.csv
```
    
There is also a more broad list of fusions from [FusionCatcher](https://github.com/ndaniel/fusioncatcher):

```
wget https://raw.githubusercontent.com/ndaniel/fusioncatcher/master/bin/generate_known.py
grep "        \['" generate_known.py | sed "s#        \['##" | sed "s#','#,#" | sed "s#'\],##" | sed "s#'\]##" | sort -u > fusioncatcher_pairs.txt
```

There are also known fusion genes in [NGC](http://ncg.kcl.ac.uk/), one of the sources for UMCCR cancer gene list.

We compare the lists with `compare.R`. 

The FusionsCatcher list is too big so we assign a lower tier to a matching fusion when prioritizing.
"""

def get_known_fusion_pairs():
    return _get('fusions/knownFusionPairs.csv')

def get_known_fusion_heads():
    return _get('fusions/knownPromiscuousFive.csv')

def get_known_fusion_tails():
    return _get('fusions/knownPromiscuousThree.csv')

def get_fusioncatcher_pairs():
    return _get('fusions/fusioncatcher_pairs.txt')



