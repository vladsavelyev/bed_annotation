import gffutils
import tempfile
import os
import random
import gzip
from ngs_utils.file_utils import file_exists, file_transaction
from ngs_utils.logger import debug


def guess_infer_extent(gtf_file):
    """
    guess if we need to use the gene extent option when making a gffutils
    database by making a tiny database of 1000 lines from the original
    GTF and looking for all of the features
    """
    _, ext = os.path.splitext(gtf_file)
    tmp_out = tempfile.NamedTemporaryFile(suffix=".gtf", delete=False).name
    with open(tmp_out, "w") as out_handle:
        count = 0
        in_handle = open(gtf_file) if ext != ".gz" else gzip.open(gtf_file)
        for line in in_handle:
            if count > 1000:
                break
            out_handle.write(line)
            count += 1
        in_handle.close()
    db = gffutils.create_db(tmp_out, dbfn=":memory:", infer_gene_extent=False)
    os.remove(tmp_out)
    features = [x for x in db.featuretypes()]
    if "gene" in features and "transcript" in features:
        return False
    else:
        return True

def get_gtf_db(gtf, in_memory=False):
    """
    create a gffutils DB
    """
    db_file = gtf + '.db'
    if gtf.endswith('.gz'):
        db_file = gtf[:-3] + '.db'
    if file_exists(db_file):
        return gffutils.FeatureDB(db_file)
    db_file = ':memory:' if in_memory else db_file
    if in_memory or not file_exists(db_file):
        debug('GTF database does not exist, creating...')
        infer_extent = guess_infer_extent(gtf)
        db = gffutils.create_db(gtf, dbfn=db_file,
                                infer_gene_extent=infer_extent)
        return db
    else:
        return gffutils.FeatureDB(db_file)

def gtf_to_bed(gtf, alt_out_dir=None):
    """
    create a BED file of transcript-level features with attached gene name
    or gene ids
    """
    out_file = os.path.splitext(gtf)[0] + '.bed'
    if file_exists(out_file):
        return out_file
    if not os.access(os.path.dirname(out_file), os.W_OK | os.X_OK):
        if not alt_out_dir:
            raise IOError('Cannot write transcript BED output file %s' % out_file)
        else:
            out_file = os.path.join(alt_out_dir, os.path.basename(out_file))
    with open(out_file, "w") as out_handle:
        db = get_gtf_db(gtf)
        for feature in db.features_of_type('transcript', order_by=("seqid", "start", "end")):
            chrom = feature.chrom
            start = feature.start
            end = feature.end
            attributes = feature.attributes.keys()
            strand = feature.strand
            name = (feature['gene_name'][0] if 'gene_name' in attributes else
                    feature['gene_id'][0])
            line = "\t".join([str(x) for x in [chrom, start, end, name, ".",
                                               strand]])
            out_handle.write(line + "\n")
    return out_file

def complete_features(db):
    """
    iterator returning features which are complete (have a 'gene_id' and a
    'transcript_id')
    """
    for feature in db.all_features():
        gene_id = feature.attributes.get('gene_id', [None])[0]
        transcript_id = feature.attributes.get('transcript_id', [None])[0]
        if gene_id and transcript_id and feature.featuretype != "transcript":
            yield feature

def partition_gtf(gtf, coding=False, out_file=False):
    """
    return a GTF file of all non-coding or coding transcripts. the GTF must be annotated
    with gene_biotype = "protein_coding" or to have the source column set to the
    biotype for all coding transcripts. set coding to
    True to get only the coding, false to get only the non-coding
    """
    if out_file and file_exists(out_file):
        return out_file
    if not out_file:
        out_file = tempfile.NamedTemporaryFile(delete=False,
                                               suffix=".gtf").name

    if coding:
        pred = lambda biotype: biotype and biotype == "protein_coding"
    else:
        pred = lambda biotype: biotype and biotype != "protein_coding"

    biotype_lookup = _biotype_lookup_fn(gtf)

    db = get_gtf_db(gtf)
    with file_transaction(out_file) as tx_out_file:
        with open(tx_out_file, "w") as out_handle:
            for feature in db.all_features():
                biotype = biotype_lookup(feature)
                if pred(biotype):
                    out_handle.write(str(feature) + "\n")
    return out_file

def get_coding_noncoding_transcript_ids(gtf):
    """
    return a set of coding and non-coding transcript_ids from a GTF
    """
    coding_gtf = partition_gtf(gtf, coding=True)
    coding_db = get_gtf_db(coding_gtf)
    coding_ids = set([x['transcript_id'][0] for x in coding_db.all_features()
                  if 'transcript_id' in x.attributes])
    noncoding_gtf = partition_gtf(gtf)
    noncoding_db = get_gtf_db(noncoding_gtf)
    noncoding_ids = set([x['transcript_id'][0] for x in noncoding_db.all_features()
                     if 'transcript_id' in x.attributes])
    return coding_ids, noncoding_ids

def get_gene_source_set(gtf):
    """
    get a dictionary of the set of all sources for a gene
    """
    gene_to_source = {}
    db = get_gtf_db(gtf)
    for feature in complete_features(db):
        gene_id = feature['gene_id'][0]
        sources = gene_to_source.get(gene_id, set([])).union(set([feature.source]))
        gene_to_source[gene_id] = sources
    return gene_to_source

def get_transcript_source_set(gtf):
    """
    get a dictionary of the set of all sources of the gene for a given
    transcript
    """
    gene_to_source = get_gene_source_set(gtf)
    transcript_to_source = {}
    db = get_gtf_db(gtf)
    for feature in complete_features(db):
        gene_id = feature['gene_id'][0]
        transcript_to_source[feature['transcript_id'][0]] = gene_to_source[gene_id]
    return transcript_to_source

def _biotype_lookup_fn(gtf):
    """
    return a function that will look up the biotype of a feature
    this checks for either gene_biotype or biotype being set or for the source
    column to have biotype information
    """
    db = get_gtf_db(gtf)
    sources = set([feature.source for feature in db.all_features()])
    gene_biotypes = set([feature.attributes.get("gene_biotype", [None])[0]
                         for feature in db.all_features()])
    biotypes = set([feature.attributes.get("biotype", [None])[0]
                    for feature in db.all_features()])
    if "protein_coding" in sources:
        return lambda feature: feature.source
    elif "protein_coding" in biotypes:
        return lambda feature: feature.attributes.get("biotype", [None])[0]
    elif "protein_coding" in gene_biotypes:
        return lambda feature: feature.attributes.get("gene_biotype", [None])[0]
    else:
        return None

def tx2genefile(gtf, out_file=None):
    """
    write out a file of transcript->gene mappings.
    use the installed tx2gene.csv if it exists, else write a new one out
    """
    installed_tx2gene = os.path.join(os.path.dirname(gtf), "tx2gene.csv")
    if file_exists(installed_tx2gene):
        return installed_tx2gene
    if file_exists(out_file):
        return out_file
    with file_transaction(out_file) as tx_out_file:
        with open(tx_out_file, "w") as out_handle:
            for k, v in transcript_to_gene(gtf).items():
                out_handle.write(",".join([k, v]) + "\n")
    return out_file

def transcript_to_gene(gtf):
    """
    return a dictionary keyed by transcript_id of the associated gene_id
    """
    gene_lookup = {}
    for feature in complete_features(get_gtf_db(gtf)):
        gene_id = feature.attributes.get('gene_id', [None])[0]
        transcript_id = feature.attributes.get('transcript_id', [None])[0]
        gene_lookup[transcript_id] = gene_id
    return gene_lookup
