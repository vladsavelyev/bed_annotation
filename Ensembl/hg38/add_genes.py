import os
import gffutils
import sys
from bcbio import utils
import pybedtools
import subprocess
in_file = sys.argv[1]

def gtf_to_bed(gtf):
    """
    create a BED file of transcript-level features with attached gene name
    or gene ids
    """
    out_file = gtf.replace('.gtf.db', '.bed')
    if os.path.isfile(out_file):
        return out_file
    if not os.access(os.path.dirname(out_file), os.W_OK | os.X_OK):
        if not alt_out_dir:
            raise IOError("Cannot write transcript BED output file %s" % out_file)
        else:
            out_file = os.path.join(alt_out_dir, os.path.basename(out_file))
    with open(out_file, "w") as out_handle:
        db = gffutils.FeatureDB(gtf)
        for feature in db.features_of_type('transcript', order_by=("seqid", "start", "end")):
            chrom = feature.chrom
            start = feature.start - 1
            end = feature.end
            attributes = feature.attributes.keys()
            strand = feature.strand
            name = (feature['gene_name'][0] if 'gene_name' in attributes else
                    feature['gene_id'][0])
            line = "\t".join([str(x) for x in [chrom, start, end, name, ".",
                                               strand]])
            out_handle.write(line + "\n")
    return out_file

def add_genes(in_file, gene_file, max_distance=0):
    if gene_file and os.path.isfile(in_file):
        out_file = "%s.add_genes.bed" % utils.splitext_plus(in_file)[0]
        if not utils.file_uptodate(out_file, in_file):
            input_rec = iter(pybedtools.BedTool(in_file)).next()
            # keep everything after standard chrom/start/end, 1-based
            extra_fields = range(4, len(input_rec.fields) + 1)
            # keep the new gene annotation
            gene_index = len(input_rec.fields) + 4
            extra_fields.append(gene_index)
            columns = ",".join([str(x) for x in extra_fields])
            max_column = max(extra_fields) + 1
            ops = ",".join(["distinct"] * len(extra_fields))
            fai_file = '/Users/vlad/vagrant/TargQC/Utils/reference_data/fai/hg38.fa.fai'
            # swap over gene name to '.' if beyond maximum distance
            # cut removes the last distance column which can cause issues
            # with bedtools merge: 'ERROR: illegal character '.' found in integer conversion of string'
            distance_filter = (r"""awk -F$'\t' -v OFS='\t' '{if ($NF > %s) $%s = "."} {print}'""" %
                               (max_distance, gene_index))
            cmd = ("set -o pipefail; bedtools closest -g <(cut -f1,2 {fai_file}) "
                   "-d -t all -a {in_file} -b {gene_file} | "
                   "{distance_filter} | cut -f 1-{max_column} | "
                   "bedtools merge -i - -c {columns} -o {ops} -delim ',' > {out_file}")
            res = subprocess.call(cmd.format(**locals()), shell=True)

gene_file = gtf_to_bed('/Users/vlad/vagrant/TargQC/GeneAnnotation/Ensembl/hg38/gtf/ref-transcripts.gtf.db')
add_genes(in_file, gene_file)


# with open('test.anno.dups.bed', 'w') as out:
# 	line_counts = dict()
# 	with open('test.anno.bed') as f:
# 		for i, l in enumerate(f):
# 			fs = l.strip('\n').split('\t')
# 			reg = tuple(fs[:3])
# 			if reg not in line_counts:
# 				line_counts[reg] = 0
# 			line_counts[reg] += 1

# 	with open('test.anno.bed') as f:
# 		for i, l in enumerate(f):
# 			fs = l.strip('\n').split('\t')
# 			reg = tuple(fs[:3])
# 			if line_counts[reg] > 1:
# 				out.write('\t'.join(fs) + '\n')




'''
add_genes(distance=100)
compare add_genes with mine
try selecting most trustful
'''

