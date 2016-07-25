from os.path import isfile


# Merge RefSeq ncRNA and mRNA transcript fields
if isfile('biomart_ensembl_to_refseq.tsv'):
    with open('biomart_ensembl_to_refseq.tsv') as f, \
         open('biomart_ensembl_to_refseq_refseqtx.tsv', 'w') as out:
        for i, l in enumerate(f):
            fs = l.strip('\n').split('\t')
            if i == 0:
                fs = ['Ensembl', 'RefSeq', fs[3]] + fs[6:]
                out.write('\t'.join(fs) + '\n')
            else:
                tx = fs[1] or fs[2]
                if not tx.startswith('N'):
                    tx = ''
                fs = [fs[0], tx, fs[3]] + fs[6:]
                out.write('\t'.join(fs) + '\n')


# /Users/vlad/vagrant/venv/bin/python /Users/vlad/vagrant/TargQC/GeneAnnotation/annotate_bed.py
# /Users/vlad/vagrant/TargQC/GeneAnnotation/Ensembl/hg38/test.bed -g hg38 --debug -o /Users/vlad/vagrant/TargQC/GeneAnnotation/Ensembl/hg38/test.anno.bed
