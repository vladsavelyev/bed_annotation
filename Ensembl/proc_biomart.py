# with open('biomart_ensembl_to_refseq_refseqtx.tsv', 'w') as out, open('biomart_ensembl_to_refseq.tsv') as f:
# 	for i, l in enumerate(f):
# 		fs = l.strip('\n').split('\t')
# 		if i == 0:
# 			fs = ['Ensembl', 'RefSeq', fs[3]] + fs[6:] 
# 			out.write('\t'.join(fs) + '\n')	
# 		else:
# 			tx = fs[1] or fs[2]
# 			if not tx.startswith('N'): 
# 				tx = ''
# 			fs = [fs[0], tx, fs[3]] + fs[6:]
# 			out.write('\t'.join(fs) + '\n')

with open('biomart.tsv', 'w') as out, open('tsl_and_refseq_by_ensembl.tsv') as f:
	for i, l in enumerate(f):
		fs = l.strip('\n').split('\t')
		if i == 0:
			fs = ['Ensembl', 'RefSeq', 'TSL']
			# fs = ['Ensembl', 'RefSeq', fs[3]] + fs[6:] 
			out.write('\t'.join(fs) + '\n')	
		else:
			# tx = fs[1] or fs[2]
			# if not tx.startswith('N'): 
			# 	tx = ''
			# fs = [fs[0], tx, fs[3]] + fs[6:]
			out.write('\t'.join(fs) + '\n')


# /Users/vlad/vagrant/venv/bin/python /Users/vlad/vagrant/TargQC/GeneAnnotation/annotate_bed.py /Users/vlad/vagrant/TargQC/GeneAnnotation/Ensembl/hg38/test.bed -g hg38 --debug -o /Users/vlad/vagrant/TargQC/GeneAnnotation/Ensembl/hg38/test.anno.bed

