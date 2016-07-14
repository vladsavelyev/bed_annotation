with open('ensembl.bed') as f:
	for l in f:
		print str(len(l.split('\t'))), l.strip()
