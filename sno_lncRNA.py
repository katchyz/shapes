gtf = open('/Volumes/USELESS/DATA/genomes/GTF/danio_rerio/Danio_rerio.GRCz10.84.gtf', 'r')


r = []
sno = []
for line in gtf:
	if not line.startswith('#'):
		fields = line.split('\t')
		ginfo = fields[8].split(';')
		for gi in ginfo:
			if gi.startswith(' transcript_id'):
				tid = gi[16:34]
			if gi.startswith(' gene_biotype'):
				if gi == ' gene_biotype "rRNA"':
					r.append(tid)
				elif gi == ' gene_biotype "snoRNA"':
					sno.append(tid)


sno = set(sno)
sno = list(sno)
r = set(r)
r = list(r)


snoDf = pandas.DataFrame(sno)
rDf = pandas.DataFrame(r)

rDf.to_csv("/Volumes/USELESS/META/SHAPES/rDf.csv")
snoDf.to_csv("/Volumes/USELESS/META/SHAPES/snoDf.csv")


lnc = []
for line in gtf:
	if not line.startswith('#'):
		fields = line.split('\t')
		ginfo = fields[8].split(';')
		for gi in ginfo:
			if gi.startswith(' transcript_id'):
				tid = gi[16:34]
			if gi.startswith(' gene_biotype'):
				if gi == ' gene_biotype "lincRNA"':
					lnc.append(tid)


lnc = set(lnc)
lnc = list(lnc)

lncDf = pandas.DataFrame(lnc)
lncDf.to_csv("/Volumes/USELESS/META/SHAPES/lncDf.csv")



