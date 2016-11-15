### dinucleotide relative abundances

from stalling import get_FASTA_sequence

fasta = get_FASTA_sequence('/Volumes/USELESS/DATA/fasta/danio_rerio/GRCz10/cds/Danio_rerio.GRCz10.cds.all.fa.gz')


bases = ['T', 'C', 'A', 'G']
dinucleotides = [a+b for a in bases for b in bases]

dNTrelA = {}
diNTrelA = {}
for tr in fasta:
	dNTrelA[tr[0:18]] = {dNT: [0,0,0] for dNT in dinucleotides}
	diNTrelA[tr[0:18]] = {dNT+'_12': 0 for dNT in dinucleotides}
	diNTrelA[tr[0:18]] = {dNT+'_23': 0 for dNT in dinucleotides}
	diNTrelA[tr[0:18]] = {dNT+'_31': 0 for dNT in dinucleotides}
	l3 = float(len(fasta[tr])/3)
	l2 = float(len(fasta[tr])-1)
	# counts of bases at positions 1, 2, 3 in codon
	baseFreq = {base: [fasta[tr][0::3].count(base)/l3, fasta[tr][1::3].count(base)/l3, fasta[tr][2::3].count(base)/l3] for base in bases} ## local
	# dinucleotide frequencies
	dNTcount = {dNT: [0,0,0] for dNT in dinucleotides} ## local
	for i in range(1, len(fasta[tr])-1):
		ntX = fasta[tr][i]
		ntY = fasta[tr][i+1]
		dntXY = ntX+ntY
		if dntXY in dinucleotides:
			if i % 3 == 0: # 12
				dNTcount[dntXY][0] += 1
			elif i % 3 == 1: # 23
				dNTcount[dntXY][1] += 1
			elif i % 3 == 2: # 31
				dNTcount[dntXY][2] += 1
	# relative abundancies
	print dNTcount
	print baseFreq
	if all(c > 0 for c in sum(baseFreq.values(), [])):
		for dNT in dinucleotides:
			a12 = (dNTcount[dNT][0]/l2) / float( baseFreq[dNT[0]][0] * baseFreq[dNT[1]][1])
			a23 = (dNTcount[dNT][1]/l2) / float( baseFreq[dNT[0]][1] * baseFreq[dNT[1]][2])
			a31 = (dNTcount[dNT][2]/l2) / float( baseFreq[dNT[0]][2] * baseFreq[dNT[1]][0])
			dNTrelA[tr[0:18]][dNT] = [a12, a23, a31]
			diNTrelA[tr[0:18]][dNT+'_12'] = a12
			diNTrelA[tr[0:18]][dNT+'_23'] = a23
			diNTrelA[tr[0:18]][dNT+'_31'] = a31




### saving just in case
import cPickle as pickle
import gzip
pickle.dump(dNTrelA, gzip.open("/Volumes/USELESS/META/SHAPES/dNTrelA.p.gz", 'wb'))

import pandas
diNTrelAdf = pandas.DataFrame(diNTrelA).transpose()

diNTrelAdf.to_csv("/Volumes/USELESS/META/SHAPES/diNTrelAdf.csv")



l = []
for tr in dNTrelA:
	for dn in dNTrelA[tr]:
		l.extend(dNTrelA[tr][dn])


###################
### GC content
###################

GCcontent = {}
for tr in fasta:
	GCcontent[tr[0:18]] = (fasta[tr].count('G') + fasta[tr].count('C')) / float(len(fasta[tr]))


GCdf = pandas.DataFrame(GCcontent.items())
GCdf.to_csv("/Volumes/USELESS/META/SHAPES/GCdf.csv")


fasta_cdna = get_FASTA_sequence('/Volumes/USELESS/DATA/fasta/danio_rerio/GRCz10/cdna/Danio_rerio.GRCz10.cdna.all.fa.gz')

first5 = {bp: [0,0,0,0,0] for bp in bases}
######## first 5 bases
for tr in fasta_cdna:
	if len(fasta_cdna[tr]) > 5:
		for i in range(0,5):
			if fasta_cdna[tr][i] in bases:
				first5[fasta_cdna[tr][i]][i] += 1


total = [0,0,0,0,0]
for bp in first5:
	for i in range(0,5):
		total[i] += first5[bp][i]

freq = {bp: [0,0,0,0,0] for bp in bases}
for bp in first5:
	for i in range(0,5):
		freq[bp][i] = first5[bp][i] / float(total[i])



