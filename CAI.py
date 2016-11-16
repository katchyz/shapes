### CAI

from stalling import get_FASTA_sequence
import pandas

fasta = get_FASTA_sequence('/Volumes/USELESS/DATA/fasta/danio_rerio/GRCz10/cds/Danio_rerio.GRCz10.cds.all.fa.gz')


# codon usage table
### cut[AA][codon] = [fraction, frequency:per thousand, number]
### cut[AA][codon] = frequency(per thousand)

def get_codon_usage(filepath):
	""" Reads codon usage table """
	f = open(filepath, 'r')
	aas = 'ACDEFGHIKLMNPQRSTVWY*'
	cut = {aa: {} for aa in aas}
	for line in f:
		line = line.strip()
		fields = line.split()
		for i in range(0,len(fields),5):
			#cut[fields[i+1]][fields[i]] = [float(fields[i+2]), float(fields[i+3]), eval(fields[i+4])]
			cut[fields[i+1]][fields[i].replace('U', 'T')] = float(fields[i+3])
	return cut


cut = get_codon_usage('/Volumes/USELESS/DATA/fasta/danio_rerio/GRCz10/codon_usage_table.txt')


### calculate RSCU and W (from table)
# relative synonymous codon usage
rscu = {}
for aa in cut:
	rscu[aa] = {}
	for codon in cut[aa]:
		rscu[aa][codon] = cut[aa][codon] / ((1/float(len(cut[aa]))) * sum(cut[aa].values()))


bases = ['T', 'C', 'A', 'G']
codons = [a+b+c for a in bases for b in bases for c in bases]
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
codon_table = dict(zip(codons, amino_acids))


# w = {}
# for aa in rscu:
# 	for codon in rscu[aa]:
# 		w[codon] = rscu[aa][codon] / max(rscu[aa].values())


### calculate CAI for each gene in fasta
cai = {}
for tr in fasta:
	cai_obs = 1
	cai_max = 1
	for i in range(0,len(fasta[tr]),3):
		codon = fasta[tr][i:i+3]
		if len(codon) == 3:
			if codon in codon_table:
				cai_obs *= rscu[codon_table[codon]][codon]
				cai_max *= max(rscu[codon_table[codon]].values())
	cai[tr[0:18]] = (cai_obs / cai_max) ** (1.0/(len(fasta[tr])/3))


caiDf = pandas.DataFrame(cai.items)
caiDf.to_csv("/Volumes/USELESS/META/SHAPES/caiDf.csv")


#################
### CAI on the last 50 codons

### calculate CAI for each gene in fasta
cai_last50 = {}
for tr in fasta:
	cai_obs = 1
	cai_max = 1
	for i in range(len(fasta[tr])-150,len(fasta[tr]),3):
		codon = fasta[tr][i:i+3]
		if len(codon) == 3:
			if codon in codon_table:
				cai_obs *= rscu[codon_table[codon]][codon]
				cai_max *= max(rscu[codon_table[codon]].values())
	cai_last50[tr[0:18]] = (cai_obs / cai_max) ** (1.0/50)


caiDf_last50 = pandas.DataFrame(cai_last50.items())
caiDf_last50.to_csv("/Volumes/USELESS/META/SHAPES/caiDf_last50.csv")












