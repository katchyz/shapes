### AVOIDANCE

from RNAup_mRNAs_ncRNAs import *
from stalling import get_FASTA_sequence
from random import shuffle
from shoelaces_full_distribution import BED_get_CDS_mRNA_coord

from Bio.Seq import Seq
from Bio import SeqIO

import textwrap

# should it be CDS or CDNA?
# fasta = get_FASTA_sequence('/Volumes/USELESS/DATA/fasta/danio_rerio/GRCz10/cds/Danio_rerio.GRCz10.cds.all.fa.gz')

# RNAup_Executer(Single_RNA,Multiple_RNA,Interaction_Region)

# shuffle, preserving diNT sequence

# compare distribution of energies for 21nt fragments vs 200nt


########
s = 'somestring'
l = [s[i:i+2] for i in range(0, len(s), 2)]
shuffle(l) # now l is shuffled

##
r1 = "/Users/kasia/Desktop/r1.txt"
r2 = "/Users/kasia/Desktop/r2.txt"

RNAup_Executer(r1,r2,(1,20))

#########-----------------------------------
## exp
from RNAup_mRNAs_ncRNAs import *
f1 = "mrna2.fa"
f2 = "Danio_rerio.GRCz10.ncrna.fa"
#RNAup_Executer(f1,f2,(0,40))
RNAup_Executer_Parallel(f1,f2,(0,40),32)
## ctr
from RNAup_mRNAs_ncRNAs import *
f1 = "mrna2_shuffled.fa"
f2 = "Danio_rerio.GRCz10.ncrna.fa"
#RNAup_Executer(f1,f2,(0,40))
RNAup_Executer_Parallel(f1,f2,(0,40),32)
#########-----------------------------------


###########
# execute from command line

# fasta1_exp: get exact string (-20:20)
# fasta2_exp: get whole cdna sequence

# fasta1_ctr: shuffle fasta1_exp, preserving diNT
# fasta2_ctr: same as fasta2_exp

##
# test_ctr: shuffle 200nt from 5'

cdna = "/Volumes/USELESS/DATA/fasta/danio_rerio/GRCz10/cdna/Danio_rerio.GRCz10.cdna.all.fa.gz"
#fasta = SeqIO.parse(cdna,"fasta")
fasta = get_FASTA_sequence(cdna)


bedpath = '/Volumes/USELESS/DATA/genomes/BED/zebrafish_GRCz10.bed'
(intervals, geneToInterval, geneLength) = BED_to_intervals(bedpath)
orfs = BED_get_CDS_mRNA_coord(bedpath)

# geneToInterval['ENSDART00000124467.3'][0].strand
# orfs['ENSDART00000124467.3']

# get utr5 length (from  shoelaces_full_distribution, BED file)
utr5len = {}
for tx in orfs:
	if geneToInterval[tx][0].strand == '+':
		utr5len[tx] = orfs[tx][0][0]
	elif geneToInterval[tx][0].strand == '-':
		utr5len[tx] = geneLength[tx] - orfs[tx][0][1]


fasta1_exp = open('/Volumes/USELESS/META/SHAPES/avoidance/fasta1_exp.fa', 'w')
fasta2_exp = open('/Volumes/USELESS/META/SHAPES/avoidance/fasta2_exp.fa', 'w')
# extract fragments, write to file fasta1_exp
for tx in fasta:
	if tx in utr5len:
		if utr5len[tx] > 30:
			seq = fasta[tx][utr5len[tx]-20:utr5len[tx]+20]
			##
			#fasta1_exp.write("%s%s\n%s\n" % ('>', tx, seq))
			#fasta1_exp.write(textwrap.fill("%s%s\n%s\n" % ('>', tx, seq)), width=60)
			#fasta2_exp.write("%s%s\n%s\n" % ('>', tx, fasta[tx]))
			#fasta2_exp.write(textwrap.fill("%s%s\n%s\n" % ('>', tx, fasta[tx])), width=60)
			##
			fasta1_exp.write("%s%s\n" % ('>', tx))
			fasta1_exp.write(textwrap.fill(seq, width=60))
			fasta1_exp.write('\n')
			fasta2_exp.write("%s%s\n" % ('>', tx))
			fasta2_exp.write(textwrap.fill(fasta[tx], width=60))
			fasta2_exp.write('\n')

fasta1_exp.close()
fasta2_exp.close()

##
utr5 = {}
for name in utr5len:
	utr5[name[0:18]] = utr5len[name]

fas = {}
for name in fasta:
	fas[name[0:18]] = fasta[name]

names_mrna = open('/Volumes/USELESS/META/SHAPES/avoidance/het.txt', 'r')
mrna = open('/Volumes/USELESS/META/SHAPES/avoidance/mrna.fa', 'w')

for line in names_mrna:
	tx = line.split()[0]
	if tx in utr5:
		if utr5[tx] > 30:
			if tx in fas:
				seq = fas[tx][utr5[tx]-20:utr5[tx]+20]
				mrna.write("%s%s\n%s\n" % ('>', tx, seq))

mrna.close()

##
names_mrna = open('/Volumes/USELESS/META/SHAPES/avoidance/het.txt', 'r')
mrna_shuffled = open('/Volumes/USELESS/META/SHAPES/avoidance/mrna_shuffled.fa', 'w')

for line in names_mrna:
	tx = line.split()[0]
	if tx in utr5:
		if utr5[tx] > 30:
			if tx in fas:
				seq = fas[tx][utr5[tx]-20:utr5[tx]+20]
				l = [seq[i:i+2] for i in range(0, len(seq), 2)]
				shuffle(l)
				j = "".join(l)
				mrna_shuffled.write("%s%s\n%s\n" % ('>', tx, j))

mrna_shuffled.close()

# shuffle fragments, write to file fasta1_ctr
fasta1_ctr = open('/Volumes/USELESS/META/SHAPES/avoidance/fasta1_ctr.fa', 'w')

for tx in fasta:
	if tx in utr5len:
		if utr5len[tx] > 30:
			seq = fasta[tx][utr5len[tx]-20:utr5len[tx]+20]
			##
			l = [seq[i:i+2] for i in range(0, len(seq), 2)]
			shuffle(l)
			j = "".join(l)
			##
			fasta1_ctr.write("%s%s\n%s\n" % ('>', tx, j))


fasta1_ctr.close()


#############
## 5' ends of transcripts vs tRNAs

mrna = open('/Volumes/USELESS/META/SHAPES/avoidance/mrna_5end.fa', 'w')

for line in names_mrna:
	tx = line.split()[0]
	if tx in fas:
		seq = fas[tx][0:20]
		mrna.write("%s%s\n%s\n" % ('>', tx, seq))

mrna.close()

mrna_shuffled = open('/Volumes/USELESS/META/SHAPES/avoidance/mrna_shuffled_5end.fa', 'w')

for line in names_mrna:
	tx = line.split()[0]
	if tx in fas:
		seq = fas[tx][0:20]
		l = [seq[i:i+2] for i in range(0, len(seq), 2)]
		shuffle(l)
		j = "".join(l)
		mrna_shuffled.write("%s%s\n%s\n" % ('>', tx, j))

mrna_shuffled.close()

#########-----------------------------------
## exp
from RNAup_mRNAs_ncRNAs import *
f1 = "mrna_5end.fa" # get cdna
f2 = "cdhit.fa"
#RNAup_Executer(f1,f2,(0,40))
RNAup_Executer_Parallel(f1,f2,(0,20),32)
## ctr
from RNAup_mRNAs_ncRNAs import *
f1 = "mrna_shuffled_5end.fa"
f2 = "cdhit.fa"
#RNAup_Executer(f1,f2,(0,40))
RNAup_Executer_Parallel(f1,f2,(0,20),32)
#########-----------------------------------

#########-----------------------------------
## exp
from RNAup_mRNAs_ncRNAs import *
f1 = "mrna_5end.fa" # get cdna
f2 = "Danio_rerio.GRCz10.ncrna.fa"
#RNAup_Executer(f1,f2,(0,40))
RNAup_Executer_Parallel(f1,f2,(0,20),32)
## ctr
from RNAup_mRNAs_ncRNAs import *
f1 = "mrna_shuffled_5end.fa"
f2 = "Danio_rerio.GRCz10.ncrna.fa"
#RNAup_Executer(f1,f2,(0,40))
RNAup_Executer_Parallel(f1,f2,(0,20),32)
#########-----------------------------------


