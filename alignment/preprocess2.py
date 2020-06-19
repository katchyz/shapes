### preprocess.py read1 read2
### trim barcodes, save to file
### remove empty reads
### save trimmed fastq

import sys
from itertools import islice
from itertools import izip_longest


def grouper(iterable, n, fillvalue=None):
	args = [iter(iterable)] * n
	return izip_longest(*args, fillvalue=fillvalue)


N = 4 # process four lines at the time
trim1_3end = 15
trim2_5end = 7
trim2_3end = 50

o1 = open("read1.fastq", "w")
b = open("barcodes.txt", "w")
o2 = open("read2.fastq", "w")

with open(sys.argv[1]) as f1, open(sys.argv[2]) as f2:
	for lines in grouper(zip(f1,f2), N, ''):
		# read1
		read_name = lines[0][0].split()[0][1:]
		barcode = lines[1][0].split()[0][0:7]
		trimmed_seq1 = lines[1][0].split()[0][7:-trim1_3end]
		trimmed_qual1 = lines[3][0].split()[0][7:-trim1_3end]
		# read2
		trimmed_seq2 = lines[1][1].split()[0][trim2_5end:-trim2_3end]
		trimmed_qual2 = lines[3][1].split()[0][trim2_5end:-trim2_3end]
		# save to files if not empty
		if trimmed_seq1 and trimmed_seq2:
			# barcodes
			b.write(read_name + "\t" + barcode + "\n")
			# read1
			o1.write(lines[0][0] + trimmed_seq1 + "\n" + lines[2][0] + trimmed_qual1 + "\n")
			# read2
			o2.write(lines[0][1] + trimmed_seq2 + "\n" + lines[2][1] + trimmed_qual2 + "\n")



b.close()
o1.close()
o2.close()





