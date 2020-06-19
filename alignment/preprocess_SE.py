### preprocess.py read1 (SINGLE END)
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
trim1_3end = 1

o1 = open("read1.fastq", "w")
b = open("barcodes.txt", "w")

with open(sys.argv[1]) as f1:
	for lines in grouper(f1, N, ''):
		# read1
		read_name = lines[0].split()[0][1:]
		barcode = lines[1].split()[0][0:7]
		trimmed_seq1 = lines[1].split()[0][7:-trim1_3end]
		trimmed_qual1 = lines[3].split()[0][7:-trim1_3end]
		# save to files if not empty
		if trimmed_seq1:
			# barcodes
			b.write(read_name + "\t" + barcode + "\n")
			# read1
			o1.write(lines[0][0] + trimmed_seq1 + "\n" + lines[2][0] + trimmed_qual1 + "\n")



b.close()
o1.close()


