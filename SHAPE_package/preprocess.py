### preprocess.py read1 read2
### trim barcodes, save to file
### save trimmed fastq

import sys

f1 = open(sys.argv[1]) #read1
o1 = open("read1.fastq", "w")
b = open("barcodes.txt", "w")

for line in f1:
        if line.startswith("@"):
                read_name = line.split()[0][1:]
                nextline = next(f1)
                barcode = nextline.split()[0][0:7]
                trimmed_seq = nextline.split()[0][7:-15]
                # save
                o1.write(line)
                o1.write(trimmed_seq + "\n")
                b.write(read_name + "\t" + barcode + "\n")
        elif line.startswith("+"):
                qualline = next(f1)
                trimmed_qual = qualline.split()[0][7:-15]
                # save
                o1.write(line)
                o1.write(trimmed_qual + "\n")


f1.close
o1.close()
b.close()



f2 = open(sys.argv[2]) #read2
o2 = open("read2.fastq", "w")

for line in f2:
        if line.startswith("@"):
                nextline = next(f2)
                trimmed_seq = nextline.split()[0][15:]
                # save
                o2.write(line)
                o2.write(trimmed_seq + "\n")
        elif line.startswith("+"):
                qualline = next(f2)
                trimmed_qual = qualline.split()[0][15:]
                # save
                o2.write(line)
                o2.write(trimmed_qual + "\n")


f2.close
o2.close()


