# parse structures
import re

#p = re.compile('[\.\(\)]+\s\(')

#p = re.compile('\-\d+\.\d+')
p = re.compile('-?\d+(\.\d+)?')

fh = open("/Volumes/USELESS/DATA/fasta/danio_rerio/GRCz10/cdna/structures.txt", "r")
fout = open("/Volumes/USELESS/DATA/fasta/danio_rerio/GRCz10/cdna/min_free_energies.txt", "w")

f = 0
for line in fh:
	if f == 1:
		fout.write(tx+'\t'+minfe+'\n')
	f = 0
	if line.startswith(">"):
		tx = line[1:19]
		print tx
	elif line.startswith(".") or line.startswith("("):
		minfe = p.search(line).group()
		f = 1
		print minfe


# for the last line
fout.write(tx+'\t'+minfe+'\n')

fout.close()



###############
fh = open("/Volumes/USELESS/DATA/fasta/danio_rerio/GRCz10/cdna/structures.txt", "r")
fout = open("/Volumes/USELESS/DATA/fasta/danio_rerio/GRCz10/cdna/structures_only.txt", "w")

for line in fh:
	if line.startswith(">"):
		fout.write(line)
	elif line.startswith(".") or line.startswith("("):
		struc = line.split()[0]
		fout.write(struc+"\n")



fout.close()
