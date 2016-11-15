### cutting adapters
# 2-4 cell
cutadapt -q 20 -a AGATCGGAAGAGCACACGTCT --minimum-length 20 --paired-output tmp.2.fastq -o tmp.1.fastq R1_ctrl.fastq.gz R2_ctrl.fastq.gz
cutadapt -q 20 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 20 --paired-output R1_ctrl_cut.fastq.gz -o R2_ctrl_cut.fastq.gz tmp.2.fastq tmp.1.fastq
rm tmp.1.fastq tmp.2.fastq


cutadapt -q 20 -a AGATCGGAAGAGCACACGTCT --minimum-length 20 --paired-output tmp.2.fastq -o tmp.1.fastq R1_nai.fastq.gz R2_nai.fastq.gz
cutadapt -q 20 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 20 --paired-output R1_nai_cut.fastq.gz -o R2_nai_cut.fastq.gz tmp.2.fastq tmp.1.fastq
rm tmp.1.fastq tmp.2.fastq

# 256 cell
cutadapt -q 20 -a AGATCGGAAGAGCACACGTCT --minimum-length 20 --paired-output tmp.2.fastq -o tmp.1.fastq 256cell_ctrl_R1.fastq.gz 256cell_ctrl_R2.fastq.gz
cutadapt -q 20 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 20 --paired-output 256cell_ctrl_R1_cut.fastq.gz -o 256cell_ctrl_R2_cut.fastq.gz tmp.2.fastq tmp.1.fastq
rm tmp.1.fastq tmp.2.fastq

cutadapt -q 20 -a AGATCGGAAGAGCACACGTCT --minimum-length 20 --paired-output tmp.2.fastq -o tmp.1.fastq 256cell_NAI_R1.fastq.gz 256cell_NAI_R2.fastq.gz
cutadapt -q 20 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 20 --paired-output 256cell_NAI_R1_cut.fastq.gz -o 256cell_NAI_R2_cut.fastq.gz tmp.2.fastq tmp.1.fastq
rm tmp.1.fastq tmp.2.fastq


### preprocessing

cd /Volumes/USELESS/DATA/Shape-Seq

~/Documents/PhD/scripts/SHAPES/RNAprobBash-master/preprocessing.sh -1 256cell_NAI_R1_cut.fastq -2 256cell_NAI_R2_cut.fastq -b NNNNNNN -t 15 -o "preproc_256cell_NAI"



### mapping
# 2-4cell
r1="/Home/ii/katchyz/OUT/SHAPES_out/preproc_2-4cell_NAI/read1.fastq"
r2="/Home/ii/katchyz/OUT/SHAPES_out/preproc_2-4cell_NAI/read2.fastq"

bowtie2 -p 8 -X 700 -x /Home/ii/katchyz/DATA/genomes/index_GRCz10 -1 $r1 -2 $r2 | gzip > aln_2-4cell_NAI.sam.gz

# 256cell
r1="/Home/ii/katchyz/OUT/SHAPES_out/preproc_256cell_NAI/read1.fastq"
r2="/Home/ii/katchyz/OUT/SHAPES_out/preproc_256cell_NAI/read2.fastq"

bowtie2 -p 8 -X 700 -x /Home/ii/katchyz/DATA/genomes/index_GRCz10 -1 $r1 -2 $r2 | gzip > aln_256cell_NAI.sam.gz


### summarize barcodes >>>>>>>>>>>> run locally, gawk problem (download sam.gz files)
#barcodes="/Home/ii/katchyz/OUT/SHAPES_out/preproc_2-4cell_NAI/barcodes.txt"
#~/Documents/PhD/scripts/SHAPES/RNAprobBash-master/summarize_unique_barcodes2.sh -f aln_2-4cell_NAI.sam.gz -b /Volumes/USELESS/DATA/Shape-Seq/preproc_2-4cell_NAI/barcodes.txt -t -k -o "summarize_2-4cell_NAI"

#~/Documents/PhD/scripts/SHAPES/RNAprobBash-master/summarize_unique_barcodes2.sh -f ../aln_256cell_ctrl.sam.gz -b /Volumes/USELESS/DATA/Shape-Seq/preproc_256cell_ctrl/barcodes.txt -t -k -o "../summarize_256cell_ctrl"

#scp -r summarize_256cell_ctrl katchyz@login.ii.uib.no:/Home/ii/katchyz/OUT/SHAPES_out/

~/Documents/PhD/scripts/SHAPES/RNAprobBash-master/summarize_unique_barcodes2.sh -f ../aln_256cell_NAI.sam.gz -b /Volumes/USELESS/DATA/Shape-Seq/preproc_256cell_NAI/barcodes.txt -t -k -o "../summarize_256cell_NAI"

scp -r summarize_256cell_NAI katchyz@login.ii.uib.no:/Home/ii/katchyz/OUT/SHAPES_out/

### normalization (RNAprobR)