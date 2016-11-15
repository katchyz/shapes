#!/bin/sh
# 256 cell: 50 mM NAI
data_dir='/Volumes/USELESS/DATA/Shape-Seq'
out_dir='/Volumes/USELESS/OUT/SHAPES_out/256cell_50_mM_NAI'

bowtie_index='/Volumes/USELESS/DATA/zebrafish_GRCz10/bowtie2_index/GRCz10_bowtie'
tophat_trans='/Volumes/USELESS/DATA/zebrafish_GRCz10/tophat2_transcriptome/Danio_rerio.GRCz10.84.chr'

for lane in 1 2 3 4
do
	r1=${data_dir}/SHAPE-sample-6_S8_L00${lane}_R1_001.fastq.gz
	r2=${data_dir}/SHAPE-sample-6_S8_L00${lane}_R2_001.fastq.gz
	#
	tmp1=${r1%.fastq.gz}.TEMP.fastq.gz
	tmp2=${r2%.fastq.gz}.TEMP.fastq.gz
	#
	r1trimmed=${r1%.fastq.gz}.trimmed.fastq.gz
	r2trimmed=${r2%.fastq.gz}.trimmed.fastq.gz
	#
	# cut adapters and remove pairs shorter than 20 nt; trim 7N at 5' of R1 and 15N at 5' of R2
	cutadapt -q 20 -a AGATCGGAAGAGCACACGTCT -u7 --minimum-length 20 --paired-output $tmp2 -o $tmp1 $r1 $r2
	cutadapt -q 20 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -u 15 --minimum-length 20 --paired-output $r1trimmed -o $r2trimmed $tmp2 $tmp1
	rm $tmp1 $tmp2
	# map to transcriptome / genome
	out=${out_dir}/lane${lane}
	tophat2 -p 8 --transcriptome-index=$tophat_trans -o $out $bowtie_index $r1trimmed $r2trimmed
	#
	echo finished lane $lane
done
