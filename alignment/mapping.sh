### processing SHAPES data


# ligation adapter with random barcode: NNNNNNNAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
# PCR forward primer: AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCT
# index primer: CAAGCAGAAGACGGCATACGAGAT>>TACAAGG<<TGACTGGAGTTCAGACGTGTGCTCTTCCGATCT (barcode for index 12)
############### CAAGCAGAAGACGGCATACGAGAT  AAGCTAG  TGACTGGAGTTCAGACGTGTGCTCTTCCGATCT (index 10)

# libraries
raw="/export/valenfs/data/raw_data/SHAPES/raw"

# index_12_S4_R1_001.fastq.gz # oblong_invivo_DMSO
# index_13_S5_R1_001.fastq.gz # oblong_invivo_NAI
# index_14_S6_R1_001.fastq.gz # oblong_CHX_invivo_DMSO
# index_15_S7_R1_001.fastq.gz # oblong_CHX_invivo_NAI
# index_1_S1_R1_001.fastq.gz # cell2_4_invivo_NAI_N3
# index_6_S2_R1_001.fastq.gz # cell1_invitro_NAI_N3_nonsel
# index_9_S3_R1_001.fastq.gz # cell1_invitro_NAI_N3

## cutting adapters (cutadapt)
cut="/export/valenfs/data/raw_data/SHAPES/cut"

#cutadapt -a ADAPTER [options] [-o output.fastq] input.fastq
cutadapt -q 20 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 20 -o ${cut}/oblong_invivo_DMSO.fastq.gz ${raw}/index_12_S4_R1_001.fastq.gz
cutadapt -q 20 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 20 -o ${cut}/oblong_invivo_NAI.fastq.gz ${raw}/index_13_S5_R1_001.fastq.gz
cutadapt -q 20 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 20 -o ${cut}/oblong_CHX_invivo_DMSO.fastq.gz ${raw}/index_14_S6_R1_001.fastq.gz
cutadapt -q 20 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 20 -o ${cut}/oblong_CHX_invivo_NAI.fastq.gz ${raw}/index_15_S7_R1_001.fastq.gz
cutadapt -q 20 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 20 -o ${cut}/cell2_4_invivo_NAI_N3.fastq.gz ${raw}/index_1_S1_R1_001.fastq.gz
cutadapt -q 20 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 20 -o ${cut}/cell1_invitro_NAI_N3_nonsel.fastq.gz ${raw}/index_6_S2_R1_001.fastq.gz
cutadapt -q 20 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 20 -o ${cut}/cell1_invitro_NAI_N3.fastq.gz ${raw}/index_9_S3_R1_001.fastq.gz

### cutting extra primers etc.
raw="/export/valenfs/data/raw_data/SHAPES/raw"
trimmed="/export/valenfs/data/raw_data/SHAPES/trimmed"
# index 12
cutadapt -a CAAGCAGAAGACGGCATACGAGATTACAAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA -o ${trimmed}/1_oblong_invivo_DMSO.fastq.gz ${raw}/index_12_S4_R1_001.fastq.gz
cutadapt -a AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCT -o ${trimmed}/2_oblong_invivo_DMSO.fastq.gz ${trimmed}/1_oblong_invivo_DMSO.fastq.gz
cutadapt -q 20 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 20 -o ${trimmed}/oblong_invivo_DMSO.fastq.gz ${trimmed}/2_oblong_invivo_DMSO.fastq.gz
rm ${trimmed}/1_oblong_invivo_DMSO.fastq.gz ${trimmed}/2_oblong_invivo_DMSO.fastq.gz

# index 13
raw="/export/valenfs/data/raw_data/SHAPES/raw"
trimmed="/export/valenfs/data/raw_data/SHAPES/trimmed"
cutadapt -a CAAGCAGAAGACGGCATACGAGATTTGACTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA -o ${trimmed}/1_oblong_invivo_NAI.fastq.gz ${raw}/index_13_S5_R1_001.fastq.gz
cutadapt -a AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCT -o ${trimmed}/2_oblong_invivo_NAI.fastq.gz ${trimmed}/1_oblong_invivo_NAI.fastq.gz
cutadapt -q 20 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 20 -o ${trimmed}/oblong_invivo_NAI.fastq.gz ${trimmed}/2_oblong_invivo_NAI.fastq.gz
rm ${trimmed}/1_oblong_invivo_NAI.fastq.gz ${trimmed}/2_oblong_invivo_NAI.fastq.gz

# index 14
raw="/export/valenfs/data/raw_data/SHAPES/raw"
trimmed="/export/valenfs/data/raw_data/SHAPES/trimmed"
cutadapt -a CAAGCAGAAGACGGCATACGAGATGGAACTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA -o ${trimmed}/1_oblong_CHX_invivo_DMSO.fastq.gz ${raw}/index_14_S6_R1_001.fastq.gz
cutadapt -a AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCT -o ${trimmed}/2_oblong_CHX_invivo_DMSO.fastq.gz ${trimmed}/1_oblong_CHX_invivo_DMSO.fastq.gz
cutadapt -q 20 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 20 -o ${trimmed}/oblong_CHX_invivo_DMSO.fastq.gz ${trimmed}/2_oblong_CHX_invivo_DMSO.fastq.gz
rm ${trimmed}/1_oblong_CHX_invivo_DMSO.fastq.gz ${trimmed}/2_oblong_CHX_invivo_DMSO.fastq.gz

# index 15
raw="/export/valenfs/data/raw_data/SHAPES/raw"
trimmed="/export/valenfs/data/raw_data/SHAPES/trimmed"
cutadapt -a CAAGCAGAAGACGGCATACGAGATTGACATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA -o ${trimmed}/1_oblong_CHX_invivo_NAI.fastq.gz ${raw}/index_15_S7_R1_001.fastq.gz
cutadapt -a AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCT -o ${trimmed}/2_oblong_CHX_invivo_NAI.fastq.gz ${trimmed}/1_oblong_CHX_invivo_NAI.fastq.gz
cutadapt -q 20 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 20 -o ${trimmed}/oblong_CHX_invivo_NAI.fastq.gz ${trimmed}/2_oblong_CHX_invivo_NAI.fastq.gz
rm ${trimmed}/1_oblong_CHX_invivo_NAI.fastq.gz ${trimmed}/2_oblong_CHX_invivo_NAI.fastq.gz

# index 1
raw="/export/valenfs/data/raw_data/SHAPES/raw"
trimmed="/export/valenfs/data/raw_data/SHAPES/trimmed"
cutadapt -a CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA -o ${trimmed}/1_cell2_4_invivo_NAI_N3.fastq.gz ${raw}/index_1_S1_R1_001.fastq.gz
cutadapt -a AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCT -o ${trimmed}/2_cell2_4_invivo_NAI_N3.fastq.gz ${trimmed}/1_cell2_4_invivo_NAI_N3.fastq.gz
cutadapt -q 20 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 20 -o ${trimmed}/cell2_4_invivo_NAI_N3.fastq.gz ${trimmed}/2_cell2_4_invivo_NAI_N3.fastq.gz
rm ${trimmed}/1_cell2_4_invivo_NAI_N3.fastq.gz ${trimmed}/2_cell2_4_invivo_NAI_N3.fastq.gz

# index 6
raw="/export/valenfs/data/raw_data/SHAPES/raw"
trimmed="/export/valenfs/data/raw_data/SHAPES/trimmed"
cutadapt -a CAAGCAGAAGACGGCATACGAGATATTGGCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA -o ${trimmed}/1_cell1_invitro_NAI_N3_nonsel.fastq.gz ${raw}/index_6_S2_R1_001.fastq.gz
cutadapt -a AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCT -o ${trimmed}/2_cell1_invitro_NAI_N3_nonsel.fastq.gz ${trimmed}/1_cell1_invitro_NAI_N3_nonsel.fastq.gz
cutadapt -q 20 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 20 -o ${trimmed}/cell1_invitro_NAI_N3_nonsel.fastq.gz ${trimmed}/2_cell1_invitro_NAI_N3_nonsel.fastq.gz
rm ${trimmed}/1_cell1_invitro_NAI_N3_nonsel.fastq.gz ${trimmed}/2_cell1_invitro_NAI_N3_nonsel.fastq.gz

# index 7 (former 9)
raw="/export/valenfs/data/raw_data/SHAPES/raw"
trimmed="/export/valenfs/data/raw_data/SHAPES/trimmed"
cutadapt -a CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA -o ${trimmed}/1_cell1_invitro_NAI_N3.fastq.gz ${raw}/index_9_S3_R1_001.fastq.gz
cutadapt -a AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCT -o ${trimmed}/2_cell1_invitro_NAI_N3.fastq.gz ${trimmed}/1_cell1_invitro_NAI_N3.fastq.gz
cutadapt -q 20 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 20 -o ${trimmed}/cell1_invitro_NAI_N3.fastq.gz ${trimmed}/2_cell1_invitro_NAI_N3.fastq.gz
rm ${trimmed}/1_cell1_invitro_NAI_N3.fastq.gz ${trimmed}/2_cell1_invitro_NAI_N3.fastq.gz

## preprocessing (preprocess.sh from RNAprobBash-master)
#preprocessing.sh -1 <READ1> -2 <READ2> -b <BARCODE_SEQ> -t <TRIM_LENGTH> -o <output_dir>

~/scripts/preprocessing.sh -1 ${cut}/oblong_invivo_DMSO.fastq -b NNNNNNN -t 15 -o "preproc_oblong_invivo_DMSO"						# 12
~/scripts/preprocessing.sh -1 ${cut}/oblong_invivo_NAI.fastq -b NNNNNNN -t 15 -o "preproc_oblong_invivo_NAI"							# 13
~/scripts/preprocessing.sh -1 ${cut}/oblong_CHX_invivo_DMSO.fastq -b NNNNNNN -t 15 -o "preproc_oblong_CHX_invivo_DMSO"				# 14
~/scripts/preprocessing.sh -1 ${cut}/oblong_CHX_invivo_NAI.fastq -b NNNNNNN -t 15 -o "preproc_oblong_CHX_invivo_NAI"					# 15
~/scripts/preprocessing.sh -1 ${cut}/cell2_4_invivo_NAI_N3.fastq -b NNNNNNN -t 15 -o "preproc_cell2_4_invivo_NAI_N3"					# 1
~/scripts/preprocessing.sh -1 ${cut}/cell1_invitro_NAI_N3_nonsel.fastq -b NNNNNNN -t 15 -o "preproc_cell1_invitro_NAI_N3_nonsel"		# 6
~/scripts/preprocessing.sh -1 ${cut}/cell1_invitro_NAI_N3.fastq -b NNNNNNN -t 15 -o "preproc_cell1_invitro_NAI_N3"					# 9


## preprocessing (preprocess.sh from RNAprobBash-master)
#preprocessing.sh -1 <READ1> -2 <READ2> -b <BARCODE_SEQ> -t <TRIM_LENGTH> -o <output_dir>

trimmed="/export/valenfs/data/raw_data/SHAPES/trimmed"
~/scripts/preprocessing.sh -1 ${trimmed}/oblong_invivo_DMSO.fastq -b NNNNNNN -t 15 -o "pp_oblong_invivo_DMSO"
trimmed="/export/valenfs/data/raw_data/SHAPES/trimmed"
~/scripts/preprocessing.sh -1 ${trimmed}/oblong_invivo_NAI.fastq -b NNNNNNN -t 15 -o "pp_oblong_invivo_NAI"
trimmed="/export/valenfs/data/raw_data/SHAPES/trimmed"
~/scripts/preprocessing.sh -1 ${trimmed}/oblong_CHX_invivo_DMSO.fastq -b NNNNNNN -t 15 -o "pp_oblong_CHX_invivo_DMSO"
trimmed="/export/valenfs/data/raw_data/SHAPES/trimmed"
~/scripts/preprocessing.sh -1 ${trimmed}/oblong_CHX_invivo_NAI.fastq -b NNNNNNN -t 15 -o "pp_oblong_CHX_invivo_NAI"
trimmed="/export/valenfs/data/raw_data/SHAPES/trimmed"
~/scripts/preprocessing.sh -1 ${trimmed}/cell2_4_invivo_NAI_N3.fastq -b NNNNNNN -t 15 -o "pp_cell2_4_invivo_NAI_N3"
trimmed="/export/valenfs/data/raw_data/SHAPES/trimmed"
~/scripts/preprocessing.sh -1 ${trimmed}/cell1_invitro_NAI_N3_nonsel.fastq -b NNNNNNN -t 15 -o "pp_cell1_invitro_NAI_N3_nonsel"
trimmed="/export/valenfs/data/raw_data/SHAPES/trimmed"
~/scripts/preprocessing.sh -1 ${trimmed}/cell1_invitro_NAI_N3.fastq -b NNNNNNN -t 15 -o "pp_cell1_invitro_NAI_N3"



#id="oblong_invivo_DMSO"
shapes="/export/valenfs/data/raw_data/SHAPES"
id="256cell_ctrl"
# map to transcriptome / genome
bowtie_index='/Home/ii/katchyz/DATA/zebrafish_GRCz10/bowtie2_index/GRCz10_bowtie'
tophat_trans='/Home/ii/katchyz/DATA/zebrafish_GRCz10/tophat2_transcriptome/Danio_rerio.GRCz10.84.chr'
out=${shapes}/${id}
r1trimmed=${shapes}/preproc_${id}/read1.fastq
tophat -p 8 --transcriptome-index=$tophat_trans -o $out $bowtie_index $r1trimmed


## where do the multimappers map to? is it mitochondrial?

## mapping 

shapes="/export/valenfs/data/raw_data/SHAPES"
bowtie2 -p 8 -X 700 -x ~/DATA/genomes/index_GRCz10 -U ${shapes}/preproc_oblong_invivo_DMSO/read1.fastq | gzip > ${shapes}/aligned/oblong_invivo_DMSO.sam.gz
bowtie2 -p 8 -X 700 -x ~/DATA/genomes/index_GRCz10 -U ${shapes}/preproc_oblong_invivo_NAI/read1.fastq | gzip > ${shapes}/aligned/oblong_invivo_NAI.sam.gz
bowtie2 -p 8 -X 700 -x ~/DATA/genomes/index_GRCz10 -U ${shapes}/preproc_oblong_CHX_invivo_DMSO/read1.fastq | gzip > ${shapes}/aligned/oblong_CHX_invivo_DMSO.sam.gz
bowtie2 -p 8 -X 700 -x ~/DATA/genomes/index_GRCz10 -U ${shapes}/preproc_oblong_CHX_invivo_NAI/read1.fastq | gzip > ${shapes}/aligned/oblong_CHX_invivo_NAI.sam.gz
bowtie2 -p 8 -X 700 -x ~/DATA/genomes/index_GRCz10 -U ${shapes}/preproc_cell2_4_invivo_NAI_N3/read1.fastq | gzip > ${shapes}/aligned/cell2_4_invivo_NAI_N3.sam.gz
bowtie2 -p 8 -X 700 -x ~/DATA/genomes/index_GRCz10 -U ${shapes}/preproc_cell1_invitro_NAI_N3_nonsel/read1.fastq | gzip > ${shapes}/aligned/cell1_invitro_NAI_N3_nonsel.sam.gz
bowtie2 -p 8 -X 700 -x ~/DATA/genomes/index_GRCz10 -U ${shapes}/preproc_cell1_invitro_NAI_N3/read1.fastq | gzip > ${shapes}/aligned/cell1_invitro_NAI_N3.sam.gz



## sumarize barcodes (summarize_unique_barcodes2 from RNAprobBash-master)
~/scripts/RNAprobBash-master/summarize_unique_barcodes.sh
summarize_unique_barcodes.sh -f <SAM_file> -b <BARCODES> -p <PRIMING_POSITION> -t -k -r <R_SCRIPT_PATH>

~/Documents/PhD/scripts/SHAPES/RNAprobBash-master/summarize_unique_barcodes2.sh -f ../aln_256cell_NAI.sam.gz -b /Volumes/USELESS/DATA/Shape-Seq/preproc_256cell_NAI/barcodes.txt -t -k -o "../summarize_256cell_NAI"

# change bam to sam.gz
id="oblong_invivo_DMSO"
samtools view -h -o ${id}/${id}.sam ${id}/accepted_hits.bam
gzip ${id}/${id}.sam


# run locally?

~/scripts/RNAprobBash-master/summarize_unique_barcodes.sh -f ${id}/${id}.sam.gz -b preproc_${id}/barcodes.txt -t -k -o summarize_${id}


id="cell1_invitro_NAI_N3"
~/Documents/PhD/scripts/SHAPES/RNAprobBash-master/summarize_unique_barcodes2.sh -f ${id}.sam.gz -b preproc_${id}/barcodes.txt -t -k -o summarize_${id}


## normalization (RNAprobR)





############## 18S
id="cell1_invitro_NAI_N3"
shapes="/export/valenfs/data/raw_data/SHAPES"
bowtie2 -p 8 -X 700 -x ~/DATA/zebrafish_18S -U ${shapes}/preproc_${id}/read1.fastq | gzip > ${shapes}/18S/${id}.sam.gz


