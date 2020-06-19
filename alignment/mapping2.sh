### processing SHAPE-Seq data


# ligation adapter with random barcode: NNNNNNNAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
# PCR forward primer: AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCT
# index primer: CAAGCAGAAGACGGCATACGAGAT>>TACAAGG<<TGACTGGAGTTCAGACGTGTGCTCTTCCGATCT (barcode for index 12)
############### CAAGCAGAAGACGGCATACGAGAT  AAGCTAG  TGACTGGAGTTCAGACGTGTGCTCTTCCGATCT (index 10)

# libraries
raw="/export/valenfs/data/raw_data/SHAPE-Seq/raw"

# 1-2-1KCell-invivo-DMSO_S1_R1_001.fastq.gz # cell1K_invivo_DMSO_et
# 2-13-1Cell-invivo-DMSO_S2_R1_001.fastq.gz # cell1_invivo_DMSO_et1
# 3-14-1Cell-invivo-DMSO_S3_R1_001.fastq.gz # cell1_invivo_DMSO_et2
# 4-25-1Cell-invitro-NAIN3_S4_R1_001.fastq.gz # cell1_invitro_NAI_et1
# 5-26-1Cell-invitro-NAIN3_S5_R1_001.fastq.gz # cell1_invitro_NAI_et2
# 6-38-1Cell-invivo-NAIN3_S6_R1_001.fastq.gz # cell1_invivo_NAI_et1
# 7-39-1Cell-invivo-NAIN3_S7_R1_001.fastq.gz # cell1_invivo_NAI_et2
# 8-44-1KCell-invivo-NAIN3_S8_R1_001.fastq.gz # cell1K_invivo_NAI_et1
# 9-45-1KCell-invivo-NAIN3_S9_R1_001.fastq.gz # cell1K_invivo_NAI_et2


## cutting adapters (cutadapt)
cut="/export/valenfs/data/raw_data/SHAPE-Seq/cut"

#cutadapt -a ADAPTER [options] [-o output.fastq] input.fastq
cutadapt -q 20 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 20 -o ${cut}/cell1K_invivo_DMSO_et.fastq.gz ${raw}/1-2-1KCell-invivo-DMSO_S1_R1_001.fastq.gz
cutadapt -q 20 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 20 -o ${cut}/cell1_invivo_DMSO_et1.fastq.gz ${raw}/2-13-1Cell-invivo-DMSO_S2_R1_001.fastq.gz
cutadapt -q 20 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 20 -o ${cut}/cell1_invivo_DMSO_et2.fastq.gz ${raw}/3-14-1Cell-invivo-DMSO_S3_R1_001.fastq.gz
cutadapt -q 20 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 20 -o ${cut}/cell1_invitro_NAI_et1.fastq.gz ${raw}/4-25-1Cell-invitro-NAIN3_S4_R1_001.fastq.gz
cutadapt -q 20 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 20 -o ${cut}/cell1_invitro_NAI_et2.fastq.gz ${raw}/5-26-1Cell-invitro-NAIN3_S5_R1_001.fastq.gz
cutadapt -q 20 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 20 -o ${cut}/cell1_invivo_NAI_et1.fastq.gz ${raw}/6-38-1Cell-invivo-NAIN3_S6_R1_001.fastq.gz
cutadapt -q 20 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 20 -o ${cut}/cell1_invivo_NAI_et2.fastq.gz ${raw}/7-39-1Cell-invivo-NAIN3_S7_R1_001.fastq.gz
cutadapt -q 20 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 20 -o ${cut}/cell1K_invivo_NAI_et1.fastq.gz ${raw}/8-44-1KCell-invivo-NAIN3_S8_R1_001.fastq.gz
cutadapt -q 20 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 20 -o ${cut}/cell1K_invivo_NAI_et2.fastq.gz ${raw}/9-45-1KCell-invivo-NAIN3_S9_R1_001.fastq.gz


######
raw="/export/valenfs/data/raw_data/SHAPE-Seq/raw"
cut="/export/valenfs/data/raw_data/SHAPE-Seq/cut"
cutadapt -q 20 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 20 -o ${cut}/cell1K_invivo_NAI_et1.fastq.gz ${raw}/8-44-1KCell-invivo-NAIN3_S8_R1_001.fastq.gz
## preprocessing (preprocess.sh from RNAprobBash-master)
#preprocessing.sh -1 <READ1> -2 <READ2> -b <BARCODE_SEQ> -t <TRIM_LENGTH> -o <output_dir>
id="cell1K_invivo_NAI_et1"

gunzip ${cut}/${id}.fastq.gz

mkdir temp_${id}
cd temp_${id}
~/scripts/preprocessing.sh -1 ${cut}/${id}.fastq -b NNNNNNN -t 15 -o "preproc_${id}"
## mapping
shapes="/export/valenfs/data/raw_data/SHAPE-Seq"
bowtie_index='/Home/ii/katchyz/DATA/zebrafish_GRCz10/bowtie2_index/GRCz10_bowtie'
tophat_trans='/Home/ii/katchyz/DATA/zebrafish_GRCz10/tophat2_transcriptome/Danio_rerio.GRCz10.84.chr'
out=${shapes}/${id}
r1trimmed=${shapes}/temp_${id}/preproc_${id}/read1.fastq
tophat -p 8 --transcriptome-index=$tophat_trans -o $out $bowtie_index $r1trimmed
# change bam to sam.gz
samtools view -h -o ${shapes}/${id}/${id}.sam ${shapes}/${id}/accepted_hits.bam
gzip ${shapes}/${id}/${id}.sam

### run locally summarize_unique_barcodes2
### (% protein coding)
### normalization (RNAprobR)

#cell1K_invivo_DMSO_et
#cell1K_invivo_NAI_et1
#cell1K_invivo_NAI_et2
#cell1_invitro_NAI_et1
#cell1_invitro_NAI_et2
#cell1_invivo_DMSO_et1
#cell1_invivo_DMSO_et2
#cell1_invivo_NAI_et1
#cell1_invivo_NAI_et2

id="cell1_invivo_NAI_et2"
cd preproc_${id}
~/Documents/PhD/scripts/SHAPES/RNAprobBash-master/summarize_unique_barcodes2.sh -f ../${id}.sam.gz -b /Volumes/USELESS/test2/preproc_${id}/barcodes.txt -t -k -o "../summarize_${id}"


