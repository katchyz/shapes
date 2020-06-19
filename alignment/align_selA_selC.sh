### SHAPES, December 2017

## Sel-A (old compound)
## Sel-C (new compound)

## R1: trim 3 last; R2: trim 15 last
## cut adapters
## preprocess (save barcodes, trim 3'nts)
## align (tophat, 5 mm)

##########
id="Sel-A"
r1="Sel-A_S1_L001_R1_001.fastq.gz"
r2="Sel-A_S1_L001_R2_001.fastq.gz"

id="Sel-C"
r1="Sel-C_S2_L001_R1_001.fastq.gz"
r2="Sel-C_S2_L001_R2_001.fastq.gz"


shapes="/export/valenfs/data/raw_data/SHAPE/Shapes_Dec2017"
bowtie_index='/Home/ii/katchyz/DATA/zebrafish_GRCz10/bowtie2_index/GRCz10_bowtie'
tophat_trans='/Home/ii/katchyz/DATA/zebrafish_GRCz10/tophat2_transcriptome/Danio_rerio.GRCz10.84.chr'

cd $shapes
mkdir temp_${id}
cd temp_${id}
cutadapt -q 20 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 20 -o r1.fastq -p r2.fastq ${shapes}/raw/$r1 ${shapes}/raw/$r2
python ../preprocess.py r1.fastq r2.fastq

### five mismatches
out=${shapes}/PAIRED_${id}
tophat -p 8 --no-coverage-search -N 7 --read-edit-dist 9 --transcriptome-index=$tophat_trans -o $out $bowtie_index read1.fastq read2.fastq

## normalize (2-step): check if '-' strand is ok with Rle-vectors
### TUN_paired
### TUN_into_GR


python ../preprocess2.py r1.fastq r2.fastq

### five mismatches
out=${shapes}/PAIRED2_${id}
tophat -p 8 --no-coverage-search -N 7 --read-edit-dist 9 --transcriptome-index=$tophat_trans -o $out $bowtie_index read1_2.fastq read2_2.fastq



