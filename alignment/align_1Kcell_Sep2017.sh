### Shape-Seq August 2017


# 1Kcell, RNA control
# 1Kcell, in vitro, DMSO, RZ
# 1Kcell, in vitro, NAI, RZ
# 1Kcell, DMSO, RZ
# 1Kcell, NAI, RZ
# 1Kcell, in vitro, DMSO, polyA
# 1Kcell, in vitro, NAI, polyA
# 1Kcell, DMSO, polyA
# 1Kcell, NAI, polyA

# ligation adapter with random barcode: NNNNNNNAGATC GTCGTGTAGGGAAAGAGTGT
# PCR forward primer: AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCT
# index primer: CAAGCAGAAGACGGCATACGAGAT>>TACAAGG<<TGACTGGAGTTCAGACGTGTGCTCTTCCGATCT (barcode for index 12)
############### CAAGCAGAAGACGGCATACGAGAT  AAGCTAG  TGACTGGAGTTCAGACGTGTGCTCTTCCGATCT (index 10)

### cut adapters (why are the reads so long???)
# cutadapt -a ADAPT1 -A ADAPT2 [options] -o out1.fastq -p out2.fastq in1.fastq in2.fastq
# ligation adapter with random barcode: NNNNNNNAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
# libraries
# 1
raw="/export/valenfs/data/raw_data/Shape-Seq_Aug2017/raw"
cut="/export/valenfs/data/raw_data/Shape-Seq_Aug2017/cut"

r1="1Kcell_polyA_DMSO_R1.fastq.gz"
r2="1Kcell_polyA_DMSO_R2.fastq.gz"
cutadapt -q 20 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 20 -o ${cut}/$r1 -p ${cut}/$r2 ${raw}/$r1 ${raw}/$r2

# 2
raw="/export/valenfs/data/raw_data/Shape-Seq_Aug2017/raw"
cut="/export/valenfs/data/raw_data/Shape-Seq_Aug2017/cut"

r1="1Kcell_polyA_NAI_R1.fastq.gz"
r2="1Kcell_polyA_NAI_R2.fastq.gz"
cutadapt -q 20 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 20 -o ${cut}/$r1 -p ${cut}/$r2 ${raw}/$r1 ${raw}/$r2

# 3
raw="/export/valenfs/data/raw_data/Shape-Seq_Aug2017/raw"
cut="/export/valenfs/data/raw_data/Shape-Seq_Aug2017/cut"

r1="1Kcell_polyA_vitro_DMSO_R1.fastq.gz"
r2="1Kcell_polyA_vitro_DMSO_R2.fastq.gz"
cutadapt -q 20 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 20 -o ${cut}/$r1 -p ${cut}/$r2 ${raw}/$r1 ${raw}/$r2

# 4
raw="/export/valenfs/data/raw_data/Shape-Seq_Aug2017/raw"
cut="/export/valenfs/data/raw_data/Shape-Seq_Aug2017/cut"

r1="1Kcell_polyA_vitro_NAI_R1.fastq.gz"
r2="1Kcell_polyA_vitro_NAI_R2.fastq.gz"
cutadapt -q 20 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 20 -o ${cut}/$r1 -p ${cut}/$r2 ${raw}/$r1 ${raw}/$r2

# 5
raw="/export/valenfs/data/raw_data/Shape-Seq_Aug2017/raw"
cut="/export/valenfs/data/raw_data/Shape-Seq_Aug2017/cut"

r1="1Kcell_RZ_DMSO_R1.fastq.gz"
r2="1Kcell_RZ_DMSO_R2.fastq.gz"
cutadapt -q 20 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 20 -o ${cut}/$r1 -p ${cut}/$r2 ${raw}/$r1 ${raw}/$r2

# 6
raw="/export/valenfs/data/raw_data/Shape-Seq_Aug2017/raw"
cut="/export/valenfs/data/raw_data/Shape-Seq_Aug2017/cut"

r1="1Kcell_RZ_NAI_R1.fastq.gz"
r2="1Kcell_RZ_NAI_R2.fastq.gz"
cutadapt -q 20 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 20 -o ${cut}/$r1 -p ${cut}/$r2 ${raw}/$r1 ${raw}/$r2

# 7
raw="/export/valenfs/data/raw_data/Shape-Seq_Aug2017/raw"
cut="/export/valenfs/data/raw_data/Shape-Seq_Aug2017/cut"

r1="1Kcell_RZ_RNAcontrol_R1.fastq.gz"
r2="1Kcell_RZ_RNAcontrol_R2.fastq.gz"
cutadapt -q 20 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 20 -o ${cut}/$r1 -p ${cut}/$r2 ${raw}/$r1 ${raw}/$r2

# 8
raw="/export/valenfs/data/raw_data/Shape-Seq_Aug2017/raw"
cut="/export/valenfs/data/raw_data/Shape-Seq_Aug2017/cut"

r1="1Kcell_RZ_vitro_DMSO_R1.fastq.gz"
r2="1Kcell_RZ_vitro_DMSO_R2.fastq.gz"
cutadapt -q 20 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 20 -o ${cut}/$r1 -p ${cut}/$r2 ${raw}/$r1 ${raw}/$r2

# 9
raw="/export/valenfs/data/raw_data/Shape-Seq_Aug2017/raw"
cut="/export/valenfs/data/raw_data/Shape-Seq_Aug2017/cut"

r1="1Kcell_RZ_vitro_NAI_R1.fastq.gz"
r2="1Kcell_RZ_vitro_NAI_R2.fastq.gz"
cutadapt -q 20 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 20 -o ${cut}/$r1 -p ${cut}/$r2 ${raw}/$r1 ${raw}/$r2


###################
### preprocessing

## preprocessing (preprocess.sh from RNAprobBash-master)
#preprocessing.sh -1 <READ1> -2 <READ2> -b <BARCODE_SEQ> -t <TRIM_LENGTH> -o <output_dir>


#id="1Kcell_polyA_DMSO" 		#1 YES pp
#id="1Kcell_polyA_NAI" 			#2 YES pp
#id="1Kcell_polyA_vitro_DMSO" 	#3 YES pp
#id="1Kcell_polyA_vitro_NAI" 	#4 YES fix # screen -r 6972.pts-14.kjempetuja pp
#id="1Kcell_RZ_DMSO" 			#5 YES
#id="1Kcell_RZ_NAI" 			#6 YES pp
#id="1Kcell_RZ_RNAcontrol" 		#7 YES pp
#id="1Kcell_RZ_vitro_DMSO" 		#8 YES pp
#id="1Kcell_RZ_vitro_NAI" 		#9 YES


mkdir temp_${id}
cd temp_${id}
cp ${cut}/${id}* .
gunzip *.fastq.gz

~/scripts/preprocessing.sh -1 ${id}_R1.fastq -2 ${id}_R2.fastq -b NNNNNNN -t 15 -o preproc_${id}


###################
### mapping

## mapping
shapes="/export/valenfs/data/raw_data/Shape-Seq_Aug2017"
bowtie_index='/Home/ii/katchyz/DATA/zebrafish_GRCz10/bowtie2_index/GRCz10_bowtie'
tophat_trans='/Home/ii/katchyz/DATA/zebrafish_GRCz10/tophat2_transcriptome/Danio_rerio.GRCz10.84.chr'

out=${shapes}/${id}
r1trimmed=${shapes}/temp_${id}/preproc_${id}/read1.fastq
r2trimmed=${shapes}/temp_${id}/preproc_${id}/read2.fastq
tophat -p 8 --transcriptome-index=$tophat_trans -o $out $bowtie_index $r1trimmed $r2trimmed


### ~/bbmap/repair.sh in1=read1.fastq in2=read2.fastq out1=fixed1.fastq out2=fixed2.fastq outsingle=singletons.fastq


####### FIXED

out=${shapes}/fixed_${id}
r1trimmed=${shapes}/temp_${id}/preproc_${id}/fixed1.fastq
r2trimmed=${shapes}/temp_${id}/preproc_${id}/fixed2.fastq
tophat -p 8 --no-coverage-search --transcriptome-index=$tophat_trans -o $out $bowtie_index $r1trimmed $r2trimmed


#### preproc just 1st read
cd temp_${id}
~/scripts/preprocessing.sh -1 ${id}_R1.fastq -b NNNNNNN -t 15 -o pp_${id}


##################
##### simple cutting of barcodes

#id="1Kcell_polyA_DMSO" 			# 35097.pts-14.kjempetuja
#id="1Kcell_polyA_NAI" 				# 59755.pts-14.kjempetuja
#id="1Kcell_polyA_vitro_DMSO"		# 71558.pts-14.kjempetuja
#id="1Kcell_polyA_vitro_NAI" 		# 77033.pts-14.kjempetuja
#id="1Kcell_RZ_DMSO" 				# 79589.pts-14.kjempetuja
#id="1Kcell_RZ_NAI" 				# 80812.pts-14.kjempetuja	
#id="1Kcell_RZ_RNAcontrol" 			# 2046.pts-14.kjempetuja
#id="1Kcell_RZ_vitro_DMSO" 			# 3120.pts-14.kjempetuja
#id="1Kcell_RZ_vitro_NAI" 			# 4881.pts-14.kjempetuja

cd temp_${id}
python ../trim_barcodes.py ${id}_R1.fastq

shapes="/export/valenfs/data/raw_data/Shape-Seq_Aug2017"
bowtie_index='/Home/ii/katchyz/DATA/zebrafish_GRCz10/bowtie2_index/GRCz10_bowtie'
tophat_trans='/Home/ii/katchyz/DATA/zebrafish_GRCz10/tophat2_transcriptome/Danio_rerio.GRCz10.84.chr'

out=${shapes}/out_${id}
r1=${shapes}/temp_${id}/read1.fastq
r2=${shapes}/temp_${id}/${id}_R2.fastq
tophat -p 8 --no-coverage-search --transcriptome-index=$tophat_trans -o $out $bowtie_index $r1 $r2



####### dovetailing reads?

#id="1Kcell_polyA_DMSO" 	
#id="1Kcell_polyA_NAI" 		
#id="1Kcell_polyA_vitro_DMSO"
#id="1Kcell_polyA_vitro_NAI" 
#id="1Kcell_RZ_DMSO" 		
#id="1Kcell_RZ_NAI" 		
#id="1Kcell_RZ_RNAcontrol" 	
#id="1Kcell_RZ_vitro_DMSO" 	
#id="1Kcell_RZ_vitro_NAI" 	

shapes="/export/valenfs/data/raw_data/Shape-Seq_Aug2017"
bowtie_index='/Home/ii/katchyz/DATA/zebrafish_GRCz10/bowtie2_index/GRCz10_bowtie'
tophat_trans='/Home/ii/katchyz/DATA/zebrafish_GRCz10/tophat2_transcriptome/Danio_rerio.GRCz10.84.chr'

out=${shapes}/dovetail_${id}
r1trimmed=${shapes}/temp_${id}/preproc_${id}/read1.fastq
r2trimmed=${shapes}/temp_${id}/preproc_${id}/read2.fastq
tophat -p 8 --no-coverage-search --dovetail --transcriptome-index=$tophat_trans -o $out $bowtie_index $r1trimmed $r2trimmed

####### 

shapes="/export/valenfs/data/raw_data/Shape-Seq_Aug2017"
bowtie_index='/Home/ii/katchyz/DATA/zebrafish_GRCz10/bowtie2_index/GRCz10_bowtie'
tophat_trans='/Home/ii/katchyz/DATA/zebrafish_GRCz10/tophat2_transcriptome/Danio_rerio.GRCz10.84.chr'

out=${shapes}/tophat2_${id}
r1trimmed=${shapes}/temp_${id}/preproc_${id}/read1.fastq
r2trimmed=${shapes}/temp_${id}/preproc_${id}/read2.fastq
tophat -p 8 --no-coverage-search --transcriptome-index=$tophat_trans -o $out $bowtie_index $r1trimmed $r2trimmed




########## preprocess

python ../preprocess.py ${id}_R1.fastq ${id}_R2.fastq


out=${shapes}/FIXED_${id}
r1trimmed=${shapes}/temp_${id}/read1.fastq
r2trimmed=${shapes}/temp_${id}/read2.fastq
tophat -p 8 --no-coverage-search --transcriptome-index=$tophat_trans -o $out $bowtie_index $r1trimmed $r2trimmed


############# preprocessed again, trimming a lot of the 2nd read
cd temp_${id}
python ../preprocess.py ${id}_R1.fastq ${id}_R2.fastq


cd $shapes
out=${shapes}/OUT_${id}
r1trimmed=${shapes}/temp_${id}/read1.fastq
r2trimmed=${shapes}/temp_${id}/read2.fastq
tophat -p 8 --no-coverage-search --transcriptome-index=$tophat_trans -o $out $bowtie_index $r1trimmed $r2trimmed


####### SINGLE
out=${shapes}/SINGLE_${id}
r1trimmed=${shapes}/temp_${id}/read1.fastq
tophat -p 8 --no-coverage-search --transcriptome-index=$tophat_trans -o $out $bowtie_index $r1trimmed

###############################################################
###############################################################
cd temp_${id}
cutadapt -q 20 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 20 -o r1.fastq -p r2.fastq ${id}_R1.fastq ${id}_R2.fastq
python ../preprocess.py r1.fastq r2.fastq

### five mismatches
out=${shapes}/PAIRED_${id}
tophat -p 8 --no-coverage-search -N 5 --read-edit-dist 7 --transcriptome-index=$tophat_trans -o $out $bowtie_index read1.fastq read2.fastq


### remap 2-4 cell & 256cell
id="cell24_DMSO"
id="cell24_NAI"
id="cell256_DMSO"
id="cell256_NAI"

shapes="/export/valenfs/data/raw_data/Shape-Seq_Aug2017/remap_old"
bowtie_index='/Home/ii/katchyz/DATA/zebrafish_GRCz10/bowtie2_index/GRCz10_bowtie'
tophat_trans='/Home/ii/katchyz/DATA/zebrafish_GRCz10/tophat2_transcriptome/Danio_rerio.GRCz10.84.chr'

cd $shapes
mkdir temp_${id}
cd temp_${id}
cutadapt -q 20 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 20 -o r1.fastq -p r2.fastq ../${id}_R1.fastq ../${id}_R2.fastq
python ../preprocess.py r1.fastq r2.fastq

### five mismatches
out=${shapes}/PAIRED_${id}
tophat -p 8 --no-coverage-search -N 5 --read-edit-dist 7 --transcriptome-index=$tophat_trans -o $out $bowtie_index read1.fastq read2.fastq

