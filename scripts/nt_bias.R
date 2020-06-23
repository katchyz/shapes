## check nt bias on 5'ends of reads (after trimming untemplated nt)
library(seqinr)

## load fasta
#fasta_dna <- read.fasta("/Home/ii/katchyz/DATA/fasta/danio_rerio/GRCz10/dna/Danio_rerio.GRCz10.dna.toplevel.fa")
fasta_dna <- read.fasta("/Volumes/USELESS/DATA/fasta/danio_rerio/GRCz10/dna/Danio_rerio.GRCz10.dna.toplevel.fa")

## load norm1, for each read:
#tun_file <- "/export/valenfs/data/raw_data/SHAPE/final_results/save_rt_granges/norm1_cell1_invitro_NAIN3.Rsave"

libs_path <- "/Volumes/USELESS/DATA/SHAPES/SE"

libs <- c("norm1_cell1_invitro_NAIN3.Rsave", "norm1_cell24_invivo_DMSO.Rsave", "norm1_cell24_invivo_NAI.Rsave",
          "norm1_cell256_invivo_DMSO.Rsave", "norm1_cell256_invivo_NAI.Rsave", "norm1_oblong_CHX_invivo_DMSO.Rsave",
          "norm1_oblong_CHX_invivo_NAI.Rsave", "norm1_oblong_invivo_DMSO.Rsave", "norm1_oblong_invivo_NAI.Rsave",
          "norm1_sphere_invitro_A.Rsave", "norm1_sphere_invitro_C.Rsave")

## get start from fasta for plus strand, end for minus strand
## (multiply by number of unique_barcodes)?
## sum As, Ts, Cs, Gs...
for (i in 1:length(libs)) {
  print(libs[i])
  tun <- get(load(file = c(file.path(libs_path, libs[i]))))
  plus <- tun[strand(tun) == "+"]
  minus <- tun[strand(tun) == "-"]
  plus_start <- as.character(mapply(function(x,y){fasta_dna[[x]][y]}, as.character(seqnames(plus)), start(plus)))
  plus_start <- rep(plus_start, as.numeric(plus$unique_barcodes))
  minus_end <- as.character(mapply(function(x,y){fasta_dna[[x]][y]}, as.character(seqnames(minus)), end(minus)))
  minus_end <- rep(minus_end, as.numeric(minus$unique_barcodes))
  nt <- c(plus_start, minus_end)
  nt <- round(table(nt) / sum(table(nt)), digits=2)
  print(nt)
}


