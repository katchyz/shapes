library(RNAprobR)
library(GenomicRanges)
library(GenomicFeatures)
##### import readsamples.R #####
##### import comp.R #####

# cell1K_invivo_DMSO_et
# cell1K_invivo_NAI_et1
# cell1K_invivo_NAI_et2
# cell1_invitro_NAI_et1
# cell1_invitro_NAI_et2
# cell1_invivo_DMSO_et1
# cell1_invivo_DMSO_et2
# cell1_invivo_NAI_et1
# cell1_invivo_NAI_et2

ctr <- "cell1_invivo_DMSO_et1"
tre <- "cell1_invivo_NAI_et1"
nor <- "cell1_invivo_et1"

###########################################
###########################################

# reading in datasets
samples_path <- "/export/valenfs/data/raw_data/SHAPE-Seq"
control <- c(file.path(samples_path, paste("summarize_", ctr, sep=""), "unique_barcodes.txt"))
treated <- c(file.path(samples_path, paste("summarize_", tre, sep=""), "unique_barcodes.txt"))

##### import readsamples.R #####

control_counts <- readsamples(control, euc="counts")
treated_counts <- readsamples(treated, euc="counts")

txdb_can <- makeTxDbFromGFF("/Home/ii/katchyz/DATA/genomes/GTF/Danio_rerio.GRCz10.84.gtf", format = "gtf")
exons <- exonsBy(txdb_can, by="tx", use.names=TRUE)

control_tx <- mapToTranscripts(control_counts, exons, ignore.strand = FALSE)
control_tx$EUC <- control_counts[control_tx$xHits]$EUC
#control_tx <- split(control_tx, seqnames(control_tx))

treated_tx <- mapToTranscripts(treated_counts, exons)
treated_tx$EUC <- treated_counts[treated_tx$xHits]$EUC

##### import comp.R #####
# compiling positional data
control_comp <- comp(control_tx, cutoff=1, exons)
treated_comp <- comp(treated_tx, cutoff=1, exons)

# save...
save(control_comp, file=paste("/export/valenfs/data/raw_data/SHAPE-Seq/normalized/RAW_", ctr, ".Rsave", sep = ""))
save(treated_comp, file=paste("/export/valenfs/data/raw_data/SHAPE-Seq/normalized/RAW_", tre, ".Rsave", sep = ""))


# ???
cc <- split(control_comp, seqnames(control_comp))
tc <- split(treated_comp, seqnames(treated_comp))

## add to control ones that exist in treated (with zeros)
toadd <- unlist(tc[!(names(tc) %in% names(cc))])
toadd$TC <- rep(0, length(toadd$TC))
toadd <- split(toadd, seqnames(toadd))
toadd <- toadd[sapply(toadd, function(x){length(x) > 0})]

cc <- c(cc, toadd)
cc <- cc[names(cc) %in% names(tc)]
cc <- cc[order(names(cc))]

ucc <- unlist(cc)
ucc$TC.control <- ucc$TC
ucc$TC.treated <- unlist(tc)$TC
ucc$TC <- NULL
P <- 1
ucc$log2ratio <- log2(ucc$TC.treated+P) - log2(ucc$TC.control+P)

## set negative log2ratio to 0
ucc$log2ratio[ucc$log2ratio < 0] <- 0

shapeseq_norm <- split(ucc, seqnames(ucc))
shapeseq_norm <- shapeseq_norm[sapply(shapeseq_norm, function(x){length(x) > 0})]

##
.process_oneRNA_off <- function(rna, df){
  rna@elementMetadata <- rbind(rna@elementMetadata[2:nrow(rna@elementMetadata),], df)
  rna
}

df <- data.frame(TC.control = 0, TC.treated = 0, log2ratio = 0)
shapeseq_norm <- endoapply(shapeseq_norm, FUN=.process_oneRNA_off, df)

# save
save(shapeseq_norm, file=paste("/export/valenfs/data/raw_data/SHAPE-Seq/normalized/", nor, ".Rsave", sep = ""))

### map to genomic coordinates

gc_shape <- mapFromTranscripts(ucc, exons, ignore.strand = FALSE)
# save
save(gc_shape, file=paste("/export/valenfs/data/raw_data/SHAPE-Seq/normalized/gc_", nor, ".Rsave", sep = ""))



