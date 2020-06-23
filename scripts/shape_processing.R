library(RNAprobR)
library(GenomicRanges)

# reading in datasets
samples_path = "/export/valenfs/data/raw_data/SHAPES/"
control = c(file.path(samples_path, "summarize_256cell_ctrl", "unique_barcodes.txt"))
treated = c(file.path(samples_path, "summarize_256cell_NAI", "unique_barcodes.txt"))

control_counts = readsamples(control, euc="counts")
treated_counts = readsamples(treated, euc="counts")

txdb_can <- makeTxDbFromGFF("/Home/ii/katchyz/DATA/genomes/GTF/Danio_rerio.GRCz10.84.gtf", format = "gtf")
exons <- exonsBy(txdb_can, by="tx", use.names=TRUE)

control_tx <- mapToTranscripts(control_counts, exons, ignore.strand = FALSE)
control_tx$EUC <- control_counts[control_tx$xHits]$EUC
#control_tx <- split(control_tx, seqnames(control_tx))

treated_tx <- mapToTranscripts(treated_counts, exons)
treated_tx$EUC <- treated_counts[treated_tx$xHits]$EUC

# compiling positional data
control_comp = comp(control_tx, cutoff=1, exons)
treated_comp = comp(treated_tx, cutoff=1, exons)

# save...
save(control_comp, file="/export/valenfs/data/raw_data/SHAPES/normalized/RAW_256cell_ctrl.Rsave")
save(treated_comp, file="/export/valenfs/data/raw_data/SHAPES/normalized/RAW_256cell_NAI.Rsave")

cc <- split(control_comp, seqnames(control_comp))
tc <- split(treated_comp, seqnames(treated_comp))

## add to control ones that exist in treated (with zeros)
toadd <- unlist(tc[!(names(tc) %in% names(cc))])
toadd$TC <- rep(0, length(toadd$TC))
c
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

shapeseq_norm <- split(ucc, seqnames(ucc))
shapeseq_norm <- shapeseq_norm[sapply(shapeseq_norm, function(x){length(x) > 0})]

###########################

.process_oneRNA_off <- function(rna, df){
  rna@elementMetadata <- rbind(rna@elementMetadata[2:nrow(rna@elementMetadata),], df)
  rna
}


df <- data.frame(TC.control = 0, TC.treated = 0, log2ratio = 0)
shapeseq_norm <- endoapply(shapeseq_norm, FUN=.process_oneRNA_off, df)



# save
save(shapeseq_norm, file="/export/valenfs/data/raw_data/SHAPES/normalized/256cell.Rsave")


