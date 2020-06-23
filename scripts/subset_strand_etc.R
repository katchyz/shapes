### get transcript dtabase
txdb_can <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/danio_rerio/zebrafish_GRCz10.gtf", format = "gtf")
cds <- cdsBy(txdb_can, by="tx", use.names=TRUE)
fiveUTR <- fiveUTRsByTranscript(txdb_can, use.names=TRUE)
threeUTR <- threeUTRsByTranscript(txdb_can, use.names=TRUE)
exons <- exonsBy(txdb_can, by="tx", use.names=TRUE)
##############

load("/Volumes/USELESS/OUT/SHAPES_out/shapeseq_norm_dtcr_2-4cell.Rsave")
dtcr_2_4 <- shapeseq_norm
rm(shapeseq_norm)

exon_minus <- exons[strand(exons) == "-"]
strand(dtcr_2_4[seqnames(dtcr_2_4) %in% names(exon_minus)]) <- rep("-", length(dtcr_2_4[seqnames(dtcr_2_4) %in% names(exon_minus)]))
rm(exon_minus)

gc_dtcr_2_4 <- mapFromTranscripts(dtcr_2_4, exons, ignore.strand = FALSE)
########

###
exon_minus <- names(unlist(exons)[strand(unlist(exons)) == "-"])
exon_plus <- names(unlist(exons)[strand(unlist(exons)) == "+"])

strand(shapeseq_norm[seqnames(shapeseq_norm) %in% exon_minus]) <- rep("-", length(shapeseq_norm[seqnames(shapeseq_norm) %in% exon_minus]))
strand(shapeseq_norm[seqnames(shapeseq_norm) %in% exon_plus]) <- rep("+", length(shapeseq_norm[seqnames(shapeseq_norm) %in% exon_plus]))

shapeseq_norm <- shapeseq_norm[seqnames(shapeseq_norm) %in% names(exons)]
gc_dtcr_2_4 <- mapFromTranscripts(shapeseq_norm, exons, ignore.strand = FALSE)
gc_dtcr_2_4$dtcr <- elementMetadata(shapeseq_norm)$dtcr
gc_dtcr_2_4$trnames <- seqnames(shapeseq_norm)

#dtcr_256_list <- split(gc_dtcr_256, gc_dtcr_256$trnames)
save(gc_dtcr_2_4, "/Home/ii/katchyz/OUT/SHAPES_out/gc_dtcr_2_4.Rsave")

####
scds <- subsetByOverlaps(gc_sam, cds)
scds_list <- split(scds, as.vector(scds$trnames))

### get first 50nt in cds (5'end)
meta <- Reduce("+", lapply(scds_list, function(x){x[0:50]$dtcr}))

# check if sum is 0...
#lapply(scds_list, function(x){(x[0:50]$dtcr)/(sum(x[0:50]$dtcr)})

for (sn in s) {
  print(exons[[sn]])
  print("*************")
  sam <- dtcr_2_4[seqnames(dtcr_2_4) %in% c(sn)]
  strand(sam[seqnames(sam) %in% exon_plus]) <- rep("+", length(sam[seqnames(sam) %in% exon_plus]))
  strand(sam[seqnames(sam) %in% exon_minus]) <- rep("-", length(sam[seqnames(sam) %in% exon_minus]))
  sam <- sam[seqnames(sam) %in% names(exons)]
  print(sam)
  print("*************")
  gc_sam <- mapFromTranscripts(sam, exons, ignore.strand = FALSE)
  print(gc_sam)
  print("__________________*****************_________________")
}
