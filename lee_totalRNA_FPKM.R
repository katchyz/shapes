### calculating FPKM
### Lee 2013, total RNA, 2h & 4h, WT
library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)
library(plyr)

### get transcript database
txdb_can <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/Danio_rerio.Zv9.79.gtf", format = "gtf")
exons <- exonsBy(txdb_can, by="tx", use.names=TRUE)
exons <- exons[order(names(exons))]

txLengths <- transcriptLengths(txdb_can, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
rownames(txLengths) <- txLengths$tx_name
txLengths <- txLengths[order(rownames(txLengths)),]

### RNA-Seq

RNAseq <- readGAlignments("/Volumes/USELESS/DATA/RNA-Seq/lee2013_totalRNA_2h.bam")
read_number_mRNA <- length(RNAseq)
exons_RNA <- countOverlaps(exons, RNAseq)
exons_RNA <- exons_RNA[order(names(exons_RNA))]

exons_len <- txLengths[rownames(txLengths) %in% names(exons_RNA),]
exons_len <- exons_len[order(rownames(exons_len)),]$tx_len

lee_RNA_FPKM <- (exons_RNA / exons_len) * (10^9 / read_number_mRNA)

save(lee_RNA_FPKM, file="/Volumes/USELESS/META/SHAPES/FPKM_Lee.Rsave")
write.csv(lee_RNA_FPKM, file="/Volumes/USELESS/META/SHAPES/FPKM_Lee.csv")

# calculate on our data
load(file = "/Volumes/USELESS/META/SHAPES/FPKM_256.Rdata")

o <- FPKM_256[with(FPKM_256, order(-exons_RNA_fpkm)), ]
het <- rownames(o[1:4000,])
het <- sort(het)
het <- het[het %in% txlen$tx_name]

fileConn<-file("/Volumes/USELESS/META/SHAPES/avoidance/het.txt")
writeLines(het, fileConn)
close(fileConn)


### OTPG

# get longest transcript per gene
txlen <- arrange(txLengths, gene_id, desc(tx_len))
txlen <- txlen[!duplicated(txlen$gene_id),]
# get transcripts with at least 33 on CDS and 30 on utr3
otpg <- sort(txlen[(txlen$cds_len > 33) & (txlen$utr3_len > 30), ]$tx_name)

