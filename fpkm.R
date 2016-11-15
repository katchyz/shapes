### calculating FPKM
library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)
library(plyr)

### get transcript database
txdb_can <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/danio_rerio/zebrafish_GRCz10.gtf", format = "gtf")
cds <- cdsBy(txdb_can, by="tx", use.names=TRUE)
cds <- cds[order(names(cds))]
fiveUTR <- fiveUTRsByTranscript(txdb_can, use.names=TRUE)
fiveUTR <- fiveUTR[order(names(fiveUTR))]
threeUTR <- threeUTRsByTranscript(txdb_can, use.names=TRUE)
threeUTR <- threeUTR[order(names(threeUTR))]
exons <- exonsBy(txdb_can, by="tx", use.names=TRUE)
exons <- exons[order(names(exons))]

names(cds) <- sapply(names(cds), function(x){substr(x,1,18)})
names(fiveUTR) <- sapply(names(fiveUTR), function(x){substr(x,1,18)})
names(threeUTR) <- sapply(names(threeUTR), function(x){substr(x,1,18)})
names(exons) <- sapply(names(exons), function(x){substr(x,1,18)})

txLengths <- transcriptLengths(txdb_can, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
rownames(txLengths) <- txLengths$tx_name
rownames(txLengths) <- sapply(rownames(txLengths), function(x){substr(x,1,18)})
txLengths <- txLengths[order(rownames(txLengths)),]

### RNA-Seq

RNAseq_2to4 <- readGAlignments("/Volumes/USELESS/DATA/RNA-Seq/2to4Cell.bam")
read_number_mRNA <- length(RNAseq_2to4)
cds_RNA <- countOverlaps(cds, RNAseq_2to4)
cds_RNA <- cds_RNA[order(names(cds_RNA))]
utr5_RNA <- countOverlaps(fiveUTR, RNAseq_2to4)
utr5_RNA <- utr5_RNA[order(names(utr5_RNA))]
utr3_RNA <- countOverlaps(threeUTR, RNAseq_2to4)
utr3_RNA <- utr3_RNA[order(names(utr3_RNA))]
exons_RNA <- countOverlaps(exons, RNAseq_2to4)
exons_RNA <- exons_RNA[order(names(exons_RNA))]

cds_len <- txLengths[rownames(txLengths) %in% names(cds_RNA),]
cds_len <- cds_len[order(rownames(cds_len)),]$cds_len
utr5_len <- txLengths[rownames(txLengths) %in% names(utr5_RNA),]
utr5_len <- utr5_len[order(rownames(utr5_len)),]$utr5_len
utr3_len <- txLengths[rownames(txLengths) %in% names(utr3_RNA),]
utr3_len <- utr3_len[order(rownames(utr3_len)),]$utr3_len
exons_len <- txLengths[rownames(txLengths) %in% names(exons_RNA),]
exons_len <- exons_len[order(rownames(exons_len)),]$tx_len

cds_RNA_FPKM <- (cds_RNA / cds_len) * (10^9 / read_number_mRNA)
utr5_RNA_FPKM <- (utr5_RNA / utr5_len) * (10^9 / read_number_mRNA)
utr3_RNA_FPKM <- (utr3_RNA / utr3_len) * (10^9 / read_number_mRNA)
exons_RNA_FPKM <- (exons_RNA / exons_len) * (10^9 / read_number_mRNA)

### Ribo-Seq

# read BigWigs
# substitute strand, merge
bw24fw <- import.bw("/Volumes/USELESS/OUT/zebrafish_our_out/GRCz10/2-4Cell_fw.bw")
strand(bw24fw) <- rep("+", length(bw24fw))
bw24rv <- import.bw("/Volumes/USELESS/OUT/zebrafish_our_out/GRCz10/2-4Cell_rv.bw")
strand(bw24rv) <- rep("-", length(bw24rv))
bw24 <- c(bw24fw, bw24rv)

# map to transcripts, subset CDSs
ribo_cds_24 <- mapToTranscripts(bw24, cds, ignore.strand = FALSE)
mcols(ribo_cds_24) <- cbind(mcols(ribo_cds_24), DataFrame(bw24[ribo_cds_24$xHits]))
cds_reads <- sum(score(ribo_cds_24))
ribo_cds_24 <- split(ribo_cds_24, seqnames(ribo_cds_24))
ribo_cds_24 <- ribo_cds_24[order(names(ribo_cds_24))]

# map to transcripts, subset fiveUTRs
ribo_utr5_24 <- mapToTranscripts(bw24, fiveUTR, ignore.strand = FALSE)
mcols(ribo_utr5_24) <- cbind(mcols(ribo_utr5_24), DataFrame(bw24[ribo_utr5_24$xHits]))
utr5_reads <- sum(score(ribo_utr5_24))
ribo_utr5_24 <- split(ribo_utr5_24, seqnames(ribo_utr5_24))
ribo_utr5_24 <- ribo_utr5_24[order(names(ribo_utr5_24))]

# map to transcripts, subset threeUTRs
ribo_utr3_24 <- mapToTranscripts(bw24, threeUTR, ignore.strand = FALSE)
mcols(ribo_utr3_24) <- cbind(mcols(ribo_utr3_24), DataFrame(bw24[ribo_utr3_24$xHits]))
utr3_reads <- sum(score(ribo_utr3_24))
ribo_utr3_24 <- split(ribo_utr3_24, seqnames(ribo_utr3_24))
ribo_utr3_24 <- ribo_utr3_24[order(names(ribo_utr3_24))]

# map to transcripts, subset threeUTRs
ribo_exons_24 <- mapToTranscripts(bw24, exons, ignore.strand = FALSE)
mcols(ribo_exons_24) <- cbind(mcols(ribo_exons_24), DataFrame(bw24[ribo_exons_24$xHits]))
exons_reads <- sum(score(ribo_exons_24))
ribo_exons_24 <- split(ribo_exons_24, seqnames(ribo_exons_24))
ribo_exons_24 <- ribo_exons_24[order(names(ribo_exons_24))]

###
#footprint_number_total <- cds_reads + utr3_reads + utr5_reads
footprint_number_total <- exons_reads

ribo_cds <- sapply(ribo_cds_24, function(x){sum(x$score)})
ribo_utr5 <- sapply(ribo_utr5_24, function(x){sum(x$score)})
ribo_utr3 <- sapply(ribo_utr3_24, function(x){sum(x$score)})
ribo_exons <- sapply(ribo_exons_24, function(x){sum(x$score)})

cds_len <- txLengths[rownames(txLengths) %in% names(ribo_cds),]$cds_len
utr5_len <- txLengths[rownames(txLengths) %in% names(ribo_utr5),]$utr5_len
utr3_len <- txLengths[rownames(txLengths) %in% names(ribo_utr3),]$utr3_len
exons_len <- txLengths[rownames(txLengths) %in% names(ribo_exons),]$tx_len

cds_Ribo_FPKM <- (ribo_cds / cds_len) * (10^9 / footprint_number_total)
utr5_Ribo_FPKM <- (ribo_utr5 / utr5_len) * (10^9 / footprint_number_total)
utr3_Ribo_FPKM <- (ribo_utr3 / utr3_len) * (10^9 / footprint_number_total)
exons_Ribo_FPKM <- (ribo_exons / exons_len) * (10^9 / footprint_number_total)

### PUT IN ONE DATA FRAME
#x <- c(ens1 = 1, ens2 = 2)
#y <- c(ens2 = 3, ens3 = 4, ens4 = 5)
#z <- c(ens1 = 6, ens3 = 7, ens5 = 8)
#t(do.call(rbind.fill, list(data.frame(t(x)),data.frame(t(y)), data.frame(t(z)))))

FPKM_24 <- data.frame(t(do.call(rbind.fill, list(data.frame(t(utr5_RNA_FPKM)),data.frame(t(cds_RNA_FPKM)), data.frame(t(utr3_RNA_FPKM)),data.frame(t(exons_RNA_FPKM)), data.frame(t(utr5_Ribo_FPKM)),data.frame(t(cds_Ribo_FPKM)), data.frame(t(utr3_Ribo_FPKM)),data.frame(t(exons_Ribo_FPKM))))))
colnames(FPKM_24) <- c("utr5_RNA_fpkm", "cds_RNA_fpkm", "utr3_RNA_fpkm", "exons_RNA_fpkm", "utr5_Ribo_fpkm", "cds_Ribo_fpkm", "utr3_Ribo_fpkm", "exons_Ribo_fpkm")

save(FPKM_24, file = "/Volumes/USELESS/META/SHAPES/FPKM_24.Rdata")


######
######
#####
#####
#####
###
###
##
##
#

### RNA-Seq

RNAseq_256 <- readGAlignments("/Volumes/USELESS/DATA/RNA-Seq/256Cell.bam")
read_number_mRNA <- length(RNAseq_256)
cds_RNA <- countOverlaps(cds, RNAseq_256)
cds_RNA <- cds_RNA[order(names(cds_RNA))]
utr5_RNA <- countOverlaps(fiveUTR, RNAseq_256)
utr5_RNA <- utr5_RNA[order(names(utr5_RNA))]
utr3_RNA <- countOverlaps(threeUTR, RNAseq_256)
utr3_RNA <- utr3_RNA[order(names(utr3_RNA))]
exons_RNA <- countOverlaps(exons, RNAseq_256)
utr3_RNA <- utr3_RNA[order(names(utr3_RNA))]

cds_len <- txLengths[rownames(txLengths) %in% names(cds_RNA),]
cds_len <- cds_len[order(rownames(cds_len)),]$cds_len
utr5_len <- txLengths[rownames(txLengths) %in% names(utr5_RNA),]
utr5_len <- utr5_len[order(rownames(utr5_len)),]$utr5_len
utr3_len <- txLengths[rownames(txLengths) %in% names(utr3_RNA),]
utr3_len <- utr3_len[order(rownames(utr3_len)),]$utr3_len
exons_len <- txLengths[rownames(txLengths) %in% names(exons_RNA),]
exons_len <- exons_len[order(rownames(exons_len)),]$tx_len

cds_RNA_FPKM <- (cds_RNA / cds_len) * (10^9 / read_number_mRNA)
utr5_RNA_FPKM <- (utr5_RNA / utr5_len) * (10^9 / read_number_mRNA)
utr3_RNA_FPKM <- (utr3_RNA / utr3_len) * (10^9 / read_number_mRNA)
exons_RNA_FPKM <- (exons_RNA / exons_len) * (10^9 / read_number_mRNA)

### Ribo-Seq

# read BigWigs
# substitute strand, merge
bw256fw <- import.bw("/Volumes/USELESS/OUT/zebrafish_our_out/GRCz10/256Cell_fw.bw")
strand(bw256fw) <- rep("+", length(bw256fw))
bw256rv <- import.bw("/Volumes/USELESS/OUT/zebrafish_our_out/GRCz10/256Cell_rv.bw")
strand(bw256rv) <- rep("-", length(bw256rv))
bw256 <- c(bw256fw, bw256rv)

# map to transcripts, subset CDSs
ribo_cds_256 <- mapToTranscripts(bw256, cds, ignore.strand = FALSE)
mcols(ribo_cds_256) <- cbind(mcols(ribo_cds_256), DataFrame(bw256[ribo_cds_256$xHits]))
cds_reads <- sum(score(ribo_cds_256))
ribo_cds_256 <- split(ribo_cds_256, seqnames(ribo_cds_256))
ribo_cds_256 <- ribo_cds_256[order(names(ribo_cds_256))]

# map to transcripts, subset fiveUTRs
ribo_utr5_256 <- mapToTranscripts(bw256, fiveUTR, ignore.strand = FALSE)
mcols(ribo_utr5_256) <- cbind(mcols(ribo_utr5_256), DataFrame(bw256[ribo_utr5_256$xHits]))
utr5_reads <- sum(score(ribo_utr5_256))
ribo_utr5_256 <- split(ribo_utr5_256, seqnames(ribo_utr5_256))
ribo_utr5_256 <- ribo_utr5_256[order(names(ribo_utr5_256))]

# map to transcripts, subset threeUTRs
ribo_utr3_256 <- mapToTranscripts(bw256, threeUTR, ignore.strand = FALSE)
mcols(ribo_utr3_256) <- cbind(mcols(ribo_utr3_256), DataFrame(bw256[ribo_utr3_256$xHits]))
utr3_reads <- sum(score(ribo_utr3_256))
ribo_utr3_256 <- split(ribo_utr3_256, seqnames(ribo_utr3_256))
ribo_utr3_256 <- ribo_utr3_256[order(names(ribo_utr3_256))]

# map to transcripts, subset exons
ribo_exons_256 <- mapToTranscripts(bw256, exons, ignore.strand = FALSE)
mcols(ribo_exons_256) <- cbind(mcols(ribo_exons_256), DataFrame(bw256[ribo_exons_256$xHits]))
exons_reads <- sum(score(ribo_exons_256))
ribo_exons_256 <- split(ribo_exons_256, seqnames(ribo_exons_256))
ribo_exons_256 <- ribo_exons_256[order(names(ribo_exons_256))]

###
#footprint_number_total <- cds_reads + utr3_reads + utr5_reads
footprint_number_total <- exons_reads

ribo_cds <- sapply(ribo_cds_256, function(x){sum(x$score)})
ribo_utr5 <- sapply(ribo_utr5_256, function(x){sum(x$score)})
ribo_utr3 <- sapply(ribo_utr3_256, function(x){sum(x$score)})
ribo_exons <- sapply(ribo_exons_256, function(x){sum(x$score)})

cds_len <- txLengths[rownames(txLengths) %in% names(ribo_cds),]$cds_len
utr5_len <- txLengths[rownames(txLengths) %in% names(ribo_utr5),]$utr5_len
utr3_len <- txLengths[rownames(txLengths) %in% names(ribo_utr3),]$utr3_len
exons_len <- txLengths[rownames(txLengths) %in% names(ribo_exons),]$tx_len

cds_Ribo_FPKM <- (ribo_cds / cds_len) * (10^9 / footprint_number_total)
utr5_Ribo_FPKM <- (ribo_utr5 / utr5_len) * (10^9 / footprint_number_total)
utr3_Ribo_FPKM <- (ribo_utr3 / utr3_len) * (10^9 / footprint_number_total)
exons_Ribo_FPKM <- (ribo_exons / exons_len) * (10^9 / footprint_number_total)

### PUT IN ONE DATA FRAME
x <- c(ens1 = 1, ens2 = 2)
y <- c(ens2 = 3, ens3 = 4, ens4 = 5)
z <- c(ens1 = 6, ens3 = 7, ens5 = 8)
t(do.call(rbind.fill, list(data.frame(t(x)),data.frame(t(y)), data.frame(t(z)))))

FPKM_256 <- data.frame(t(do.call(rbind.fill, list(data.frame(t(utr5_RNA_FPKM)),data.frame(t(cds_RNA_FPKM)), data.frame(t(utr3_RNA_FPKM)),data.frame(t(exons_RNA_FPKM)), data.frame(t(utr5_Ribo_FPKM)),data.frame(t(cds_Ribo_FPKM)), data.frame(t(utr3_Ribo_FPKM)),data.frame(t(exons_Ribo_FPKM))))))
colnames(FPKM_256) <- c("utr5_RNA_fpkm", "cds_RNA_fpkm", "utr3_RNA_fpkm", "exons_RNA_fpkm", "utr5_Ribo_fpkm", "cds_Ribo_fpkm", "utr3_Ribo_fpkm", "exons_Ribo_fpkm")

save(FPKM_256, file = "/Volumes/USELESS/META/SHAPES/FPKM_256.Rdata")

