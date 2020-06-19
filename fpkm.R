### calculating FPKM
library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)
library(plyr)

### get transcript database
txdb_can <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/danio_rerio/zebrafish_GRCz10.gtf", format = "gtf")
txdb_can <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/danio_rerio/Danio_rerio.GRCz10.84.gtf", format = "gtf")
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


###################
### 1K cell #######

### RNA-Seq

RNAseq_1K <- readGAlignments("/Volumes/USELESS/DATA/RNA-Seq/1Kcell_mRNA.bam")
read_number_mRNA <- length(RNAseq_1K)
cds_RNA <- countOverlaps(cds, RNAseq_1K)
cds_RNA <- cds_RNA[order(names(cds_RNA))]
utr5_RNA <- countOverlaps(fiveUTR, RNAseq_1K)
utr5_RNA <- utr5_RNA[order(names(utr5_RNA))]
utr3_RNA <- countOverlaps(threeUTR, RNAseq_1K)
utr3_RNA <- utr3_RNA[order(names(utr3_RNA))]
exons_RNA <- countOverlaps(exons, RNAseq_1K)
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
bw1Kfw <- import.bw("/Volumes/USELESS/OUT/zebrafish_our_out/GRCz10/1KCell_fw.bw")
strand(bw1Kfw) <- rep("+", length(bw1Kfw))
bw1Krv <- import.bw("/Volumes/USELESS/OUT/zebrafish_our_out/GRCz10/1KCell_rv.bw")
strand(bw1Krv) <- rep("-", length(bw1Krv))
bw1K <- c(bw1Kfw, bw1Krv)

# map to transcripts, subset CDSs
ribo_cds_1K <- mapToTranscripts(bw1K, cds, ignore.strand = FALSE)
mcols(ribo_cds_1K) <- cbind(mcols(ribo_cds_1K), DataFrame(bw1K[ribo_cds_1K$xHits]))
cds_reads <- sum(score(ribo_cds_1K))
ribo_cds_1K <- split(ribo_cds_1K, seqnames(ribo_cds_1K))
ribo_cds_1K <- ribo_cds_1K[order(names(ribo_cds_1K))]

# map to transcripts, subset fiveUTRs
ribo_utr5_1K <- mapToTranscripts(bw1K, fiveUTR, ignore.strand = FALSE)
mcols(ribo_utr5_1K) <- cbind(mcols(ribo_utr5_1K), DataFrame(bw1K[ribo_utr5_1K$xHits]))
utr5_reads <- sum(score(ribo_utr5_K))
ribo_utr5_1K <- split(ribo_utr5_1K, seqnames(ribo_utr5_1K))
ribo_utr5_1K <- ribo_utr5_1K[order(names(ribo_utr5_1K))]

# map to transcripts, subset threeUTRs
ribo_utr3_1K <- mapToTranscripts(bw1K, threeUTR, ignore.strand = FALSE)
mcols(ribo_utr3_1K) <- cbind(mcols(ribo_utr3_1K), DataFrame(bw1K[ribo_utr3_1K$xHits]))
utr3_reads <- sum(score(ribo_utr3_1K))
ribo_utr3_1K <- split(ribo_utr3_1K, seqnames(ribo_utr3_1K))
ribo_utr3_1K <- ribo_utr3_1K[order(names(ribo_utr3_1K))]

# map to transcripts, subset exons
ribo_exons_1K <- mapToTranscripts(bw1K, exons, ignore.strand = FALSE)
mcols(ribo_exons_1K) <- cbind(mcols(ribo_exons_1K), DataFrame(bw1K[ribo_exons_1K$xHits]))
exons_reads <- sum(score(ribo_exons_1K))
ribo_exons_1K <- split(ribo_exons_1K, seqnames(ribo_exons_1K))
ribo_exons_1K <- ribo_exons_1K[order(names(ribo_exons_1K))]

###
#footprint_number_total <- cds_reads + utr3_reads + utr5_reads
footprint_number_total <- exons_reads

ribo_cds <- sapply(ribo_cds_1K, function(x){sum(x$score)})
ribo_utr5 <- sapply(ribo_utr5_1K, function(x){sum(x$score)})
ribo_utr3 <- sapply(ribo_utr3_1K, function(x){sum(x$score)})
ribo_exons <- sapply(ribo_exons_1K, function(x){sum(x$score)})

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

FPKM_1K <- data.frame(t(do.call(rbind.fill, list(data.frame(t(utr5_RNA_FPKM)),data.frame(t(cds_RNA_FPKM)), data.frame(t(utr3_RNA_FPKM)),data.frame(t(exons_RNA_FPKM)), data.frame(t(utr5_Ribo_FPKM)),data.frame(t(cds_Ribo_FPKM)), data.frame(t(utr3_Ribo_FPKM)),data.frame(t(exons_Ribo_FPKM))))))
colnames(FPKM_1K) <- c("utr5_RNA_fpkm", "cds_RNA_fpkm", "utr3_RNA_fpkm", "exons_RNA_fpkm", "utr5_Ribo_fpkm", "cds_Ribo_fpkm", "utr3_Ribo_fpkm", "exons_Ribo_fpkm")

save(FPKM_1K, file = "/Volumes/USELESS/META/SHAPES/FPKM_1K.Rdata")


############
###################
### 1 cell #######

### RNA-Seq

RNAseq_1cell <- readGAlignments("/Volumes/USELESS/DATA/RNA-Seq/1cell_mRNA.bam")
read_number_mRNA <- length(RNAseq_1cell)
cds_RNA <- countOverlaps(cds, RNAseq_1cell)
cds_RNA <- cds_RNA[order(names(cds_RNA))]
utr5_RNA <- countOverlaps(fiveUTR, RNAseq_1cell)
utr5_RNA <- utr5_RNA[order(names(utr5_RNA))]
utr3_RNA <- countOverlaps(threeUTR, RNAseq_1cell)
utr3_RNA <- utr3_RNA[order(names(utr3_RNA))]
exons_RNA <- countOverlaps(exons, RNAseq_1cell)
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

FPKM_1cell <- data.frame(t(do.call(rbind.fill, list(data.frame(t(utr5_RNA_FPKM)),data.frame(t(cds_RNA_FPKM)), data.frame(t(utr3_RNA_FPKM)),data.frame(t(exons_RNA_FPKM))))))
colnames(FPKM_1cell) <- c("utr5_RNA_fpkm", "cds_RNA_fpkm", "utr3_RNA_fpkm", "exons_RNA_fpkm")
FPKM_1cell$n <- rownames(FPKM_1cell)

save(FPKM_1cell, file = "/Volumes/USELESS/META/SHAPES/FPKM_1cell.Rdata")

cell1 <- data.frame(n = FPKM_1cell$n, rna_cell1 = FPKM_1cell$exons_RNA_fpkm)
df <- merge(df, cell1, by="n", all=TRUE)



#######
### RNA-Seq

RNAseq_1cell <- readGAlignments("/Volumes/USELESS/DATA/RNA-Seq/cell1_aanes.bam")
read_number_mRNA <- length(RNAseq_1cell)
cds_RNA <- countOverlaps(cds, RNAseq_1cell)
cds_RNA <- cds_RNA[order(names(cds_RNA))]
utr5_RNA <- countOverlaps(fiveUTR, RNAseq_1cell)
utr5_RNA <- utr5_RNA[order(names(utr5_RNA))]
utr3_RNA <- countOverlaps(threeUTR, RNAseq_1cell)
utr3_RNA <- utr3_RNA[order(names(utr3_RNA))]
exons_RNA <- countOverlaps(exons, RNAseq_1cell)
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

FPKM_1cell_aanes <- data.frame(t(do.call(rbind.fill, list(data.frame(t(utr5_RNA_FPKM)),data.frame(t(cds_RNA_FPKM)), data.frame(t(utr3_RNA_FPKM)),data.frame(t(exons_RNA_FPKM))))))
colnames(FPKM_1cell_aanes) <- c("utr5_RNA_fpkm", "cds_RNA_fpkm", "utr3_RNA_fpkm", "exons_RNA_fpkm")
FPKM_1cell_aanes$n <- rownames(FPKM_1cell_aanes)

save(FPKM_1cell_aanes, file = "/Volumes/USELESS/META/SHAPES/FPKM_1cell_aanes.Rdata")

cell1 <- data.frame(n = FPKM_1cell_aanes$n, rna_cell1 = FPKM_1cell_aanes$exons_RNA_fpkm)
df <- merge(df, cell1, by="n", all=TRUE)

##
RNAseq_3_5h <- readGAlignments("/Volumes/USELESS/DATA/RNA-Seq/aanes_3_5h.bam")
read_number_mRNA <- length(RNAseq_3_5h)
cds_RNA <- countOverlaps(cds, RNAseq_3_5h)
cds_RNA <- cds_RNA[order(names(cds_RNA))]
utr5_RNA <- countOverlaps(fiveUTR, RNAseq_3_5h)
utr5_RNA <- utr5_RNA[order(names(utr5_RNA))]
utr3_RNA <- countOverlaps(threeUTR, RNAseq_3_5h)
utr3_RNA <- utr3_RNA[order(names(utr3_RNA))]
exons_RNA <- countOverlaps(exons, RNAseq_3_5h)
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

FPKM_3_5h_aanes <- data.frame(t(do.call(rbind.fill, list(data.frame(t(utr5_RNA_FPKM)),data.frame(t(cds_RNA_FPKM)), data.frame(t(utr3_RNA_FPKM)),data.frame(t(exons_RNA_FPKM))))))
colnames(FPKM_3_5h_aanes) <- c("utr5_RNA_fpkm", "cds_RNA_fpkm", "utr3_RNA_fpkm", "exons_RNA_fpkm")
FPKM_3_5h_aanes$n <- rownames(FPKM_3_5h_aanes)

save(FPKM_3_5h_aanes, file = "/Volumes/USELESS/META/SHAPES/FPKM_3_5h_aanes.Rdata")

h35 <- data.frame(n = FPKM_3_5h_aanes$n, rna_3_5h = FPKM_3_5h_aanes$exons_RNA_fpkm)
df <- merge(df, h35, by="n", all=TRUE)


########### cell1 total RNA (Vesterlund)
RNAseq_1cell <- readGAlignments("/Volumes/USELESS/DATA/RNA-Seq/cell1_vesterlund.bam")
read_number_mRNA <- length(RNAseq_1cell)
cds_RNA <- countOverlaps(cds, RNAseq_1cell)
cds_RNA <- cds_RNA[order(names(cds_RNA))]
utr5_RNA <- countOverlaps(fiveUTR, RNAseq_1cell)
utr5_RNA <- utr5_RNA[order(names(utr5_RNA))]
utr3_RNA <- countOverlaps(threeUTR, RNAseq_1cell)
utr3_RNA <- utr3_RNA[order(names(utr3_RNA))]
exons_RNA <- countOverlaps(exons, RNAseq_1cell)
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

FPKM_1cell_vesterlund <- data.frame(t(do.call(rbind.fill, list(data.frame(t(utr5_RNA_FPKM)),data.frame(t(cds_RNA_FPKM)), data.frame(t(utr3_RNA_FPKM)),data.frame(t(exons_RNA_FPKM))))))
colnames(FPKM_1cell_vesterlund) <- c("utr5_RNA_fpkm", "cds_RNA_fpkm", "utr3_RNA_fpkm", "exons_RNA_fpkm")
FPKM_1cell_vesterlund$n <- rownames(FPKM_1cell_vesterlund)

save(FPKM_1cell_vesterlund, file = "/Volumes/USELESS/META/SHAPES/FPKM_1cell_vesterlund.Rdata")

cell1 <- data.frame(n = FPKM_1cell_vesterlund$n, rna_cell1_total = FPKM_1cell_vesterlund$exons_RNA_fpkm)
df <- merge(df, cell1, by="n", all=TRUE)


## bazzini 4h
txdb_can <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/Danio_rerio.Zv9.79.gtf", format = "gtf")
exons <- exonsBy(txdb_can, by="tx", use.names=TRUE)
txLengths <- transcriptLengths(txdb_can, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
rownames(txLengths) <- txLengths$tx_name
txLengths <- txLengths[order(rownames(txLengths)),]

# read BigWigs
# substitute strand, merge
bw1Kfw <- import.bw("/Volumes/USELESS/DATA/Ribo-Seq/bazzini_4h_fwd.bw")
strand(bw1Kfw) <- rep("+", length(bw1Kfw))
bw1Krv <- import.bw("/Volumes/USELESS/DATA/Ribo-Seq/bazzini_4h_fwd.bw")
strand(bw1Krv) <- rep("-", length(bw1Krv))
bw1K <- c(bw1Kfw, bw1Krv)

# map to transcripts, subset exons
ribo_exons_1K <- mapToTranscripts(bw1K, exons, ignore.strand = FALSE)
mcols(ribo_exons_1K) <- cbind(mcols(ribo_exons_1K), DataFrame(bw1K[ribo_exons_1K$xHits]))
exons_reads <- sum(score(ribo_exons_1K))
ribo_exons_1K <- split(ribo_exons_1K, seqnames(ribo_exons_1K))
ribo_exons_1K <- ribo_exons_1K[order(names(ribo_exons_1K))]

###
#footprint_number_total <- cds_reads + utr3_reads + utr5_reads
footprint_number_total <- exons_reads
ribo_exons <- sapply(ribo_exons_1K, function(x){sum(x$score)})
exons_len <- txLengths[rownames(txLengths) %in% names(ribo_exons),]$tx_len
exons_Ribo_FPKM <- (ribo_exons / exons_len) * (10^9 / footprint_number_total)

FPKM_4h <- data.frame(n = names(exons_Ribo_FPKM), ribo_4h = as.numeric(exons_Ribo_FPKM))

save(FPKM_4h, file = "/Volumes/USELESS/META/SHAPES_OLD/FPKM_4h_ribo.Rdata")


