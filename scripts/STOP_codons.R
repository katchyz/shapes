##### STOP codons

## split by stop codon (UAA/UAG), bar plots

## shape vs CAI of the last 50 codons on CDS
## shape vs CAI of the last 50 codons on 3'UTR

## avg shape in windows (e.g. first/last 50) vs TE

#############
library(GenomicFeatures)
library(ggplot2)
library(GeneCycle)
library(reshape2)
library(ggplot2)

txdb_can <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/danio_rerio/zebrafish_GRCz10.gtf", format = "gtf")
cds <- cdsBy(txdb_can, by="tx", use.names=TRUE)
fiveUTR <- fiveUTRsByTranscript(txdb_can, use.names=TRUE)
threeUTR <- threeUTRsByTranscript(txdb_can, use.names=TRUE)

names(cds) <- sapply(names(cds), function(x){substr(x,1,18)})
names(fiveUTR) <- sapply(names(fiveUTR), function(x){substr(x,1,18)})
names(threeUTR) <- sapply(names(fiveUTR), function(x){substr(x,1,18)})

cds_2_4 <- subsetByOverlaps(gc_unsmoothed_2_4, cds)
cds_2_4 <- split(cds_2_4, cds_2_4$trnames)
cds_2_4 <- cds_2_4[lengths(cds_2_4) > 50]
names(cds_2_4) <- sapply(names(cds_2_4), function(x){substr(x,1,18)})

cds_256 <- subsetByOverlaps(gc_unsmoothed_256, cds)
cds_256 <- split(cds_256, cds_256$trnames)
cds_256 <- cds_256[lengths(cds_256) > 50]
names(cds_256) <- sapply(names(cds_256), function(x){substr(x,1,18)})

five_2_4 <- subsetByOverlaps(gc_unsmoothed_2_4, fiveUTR)
five_2_4 <- split(five_2_4, five_2_4$trnames)
five_2_4 <- five_2_4[lengths(five_2_4) > 50]

three_2_4 <- subsetByOverlaps(gc_unsmoothed_2_4, threeUTR)
three_2_4 <- split(three_2_4, three_2_4$trnames)
three_2_4 <- three_2_4[lengths(three_2_4) > 50]
names(three_2_4) <- sapply(names(three_2_4), function(x){substr(x,1,18)})

five_256 <- subsetByOverlaps(gc_unsmoothed_256, fiveUTR)
five_256 <- split(five_256, five_256$trnames)
five_256 <- five_256[lengths(five_256) > 50]

three_256 <- subsetByOverlaps(gc_unsmoothed_256, threeUTR)
three_256 <- split(three_256, three_256$trnames)
three_256 <- three_256[lengths(three_256) > 50]
names(three_256) <- sapply(names(three_256), function(x){substr(x,1,18)})

cds_2_4 <- cds_2_4[names(cds_2_4) %in% names(three_2_4)]
three_2_4 <- three_2_4[names(three_2_4) %in% names(cds_2_4)]

cds_256 <- cds_256[names(cds_256) %in% names(three_256)]
three_256 <- three_256[names(three_256) %in% names(cds_256)]

#### !!!!!!!!!!!!!!!
#### merge the two granges (-32, 30)

cds_2_4_last33 <- lapply(cds_2_4, function(x){x[(length(x)-32):length(x)]$dtcr})
three_2_4_fi30 <- lapply(three_2_4, function(x){x[1:30]$dtcr})

cds_three_2_4 <- t(mapply(c, cds_2_4_last33, three_2_4_fi30))

cds_256_last33 <- lapply(cds_256, function(x){x[(length(x)-32):length(x)]$dtcr})
three_256_fi30 <- lapply(three_256, function(x){x[1:30]$dtcr})

cds_three_256 <- t(mapply(c, cds_256_last33, three_256_fi30))

### normalize by sum of both CDS and 3'UTR fragments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

cds_three_2_4 <- cds_three_2_4[rowSums(cds_three_2_4) > 0,]
cds_three_256 <- cds_three_256[rowSums(cds_three_256) > 0,]

meta_2_4 <- colSums(cds_three_2_4 / rowSums(cds_three_2_4))
meta_256 <- colSums(cds_three_256 / rowSums(cds_three_256))
###### stop ######
scale <- c(-33:30)[-34]
df <- data.frame(scale,meta_2_4)
p_2_4 <- ggplot(df, aes(x=scale, y=meta_2_4)) + geom_bar(stat = "identity") + xlab("position from stop codon") + ylab("sum(reactivities)") + ggtitle("2-4cell")
ggsave(file = "/Volumes/USELESS/META/SHAPES/stop/p_2_4.png", plot = p_2_4)

##
df <- data.frame(scale,meta_256)
p_256 <- ggplot(df, aes(x=scale, y=meta_256)) + geom_bar(stat = "identity") + xlab("position from stop codon") + ylab("sum(reactivities)") + ggtitle("256cell")
ggsave(file = "/Volumes/USELESS/META/SHAPES/stop/p_256.png", plot = p_256)

##############
##############
library(seqinr)

fasta_cds <- read.fasta("/Volumes/USELESS/DATA/fasta/danio_rerio/GRCz10/cds/Danio_rerio.GRCz10.cds.all.fa")
names(fasta_cds) <- sapply(names(fasta_cds), function(x){substr(x,1,18)})

fasta_cds <- fasta_cds[unlist(lapply(fasta_cds, function(x){length(x) %% 3 == 0}))]

taa <- names(fasta_cds[unlist(lapply(fasta_cds, function(x){paste(x[(length(x)-2):length(x)], collapse = "") == "taa"}))])
tga <- names(fasta_cds[unlist(lapply(fasta_cds, function(x){paste(x[(length(x)-2):length(x)], collapse = "") == "tga"}))])

tag <- names(fasta_cds[unlist(lapply(fasta_cds, function(x){paste(x[(length(x)-2):length(x)], collapse = "") == "tag"}))])

## 2-4cell
cds_three_2_4_taa <- cds_three_2_4[rownames(cds_three_2_4) %in% taa,]
meta_2_4_taa <- colSums(cds_three_2_4_taa / rowSums(cds_three_2_4_taa))
df <- data.frame(scale,meta_2_4_taa)
p_2_4_taa <- ggplot(df, aes(x=scale, y=meta_2_4_taa)) + geom_bar(stat = "identity") + xlab("position from stop codon") + ylab("sum(reactivities)") + ggtitle("2-4cell, TAA")
ggsave(file = "/Volumes/USELESS/META/SHAPES/stop/p_2_4_taa.png", plot = p_2_4_taa)

cds_three_2_4_tga <- cds_three_2_4[rownames(cds_three_2_4) %in% tga,]
meta_2_4_tga <- colSums(cds_three_2_4_tga / rowSums(cds_three_2_4_tga))
df <- data.frame(scale,meta_2_4_tga)
p_2_4_tga <- ggplot(df, aes(x=scale, y=meta_2_4_tga)) + geom_bar(stat = "identity") + xlab("position from stop codon") + ylab("sum(reactivities)") + ggtitle("2-4cell, TGA")
ggsave(file = "/Volumes/USELESS/META/SHAPES/stop/p_2_4_tga.png", plot = p_2_4_tga)

cds_three_2_4_tag <- cds_three_2_4[rownames(cds_three_2_4) %in% tag,]
meta_2_4_tag <- colSums(cds_three_2_4_tag / rowSums(cds_three_2_4_tag))
df <- data.frame(scale,meta_2_4_tag)
p_2_4_tag <- ggplot(df, aes(x=scale, y=meta_2_4_tag)) + geom_bar(stat = "identity") + xlab("position from stop codon") + ylab("sum(reactivities)") + ggtitle("2-4cell, TAG")
ggsave(file = "/Volumes/USELESS/META/SHAPES/stop/p_2_4_tag.png", plot = p_2_4_tag)


## 256cell
cds_three_256_taa <- cds_three_256[rownames(cds_three_256) %in% taa,]
meta_256_taa <- colSums(cds_three_256_taa / rowSums(cds_three_256_taa))
df <- data.frame(scale,meta_256_taa)
p_256_taa <- ggplot(df, aes(x=scale, y=meta_256_taa)) + geom_bar(stat = "identity") + xlab("position from stop codon") + ylab("sum(reactivities)") + ggtitle("256cell, TAA")
ggsave(file = "/Volumes/USELESS/META/SHAPES/stop/p_256_taa.png", plot = p_256_taa)

cds_three_256_tga <- cds_three_256[rownames(cds_three_256) %in% tga,]
meta_256_tga <- colSums(cds_three_256_tga / rowSums(cds_three_256_tga))
df <- data.frame(scale,meta_256_tga)
p_256_tga <- ggplot(df, aes(x=scale, y=meta_256_tga)) + geom_bar(stat = "identity") + xlab("position from stop codon") + ylab("sum(reactivities)") + ggtitle("256cell, TGA")
ggsave(file = "/Volumes/USELESS/META/SHAPES/stop/p_256_tga.png", plot = p_256_tga)

cds_three_256_tag <- cds_three_256[rownames(cds_three_256) %in% tag,]
meta_256_tag <- colSums(cds_three_256_tag / rowSums(cds_three_256_tag))
df <- data.frame(scale,meta_256_tag)
p_256_tag <- ggplot(df, aes(x=scale, y=meta_256_tag)) + geom_bar(stat = "identity") + xlab("position from stop codon") + ylab("sum(reactivities)") + ggtitle("256cell, TAG")
ggsave(file = "/Volumes/USELESS/META/SHAPES/stop/p_256_tag.png", plot = p_256_tag)

#####
# Ribo-seq in 3'UTRs for TAA, TGA, TAG

# get longest transcript per gene
library(plyr)
txlen <- arrange(txLengths, gene_id, desc(tx_len))
txlen <- txlen[!duplicated(txlen$gene_id),]
# get transcripts with at least 33 on CDS and 30 on utr3
otpg <- sort(txlen[(txlen$cds_len > 33) & (txlen$utr3_len > 30), ]$tx_name)

# cov_cds_24, cov_utr3_24 --- from whole_tx_ribo.R
# names: taa, tga, (tag)

taa_common <- Reduce(intersect, list(taa, names(cov_cds_24), names(cov_utr3_24), otpg))
tga_common <- Reduce(intersect, list(tga, names(cov_cds_24), names(cov_utr3_24), otpg))
tag_common <- Reduce(intersect, list(tag, names(cov_cds_24), names(cov_utr3_24), otpg))

# sample the same number of tx ?
scale <- c(-33:30)[-34]

# 2-4cell
# TAA
meta <- t(rbind(sapply(cov_cds_24[names(cov_cds_24) %in% taa_common], function(x){x[(length(x)-32):length(x)]}), sapply(cov_utr3_24[names(cov_utr3_24) %in% taa_common], function(x){x[1:30]})))
meta[is.na(meta)] <- 0
# normalize
meta <- meta / rowSums(meta)
meta[is.nan(meta)] <- 0
meta <- colSums(meta)
df <- data.frame(scale,meta)
p_24taa <- ggplot(df, aes(x=scale, y=meta)) + geom_bar(stat = "identity") + xlab("position from stop codon") + ylab("sum(Ribo-seq)") + ggtitle("2-4cell, TAA")
ggsave(file = "/Volumes/USELESS/META/SHAPES/stop/ribo_taa_24.png", plot = p_24taa)

# TGA
meta <- t(rbind(sapply(cov_cds_24[names(cov_cds_24) %in% tga_common], function(x){x[(length(x)-32):length(x)]}), sapply(cov_utr3_24[names(cov_utr3_24) %in% tga_common], function(x){x[1:30]})))
meta[is.na(meta)] <- 0
# normalize
meta <- meta / rowSums(meta)
meta[is.nan(meta)] <- 0
meta <- colSums(meta)
df <- data.frame(scale,meta)
p_24tga <- ggplot(df, aes(x=scale, y=meta)) + geom_bar(stat = "identity") + xlab("position from stop codon") + ylab("sum(Ribo-seq)") + ggtitle("2-4cell, TGA")
ggsave(file = "/Volumes/USELESS/META/SHAPES/stop/ribo_tga_24.png", plot = p_24tga)

# TAG
meta <- t(rbind(sapply(cov_cds_24[names(cov_cds_24) %in% tag_common], function(x){x[(length(x)-32):length(x)]}), sapply(cov_utr3_24[names(cov_utr3_24) %in% tag_common], function(x){x[1:30]})))
meta[is.na(meta)] <- 0
# normalize
meta <- meta / rowSums(meta)
meta[is.nan(meta)] <- 0
meta <- colSums(meta)
df <- data.frame(scale,meta)
p_24tag <- ggplot(df, aes(x=scale, y=meta)) + geom_bar(stat = "identity") + xlab("position from stop codon") + ylab("sum(Ribo-seq)") + ggtitle("2-4cell, TAG")
ggsave(file = "/Volumes/USELESS/META/SHAPES/stop/ribo_tag_24.png", plot = p_24tag)


# 256cell

taa_common <- Reduce(intersect, list(taa, names(cov_cds_256), names(cov_utr3_256), otpg))
tga_common <- Reduce(intersect, list(tga, names(cov_cds_256), names(cov_utr3_256), otpg))
tag_common <- Reduce(intersect, list(tag, names(cov_cds_256), names(cov_utr3_256), otpg))

# TAA
meta <- t(rbind(sapply(cov_cds_256[names(cov_cds_256) %in% taa_common], function(x){x[(length(x)-32):length(x)]}), sapply(cov_utr3_256[names(cov_utr3_256) %in% taa_common], function(x){x[1:30]})))
meta[is.na(meta)] <- 0
# normalize
meta <- meta / rowSums(meta)
meta[is.nan(meta)] <- 0
meta <- colSums(meta)
df <- data.frame(scale,meta)
p_256taa <- ggplot(df, aes(x=scale, y=meta)) + geom_bar(stat = "identity") + xlab("position from stop codon") + ylab("sum(Ribo-seq)") + ggtitle("256cell, TAA")
ggsave(file = "/Volumes/USELESS/META/SHAPES/stop/ribo_taa_256.png", plot = p_256taa)

# TGA
meta <- t(rbind(sapply(cov_cds_256[names(cov_cds_256) %in% tga_common], function(x){x[(length(x)-32):length(x)]}), sapply(cov_utr3_256[names(cov_utr3_256) %in% tga_common], function(x){x[1:30]})))
meta[is.na(meta)] <- 0
# normalize
meta <- meta / rowSums(meta)
meta[is.nan(meta)] <- 0
meta <- colSums(meta)
df <- data.frame(scale,meta)
p_256tga <- ggplot(df, aes(x=scale, y=meta)) + geom_bar(stat = "identity") + xlab("position from stop codon") + ylab("sum(Ribo-seq)") + ggtitle("256cell, TGA")
ggsave(file = "/Volumes/USELESS/META/SHAPES/stop/ribo_tga_256.png", plot = p_256tga)

# TAG
meta <- t(rbind(sapply(cov_cds_256[names(cov_cds_256) %in% tag_common], function(x){x[(length(x)-32):length(x)]}), sapply(cov_utr3_256[names(cov_utr3_256) %in% tag_common], function(x){x[1:30]})))
meta[is.na(meta)] <- 0
# normalize
meta <- meta / rowSums(meta)
meta[is.nan(meta)] <- 0
meta <- colSums(meta)
df <- data.frame(scale,meta)
p_256tag <- ggplot(df, aes(x=scale, y=meta)) + geom_bar(stat = "identity") + xlab("position from stop codon") + ylab("sum(Ribo-seq)") + ggtitle("256cell, TAG")
ggsave(file = "/Volumes/USELESS/META/SHAPES/stop/ribo_tag_256.png", plot = p_256tag)

###############################
######## LOGO #################
###############################
library(seqinr)
fasta_cdna <- read.fasta("/Volumes/USELESS/DATA/fasta/danio_rerio/GRCz10/cdna/Danio_rerio.GRCz10.cdna.all.fa")
names(fasta_cdna) <- sapply(names(fasta_cdna), function(x){substr(x,1,18)})

fasta <- fasta_cdna[names(fasta_cdna) %in% names(cds)]

logo <- lapply(names(fasta), function(x){fasta[[x]][(txLengths[x,]$utr5_len+txLengths[x,]$cds_len-29):(txLengths[x,]$utr5_len+txLengths[x,]$cds_len+33)]})

logo_taa <- logo[sapply(logo, function(x){paste(x[28:30], collapse = "") == "taa"})]
logo_taa <- sapply(logo_taa, function(x){paste(x, collapse = "")})
logo_tga <- logo[sapply(logo, function(x){paste(x[28:30], collapse = "") == "tga"})]
logo_tga <- sapply(logo_tga, function(x){paste(x, collapse = "")})
logo_tag <- logo[sapply(logo, function(x){paste(x[28:30], collapse = "") == "tag"})]
logo_tag <- sapply(logo_tag, function(x){paste(x, collapse = "")})

fileConn<-file("/Volumes/USELESS/META/SHAPES/codon/stop_codon/logo_taa.txt")
writeLines(logo_taa, fileConn)
close(fileConn)

fileConn<-file("/Volumes/USELESS/META/SHAPES/codon/stop_codon/logo_tga.txt")
writeLines(logo_tga, fileConn)
close(fileConn)

fileConn<-file("/Volumes/USELESS/META/SHAPES/codon/stop_codon/logo_tag.txt")
writeLines(logo_tag, fileConn)
close(fileConn)

##########
# TAA/TGA/TAG changes from 2-4cell to 256cell (TE24 vs TE256, facet on stop codon)

load(file = "/Volumes/USELESS/META/SHAPES/FPKM_24.Rdata")
load(file = "/Volumes/USELESS/META/SHAPES/FPKM_256.Rdata")

TE <- data.frame(te24 = FPKM_24$exons_Ribo_fpkm / FPKM_24$exons_RNA_fpkm, te256 = FPKM_256$exons_Ribo_fpkm / FPKM_256$exons_RNA_fpkm)
rownames(TE) <- rownames(FPKM_24)
TE <- TE[complete.cases(TE),]

te_taa <- TE[rownames(TE) %in% taa,]
te_taa$stop_codon <- rep("taa", nrow(te_taa))
te_tga <- TE[rownames(TE) %in% tga,]
te_tga$stop_codon <- rep("tga", nrow(te_tga))
te_tag <- TE[rownames(TE) %in% tag,]
te_tag$stop_codon <- rep("tag", nrow(te_tag))

te <- rbind(te_taa, te_tga, te_tag)

ggplot(te, aes(x = te24, y = te256)) + geom_point() + facet_wrap(~ stop_codon) + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES/stop/TE24_vs_TE256_for_TAA_TAG_TGA.png")

