### stability
library(GenomicFeatures)
library(plyr)
library(rtracklayer)

mishima <- read.table(file = "/Volumes/USELESS/DATA/Mishima_2016/GSE71609_2015_MZT.txt", header = TRUE, sep = "\t")

# txdb <- makeTxDbFromGFF("/Users/kasia/Downloads/mishima.gtf", format = "gtf")
txdb <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/Danio_rerio.Zv9.74.gtf.gz", format = "gtf")
#cds <- cdsBy(txdb, by="tx", use.names=TRUE)
#fiveUTR <- fiveUTRsByTranscript(txdb, use.names=TRUE)
#threeUTR <- threeUTRsByTranscript(txdb, use.names=TRUE)

# >10 nt 5 UTR
# >100 nt ORF
# >50 nt 3 UTR
# (13471 genes)

gtf_annot <- import("/Volumes/USELESS/DATA/genomes/GTF/Danio_rerio.Zv9.74.gtf.gz", format = "gtf")
protein_coding <- unique(gtf_annot[gtf_annot$gene_biotype == "protein_coding"]$gene_id)

# otpg
txLengths <- transcriptLengths(txdb, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
rownames(txLengths) <- txLengths$tx_name
txLengths <- txLengths[order(rownames(txLengths)),]
# get longest transcript per gene
txlen <- arrange(txLengths, gene_id, desc(tx_len))
txlen <- txlen[!duplicated(txlen$gene_id),]

#nrow(txlen[txlen$utr5_len > 10 & txlen$cds_len > 100 & txlen$utr3_len > 50,])
#reliableGenes <- txlen[txlen$utr5_len > 10 & txlen$cds_len > 100 & txlen$utr3_len > 50,]$gene_id

nrow(txLengths[txLengths$utr5_len > 10 & txLengths$cds_len > 100 & txLengths$utr3_len > 50,])
rg <- txLengths[txLengths$utr5_len > 10 & txLengths$cds_len > 100 & txLengths$utr3_len > 50,]$gene_id

rg <- arrange(rg, gene_id, desc(tx_len))
rg <- rg[!duplicated(rg$gene_id),] # 13525

#ok5 <- fiveUTR[sapply(fiveUTR, function(x){sum(width(ranges(x))) > 10})]
#ok3 <- threeUTR[names(threeUTR) %in% names(ok5)]
#ok3 <- ok3[sapply(ok3, function(x){sum(width(ranges(x))) > 50})]
#okCDS <- cds[names(cds) %in% names(ok3)]
#okCDS <- okCDS[sapply(okCDS, function(x){sum(width(ranges(x))) > 100})]

# change ENSDART to ENSDARG      # names(okCDS)
#gtp <- read.table(file = "/Volumes/USELESS/DATA/genomes/gtp_mishima.txt", header = FALSE, sep = "\t")
#colnames(gtp) <- c("gene", "transcript", "protein")

#gtp <- gtp[gtp$transcript %in% names(okCDS),]

#reliableGenes <- as.character.factor(unique(gtp$gene))

### RPKM cutoff
rel <- mishima[mishima$gene_id %in% reliableGenes,]
rel <- rel[rel$value_1 > 1,] # 8170 (should be 8701)

as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
rel$log2.fold_change. <- as.numeric.factor(rel$log2.fold_change.)
rel$log2.fold_change..1 <- as.numeric.factor(rel$log2.fold_change..1)
rel$log2.fold_change..2 <- as.numeric.factor(rel$log2.fold_change..2)
rel$log2.fold_change..3 <- as.numeric.factor(rel$log2.fold_change..3)
rel$log2.fold_change..4 <- as.numeric.factor(rel$log2.fold_change..4)

rel <- rel[complete.cases(stable$value_1),]
rel <- rel[complete.cases(stable$value_2),]

stable <- rel[rel$log2.fold_change. > -0.3,]
stable <- stable[stable$log2.fold_change. < 0.3,]

# stable <- stable[stable$log2.fold_change..3 > -0.3,]
# stable <- stable[stable$log2.fold_change..3 < 0.3,]

# stable <- stable[stable$log2.fold_change..4 > -0.3,]
# stable <- stable[stable$log2.fold_change..4 < 0.3,]

stable <- stable[stable$status == "OK",]
stable <- stable[complete.cases(stable$status),]

########## screw Mishima
########## compare FPKM of 2-4cell and 256 cell

fpkm_2_4 <- read.csv("/Volumes/USELESS/META/SHAPES/fpkm/Pauli_rna/0_2to4Cell_fpkm.csv")
fpkm_256 <- read.csv("/Volumes/USELESS/META/SHAPES/fpkm/Pauli_rna/1_256Cell_fpkm.csv")

a <- fpkm_2_4$rkpm
b <- fpkm_256$rkpm
fpkm <- data.frame(a,b)
rownames(fpkm) <- fpkm_2_4$tx_id
colnames(fpkm) <- c("rpkm_2_4", "rpkm_256")
fpkm$log2fold_change <- log2(fpkm$rpkm_2_4 / fpkm$rpkm_256)

# cut-off on rpkm_256 > 1
fpkm <- fpkm[fpkm$rpkm_256 > 1,]

stable <- fpkm[fpkm$log2fold_change > -0.3,]
stable <- stable[stable$log2fold_change < 0.3,]

un <- rownames(fpkm[fpkm$log2fold_change < -1.5,])
un <- c(un, rownames(fpkm[fpkm$log2fold_change > 1.5,]))
un <- sort(un)

unstable <- fpkm[rownames(fpkm) %in% un,]

#### get SHAPE

stable_2_4 <- shape_2_4[rownames(shape_2_4) %in% rownames(stable),]
stable_2_4$type <- rep("stable_2_4", nrow(stable_2_4))
stable_256 <- shape_256[rownames(shape_256) %in% rownames(stable),]
stable_256$type <- rep("stable_256", nrow(stable_256))

unstable_2_4 <- shape_2_4[rownames(shape_2_4) %in% rownames(unstable),]
unstable_2_4$type <- rep("unstable_2_4", nrow(unstable_2_4))
unstable_256 <- shape_256[rownames(shape_256) %in% rownames(unstable),]
unstable_256$type <- rep("unstable_256", nrow(unstable_256))

su_2_4 <- rbind.data.frame(stable_2_4, unstable_2_4)
su_256 <- rbind.data.frame(stable_256, unstable_256)

stun_2_4 <- melt(su_2_4, id=c("type","cai", "te_cds"))
stun_2_4$variable <- as.character.factor(stun_2_4$variable)

ggplot(stun_2_4, aes(x = type, y = value, fill = type)) + geom_boxplot(outlier.size = NA) + scale_y_log10() + guides(fill=FALSE) + facet_wrap(~ variable, scales = "free_y") + ylab("avg shape reactivities") + xlab(NULL)
ggsave("/Volumes/USELESS/META/SHAPES/stability/stable_unstable_2_4.png")

stun_256 <- melt(su_256, id=c("type","cai", "te_cds"))
stun_256$variable <- as.character.factor(stun_256$variable)

ggplot(stun_256, aes(x = type, y = value, fill = type)) + geom_boxplot(outlier.size = NA) + scale_y_log10() + guides(fill=FALSE) + facet_wrap(~ variable, scales = "free_y") + ylab("avg shape reactivities") + xlab(NULL)
ggsave("/Volumes/USELESS/META/SHAPES/stability/stable_unstable_256.png")


##### periodicity
stable_2_4_list <- cds_2_4_list[names(cds_2_4_list) %in% rownames(stable_2_4)]
stable_2_4_list <- stable_2_4_list[lengths(stable_2_4_list) > 150]
meta <- lapply(stable_2_4_list, function(x){(x[1:150]$dtcr)/sum(x[1:150]$dtcr)})
meta <- rapply(meta, f=function(x) ifelse(is.na(x),0,x), how="replace")
meta <- Reduce("+", meta)
amplitudes <- abs(fft(meta))
amp <- amplitudes[2:(length(amplitudes)/2+1)]
periods <- 1/periodogram(meta, method = "clone")$freq
df <- data.frame(periods, amp)
p1 <- ggplot(df, aes(x=periods, y=amp)) + geom_line() + ggtitle("stable, 2_4cell") + xlim(1, 8)
ggsave(file = "/Volumes/USELESS/META/SHAPES/stability/FFT_stable_2_4.png", plot = p1)

stable_256_list <- cds_256_list[names(cds_256_list) %in% rownames(stable_256)]
stable_256_list <- stable_256_list[lengths(stable_256_list) > 150]
meta <- lapply(stable_256_list, function(x){(x[1:150]$dtcr)/sum(x[1:150]$dtcr)})
meta <- rapply(meta, f=function(x) ifelse(is.na(x),0,x), how="replace")
meta <- Reduce("+", meta)
amplitudes <- abs(fft(meta))
amp <- amplitudes[2:(length(amplitudes)/2+1)]
periods <- 1/periodogram(meta, method = "clone")$freq
df <- data.frame(periods, amp)
p2 <- ggplot(df, aes(x=periods, y=amp)) + geom_line() + ggtitle("stable, 256cell") + xlim(1, 8)
ggsave(file = "/Volumes/USELESS/META/SHAPES/stability/FFT_stable_256.png", plot = p2)

##
unstable_2_4_list <- cds_2_4_list[names(cds_2_4_list) %in% rownames(unstable_2_4)]
unstable_2_4_list <- unstable_2_4_list[lengths(unstable_2_4_list) > 150]
meta <- lapply(unstable_2_4_list, function(x){(x[1:150]$dtcr)/sum(x[1:150]$dtcr)})
meta <- rapply(meta, f=function(x) ifelse(is.na(x),0,x), how="replace")
meta <- Reduce("+", meta)
amplitudes <- abs(fft(meta))
amp <- amplitudes[2:(length(amplitudes)/2+1)]
periods <- 1/periodogram(meta, method = "clone")$freq
df <- data.frame(periods, amp)
p3 <- ggplot(df, aes(x=periods, y=amp)) + geom_line() + ggtitle("unstable, 2_4cell") + xlim(1, 8)
ggsave(file = "/Volumes/USELESS/META/SHAPES/stability/FFT_unstable_2_4.png", plot = p3)

unstable_256_list <- cds_256_list[names(cds_256_list) %in% rownames(unstable_256)]
unstable_256_list <- unstable_256_list[lengths(unstable_256_list) > 150]
meta <- lapply(unstable_256_list, function(x){(x[1:150]$dtcr)/sum(x[1:150]$dtcr)})
meta <- rapply(meta, f=function(x) ifelse(is.na(x),0,x), how="replace")
meta <- Reduce("+", meta)
amplitudes <- abs(fft(meta))
amp <- amplitudes[2:(length(amplitudes)/2+1)]
periods <- 1/periodogram(meta, method = "clone")$freq
df <- data.frame(periods, amp)
p4 <- ggplot(df, aes(x=periods, y=amp)) + geom_line() + ggtitle("unstable, 256cell") + xlim(1, 8)
ggsave(file = "/Volumes/USELESS/META/SHAPES/stability/FFT_unstable_256.png", plot = p4)


