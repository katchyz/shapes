### make plots as in papers
### -50/+100 start, -100/+50 stop; line plot
library(GenomicFeatures)
library(plyr)
library(ggplot2)
library(dplyr)

load(file = "/Volumes/USELESS/META/SHAPES/SHAPE_FPKM.Rdata")

### get transcript dtabase
txdb_can <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/danio_rerio/zebrafish_GRCz10.gtf", format = "gtf")
exons <- exonsBy(txdb_can, by="tx", use.names=TRUE)
exons <- exons[order(names(exons))]

txLengths <- transcriptLengths(txdb_can, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
rownames(txLengths) <- txLengths$tx_name
txLengths <- txLengths[order(rownames(txLengths)),]
rownames(txLengths) <- sapply(rownames(txLengths), function(x){substr(x,1,18)})

exons24 <- subsetByOverlaps(gc_unsmoothed_2_4, exons)
exons24 <- split(exons24, exons24$trnames)
names(exons24) <- sapply(names(exons24), function(x){substr(x,1,18)})

exons256 <- subsetByOverlaps(gc_unsmoothed_256, exons)
exons256 <- split(exons256, exons256$trnames)
names(exons256) <- sapply(names(exons256), function(x){substr(x,1,18)})

##
# subset on shape fpkm
wes <- rownames(SHAPE_FPKM[SHAPE_FPKM$stage_256cell > 1 & SHAPE_FPKM$stage_2_4cell > 1,])
# subset on utr5 and cds length
txl <- txLengths[rownames(txLengths) %in% wes,]
txl <- rownames(txl[txl$utr5_len > 50 & txl$cds_len > 100,])

# pick the same transcripts in 2-4cell and 256cell
t24 <- exons24[names(exons24) %in% txl]
t256 <- exons256[names(exons256) %in% txl]

t24 <- t24[sapply(t24, function(x){sum(x$dtcr) > 0})]
t256 <- t256[sapply(t256, function(x){sum(x$dtcr) > 0})]

### OTPG
txdb_can <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/danio_rerio/Danio_rerio.GRCz10.84.gtf", format = "gtf")
txLengths <- transcriptLengths(txdb_can, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
rownames(txLengths) <- txLengths$tx_name
txLengths <- txLengths[order(rownames(txLengths)),]
rownames(txLengths) <- sapply(rownames(txLengths), function(x){substr(x,1,18)})
# get longest transcript per gene
txlen <- arrange(txLengths, gene_id, desc(tx_len))
txlen <- txlen[!duplicated(txlen$gene_id),]
# get transcripts with at least 33 on CDS and 30 on utr3
otpg <- sort(txlen[(txlen$cds_len > 100) & (txlen$utr3_len > 50) & (txlen$utr5_len > 50), ]$tx_name)

###
t24 <- t24[names(t24) %in% otpg]
t256 <- t256[names(t256) %in% otpg]

t24 <- t24[names(t24) %in% names(t256)]
t256 <- t256[names(t256) %in% names(t24)]

# plot START
meta24 <- Reduce("+", lapply(names(t24), function(x){t24[[x]][(txLengths[x,]$utr5_len-49):(txLengths[x,]$utr5_len+100)]$dtcr / sum(t24[[x]]$dtcr)}))

t256[["ENSDART00000054202"]] <- NULL

meta256 <- Reduce("+", lapply(names(t256), function(x){t256[[x]][(txLengths[x,]$utr5_len-49):(txLengths[x,]$utr5_len+100)]$dtcr / sum(t256[[x]]$dtcr)}))

# plot STOP

n <- c()
for (name in names(t24)) {
  if (txLengths[name,]$utr5_len+txLengths[name,]$cds_len+50 > length(t24[[name]])) {
    print(name)
    n <- c(n, name)
  }
}

for (name in n) {
  t24[[name]] <- NULL
}

stop24 <- Reduce("+", lapply(names(t24), function(x){t24[[x]][(txLengths[x,]$utr5_len+txLengths[x,]$cds_len-99):(txLengths[x,]$utr5_len+txLengths[x,]$cds_len+50)]$dtcr / sum(t24[[x]]$dtcr)}))

for (name in n) {
  t256[[name]] <- NULL
}

stop256 <- Reduce("+", lapply(names(t256), function(x){t256[[x]][(txLengths[x,]$utr5_len+txLengths[x,]$cds_len-99):(txLengths[x,]$utr5_len+txLengths[x,]$cds_len+50)]$dtcr / sum(t256[[x]]$dtcr)}))


### plot (line)
m24 <- c(meta24, rep(NA,10), stop24)
m256 <- c(meta256, rep(NA,10), stop256)
scale <- c(c(-50:-1), c(1:100), rep(NA, 10), c(-100:-1), c(1:50))

### 2-4cell
#start
df24 <- data.frame(scale = c(-50:100)[-51], meta24)
p24 <- ggplot(df24, aes(x=scale, y=meta24)) + geom_line() + xlab("position from start codon") + ylab("SHAPE") + ggtitle("2-4cell") + ylim(1,11)
rect <- data.frame(xmin=1, xmax=3, ymin=-Inf, ymax=Inf)
p24 <- p24 + geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
              color="grey20",
              alpha=0.5,
              inherit.aes = FALSE)
ggsave(file = "/Volumes/USELESS/META/SHAPES/start/line24start.png", plot = p24)
#stop
df24 <- data.frame(scale = c(-100:50)[-101], stop24)
p24 <- ggplot(df24, aes(x=scale, y=stop24)) + geom_line() + xlab("position from stop codon") + ylab("SHAPE") + ggtitle("2-4cell") + ylim(1,11)
rect <- data.frame(xmin=-3, xmax=-1, ymin=-Inf, ymax=Inf)
p24 <- p24 + geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                       color="grey20",
                       alpha=0.5,
                       inherit.aes = FALSE)
ggsave(file = "/Volumes/USELESS/META/SHAPES/start/line24stop.png", plot = p24)

### 256cell
#start
df256 <- data.frame(scale = c(-50:100)[-51], meta256)
p256 <- ggplot(df256, aes(x=scale, y=meta256)) + geom_line() + xlab("position from start codon") + ylab("SHAPE") + ggtitle("256cell") + ylim(1,11)
rect <- data.frame(xmin=1, xmax=3, ymin=-Inf, ymax=Inf)
p256 <- p256 + geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                       color="grey20",
                       alpha=0.5,
                       inherit.aes = FALSE)
ggsave(file = "/Volumes/USELESS/META/SHAPES/start/line256start.png", plot = p256)
#stop
df256 <- data.frame(scale = c(-100:50)[-101], stop256)
p256 <- ggplot(df256, aes(x=scale, y=stop256)) + geom_line() + xlab("position from stop codon") + ylab("SHAPE") + ggtitle("256cell") + ylim(1,11)
rect <- data.frame(xmin=-3, xmax=-1, ymin=-Inf, ymax=Inf)
p256 <- p256 + geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                       color="grey20",
                       alpha=0.5,
                       inherit.aes = FALSE)
ggsave(file = "/Volumes/USELESS/META/SHAPES/start/line256stop.png", plot = p256)


###########################
###########################
### > structure at start correlating with TE (split by TE, make bar/line plots)

load(file = "/Volumes/USELESS/META/SHAPES/FPKM_24.Rdata")
load(file = "/Volumes/USELESS/META/SHAPES/FPKM_256.Rdata")

fpkm24 <- FPKM_24[rownames(FPKM_24) %in% names(t24),]
fpkm256 <- FPKM_256[rownames(FPKM_256) %in% names(t256),]
fpkm24 <- fpkm24[rownames(fpkm24) %in% rownames(fpkm256),]
te <- data.frame(te24 = fpkm24$exons_Ribo_fpkm/fpkm24$exons_RNA_fpkm, te256 = fpkm256$exons_Ribo_fpkm/fpkm256$exons_RNA_fpkm)
rownames(te) <- rownames(fpkm24)
te <- te[complete.cases(te),]

# split t24 into bins based on TE
nbins = 3
te$bin24 <- as.factor(ntile(te$te24, nbins))
te$bin256 <- as.factor(ntile(te$te256, nbins))

t24 <- t24[names(t24) %in% rownames(te)]
t256 <- t256[names(t256) %in% rownames(te)]

df <- rbind(t(te), sapply(names(t24), function(x){as.numeric(t24[[x]][(txLengths[x,]$utr5_len-49):(txLengths[x,]$utr5_len+100)]$dtcr / sum(t24[[x]]$dtcr))}))
df <- data.frame(t(df))
for (i in c(5:154)) {
  df[,i] <- as.numeric(df[,i])
}

metas24 <- by(df, df$bin24, function(x){as.numeric(colSums(x[,5:154]))})
metas256 <- by(df, df$bin256, function(x){as.numeric(colSums(x[,5:154]))})
### bin1 - lowest TE, bin3 - highest TE

df <- ldply(metas24, data.frame)
colnames(df) <- c("bin", "reactivities")
df$scale <- rep(c(-50:100)[-51],nbins)
p <- ggplot(df, aes(x = scale, y = reactivities, fill = bin)) + geom_line() + facet_wrap(~ bin) + xlab("position from start codon") + guides(fill=FALSE) + ggtitle("2-4cell, 1-lowest, 3-highest TE")
rect <- data.frame(xmin=1, xmax=3, ymin=-Inf, ymax=Inf)
p <- p + geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                         color="red",
                         alpha=0.5,
                         inherit.aes = FALSE)
ggsave(file = "/Volumes/USELESS/META/SHAPES/start/startByTE24.png", plot = p)

df <- ldply(metas256, data.frame)
colnames(df) <- c("bin", "reactivities")
df$scale <- rep(c(-50:100)[-51],nbins)
p <- ggplot(df, aes(x = scale, y = reactivities, fill = bin)) + geom_line() + facet_wrap(~ bin) + xlab("position from start codon") + guides(fill=FALSE) + ggtitle("256cell, 1-lowest, 3-highest TE")
rect <- data.frame(xmin=1, xmax=3, ymin=-Inf, ymax=Inf)
p <- p + geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                   color="red",
                   alpha=0.5,
                   inherit.aes = FALSE)
ggsave(file = "/Volumes/USELESS/META/SHAPES/start/startByTE256.png", plot = p)

