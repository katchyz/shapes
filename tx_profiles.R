### get 1000 tx with decent fpkm (otpg), plot to unit length: utr5, cds, utr3

## SHAPE
## Ribo-Seq
## CAI

library(data.table)
library(ggplot2)
library(GenomicRanges)
library(GenomicFeatures)
library(plyr)
library(GenomicAlignments)
library(rtracklayer)
library(seqinr)

### get transcript database
txdb_can <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/danio_rerio/Danio_rerio.GRCz10.84.gtf", format = "gtf")
txLengths <- transcriptLengths(txdb_can, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
rownames(txLengths) <- txLengths$tx_name
txLengths <- txLengths[order(rownames(txLengths)),]

# get longest transcript per gene
txlen <- arrange(txLengths, gene_id, desc(tx_len))
txlen <- txlen[!duplicated(txlen$gene_id),]

exons <- exonsBy(txdb_can, by="tx", use.names=TRUE)
exons <- exons[order(names(exons))]
names(exons) <- sapply(names(exons), function(x){substr(x,1,18)})

### load data
load(file="/Volumes/USELESS/META/SHAPES_NEW/general/clouds/cloud_df.Rsave")
# shape
load(file="/Volumes/USELESS/DATA/Shape-Seq/shape_cell1.Rsave")
load(file="/Volumes/USELESS/DATA/Shape-Seq/shape_cell24.Rsave")
load(file="/Volumes/USELESS/DATA/Shape-Seq/shape_cell256.Rsave")
load(file="/Volumes/USELESS/DATA/Shape-Seq/shape_cell1K.Rsave")
load(file="/Volumes/USELESS/DATA/Shape-Seq/shape_oblong.Rsave")
load(file="/Volumes/USELESS/DATA/Shape-Seq/shape_oblong_CHX.Rsave")

# shapes
shapes_cell1_invitro <- get(load(file = "/Volumes/USELESS/test/normalized/cell1_invitro_NAI_N3.Rsave"))
shapes_cell24 <- get(load(file = "/Volumes/USELESS/test/normalized/cell2_4_invivo_NAI_N3.Rsave"))
rm(shapes_norm)

### get features

# pick tx with utr5, cds, utr3
# subset on rna_fpkm (or SHAPE?)
lct <- txlen[(txlen$utr5_len > 30) & (txlen$cds_len > 100) & (txlen$utr3_len > 30), ]

lct <- lct[lct$tx_name %in% df[df$rna_cell24 > 1 & df$rna_cell256 > 1 & df$rna_cell1K > 1,]$n,]
lct <- lct[lct$tx_name %in% df[df$shape_cell24_log2ratio_internal > 0.01 
                               & df$shape_cell256_log2ratio_internal > 0.01
                               & df$shape_cell1K_log2ratio_internal > 0.01,]$n,]

sam <- sample(lct$tx_name, 500)



################
normalizeTx <- function(shape, txlen, sam){
  # subset shape GRanges
  shape <- shape[names(shape) %in% sam]
  shape <- shape[order(names(shape))]
  len <- txlen[txlen$tx_name %in% names(shape),]
  len <- len[order(len$tx_name),]
  
  leader <- mapply(function(x,y){x$log2ratio[1:y]}, shape, len$utr5_len)
  shape <- mapply(function(x,y){x$log2ratio[(y+1):length(x)]}, shape, len$utr5_len)
  cds <- mapply(function(x,y){x[1:y]}, shape, len$cds_len)
  trailer <- mapply(function(x,y){x[(y+1):length(x)]}, shape, len$cds_len)
  
  l <- leader[[1]]
  c <- cds[[1]]
  t <- trailer[[1]]
  
  # leader
  sev <- data.table(
    start = head(seq(0,1,(1/length(l))), -1),
    end = tail(seq(0,1,(1/length(l))), -1),
    value = l / (sum(l) + sum(c) + sum(t))  # values, normalized
  )
  setkey(sev, start, end)
  for (name in tail(names(leader), -1)) {
    lnew <- leader[[name]]
    cnew <- cds[[name]]
    tnew <- trailer[[name]]
    sevnew <- data.table(
      start = head(seq(0,1,(1/length(lnew))), -1),
      end = tail(seq(0,1,(1/length(lnew))), -1),
      value = lnew / (sum(lnew) + sum(cnew) + sum(tnew))  # values, normalized
    )
    setkey(sevnew, start, end)
    over <- foverlaps(sev, sevnew)
    over2 <- data.table(
      start = over[, ifelse(start > i.start, start, i.start)],
      end = over[, ifelse(end < i.end, end, i.end)],
      value = over[, value + i.value])
    # merge some intervals ? (smooth over, so to avoid running out of significant digits)
    sev <- over2
    sev <- sev[sev$start < sev$end]
    setkey(sev, start, end)
  }
  sevL <- sev
  # CDS
  sev <- data.table(
    start = head(seq(0,1,(1/length(c))), -1),
    end = tail(seq(0,1,(1/length(c))), -1),
    value = c / (sum(l) + sum(c) + sum(t))  # values, normalized
  )
  setkey(sev, start, end)
  for (name in tail(names(cds), -1)) {
    lnew <- leader[[name]]
    cnew <- cds[[name]]
    tnew <- trailer[[name]]
    sevnew <- data.table(
      start = head(seq(0,1,(1/length(cnew))), -1),
      end = tail(seq(0,1,(1/length(cnew))), -1),
      value = cnew / (sum(lnew) + sum(cnew) + sum(tnew))  # values, normalized
    )
    setkey(sevnew, start, end)
    over <- foverlaps(sev, sevnew)
    over2 <- data.table(
      start = over[, ifelse(start > i.start, start, i.start)],
      end = over[, ifelse(end < i.end, end, i.end)],
      value = over[, value + i.value])
    # merge some intervals ? (smooth over, so to avoid running out of significant digits)
    sev <- over2
    sev <- sev[sev$start < sev$end]
    setkey(sev, start, end)
  }
  sevC <- sev
  sevC$start <- sevC$start + 1
  sevC$end <- sevC$end + 1
  # trailer
  sev <- data.table(
    start = head(seq(0,1,(1/length(t))), -1),
    end = tail(seq(0,1,(1/length(t))), -1),
    value = t / (sum(l) + sum(c) + sum(t))  # values, normalized
  )
  setkey(sev, start, end)
  for (name in tail(names(trailer), -1)) {
    lnew <- leader[[name]]
    cnew <- cds[[name]]
    tnew <- trailer[[name]]
    sevnew <- data.table(
      start = head(seq(0,1,(1/length(tnew))), -1),
      end = tail(seq(0,1,(1/length(tnew))), -1),
      value = tnew / (sum(lnew) + sum(cnew) + sum(tnew))  # values, normalized
    )
    setkey(sevnew, start, end)
    over <- foverlaps(sev, sevnew)
    over2 <- data.table(
      start = over[, ifelse(start > i.start, start, i.start)],
      end = over[, ifelse(end < i.end, end, i.end)],
      value = over[, value + i.value])
    # merge some intervals ? (smooth over, so to avoid running out of significant digits)
    sev <- over2
    sev <- sev[sev$start < sev$end]
    setkey(sev, start, end)
  }
  sevT <- sev
  sevT$start <- sevT$start + 2
  sevT$end <- sevT$end + 2
  sev <- rbind(sevL, sevC, sevT)
  return(sev)
}

sam <- sample(lct$tx_name, 1000)
sam <- lct$tx_name

sev <- normalizeTx(shape_cell1, txlen, sam)
ggplot() + geom_step(data = sev, aes(start, value), size=0.1) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) + xlab("leader, CDS, trailer") + ylab("sum(SHAPE)") + ggtitle("1Cell")
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/tx_profiles/shape1.png")

sev <- normalizeTx(shape_cell24, txlen, sam)
ggplot() + geom_step(data = sev, aes(start, value), size=0.1) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) + xlab("leader, CDS, trailer") + ylab("sum(SHAPE)") + ggtitle("2-4Cell")
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/tx_profiles/shape24.png")

sev <- normalizeTx(shape_cell256, txlen, sam)
ggplot() + geom_step(data = sev, aes(start, value), size=0.1) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) + xlab("leader, CDS, trailer") + ylab("sum(SHAPE)") + ggtitle("256Cell")
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/tx_profiles/shape256.png")

sev <- normalizeTx(shape_cell1K, txlen, sam)
ggplot() + geom_step(data = sev, aes(start, value), size=0.1) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) + xlab("leader, CDS, trailer") + ylab("sum(SHAPE)") + ggtitle("1KCell")
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/tx_profiles/shape1K.png")

sev <- normalizeTx(shape_oblong, txlen, sam)
ggplot() + geom_step(data = sev, aes(start, value), size=0.1) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) + xlab("leader, CDS, trailer") + ylab("sum(SHAPE)") + ggtitle("oblong")
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/tx_profiles/shape_oblong.png")

sev <- normalizeTx(shape_oblong_CHX, txlen, sam)
ggplot() + geom_step(data = sev, aes(start, value), size=0.1) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) + xlab("leader, CDS, trailer") + ylab("sum(SHAPE)") + ggtitle("oblong CHX")
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/tx_profiles/shape_oblong_CHX.png")



#########
txdb_can <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/danio_rerio/zebrafish_GRCz10.gtf", format = "gtf")
exons <- exonsBy(txdb_can, by="tx", use.names=TRUE)
exons <- exons[order(names(exons))]
names(exons) <- sapply(names(exons), function(x){substr(x,1,18)})

### plot Ribo-seq

## get data

### cell24
# read BigWigs
# substitute strand, merge
bw24fw <- import.bw("/Volumes/USELESS/OUT/zebrafish_our_out/GRCz10/2-4Cell_fw.bw")
strand(bw24fw) <- rep("+", length(bw24fw))
bw24rv <- import.bw("/Volumes/USELESS/OUT/zebrafish_our_out/GRCz10/2-4Cell_rv.bw")
strand(bw24rv) <- rep("-", length(bw24rv))
bw24 <- c(bw24fw, bw24rv)

# map to transcripts, subset exons
ribo24 <- mapToTranscripts(bw24, exons, ignore.strand = FALSE)
mcols(ribo24) <- cbind(mcols(ribo24), DataFrame(bw24[ribo24$xHits]))
ribo24 <- split(ribo24, seqnames(ribo24))
ribo24 <- ribo24[names(ribo24) %in% rownames(txLengths)]
ribo24 <- ribo24[order(names(ribo24))]

sl <- txLengths[rownames(txLengths) %in% names(ribo24),]
sl <- sl[order(match(rownames(sl), names(seqlengths(ribo24)))),]

seqlevels(ribo24) <- seqlevels(ribo24)[seqlevels(ribo24) %in% rownames(sl)]
seqinfo(ribo24) <- seqinfo(ribo24)[rownames(sl)]
seqlengths(ribo24) <- sl$tx_len


# exons
cov <- coverage(ribo24)
cov <- cov[order(names(cov))]
plus <- ribo24[strand(ribo24) == "+"]
plus <- plus[sapply(plus, function(x){length(x) > 0})]
minus <- ribo24[strand(ribo24) == "-"]
minus <- minus[sapply(minus, function(x){length(x) > 0})]
# whole vector
# for - strand:rev()
#as.vector(cov[["ENSDART00000162928"]])
cov_plus <- sapply(cov[names(cov) %in% names(plus)], function(x){as.vector(x)})
cov_minus <- sapply(cov[names(cov) %in% names(minus)], function(x){rev(as.vector(x))})
for (name in names(cov_plus)) {
  cov_plus[[name]][cov_plus[[name]] == 1] <- ribo24[[name]]$score
}
for (name in names(cov_minus)) {
  cov_minus[[name]][cov_minus[[name]] == 1] <- rev(ribo24[[name]]$score)
}
cov <- c(cov_plus, cov_minus)
cov <- cov[order(names(cov))]
ribo24 <- cov

save(ribo24, file="/Volumes/USELESS/DATA/Ribo-Seq/ribo24.Rsave")



### cell256
# read BigWigs
# substitute strand, merge
bw256fw <- import.bw("/Volumes/USELESS/OUT/zebrafish_our_out/GRCz10/256Cell_fw.bw")
strand(bw256fw) <- rep("+", length(bw256fw))
bw256rv <- import.bw("/Volumes/USELESS/OUT/zebrafish_our_out/GRCz10/256Cell_rv.bw")
strand(bw256rv) <- rep("-", length(bw256rv))
bw256 <- c(bw256fw, bw256rv)

# map to transcripts, subset exons
ribo256 <- mapToTranscripts(bw256, exons, ignore.strand = FALSE)
mcols(ribo256) <- cbind(mcols(ribo256), DataFrame(bw256[ribo256$xHits]))
ribo256 <- split(ribo256, seqnames(ribo256))
ribo256 <- ribo256[names(ribo256) %in% rownames(txLengths)]
ribo256 <- ribo256[order(names(ribo256))]

sl <- txLengths[rownames(txLengths) %in% names(ribo256),]
sl <- sl[order(match(rownames(sl), names(seqlengths(ribo256)))),]

seqlevels(ribo256) <- seqlevels(ribo256)[seqlevels(ribo256) %in% rownames(sl)]
seqinfo(ribo256) <- seqinfo(ribo256)[rownames(sl)]
seqlengths(ribo256) <- sl$tx_len

# exons
cov <- coverage(ribo256)
cov <- cov[order(names(cov))]
plus <- ribo256[strand(ribo256) == "+"]
plus <- plus[sapply(plus, function(x){length(x) > 0})]
minus <- ribo256[strand(ribo256) == "-"]
minus <- minus[sapply(minus, function(x){length(x) > 0})]
# whole vector
# for - strand:rev()
#as.vector(cov[["ENSDART00000162928"]])
cov_plus <- sapply(cov[names(cov) %in% names(plus)], function(x){as.vector(x)})
cov_minus <- sapply(cov[names(cov) %in% names(minus)], function(x){rev(as.vector(x))})
for (name in names(cov_plus)) {
  cov_plus[[name]][cov_plus[[name]] == 1] <- ribo256[[name]]$score
}
for (name in names(cov_minus)) {
  cov_minus[[name]][cov_minus[[name]] == 1] <- rev(ribo256[[name]]$score)
}
cov <- c(cov_plus, cov_minus)
cov <- cov[order(names(cov))]
ribo256 <- cov

save(ribo256, file="/Volumes/USELESS/DATA/Ribo-Seq/ribo256.Rsave")



### cell1K
# read BigWigs
# substitute strand, merge
bw1Kfw <- import.bw("/Volumes/USELESS/OUT/zebrafish_our_out/GRCz10/1KCell_fw.bw")
strand(bw1Kfw) <- rep("+", length(bw1Kfw))
bw1Krv <- import.bw("/Volumes/USELESS/OUT/zebrafish_our_out/GRCz10/1KCell_rv.bw")
strand(bw1Krv) <- rep("-", length(bw1Krv))
bw1K <- c(bw1Kfw, bw1Krv)

# map to transcripts, subset exons
ribo1K <- mapToTranscripts(bw1K, exons, ignore.strand = FALSE)
mcols(ribo1K) <- cbind(mcols(ribo1K), DataFrame(bw1K[ribo1K$xHits]))
ribo1K <- split(ribo1K, seqnames(ribo1K))
ribo1K <- ribo1K[names(ribo1K) %in% rownames(txLengths)]
ribo1K <- ribo1K[order(names(ribo1K))]

sl <- txLengths[rownames(txLengths) %in% names(ribo1K),]
sl <- sl[order(match(rownames(sl), names(seqlengths(ribo1K)))),]

seqlevels(ribo1K) <- seqlevels(ribo1K)[seqlevels(ribo1K) %in% rownames(sl)]
seqinfo(ribo1K) <- seqinfo(ribo1K)[rownames(sl)]
seqlengths(ribo1K) <- sl$tx_len

# exons
cov <- coverage(ribo1K)
cov <- cov[order(names(cov))]
plus <- ribo1K[strand(ribo1K) == "+"]
plus <- plus[sapply(plus, function(x){length(x) > 0})]
minus <- ribo1K[strand(ribo1K) == "-"]
minus <- minus[sapply(minus, function(x){length(x) > 0})]
# whole vector
# for - strand:rev()
#as.vector(cov[["ENSDART00000162928"]])
cov_plus <- sapply(cov[names(cov) %in% names(plus)], function(x){as.vector(x)})
cov_minus <- sapply(cov[names(cov) %in% names(minus)], function(x){rev(as.vector(x))})
for (name in names(cov_plus)) {
  cov_plus[[name]][cov_plus[[name]] == 1] <- ribo1K[[name]]$score
}
for (name in names(cov_minus)) {
  cov_minus[[name]][cov_minus[[name]] == 1] <- rev(ribo1K[[name]]$score)
}
cov <- c(cov_plus, cov_minus)
cov <- cov[order(names(cov))]
ribo1K <- cov

save(ribo1K, file="/Volumes/USELESS/DATA/Ribo-Seq/ribo1K.Rsave")

#####################
#####################
#####################

normalizeTxRibo <- function(ribo, txlen, sam){
  # subset shape GRanges
  ribo <- ribo[names(ribo) %in% sam]
  ribo <- ribo[order(names(ribo))]
  len <- txlen[txlen$tx_name %in% names(ribo),]
  len <- len[order(len$tx_name),]
  
  leader <- mapply(function(x,y){x[1:y]}, ribo, len$utr5_len)
  ribo <- mapply(function(x,y){x[(y+1):length(x)]}, ribo, len$utr5_len)
  cds <- mapply(function(x,y){x[1:y]}, ribo, len$cds_len)
  trailer <- mapply(function(x,y){x[(y+1):length(x)]}, ribo, len$cds_len)
  
  l <- leader[[1]]
  c <- cds[[1]]
  t <- trailer[[1]]
  
  # leader
  sev <- data.table(
    start = head(seq(0,1,(1/length(l))), -1),
    end = tail(seq(0,1,(1/length(l))), -1),
    value = l / (sum(l) + sum(c) + sum(t))  # values, normalized
  )
  setkey(sev, start, end)
  for (name in tail(names(leader), -1)) {
    lnew <- leader[[name]]
    cnew <- cds[[name]]
    tnew <- trailer[[name]]
    sevnew <- data.table(
      start = head(seq(0,1,(1/length(lnew))), -1),
      end = tail(seq(0,1,(1/length(lnew))), -1),
      value = lnew / (sum(lnew) + sum(cnew) + sum(tnew))  # values, normalized
    )
    setkey(sevnew, start, end)
    over <- foverlaps(sev, sevnew)
    over2 <- data.table(
      start = over[, ifelse(start > i.start, start, i.start)],
      end = over[, ifelse(end < i.end, end, i.end)],
      value = over[, value + i.value])
    # merge some intervals ? (smooth over, so to avoid running out of significant digits)
    sev <- over2
    sev <- sev[sev$start < sev$end]
    setkey(sev, start, end)
  }
  sevL <- sev
  # CDS
  sev <- data.table(
    start = head(seq(0,1,(1/length(c))), -1),
    end = tail(seq(0,1,(1/length(c))), -1),
    value = c / (sum(l) + sum(c) + sum(t))  # values, normalized
  )
  setkey(sev, start, end)
  for (name in tail(names(cds), -1)) {
    lnew <- leader[[name]]
    cnew <- cds[[name]]
    tnew <- trailer[[name]]
    sevnew <- data.table(
      start = head(seq(0,1,(1/length(cnew))), -1),
      end = tail(seq(0,1,(1/length(cnew))), -1),
      value = cnew / (sum(lnew) + sum(cnew) + sum(tnew))  # values, normalized
    )
    setkey(sevnew, start, end)
    over <- foverlaps(sev, sevnew)
    over2 <- data.table(
      start = over[, ifelse(start > i.start, start, i.start)],
      end = over[, ifelse(end < i.end, end, i.end)],
      value = over[, value + i.value])
    # merge some intervals ? (smooth over, so to avoid running out of significant digits)
    sev <- over2
    sev <- sev[sev$start < sev$end]
    setkey(sev, start, end)
  }
  sevC <- sev
  sevC$start <- sevC$start + 1
  sevC$end <- sevC$end + 1
  # trailer
  sev <- data.table(
    start = head(seq(0,1,(1/length(t))), -1),
    end = tail(seq(0,1,(1/length(t))), -1),
    value = t / (sum(l) + sum(c) + sum(t))  # values, normalized
  )
  setkey(sev, start, end)
  for (name in tail(names(trailer), -1)) {
    lnew <- leader[[name]]
    cnew <- cds[[name]]
    tnew <- trailer[[name]]
    sevnew <- data.table(
      start = head(seq(0,1,(1/length(tnew))), -1),
      end = tail(seq(0,1,(1/length(tnew))), -1),
      value = tnew / (sum(lnew) + sum(cnew) + sum(tnew))  # values, normalized
    )
    setkey(sevnew, start, end)
    over <- foverlaps(sev, sevnew)
    over2 <- data.table(
      start = over[, ifelse(start > i.start, start, i.start)],
      end = over[, ifelse(end < i.end, end, i.end)],
      value = over[, value + i.value])
    # merge some intervals ? (smooth over, so to avoid running out of significant digits)
    sev <- over2
    sev <- sev[sev$start < sev$end]
    setkey(sev, start, end)
  }
  sevT <- sev
  sevT$start <- sevT$start + 2
  sevT$end <- sevT$end + 2
  sev <- rbind(sevL, sevC, sevT)
  return(sev)
}

sam <- sample(lct$tx_name, 1000)
sam <- lct$tx_name

sev <- normalizeTxRibo(ribo24, txlen, sam)
ggplot() + geom_step(data = sev, aes(start, value), size=0.1) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) + xlab("leader, CDS, trailer") + ylab("sum(Ribo)") + ggtitle("2-4Cell")
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/tx_profiles/ribo/ribo24.png")

sev <- normalizeTxRibo(ribo256, txlen, sam)
ggplot() + geom_step(data = sev, aes(start, value), size=0.1) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) + xlab("leader, CDS, trailer") + ylab("sum(Ribo)") + ggtitle("256Cell")
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/tx_profiles/ribo/ribo256.png")

sev <- normalizeTxRibo(ribo1K, txlen, sam)
ggplot() + geom_step(data = sev, aes(start, value), size=0.1) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) + xlab("leader, CDS, trailer") + ylab("sum(Ribo)") + ggtitle("1KCell")
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/tx_profiles/ribo/ribo1K.png")


####### profile on RNAfold
fasta_vienna <- read.fasta("/Volumes/USELESS/DATA/fasta/danio_rerio/GRCz10/cdna/structures_only.txt")
fasta_vienna <- fasta_vienna[names(fasta_vienna) %in% sam]
fasta_vienna <- sapply(fasta_vienna, function(x){recode(getSequence(x),"'.'=1;'('=0;')'=0")})


sev <- normalizeTxRibo(fasta_vienna, txlen, sam)
ggplot() + geom_step(data = sev, aes(start, value), size=0.1) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) + xlab("leader, CDS, trailer") + ylab("sum(RNAfold)") + ggtitle("RNAfold ('.'=1 '(',')'=0)")
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/tx_profiles/RNAfold.png")



####### clouds vs the rest
##########
### plot cloud and main population separately
tx <- df[df$cloud_cell24 == "cloud" & !is.na(df$cloud_cell24),]$n
tx <- tx[tx %in% txlen$tx_name]
s <- shape_cell24[names(shape_cell24) %in% tx]
s <- s[sapply(s, function(x){sum(x$log2ratio) > 0})]
sam <- names(s)
sev <- normalizeTx(s, txlen, sam)
ggplot() + geom_step(data = sev, aes(start, value), size=0.1) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) + xlab("leader, CDS, trailer") + ylab("sum(SHAPE)") + ggtitle("cloud, cell24")
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/tx_profiles/cloud/cloud_cell24.png")

#ribo <- ribo[sapply(ribo, function(x){sum(x) > 0})]
sev <- normalizeTxRibo(ribo24, txlen, sam)
ggplot() + geom_step(data = sev, aes(start, value), size=0.1) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) + xlab("leader, CDS, trailer") + ylab("sum(Ribo)") + ggtitle("2-4Cell")
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/tx_profiles/ribo/ribo24.png")