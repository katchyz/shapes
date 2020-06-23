### get 1000 tx with decent fpkm (otpg), plot to unit length: utr5, cds, utr3

## SHAPE
## Ribo-Seq
## CAI

library(data.table)
library(ggplot2)
library(GenomicRanges)
library(GenomicFeatures)
library(plyr)

### get transcript database
txdb_can <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/danio_rerio/Danio_rerio.GRCz10.84.gtf", format = "gtf")
txLengths <- transcriptLengths(txdb_can, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
rownames(txLengths) <- txLengths$tx_name
txLengths <- txLengths[order(rownames(txLengths)),]

# get longest transcript per gene
txlen <- arrange(txLengths, gene_id, desc(tx_len))
txlen <- txlen[!duplicated(txlen$gene_id),]

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

