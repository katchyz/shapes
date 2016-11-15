library(data.table)
library(ggplot2)
library(GenomicRanges)
library(GenomicFeatures)

### get transcript database
txdb_can <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/danio_rerio/zebrafish_GRCz10.gtf", format = "gtf")
cds <- cdsBy(txdb_can, by="tx", use.names=TRUE)
fiveUTR <- fiveUTRsByTranscript(txdb_can, use.names=TRUE)
threeUTR <- threeUTRsByTranscript(txdb_can, use.names=TRUE)
exons <- exonsBy(txdb_can, by="tx", use.names=TRUE)

names(cds) <- sapply(names(cds), function(x){substr(x,1,18)})
names(fiveUTR) <- sapply(names(fiveUTR), function(x){substr(x,1,18)})
names(threeUTR) <- sapply(names(threeUTR), function(x){substr(x,1,18)})
names(exons) <- sapply(names(exons), function(x){substr(x,1,18)})

txLengths <- transcriptLengths(txdb_can, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
rownames(txLengths) <- txLengths$tx_name
rownames(txLengths) <- sapply(rownames(txLengths), function(x){substr(x,1,18)})

### features

leader24 <- subsetByOverlaps(gc_unsmoothed_2_4, fiveUTR)
leader24 <- split(leader24, leader24$trnames)
names(leader24) <- sapply(names(leader24), function(x){substr(x,1,18)})

cds24 <- subsetByOverlaps(gc_unsmoothed_2_4, cds)
cds24 <- split(cds24, cds24$trnames)
names(cds24) <- sapply(names(cds24), function(x){substr(x,1,18)})

trailer24 <- subsetByOverlaps(gc_unsmoothed_2_4, threeUTR)
trailer24 <- split(trailer24, trailer24$trnames)
names(trailer24) <- sapply(names(trailer24), function(x){substr(x,1,18)})

exons24 <- subsetByOverlaps(gc_unsmoothed_2_4, exons)
exons24 <- split(exons24, exons24$trnames)
names(exons24) <- sapply(names(exons24), function(x){substr(x,1,18)})

leader256 <- subsetByOverlaps(gc_unsmoothed_256, fiveUTR)
leader256 <- split(leader256, leader256$trnames)
names(leader256) <- sapply(names(leader256), function(x){substr(x,1,18)})

cds256 <- subsetByOverlaps(gc_unsmoothed_256, cds)
cds256 <- split(cds256, cds256$trnames)
names(cds256) <- sapply(names(cds256), function(x){substr(x,1,18)})

trailer256 <- subsetByOverlaps(gc_unsmoothed_256, threeUTR)
trailer256 <- split(trailer256, trailer256$trnames)
names(trailer256) <- sapply(names(trailer256), function(x){substr(x,1,18)})

exons256 <- subsetByOverlaps(gc_unsmoothed_256, exons)
exons256 <- split(exons256, exons256$trnames)
names(exons256) <- sapply(names(exons256), function(x){substr(x,1,18)})

## DATA
cds24 <- cds_2_4[sapply(cds_2_4, function(x){length(x) > 100})]
cds24 <- cds24[sapply(cds24, function(x){sum(x$dtcr) > 0})]
#sam24 <- cds24[sapply(cds24, function(x){median(x$dtcr) > 0})]

cds256 <- cds_256[sapply(cds_256, function(x){length(x) > 100})]
cds256 <- cds256[sapply(cds256, function(x){sum(x$dtcr) > 0})]
#sam256 <- cds256[sapply(cds256, function(x){median(x$dtcr) > 0})]

##
# FPKM_24
load(file = "/Volumes/USELESS/META/SHAPES/FPKM_24.Rdata")
# FPKM_256
load(file = "/Volumes/USELESS/META/SHAPES/FPKM_256.Rdata")

# leaders, CDSs, trailers (all present)
lct <- rownames(txLengths[(txLengths$utr5_len > 30) & (txLengths$cds_len > 100) & (txLengths$utr3_len > 30), ])

# FPKM cut-off
fpkm <- FPKM_24[rownames(FPKM_24) %in% lct,]
fpkm24 <- rownames(fpkm[(fpkm$utr5_RNA_fpkm > 1) & (fpkm$cds_RNA_fpkm > 1) & (fpkm$utr3_RNA_fpkm > 1),])
# choose one transcript per gene ?
fpkm <- FPKM_256[rownames(FPKM_256) %in% lct,]
fpkm256 <- rownames(fpkm[(fpkm$utr5_RNA_fpkm > 1) & (fpkm$cds_RNA_fpkm > 1) & (fpkm$utr3_RNA_fpkm > 1),])


### get sample for leaders, CDSs and trailers
leader24_fpkm <- leader24[names(leader24) %in% fpkm24]
leader24_fpkm <- leader24_fpkm[sapply(leader24_fpkm, function(x){length(x) > 30})]
leader24_fpkm <- leader24_fpkm[sapply(leader24_fpkm, function(x){sum(as.numeric(x$dtcr)) > 0})]
cds24_fpkm <- cds24[names(cds24) %in% fpkm24]
cds24_fpkm <- cds24_fpkm[sapply(cds24_fpkm, function(x){length(x) > 100})]
cds24_fpkm <- cds24_fpkm[sapply(cds24_fpkm, function(x){sum(as.numeric(x$dtcr)) > 0})]
trailer24_fpkm <- trailer24[names(trailer24) %in% fpkm24]
trailer24_fpkm <- trailer24_fpkm[sapply(trailer24_fpkm, function(x){length(x) > 30})]
trailer24_fpkm <- trailer24_fpkm[sapply(trailer24_fpkm, function(x){sum(as.numeric(x$dtcr)) > 0})]

leader256_fpkm <- leader256[names(leader256) %in% fpkm256]
leader256_fpkm <- leader256_fpkm[sapply(leader256_fpkm, function(x){length(x) > 30})]
leader256_fpkm <- leader256_fpkm[sapply(leader256_fpkm, function(x){sum(as.numeric(x$dtcr)) > 0})]
cds256_fpkm <- cds256[names(cds256) %in% fpkm256]
cds256_fpkm <- cds256_fpkm[sapply(cds256_fpkm, function(x){length(x) > 100})]
cds256_fpkm <- cds256_fpkm[sapply(cds256_fpkm, function(x){sum(as.numeric(x$dtcr)) > 0})]
trailer256_fpkm <- trailer256[names(trailer256) %in% fpkm256]
trailer256_fpkm <- trailer256_fpkm[sapply(trailer256_fpkm, function(x){length(x) > 30})]
trailer256_fpkm <- trailer256_fpkm[sapply(trailer256_fpkm, function(x){sum(as.numeric(x$dtcr)) > 0})]

### subsample some of the transcripts
common <- Reduce(intersect, list(names(leader24_fpkm), names(cds24_fpkm), names(trailer24_fpkm), names(leader256_fpkm), names(cds256_fpkm), names(trailer256_fpkm)))

sam <- sample(common, 1000)

### run algorithm on each, plot together


### !!!!!! merge intervals with the same value (zeros)
# multiply by a BIIIIG number ? (significant digits problem...)
sam24 <- cds24_fpkm[names(cds24_fpkm) %in% sam]
sam24 <- leader24_fpkm[names(leader24_fpkm) %in% sam]
#sam24 <- sam24[sapply(sam24, function(x){length(x) > 0})]

s <- cds24[[1]]$dtcr
sev <- data.table(
  start = head(seq(0,1,(1/length(s))), -1),
  end = tail(seq(0,1,(1/length(s))), -1),
  value = s / sum(s) # values, normalized
)
setkey(sev, start, end)

for (name in tail(names(cds24), -1)) {
  snew <- cds24[[name]]$dtcr
  sevnew <- data.table(
    start = head(seq(0,1,(1/length(snew))), -1),
    end = tail(seq(0,1,(1/length(snew))), -1),
    value = snew / sum(snew) # values, normalized
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

#### PLOTTING
### add a row to over2 for start = 1.0, end = 1.0, value = (as in last row)
p <- ggplot() + geom_step(data = over2, aes(start, value)) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) + xlab("CDS") + ylab("sum(SHAPE)") + ggtitle("256Cell")

ggsave(file = "/Volumes/USELESS/META/SHAPES/whole_cds/sam256.png", plot = p)


#########
#########
#########

normalizeTx <- function(leader, cds, trailer, sam){
  leader <- leader[names(leader) %in% sam]
  cds <- cds[names(cds) %in% sam]
  trailer <- trailer[names(trailer) %in% sam]
  l <- leader[[1]]$dtcr
  c <- cds[[1]]$dtcr
  t <- trailer[[1]]$dtcr
  # leader
  sev <- data.table(
    start = head(seq(0,1,(1/length(l))), -1),
    end = tail(seq(0,1,(1/length(l))), -1),
    value = l / (sum(l) + sum(c) + sum(t))  # values, normalized
  )
  setkey(sev, start, end)
  for (name in tail(names(leader), -1)) {
    lnew <- leader[[name]]$dtcr
    cnew <- cds[[name]]$dtcr
    tnew <- trailer[[name]]$dtcr
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
    lnew <- leader[[name]]$dtcr
    cnew <- cds[[name]]$dtcr
    tnew <- trailer[[name]]$dtcr
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
    lnew <- leader[[name]]$dtcr
    cnew <- cds[[name]]$dtcr
    tnew <- trailer[[name]]$dtcr
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


sam <- sample(common, 1000)

sev24 <- normalizeTx(leader24_fpkm, cds24_fpkm, trailer24_fpkm, sam)

p24 <- ggplot() + geom_step(data = sev24, aes(start, value)) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) + xlab("leader, CDS, trailer") + ylab("sum(SHAPE)") + ggtitle("2-4Cell")
ggsave(file = "/Volumes/USELESS/META/SHAPES/whole_cds/lct24_sam.png", plot = p24)

sev256 <- normalizeTx(leader256_fpkm, cds256_fpkm, trailer256_fpkm, sam)

p256 <- ggplot() + geom_step(data = sev256, aes(start, value)) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) + xlab("leader, CDS, trailer") + ylab("sum(SHAPE)") + ggtitle("256Cell")
ggsave(file = "/Volumes/USELESS/META/SHAPES/whole_cds/lct256_sam.png", plot = p256)

###################

normalizeCDS <- function(cds, sam) {
  cds <- cds[names(cds) %in% sam]
  s <- cds[[1]]$dtcr
  sev <- data.table(
    start = head(seq(0,1,(1/length(s))), -1),
    end = tail(seq(0,1,(1/length(s))), -1),
    value = s / sum(s) # values, normalized
  )
  setkey(sev, start, end)
  for (name in tail(names(cds), -1)) {
    snew <- cds[[name]]$dtcr
    sevnew <- data.table(
      start = head(seq(0,1,(1/length(snew))), -1),
      end = tail(seq(0,1,(1/length(snew))), -1),
      value = snew / sum(snew) # values, normalized
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
  return(sev)
}

sevSH <- normalizeCDS(cds256_fpkm, sam)
pSH <- ggplot() + geom_step(data = sevSH, aes(start, value)) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) + xlab("CDS") + ylab("sum(SHAPE)") + ggtitle("256Cell")
ggsave(file = "/Volumes/USELESS/META/SHAPES/whole_cds/shape_CDS_256_same_as_rac2.png", plot = pSH)
