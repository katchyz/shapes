library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)

### get transcript database
txdb_can <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/danio_rerio/Danio_rerio.GRCz10.84.gtf", format = "gtf")
cds <- cdsBy(txdb_can, by="tx", use.names=TRUE)
cds <- cds[order(names(cds))]
fiveUTR <- fiveUTRsByTranscript(txdb_can, use.names=TRUE)
fiveUTR <- fiveUTR[order(names(fiveUTR))]
threeUTR <- threeUTRsByTranscript(txdb_can, use.names=TRUE)
threeUTR <- threeUTR[order(names(threeUTR))]

names(cds) <- sapply(names(cds), function(x){substr(x,1,18)})
names(fiveUTR) <- sapply(names(fiveUTR), function(x){substr(x,1,18)})
names(threeUTR) <- sapply(names(threeUTR), function(x){substr(x,1,18)})

txLengths <- transcriptLengths(txdb_can, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
rownames(txLengths) <- txLengths$tx_name
rownames(txLengths) <- sapply(rownames(txLengths), function(x){substr(x,1,18)})
txLengths <- txLengths[order(rownames(txLengths)),]

### Ribo-Seq

##################################### 2-4Cell #################################################

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

####
# make granges with 0 score of transcript length, one by one
# substitute values

seqlengths(ribo_cds_24) <- txLengths[rownames(txLengths) %in% names(ribo_cds_24),]$cds_len
seqlengths(ribo_utr5_24) <- txLengths[rownames(txLengths) %in% names(ribo_utr5_24),]$utr5_len
seqlengths(ribo_utr3_24) <- txLengths[rownames(txLengths) %in% names(ribo_utr3_24),]$utr3_len

# cds
cov <- coverage(ribo_cds_24)
cov <- cov[order(names(cov))]
plus <- ribo_cds_24[strand(ribo_cds_24) == "+"]
plus <- plus[sapply(plus, function(x){length(x) > 0})]
minus <- ribo_cds_24[strand(ribo_cds_24) == "-"]
minus <- minus[sapply(minus, function(x){length(x) > 0})]
# whole vector
# for - strand:rev()
#as.vector(cov[["ENSDART00000162928"]])
cov_plus <- sapply(cov[names(cov) %in% names(plus)], function(x){as.vector(x)})
cov_minus <- sapply(cov[names(cov) %in% names(minus)], function(x){rev(as.vector(x))})
for (name in names(cov_plus)) {
  cov_plus[[name]][cov_plus[[name]] == 1] <- ribo_cds_24[[name]]$score
}
for (name in names(cov_minus)) {
  cov_minus[[name]][cov_minus[[name]] == 1] <- rev(ribo_cds_24[[name]]$score)
}
cov <- c(cov_plus, cov_minus)
cov <- cov[order(names(cov))]
cov_cds_24 <- cov

# utr5
cov <- coverage(ribo_utr5_24)
cov <- cov[order(names(cov))]
plus <- ribo_utr5_24[strand(ribo_utr5_24) == "+"]
plus <- plus[sapply(plus, function(x){length(x) > 0})]
minus <- ribo_utr5_24[strand(ribo_utr5_24) == "-"]
minus <- minus[sapply(minus, function(x){length(x) > 0})]
# whole vector
# for - strand:rev()
#as.vector(cov[["ENSDART00000162928"]])
cov_plus <- sapply(cov[names(cov) %in% names(plus)], function(x){as.vector(x)})
cov_minus <- sapply(cov[names(cov) %in% names(minus)], function(x){rev(as.vector(x))})
for (name in names(cov_plus)) {
  cov_plus[[name]][cov_plus[[name]] == 1] <- ribo_utr5_24[[name]]$score
}
for (name in names(cov_minus)) {
  cov_minus[[name]][cov_minus[[name]] == 1] <- rev(ribo_utr5_24[[name]]$score)
}
cov <- c(cov_plus, cov_minus)
cov <- cov[order(names(cov))]
cov_utr5_24 <- cov

# utr3
cov <- coverage(ribo_utr3_24)
cov <- cov[order(names(cov))]
plus <- ribo_utr3_24[strand(ribo_utr3_24) == "+"]
plus <- plus[sapply(plus, function(x){length(x) > 0})]
minus <- ribo_utr3_24[strand(ribo_utr3_24) == "-"]
minus <- minus[sapply(minus, function(x){length(x) > 0})]
# whole vector
# for - strand:rev()
#as.vector(cov[["ENSDART00000162928"]])
cov_plus <- sapply(cov[names(cov) %in% names(plus)], function(x){as.vector(x)})
cov_minus <- sapply(cov[names(cov) %in% names(minus)], function(x){rev(as.vector(x))})
for (name in names(cov_plus)) {
  cov_plus[[name]][cov_plus[[name]] == 1] <- ribo_utr3_24[[name]]$score
}
for (name in names(cov_minus)) {
  cov_minus[[name]][cov_minus[[name]] == 1] <- rev(ribo_utr3_24[[name]]$score)
}
cov <- c(cov_plus, cov_minus)
cov <- cov[order(names(cov))]
cov_utr3_24 <- cov


# subset on Ribo-FPKM
load(file = "/Volumes/USELESS/META/SHAPES/FPKM_24.Rdata")
complete24 <- rownames(FPKM_24[complete.cases(FPKM_24),])

##################################### 256Cell #################################################

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

####
# make granges with 0 score of transcript length, one by one
# substitute values

seqlengths(ribo_cds_256) <- txLengths[rownames(txLengths) %in% names(ribo_cds_256),]$cds_len
seqlengths(ribo_utr5_256) <- txLengths[rownames(txLengths) %in% names(ribo_utr5_256),]$utr5_len
seqlengths(ribo_utr3_256) <- txLengths[rownames(txLengths) %in% names(ribo_utr3_256),]$utr3_len

# cds
cov <- coverage(ribo_cds_256)
cov <- cov[order(names(cov))]
plus <- ribo_cds_256[strand(ribo_cds_256) == "+"]
plus <- plus[sapply(plus, function(x){length(x) > 0})]
minus <- ribo_cds_256[strand(ribo_cds_256) == "-"]
minus <- minus[sapply(minus, function(x){length(x) > 0})]
# whole vector
# for - strand:rev()
#as.vector(cov[["ENSDART00000162928"]])
cov_plus <- sapply(cov[names(cov) %in% names(plus)], function(x){as.vector(x)})
cov_minus <- sapply(cov[names(cov) %in% names(minus)], function(x){rev(as.vector(x))})
for (name in names(cov_plus)) {
  cov_plus[[name]][cov_plus[[name]] == 1] <- ribo_cds_256[[name]]$score
}
for (name in names(cov_minus)) {
  cov_minus[[name]][cov_minus[[name]] == 1] <- rev(ribo_cds_256[[name]]$score)
}
cov <- c(cov_plus, cov_minus)
cov <- cov[order(names(cov))]
cov_cds_256 <- cov

# utr5
cov <- coverage(ribo_utr5_256)
cov <- cov[order(names(cov))]
plus <- ribo_utr5_256[strand(ribo_utr5_256) == "+"]
plus <- plus[sapply(plus, function(x){length(x) > 0})]
minus <- ribo_utr5_256[strand(ribo_utr5_256) == "-"]
minus <- minus[sapply(minus, function(x){length(x) > 0})]
# whole vector
# for - strand:rev()
#as.vector(cov[["ENSDART00000162928"]])
cov_plus <- sapply(cov[names(cov) %in% names(plus)], function(x){as.vector(x)})
cov_minus <- sapply(cov[names(cov) %in% names(minus)], function(x){rev(as.vector(x))})
for (name in names(cov_plus)) {
  cov_plus[[name]][cov_plus[[name]] == 1] <- ribo_utr5_256[[name]]$score
}
for (name in names(cov_minus)) {
  cov_minus[[name]][cov_minus[[name]] == 1] <- rev(ribo_utr5_256[[name]]$score)
}
cov <- c(cov_plus, cov_minus)
cov <- cov[order(names(cov))]
cov_utr5_256 <- cov

# utr3
cov <- coverage(ribo_utr3_256)
cov <- cov[order(names(cov))]
plus <- ribo_utr3_256[strand(ribo_utr3_256) == "+"]
plus <- plus[sapply(plus, function(x){length(x) > 0})]
minus <- ribo_utr3_256[strand(ribo_utr3_256) == "-"]
minus <- minus[sapply(minus, function(x){length(x) > 0})]
# whole vector
# for - strand:rev()
#as.vector(cov[["ENSDART00000162928"]])
cov_plus <- sapply(cov[names(cov) %in% names(plus)], function(x){as.vector(x)})
cov_minus <- sapply(cov[names(cov) %in% names(minus)], function(x){rev(as.vector(x))})
for (name in names(cov_plus)) {
  cov_plus[[name]][cov_plus[[name]] == 1] <- ribo_utr3_256[[name]]$score
}
for (name in names(cov_minus)) {
  cov_minus[[name]][cov_minus[[name]] == 1] <- rev(ribo_utr3_256[[name]]$score)
}
cov <- c(cov_plus, cov_minus)
cov <- cov[order(names(cov))]
cov_utr3_256 <- cov


# subset on Ribo-FPKM
load(file = "/Volumes/USELESS/META/SHAPES/FPKM_256.Rdata")
complete256 <- rownames(FPKM_256[complete.cases(FPKM_256),])

# common in 2-4 and 256
common <- Reduce(intersect, list(complete24, complete256))

#####
normalizeTxRibo <- function(leader, cds, trailer, sam){
  leader <- leader[names(leader) %in% sam]
  cds <- cds[names(cds) %in% sam]
  trailer <- trailer[names(trailer) %in% sam]
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

#
sev24 <- normalizeTxRibo(cov_utr5_24, cov_cds_24, cov_utr3_24, common)

p24 <- ggplot() + geom_step(data = sev24, aes(start, value)) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) + xlab("leader, CDS, trailer") + ylab("sum(Ribo)") + ggtitle("2-4Cell")
ggsave(file = "/Volumes/USELESS/META/SHAPES/whole_cds/lct24_ribo.png", plot = p24)

sev256 <- normalizeTxRibo(cov_utr5_256, cov_cds_256, cov_utr3_256, common)

p256 <- ggplot() + geom_step(data = sev256, aes(start, value)) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) + xlab("leader, CDS, trailer") + ylab("sum(Ribo)") + ggtitle("256Cell")
ggsave(file = "/Volumes/USELESS/META/SHAPES/whole_cds/lct256_ribo.png", plot = p256)

###
# split by footprint lengths (group on short, avg and long)


##### group tx into genes

transcriptsBy(txdb_can, by = "gene")
txdb_can <- makeTxDbFromGFF("/home/ai/Projects/data/genomes/GRCh38/Homo_sapiens.GRCh38.82.chr.gtf")
library(GenomicFeatures)


###############
normalizeCDSRibo <- function(cds, sam){
  #cds <- cds[names(cds) %in% sam]
  c <- cds[[1]]
  # CDS
  sev <- data.table(
    start = head(seq(0,1,(1/length(c))), -1),
    end = tail(seq(0,1,(1/length(c))), -1),
    value = c / sum(c)  # values, normalized
  )
  setkey(sev, start, end)
  for (name in tail(names(cds), -1)) {
    cnew <- cds[[name]]
    sevnew <- data.table(
      start = head(seq(0,1,(1/length(cnew))), -1),
      end = tail(seq(0,1,(1/length(cnew))), -1),
      value = cnew / sum(cnew)  # values, normalized
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

#
sevRIBO <- normalizeCDSRibo(cdsX)

pRIBO <- ggplot() + geom_step(data = sevRIBO, aes(start, value)) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) + xlab("CDS") + ylab("sum(Ribo)") + ggtitle("25Cell, ribo")
ggsave(file = "/Volumes/USELESS/META/SHAPES/whole_cds/cds_ribo_same_as_rca2.png", plot = pRIBO)
