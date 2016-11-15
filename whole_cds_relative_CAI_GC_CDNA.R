# whole CDS plotting, relative adaptiveness of codons
library(seqinr)

# relative adaptiveness of a codon
w <- read.csv("/Volumes/USELESS/META/SHAPES/w.csv", header = TRUE, row.names = 1)
# w['AAA',]

# fasta
fasta_cds <- read.fasta("/Volumes/USELESS/DATA/fasta/danio_rerio/GRCz10/cds/Danio_rerio.GRCz10.cds.all.fa")
names(fasta_cds) <- sapply(names(fasta_cds), function(x){substr(x,1,18)})
fasta_cds <- fasta_cds[unlist(lapply(fasta_cds, function(x){length(x) %% 3 == 0}))]

# for each CDS sequence, put values from w in place of codons
f <- fasta_cds[[1]]
tapply(x, rep(1:(length(x)/3), each = 3), paste, collapse = "")

sapply(fasta_cds, function(x){tapply(x, rep(1:(length(x)/3), each = 3), paste, collapse = "")})
# pick the same genes as for shape & ribo wholeTx plots
# normalizeCDS

#sam <- fasta_cds[1:10]
sam <- sample(common, 1000)
samCAI <- fasta_cds[names(fasta_cds) %in% sam]

temp <- sapply(samCAI, function(x){tapply(x, rep(1:(length(x)/3), each = 3), paste, collapse = "")})
codon_cai <- sapply(temp, function(t){as.vector(sapply(t, function(x){w[toupper(x[[1]]),]}))})

#sapply(sam, function(x){split(x, ceiling(seq_along(x)/3))})


###########
###########

normalizeCdsCAI <- function(cds){
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
sev <- normalizeCdsCAI(codon_cai)

p <- ggplot() + geom_step(data = sev, aes(start, value)) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) + xlab("CDS") + ylab("sum(w)") + ggtitle("relative adaptiveness of codons")
ggsave(file = "/Volumes/USELESS/META/SHAPES/whole_cds/relative_adaptiveness_of_codons2.png", plot = p)


##############
### GC content
sam <- sample(common, 1000)
samGC <- fasta_cds[names(fasta_cds) %in% sam]
gc <- data.table(a = 2, t = 2, g = 3, c = 3)

gc_score <- sapply(samGC, function(s){as.vector(sapply(s, function(x){gc[[x]]}))})


normalizeGC <- function(cds){
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
sevGC <- normalizeGC(gc_score)

pGC <- ggplot() + geom_step(data = sevGC, aes(start, value)) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) + xlab("CDS") + ylab("sum(gc)") + ggtitle("GC bias (G,C = 3; A,T = 2)")
ggsave(file = "/Volumes/USELESS/META/SHAPES/whole_cds/GC_bias_CDS.png", plot = pGC)


ix <- 0
for (i in gc3_score) {
  ix <- ix + 1
  print(sum(i))
}

## GC3

gc3_score <- sapply(samGC, function(s){as.vector(sapply(s[seq(3, length(s), 3)], function(x){gc[[x]]}))})

sevGC3 <- normalizeGC(gc3_score)

pGC3 <- ggplot() + geom_step(data = sevGC3, aes(start, value)) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) + xlab("CDS") + ylab("sum(gc3)") + ggtitle("GC bias at 3rd nt in codon (G,C = 3; A,T = 2)")
ggsave(file = "/Volumes/USELESS/META/SHAPES/whole_cds/GC3_bias_CDS.png", plot = pGC3)

pGC3 <- ggplot() + geom_step(data = sevGC3, aes(start, value)) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) + xlab("CDS") + ylab("sum(gc3)") + ggtitle("GC bias at 3rd nt in codon (G,C = 5; A,T = 1)")

###
## GC over tx
# get tx fasta
fasta_cdna <- read.fasta("/Volumes/USELESS/DATA/fasta/danio_rerio/GRCz10/cdna/Danio_rerio.GRCz10.cdna.all.fa")
names(fasta_cdna) <- sapply(names(fasta_cdna), function(x){substr(x,1,18)})

# pick sample
samCDNA <- fasta_cdna[names(fasta_cdna) %in% sam]
gc <- data.table(a = 2, t = 2, g = 3, c = 3)
gc_cdna <- sapply(samCDNA, function(s){as.vector(sapply(s, function(x){gc[[x]]}))})
# divide - length of fiveUTR, cds (txLengths)

normalizeCDNA <- function(cdna, txLengths){
  n <- names(cdna[1])
  l <- cdna[[n]][1:txLengths[n,]$utr5_len]
  c <- cdna[[n]][(txLengths[n,]$utr5_len+1):(txLengths[n,]$utr5_len+txLengths[n,]$cds_len)]
  t <- cdna[[n]][(txLengths[n,]$utr5_len+txLengths[n,]$cds_len+1):length(cdna[[n]])]
  # leader
  sev <- data.table(
    start = head(seq(0,1,(1/length(l))), -1),
    end = tail(seq(0,1,(1/length(l))), -1),
    value = l / (sum(l) + sum(c) + sum(t))  # values, normalized
  )
  setkey(sev, start, end)
  for (name in tail(names(cdna), -1)) {
    lnew <- cdna[[name]][1:txLengths[name,]$utr5_len]
    cnew <- cdna[[name]][(txLengths[name,]$utr5_len+1):(txLengths[name,]$utr5_len+txLengths[name,]$cds_len)]
    tnew <- cdna[[name]][(txLengths[name,]$utr5_len+txLengths[name,]$cds_len+1):length(cdna[[name]])]
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
  for (name in tail(names(cdna), -1)) {
    lnew <- cdna[[name]][1:txLengths[name,]$utr5_len]
    cnew <- cdna[[name]][(txLengths[name,]$utr5_len+1):(txLengths[name,]$utr5_len+txLengths[name,]$cds_len)]
    tnew <- cdna[[name]][(txLengths[name,]$utr5_len+txLengths[name,]$cds_len+1):length(cdna[[name]])]
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
  for (name in tail(names(cdna), -1)) {
    lnew <- cdna[[name]][1:txLengths[name,]$utr5_len]
    cnew <- cdna[[name]][(txLengths[name,]$utr5_len+1):(txLengths[name,]$utr5_len+txLengths[name,]$cds_len)]
    tnew <- cdna[[name]][(txLengths[name,]$utr5_len+txLengths[name,]$cds_len+1):length(cdna[[name]])]
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

sevCDNA <- normalizeCDNA(gc_cdna, txLengths)

pCDNA <- ggplot() + geom_step(data = sevCDNA, aes(start, value)) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) + xlab("leader, CDS, trailer") + ylab("sum(gc)") + ggtitle("GC bias (G,C = 3; A,T = 2)")
ggsave(file = "/Volumes/USELESS/META/SHAPES/whole_cds/GC_bias_CDNA.png", plot = pCDNA)
