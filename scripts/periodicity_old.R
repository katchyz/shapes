### periodicity in CDSs
library(GenomicFeatures)
library(ggplot2)
library(GeneCycle)
txdb_can <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/danio_rerio/Danio_rerio.GRCz10.84.gtf", format = "gtf")
cds <- cdsBy(txdb_can, by="tx", use.names=TRUE)
fiveUTR <- fiveUTRsByTranscript(txdb_can, use.names=TRUE)

#trnames <- seqnames(seqinfo(dtcr_2_4))
#trnames <- sapply(trnames, function(x){substr(x,1,18)})
### seqnames(seqinfo(dtcr_2_4)) <- sapply(seqnames(seqinfo(dtcr_2_4)), function(x){substr(x,1,18)})
#exons[names(exons) %in% trnames]

scale <- c(1:150)
## 2-4cell, dtcr, normalized to 1
meta = rep(0,150)
for (tr in seqnames(seqinfo(dtcr_2_4))) {
  slr <- elementMetadata(dtcr_2_4)$dtcr[as.vector(seqnames(dtcr_2_4)) == tr]
  tr_name <- substr(tr,1,18)
  if (length(fiveUTR[[tr_name]]) > 0) { # if fiveUTR exists (len_f is not NULL)
    len_f <- sum(width(ranges(fiveUTR[[tr_name]])))
    meta_tr <- slr[(len_f+1):(len_f+150)]
    if (sum(meta_tr, na.rm=TRUE) > 0) {
      meta_tr <- meta_tr/sum(meta_tr)
      meta <- meta + meta_tr
    }
  } else {
    meta_tr <- slr[1:150]
    if (sum(meta_tr, na.rm=TRUE) > 0) {
      meta_tr <- meta_tr/sum(meta_tr)
      meta <- meta + meta_tr
    }
  }
}

# periodicity
amplitudes <- abs(fft(meta))
amp <- amplitudes[2:(length(amplitudes)/2+1)]
periods <- 1/periodogram(meta, method = "clone")$freq

df <- data.frame(periods, amp)
p1 <- ggplot(df, aes(x=periods, y=amp)) + geom_line() + ggtitle("FFT, 2-4cell") + xlim(0, 10)


## 256cell, dtcr
meta = rep(0,50)
for (tr in seqnames(seqinfo(dtcr_256))) {
  slr <- elementMetadata(dtcr_256)$dtcr[as.vector(seqnames(dtcr_256)) == tr]
  tr_name <- substr(tr,1,18)
  if (length(fiveUTR[[tr_name]]) > 0) { # if fiveUTR exists (len_f is not NULL)
    len_f <- sum(width(ranges(fiveUTR[[tr_name]])))
    meta_tr <- slr[(len_f+1):(len_f+150)]
    if (sum(meta_tr, na.rm=TRUE) > 0) {
      meta_tr <- meta_tr/sum(meta_tr)
      meta <- meta + meta_tr
    }
  } else {
    meta_tr <- slr[1:150]
    if (sum(meta_tr, na.rm=TRUE) > 0) {
      meta_tr <- meta_tr/sum(meta_tr)
      meta <- meta + meta_tr
    }
  }
}

# periodicity
amplitudes <- abs(fft(meta))
amp <- amplitudes[2:(length(amplitudes)/2+1)]
periods <- 1/periodogram(meta, method = "clone")$freq

df <- data.frame(periods, amp)
p2 <- ggplot(df, aes(x=periods, y=amp)) + geom_line() + ggtitle("FFT, 256cell") + xlim(0, 10)


ggsave(file = "/Volumes/USELESS/META/SHAPES/periodicity/p1.png", plot = p1)
ggsave(file = "/Volumes/USELESS/META/SHAPES/periodicity/p2.png", plot = p2)

#convert p1.png p2.png -append periodicity.png
