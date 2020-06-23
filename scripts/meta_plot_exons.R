### 5'ends of exons
library(GenomicFeatures)
library(ggplot2)
txdb_can <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/danio_rerio/Danio_rerio.GRCz10.84.gtf", format = "gtf")
exons <- exonsBy(txdb_can, by="tx", use.names=TRUE)

trnames <- seqnames(seqinfo(dtcr_2_4))
trnames <- sapply(trnames, function(x){substr(x,1,18)})
exons[names(exons) %in% trnames]

scale <- c(1:40)
## 2-4cell, dtcr
meta = rep(0,40)
for (tr in seqnames(seqinfo(dtcr_2_4))) {
  slr <- elementMetadata(dtcr_2_4)$dtcr[as.vector(seqnames(dtcr_2_4)) == tr]
  tr_name <- substr(tr,1,18)
  exs <- c(0, width(ranges(exons[[tr_name]])))
  len <- 0
  for (w in exs[2:(length(exs)-1)]) {
    len <- len + w
    meta_tr <- slr[(len+1):(len+40)]
    if (sum(meta_tr, na.rm=TRUE) > 0) {
      meta_tr <- meta_tr/sum(meta_tr)
      meta <- meta + meta_tr
    }
  }
}

df <- data.frame(scale,meta)
p1 <- ggplot(df, aes(x=scale, y=meta)) + geom_point() + xlab("position from 5'end of exon") + ylab("sum(reactivities)") + ggtitle("2-4cell, dtcr")


## 256cell, dtcr
meta = rep(0,50)
for (tr in seqnames(seqinfo(dtcr_256))) {
  slr <- elementMetadata(dtcr_256)$dtcr[as.vector(seqnames(dtcr_256)) == tr]
  tr_name <- substr(tr,1,18)
  exs <- c(0, width(ranges(exons[[tr_name]])))
  len <- 0
  for (w in exs[2:(length(exs)-1)]) {
    len <- len + w
    meta_tr <- slr[(len+1):(len+40)]
    if (sum(meta_tr, na.rm=TRUE) > 0) {
      meta_tr <- meta_tr/sum(meta_tr)
      meta <- meta + meta_tr
    }
  }
}

df <- data.frame(scale,meta)
p2 <- ggplot(df, aes(x=scale, y=meta)) + geom_point() + xlab("position from 5'end of exon") + ylab("sum(reactivities)") + ggtitle("256cell, dtcr")



ggsave(file = "/Volumes/USELESS/META/SHAPES/exons/p1.png", plot = p1)
ggsave(file = "/Volumes/USELESS/META/SHAPES/exons/p2.png", plot = p2)

#convert p1.png p2.png -append exons.png
