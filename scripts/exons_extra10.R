# extend exons 10nt upstream
library(GenomicFeatures)

txdb_can <- makeTxDbFromGFF("/Home/ii/katchyz/DATA/genomes/GTF/Danio_rerio.GRCz10.84.gtf", format = "gtf")
exons <- exonsBy(txdb_can, by="tx", use.names=TRUE)

### extend transcripts by 10nt upstream
# for plus: start of first exon -10
# for minus: end of first exon +10


for (i in 1:length(exons)) {
  if (as.character(strand(exons[[i]]))[1] == "+") {
    e1df <- data.frame(exons[[i]][1])
    e1df$start <- e1df$start - 10
    # end stays the same
    e1df$width <- e1df$width + 10
    exons[[i]][1] <- GRanges(e1df)
  } else {
    e1df <- data.frame(exons[[i]][1])
    e1df$width <- e1df$width + 10
    e1df$end <- e1df$end + 10
    exons[[i]][1] <- GRanges(e1df)
  }
}

save(exons, file = "/Home/ii/katchyz/DATA/genomes/GTF/exons_10extra.Rsave")

