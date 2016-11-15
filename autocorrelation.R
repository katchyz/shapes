### AUTOCORRELATION
library(rtracklayer)
library(GenomicFeatures)

txdb_can <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/danio_rerio/zebrafish_GRCz10.gtf", format = "gtf")
#exons <- exonsBy(txdb_can, by="tx", use.names=TRUE)
cds <- cdsBy(txdb_can, by="tx", use.names=TRUE)
names(cds) <- sapply(names(cds), function(x){substr(x,1,18)})

a <- c(117,80,32,190,120,31,108,68,52)
b <- c(201,144,32,175,40,15,89,70,20)

ab <- ccf(a,b)
#ab <- ccf(a,b, lag.max = 3)

ab[3] # autocorrelation at lag of 3

ab$acf # list of autocorrelation values [-lag.max, lag.max]

#####
## Ribo-Seq data
# read Wigs
wig_2_4_fw <- import.wig("/Volumes/USELESS/OUT/zebrafish_our_out/GRCz10/2-4Cell_forward.wig")
wig_2_4_rv <- import.wig("/Volumes/USELESS/OUT/zebrafish_our_out/GRCz10/2-4Cell_reverse.wig")
wig_256_fw <- import.wig("/Volumes/USELESS/OUT/zebrafish_our_out/GRCz10/256Cell_forward.wig")
wig_256_rv <- import.wig("/Volumes/USELESS/OUT/zebrafish_our_out/GRCz10/256Cell_reverse.wig")

# substitute strand, merge
strand(wig_2_4_fw) <- rep("+", length(wig_2_4_fw))
strand(wig_2_4_rv) <- rep("-", length(wig_2_4_rv))
wig_2_4 <- c(wig_2_4_fw, wig_2_4_rv)

strand(wig_256_fw) <- rep("+", length(wig_256_fw))
strand(wig_256_rv) <- rep("-", length(wig_256_rv))
wig_256 <- c(wig_256_fw, wig_256_rv)

# map to transcripts, subset CDSs
ribo_2_4 <- mapToTranscripts(wig_2_4, cds, ignore.strand = FALSE)
mcols(ribo_2_4) <- cbind(mcols(ribo_2_4), DataFrame(wig_2_4[ribo_2_4$xHits]))
ribo_2_4 <- split(ribo_2_4, seqnames(ribo_2_4))

ribo_256 <- mapToTranscripts(wig_256, cds, ignore.strand = FALSE)
mcols(ribo_256) <- cbind(mcols(ribo_256), DataFrame(wig_256[ribo_256$xHits]))
ribo_256 <- split(ribo_256, seqnames(ribo_256))

# enumerate ribo-seq 1 by 1 ?

sum(width(ranges(cds))) # CDS lengths

sum(width(ranges(cds[["ENSDART00000124088"]])))

riboseq_2_4 <- list()
for (name in names(ribo_2_4)) {
  tx <- rep(0, sum(width(ranges(cds[[name]]))))
  for (i in 1:length(ranges(ribo_2_4[[name]]))) {
    tx[start(ranges(ribo_2_4[[name]])[i]):start(ranges(ribo_2_4[[name]])[i])] <- score(ribo_2_4[[name]][i])
  }
  riboseq_2_4[[name]] <- tx
}

riboseq_256 <- list()
for (name in names(ribo_256)) {
  tx <- rep(0, sum(width(ranges(cds[[name]]))))
  for (i in 1:length(ranges(ribo_256[[name]]))) {
    tx[start(ranges(ribo_256[[name]])[i]):start(ranges(ribo_256[[name]])[i])] <- score(ribo_256[[name]][i])
  }
  riboseq_256[[name]] <- tx
}


load("/Volumes/USELESS/META/SHAPES/riboseq_2_4.Rdata")
load("/Volumes/USELESS/META/SHAPES/riboseq_256.Rdata")

# cut-off 
we_2_4 <- riboseq_2_4[sapply(riboseq_2_4, function(x){median(unname(tapply(x, (seq_along(x)-1) %/% 3, sum))) > 0})]
we_256 <- riboseq_256[sapply(riboseq_256, function(x){median(unname(tapply(x, (seq_along(x)-1) %/% 3, sum))) > 0})]

#lapply(we_256, function(x){acf(x, lag.max = 3, plot = FALSE)[3]})

cds_we_2_4 <- cds_2_4[names(cds_2_4) %in% names(we_2_4)]
ribo_we_2_4 <- we_2_4[names(we_2_4) %in% names(cds_we_2_4)]
ribo_we_2_4 <- ribo_we_2_4[order(names(ribo_we_2_4))]

cds_we_256 <- cds_256[names(cds_256) %in% names(we_256)]
ribo_we_256 <- we_256[names(we_256) %in% names(cds_we_256)]
ribo_we_256 <- ribo_we_256[order(names(ribo_we_256))]

ccf0_256 <- mapply(function(x,y){ccf(x$dtcr,y, plot = FALSE)[0]$acf[1]}, x = cds_we_256, y = ribo_we_256)
ccf3_256 <- mapply(function(x,y){ccf(x$dtcr,y, plot = FALSE)[3]$acf[1]}, x = cds_we_256, y = ribo_we_256)

ccf3_256[which.max(ccf3_256)]

########
# CROSS-CORRELATION 2-4cell to 256cell
ccf(cds24_fpkm[["ENSDART00000124467"]]$dtcr, cds256_fpkm[["ENSDART00000124467"]]$dtcr)[0]$acf[1]

c1 <- cds24_fpkm[names(cds24_fpkm) %in% names(cds256_fpkm)]
c2 <- cds256_fpkm[names(cds256_fpkm) %in% names(c1)]

# mapply ?
cc0 <- c()
for (name in names(c1)) {
  #cc0 <- c(cc0, name = ccf(c1[[name]]$dtcr, c2[[name]]$dtcr, plot = FALSE)[0]$acf[1])
  cc0 <- c(cc0, setNames(list(ccf(c1[[name]]$dtcr, c2[[name]]$dtcr, plot = FALSE)[0]$acf[1]), paste(name)))
}
temp <- data.frame(unname(unlist(cc0)))
colnames(temp) <- "cc0"
ggplot() + geom_histogram(data = temp, aes(cc0)) + ggtitle("cross-correlation at lag 0 on CDSs (2-4cell & 256cell)")
ggsave(file = "/Volumes/USELESS/META/SHAPES/cross_correlation/cc0_CDS.png")

# cut-off on RNA fpkm

