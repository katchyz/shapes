library(GenomicFeatures)
library(ggplot2)
txdb_can <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/danio_rerio/Danio_rerio.GRCz10.84.gtf", format = "gtf")

fiveUTR <- fiveUTRsByTranscript(txdb_can, use.names=TRUE)
threeUTR <- threeUTRsByTranscript(txdb_can, use.names=TRUE)
cds <- cdsBy(txdb_can, by="tx", use.names=TRUE)

##############
# access slograt values for transcript
# seqnames(seqinfo(shapeseq_norm)) # 310 transcripts
slr <- elementMetadata(shapeseq_norm)$slograt[as.vector(seqnames(shapeseq_norm)) == "ENSDART00000000992.7"]
#slr <- elementMetadata(shapeseq_norm)$slograt[seqnames(seqinfo(shapeseq_norm)) == "ENSDART00000000992.7"] # longer??? slograt values are different

c <- cds[["ENSDART00000000992"]]
f <- fiveUTR[["ENSDART00000000992"]]
t <- threeUTR[["ENSDART00000000992"]]

y <- cds[[ENSDART00000157370.1]]

length(slr) # = 1032 (is one nt shorter than should be!)
# (if f and t exist = are not NULL)
len_f <- sum(width(ranges(f))) # = 173
len_c <- sum(width(ranges(c))) # = 549
len_t <- sum(width(ranges(t))) # = 311

# sum(width(ranges(f))) + sum(width(ranges(c))) + sum(width(ranges(t)))

# names(threeUTR) %in% names(cds)
################

#### only with leaders
meta = rep(0,33)
utr <- 0
### fix leaders shorter than 15 nt!!!
# for each name in 310
for (tr in seqnames(seqinfo(shapeseq_norm))) {
  slr <- elementMetadata(shapeseq_norm)$slograt[as.vector(seqnames(shapeseq_norm)) == tr]
  tr_name <- substr(tr,1,18)
  if (length(fiveUTR[[tr_name]]) > 0) { # if fiveUTR exists (len_f is not NULL)
    len_f <- sum(width(ranges(fiveUTR[[tr_name]])))
    meta[1:15] <- meta[1:15] + slr[(len_f-15):len_f]
    meta[16:33] <- meta[16:33] + slr[(len_f+1):(len_f+18)]
    utr <- utr + 1
  }
}

scale <- c(-15:18)[-16]
df <- data.frame(scale,meta)
ggplot(df, aes(x=scale, y=meta)) + geom_point() + xlab("position from start codon") + ylab("sum(reactivities)") + ggtitle("2-4cell, slograt, 261 transcripts")

ggsave(file = "/Volumes/USELESS/META/SHAPES/NEW/2-4cell_slograt_subset.png")



# Kornel's package ORFik
#cdsStarts <- extract_START_sites_from_CDSs(cds)
#cdsEnds <- extract_STOP_sites_from_CDSs(cds)


meta = rep(0,33)
utr <- 0
# for each name in 310
for (tr in seqnames(seqinfo(shapeseq_norm))) {
  slr <- elementMetadata(shapeseq_norm)$dtcr[as.vector(seqnames(shapeseq_norm)) == tr]
  tr_name <- substr(tr,1,18)
  if (length(fiveUTR[[tr_name]]) > 0) { # if fiveUTR exists (len_f is not NULL)
    len_f <- sum(width(ranges(fiveUTR[[tr_name]])))
    print(len_f)
    meta[1:15] <- meta[1:15] + slr[(len_f-15):len_f]
    meta[16:33] <- meta[16:33] + slr[(len_f+1):(len_f+18)]
    utr <- utr + 1
  }
}


scale <- c(-15:18)[-16]
df <- data.frame(scale,meta)
ggplot(df, aes(x=scale, y=meta)) + geom_point() + xlab("position from start codon") + ylab("sum(reactivities)") + ggtitle("2-4cell, dtcr, 261 transcripts")

ggsave(file = "/Volumes/USELESS/META/SHAPES/2-4cell_dtcr_subset.png")

