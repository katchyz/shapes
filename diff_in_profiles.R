####### difference in profiles
library(reshape2)

sum_subtr <- data.frame(transcript = character(0), value=numeric(0), stringsAsFactors=F)
for (tr in seqnames(seqinfo(dtcr_2_4))) {
  slr2_4 <- elementMetadata(dtcr_2_4)$dtcr[as.vector(seqnames(dtcr_2_4)) == tr]
  slr256 <- elementMetadata(dtcr_256)$dtcr[as.vector(seqnames(dtcr_256)) == tr]
  sum_subtr[nrow(sum_subtr)+1,] <- c(tr, sum(slr2_4-slr256))
}

sum_subtr$value <- as.numeric(sum_subtr$value)

ggplot(sum_subtr, aes(x=value)) + geom_histogram(bins = 50) + xlab("sum(react_2-4cell - react_256cell)")
