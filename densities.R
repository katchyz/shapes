### leaders, cds, trailers
# (normalized by length)
leaders <- c()
cdss <- c()
trailers <- c()

## 2-4cell, slograt
for (tr in seqnames(seqinfo(slograt_2_4))) {
  slr <- elementMetadata(slograt_2_4)$slograt[as.vector(seqnames(slograt_2_4)) == tr]
  tr_name <- substr(tr,1,18)
  if (length(fiveUTR[[tr_name]]) > 0) {
    if (length(threeUTR[[tr_name]]) > 0) { # if fiveUTR AND threeUTR exist
      len_f <- sum(width(ranges(fiveUTR[[tr_name]])))
      len_c <- sum(width(ranges(cds[[tr_name]])))
      len_t <- sum(width(ranges(threeUTR[[tr_name]])))
      leader <- slr[1:len_f]
      coding <- slr[(len_f+1):(len_f+len_c)]
      trailer <- slr[(len_f+len_c+1):length(slr)]
      leaders <- c(leaders, sum(leader)/len_f)
      cdss <- c(cdss, sum(coding)/len_c)
      trailers <- c(trailers, sum(trailer)/len_t)
    }
  }
}

df <- data.frame(leaders, cdss, trailers)
df$scale <- c(1:247)
d <- melt(df, id.var="scale")
p1 <- ggplot(d, aes(x=value, fill=variable)) + geom_density(alpha=.3) + ggtitle("2-4cell, slograt")


leaders <- c()
cdss <- c()
trailers <- c()

## 2-4cell, dtcr
for (tr in seqnames(seqinfo(dtcr_2_4))) {
  slr <- elementMetadata(dtcr_2_4)$dtcr[as.vector(seqnames(dtcr_2_4)) == tr]
  tr_name <- substr(tr,1,18)
  if (length(fiveUTR[[tr_name]]) > 0) {
    if (length(threeUTR[[tr_name]]) > 0) { # if fiveUTR AND threeUTR exist
      len_f <- sum(width(ranges(fiveUTR[[tr_name]])))
      len_c <- sum(width(ranges(cds[[tr_name]])))
      len_t <- sum(width(ranges(threeUTR[[tr_name]])))
      leader <- slr[1:len_f]
      coding <- slr[(len_f+1):(len_f+len_c)]
      trailer <- slr[(len_f+len_c+1):length(slr)]
      leaders <- c(leaders, sum(leader)/len_f)
      cdss <- c(cdss, sum(coding)/len_c)
      trailers <- c(trailers, sum(trailer)/len_t)
    }
  }
}

df <- data.frame(leaders, cdss, trailers)
df$scale <- c(1:247)
d <- melt(df, id.var="scale")
p2 <- ggplot(d, aes(x=value, fill=variable)) + geom_density(alpha=.3) + ggtitle("2-4cell, dtcr") + xlab("average reactivity")

leaders <- c()
cdss <- c()
trailers <- c()

## 256cell, slograt
for (tr in seqnames(seqinfo(slograt_256))) {
  slr <- elementMetadata(slograt_256)$slograt[as.vector(seqnames(slograt_256)) == tr]
  tr_name <- substr(tr,1,18)
  if (length(fiveUTR[[tr_name]]) > 0) {
    if (length(threeUTR[[tr_name]]) > 0) { # if fiveUTR AND threeUTR exist
      len_f <- sum(width(ranges(fiveUTR[[tr_name]])))
      len_c <- sum(width(ranges(cds[[tr_name]])))
      len_t <- sum(width(ranges(threeUTR[[tr_name]])))
      leader <- slr[1:len_f]
      coding <- slr[(len_f+1):(len_f+len_c)]
      trailer <- slr[(len_f+len_c+1):length(slr)]
      leaders <- c(leaders, sum(leader)/len_f)
      cdss <- c(cdss, sum(coding)/len_c)
      trailers <- c(trailers, sum(trailer)/len_t)
    }
  }
}

df <- data.frame(leaders, cdss, trailers)
df$scale <- c(1:247)
d <- melt(df, id.var="scale")
p3 <- ggplot(d, aes(x=value, fill=variable)) + geom_density(alpha=.3) + ggtitle("256cell, slograt")

leaders <- c()
cdss <- c()
trailers <- c()

## 256cell, dtcr
for (tr in seqnames(seqinfo(dtcr_256))) {
  slr <- elementMetadata(dtcr_256)$dtcr[as.vector(seqnames(dtcr_256)) == tr]
  tr_name <- substr(tr,1,18)
  if (length(fiveUTR[[tr_name]]) > 0) {
    if (length(threeUTR[[tr_name]]) > 0) { # if fiveUTR AND threeUTR exist
      len_f <- sum(width(ranges(fiveUTR[[tr_name]])))
      len_c <- sum(width(ranges(cds[[tr_name]])))
      len_t <- sum(width(ranges(threeUTR[[tr_name]])))
      leader <- slr[1:len_f]
      coding <- slr[(len_f+1):(len_f+len_c)]
      trailer <- slr[(len_f+len_c+1):length(slr)]
      leaders <- c(leaders, sum(leader)/len_f)
      cdss <- c(cdss, sum(coding)/len_c)
      trailers <- c(trailers, sum(trailer)/len_t)
    }
  }
}

df <- data.frame(leaders, cdss, trailers)
df$scale <- c(1:247)
d <- melt(df, id.var="scale")
p4 <- ggplot(d, aes(x=value, fill=variable)) + geom_density(alpha=.3) + ggtitle("256cell, dtcr") + xlab("average reactivity")


multiplot(p1, p2, p3, p4, cols=2)

ggsave(file = "/Volumes/USELESS/META/SHAPES/distributions/p1norm.png", plot = p1)
ggsave(file = "/Volumes/USELESS/META/SHAPES/distributions/p2norm.png", plot = p2)
ggsave(file = "/Volumes/USELESS/META/SHAPES/distributions/p3norm.png", plot = p3)
ggsave(file = "/Volumes/USELESS/META/SHAPES/distributions/p4norm.png", plot = p4)

#convert p1norm.png p3norm.png -append slograt.png
#convert p2norm.png p4norm.png -append dtcr.png
#convert slograt.png dtcr.png +append normalized.png
