### 3'ends of exons


scale <- c(-40:-1)
## 2-4cell, dtcr
meta = rep(0,40)
for (tr in seqnames(seqinfo(dtcr_2_4))) {
  slr <- elementMetadata(dtcr_2_4)$dtcr[as.vector(seqnames(dtcr_2_4)) == tr]
  tr_name <- substr(tr,1,18)
  exs <- c(0, width(ranges(exons[[tr_name]])))
  len <- 0
  for (w in exs[2:length(exs)]) {
    len <- len + w
    if (w > 40) {
      meta_tr <- slr[(len-39):len]
      if (sum(meta_tr, na.rm=TRUE) > 0) {
        if (!any(is.na(meta_tr))) {
          meta_tr <- meta_tr/sum(meta_tr)
          meta <- meta + meta_tr
        }
      }
    }
  }
}

df <- data.frame(scale,meta)
p1 <- ggplot(df, aes(x=scale, y=meta)) + geom_point() + xlab("position from 3'end of exon") + ylab("sum(reactivities)") + ggtitle("2-4cell, dtcr")


## 256cell, dtcr
meta = rep(0,40)
for (tr in seqnames(seqinfo(dtcr_256))) {
  slr <- elementMetadata(dtcr_256)$dtcr[as.vector(seqnames(dtcr_256)) == tr]
  tr_name <- substr(tr,1,18)
  exs <- c(0, width(ranges(exons[[tr_name]])))
  len <- 0
  for (w in exs[2:length(exs)]) {
    len <- len + w
    if (w > 40) {
      meta_tr <- slr[(len-39):len]
      if (sum(meta_tr, na.rm=TRUE) > 0) {
        if (!any(is.na(meta_tr))) {
          meta_tr <- meta_tr/sum(meta_tr)
          meta <- meta + meta_tr
        }
      }
    }
  }
}

df <- data.frame(scale,meta)
p2 <- ggplot(df, aes(x=scale, y=meta)) + geom_point() + xlab("position from 3'end of exon") + ylab("sum(reactivities)") + ggtitle("256cell, dtcr")



ggsave(file = "/Volumes/USELESS/META/SHAPES/exons/three/p1.png", plot = p1)
ggsave(file = "/Volumes/USELESS/META/SHAPES/exons/three/p2.png", plot = p2)

#convert p1.png p2.png -append exons_ends.png
