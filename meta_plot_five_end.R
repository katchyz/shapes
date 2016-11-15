### 5'end of transcript
scale <- c(1:50)

## 2-4cell, slograt
meta = rep(0,50)
for (tr in seqnames(seqinfo(slograt_2_4))) {
  slr <- elementMetadata(slograt_2_4)$slograt[as.vector(seqnames(slograt_2_4)) == tr]
  tr_name <- substr(tr,1,18)
  meta_tr <- slr[0:50]
  if (sum(meta_tr) > 0) {
    meta_tr <- meta_tr/sum(meta_tr)
    meta <- meta + meta_tr
  }
}

df <- data.frame(scale,meta)
p1 <- ggplot(df, aes(x=scale, y=meta)) + geom_point() + xlab("position from 5'end of transcript") + ylab("sum(reactivities)") + ggtitle("2-4cell, slograt")

## 2-4cell, dtcr
meta = rep(0,50)
for (tr in seqnames(seqinfo(dtcr_2_4))) {
  slr <- elementMetadata(dtcr_2_4)$dtcr[as.vector(seqnames(dtcr_2_4)) == tr]
  tr_name <- substr(tr,1,18)
  meta_tr <- slr[0:50]
  if (sum(meta_tr) > 0) {
    meta_tr <- meta_tr/sum(meta_tr)
    meta <- meta + meta_tr
  }
}

df <- data.frame(scale,meta)
p2 <- ggplot(df, aes(x=scale, y=meta)) + geom_point() + xlab("position from 5'end of transcript") + ylab("sum(reactivities)") + ggtitle("2-4cell, dtcr")

## 256cell, slograt
meta = rep(0,50)
for (tr in seqnames(seqinfo(slograt_256))) {
  slr <- elementMetadata(slograt_256)$slograt[as.vector(seqnames(slograt_256)) == tr]
  tr_name <- substr(tr,1,18)
  meta_tr <- slr[0:50]
  if (sum(meta_tr) > 0) {
    meta_tr <- meta_tr/sum(meta_tr)
    meta <- meta + meta_tr
  }
}

df <- data.frame(scale,meta)
p3 <- ggplot(df, aes(x=scale, y=meta)) + geom_point() + xlab("position from 5'end of transcript") + ylab("sum(reactivities)") + ggtitle("256cell, slograt")

## 256cell, dtcr
meta = rep(0,50)
for (tr in seqnames(seqinfo(dtcr_256))) {
  slr <- elementMetadata(dtcr_256)$dtcr[as.vector(seqnames(dtcr_256)) == tr]
  tr_name <- substr(tr,1,18)
  meta_tr <- slr[0:50]
  if (sum(meta_tr) > 0) {
    meta_tr <- meta_tr/sum(meta_tr)
    meta <- meta + meta_tr
  }
}

df <- data.frame(scale,meta)
p4 <- ggplot(df, aes(x=scale, y=meta)) + geom_point() + xlab("position from 5'end of transcript") + ylab("sum(reactivities)") + ggtitle("256cell, dtcr")


multiplot(p1, p2, p3, p4, cols=2)

ggsave(file = "/Volumes/USELESS/META/SHAPES/five_end/p1norm.png", plot = p1)
ggsave(file = "/Volumes/USELESS/META/SHAPES/five_end/p2norm.png", plot = p2)
ggsave(file = "/Volumes/USELESS/META/SHAPES/five_end/p3norm.png", plot = p3)
ggsave(file = "/Volumes/USELESS/META/SHAPES/five_end/p4norm.png", plot = p4)

#convert p1norm.png p3norm.png -append slograt.png
#convert p2norm.png p4norm.png -append dtcr.png
#convert slograt.png dtcr.png +append normalized.png
