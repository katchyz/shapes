# stop codon
scale <- c(-30:33)[-31]

## 2-4cell, slograt
meta = rep(0,63)
### fix leaders shorter than 15 nt!!!
for (tr in seqnames(seqinfo(slograt_2_4))) {
  slr <- elementMetadata(slograt_2_4)$slograt[as.vector(seqnames(slograt_2_4)) == tr]
  tr_name <- substr(tr,1,18)
  if (length(fiveUTR[[tr_name]]) > 0) {
    if (length(threeUTR[[tr_name]]) > 0) { # if fiveUTR AND threeUTR exist
      len_f <- sum(width(ranges(fiveUTR[[tr_name]])))
      len_c <- sum(width(ranges(cds[[tr_name]])))
      len_t <- sum(width(ranges(threeUTR[[tr_name]])))
      meta_tr = rep(0,63)
      meta_tr[1:33] <- slr[(len_f+len_c-33):(len_f+len_c)]
      meta_tr[34:63] <- slr[(len_f+len_c+1):(len_f+len_c+30)]
      if (sum(meta_tr) > 0) {
        meta_tr <- meta_tr/sum(meta_tr)
        meta <- meta + meta_tr
      } 
    }
  } else {
    if (length(threeUTR[[tr_name]]) > 0) { # no fiveUTR AND three UTR
      len_c <- sum(width(ranges(cds[[tr_name]])))
      len_t <- sum(width(ranges(threeUTR[[tr_name]])))
      meta_tr = rep(0,63)
      meta_tr[1:33] <- slr[(len_c-33):(len_c)]
      meta_tr[34:63] <- slr[(len_c+1):(len_c+30)]
      if (sum(meta_tr) > 0) {
        meta_tr <- meta_tr/sum(meta_tr)
        meta <- meta + meta_tr
      }
    }
  }
}

df <- data.frame(scale,meta)
p1 <- ggplot(df, aes(x=scale, y=meta)) + geom_point() + xlab("position from stop codon") + ylab("sum(reactivities)") + ggtitle("2-4cell, slograt")

## 2-4cell, dtcr
meta = rep(0,63)
### fix leaders shorter than 15 nt!!!
for (tr in seqnames(seqinfo(dtcr_2_4))) {
  slr <- elementMetadata(dtcr_2_4)$dtcr[as.vector(seqnames(dtcr_2_4)) == tr]
  tr_name <- substr(tr,1,18)
  if (length(fiveUTR[[tr_name]]) > 0) {
    if (length(threeUTR[[tr_name]]) > 0) { # if fiveUTR AND threeUTR exist
      len_f <- sum(width(ranges(fiveUTR[[tr_name]])))
      len_c <- sum(width(ranges(cds[[tr_name]])))
      len_t <- sum(width(ranges(threeUTR[[tr_name]])))
      meta_tr = rep(0,63)
      meta_tr[1:33] <- slr[(len_f+len_c-33):(len_f+len_c)]
      meta_tr[34:63] <- slr[(len_f+len_c+1):(len_f+len_c+30)]
      if (sum(meta_tr) > 0) {
        meta_tr <- meta_tr/sum(meta_tr)
        meta <- meta + meta_tr
      } 
    }
  } else {
    if (length(threeUTR[[tr_name]]) > 0) { # no fiveUTR AND three UTR
      len_c <- sum(width(ranges(cds[[tr_name]])))
      len_t <- sum(width(ranges(threeUTR[[tr_name]])))
      meta_tr = rep(0,63)
      meta_tr[1:33] <- slr[(len_c-33):(len_c)]
      meta_tr[34:63] <- slr[(len_c+1):(len_c+30)]
      if (sum(meta_tr) > 0) {
        meta_tr <- meta_tr/sum(meta_tr)
        meta <- meta + meta_tr
      }
    }
  }
}

df <- data.frame(scale,meta)
p2 <- ggplot(df, aes(x=scale, y=meta)) + geom_point() + xlab("position from stop codon") + ylab("sum(reactivities)") + ggtitle("2-4cell, dtcr")

## 256cell, slograt
meta = rep(0,63)
### fix leaders shorter than 15 nt!!!
for (tr in seqnames(seqinfo(slograt_256))) {
  slr <- elementMetadata(slograt_256)$slograt[as.vector(seqnames(slograt_256)) == tr]
  tr_name <- substr(tr,1,18)
  if (length(fiveUTR[[tr_name]]) > 0) {
    if (length(threeUTR[[tr_name]]) > 0) { # if fiveUTR AND threeUTR exist
      len_f <- sum(width(ranges(fiveUTR[[tr_name]])))
      len_c <- sum(width(ranges(cds[[tr_name]])))
      len_t <- sum(width(ranges(threeUTR[[tr_name]])))
      meta_tr = rep(0,63)
      meta_tr[1:33] <- slr[(len_f+len_c-33):(len_f+len_c)]
      meta_tr[34:63] <- slr[(len_f+len_c+1):(len_f+len_c+30)]
      if (sum(meta_tr) > 0) {
        meta_tr <- meta_tr/sum(meta_tr)
        meta <- meta + meta_tr
      } 
    }
  } else {
    if (length(threeUTR[[tr_name]]) > 0) { # no fiveUTR AND three UTR
      len_c <- sum(width(ranges(cds[[tr_name]])))
      len_t <- sum(width(ranges(threeUTR[[tr_name]])))
      meta_tr = rep(0,63)
      meta_tr[1:33] <- slr[(len_c-33):(len_c)]
      meta_tr[34:63] <- slr[(len_c+1):(len_c+30)]
      if (sum(meta_tr) > 0) {
        meta_tr <- meta_tr/sum(meta_tr)
        meta <- meta + meta_tr
      }
    }
  }
}

df <- data.frame(scale,meta)
p3 <- ggplot(df, aes(x=scale, y=meta)) + geom_point() + xlab("position from stop codon") + ylab("sum(reactivities)") + ggtitle("256cell, slograt")

## 256cell, dtcr
meta = rep(0,63)
### fix leaders shorter than 15 nt!!!
for (tr in seqnames(seqinfo(dtcr_256))) {
  slr <- elementMetadata(dtcr_256)$dtcr[as.vector(seqnames(dtcr_256)) == tr]
  tr_name <- substr(tr,1,18)
  if (length(fiveUTR[[tr_name]]) > 0) {
    if (length(threeUTR[[tr_name]]) > 0) { # if fiveUTR AND threeUTR exist
      len_f <- sum(width(ranges(fiveUTR[[tr_name]])))
      len_c <- sum(width(ranges(cds[[tr_name]])))
      len_t <- sum(width(ranges(threeUTR[[tr_name]])))
      meta_tr = rep(0,63)
      meta_tr[1:33] <- slr[(len_f+len_c-33):(len_f+len_c)]
      meta_tr[34:63] <- slr[(len_f+len_c+1):(len_f+len_c+30)]
      if (sum(meta_tr) > 0) {
        meta_tr <- meta_tr/sum(meta_tr)
        meta <- meta + meta_tr
      } 
    }
  } else {
    if (length(threeUTR[[tr_name]]) > 0) { # no fiveUTR AND three UTR
      len_c <- sum(width(ranges(cds[[tr_name]])))
      len_t <- sum(width(ranges(threeUTR[[tr_name]])))
      meta_tr = rep(0,63)
      meta_tr[1:33] <- slr[(len_c-33):(len_c)]
      meta_tr[34:63] <- slr[(len_c+1):(len_c+30)]
      if (sum(meta_tr) > 0) {
        meta_tr <- meta_tr/sum(meta_tr)
        meta <- meta + meta_tr
      }
    }
  }
}

df <- data.frame(scale,meta)
p4 <- ggplot(df, aes(x=scale, y=meta)) + geom_point() + xlab("position from stop codon") + ylab("sum(reactivities)") + ggtitle("256cell, dtcr")


multiplot(p1, p2, p3, p4, cols=2)

ggsave(file = "/Volumes/USELESS/META/SHAPES/stop/p1norm.png", plot = p1)
ggsave(file = "/Volumes/USELESS/META/SHAPES/stop/p2norm.png", plot = p2)
ggsave(file = "/Volumes/USELESS/META/SHAPES/stop/p3norm.png", plot = p3)
ggsave(file = "/Volumes/USELESS/META/SHAPES/stop/p4norm.png", plot = p4)

#convert p1norm.png p3norm.png -append slograt.png
#convert p2norm.png p4norm.png -append dtcr.png
#convert slograt.png dtcr.png +append normalized.png
