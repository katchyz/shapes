###
# plot sum treated vs sum control
c256 <- ctrl_256[names(ctrl_256) %in% names(trea_256)]
c256 <- sapply(c256, function(x){sum(x$TC)})
t24 <- trea_2_4[names(trea_2_4) %in% names(ctrl_2_4)]
t24 <- sapply(t24, function(x){sum(x$TC)})

df <- data.frame(control = as.vector(c24), treated = as.vector(t24))
ggplot(df, aes(x = control, y = treated)) + geom_point() + scale_x_log10() + scale_y_log10() + ggtitle("2-4cell")
ggsave(file = "/Volumes/USELESS/META/SHAPES/general/sum_control_vs_sum_treated24.png")

# sum normalized correlates with sum raw (treated?)
e24 <- exons24[names(exons24) %in% names(trea_2_4)]
e24 <- sapply(e24, function(x){sum(x$dtcr)})
t24 <- trea_2_4[names(trea_2_4) %in% names(exons24)]
t24 <- sapply(t24, function(x){sum(x$TC)})

df <- data.frame(normalized = as.vector(e24), treated = as.vector(t24))
ggplot(df, aes(x = normalized, y = treated)) + geom_point() + scale_x_log10() + scale_y_log10() + ggtitle("2-4cell")
ggsave(file = "/Volumes/USELESS/META/SHAPES/general/sum_normalized_vs_sum_treated24.png")

# exons - calculate shape fpkm, cut-off (on treated)
sum_tx <- sapply(trea_2_4, function(x){sum(x$TC)})
total_sum <- sum(as.vector(sum_tx))
len_tx <- txLengths[rownames(txLengths) %in% names(sum_tx),]
len_tx <- len_tx[order(rownames(len_tx)),]$tx_len

SHAPE_FPKM_24 <- (sum_tx / len_tx) * (10^9 / total_sum)
save(SHAPE_FPKM_24, file = "/Volumes/USELESS/META/SHAPES/SHAPE_FPKM_24.Rdata")

#
sum_tx <- sapply(trea_256, function(x){sum(x$TC)})
total_sum <- sum(as.vector(sum_tx))
len_tx <- txLengths[rownames(txLengths) %in% names(sum_tx),]
len_tx <- len_tx[order(rownames(len_tx)),]$tx_len

SHAPE_FPKM_256 <- (sum_tx / len_tx) * (10^9 / total_sum)
save(SHAPE_FPKM_256, file = "/Volumes/USELESS/META/SHAPES/SHAPE_FPKM_256.Rdata")

###
# SHAPE_fpkm vs Ribo_fpkm
s <- SHAPE_FPKM_24[names(SHAPE_FPKM_24) %in% rownames(FPKM_24)]
r <- FPKM_24[rownames(FPKM_24) %in% names(SHAPE_FPKM_24),]
df <- data.frame(shape_fpkm = as.vector(s), rna_fpkm = r$exons_RNA_fpkm)
ggplot(df, aes(x = shape_fpkm, y = rna_fpkm)) + geom_point() + scale_x_log10() + scale_y_log10() + ggtitle("2-4cell")
ggsave(file = "/Volumes/USELESS/META/SHAPES/general/shape_fpkm_vs_rna_fpkm24.png")

# get long enough (utr5 > 30, cds > 30...)
lc <- rownames(txLengths[(txLengths$utr5_len > 30) & (txLengths$cds_len > 50), ])
# cut-off on shape fpkm
s24 <- SHAPE_FPKM_24[SHAPE_FPKM_24 > 1]
s24 <- names(s24[names(s24) %in% lc])

s256 <- SHAPE_FPKM_256[SHAPE_FPKM_256 > 1]
s256 <- names(s256[names(s256) %in% lc])

# txLengths - get fragment around start codon
scale <- c(-30:33)[-31]
## 2-4cell
meta = rep(0,63)
for (tr in s24) {
  if (sum(exons24[[tr]]$dtcr) > 0 & (txLengths[tr,]$utr5_len+33) < length(exons24[[tr]])) {
    #print(exons24[[tr]]$dtcr[(txLengths[tr,]$utr5_len-14):(txLengths[tr,]$utr5_len+18)])
    meta_tr <- exons24[[tr]]$dtcr[(txLengths[tr,]$utr5_len-29):(txLengths[tr,]$utr5_len+33)] / sum(exons24[[tr]]$dtcr)
    if (length(meta_tr) == 63) {
      meta <- meta + meta_tr
      print(tr)
      print(meta)
    }
  }
}

# exons24[["ENSDART00000083470"]] <- NULL

# plot
df <- data.frame(scale,meta)
ggplot(df, aes(x=scale, y=meta)) + geom_bar(stat = "identity") + xlab("position from start codon") + ylab("sum(reactivities)") + ggtitle("2-4cell")
ggsave(file = "/Volumes/USELESS/META/SHAPES/start/start24.png")

## 256cell
meta = rep(0,63)
for (tr in s256) {
  if (sum(exons256[[tr]]$dtcr) > 0 & (txLengths[tr,]$utr5_len+33) < length(exons256[[tr]])) {
    meta_tr <- exons256[[tr]]$dtcr[(txLengths[tr,]$utr5_len-29):(txLengths[tr,]$utr5_len+33)] / sum(exons256[[tr]]$dtcr)
    if (length(meta_tr) == 63) {
      meta <- meta + meta_tr
    }
  }
}

# plot
df <- data.frame(scale,meta)
ggplot(df, aes(x=scale, y=meta)) + geom_bar(stat = "identity") + xlab("position from start codon") + ylab("sum(reactivities)") + ggtitle("256cell")
ggsave(file = "/Volumes/USELESS/META/SHAPES/start/start256.png")

##### common
common <- Reduce(intersect, list(s24,s256))

##
## 2-4cell, common
meta = rep(0,63)
for (tr in common) {
  if (sum(exons24[[tr]]$dtcr) > 0 & (txLengths[tr,]$utr5_len+33) < length(exons24[[tr]])) {
    meta_tr <- exons24[[tr]]$dtcr[(txLengths[tr,]$utr5_len-29):(txLengths[tr,]$utr5_len+33)] / sum(exons24[[tr]]$dtcr)
    if (length(meta_tr) == 63) {
      meta <- meta + meta_tr
    }
  }
}
# plot
df <- data.frame(scale,meta)
ggplot(df, aes(x=scale, y=meta)) + geom_bar(stat = "identity") + xlab("position from start codon") + ylab("sum(reactivities)") + ggtitle("2-4cell")
ggsave(file = "/Volumes/USELESS/META/SHAPES/start/start_common24.png")

## 256cell, common
meta = rep(0,63)
for (tr in common) {
  if (sum(exons256[[tr]]$dtcr) > 0 & (txLengths[tr,]$utr5_len+33) < length(exons256[[tr]])) {
    meta_tr <- exons256[[tr]]$dtcr[(txLengths[tr,]$utr5_len-29):(txLengths[tr,]$utr5_len+33)] / sum(exons256[[tr]]$dtcr)
    if (length(meta_tr) == 63) {
      meta <- meta + meta_tr
    }
  }
}
# plot
df <- data.frame(scale,meta)
ggplot(df, aes(x=scale, y=meta)) + geom_bar(stat = "identity") + xlab("position from start codon") + ylab("sum(reactivities)") + ggtitle("256cell")
ggsave(file = "/Volumes/USELESS/META/SHAPES/start/start_common256.png")

###
SHAPE_FPKM <- t(do.call(rbind.fill, list(data.frame(t(SHAPE_FPKM_24)),data.frame(t(SHAPE_FPKM_256)))))
SHAPE_FPKM <- data.frame(SHAPE_FPKM)
colnames(SHAPE_FPKM) <- c("stage_2_4cell", "stage_256cell")
save(SHAPE_FPKM, file = "/Volumes/USELESS/META/SHAPES/SHAPE_FPKM.Rdata")

# log2 fold change: minus - increase, plus - decrease
SHAPE_FPKM$log2_fold_change <- log2(SHAPE_FPKM$stage_2_4cell / SHAPE_FPKM$stage_256cell)
SHAPE_FPKM <- SHAPE_FPKM[complete.cases(SHAPE_FPKM),]

# log2_fold_change > 1.5 high_low
high_low <- rownames(SHAPE_FPKM[with(SHAPE_FPKM, log2_fold_change > 1.5),])

# log2_fold_change < -1.5 low_high
low_high <- rownames(SHAPE_FPKM[with(SHAPE_FPKM, log2_fold_change < -1.5),])

# -0.5 < log2_fold_change < 0.5 stable
stable <- rownames(SHAPE_FPKM[with(SHAPE_FPKM, log2_fold_change > -0.2 & log2_fold_change < 0.2),])

# (plot each for 2-4 and 256)
##########
## 2-4cell, high_low
meta = rep(0,63)
for (tr in high_low) {
  if (txLengths[tr,]$utr5_len > 30 & sum(exons24[[tr]]$dtcr) > 0 & (txLengths[tr,]$utr5_len+33) < length(exons24[[tr]])) {
    meta_tr <- exons24[[tr]]$dtcr[(txLengths[tr,]$utr5_len-29):(txLengths[tr,]$utr5_len+33)] / sum(exons24[[tr]]$dtcr)
    if (length(meta_tr) == 63) {
      meta <- meta + meta_tr
    }
  }
}
# plot
df <- data.frame(scale,meta)
ggplot(df, aes(x=scale, y=meta)) + geom_bar(stat = "identity") + xlab("position from start codon") + ylab("sum(reactivities)") + ggtitle("2-4cell, high_low")
ggsave(file = "/Volumes/USELESS/META/SHAPES/start/fold_change/high_low24.png")

## 256cell, high_low
meta = rep(0,63)
for (tr in high_low) {
  if (txLengths[tr,]$utr5_len > 30 & sum(exons256[[tr]]$dtcr) > 0 & (txLengths[tr,]$utr5_len+33) < length(exons256[[tr]])) {
    meta_tr <- exons256[[tr]]$dtcr[(txLengths[tr,]$utr5_len-29):(txLengths[tr,]$utr5_len+33)] / sum(exons256[[tr]]$dtcr)
    if (length(meta_tr) == 63) {
      meta <- meta + meta_tr
    }
  }
}
# plot
df <- data.frame(scale,meta)
ggplot(df, aes(x=scale, y=meta)) + geom_bar(stat = "identity") + xlab("position from start codon") + ylab("sum(reactivities)") + ggtitle("256cell, high_low")
ggsave(file = "/Volumes/USELESS/META/SHAPES/start/fold_change/high_low256.png")

#########
## 2-4cell, low_high
meta = rep(0,63)
for (tr in low_high) {
  if (txLengths[tr,]$utr5_len > 30 & sum(exons24[[tr]]$dtcr) > 0 & (txLengths[tr,]$utr5_len+33) < length(exons24[[tr]])) {
    meta_tr <- exons24[[tr]]$dtcr[(txLengths[tr,]$utr5_len-29):(txLengths[tr,]$utr5_len+33)] / sum(exons24[[tr]]$dtcr)
    if (length(meta_tr) == 63) {
      meta <- meta + meta_tr
    }
  }
}
# plot
df <- data.frame(scale,meta)
ggplot(df, aes(x=scale, y=meta)) + geom_bar(stat = "identity") + xlab("position from start codon") + ylab("sum(reactivities)") + ggtitle("2-4cell, low_high")
ggsave(file = "/Volumes/USELESS/META/SHAPES/start/fold_change/low_high24.png")

## 256cell, low_high
meta = rep(0,63)
for (tr in low_high) {
  if (txLengths[tr,]$utr5_len > 30 & sum(exons256[[tr]]$dtcr) > 0 & (txLengths[tr,]$utr5_len+33) < length(exons256[[tr]])) {
    meta_tr <- exons256[[tr]]$dtcr[(txLengths[tr,]$utr5_len-29):(txLengths[tr,]$utr5_len+33)] / sum(exons256[[tr]]$dtcr)
    if (length(meta_tr) == 63) {
      meta <- meta + meta_tr
    }
  }
}
# plot
df <- data.frame(scale,meta)
ggplot(df, aes(x=scale, y=meta)) + geom_bar(stat = "identity") + xlab("position from start codon") + ylab("sum(reactivities)") + ggtitle("256cell, low_high")
ggsave(file = "/Volumes/USELESS/META/SHAPES/start/fold_change/low_high256.png")

#########
## 2-4cell, stable
meta = rep(0,63)
for (tr in stable) {
  if (txLengths[tr,]$utr5_len > 30 & sum(exons24[[tr]]$dtcr) > 0 & (txLengths[tr,]$utr5_len+33) < length(exons24[[tr]])) {
    meta_tr <- exons24[[tr]]$dtcr[(txLengths[tr,]$utr5_len-29):(txLengths[tr,]$utr5_len+33)] / sum(exons24[[tr]]$dtcr)
    if (length(meta_tr) == 63) {
      meta <- meta + meta_tr
    }
  }
}
# plot
df <- data.frame(scale,meta)
ggplot(df, aes(x=scale, y=meta)) + geom_bar(stat = "identity") + xlab("position from start codon") + ylab("sum(reactivities)") + ggtitle("2-4cell, stable")
ggsave(file = "/Volumes/USELESS/META/SHAPES/start/fold_change/stable24.png")

## 256cell, stable
meta = rep(0,63)
for (tr in stable) {
  if (txLengths[tr,]$utr5_len > 30 & sum(exons256[[tr]]$dtcr) > 0 & (txLengths[tr,]$utr5_len+33) < length(exons256[[tr]])) {
    meta_tr <- exons256[[tr]]$dtcr[(txLengths[tr,]$utr5_len-29):(txLengths[tr,]$utr5_len+33)] / sum(exons256[[tr]]$dtcr)
    if (length(meta_tr) == 63) {
      meta <- meta + meta_tr
    }
  }
}
# plot
df <- data.frame(scale,meta)
ggplot(df, aes(x=scale, y=meta)) + geom_bar(stat = "identity") + xlab("position from start codon") + ylab("sum(reactivities)") + ggtitle("256cell, stable")
ggsave(file = "/Volumes/USELESS/META/SHAPES/start/fold_change/stable256.png")

########
########

normalizeTx <- function(leader, cds, trailer, sam){
  leader <- leader[names(leader) %in% sam]
  cds <- cds[names(cds) %in% sam]
  trailer <- trailer[names(trailer) %in% sam]
  l <- leader[[1]]$dtcr
  c <- cds[[1]]$dtcr
  t <- trailer[[1]]$dtcr
  # leader
  sev <- data.table(
    start = head(seq(0,1,(1/length(l))), -1),
    end = tail(seq(0,1,(1/length(l))), -1),
    value = l / (sum(l) + sum(c) + sum(t))  # values, normalized
  )
  setkey(sev, start, end)
  for (name in tail(names(leader), -1)) {
    lnew <- leader[[name]]$dtcr
    cnew <- cds[[name]]$dtcr
    tnew <- trailer[[name]]$dtcr
    sevnew <- data.table(
      start = head(seq(0,1,(1/length(lnew))), -1),
      end = tail(seq(0,1,(1/length(lnew))), -1),
      value = lnew / (sum(lnew) + sum(cnew) + sum(tnew))  # values, normalized
    )
    setkey(sevnew, start, end)
    print(name)
    print(sev)
    print(sevnew)
    over <- foverlaps(sev, sevnew)
    over2 <- data.table(
      start = over[, ifelse(start > i.start, start, i.start)],
      end = over[, ifelse(end < i.end, end, i.end)],
      value = over[, value + i.value])
    # merge some intervals ? (smooth over, so to avoid running out of significant digits)
    sev <- over2
    sev <- sev[sev$start < sev$end]
    setkey(sev, start, end)
  }
  sevL <- sev
  # CDS
  sev <- data.table(
    start = head(seq(0,1,(1/length(c))), -1),
    end = tail(seq(0,1,(1/length(c))), -1),
    value = c / (sum(l) + sum(c) + sum(t))  # values, normalized
  )
  setkey(sev, start, end)
  for (name in tail(names(cds), -1)) {
    lnew <- leader[[name]]$dtcr
    cnew <- cds[[name]]$dtcr
    tnew <- trailer[[name]]$dtcr
    sevnew <- data.table(
      start = head(seq(0,1,(1/length(cnew))), -1),
      end = tail(seq(0,1,(1/length(cnew))), -1),
      value = cnew / (sum(lnew) + sum(cnew) + sum(tnew))  # values, normalized
    )
    setkey(sevnew, start, end)
    over <- foverlaps(sev, sevnew)
    over2 <- data.table(
      start = over[, ifelse(start > i.start, start, i.start)],
      end = over[, ifelse(end < i.end, end, i.end)],
      value = over[, value + i.value])
    # merge some intervals ? (smooth over, so to avoid running out of significant digits)
    sev <- over2
    sev <- sev[sev$start < sev$end]
    setkey(sev, start, end)
  }
  sevC <- sev
  sevC$start <- sevC$start + 1
  sevC$end <- sevC$end + 1
  # trailer
  sev <- data.table(
    start = head(seq(0,1,(1/length(t))), -1),
    end = tail(seq(0,1,(1/length(t))), -1),
    value = t / (sum(l) + sum(c) + sum(t))  # values, normalized
  )
  setkey(sev, start, end)
  for (name in tail(names(trailer), -1)) {
    lnew <- leader[[name]]$dtcr
    cnew <- cds[[name]]$dtcr
    tnew <- trailer[[name]]$dtcr
    sevnew <- data.table(
      start = head(seq(0,1,(1/length(tnew))), -1),
      end = tail(seq(0,1,(1/length(tnew))), -1),
      value = tnew / (sum(lnew) + sum(cnew) + sum(tnew))  # values, normalized
    )
    setkey(sevnew, start, end)
    over <- foverlaps(sev, sevnew)
    over2 <- data.table(
      start = over[, ifelse(start > i.start, start, i.start)],
      end = over[, ifelse(end < i.end, end, i.end)],
      value = over[, value + i.value])
    # merge some intervals ? (smooth over, so to avoid running out of significant digits)
    sev <- over2
    sev <- sev[sev$start < sev$end]
    setkey(sev, start, end)
  }
  sevT <- sev
  sevT$start <- sevT$start + 2
  sevT$end <- sevT$end + 2
  sev <- rbind(sevL, sevC, sevT)
  return(sev)
}

##
high_low <- high_low[txLengths[rownames(txLengths) %in% high_low,]$utr5_len > 30]
high_low <- high_low[sapply(leader24[names(leader24) %in% high_low], function(x){length(x) > 30})]
high_low <- high_low[sapply(trailer24[names(trailer24) %in% high_low], function(x){length(x) > 30})]
high_low <- high_low[sapply(cds24[names(cds24) %in% high_low], function(x){length(x) > 30})]

low_high <- low_high[txLengths[rownames(txLengths) %in% low_high,]$utr5_len > 30]
low_high <- low_high[sapply(leader24[names(leader24) %in% low_high], function(x){length(x) > 30})]
low_high <- low_high[sapply(trailer24[names(trailer24) %in% low_high], function(x){length(x) > 30})]
low_high <- low_high[sapply(cds24[names(cds24) %in% low_high], function(x){length(x) > 30})]

stable <- stable[txLengths[rownames(txLengths) %in% stable,]$utr5_len > 30]
stable <-stable[sapply(leader24[names(leader24) %in% stable], function(x){length(x) > 30})]
stable <-stable[sapply(trailer24[names(trailer24) %in% stable], function(x){length(x) > 30})]
stable <-stable[sapply(cds24[names(cds24) %in% stable], function(x){length(x) > 30})]



hl <- normalizeTx(leader24, cds24, trailer24, high_low)
  