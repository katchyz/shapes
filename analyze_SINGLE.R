### on SINGLE:
### start/stop on raw
### start/stop smoothed
library(GenomicFeatures)
library(ggplot2)
library(plyr)
library(data.table)
library(reshape2)

# libs:
# cell1, cell1K, oblong, oblong_CHX
# only TC: cell1_vitro

libs_path <- "/Volumes/USELESS/DATA/Shape-Seq"

libs <- c("shape_cell1.Rsave", "shape_cell1K.Rsave", "shape_oblong.Rsave", "shape_oblong_CHX.Rsave")
# "shapes_cell1_invitro.Rsave"


######### START
txdb_can <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/danio_rerio/Danio_rerio.GRCz10.84.gtf", format = "gtf")
txLengths <- transcriptLengths(txdb_can, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
rownames(txLengths) <- txLengths$tx_name
txLengths <- txLengths[order(rownames(txLengths)),]

# get longest transcript per gene
txlen <- arrange(txLengths, gene_id, desc(tx_len))
txlen <- txlen[!duplicated(txlen$gene_id),]

### tx with utr5 and cds longer than 100
txlen <- txlen[txlen$utr5_len > 100,]
txlen <- txlen[txlen$cds_len > 100,]
rownames(txlen) <- txlen$tx_name
txlen <- txlen[order(rownames(txlen)),]

start <- GRanges(Rle(rownames(txlen), rep(1, nrow(txlen))), IRanges(txlen$utr5_len-99, width=rep(200, nrow(txlen))))
start <- split(start, seqnames(start))

scale <- c(-100:100)[-101]


### 
save_path_start_log2ratio <- "/Volumes/USELESS/META/SHAPES_NEW/SINGLE/start/log2ratio"
save_path_start_TC <- "/Volumes/USELESS/META/SHAPES_NEW/SINGLE/start/TC"

for (i in 1:length(libs)) {
  u <- get(load(file = c(file.path(libs_path, libs[i]))))
  gr <- GRanges(seqnames=Rle(names(u)),
                IRanges(start(u), width=1),
                strand='*')
  gr$log2ratio <- u$log2ratio
  gr$TC.control <- u$TC.control
  gr$TC.treated <- u$TC.treated
  meta <- subsetByOverlaps(gr, start, ignore.strand=TRUE)
  meta <- split(meta, seqnames(meta))
  meta <- meta[sapply(meta, function(x){length(x) > 0})]
  metad <- meta[sapply(meta, function(x){sum(x$log2ratio) > 0})]
  metad <- Reduce("+", lapply(metad, function(x){x$log2ratio/sum(x$log2ratio)}))
  d <- data.frame(scale = scale, shape = metad)
  
  t <- paste("START, ", substr(libs[i], 7, (nchar(libs[i])-6)), ", log2ratio", sep="")
  p_start <- ggplot(d, aes(x = scale, y = shape)) + geom_bar(stat = "identity") + ggtitle(t)
  fstart = c(file.path(save_path_start_log2ratio, paste("start_", substr(libs[i], 7, (nchar(libs[i])-6)), ".png", sep="")))
  ggsave(file = fstart, p_start)
  
  metac <- meta[sapply(meta, function(x){sum(x$TC.control) > 0})]
  metac <- Reduce("+", lapply(metac, function(x){x$TC.control/sum(x$TC.control)}))
  metat <- meta[sapply(meta, function(x){sum(x$TC.treated) > 0})]
  metat <- Reduce("+", lapply(metat, function(x){x$TC.treated/sum(x$TC.treated)}))
  
  d <- data.frame(scale = scale, shape = metac)
  t <- paste("START, ", substr(libs[i], 7, (nchar(libs[i])-6)), ", TC.control", sep="")
  p_start <- ggplot(d, aes(x = scale, y = shape)) + geom_bar(stat = "identity") + ggtitle(t)
  fstart = c(file.path(save_path_start_TC, paste("start_", substr(libs[i], 7, (nchar(libs[i])-6)), "_CONTROL.png", sep="")))
  ggsave(file = fstart, p_start)
  
  d <- data.frame(scale = scale, shape = metat)
  t <- paste("START, ", substr(libs[i], 7, (nchar(libs[i])-6)), ", TC.treated", sep="")
  p_start <- ggplot(d, aes(x = scale, y = shape)) + geom_bar(stat = "identity") + ggtitle(t)
  fstart = c(file.path(save_path_start_TC, paste("start_", substr(libs[i], 7, (nchar(libs[i])-6)), "_TREATED.png", sep="")))
  ggsave(file = fstart, p_start)
  
}


####################################################################
####################################################################
######### STOP
# get longest transcript per gene
txlen <- arrange(txLengths, gene_id, desc(tx_len))
txlen <- txlen[!duplicated(txlen$gene_id),]

### tx with utr3 and cds longer than 100
txlen <- txlen[txlen$utr3_len > 100,]
txlen <- txlen[txlen$cds_len > 100,]
rownames(txlen) <- txlen$tx_name
txlen <- txlen[order(rownames(txlen)),]

stop <- GRanges(Rle(rownames(txlen), rep(1, nrow(txlen))), IRanges(txlen$utr5_len+txlen$cds_len-99, width=rep(200, nrow(txlen))))
stop <- split(stop, seqnames(stop))

scale <- c(-100:100)[-101]

### dCTR, TC, TCR
save_path_stop_log2ratio <- "/Volumes/USELESS/META/SHAPES_NEW/SINGLE/stop/log2ratio"
save_path_stop_TC <- "/Volumes/USELESS/META/SHAPES_NEW/SINGLE/stop/TC"

for (i in 1:length(libs)) {
  u <- get(load(file = c(file.path(libs_path, libs[i]))))
  gr <- GRanges(seqnames=Rle(names(u)),
                IRanges(start(u), width=1),
                strand='*')
  gr$log2ratio <- u$log2ratio
  gr$TC.control <- u$TC.control
  gr$TC.treated <- u$TC.treated
  
  meta <- subsetByOverlaps(gr, stop, ignore.strand=TRUE)
  meta <- split(meta, seqnames(meta))
  meta <- meta[sapply(meta, function(x){length(x) > 0})]
  
  # log2ratio
  metad <- meta[sapply(meta, function(x){sum(x$log2ratio) > 0})]
  metad <- Reduce("+", lapply(metad, function(x){x$log2ratio/sum(x$log2ratio)}))
  d <- data.frame(scale = scale, shape = metad)
  t <- paste("STOP, ", substr(libs[i], 7, (nchar(libs[i])-6)), ", log2ratio", sep="")
  p_stop <- ggplot(d, aes(x = scale, y = shape)) + geom_bar(stat = "identity") + ggtitle(t)
  fstop = c(file.path(save_path_stop_log2ratio, paste("stop_", substr(libs[i], 7, (nchar(libs[i])-6)), ".png", sep="")))
  ggsave(file = fstop, p_stop)
  
  # TC
  metac <- meta[sapply(meta, function(x){sum(x$TC.control) > 0})]
  metac <- Reduce("+", lapply(metac, function(x){x$TC.control/sum(x$TC.control)}))
  metat <- meta[sapply(meta, function(x){sum(x$TC.treated) > 0})]
  metat <- Reduce("+", lapply(metat, function(x){x$TC.treated/sum(x$TC.treated)}))
  d <- data.frame(scale = scale, shape = metac)
  t <- paste("STOP, ", substr(libs[i], 7, (nchar(libs[i])-6)), ", TC.control", sep="")
  p_stop <- ggplot(d, aes(x = scale, y = shape)) + geom_bar(stat = "identity") + ggtitle(t)
  fstop = c(file.path(save_path_stop_TC, paste("stop_", substr(libs[i], 7, (nchar(libs[i])-6)), "_CONTROL.png", sep="")))
  ggsave(file = fstop, p_stop)
  d <- data.frame(scale = scale, shape = metat)
  t <- paste("STOP, ", substr(libs[i], 7, (nchar(libs[i])-6)), ", TC.treated", sep="")
  p_stop <- ggplot(d, aes(x = scale, y = shape)) + geom_bar(stat = "identity") + ggtitle(t)
  fstop = c(file.path(save_path_stop_TC, paste("stop_", substr(libs[i], 7, (nchar(libs[i])-6)), "_TREATED.png", sep="")))
  ggsave(file = fstop, p_stop)
  
}


####################################################################
####################################################################
######### tx_profiles

load(file="/Volumes/USELESS/META/SHAPES_NEW/general/clouds/cloud_df.Rsave")

lct <- txlen[(txlen$utr5_len > 100) & (txlen$cds_len > 100) & (txlen$utr3_len > 100), ]

lct <- lct[lct$tx_name %in% df[df$rna_cell1K > 1,]$n,]
lct <- lct[lct$tx_name %in% df[df$shape_cell1K_RZ_dTCR_internal > 0.01,]$n,]

sam <- sample(lct$tx_name, 500)

#### normalizeTx (from tx_profiles.R)
normalizeTx <- function(shape, txlen, sam, norm="dTCR"){
  # subset shape GRanges
  shape <- shape[names(shape) %in% sam]
  shape <- shape[order(names(shape))]
  len <- txlen[txlen$tx_name %in% names(shape),]
  len <- len[order(len$tx_name),]
  
  leader <- mapply(function(x,y){mcols(x)[[norm]][1:y]}, shape, len$utr5_len)
  shape <- mapply(function(x,y){mcols(x)[[norm]][(y+1):length(x)]}, shape, len$utr5_len)
  cds <- mapply(function(x,y){x[1:y]}, shape, len$cds_len)
  trailer <- mapply(function(x,y){x[(y+1):length(x)]}, shape, len$cds_len)
  
  l <- leader[[1]]
  c <- cds[[1]]
  t <- trailer[[1]]
  
  # leader
  sev <- data.table(
    start = head(seq(0,1,(1/length(l))), -1),
    end = tail(seq(0,1,(1/length(l))), -1),
    value = l / (sum(l) + sum(c) + sum(t))  # values, normalized
  )
  setkey(sev, start, end)
  for (name in tail(names(leader), -1)) {
    lnew <- leader[[name]]
    cnew <- cds[[name]]
    tnew <- trailer[[name]]
    sevnew <- data.table(
      start = head(seq(0,1,(1/length(lnew))), -1),
      end = tail(seq(0,1,(1/length(lnew))), -1),
      value = lnew / (sum(lnew) + sum(cnew) + sum(tnew))  # values, normalized
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
  sevL <- sev
  # CDS
  sev <- data.table(
    start = head(seq(0,1,(1/length(c))), -1),
    end = tail(seq(0,1,(1/length(c))), -1),
    value = c / (sum(l) + sum(c) + sum(t))  # values, normalized
  )
  setkey(sev, start, end)
  for (name in tail(names(cds), -1)) {
    lnew <- leader[[name]]
    cnew <- cds[[name]]
    tnew <- trailer[[name]]
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
    lnew <- leader[[name]]
    cnew <- cds[[name]]
    tnew <- trailer[[name]]
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



save_path_tx_log2ratio <- "/Volumes/USELESS/META/SHAPES_NEW/SINGLE/tx_profiles/utr5_cds_utr3/log2ratio"
save_path_tx_TC <- "/Volumes/USELESS/META/SHAPES_NEW/SINGLE/tx_profiles/utr5_cds_utr3/TC"

for (i in 1:length(libs)) {
  u <- get(load(file = c(file.path(libs_path, libs[i]))))
  gr <- GRanges(seqnames=Rle(names(u)),
                IRanges(start(u), width=1),
                strand='*')
  gr$log2ratio <- u$log2ratio
  gr$TC.control <- u$TC.control
  gr$TC.treated <- u$TC.treated
  gr <- split(gr, seqnames(gr))
  
  sev <- normalizeTx(gr, txlen, sam, norm="log2ratio")
  t <- paste("tx_profile, ", substr(libs[i], 7, (nchar(libs[i])-6)), ", log2ratio", sep="")
  p_tx <- ggplot() + geom_step(data = sev, aes(start, value), size=0.1) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) +
    xlab("leader, CDS, trailer") + ylab("sum(SHAPE)") + ggtitle(t)
  fp = c(file.path(save_path_tx_log2ratio, paste("tx_", substr(libs[i], 7, (nchar(libs[i])-6)), ".png", sep="")))
  ggsave(file = fp, p_tx)
  
  sev <- normalizeTx(gr, txlen, sam, norm="TC.control")
  t <- paste("tx_profile, ", substr(libs[i], 7, (nchar(libs[i])-6)), ", TC.control", sep="")
  p_tx <- ggplot() + geom_step(data = sev, aes(start, value), size=0.1) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) +
    xlab("leader, CDS, trailer") + ylab("sum(SHAPE)") + ggtitle(t)
  fp = c(file.path(save_path_tx_TC, paste("tx_", substr(libs[i], 7, (nchar(libs[i])-6)), "_CONTROL.png", sep="")))
  ggsave(file = fp, p_tx)
  
  sev <- normalizeTx(gr, txlen, sam, norm="TC.treated")
  t <- paste("tx_profile, ", substr(libs[i], 7, (nchar(libs[i])-6)), ", TC.treated", sep="")
  p_tx <- ggplot() + geom_step(data = sev, aes(start, value), size=0.1) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) +
    xlab("leader, CDS, trailer") + ylab("sum(SHAPE)") + ggtitle(t)
  fp = c(file.path(save_path_tx_TC, paste("tx_", substr(libs[i], 7, (nchar(libs[i])-6)), "_TREATED.png", sep="")))
  ggsave(file = fp, p_tx)
  
}


###################################################
###################################################
######### tx_profiles: whole tx to unit length

######### tx_profiles

load(file="/Volumes/USELESS/META/SHAPES_NEW/general/clouds/cloud_df.Rsave")

lct <- txlen[(txlen$utr5_len > 50) & (txlen$cds_len > 100) & (txlen$utr3_len > 50), ]

lct <- lct[lct$tx_name %in% df[df$rna_cell1K > 1,]$n,]
lct <- lct[lct$tx_name %in% df[df$shape_cell1K_RZ_dTCR_internal > 0.01,]$n,]

sam <- sample(lct$tx_name, 1000)

#### normalizeTx (from tx_profiles.R)
normalizeTxUL <- function(shape, sam, norm="dTCR"){
  # subset shape GRanges
  shape <- shape[names(shape) %in% sam]
  shape <- shape[order(names(shape))]
  
  tx <- sapply(shape, function(x){mcols(x)[[norm]]})
  t <- tx[[1]]
  # tx
  sev <- data.table(
    start = head(seq(0,1,(1/length(t))), -1),
    end = tail(seq(0,1,(1/length(t))), -1),
    value = t / sum(t)  # values, normalized
  )
  setkey(sev, start, end)
  for (name in tail(names(tx), -1)) {
    tnew <- tx[[name]]
    sevnew <- data.table(
      start = head(seq(0,1,(1/length(tnew))), -1),
      end = tail(seq(0,1,(1/length(tnew))), -1),
      value = tnew / sum(tnew)  # values, normalized
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
  return(sev)
}


save_path_tx_log2ratio <- "/Volumes/USELESS/META/SHAPES_NEW/SINGLE/tx_profiles/tx_unit_length/log2ratio"
save_path_tx_TC <- "/Volumes/USELESS/META/SHAPES_NEW/SINGLE/tx_profiles/tx_unit_length/TC"

for (i in 1:length(libs)) {
  u <- get(load(file = c(file.path(libs_path, libs[i]))))
  gr <- GRanges(seqnames=Rle(names(u)),
                IRanges(start(u), width=1),
                strand='*')
  gr$log2ratio <- u$log2ratio
  gr$TC.control <- u$TC.control
  gr$TC.treated <- u$TC.treated
  gr <- split(gr, seqnames(gr))
  
  sev <- normalizeTxUL(gr, sam, norm="log2ratio")
  t <- paste("tx_profile, ", substr(libs[i], 7, (nchar(libs[i])-6)), ", log2ratio", sep="")
  p_tx <- ggplot() + geom_step(data = sev, aes(start, value), size=0.1) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) +
    xlab("leader, CDS, trailer") + ylab("sum(SHAPE)") + ggtitle(t)
  fp = c(file.path(save_path_tx_log2ratio, paste("tx_", substr(libs[i], 7, (nchar(libs[i])-6)), ".png", sep="")))
  ggsave(file = fp, p_tx)
  
  sev <- normalizeTxUL(gr, sam, norm="TC.control")
  t <- paste("tx_profile, ", substr(libs[i], 7, (nchar(libs[i])-6)), ", TC.control", sep="")
  p_tx <- ggplot() + geom_step(data = sev, aes(start, value), size=0.1) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) +
    xlab("leader, CDS, trailer") + ylab("sum(SHAPE)") + ggtitle(t)
  fp = c(file.path(save_path_tx_TC, paste("tx_", substr(libs[i], 7, (nchar(libs[i])-6)), "_CONTROL.png", sep="")))
  ggsave(file = fp, p_tx)
  
  sev <- normalizeTxUL(gr, sam, norm="TC.treated")
  t <- paste("tx_profile, ", substr(libs[i], 7, (nchar(libs[i])-6)), ", TC.treated", sep="")
  p_tx <- ggplot() + geom_step(data = sev, aes(start, value), size=0.1) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) +
    xlab("leader, CDS, trailer") + ylab("sum(SHAPE)") + ggtitle(t)
  fp = c(file.path(save_path_tx_TC, paste("tx_", substr(libs[i], 7, (nchar(libs[i])-6)), "_TREATED.png", sep="")))
  ggsave(file = fp, p_tx)
  
}

#####
#### normalizeTx (from tx_profiles.R)
normalizeTxUL <- function(shape, sam, norm="dTCR"){
  # subset shape GRanges
  shape <- shape[names(shape) %in% sam]
  shape <- shape[order(names(shape))]
  
  tx <- sapply(shape, function(x){mcols(x)[[norm]]})
  t <- tx[[1]]
  if (sum(t) > 0) {
    val = t / sum(t)
  } else {
    val = t
  }
  # tx
  sev <- data.table(
    start = head(seq(0,1,(1/length(t))), -1),
    end = tail(seq(0,1,(1/length(t))), -1),
    value = val  # values, normalized
  )
  setkey(sev, start, end)
  for (name in tail(names(tx), -1)) {
    tnew <- tx[[name]]
    if (sum(tnew) > 0) {
      val = tnew / sum(tnew)
    } else {
      val = tnew
    }
    sevnew <- data.table(
      start = head(seq(0,1,(1/length(tnew))), -1),
      end = tail(seq(0,1,(1/length(tnew))), -1),
      value = val  # values, normalized
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
  return(sev)
}


save_path_tx_TC <- "/Volumes/USELESS/META/SHAPES_NEW/SINGLE/tx_profiles/tx_unit_length/TC"

for (i in 1:length(libs)) {
  u <- get(load(file = c(file.path(libs_path, libs[i]))))
  gr <- GRanges(seqnames=Rle(names(u)),
                IRanges(start(u), width=1),
                strand='*')
  gr$log2ratio <- u$log2ratio
  gr$TC.control <- u$TC.control
  gr$TC.treated <- u$TC.treated
  gr <- split(gr, seqnames(gr))
  
  sev <- normalizeTxUL(gr, sam, norm="TC.control")
  t <- paste("tx_profile, ", substr(libs[i], 7, (nchar(libs[i])-6)), ", TC.control", sep="")
  p_tx <- ggplot() + geom_step(data = sev, aes(start, value), size=0.1) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) +
    xlab("leader, CDS, trailer") + ylab("sum(SHAPE)") + ggtitle(t)
  fp = c(file.path(save_path_tx_TC, paste("tx_", substr(libs[i], 7, (nchar(libs[i])-6)), "_CONTROL.png", sep="")))
  ggsave(file = fp, p_tx)
  
}


# shapes nonsel

for (i in 1:length(libs)) {
  u <- get(load(file = c(file.path(libs_path, libs[i]))))
  gr <- GRanges(seqnames=Rle(names(u)),
                IRanges(start(u), width=1),
                strand='*')
  gr$TC <- u$TC
  
  meta <- subsetByOverlaps(gr, start, ignore.strand=TRUE)
  meta <- split(meta, seqnames(meta))
  meta <- meta[sapply(meta, function(x){length(x) > 0})]
  
  # TC
  metac <- meta[sapply(meta, function(x){sum(x$TC) > 0})]
  metac <- Reduce("+", lapply(metac, function(x){x$TC/sum(x$TC)}))
  d <- data.frame(scale = scale, shape = metac)
  p_start <- ggplot(d, aes(x = scale, y = shape)) + geom_bar(stat = "identity") + ggtitle("shapes_cell1_invitro_nonsel")
  fstart = c(file.path(save_path_start_TC, "shapes_cell1_invitro_nonsel.png"))
  ggsave(file = fstart, p_start)
  
}



save_path_tx_TC <- "/Volumes/USELESS/META/SHAPES_NEW/SINGLE/tx_profiles/tx_unit_length/TC"


for (i in 1:length(libs)) {
  u <- get(load(file = c(file.path(libs_path, libs[i]))))
  gr <- GRanges(seqnames=Rle(names(u)),
                IRanges(start(u), width=1),
                strand='*')
  gr$TC <- u$TC
  gr <- split(gr, seqnames(gr))
  
  sev <- normalizeTxUL(gr, sam, norm="TC")
  p_tx <- ggplot() + geom_step(data = sev, aes(start, value), size=0.1) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) +
    xlab("leader, CDS, trailer") + ylab("sum(SHAPE)") + ggtitle("shapes_cell1_invitro_nonsel")
  fp = c(file.path(save_path_tx_TC, "shapes_cell1_invitro_nonsel.png"))
  ggsave(file = fp, p_tx)
  
}


save_path_tx_TC <- "/Volumes/USELESS/META/SHAPES_NEW/SINGLE/tx_profiles/utr5_cds_utr3/TC"
libs <- c("shapes_cell1_invitro.Rsave", "shapes_cell1_invitro_nonsel.Rsave")

for (i in 1:length(libs)) {
  u <- get(load(file = c(file.path(libs_path, libs[i]))))
  gr <- GRanges(seqnames=Rle(names(u)),
                IRanges(start(u), width=1),
                strand='*')
  gr$TC <- u$TC
  gr <- split(gr, seqnames(gr))
  gr <- gr[sapply(gr, function(x){sum(x$TC) > 0})]
  
  sev <- normalizeTx(gr, txlen, sam, norm="TC")
  t <- substr(libs[i], 1, (nchar(libs[i])-6))
  p_tx <- ggplot() + geom_step(data = sev, aes(start, value), size=0.1) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) +
    xlab("leader, CDS, trailer") + ylab("sum(SHAPE)") + ggtitle(t)
  fp = c(file.path(save_path_tx_TC, paste(substr(libs[i], 1, (nchar(libs[i])-6)), ".png", sep="")))
  ggsave(file = fp, p_tx)
  
}



libs <- c("shape_cell1K.Rsave")

for (i in 1:length(libs)) {
  u <- get(load(file = c(file.path(libs_path, libs[i]))))
  gr <- GRanges(seqnames=Rle(names(u)),
                IRanges(start(u), width=1),
                strand='*')
  gr$log2ratio <- u$log2ratio
  gr <- split(gr, seqnames(gr))
  gr <- gr[sapply(gr, function(x){sum(x$log2ratio) > 0})]
  
  sev <- normalizeTx(gr, txlen, sam, norm="log2ratio")
  t <- paste("tx_profile, ", substr(libs[i], 7, (nchar(libs[i])-6)), ", log2ratio", sep="")
  p_tx <- ggplot() + geom_step(data = sev, aes(start, value), size=0.1) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) +
    xlab("leader, CDS, trailer") + ylab("sum(SHAPE)") + ggtitle(t)
  fp = c(file.path(save_path_tx_log2ratio, paste("tx_", substr(libs[i], 7, (nchar(libs[i])-6)), ".png", sep="")))
  ggsave(file = fp, p_tx)
  
}




libs <- c("shape_cell1.Rsave", "shape_cell1K.Rsave", "shape_oblong.Rsave")

for (i in 1:length(libs)) {
  u <- get(load(file = c(file.path(libs_path, libs[i]))))
  gr <- GRanges(seqnames=Rle(names(u)),
                IRanges(start(u), width=1),
                strand='*')
  gr$log2ratio <- u$log2ratio
  gr <- split(gr, seqnames(gr))
  
  sev <- normalizeTxUL(gr, sam, norm="log2ratio")
  t <- paste("tx_profile, ", substr(libs[i], 7, (nchar(libs[i])-6)), ", log2ratio", sep="")
  p_tx <- ggplot() + geom_step(data = sev, aes(start, value), size=0.1) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) +
    xlab("leader, CDS, trailer") + ylab("sum(SHAPE)") + ggtitle(t)
  fp = c(file.path(save_path_tx_log2ratio, paste("tx_", substr(libs[i], 7, (nchar(libs[i])-6)), ".png", sep="")))
  ggsave(file = fp, p_tx)
  
}
