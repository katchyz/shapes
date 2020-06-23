### on dTCR:
### start/stop on raw
### start/stop smoothed
library(GenomicFeatures)
library(ggplot2)
library(plyr)
library(data.table)
library(reshape2)

# libs:
# 1Kcell_polyA, 1Kcell_RZ, vivo
# 1Kcell_polyA_vitro, 1Kcell_RZ_vitro, vitro
# 2-4cell
# 256cell

libs_path <- "/Volumes/USELESS/DATA/Shape-Seq/dTCR"

cell24 <- get(load(file = c(file.path(libs_path, "PAIR_GC_cell24.Rsave"))))
cell256 <- get(load(file = c(file.path(libs_path, "PAIR_GC_cell256.Rsave"))))
cell1K <- get(load(file = c(file.path(libs_path, "1Kcell_vivo_dTCR.Rsave"))))
cell1K_vitro <- get(load(file = c(file.path(libs_path, "1Kcell_vitro_dTCR.Rsave"))))

rm(GR)

cell1K_polyA <- load(file = c(file.path(libs_path, "PAIR_GC_1Kcell_polyA.Rsave")))
cell1K_RZ <- load(file = c(file.path(libs_path, "PAIR_GC_1Kcell_RZ.Rsave")))
cell1K_polyA_vitro <- load(file = c(file.path(libs_path, "PAIR_GC_1Kcell_polyA_vitro.Rsave")))
cell1K_RZ_vitro <- load(file = c(file.path(libs_path, "PAIR_GC_1Kcell_RZ_vitro.Rsave")))


# 1Kcell_vitro_dTCR.Rsave             vitro
# 1Kcell_vivo_dTCR.Rsave              vivo
# PAIR_GC_1Kcell_RZ.Rsave             GR
# PAIR_GC_1Kcell_RZ_vitro.Rsave       GR
# PAIR_GC_1Kcell_polyA.Rsave          GR
# PAIR_GC_1Kcell_polyA_vitro.Rsave    GR
# PAIR_GC_cell24.Rsave                GR
# PAIR_GC_cell256.Rsave               GR

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

############ subset by overlaps

#libs <- c("PAIR_GC_cell24.Rsave", "PAIR_GC_cell256.Rsave", "1Kcell_vivo_dTCR.Rsave",
          "1Kcell_vitro_dTCR.Rsave")

libs <- c("PAIR_GC_cell24.Rsave", "PAIR_GC_cell256.Rsave", "1Kcell_vivo_dTCR.Rsave",
          "1Kcell_vitro_dTCR.Rsave", "PAIR_GC_1Kcell_RZ.Rsave",
          "PAIR_GC_1Kcell_RZ_vitro.Rsave", "PAIR_GC_1Kcell_polyA.Rsave",
          "PAIR_GC_1Kcell_polyA_vitro.Rsave")


### dCTR
save_path_start <- "/Volumes/USELESS/META/SHAPES_NEW/PAIRED/start/dTCR"

for (i in 1:length(libs)) {
  u <- get(load(file = c(file.path(libs_path, libs[i]))))
  gr <- GRanges(seqnames=Rle(names(u)),
                IRanges(start(u), width=1),
                strand='*')
  gr$dTCR <- u$dTCR
  meta <- subsetByOverlaps(gr, start, ignore.strand=TRUE)
  meta <- split(meta, seqnames(meta))
  meta <- meta[sapply(meta, function(x){length(x) > 0})]
  metad <- meta[sapply(meta, function(x){sum(x$dTCR) > 0})]
  metad <- Reduce("+", lapply(metad, function(x){x$dTCR/sum(x$dTCR)}))
  d <- data.frame(scale = scale, shape = metad)
  
  t <- paste("START, ", substr(libs[i], 9, (nchar(libs[i])-6)), ", dTCR", sep="")
  p_start <- ggplot(d, aes(x = scale, y = shape)) + geom_bar(stat = "identity") + ggtitle(t)
  fstart = c(file.path(save_path_start, paste("start_", substr(libs[i], 9, (nchar(libs[i])-6)), ".png", sep="")))
  ggsave(file = fstart, p_start)
}

### TC
save_path_start <- "/Volumes/USELESS/META/SHAPES_NEW/PAIRED/start/TC"

for (i in 1:length(libs)) {
  u <- get(load(file = c(file.path(libs_path, libs[i]))))
  gr <- GRanges(seqnames=Rle(names(u)),
                IRanges(start(u), width=1),
                strand='*')
  gr$TC.control <- u$TC.control
  gr$TC.treated <- u$TC.treated
  meta <- subsetByOverlaps(gr, start, ignore.strand=TRUE)
  meta <- split(meta, seqnames(meta))
  meta <- meta[sapply(meta, function(x){length(x) > 0})]
  
  metac <- meta[sapply(meta, function(x){sum(x$TC.control) > 0})]
  metac <- Reduce("+", lapply(metac, function(x){x$TC.control/sum(x$TC.control)}))
  metat <- meta[sapply(meta, function(x){sum(x$TC.treated) > 0})]
  metat <- Reduce("+", lapply(metat, function(x){x$TC.treated/sum(x$TC.treated)}))
  
  d <- data.frame(scale = scale, shape = metac)
  t <- paste("START, ", substr(libs[i], 9, (nchar(libs[i])-6)), ", TC.control", sep="")
  p_start <- ggplot(d, aes(x = scale, y = shape)) + geom_bar(stat = "identity") + ggtitle(t)
  fstart = c(file.path(save_path_start, paste("start_", substr(libs[i], 9, (nchar(libs[i])-6)), "_CONTROL.png", sep="")))
  ggsave(file = fstart, p_start)
  
  d <- data.frame(scale = scale, shape = metat)
  t <- paste("START, ", substr(libs[i], 9, (nchar(libs[i])-6)), ", TC.treated", sep="")
  p_start <- ggplot(d, aes(x = scale, y = shape)) + geom_bar(stat = "identity") + ggtitle(t)
  fstart = c(file.path(save_path_start, paste("start_", substr(libs[i], 9, (nchar(libs[i])-6)), "_TREATED.png", sep="")))
  ggsave(file = fstart, p_start)
}

### TCR
save_path_start <- "/Volumes/USELESS/META/SHAPES_NEW/PAIRED/start/TCR"

for (i in 1:length(libs)) {
  u <- get(load(file = c(file.path(libs_path, libs[i]))))
  gr <- GRanges(seqnames=Rle(names(u)),
                IRanges(start(u), width=1),
                strand='*')
  gr$TCR.control <- u$TCR.control
  gr$TCR.treated <- u$TCR.treated
  meta <- subsetByOverlaps(gr, start, ignore.strand=TRUE)
  meta <- split(meta, seqnames(meta))
  meta <- meta[sapply(meta, function(x){length(x) > 0})]
  
  metac <- meta[sapply(meta, function(x){sum(x$TCR.control) > 0})]
  metac <- Reduce("+", lapply(metac, function(x){x$TCR.control/sum(x$TCR.control)}))
  metat <- meta[sapply(meta, function(x){sum(x$TCR.treated) > 0})]
  metat <- Reduce("+", lapply(metat, function(x){x$TCR.treated/sum(x$TCR.treated)}))
  
  d <- data.frame(scale = scale, shape = metac)
  t <- paste("START, ", substr(libs[i], 9, (nchar(libs[i])-6)), ", TCR.control", sep="")
  p_start <- ggplot(d, aes(x = scale, y = shape)) + geom_bar(stat = "identity") + ggtitle(t)
  fstart = c(file.path(save_path_start, paste("start_", substr(libs[i], 9, (nchar(libs[i])-6)), "_CONTROL.png", sep="")))
  ggsave(file = fstart, p_start)
  
  d <- data.frame(scale = scale, shape = metat)
  t <- paste("START, ", substr(libs[i], 9, (nchar(libs[i])-6)), ", TCR.treated", sep="")
  p_start <- ggplot(d, aes(x = scale, y = shape)) + geom_bar(stat = "identity") + ggtitle(t)
  fstart = c(file.path(save_path_start, paste("start_", substr(libs[i], 9, (nchar(libs[i])-6)), "_TREATED.png", sep="")))
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
save_path_stop_dTCR <- "/Volumes/USELESS/META/SHAPES_NEW/PAIRED/stop/dTCR"
save_path_stop_TC <- "/Volumes/USELESS/META/SHAPES_NEW/PAIRED/stop/TC"
save_path_stop_TCR <- "/Volumes/USELESS/META/SHAPES_NEW/PAIRED/stop/TCR"

for (i in 1:length(libs)) {
  u <- get(load(file = c(file.path(libs_path, libs[i]))))
  gr <- GRanges(seqnames=Rle(names(u)),
                IRanges(start(u), width=1),
                strand='*')
  gr$dTCR <- u$dTCR
  gr$TC.control <- u$TC.control
  gr$TC.treated <- u$TC.treated
  gr$TCR.control <- u$TCR.control
  gr$TCR.treated <- u$TCR.treated
  
  meta <- subsetByOverlaps(gr, stop, ignore.strand=TRUE)
  meta <- split(meta, seqnames(meta))
  meta <- meta[sapply(meta, function(x){length(x) > 0})]
  
  # dTCR
  metad <- meta[sapply(meta, function(x){sum(x$dTCR) > 0})]
  metad <- Reduce("+", lapply(metad, function(x){x$dTCR/sum(x$dTCR)}))
  d <- data.frame(scale = scale, shape = metad)
  t <- paste("STOP, ", substr(libs[i], 9, (nchar(libs[i])-6)), ", dTCR", sep="")
  p_stop <- ggplot(d, aes(x = scale, y = shape)) + geom_bar(stat = "identity") + ggtitle(t)
  fstop = c(file.path(save_path_stop_dTCR, paste("stop_", substr(libs[i], 9, (nchar(libs[i])-6)), ".png", sep="")))
  ggsave(file = fstop, p_stop)
  
  # TC
  metac <- meta[sapply(meta, function(x){sum(x$TC.control) > 0})]
  metac <- Reduce("+", lapply(metac, function(x){x$TC.control/sum(x$TC.control)}))
  metat <- meta[sapply(meta, function(x){sum(x$TC.treated) > 0})]
  metat <- Reduce("+", lapply(metat, function(x){x$TC.treated/sum(x$TC.treated)}))
  d <- data.frame(scale = scale, shape = metac)
  t <- paste("STOP, ", substr(libs[i], 9, (nchar(libs[i])-6)), ", TC.control", sep="")
  p_stop <- ggplot(d, aes(x = scale, y = shape)) + geom_bar(stat = "identity") + ggtitle(t)
  fstop = c(file.path(save_path_stop_TC, paste("stop_", substr(libs[i], 9, (nchar(libs[i])-6)), "_CONTROL.png", sep="")))
  ggsave(file = fstop, p_stop)
  d <- data.frame(scale = scale, shape = metat)
  t <- paste("STOP, ", substr(libs[i], 9, (nchar(libs[i])-6)), ", TC.treated", sep="")
  p_stop <- ggplot(d, aes(x = scale, y = shape)) + geom_bar(stat = "identity") + ggtitle(t)
  fstop = c(file.path(save_path_stop_TC, paste("stop_", substr(libs[i], 9, (nchar(libs[i])-6)), "_TREATED.png", sep="")))
  ggsave(file = fstop, p_stop)
  
  # TCR
  metac <- meta[sapply(meta, function(x){sum(x$TC.control) > 0})]
  metac <- Reduce("+", lapply(metac, function(x){x$TC.control/sum(x$TC.control)}))
  metat <- meta[sapply(meta, function(x){sum(x$TC.treated) > 0})]
  metat <- Reduce("+", lapply(metat, function(x){x$TC.treated/sum(x$TC.treated)}))
  d <- data.frame(scale = scale, shape = metac)
  t <- paste("STOP, ", substr(libs[i], 9, (nchar(libs[i])-6)), ", TCR.control", sep="")
  p_stop <- ggplot(d, aes(x = scale, y = shape)) + geom_bar(stat = "identity") + ggtitle(t)
  fstop = c(file.path(save_path_stop_TCR, paste("stop_", substr(libs[i], 9, (nchar(libs[i])-6)), "_CONTROL.png", sep="")))
  ggsave(file = fstop, p_stop)
  d <- data.frame(scale = scale, shape = metat)
  t <- paste("STOP, ", substr(libs[i], 9, (nchar(libs[i])-6)), ", TCR.treated", sep="")
  p_stop <- ggplot(d, aes(x = scale, y = shape)) + geom_bar(stat = "identity") + ggtitle(t)
  fstop = c(file.path(save_path_stop_TCR, paste("stop_", substr(libs[i], 9, (nchar(libs[i])-6)), "_TREATED.png", sep="")))
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


save_path_tx_dTCR <- "/Volumes/USELESS/META/SHAPES_NEW/PAIRED/tx_profiles/utr5_cds_utr3/dTCR"
save_path_tx_TC <- "/Volumes/USELESS/META/SHAPES_NEW/PAIRED/tx_profiles/utr5_cds_utr3/TC"
save_path_tx_TCR <- "/Volumes/USELESS/META/SHAPES_NEW/PAIRED/tx_profiles/utr5_cds_utr3/TCR"

for (i in 1:length(libs)) {
  u <- get(load(file = c(file.path(libs_path, libs[i]))))
  gr <- GRanges(seqnames=Rle(names(u)),
                IRanges(start(u), width=1),
                strand='*')
  gr$dTCR <- u$dTCR
  gr$TC.control <- u$TC.control
  gr$TC.treated <- u$TC.treated
  gr$TCR.control <- u$TCR.control
  gr$TCR.treated <- u$TCR.treated
  gr <- split(gr, seqnames(gr))
  
  sev <- normalizeTx(gr, txlen, sam, norm="dTCR")
  t <- paste("tx_profile, ", substr(libs[i], 9, (nchar(libs[i])-6)), ", dTCR", sep="")
  p_tx <- ggplot() + geom_step(data = sev, aes(start, value), size=0.1) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) +
    xlab("leader, CDS, trailer") + ylab("sum(SHAPE)") + ggtitle(t)
  fp = c(file.path(save_path_tx_dTCR, paste("tx_", substr(libs[i], 9, (nchar(libs[i])-6)), ".png", sep="")))
  ggsave(file = fp, p_tx)
  
  sev <- normalizeTx(gr, txlen, sam, norm="TC.control")
  t <- paste("tx_profile, ", substr(libs[i], 9, (nchar(libs[i])-6)), ", TC.control", sep="")
  p_tx <- ggplot() + geom_step(data = sev, aes(start, value), size=0.1) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) +
    xlab("leader, CDS, trailer") + ylab("sum(SHAPE)") + ggtitle(t)
  fp = c(file.path(save_path_tx_TC, paste("tx_", substr(libs[i], 9, (nchar(libs[i])-6)), "_CONTROL.png", sep="")))
  ggsave(file = fp, p_tx)
  
  sev <- normalizeTx(gr, txlen, sam, norm="TC.treated")
  t <- paste("tx_profile, ", substr(libs[i], 9, (nchar(libs[i])-6)), ", TC.treated", sep="")
  p_tx <- ggplot() + geom_step(data = sev, aes(start, value), size=0.1) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) +
    xlab("leader, CDS, trailer") + ylab("sum(SHAPE)") + ggtitle(t)
  fp = c(file.path(save_path_tx_TC, paste("tx_", substr(libs[i], 9, (nchar(libs[i])-6)), "_TREATED.png", sep="")))
  ggsave(file = fp, p_tx)
  
  sev <- normalizeTx(gr, txlen, sam, norm="TCR.control")
  t <- paste("tx_profile, ", substr(libs[i], 9, (nchar(libs[i])-6)), ", TCR.control", sep="")
  p_tx <- ggplot() + geom_step(data = sev, aes(start, value), size=0.1) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) +
    xlab("leader, CDS, trailer") + ylab("sum(SHAPE)") + ggtitle(t)
  fp = c(file.path(save_path_tx_TCR, paste("tx_", substr(libs[i], 9, (nchar(libs[i])-6)), "_CONTROL.png", sep="")))
  ggsave(file = fp, p_tx)
  
  sev <- normalizeTx(gr, txlen, sam, norm="TCR.treated")
  t <- paste("tx_profile, ", substr(libs[i], 9, (nchar(libs[i])-6)), ", TCR.treated", sep="")
  p_tx <- ggplot() + geom_step(data = sev, aes(start, value), size=0.1) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) +
    xlab("leader, CDS, trailer") + ylab("sum(SHAPE)") + ggtitle(t)
  fp = c(file.path(save_path_tx_TCR, paste("tx_", substr(libs[i], 9, (nchar(libs[i])-6)), "_TREATED.png", sep="")))
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


save_path_tx_dTCR <- "/Volumes/USELESS/META/SHAPES_NEW/PAIRED/tx_profiles/tx_unit_length/dTCR"
save_path_tx_TC <- "/Volumes/USELESS/META/SHAPES_NEW/PAIRED/tx_profiles/tx_unit_length/TC"
save_path_tx_TCR <- "/Volumes/USELESS/META/SHAPES_NEW/PAIRED/tx_profiles/tx_unit_length/TCR"

for (i in 1:length(libs)) {
  u <- get(load(file = c(file.path(libs_path, libs[i]))))
  gr <- GRanges(seqnames=Rle(names(u)),
                IRanges(start(u), width=1),
                strand='*')
  gr$dTCR <- u$dTCR
  gr$TC.control <- u$TC.control
  gr$TC.treated <- u$TC.treated
  gr$TCR.control <- u$TCR.control
  gr$TCR.treated <- u$TCR.treated
  gr <- split(gr, seqnames(gr))
  
  sev <- normalizeTxUL(gr, sam, norm="dTCR")
  t <- paste("tx_profile, ", substr(libs[i], 9, (nchar(libs[i])-6)), ", dTCR", sep="")
  p_tx <- ggplot() + geom_step(data = sev, aes(start, value), size=0.1) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) +
    xlab("leader, CDS, trailer") + ylab("sum(SHAPE)") + ggtitle(t)
  fp = c(file.path(save_path_tx_dTCR, paste("tx_", substr(libs[i], 9, (nchar(libs[i])-6)), ".png", sep="")))
  ggsave(file = fp, p_tx)
  
  sev <- normalizeTxUL(gr, sam, norm="TC.control")
  t <- paste("tx_profile, ", substr(libs[i], 9, (nchar(libs[i])-6)), ", TC.control", sep="")
  p_tx <- ggplot() + geom_step(data = sev, aes(start, value), size=0.1) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) +
    xlab("leader, CDS, trailer") + ylab("sum(SHAPE)") + ggtitle(t)
  fp = c(file.path(save_path_tx_TC, paste("tx_", substr(libs[i], 9, (nchar(libs[i])-6)), "_CONTROL.png", sep="")))
  ggsave(file = fp, p_tx)
  
  sev <- normalizeTxUL(gr, sam, norm="TC.treated")
  t <- paste("tx_profile, ", substr(libs[i], 9, (nchar(libs[i])-6)), ", TC.treated", sep="")
  p_tx <- ggplot() + geom_step(data = sev, aes(start, value), size=0.1) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) +
    xlab("leader, CDS, trailer") + ylab("sum(SHAPE)") + ggtitle(t)
  fp = c(file.path(save_path_tx_TC, paste("tx_", substr(libs[i], 9, (nchar(libs[i])-6)), "_TREATED.png", sep="")))
  ggsave(file = fp, p_tx)
  
  sev <- normalizeTxUL(gr, sam, norm="TCR.control")
  t <- paste("tx_profile, ", substr(libs[i], 9, (nchar(libs[i])-6)), ", TCR.control", sep="")
  p_tx <- ggplot() + geom_step(data = sev, aes(start, value), size=0.1) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) +
    xlab("leader, CDS, trailer") + ylab("sum(SHAPE)") + ggtitle(t)
  fp = c(file.path(save_path_tx_TCR, paste("tx_", substr(libs[i], 9, (nchar(libs[i])-6)), "_CONTROL.png", sep="")))
  ggsave(file = fp, p_tx)
  
  sev <- normalizeTxUL(gr, sam, norm="TCR.treated")
  t <- paste("tx_profile, ", substr(libs[i], 9, (nchar(libs[i])-6)), ", TCR.treated", sep="")
  p_tx <- ggplot() + geom_step(data = sev, aes(start, value), size=0.1) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) +
    xlab("leader, CDS, trailer") + ylab("sum(SHAPE)") + ggtitle(t)
  fp = c(file.path(save_path_tx_TCR, paste("tx_", substr(libs[i], 9, (nchar(libs[i])-6)), "_TREATED.png", sep="")))
  ggsave(file = fp, p_tx)
  
}

