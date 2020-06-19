### analyze Sel-A and Sel-C
library(GenomicFeatures)
library(ggplot2)
library(plyr)
library(data.table)
library(reshape2)
library(GeneCycle)
library(seqinr)


libs_path <- "/Volumes/USELESS/DATA/Shape-Seq/shapes_Dec2017"

selA <- get(load(file = c(file.path(libs_path, "PAIR_GC_SelA.Rsave"))))
selC <- get(load(file = c(file.path(libs_path, "PAIR_GC_SelC.Rsave"))))
rm(GR_treated)


###################################################
######### tx_profiles: whole tx to unit length

######### tx_profiles

load(file="/Volumes/USELESS/META/SHAPES_NEW/general/clouds/cloud_df.Rsave")

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

libs <- c("PAIR_GC_SelA.Rsave", "PAIR_GC_SelC.Rsave")

save_path_tx_TC <- "/Volumes/USELESS/META/SHAPES_NEW/PAIRED/tx_profiles/tx_unit_length/TC"
save_path_tx_TCR <- "/Volumes/USELESS/META/SHAPES_NEW/PAIRED/tx_profiles/tx_unit_length/TCR"

for (i in 1:length(libs)) {
  u <- get(load(file = c(file.path(libs_path, libs[i]))))
  u <- unlist(u)
  names(u) <- substring(names(u), 1, 18)
  gr <- GRanges(seqnames=Rle(names(u)),
                IRanges(start(u), width=1),
                strand='*')
  gr$TC.treated <- u$TC.treated
  #gr$TCR.treated <- u$TCR.treated
  gr <- split(gr, seqnames(gr))
  
  sev <- normalizeTxUL(gr, sam, norm="TC.treated")
  t <- paste("tx_profile, ", substr(libs[i], 9, (nchar(libs[i])-6)), ", TC.treated", sep="")
  p_tx <- ggplot() + geom_step(data = sev, aes(start, value), size=0.1) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) +
    xlab("leader, CDS, trailer") + ylab("sum(SHAPE)") + ggtitle(t)
  fp = c(file.path(save_path_tx_TC, paste("tx_", substr(libs[i], 9, (nchar(libs[i])-6)), "_TREATED.png", sep="")))
  ggsave(file = fp, p_tx)
  
  #sev <- normalizeTxUL(gr, sam, norm="TCR.treated")
  #t <- paste("tx_profile, ", substr(libs[i], 9, (nchar(libs[i])-6)), ", TCR.treated", sep="")
  #p_tx <- ggplot() + geom_step(data = sev, aes(start, value), size=0.1) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) +
  #  xlab("leader, CDS, trailer") + ylab("sum(SHAPE)") + ggtitle(t)
  #fp = c(file.path(save_path_tx_TCR, paste("tx_", substr(libs[i], 9, (nchar(libs[i])-6)), "_TREATED.png", sep="")))
  #ggsave(file = fp, p_tx)
  
}


######## RNA

df_selA <- data.frame(selA_TC.treated = sapply(selA, function(x){mean(x$TC.treated)}), n = names(selA))
df_selC <- data.frame(selC_TC.treated = sapply(selC, function(x){mean(x$TC.treated)}), n = names(selC))

df <- merge(df, df_selA, by="n", all=TRUE)
df <- merge(df, df_selC, by="n", all=TRUE)

# RNA
ggplot(df, aes(x = selA_TC.treated, y = rna_4h_total)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave("/Volumes/USELESS/META/SHAPES_NEW/PAIRED/general/shape_vs_rna/selA_sphere.png")
ggplot(df, aes(x = selA_TC.treated, y = rna_3_5h)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggplot(df, aes(x = selA_TC.treated, y = rna_cell1K)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave("/Volumes/USELESS/META/SHAPES_NEW/PAIRED/general/shape_vs_rna/selA.png")

ggplot(df, aes(x = selC_TC.treated, y = rna_4h_total)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave("/Volumes/USELESS/META/SHAPES_NEW/PAIRED/general/shape_vs_rna/selC_sphere.png")
ggplot(df, aes(x = selC_TC.treated, y = rna_3_5h)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggplot(df, aes(x = selC_TC.treated, y = rna_cell1K)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave("/Volumes/USELESS/META/SHAPES_NEW/PAIRED/general/shape_vs_rna/selC.png")

# Ribo
ggplot(df, aes(x = selA_TC.treated, y = ribo_cell1K)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave("/Volumes/USELESS/META/SHAPES_NEW/PAIRED/general/shape_vs_ribo/selA.png")
ggplot(df, aes(x = selC_TC.treated, y = ribo_cell1K)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave("/Volumes/USELESS/META/SHAPES_NEW/PAIRED/general/shape_vs_ribo/selC.png")


### periodicity
txdb_can <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/danio_rerio/Danio_rerio.GRCz10.84.gtf", format = "gtf")
txLengths <- transcriptLengths(txdb_can, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
rownames(txLengths) <- txLengths$tx_name
txLengths <- txLengths[order(rownames(txLengths)),]

cds150 <- GRanges(Rle(rownames(txLengths), rep(1, nrow(txLengths))), IRanges(txLengths$utr5_len+1, width=rep(150, nrow(txLengths))))
cds150 <- split(cds150, seqnames(cds150))

save_path_period <- "/Volumes/USELESS/META/SHAPES_NEW/PAIRED/periodicity"

for (i in 1:length(libs)) {
  u <- get(load(file = c(file.path(libs_path, libs[i]))))
  u <- unlist(u)
  names(u) <- substring(names(u), 1, 18)
  gr <- GRanges(seqnames=Rle(names(u)),
                IRanges(start(u), width=1),
                strand='*')
  gr$TC.treated <- u$TC.treated
  
  f150 <- subsetByOverlaps(gr, cds150, ignore.strand=TRUE)
  f150 <- split(f150, seqnames(f150))
  meta <- f150[sapply(f150, function(x){length(x) > 0})]
  
  metat <- meta[sapply(meta, function(x){sum(x$TC.treated) > 0})]
  metat <- Reduce("+", lapply(metat, function(x){x$TC.treated/sum(x$TC.treated)}))
  
  amplitudes <- abs(fft(metat))
  amp <- amplitudes[2:(length(amplitudes)/2+1)]
  periods <- 1/periodogram(metat, method = "clone")$freq
  df <- data.frame(periods, amp)
  t <- paste("FFT, ", substr(libs[i], 9, (nchar(libs[i])-6)), ", TC.treated", sep="")
  p2 <- ggplot(df, aes(x=periods, y=amp)) + geom_line() + ggtitle(t) + xlim(0, 8)
  fp = c(file.path(save_path_period, paste("t_", substr(libs[i], 9, (nchar(libs[i])-6)), ".png", sep="")))
  ggsave(file = fp, plot = p2)
  
}


###################################################### START ################################################################

# get longest transcript per gene
txlen <- arrange(txLengths, gene_id, desc(tx_len))
txlen <- txlen[!duplicated(txlen$gene_id),]

### tx with utr5 longer than 30
txlen <- txlen[txlen$utr5_len > 100,]
rownames(txlen) <- txlen$tx_name
txlen <- txlen[order(rownames(txlen)),]

start <- GRanges(Rle(rownames(txlen), rep(1, nrow(txlen))), IRanges(txlen$utr5_len-99, width=rep(200, nrow(txlen))))
start <- split(start, seqnames(start))

scale <- c(-100:100)[-101]

###

save_path_start <- "/Volumes/USELESS/META/SHAPES_NEW/PAIRED/start/TC"

for (i in 1:length(libs)) {
  u <- get(load(file = c(file.path(libs_path, libs[i]))))
  u <- unlist(u)
  names(u) <- substring(names(u), 1, 18)
  gr <- GRanges(seqnames=Rle(names(u)),
                IRanges(start(u), width=1),
                strand='*')
  gr$TC.treated <- u$TC.treated
  meta <- subsetByOverlaps(gr, start, ignore.strand=TRUE)
  meta <- split(meta, seqnames(meta))
  meta <- meta[sapply(meta, function(x){length(x) > 0})]
  metad <- meta[sapply(meta, function(x){sum(x$TC.treated) > 0})]
  metad <- Reduce("+", lapply(metad, function(x){x$TC.treated/sum(x$TC.treated)}))
  d <- data.frame(scale = scale, shape = metad)
  
  t <- paste("START, ", substr(libs[i], 9, (nchar(libs[i])-6)), ", TC.treated", sep="")
  p_start <- ggplot(d, aes(x = scale, y = shape)) + geom_bar(stat = "identity") + ggtitle(t)
  fstart = c(file.path(save_path_start, paste("start100_", substr(libs[i], 9, (nchar(libs[i])-6)), ".png", sep="")))
  ggsave(file = fstart, p_start)
}




###################################################### STOP #################################################################
txlen <- txLengths[txLengths$utr3_len > 100,]
stop <- GRanges(Rle(rownames(txlen), rep(1, nrow(txlen))), IRanges(txlen$utr5_len+txlen$cds_len-99, width=rep(200, nrow(txlen))))
stop <- split(stop, seqnames(stop))

scale <- c(-100:100)[-101]

### stop codon
fasta_cds <- read.fasta("/Volumes/USELESS/DATA/fasta/danio_rerio/GRCz10/cds/Danio_rerio.GRCz10.cds.all.fa")
names(fasta_cds) <- sapply(names(fasta_cds), function(x){substr(x,1,18)})
fasta_cds <- fasta_cds[unlist(lapply(fasta_cds, function(x){length(x) %% 3 == 0}))]

taa <- sort(names(fasta_cds[unlist(lapply(fasta_cds, function(x){paste(x[(length(x)-2):length(x)], collapse = "") == "taa"}))]))
tga <- sort(names(fasta_cds[unlist(lapply(fasta_cds, function(x){paste(x[(length(x)-2):length(x)], collapse = "") == "tga"}))]))
tag <- sort(names(fasta_cds[unlist(lapply(fasta_cds, function(x){paste(x[(length(x)-2):length(x)], collapse = "") == "tag"}))]))

###

save_path_stop <- "/Volumes/USELESS/META/SHAPES_NEW/PAIRED/stop/TC"


for (i in 1:length(libs)) {
  u <- get(load(file = c(file.path(libs_path, libs[i]))))
  u <- unlist(u)
  names(u) <- substring(names(u), 1, 18)
  gr <- GRanges(seqnames=Rle(names(u)),
                IRanges(start(u), width=1),
                strand='*')
  gr$TC.treated <- u$TC.treated
  meta <- subsetByOverlaps(gr, stop, ignore.strand=TRUE)
  meta <- split(meta, seqnames(meta))
  meta <- meta[sapply(meta, function(x){length(x) > 0})]
  metad <- meta[sapply(meta, function(x){sum(x$TC.treated) > 0})]
  
  shape_taa <- metad[names(metad) %in% taa]
  shape_tga <- metad[names(metad) %in% tga]
  shape_tag <- metad[names(metad) %in% tag]
  
  df_shape <- data.frame(scale = scale,
                         shape_taa = Reduce("+", lapply(shape_taa, function(x){x$TC.treated/sum(x$TC.treated)})),
                         shape_tga = Reduce("+", lapply(shape_tga, function(x){x$TC.treated/sum(x$TC.treated)})),
                         shape_tag = Reduce("+", lapply(shape_tag, function(x){x$TC.treated/sum(x$TC.treated)})))
  
  t <- paste("STOP, ", substr(libs[i], 9, (nchar(libs[i])-6)), ", TC.treated, TAA", sep="")
  p_taa <- ggplot(df_shape, aes(x=scale, y=shape_taa)) + geom_bar(stat = "identity") + ggtitle(t)
  fstop_taa = c(file.path(save_path_stop, paste("stop100_", substr(libs[i], 9, (nchar(libs[i])-6)), "_TAA.png", sep="")))
  ggsave(file = fstop_taa, p_taa)
  
  t <- paste("STOP, ", substr(libs[i], 9, (nchar(libs[i])-6)), ", TC.treated, TGA", sep="")
  p_tga <- ggplot(df_shape, aes(x=scale, y=shape_tga)) + geom_bar(stat = "identity") + ggtitle(t)
  fstop_tga = c(file.path(save_path_stop, paste("stop100_", substr(libs[i], 9, (nchar(libs[i])-6)), "_TGA.png", sep="")))
  ggsave(file = fstop_tga, p_tga)
  
  t <- paste("STOP, ", substr(libs[i], 9, (nchar(libs[i])-6)), ", TC.treated, TAG", sep="")
  p_tag <- ggplot(df_shape, aes(x=scale, y=shape_tag)) + geom_bar(stat = "identity") + ggtitle(t)
  fstop_tag = c(file.path(save_path_stop, paste("stop100_", substr(libs[i], 9, (nchar(libs[i])-6)), "_TAG.png", sep="")))
  ggsave(file = fstop_tag, p_tag)
}


for (i in 1:length(libs)) {
  u <- get(load(file = c(file.path(libs_path, libs[i]))))
  u <- unlist(u)
  names(u) <- substring(names(u), 1, 18)
  gr <- GRanges(seqnames=Rle(names(u)),
                IRanges(start(u), width=1),
                strand='*')
  gr$TC.treated <- u$TC.treated
  meta <- subsetByOverlaps(gr, stop, ignore.strand=TRUE)
  meta <- split(meta, seqnames(meta))
  meta <- meta[sapply(meta, function(x){length(x) > 0})]
  metad <- meta[sapply(meta, function(x){sum(x$TC.treated) > 0})]
  shape <- metad
  
  df_shape <- data.frame(scale = scale,
                         shape = Reduce("+", lapply(shape, function(x){x$TC.treated/sum(x$TC.treated)})))
  
  t <- paste("STOP, ", substr(libs[i], 9, (nchar(libs[i])-6)), ", TC.treated", sep="")
  p_taa <- ggplot(df_shape, aes(x=scale, y=shape)) + geom_bar(stat = "identity") + ggtitle(t)
  fstop_taa = c(file.path(save_path_stop, paste("stop100_", substr(libs[i], 9, (nchar(libs[i])-6)), ".png", sep="")))
  ggsave(file = fstop_taa, p_taa)
}


################ tx_profiles: utr5, cds, utr3 #######################
# get longest transcript per gene
txlen <- arrange(txLengths, gene_id, desc(tx_len))
txlen <- txlen[!duplicated(txlen$gene_id),]

lct <- txlen[(txlen$utr5_len > 100) & (txlen$cds_len > 100) & (txlen$utr3_len > 100), ]

load(file="/Volumes/USELESS/META/SHAPES_NEW/general/clouds/cloud_df.Rsave")
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


save_path_tx_TC <- "/Volumes/USELESS/META/SHAPES_NEW/PAIRED/tx_profiles/utr5_cds_utr3/TC"
save_path_tx_TCR <- "/Volumes/USELESS/META/SHAPES_NEW/PAIRED/tx_profiles/utr5_cds_utr3/TCR"

for (i in 1:length(libs)) {
  u <- get(load(file = c(file.path(libs_path, libs[i]))))
  u <- unlist(u)
  names(u) <- substring(names(u), 1, 18)
  gr <- GRanges(seqnames=Rle(names(u)),
                IRanges(start(u), width=1),
                strand='*')
  gr$TC.treated <- u$TC.treated
  #gr$TCR.treated <- u$TCR.treated
  gr <- split(gr, seqnames(gr))
  
  sev <- normalizeTx(gr, txlen, sam, norm="TC.treated")
  t <- paste("tx_profile, ", substr(libs[i], 9, (nchar(libs[i])-6)), ", TC.treated", sep="")
  p_tx <- ggplot() + geom_step(data = sev, aes(start, value), size=0.1) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) +
    xlab("leader, CDS, trailer") + ylab("sum(SHAPE)") + ggtitle(t)
  fp = c(file.path(save_path_tx_TC, paste("tx_", substr(libs[i], 9, (nchar(libs[i])-6)), "_TREATED.png", sep="")))
  ggsave(file = fp, p_tx)
  
  #sev <- normalizeTx(gr, txlen, sam, norm="TCR.treated")
  #t <- paste("tx_profile, ", substr(libs[i], 9, (nchar(libs[i])-6)), ", TCR.treated", sep="")
  #p_tx <- ggplot() + geom_step(data = sev, aes(start, value), size=0.1) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) +
  #  xlab("leader, CDS, trailer") + ylab("sum(SHAPE)") + ggtitle(t)
  #fp = c(file.path(save_path_tx_TCR, paste("tx_", substr(libs[i], 9, (nchar(libs[i])-6)), "_TREATED.png", sep="")))
  #ggsave(file = fp, p_tx)
  
}



################ cloud

#selA
ggplot(df, aes(x = log2(selA_TC.treated), y = log2(rna_4h_total))) + geom_point(size=0.1)
f <- function(x) {0.7*x + 3}
ggplot(df, aes(x = log2(selA_TC.treated), y = log2(rna_4h_total))) + geom_point(size=0.1) + geom_point(size=0.1) + stat_function(fun = f, colour = "red")

piggyback <- df$n[log2(df$rna_4h_total) > (0.7*log2(df$selA_TC.treated) + 3) & log2(df$selA_TC.treated) < -5 & log2(df$selA_TC.treated) > -13]
piggyback <- piggyback[!is.na(piggyback)]
piggyback <- data.frame(n = piggyback, cloud_selA = rep("cloud", length(piggyback)))
df <- merge(df, piggyback, by="n", all=TRUE)

#selC
ggplot(df, aes(x = log2(selC_TC.treated), y = log2(rna_4h_total))) + geom_point(size=0.1)
f <- function(x) {0.7*x + 3}
ggplot(df, aes(x = log2(selC_TC.treated), y = log2(rna_4h_total))) + geom_point(size=0.1) + geom_point(size=0.1) + stat_function(fun = f, colour = "red")

piggyback <- df$n[log2(df$rna_4h_total) > (0.7*log2(df$selC_TC.treated) + 3) & log2(df$selC_TC.treated) < -5 & log2(df$selC_TC.treated) > -13]
piggyback <- piggyback[!is.na(piggyback)]
piggyback <- data.frame(n = piggyback, cloud_selC = rep("cloud", length(piggyback)))
df <- merge(df, piggyback, by="n", all=TRUE)


cA <- as.character(df[!is.na(df$cloud_selA),]$gene_id)
cC <- as.character(df[!is.na(df$cloud_selC),]$gene_id)
csel <- unique(cA[cA %in% cC])

fileConn <- file("/Users/kasia/Desktop/cloud_sel.txt")    
writeLines(csel, fileConn)    
close(fileConn)

all <- unique(as.character(df[!is.na(df$selA_TC.treated),]$gene_id))

fileConn <- file("/Users/kasia/Desktop/all_sel.txt")    
writeLines(all, fileConn)    
close(fileConn)

## add gene_id to df
dupa <- data.frame(n = txLengths$tx_name, gene_id = txLengths$gene_id)
df <- merge(df, dupa, by="n", all=TRUE)



