### ANALYZE 1Kcell
library(ggplot2)
library(GenomicFeatures)
library(GeneCycle)
library(plyr)
library(seqinr)
library(reshape2)
library(data.table)

libs_path <- "/Volumes/USELESS/DATA/Shape-Seq/1Kcell"
## log2ratio
single_RZ <- get(load(file = c(file.path(libs_path, "GC_1Kcell_RZ.Rsave"))))
single_RZ_vitro <- get(load(file = c(file.path(libs_path, "GC_1Kcell_RZ_vitro.Rsave"))))
single_polyA <- get(load(file = c(file.path(libs_path, "GC_1Kcell_polyA.Rsave"))))
single_polyA_vitro <- get(load(file = c(file.path(libs_path, "GC_1Kcell_polyA_vitro.Rsave"))))


## dTCR
paired_RZ <- get(load(file = c(file.path(libs_path, "PAIR_GC_1Kcell_RZ.Rsave"))))
paired_RZ_vitro <- get(load(file = c(file.path(libs_path, "PAIR_GC_1Kcell_RZ_vitro.Rsave"))))
paired_polyA <- get(load(file = c(file.path(libs_path, "PAIR_GC_1Kcell_polyA.Rsave"))))
paired_polyA_vitro <- get(load(file = c(file.path(libs_path, "PAIR_GC_1Kcell_polyA_vitro.Rsave"))))

rm(GR)


### general: 
# shape vs RNA;
# shape vs Ribo;
# periodicity;
# START;
# STOP;
# tx profiles

# load(file="/Volumes/USELESS/META/SHAPES_NEW/general/clouds/cloud_df.Rsave")

single_RZ$log2ratio[single_RZ$log2ratio < 0] <- 0
single_RZ_vitro$log2ratio[single_RZ_vitro$log2ratio < 0] <- 0
single_polyA$log2ratio[single_polyA$log2ratio < 0] <- 0
single_polyA_vitro$log2ratio[single_polyA_vitro$log2ratio < 0] <- 0

single_RZ <- split(single_RZ, names(single_RZ))
single_RZ_vitro <- split(single_RZ_vitro, names(single_RZ_vitro))
single_polyA <- split(single_polyA, names(single_polyA))
single_polyA_vitro <- split(single_polyA_vitro, names(single_polyA_vitro))

#dtcr <- (paired_RZ$TCR.treated - paired_RZ$TCR.control) / (1 - paired_RZ$TCR.control)
#dtcr[dtcr < 0] <- 0
#paired_RZ$dTCR <- dtcr
#paired_RZ <- split(paired_RZ, names(paired_RZ))

###################################################### SHAPE vs RNA #########################################################
###################################################### SHAPE vs Ribo ########################################################

df_single_cell1K_RZ <- data.frame(shape_cell1K_RZ_log2ratio_internal = sapply(single_RZ, function(x){mean(x$log2ratio[2:length(x)])}),
                             n = names(single_RZ))
df_single_cell1K_RZ_vitro <- data.frame(shape_cell1K_RZ_vitro_log2ratio_internal = sapply(single_RZ_vitro,
                                                                                          function(x){mean(x$log2ratio[2:length(x)])}),
                                        n = names(single_RZ_vitro))
df_single_cell1K_polyA <- data.frame(shape_cell1K_polyA_log2ratio_internal = sapply(single_polyA,
                                                                                    function(x){mean(x$log2ratio[2:length(x)])}),
                                  n = names(single_polyA))
df_single_cell1K_polyA_vitro <- data.frame(shape_cell1K_polyA_vitro_log2ratio_internal = sapply(single_polyA_vitro,
                                                                                    function(x){mean(x$log2ratio[2:length(x)])}),
                                     n = names(single_polyA_vitro))


df <- merge(df, df_single_cell1K_RZ, by="n", all=TRUE)
df <- merge(df, df_single_cell1K_RZ_vitro, by="n", all=TRUE)
df <- merge(df, df_single_cell1K_polyA, by="n", all=TRUE)
df <- merge(df, df_single_cell1K_polyA_vitro, by="n", all=TRUE)

#ggplot(df, aes(x = log2(shape_cell1K_RZ_log2ratio_internal), y = log2(rna_cell1K))) + geom_point(size=0.1)
#ggplot(df, aes(x = log2(shape_cell1K_RZ_log2ratio_internal), y = log2(ribo_cell1K))) + geom_point(size=0.1)

Xaes = c("shape_cell1K_RZ_log2ratio_internal", "shape_cell1K_RZ_vitro_log2ratio_internal", "shape_cell1K_polyA_log2ratio_internal",
         "shape_cell1K_polyA_vitro_log2ratio_internal")


for (i in 1:4) {
  ggplot(df, aes(x = log2(df[Xaes[i]]), y = log2(df$rna_cell1K))) + geom_point(size=0.1) + xlab(Xaes[i]) + ylab("log2(rna_cell1K)")
  fp <- paste0("/Volumes/USELESS/META/SHAPES_NEW/1Kcell/general/shape_rna", i, ".png", collapse = NULL)
  ggsave(fp)
}

for (i in 1:4) {
  ggplot(df, aes(x = log2(df[Xaes[i]]), y = log2(df$ribo_cell1K))) + geom_point(size=0.1) + xlab(Xaes[i]) + ylab("log2(ribo_cell1K)")
  fp <- paste0("/Volumes/USELESS/META/SHAPES_NEW/1Kcell/general/shape_ribo", i, ".png", collapse = NULL)
  ggsave(fp)
}

### mini-cloud on ribo
ggplot(df, aes(x = log2(shape_cell1K_polyA_log2ratio_internal), y = log2(ribo_cell1K))) + geom_point(size=0.1)
f <- function(x) {1.3*x + 15}
ggplot(df, aes(x = log2(shape_cell1K_polyA_log2ratio_internal), y = log2(ribo_cell1K))) + geom_point(size=0.1) +
  stat_function(fun = f, colour = "red")

piggyback <- df$n[log2(df$ribo_cell1K) > (1.3*log2(df$shape_cell1K_polyA_log2ratio_internal) + 15) & log2(df$shape_cell1K_polyA_log2ratio_internal) > -10 & log2(df$ribo_cell1K) > 5]
piggyback <- piggyback[!is.na(piggyback)]
piggyback <- data.frame(n = piggyback, cloud_cell1K_polyA_ribo = rep("cloud", length(piggyback)))
df <- merge(df, piggyback, by="n", all=TRUE)



####### paired
paired_RZ <- split(paired_RZ, names(paired_RZ))
paired_RZ_vitro <- split(paired_RZ_vitro, names(paired_RZ_vitro))
paired_polyA <- split(paired_polyA, names(paired_polyA))
paired_polyA_vitro <- split(paired_polyA_vitro, names(paired_polyA_vitro))

df_paired_cell1K_RZ <- data.frame(shape_cell1K_RZ_dTCR_internal = sapply(paired_RZ, function(x){mean(x$dTCR[2:length(x)])}),
                                  n = names(paired_RZ))
df_paired_cell1K_RZ_vitro <- data.frame(shape_cell1K_RZ_vitro_dTCR_internal = sapply(paired_RZ_vitro,
                                                                                          function(x){mean(x$dTCR[2:length(x)])}),
                                        n = names(paired_RZ_vitro))
df_paired_cell1K_polyA <- data.frame(shape_cell1K_polyA_dTCR_internal = sapply(paired_polyA,
                                                                                    function(x){mean(x$dTCR[2:length(x)])}),
                                     n = names(paired_polyA))
df_paired_cell1K_polyA_vitro <- data.frame(shape_cell1K_polyA_vitro_dTCR_internal = sapply(paired_polyA_vitro,
                                                                                                function(x){mean(x$dTCR[2:length(x)])}),
                                           n = names(paired_polyA_vitro))

# load(file="/Volumes/USELESS/META/SHAPES_NEW/general/clouds/cloud_df.Rsave")
df <- merge(df, df_paired_cell1K_RZ, by="n", all=TRUE)
df <- merge(df, df_paired_cell1K_RZ_vitro, by="n", all=TRUE)
df <- merge(df, df_paired_cell1K_polyA, by="n", all=TRUE)
df <- merge(df, df_paired_cell1K_polyA_vitro, by="n", all=TRUE)

# save(df, file="/Volumes/USELESS/META/SHAPES_NEW/general/clouds/cloud_df.Rsave")

###################################################### periodicity ##########################################################
txdb_can <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/danio_rerio/Danio_rerio.GRCz10.84.gtf", format = "gtf")
txLengths <- transcriptLengths(txdb_can, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
rownames(txLengths) <- txLengths$tx_name
txLengths <- txLengths[order(rownames(txLengths)),]

cds150 <- GRanges(Rle(rownames(txLengths), rep(1, nrow(txLengths))), IRanges(txLengths$utr5_len+1, width=rep(150, nrow(txLengths))))
cds150 <- split(cds150, seqnames(cds150))

###
save_path <- "/Volumes/USELESS/META/SHAPES_NEW/1Kcell/periodicity"
libs <- c("GC_1Kcell_RZ.Rsave", "GC_1Kcell_RZ_vitro.Rsave", "GC_1Kcell_polyA.Rsave", "GC_1Kcell_polyA_vitro.Rsave")

for (i in 1:length(libs)) {
  u <- get(load(file = c(file.path(libs_path, libs[i]))))
  u$log2ratio[u$log2ratio < 0] <- 0
  
  gr <- GRanges(seqnames=Rle(names(u)),
                IRanges(start(u), width=1),
                strand='*')
  gr$TC.control <- u$TC.control
  gr$TC.treated <- u$TC.treated
  gr$log2ratio <- u$log2ratio
  
  f150 <- subsetByOverlaps(gr, cds150, ignore.strand=TRUE)
  f150 <- split(f150, seqnames(f150))
  meta <- f150[sapply(f150, function(x){length(x) > 0})]
  
  metac <- meta[sapply(meta, function(x){sum(x$TC.control) > 0})]
  metat <- meta[sapply(meta, function(x){sum(x$TC.treated) > 0})]
  metal <- meta[sapply(meta, function(x){sum(x$log2ratio) > 0})]
  
  metac <- Reduce("+", lapply(metac, function(x){x$TC.control/sum(x$TC.control)}))
  metat <- Reduce("+", lapply(metat, function(x){x$TC.treated/sum(x$TC.treated)}))
  metal <- Reduce("+", lapply(metal, function(x){x$log2ratio/sum(x$log2ratio)}))
  
  amplitudes <- abs(fft(metac))
  amp <- amplitudes[2:(length(amplitudes)/2+1)]
  periods <- 1/periodogram(metac, method = "clone")$freq
  df <- data.frame(periods, amp)
  t <- paste("FFT, ", substr(libs[i], 4, (nchar(libs[i])-6)), ", TC.control", sep="")
  p1 <- ggplot(df, aes(x=periods, y=amp)) + geom_line() + ggtitle(t) + xlim(0, 8)
  fp = c(file.path(save_path, paste("c_", substr(libs[i], 4, (nchar(libs[i])-6)), ".png", sep="")))
  ggsave(file = fp, plot = p1)
  
  amplitudes <- abs(fft(metat))
  amp <- amplitudes[2:(length(amplitudes)/2+1)]
  periods <- 1/periodogram(metat, method = "clone")$freq
  df <- data.frame(periods, amp)
  t <- paste("FFT, ", substr(libs[i], 4, (nchar(libs[i])-6)), ", TC.treated", sep="")
  p2 <- ggplot(df, aes(x=periods, y=amp)) + geom_line() + ggtitle(t) + xlim(0, 8)
  fp = c(file.path(save_path, paste("t_", substr(libs[i], 4, (nchar(libs[i])-6)), ".png", sep="")))
  ggsave(file = fp, plot = p2)
  
  amplitudes <- abs(fft(metal))
  amp <- amplitudes[2:(length(amplitudes)/2+1)]
  periods <- 1/periodogram(metal, method = "clone")$freq
  df <- data.frame(periods, amp)
  t <- paste("FFT, ", substr(libs[i], 4, (nchar(libs[i])-6)), ", log2ratio", sep="")
  p3 <- ggplot(df, aes(x=periods, y=amp)) + geom_line() + ggtitle(t) + xlim(0, 8)
  fp = c(file.path(save_path, paste("l_", substr(libs[i], 4, (nchar(libs[i])-6)), ".png", sep="")))
  ggsave(file = fp, plot = p3)

}

## dtcr
libs <- c("PAIR_GC_1Kcell_RZ.Rsave", "PAIR_GC_1Kcell_RZ_vitro.Rsave", "PAIR_GC_1Kcell_polyA.Rsave",
                 "PAIR_GC_1Kcell_polyA_vitro.Rsave")

for (i in 1:length(libs)) {
  u <- get(load(file = c(file.path(libs_path, libs[i]))))
  
  gr <- GRanges(seqnames=Rle(names(u)),
                IRanges(start(u), width=1),
                strand='*')
  gr$dTCR <- u$dTCR
  
  f150 <- subsetByOverlaps(gr, cds150, ignore.strand=TRUE)
  f150 <- split(f150, seqnames(f150))
  meta <- f150[sapply(f150, function(x){length(x) == 150})]
  m <- unlist(meta)
  m$dTCR[is.na(m$dTCR)] <- 0
  m$dTCR[m$dTCR == Inf] <- 0
  meta <- split(m, seqnames(m))
  
  metad <- meta[sapply(meta, function(x){sum(x$dTCR) > 0})]
  metad <- Reduce("+", lapply(metad, function(x){x$dTCR/sum(x$dTCR)}))
  
  amplitudes <- abs(fft(metad))
  amp <- amplitudes[2:(length(amplitudes)/2+1)]
  periods <- 1/periodogram(metad, method = "clone")$freq
  df <- data.frame(periods, amp)
  t <- paste("FFT, ", substr(libs[i], 9, (nchar(libs[i])-6)), ", dTCR", sep="")
  p1 <- ggplot(df, aes(x=periods, y=amp)) + geom_line() + ggtitle(t) + xlim(0, 8)
  fp = c(file.path(save_path, paste("d_", substr(libs[i], 9, (nchar(libs[i])-6)), ".png", sep="")))
  ggsave(file = fp, plot = p1)
  
}

###################################################### START ################################################################

# get longest transcript per gene
txlen <- arrange(txLengths, gene_id, desc(tx_len))
txlen <- txlen[!duplicated(txlen$gene_id),]

### tx with utr5 longer than 30
txlen <- txlen[txlen$utr5_len > 30,]
rownames(txlen) <- txlen$tx_name
txlen <- txlen[order(rownames(txlen)),]

start <- GRanges(Rle(rownames(txlen), rep(1, nrow(txlen))), IRanges(txlen$utr5_len-29, width=rep(60, nrow(txlen))))
start <- split(start, seqnames(start))

scale <- c(-30:30)[-31]

###
libs <- c("PAIR_GC_1Kcell_RZ.Rsave", "PAIR_GC_1Kcell_RZ_vitro.Rsave", "PAIR_GC_1Kcell_polyA.Rsave",
          "PAIR_GC_1Kcell_polyA_vitro.Rsave")

save_path_start <- "/Volumes/USELESS/META/SHAPES_NEW/1Kcell/start"

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

###################################################### STOP #################################################################
txLengths <- txLengths[txLengths$utr3_len > 30,]
stop <- GRanges(Rle(rownames(txLengths), rep(1, nrow(txLengths))), IRanges(txLengths$utr5_len+txLengths$cds_len-29, width=rep(60, nrow(txLengths))))
stop <- split(stop, seqnames(stop))

scale <- c(-30:30)[-31]

### stop codon
fasta_cds <- read.fasta("/Volumes/USELESS/DATA/fasta/danio_rerio/GRCz10/cds/Danio_rerio.GRCz10.cds.all.fa")
names(fasta_cds) <- sapply(names(fasta_cds), function(x){substr(x,1,18)})
fasta_cds <- fasta_cds[unlist(lapply(fasta_cds, function(x){length(x) %% 3 == 0}))]

taa <- sort(names(fasta_cds[unlist(lapply(fasta_cds, function(x){paste(x[(length(x)-2):length(x)], collapse = "") == "taa"}))]))
tga <- sort(names(fasta_cds[unlist(lapply(fasta_cds, function(x){paste(x[(length(x)-2):length(x)], collapse = "") == "tga"}))]))
tag <- sort(names(fasta_cds[unlist(lapply(fasta_cds, function(x){paste(x[(length(x)-2):length(x)], collapse = "") == "tag"}))]))

###
libs <- c("PAIR_GC_1Kcell_RZ.Rsave", "PAIR_GC_1Kcell_RZ_vitro.Rsave", "PAIR_GC_1Kcell_polyA.Rsave",
          "PAIR_GC_1Kcell_polyA_vitro.Rsave")

save_path_stop <- "/Volumes/USELESS/META/SHAPES_NEW/1Kcell/stop"


for (i in 1:length(libs)) {
  u <- get(load(file = c(file.path(libs_path, libs[i]))))
  gr <- GRanges(seqnames=Rle(names(u)),
                IRanges(start(u), width=1),
                strand='*')
  gr$dTCR <- u$dTCR
  meta <- subsetByOverlaps(gr, stop, ignore.strand=TRUE)
  meta <- split(meta, seqnames(meta))
  meta <- meta[sapply(meta, function(x){length(x) > 0})]
  metad <- meta[sapply(meta, function(x){sum(x$dTCR) > 0})]

  shape_taa <- metad[names(metad) %in% taa]
  shape_tga <- metad[names(metad) %in% tga]
  shape_tag <- metad[names(metad) %in% tag]

  df_shape <- data.frame(scale = scale,
                        shape_taa = Reduce("+", lapply(shape_taa, function(x){x$dTCR/sum(x$dTCR)})),
                        shape_tga = Reduce("+", lapply(shape_tga, function(x){x$dTCR/sum(x$dTCR)})),
                        shape_tag = Reduce("+", lapply(shape_tag, function(x){x$dTCR/sum(x$dTCR)})))

  t <- paste("STOP, ", substr(libs[i], 9, (nchar(libs[i])-6)), ", dTCR, TAA", sep="")
  p_taa <- ggplot(df_shape, aes(x=scale, y=shape_taa)) + geom_bar(stat = "identity") + ggtitle(t)
  fstop_taa = c(file.path(save_path_stop, paste("start_", substr(libs[i], 9, (nchar(libs[i])-6)), "_TAA.png", sep="")))
  ggsave(file = fstop_taa, p_taa)
  
  t <- paste("STOP, ", substr(libs[i], 9, (nchar(libs[i])-6)), ", dTCR, TGA", sep="")
  p_tga <- ggplot(df_shape, aes(x=scale, y=shape_tga)) + geom_bar(stat = "identity") + ggtitle(t)
  fstop_tga = c(file.path(save_path_stop, paste("start_", substr(libs[i], 9, (nchar(libs[i])-6)), "_TGA.png", sep="")))
  ggsave(file = fstop_tga, p_tga)
  
  t <- paste("STOP, ", substr(libs[i], 9, (nchar(libs[i])-6)), ", dTCR, TAG", sep="")
  p_tag <- ggplot(df_shape, aes(x=scale, y=shape_tag)) + geom_bar(stat = "identity") + ggtitle(t)
  fstop_tag = c(file.path(save_path_stop, paste("start_", substr(libs[i], 9, (nchar(libs[i])-6)), "_TAG.png", sep="")))
  ggsave(file = fstop_tag, p_tag)
}

###################################################### tx profiles ##########################################################
### get transcript database
txdb_can <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/danio_rerio/Danio_rerio.GRCz10.84.gtf", format = "gtf")
txLengths <- transcriptLengths(txdb_can, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
rownames(txLengths) <- txLengths$tx_name
txLengths <- txLengths[order(rownames(txLengths)),]

# get longest transcript per gene
txlen <- arrange(txLengths, gene_id, desc(tx_len))
txlen <- txlen[!duplicated(txlen$gene_id),]

exons <- exonsBy(txdb_can, by="tx", use.names=TRUE)
exons <- exons[order(names(exons))]
names(exons) <- sapply(names(exons), function(x){substr(x,1,18)})


# pick tx with utr5, cds, utr3
# subset on rna_fpkm (or SHAPE?)
lct <- txlen[(txlen$utr5_len > 30) & (txlen$cds_len > 100) & (txlen$utr3_len > 30), ]

lct <- lct[lct$tx_name %in% df[df$rna_cell1K > 1,]$n,]
lct <- lct[lct$tx_name %in% df[df$shape_cell1K_RZ_dTCR_internal > 0.01 
                               & df$shape_cell1K_RZ_vitro_dTCR_internal > 0.01
                               & df$shape_cell1K_polyA_dTCR_internal > 0.01
                               & df$shape_cell1K_polyA_vitro_dTCR_internal > 0.01,]$n,]

sam <- sample(lct$tx_name, 500)

#### normalizeTx (from tx_profiles.R)
normalizeTx <- function(shape, txlen, sam){
  # subset shape GRanges
  shape <- shape[names(shape) %in% sam]
  shape <- shape[order(names(shape))]
  len <- txlen[txlen$tx_name %in% names(shape),]
  len <- len[order(len$tx_name),]
  
  leader <- mapply(function(x,y){x$dTCR[1:y]}, shape, len$utr5_len)
  shape <- mapply(function(x,y){x$dTCR[(y+1):length(x)]}, shape, len$utr5_len)
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

libs <- c("PAIR_GC_1Kcell_RZ.Rsave", "PAIR_GC_1Kcell_RZ_vitro.Rsave", "PAIR_GC_1Kcell_polyA.Rsave",
          "PAIR_GC_1Kcell_polyA_vitro.Rsave")

save_path_tx <- "/Volumes/USELESS/META/SHAPES_NEW/1Kcell/tx_profiles"

for (i in 1:length(libs)) {
  u <- get(load(file = c(file.path(libs_path, libs[i]))))
  gr <- GRanges(seqnames=Rle(names(u)),
                IRanges(start(u), width=1),
                strand='*')
  gr$dTCR <- u$dTCR
  gr <- split(gr, seqnames(gr))

  sev <- normalizeTx(gr, txlen, sam)
  
  t <- paste("tx_profile, ", substr(libs[i], 9, (nchar(libs[i])-6)), ", dTCR", sep="")
  p_tx <- ggplot() + geom_step(data = sev, aes(start, value), size=0.1) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) +
    xlab("leader, CDS, trailer") + ylab("sum(SHAPE)") + ggtitle(t)
  fp = c(file.path(save_path_tx, paste("tx_", substr(libs[i], 9, (nchar(libs[i])-6)), ".png", sep="")))
  ggsave(file = fp, p_tx)
}

