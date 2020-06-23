### analyze all shape(s) libraries
library(GenomicFeatures)
library(ggplot2)
library(plyr)
library(data.table)
library(reshape2)
library(GeneCycle)
library(seqinr)

libs_path <- "/Volumes/USELESS/DATA/SHAPES/SE"
libs <- c("norm2_cell1_invitro_NAIN3.Rsave", "norm2_sphere_invitro_A.Rsave", "norm2_sphere_invitro_C.Rsave",
          "norm2_cell24_invivo.Rsave", "norm2_cell256_invivo.Rsave", "norm2_oblong_CHX_invivo.Rsave", "norm2_oblong_invivo.Rsave")

save_path <- "/Volumes/USELESS/META/SHAPES/SE"

## RNA and Ribo, TE, accessibility
## periodicity
## tx profile: unit length
## tx profile: utr5, cds, utr3
## start
## stop
## stop per stop codon
## splice sites (exons)
## shape per codon
## ---------------
## stability (RNA level before/after MZT)
## highly / lowly translated


############################################################################################################################
load(file = "/Volumes/USELESS/META/SHAPES/matrix_of_everything.Rdata")

### RNA
ggplot(matrix_of_everything, aes(x = log2(shapes_cell1_invitro_TC), y = log2(rna_cell1))) + geom_point(size=0.1)
ggsave(file = "/Volumes/USELESS/META/SHAPES/SE/general/shape_vs_rna/cell1.png")

ggplot(matrix_of_everything, aes(x = log2(shape_cell24_invivo_log2ratio), y = log2(rna_cell24))) + geom_point(size=0.1)
ggsave(file = "/Volumes/USELESS/META/SHAPES/SE/general/shape_vs_rna/cell24.png")

ggplot(matrix_of_everything, aes(x = log2(shape_cell256_invivo_log2ratio), y = log2(rna_cell256))) + geom_point(size=0.1)
ggsave(file = "/Volumes/USELESS/META/SHAPES/SE/general/shape_vs_rna/cell256.png")

ggplot(matrix_of_everything, aes(x = log2(shape_oblong_invivo_log2ratio), y = log2(rna_3_5h))) + geom_point(size=0.1)
ggsave(file = "/Volumes/USELESS/META/SHAPES/SE/general/shape_vs_rna/oblong.png")

ggplot(matrix_of_everything, aes(x = log2(shape_oblong_CHX_invivo_log2ratio), y = log2(rna_3_5h))) + geom_point(size=0.1)
ggsave(file = "/Volumes/USELESS/META/SHAPES/SE/general/shape_vs_rna/oblong_CHX.png")

ggplot(matrix_of_everything, aes(x = log2(shapes_sphere_invitro_A_TC), y = log2(rna_4h_total))) + geom_point(size=0.1)
ggsave(file = "/Volumes/USELESS/META/SHAPES/SE/general/shape_vs_rna/sphere_A.png")

ggplot(matrix_of_everything, aes(x = log2(shapes_sphere_invitro_C_TC), y = log2(rna_4h_total))) + geom_point(size=0.1)
ggsave(file = "/Volumes/USELESS/META/SHAPES/SE/general/shape_vs_rna/sphere_C.png")

### Ribo
ggplot(matrix_of_everything, aes(x = log2(shape_cell24_invivo_log2ratio), y = log2(ribo_cell24))) + geom_point(size=0.1)
ggsave(file = "/Volumes/USELESS/META/SHAPES/SE/general/shape_vs_ribo/cell24.png")

ggplot(matrix_of_everything, aes(x = log2(shape_cell256_invivo_log2ratio), y = log2(ribo_cell256))) + geom_point(size=0.1)
ggsave(file = "/Volumes/USELESS/META/SHAPES/SE/general/shape_vs_ribo/cell256.png")

ggplot(matrix_of_everything, aes(x = log2(shapes_sphere_invitro_A_TC), y = log2(ribo_4h))) + geom_point(size=0.1)
ggsave(file = "/Volumes/USELESS/META/SHAPES/SE/general/shape_vs_ribo/sphere_A.png")

ggplot(matrix_of_everything, aes(x = log2(shapes_sphere_invitro_C_TC), y = log2(ribo_4h))) + geom_point(size=0.1)
ggsave(file = "/Volumes/USELESS/META/SHAPES/SE/general/shape_vs_ribo/sphere_C.png")


### TE
ggplot(matrix_of_everything, aes(x = log2(shape_cell24_invivo_log2ratio), y = log2(ribo_cell24/rna_cell24))) + geom_point(size=0.1)
ggsave(file = "/Volumes/USELESS/META/SHAPES/SE/general/shape_vs_te/cell24.png")

ggplot(matrix_of_everything, aes(x = log2(shape_cell256_invivo_log2ratio), y = log2(ribo_cell256/rna_cell256))) + geom_point(size=0.1)
ggsave(file = "/Volumes/USELESS/META/SHAPES/SE/general/shape_vs_te/cell256.png")

ggplot(matrix_of_everything, aes(x = log2(shapes_sphere_invitro_A_TC), y = log2(ribo_4h/rna_4h_total))) + geom_point(size=0.1)
ggsave(file = "/Volumes/USELESS/META/SHAPES/SE/general/shape_vs_te/sphere_A.png")

ggplot(matrix_of_everything, aes(x = log2(shapes_sphere_invitro_C_TC), y = log2(ribo_4h/rna_4h_total))) + geom_point(size=0.1)
ggsave(file = "/Volumes/USELESS/META/SHAPES/SE/general/shape_vs_te/sphere_C.png")


####################

txdb_can <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/danio_rerio/Danio_rerio.GRCz10.84.gtf", format = "gtf")
txLengths <- transcriptLengths(txdb_can, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
rownames(txLengths) <- txLengths$tx_name
txLengths <- txLengths[order(rownames(txLengths)),]

# get longest transcript per gene
txlen <- arrange(txLengths, gene_id, desc(tx_len))
txlen <- txlen[!duplicated(txlen$gene_id),]
rownames(txlen) <- txlen$tx_name

cds150 <- GRanges(Rle(rownames(txlen), rep(1, nrow(txlen))), IRanges(txlen$utr5_len+1, width=rep(150, nrow(txlen))))
cds150 <- split(cds150, seqnames(cds150))

save_path_period <- "/Volumes/USELESS/META/SHAPES/SE/periodicity"

for (i in 1:length(libs)) {
  u <- get(load(file = c(file.path(libs_path, libs[i]))))
  u <- unlist(u)
  names(u) <- substring(names(u), 1, 18)
  gr <- GRanges(seqnames=Rle(names(u)),
                IRanges(start(u), width=1),
                strand='*')
  gr$TC.treated <- u$TC.treated
  if (i > 3) {
    gr$TC.control <- u$TC.control
    gr$log2ratio <- u$log2ratio
  }

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
  p <- ggplot(df, aes(x=periods, y=amp)) + geom_line() + ggtitle(t) + xlim(0, 8)
  fp = c(file.path(save_path_period, paste("t_", substr(libs[i], 7, (nchar(libs[i])-6)), ".png", sep="")))
  ggsave(file = fp, plot = p)
  
  if (i > 3) {
    metac <- meta[sapply(meta, function(x){sum(x$TC.control) > 0})]
    metac <- Reduce("+", lapply(metac, function(x){x$TC.control/sum(x$TC.control)}))
    amplitudes <- abs(fft(metac))
    amp <- amplitudes[2:(length(amplitudes)/2+1)]
    periods <- 1/periodogram(metac, method = "clone")$freq
    df <- data.frame(periods, amp)
    t <- paste("FFT, ", substr(libs[i], 9, (nchar(libs[i])-6)), ", TC.control", sep="")
    p <- ggplot(df, aes(x=periods, y=amp)) + geom_line() + ggtitle(t) + xlim(0, 8)
    fp = c(file.path(save_path_period, paste("c_", substr(libs[i], 7, (nchar(libs[i])-6)), ".png", sep="")))
    ggsave(file = fp, plot = p)
    
    metal <- meta[sapply(meta, function(x){sum(x$log2ratio) > 0})]
    metal <- Reduce("+", lapply(metal, function(x){x$log2ratio/sum(x$log2ratio)}))
    amplitudes <- abs(fft(metal))
    amp <- amplitudes[2:(length(amplitudes)/2+1)]
    periods <- 1/periodogram(metal, method = "clone")$freq
    df <- data.frame(periods, amp)
    t <- paste("FFT, ", substr(libs[i], 9, (nchar(libs[i])-6)), ", log2ratio", sep="")
    p <- ggplot(df, aes(x=periods, y=amp)) + geom_line() + ggtitle(t) + xlim(0, 8)
    fp = c(file.path(save_path_period, paste("l_", substr(libs[i], 7, (nchar(libs[i])-6)), ".png", sep="")))
    ggsave(file = fp, plot = p)
  }
  
}



###################################################### START ################################################################

### tx with utr5 longer than 30
txlen_start <- txlen[txlen$utr5_len > 100,]
rownames(txlen_start) <- txlen_start$tx_name
txlen_start <- txlen_start[order(rownames(txlen_start)),]

start <- GRanges(Rle(rownames(txlen_start), rep(1, nrow(txlen_start))), IRanges(txlen_start$utr5_len-99, width=rep(200, nrow(txlen_start))))
start <- split(start, seqnames(start))

scale <- c(-100:100)[-101]

###

save_path_start <- "/Volumes/USELESS/META/SHAPES/SE/start"

for (i in 1:length(libs)) {
  u <- get(load(file = c(file.path(libs_path, libs[i]))))
  u <- unlist(u)
  names(u) <- substring(names(u), 1, 18)
  gr <- GRanges(seqnames=Rle(names(u)),
                IRanges(start(u), width=1),
                strand='*')
  gr$TC.treated <- u$TC.treated
  if (i > 3) {
    gr$TC.control <- u$TC.control
    gr$log2ratio <- u$log2ratio
  }
  
  meta <- subsetByOverlaps(gr, start, ignore.strand=TRUE)
  meta <- split(meta, seqnames(meta))
  meta <- meta[sapply(meta, function(x){length(x) > 0})]
  
  metat <- meta[sapply(meta, function(x){sum(x$TC.treated) > 0})]
  metat <- Reduce("+", lapply(metat, function(x){x$TC.treated/sum(x$TC.treated)}))
  d <- data.frame(scale = scale, shape = metat)
  t <- paste("START, ", substr(libs[i], 7, (nchar(libs[i])-6)), ", TC.treated", sep="")
  p_start <- ggplot(d, aes(x = scale, y = shape)) + geom_bar(stat = "identity") + ggtitle(t)
  fstart = c(file.path(save_path_start, paste("t_", substr(libs[i], 7, (nchar(libs[i])-6)), ".png", sep="")))
  ggsave(file = fstart, p_start)
  
  if (i > 3) {
    metac <- meta[sapply(meta, function(x){sum(x$TC.control) > 0})]
    metac <- Reduce("+", lapply(metac, function(x){x$TC.control/sum(x$TC.control)}))
    d <- data.frame(scale = scale, shape = metac)
    t <- paste("START, ", substr(libs[i], 7, (nchar(libs[i])-6)), ", TC.control", sep="")
    p_start <- ggplot(d, aes(x = scale, y = shape)) + geom_bar(stat = "identity") + ggtitle(t)
    fstart = c(file.path(save_path_start, paste("c_", substr(libs[i], 7, (nchar(libs[i])-6)), ".png", sep="")))
    ggsave(file = fstart, p_start)
    
    metal <- meta[sapply(meta, function(x){sum(x$log2ratio) > 0})]
    metal <- Reduce("+", lapply(metal, function(x){x$log2ratio/sum(x$log2ratio)}))
    d <- data.frame(scale = scale, shape = metal)
    t <- paste("START, ", substr(libs[i], 7, (nchar(libs[i])-6)), ", log2ratio", sep="")
    p_start <- ggplot(d, aes(x = scale, y = shape)) + geom_bar(stat = "identity") + ggtitle(t)
    fstart = c(file.path(save_path_start, paste("l_", substr(libs[i], 7, (nchar(libs[i])-6)), ".png", sep="")))
    ggsave(file = fstart, p_start)
  }
}


###################################################### STOP #################################################################
txlen_stop <- txlen[txlen$utr3_len > 100,]
stop <- GRanges(Rle(rownames(txlen_stop), rep(1, nrow(txlen_stop))), IRanges(txlen_stop$utr5_len+txlen_stop$cds_len-99, width=rep(200, nrow(txlen_stop))))
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

save_path_stop <- "/Volumes/USELESS/META/SHAPES/SE/stop"


for (i in 1:length(libs)) {
  u <- get(load(file = c(file.path(libs_path, libs[i]))))
  u <- unlist(u)
  names(u) <- substring(names(u), 1, 18)
  gr <- GRanges(seqnames=Rle(names(u)),
                IRanges(start(u), width=1),
                strand='*')
  gr$TC.treated <- u$TC.treated
  if (i > 3) {
    gr$TC.control <- u$TC.control
    gr$log2ratio <- u$log2ratio
  }
  
  meta <- subsetByOverlaps(gr, stop, ignore.strand=TRUE)
  meta <- split(meta, seqnames(meta))
  meta <- meta[sapply(meta, function(x){length(x) > 0})]
  
  metad <- meta[sapply(meta, function(x){sum(x$TC.treated) > 0})]
  shape <- metad
  df_shape <- data.frame(scale = scale,
                         shape = Reduce("+", lapply(shape, function(x){x$TC.treated/sum(x$TC.treated)})))
  t <- paste("STOP, ", substr(libs[i], 7, (nchar(libs[i])-6)), ", TC.treated", sep="")
  p <- ggplot(df_shape, aes(x=scale, y=shape)) + geom_bar(stat = "identity") + ggtitle(t)
  fstop = c(file.path(save_path_stop, paste("t_", substr(libs[i], 7, (nchar(libs[i])-6)), ".png", sep="")))
  ggsave(file = fstop, p)
  
  shape_taa <- metad[names(metad) %in% taa]
  shape_tga <- metad[names(metad) %in% tga]
  shape_tag <- metad[names(metad) %in% tag]
  df_shape <- data.frame(scale = scale,
                         shape_taa = Reduce("+", lapply(shape_taa, function(x){x$TC.treated/sum(x$TC.treated)})),
                         shape_tga = Reduce("+", lapply(shape_tga, function(x){x$TC.treated/sum(x$TC.treated)})),
                         shape_tag = Reduce("+", lapply(shape_tag, function(x){x$TC.treated/sum(x$TC.treated)})))
  t <- paste("STOP, ", substr(libs[i], 7, (nchar(libs[i])-6)), ", TC.treated, TAA", sep="")
  p_taa <- ggplot(df_shape, aes(x=scale, y=shape_taa)) + geom_bar(stat = "identity") + ggtitle(t)
  fstop_taa = c(file.path(save_path_stop, paste("t_", substr(libs[i], 7, (nchar(libs[i])-6)), "_TAA.png", sep="")))
  ggsave(file = fstop_taa, p_taa)
  t <- paste("STOP, ", substr(libs[i], 7, (nchar(libs[i])-6)), ", TC.treated, TGA", sep="")
  p_tga <- ggplot(df_shape, aes(x=scale, y=shape_tga)) + geom_bar(stat = "identity") + ggtitle(t)
  fstop_tga = c(file.path(save_path_stop, paste("t_", substr(libs[i], 7, (nchar(libs[i])-6)), "_TGA.png", sep="")))
  ggsave(file = fstop_tga, p_tga)
  t <- paste("STOP, ", substr(libs[i], 7, (nchar(libs[i])-6)), ", TC.treated, TAG", sep="")
  p_tag <- ggplot(df_shape, aes(x=scale, y=shape_tag)) + geom_bar(stat = "identity") + ggtitle(t)
  fstop_tag = c(file.path(save_path_stop, paste("t_", substr(libs[i], 7, (nchar(libs[i])-6)), "_TAG.png", sep="")))
  ggsave(file = fstop_tag, p_tag)
  
  if (i > 3) {
    metac <- meta[sapply(meta, function(x){sum(x$TC.control) > 0})]
    shape <- metac
    df_shape <- data.frame(scale = scale,
                           shape = Reduce("+", lapply(shape, function(x){x$TC.control/sum(x$TC.control)})))
    t <- paste("STOP, ", substr(libs[i], 7, (nchar(libs[i])-6)), ", TC.control", sep="")
    p <- ggplot(df_shape, aes(x=scale, y=shape)) + geom_bar(stat = "identity") + ggtitle(t)
    fstop = c(file.path(save_path_stop, paste("c_", substr(libs[i], 7, (nchar(libs[i])-6)), ".png", sep="")))
    ggsave(file = fstop, p)
    
    shape_taa <- metac[names(metac) %in% taa]
    shape_tga <- metac[names(metac) %in% tga]
    shape_tag <- metac[names(metac) %in% tag]
    df_shape <- data.frame(scale = scale,
                           shape_taa = Reduce("+", lapply(shape_taa, function(x){x$TC.control/sum(x$TC.control)})),
                           shape_tga = Reduce("+", lapply(shape_tga, function(x){x$TC.control/sum(x$TC.control)})),
                           shape_tag = Reduce("+", lapply(shape_tag, function(x){x$TC.control/sum(x$TC.control)})))
    t <- paste("STOP, ", substr(libs[i], 7, (nchar(libs[i])-6)), ", TC.control, TAA", sep="")
    p_taa <- ggplot(df_shape, aes(x=scale, y=shape_taa)) + geom_bar(stat = "identity") + ggtitle(t)
    fstop_taa = c(file.path(save_path_stop, paste("c_", substr(libs[i], 7, (nchar(libs[i])-6)), "_TAA.png", sep="")))
    ggsave(file = fstop_taa, p_taa)
    t <- paste("STOP, ", substr(libs[i], 7, (nchar(libs[i])-6)), ", TC.control, TGA", sep="")
    p_tga <- ggplot(df_shape, aes(x=scale, y=shape_tga)) + geom_bar(stat = "identity") + ggtitle(t)
    fstop_tga = c(file.path(save_path_stop, paste("c_", substr(libs[i], 7, (nchar(libs[i])-6)), "_TGA.png", sep="")))
    ggsave(file = fstop_tga, p_tga)
    t <- paste("STOP, ", substr(libs[i], 7, (nchar(libs[i])-6)), ", TC.control, TAG", sep="")
    p_tag <- ggplot(df_shape, aes(x=scale, y=shape_tag)) + geom_bar(stat = "identity") + ggtitle(t)
    fstop_tag = c(file.path(save_path_stop, paste("c_", substr(libs[i], 7, (nchar(libs[i])-6)), "_TAG.png", sep="")))
    ggsave(file = fstop_tag, p_tag)
    
    
    metal <- meta[sapply(meta, function(x){sum(x$log2ratio) > 0})]
    shape <- metal
    df_shape <- data.frame(scale = scale,
                           shape = Reduce("+", lapply(shape, function(x){x$log2ratio/sum(x$log2ratio)})))
    t <- paste("STOP, ", substr(libs[i], 7, (nchar(libs[i])-6)), ", log2ratio", sep="")
    p <- ggplot(df_shape, aes(x=scale, y=shape)) + geom_bar(stat = "identity") + ggtitle(t)
    fstop = c(file.path(save_path_stop, paste("l_", substr(libs[i], 7, (nchar(libs[i])-6)), ".png", sep="")))
    ggsave(file = fstop, p)
    
    shape_taa <- metal[names(metal) %in% taa]
    shape_tga <- metal[names(metal) %in% tga]
    shape_tag <- metal[names(metal) %in% tag]
    df_shape <- data.frame(scale = scale,
                           shape_taa = Reduce("+", lapply(shape_taa, function(x){x$log2ratio/sum(x$log2ratio)})),
                           shape_tga = Reduce("+", lapply(shape_tga, function(x){x$log2ratio/sum(x$log2ratio)})),
                           shape_tag = Reduce("+", lapply(shape_tag, function(x){x$log2ratio/sum(x$log2ratio)})))
    t <- paste("STOP, ", substr(libs[i], 7, (nchar(libs[i])-6)), ", log2ratio, TAA", sep="")
    p_taa <- ggplot(df_shape, aes(x=scale, y=shape_taa)) + geom_bar(stat = "identity") + ggtitle(t)
    fstop_taa = c(file.path(save_path_stop, paste("l_", substr(libs[i], 7, (nchar(libs[i])-6)), "_TAA.png", sep="")))
    ggsave(file = fstop_taa, p_taa)
    t <- paste("STOP, ", substr(libs[i], 7, (nchar(libs[i])-6)), ", log2ratio, TGA", sep="")
    p_tga <- ggplot(df_shape, aes(x=scale, y=shape_tga)) + geom_bar(stat = "identity") + ggtitle(t)
    fstop_tga = c(file.path(save_path_stop, paste("l_", substr(libs[i], 7, (nchar(libs[i])-6)), "_TGA.png", sep="")))
    ggsave(file = fstop_tga, p_tga)
    t <- paste("STOP, ", substr(libs[i], 7, (nchar(libs[i])-6)), ", log2ratio, TAG", sep="")
    p_tag <- ggplot(df_shape, aes(x=scale, y=shape_tag)) + geom_bar(stat = "identity") + ggtitle(t)
    fstop_tag = c(file.path(save_path_stop, paste("l_", substr(libs[i], 7, (nchar(libs[i])-6)), "_TAG.png", sep="")))
    ggsave(file = fstop_tag, p_tag)
  }
}




#############################
######### tx_profiles

### tx with utr5 and cds longer than 100
txlen_txp <- txlen

lct <- txlen_txp[(txlen_txp$utr5_len > 50) & (txlen_txp$cds_len > 100) & (txlen_txp$utr3_len > 50), ]

lct <- lct[lct$tx_name %in% matrix_of_everything[matrix_of_everything$rna_cell1K > 1,]$n,]
lct <- lct[lct$tx_name %in% matrix_of_everything[matrix_of_everything$shape_cell256_invivo_log2ratio > 0.5,]$n,]

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

save_path_tx_unit <- "/Volumes/USELESS/META/SHAPES/SE/tx_profiles/tx_unit_length"

for (i in 1:length(libs)) {
  u <- get(load(file = c(file.path(libs_path, libs[i]))))
  u <- unlist(u)
  names(u) <- substring(names(u), 1, 18)
  gr <- GRanges(seqnames=Rle(names(u)),
                IRanges(start(u), width=1),
                strand='*')
  gr$TC.treated <- u$TC.treated
  if (i > 3) {
    gr$TC.control <- u$TC.control
    gr$log2ratio <- u$log2ratio
  }
  gr <- split(gr, seqnames(gr))
  
  sev <- normalizeTxUL(gr, sam, norm="TC.treated")
  t <- paste("tx_profile, ", substr(libs[i], 7, (nchar(libs[i])-6)), ", TC.treated", sep="")
  p_tx <- ggplot() + geom_step(data = sev, aes(start, value), size=0.1) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) +
    xlab("tx_unit_length") + ylab("sum(SHAPE)") + ggtitle(t)
  fp = c(file.path(save_path_tx_unit, paste("tx_", substr(libs[i], 7, (nchar(libs[i])-6)), "_TREATED.png", sep="")))
  ggsave(file = fp, p_tx)
  
  if (i > 3) {
    sev <- normalizeTxUL(gr, sam, norm="TC.control")
    t <- paste("tx_profile, ", substr(libs[i], 7, (nchar(libs[i])-6)), ", TC.control", sep="")
    p_tx <- ggplot() + geom_step(data = sev, aes(start, value), size=0.1) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) +
      xlab("tx_unit_length") + ylab("sum(SHAPE)") + ggtitle(t)
    fp = c(file.path(save_path_tx_unit, paste("tx_", substr(libs[i], 7, (nchar(libs[i])-6)), "_CONTROL.png", sep="")))
    ggsave(file = fp, p_tx)
    
    sev <- normalizeTxUL(gr, sam, norm="log2ratio")
    t <- paste("tx_profile, ", substr(libs[i], 7, (nchar(libs[i])-6)), ", log2ratio", sep="")
    p_tx <- ggplot() + geom_step(data = sev, aes(start, value), size=0.1) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) +
      xlab("tx_unit_length") + ylab("sum(SHAPE)") + ggtitle(t)
    fp = c(file.path(save_path_tx_unit, paste("tx_", substr(libs[i], 7, (nchar(libs[i])-6)), "_LOG2RATIO.png", sep="")))
    ggsave(file = fp, p_tx)
  }
  
}


#######################################
################ tx_profiles: utr5, cds, utr3 #######################

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


save_path_tx_ucu <- "/Volumes/USELESS/META/SHAPES/SE/tx_profiles/utr5_cds_utr3"

for (i in 1:length(libs)) {
  u <- get(load(file = c(file.path(libs_path, libs[i]))))
  u <- unlist(u)
  names(u) <- substring(names(u), 1, 18)
  gr <- GRanges(seqnames=Rle(names(u)),
                IRanges(start(u), width=1),
                strand='*')
  gr$TC.treated <- u$TC.treated
  if (i > 3) {
    gr$TC.control <- u$TC.control
    gr$log2ratio <- u$log2ratio
  }
  gr <- split(gr, seqnames(gr))
  
  sev <- normalizeTx(gr, txlen, sam, norm="TC.treated")
  t <- paste("tx_profile, ", substr(libs[i], 7, (nchar(libs[i])-6)), ", TC.treated", sep="")
  p_tx <- ggplot() + geom_step(data = sev, aes(start, value), size=0.1) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) +
    xlab("leader, CDS, trailer") + ylab("sum(SHAPE)") + ggtitle(t)
  fp = c(file.path(save_path_tx_ucu, paste("tx_", substr(libs[i], 7, (nchar(libs[i])-6)), "_TREATED.png", sep="")))
  ggsave(file = fp, p_tx)
  
  if (i > 3) {
    sev <- normalizeTx(gr, txlen, sam, norm="TC.control")
    t <- paste("tx_profile, ", substr(libs[i], 7, (nchar(libs[i])-6)), ", TC.control", sep="")
    p_tx <- ggplot() + geom_step(data = sev, aes(start, value), size=0.1) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) +
      xlab("leader, CDS, trailer") + ylab("sum(SHAPE)") + ggtitle(t)
    fp = c(file.path(save_path_tx_ucu, paste("tx_", substr(libs[i], 7, (nchar(libs[i])-6)), "_CONTROL.png", sep="")))
    ggsave(file = fp, p_tx)
    
    sev <- normalizeTx(gr, txlen, sam, norm="log2ratio")
    t <- paste("tx_profile, ", substr(libs[i], 7, (nchar(libs[i])-6)), ", log2ratio", sep="")
    p_tx <- ggplot() + geom_step(data = sev, aes(start, value), size=0.1) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) +
      xlab("leader, CDS, trailer") + ylab("sum(SHAPE)") + ggtitle(t)
    fp = c(file.path(save_path_tx_ucu, paste("tx_", substr(libs[i], 7, (nchar(libs[i])-6)), "_LOG2RATIO.png", sep="")))
    ggsave(file = fp, p_tx)
  }
  
}


######################################
######### EXONS (splice sites)

exons <- exonsBy(txdb_can, by="tx", use.names=TRUE)
exons <- exons[names(exons) %in% txlen$tx_name]
we <- width(exons)
exons <- cumsum(width(exons))
exons <- exons[we > 60]
exons <- exons[sapply(exons, function(x){length(x) > 3})]
exons <- sapply(exons, function(x){x[2:(length(x)-2)]})
exons <- unlist(exons)
exons <- exons[exons > 30]

ex30 <- GRanges(Rle(names(exons), rep(1, length(exons))), IRanges(exons-29, width=rep(60, length(exons))))
ex30 <- split(ex30, seqnames(ex30))
names(ex30) <- substring(names(ex30), 1, 18)
ex30 <- GRanges(seqnames=Rle(names(ex30)),
                IRanges(as.numeric(start(ex30)), width=60),
                strand='*')

scale <- c(-30:30)[-31]

save_path_exons <- "/Volumes/USELESS/META/SHAPES/SE/exons"

for (i in 1:length(libs)) {
  u <- get(load(file = c(file.path(libs_path, libs[i]))))
  u <- unlist(u)
  names(u) <- substring(names(u), 1, 18)
  gr <- GRanges(seqnames=Rle(names(u)),
                IRanges(start(u), width=1),
                strand='*')
  gr$TC.treated <- u$TC.treated
  if (i > 3) {
    gr$TC.control <- u$TC.control
    gr$log2ratio <- u$log2ratio
  }
  
  meta <- subsetByOverlaps(gr, ex30, ignore.strand=TRUE)
  # split every 60
  meta <- split(meta, ceiling(seq_along(meta)/60))
  
  metat <- meta[sapply(meta, function(x){sum(x$TC.treated) > 0})]
  metat <- Reduce("+", lapply(metat, function(x){x$TC.treated/sum(x$TC.treated)}))
  d <- data.frame(scale = scale, shape = metat)
  t <- paste("EXONS, ", substr(libs[i], 7, (nchar(libs[i])-6)), ", TC.treated", sep="")
  p_start <- ggplot(d, aes(x = scale, y = shape)) + geom_bar(stat = "identity") + ggtitle(t)
  fstart = c(file.path(save_path_exons, paste("t_", substr(libs[i], 7, (nchar(libs[i])-6)), ".png", sep="")))
  ggsave(file = fstart, p_start)
  
  if (i > 3) {
    metac <- meta[sapply(meta, function(x){sum(x$TC.control) > 0})]
    metac <- Reduce("+", lapply(metac, function(x){x$TC.control/sum(x$TC.control)}))
    d <- data.frame(scale = scale, shape = metac)
    t <- paste("EXONS, ", substr(libs[i], 7, (nchar(libs[i])-6)), ", TC.control", sep="")
    p_start <- ggplot(d, aes(x = scale, y = shape)) + geom_bar(stat = "identity") + ggtitle(t)
    fstart = c(file.path(save_path_exons, paste("c_", substr(libs[i], 7, (nchar(libs[i])-6)), ".png", sep="")))
    ggsave(file = fstart, p_start)
    
    metal <- meta[sapply(meta, function(x){sum(x$log2ratio) > 0})]
    metal <- Reduce("+", lapply(metal, function(x){x$log2ratio/sum(x$log2ratio)}))
    d <- data.frame(scale = scale, shape = metal)
    t <- paste("EXONS, ", substr(libs[i], 7, (nchar(libs[i])-6)), ", log2ratio", sep="")
    p_start <- ggplot(d, aes(x = scale, y = shape)) + geom_bar(stat = "identity") + ggtitle(t)
    fstart = c(file.path(save_path_exons, paste("l_", substr(libs[i], 7, (nchar(libs[i])-6)), ".png", sep="")))
    ggsave(file = fstart, p_start)
  }
}



