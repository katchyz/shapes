### analyze SHAPES libraries (June 2018)
library(GenomicFeatures)
library(ggplot2)
library(plyr)
library(data.table)
library(reshape2)
library(GeneCycle)
library(seqinr)


libs_path <- "/Volumes/USELESS/DATA/SHAPES/SE/June2018"
libs <- list.files(path = libs_path)
libs <- libs[10:18]

cell256_HP <- get(load(file.path(libs_path, libs[1])))
cell256_old <- get(load(file.path(libs_path, libs[2])))
cell256_P <- get(load(file.path(libs_path, libs[3])))

oblong_invitro <- get(load(file.path(libs_path, libs[4])))
oblong <- get(load(file.path(libs_path, libs[5])))
oblong_old <- get(load(file.path(libs_path, libs[6])))
oblong_invitro_old <- get(load(file.path(libs_path, libs[7])))

pool1 <- get(load(file.path(libs_path, libs[8])))
pool2 <- get(load(file.path(libs_path, libs[9])))


#################################################### RNA / Ribo ##########################################################
load(file="/Volumes/USELESS/META/SHAPES/matrix_of_everything.Rdata")


df_cell256_HP <- data.frame(shapes_cell256_HP_TC = sapply(cell256_HP, function(x){mean(x$TC.treated)}), n = names(cell256_HP))
df_cell256_old <- data.frame(shapes_cell256_old_TC = sapply(cell256_old, function(x){mean(x$TC.treated)}), n = names(cell256_old))
df_cell256_P <- data.frame(shapes_cell256_P_TC = sapply(cell256_P, function(x){mean(x$TC.treated)}), n = names(cell256_P))
matrix_of_everything <- merge(matrix_of_everything, df_cell256_HP, by="n", all=TRUE)
matrix_of_everything <- merge(matrix_of_everything, df_cell256_old, by="n", all=TRUE)
matrix_of_everything <- merge(matrix_of_everything, df_cell256_P, by="n", all=TRUE)

df_oblong <- data.frame(shapes_oblong_TC = sapply(oblong, function(x){mean(x$TC.treated)}), n = names(oblong))
df_oblong_old <- data.frame(shapes_oblong_old_TC = sapply(oblong_old, function(x){mean(x$TC.treated)}), n = names(oblong_old))
df_oblong_invitro <- data.frame(shapes_oblong_invitro_TC = sapply(oblong_invitro, function(x){mean(x$TC.treated)}), n = names(oblong_invitro))
df_oblong_invitro_old <- data.frame(shapes_oblong_invitro_old_TC = sapply(oblong_invitro_old, function(x){mean(x$TC.treated)}), n = names(oblong_invitro_old))
matrix_of_everything <- merge(matrix_of_everything, df_oblong, by="n", all=TRUE)
matrix_of_everything <- merge(matrix_of_everything, df_oblong_old, by="n", all=TRUE)
matrix_of_everything <- merge(matrix_of_everything, df_oblong_invitro, by="n", all=TRUE)
matrix_of_everything <- merge(matrix_of_everything, df_oblong_invitro_old, by="n", all=TRUE)

df_pool1 <- data.frame(shapes_pool1_TC = sapply(pool1, function(x){mean(x$TC.treated)}), n = names(pool1))
df_pool2 <- data.frame(shapes_pool2_TC = sapply(pool2, function(x){mean(x$TC.treated)}), n = names(pool2))
matrix_of_everything <- merge(matrix_of_everything, df_pool1, by="n", all=TRUE)
matrix_of_everything <- merge(matrix_of_everything, df_pool2, by="n", all=TRUE)


#### check if they correlate with each other - if so, POOL!
ggplot(matrix_of_everything, aes(x = shapes_cell256_HP_TC, y = shapes_cell256_old_TC)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggplot(matrix_of_everything, aes(x = shapes_cell256_P_TC, y = shapes_cell256_old_TC)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggplot(matrix_of_everything, aes(x = shapes_cell256_HP_TC, y = shapes_cell256_P_TC)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggplot(matrix_of_everything, aes(x = shapes_oblong_TC, y = shapes_oblong_old_TC)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggplot(matrix_of_everything, aes(x = shapes_oblong_invitro_TC, y = shapes_oblong_invitro_old_TC)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggplot(matrix_of_everything, aes(x = shapes_pool1_TC, y = shapes_pool2_TC)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()

#### check if they correlate with SHAPE libs
ggplot(matrix_of_everything, aes(x = shapes_cell256_HP_TC, y = shape_cell256_invivo_log2ratio)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggplot(matrix_of_everything, aes(x = shapes_oblong_TC, y = shape_oblong_invivo_log2ratio)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()

## they do!!


## plot
# RNA
ggplot(matrix_of_everything, aes(x = shapes_cell256_HP_TC, y = rna_cell256)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave("/Volumes/USELESS/META/SHAPES/June2018/general/shape_vs_rna/cell256_HP.png")
ggplot(matrix_of_everything, aes(x = shapes_cell256_P_TC, y = rna_cell256)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave("/Volumes/USELESS/META/SHAPES/June2018/general/shape_vs_rna/cell256_P.png")
ggplot(matrix_of_everything, aes(x = shapes_cell256_old_TC, y = rna_cell256)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave("/Volumes/USELESS/META/SHAPES/June2018/general/shape_vs_rna/cell256_old.png")

ggplot(matrix_of_everything, aes(x = shapes_oblong_TC, y = rna_3_5h)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave("/Volumes/USELESS/META/SHAPES/June2018/general/shape_vs_rna/oblong.png")
ggplot(matrix_of_everything, aes(x = shapes_oblong_old_TC, y = rna_3_5h)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave("/Volumes/USELESS/META/SHAPES/June2018/general/shape_vs_rna/oblong_old.png")
ggplot(matrix_of_everything, aes(x = shapes_oblong_invitro_TC, y = rna_3_5h)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave("/Volumes/USELESS/META/SHAPES/June2018/general/shape_vs_rna/oblong_invitro.png")
ggplot(matrix_of_everything, aes(x = shapes_oblong_invitro_old_TC, y = rna_3_5h)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave("/Volumes/USELESS/META/SHAPES/June2018/general/shape_vs_rna/oblong_invitro_old.png")

# Ribo
ggplot(matrix_of_everything, aes(x = shapes_cell256_HP_TC, y = ribo_cell256)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave("/Volumes/USELESS/META/SHAPES/June2018/general/shape_vs_ribo/cell256_HP.png")
ggplot(matrix_of_everything, aes(x = shapes_cell256_P_TC, y = ribo_cell256)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave("/Volumes/USELESS/META/SHAPES/June2018/general/shape_vs_ribo/cell256_P.png")
ggplot(matrix_of_everything, aes(x = shapes_cell256_old_TC, y = ribo_cell256)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave("/Volumes/USELESS/META/SHAPES/June2018/general/shape_vs_ribo/cell256_old.png")

ggplot(matrix_of_everything, aes(x = shapes_oblong_TC, y = ribo_4h)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave("/Volumes/USELESS/META/SHAPES/June2018/general/shape_vs_ribo/oblong.png")
ggplot(matrix_of_everything, aes(x = shapes_oblong_old_TC, y = ribo_4h)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave("/Volumes/USELESS/META/SHAPES/June2018/general/shape_vs_ribo/oblong_old.png")
ggplot(matrix_of_everything, aes(x = shapes_oblong_invitro_TC, y = ribo_4h)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave("/Volumes/USELESS/META/SHAPES/June2018/general/shape_vs_ribo/oblong_invitro.png")
ggplot(matrix_of_everything, aes(x = shapes_oblong_invitro_old_TC, y = ribo_4h)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave("/Volumes/USELESS/META/SHAPES/June2018/general/shape_vs_ribo/oblong_invitro_old.png")

# TE
ggplot(matrix_of_everything, aes(x = shapes_cell256_HP_TC, y = ribo_cell256/rna_cell256)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave("/Volumes/USELESS/META/SHAPES/June2018/general/shape_vs_te/cell256_HP.png")
ggplot(matrix_of_everything, aes(x = shapes_cell256_P_TC, y = ribo_cell256/rna_cell256)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave("/Volumes/USELESS/META/SHAPES/June2018/general/shape_vs_te/cell256_P.png")
ggplot(matrix_of_everything, aes(x = shapes_cell256_old_TC, y = ribo_cell256/rna_cell256)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave("/Volumes/USELESS/META/SHAPES/June2018/general/shape_vs_te/cell256_old.png")

ggplot(matrix_of_everything, aes(x = shapes_oblong_TC, y = ribo_4h/rna_3_5h)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave("/Volumes/USELESS/META/SHAPES/June2018/general/shape_vs_te/oblong.png")
ggplot(matrix_of_everything, aes(x = shapes_oblong_old_TC, y = ribo_4h/rna_3_5h)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave("/Volumes/USELESS/META/SHAPES/June2018/general/shape_vs_te/oblong_old.png")
ggplot(matrix_of_everything, aes(x = shapes_oblong_invitro_TC, y = ribo_4h/rna_3_5h)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave("/Volumes/USELESS/META/SHAPES/June2018/general/shape_vs_te/oblong_invitro.png")
ggplot(matrix_of_everything, aes(x = shapes_oblong_invitro_old_TC, y = ribo_4h/rna_3_5h)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave("/Volumes/USELESS/META/SHAPES/June2018/general/shape_vs_te/oblong_invitro_old.png")



##################################################### periodicity #########################################################
txdb_can <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/danio_rerio/Danio_rerio.GRCz10.84.gtf", format = "gtf")
txLengths <- transcriptLengths(txdb_can, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
rownames(txLengths) <- txLengths$tx_name
txLengths <- txLengths[order(rownames(txLengths)),]

cds150 <- GRanges(Rle(rownames(txLengths), rep(1, nrow(txLengths))), IRanges(txLengths$utr5_len+1, width=rep(150, nrow(txLengths))))
cds150 <- split(cds150, seqnames(cds150))

save_path_period <- "/Volumes/USELESS/META/SHAPES/June2018/periodicity"

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
  t <- paste("FFT, ", substr(libs[i], 1, (nchar(libs[i])-6)), ", TC.treated", sep="")
  p2 <- ggplot(df, aes(x=periods, y=amp)) + geom_line() + ggtitle(t) + xlim(0, 8)
  fp = c(file.path(save_path_period, paste("t_", substr(libs[i], 1, (nchar(libs[i])-6)), ".png", sep="")))
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

save_path_start <- "/Volumes/USELESS/META/SHAPES/June2018/start"

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
  
  t <- paste("START, ", substr(libs[i], 1, (nchar(libs[i])-6)), ", TC.treated", sep="")
  p_start <- ggplot(d, aes(x = scale, y = shape)) + geom_bar(stat = "identity") + ggtitle(t)
  fstart = c(file.path(save_path_start, paste("start100_", substr(libs[i], 1, (nchar(libs[i])-6)), ".png", sep="")))
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

save_path_stop <- "/Volumes/USELESS/META/SHAPES/June2018/stop"


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
  
  t <- paste("STOP, ", substr(libs[i], 1, (nchar(libs[i])-6)), ", TC.treated, TAA", sep="")
  p_taa <- ggplot(df_shape, aes(x=scale, y=shape_taa)) + geom_bar(stat = "identity") + ggtitle(t)
  fstop_taa = c(file.path(save_path_stop, paste("stop100_", substr(libs[i], 1, (nchar(libs[i])-6)), "_TAA.png", sep="")))
  ggsave(file = fstop_taa, p_taa)
  
  t <- paste("STOP, ", substr(libs[i], 1, (nchar(libs[i])-6)), ", TC.treated, TGA", sep="")
  p_tga <- ggplot(df_shape, aes(x=scale, y=shape_tga)) + geom_bar(stat = "identity") + ggtitle(t)
  fstop_tga = c(file.path(save_path_stop, paste("stop100_", substr(libs[i], 1, (nchar(libs[i])-6)), "_TGA.png", sep="")))
  ggsave(file = fstop_tga, p_tga)
  
  t <- paste("STOP, ", substr(libs[i], 1, (nchar(libs[i])-6)), ", TC.treated, TAG", sep="")
  p_tag <- ggplot(df_shape, aes(x=scale, y=shape_tag)) + geom_bar(stat = "identity") + ggtitle(t)
  fstop_tag = c(file.path(save_path_stop, paste("stop100_", substr(libs[i], 1, (nchar(libs[i])-6)), "_TAG.png", sep="")))
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
  
  t <- paste("STOP, ", substr(libs[i], 1, (nchar(libs[i])-6)), ", TC.treated", sep="")
  p_taa <- ggplot(df_shape, aes(x=scale, y=shape)) + geom_bar(stat = "identity") + ggtitle(t)
  fstop_taa = c(file.path(save_path_stop, paste("stop100_", substr(libs[i], 1, (nchar(libs[i])-6)), ".png", sep="")))
  ggsave(file = fstop_taa, p_taa)
}


#########################################################################################################################
########################################### pool libs together ##########################################################
############################ cell256, oblong, oblong_invitro ############################################################

cell256 <- names(cell256_HP)[names(cell256_HP) %in% names(cell256_P)]
cell256 <- cell256[cell256 %in% names(cell256_old)] # common in all 3 (HP, P, old)
hp_p <- names(cell256_HP)[names(cell256_HP) %in% names(cell256_P)]
hp_p <- hp_p[!hp_p %in% cell256] # common in HP and P
hp_old <- names(cell256_HP)[names(cell256_HP) %in% names(cell256_old)]
hp_old <- hp_old[!hp_old %in% cell256] # common in HP and old
p_old <- names(cell256_P)[names(cell256_P) %in% names(cell256_old)]
p_old <- p_old[!p_old %in% cell256] # common in P and old
hp <- names(cell256_HP)[!names(cell256_HP) %in% names(cell256_P)]
hp <- hp[!hp %in% names(cell256_old)] # unique HP
p <- names(cell256_P)[!names(cell256_P) %in% names(cell256_HP)]
p <- p[!p %in% names(cell256_old)] # unique P
old <- names(cell256_old)[!names(cell256_old) %in% names(cell256_P)]
old <- old[!old %in% names(cell256_HP)] # unique old

## add up
shape_cell256 <- unlist(cell256_HP[names(cell256_HP) %in% cell256])
shape_cell256$TC.treated <- shape_cell256$TC.treated +
  unlist(cell256_P[names(cell256_P) %in% cell256])$TC.treated +
  unlist(cell256_old[names(cell256_old) %in% cell256])$TC.treated
shape_cell256 <- split(shape_cell256, names(shape_cell256))

shape_hp_p <- unlist(cell256_HP[names(cell256_HP) %in% hp_p])
shape_hp_p$TC.treated <- shape_hp_p$TC.treated +
  unlist(cell256_P[names(cell256_P) %in% hp_p])$TC.treated
shape_hp_p <- split(shape_hp_p, names(shape_hp_p))

shape_hp_old <- unlist(cell256_HP[names(cell256_HP) %in% hp_old])
shape_hp_old$TC.treated <- shape_hp_old$TC.treated +
  unlist(cell256_old[names(cell256_old) %in% hp_old])$TC.treated
shape_hp_old <- split(shape_hp_old, names(shape_hp_old))

shape_p_old <- unlist(cell256_P[names(cell256_P) %in% p_old])
shape_p_old$TC.treated <- shape_p_old$TC.treated +
  unlist(cell256_old[names(cell256_old) %in% p_old])$TC.treated
shape_p_old <- split(shape_p_old, names(shape_p_old))

cell256_pooled <- c(shape_cell256, shape_hp_p, shape_hp_old, shape_p_old,
                    cell256_HP[names(cell256_HP) %in% hp],
                    cell256_P[names(cell256_P) %in% p],
                    cell256_old[names(cell256_old) %in% old])
cell256_pooled <- cell256_pooled[order(names(cell256_pooled))]

# save(cell256_pooled, file="/Volumes/USELESS/DATA/SHAPES/SE/June2018/cell256_pooled.Rsave")

oblong_pooled <- unlist(oblong[names(oblong) %in% names(oblong_old)])
oblong_pooled$TC.treated <- oblong_pooled$TC.treated + unlist(oblong_old[names(oblong_old) %in% names(oblong)])$TC.treated
oblong_pooled <- split(oblong_pooled, names(oblong_pooled))
oblong_pooled <- c(oblong_pooled, oblong[!names(oblong) %in% names(oblong_old)], oblong_old[!names(oblong_old) %in% names(oblong)])
oblong_pooled <- oblong_pooled[order(names(oblong_pooled))]

# save(oblong_pooled, file="/Volumes/USELESS/DATA/SHAPES/SE/June2018/oblong_pooled.Rsave")

oblong_invitro_pooled <- unlist(oblong_invitro[names(oblong_invitro) %in% names(oblong_invitro_old)])
oblong_invitro_pooled$TC.treated <- oblong_invitro_pooled$TC.treated + unlist(oblong_invitro_old[names(oblong_invitro_old) %in% names(oblong_invitro)])$TC.treated
oblong_invitro_pooled <- split(oblong_invitro_pooled, names(oblong_invitro_pooled))
oblong_invitro_pooled <- c(oblong_invitro_pooled, oblong_invitro[!names(oblong_invitro) %in% names(oblong_invitro_old)], oblong_invitro_old[!names(oblong_invitro_old) %in% names(oblong_invitro)])
oblong_invitro_pooled <- oblong_invitro_pooled[order(names(oblong_invitro_pooled))]

# save(oblong_invitro_pooled, file="/Volumes/USELESS/DATA/SHAPES/SE/June2018/oblong_invitro_pooled.Rsave")


pooled_libs <- c("cell256_pooled.Rsave", "oblong_pooled.Rsave", "oblong_invitro_pooled.Rsave")
libs <- pooled_libs
libs <- c("cell256_pooled.Rsave", "oblong_pooled.Rsave", "oblong_invitro_pooled.Rsave", "norm2_pool1.Rsave", "norm2_pool2.Rsave")

###################################
###################################
### repeat on POOLED
###################################
###################################

load(file="/Volumes/USELESS/META/SHAPES/matrix_of_everything.Rdata")
cell256 <- get(load(file = c(file.path(libs_path, libs[1]))))
oblong <- get(load(file = c(file.path(libs_path, libs[2]))))
oblong_invitro <- get(load(file = c(file.path(libs_path, libs[3]))))
pool1 <- get(load(file.path(libs_path, "norm2_pool1.Rsave")))
pool2 <- get(load(file.path(libs_path, "norm2_pool2.Rsave")))

df_cell256 <- data.frame(shapes_cell256_TC = sapply(cell256, function(x){mean(x$TC.treated)}), n = names(cell256))
matrix_of_everything <- merge(matrix_of_everything, df_cell256, by="n", all=TRUE)

df_oblong <- data.frame(shapes_oblong_TC = sapply(oblong, function(x){mean(x$TC.treated)}), n = names(oblong))
df_oblong_invitro <- data.frame(shapes_oblong_invitro_TC = sapply(oblong_invitro, function(x){mean(x$TC.treated)}), n = names(oblong_invitro))
matrix_of_everything <- merge(matrix_of_everything, df_oblong, by="n", all=TRUE)
matrix_of_everything <- merge(matrix_of_everything, df_oblong_invitro, by="n", all=TRUE)

df_pool1 <- data.frame(shapes_pool1_TC = sapply(pool1, function(x){mean(x$TC.treated)}), n = names(pool1))
df_pool2 <- data.frame(shapes_pool2_TC = sapply(pool2, function(x){mean(x$TC.treated)}), n = names(pool2))
matrix_of_everything <- merge(matrix_of_everything, df_pool1, by="n", all=TRUE)
matrix_of_everything <- merge(matrix_of_everything, df_pool2, by="n", all=TRUE)

# save(matrix_of_everything, file="/Volumes/USELESS/META/SHAPES/matrix_of_everything.Rdata")


### STOP: for control, align at off-frame stop codons

############################################## tx profiles ##################################################
##################################### whole tx to unit length ###############################################

# get longest transcript per gene
txlen <- arrange(txLengths, gene_id, desc(tx_len))
txlen <- txlen[!duplicated(txlen$gene_id),]

rownames(txlen) <- txlen$tx_name
txlen <- txlen[order(rownames(txlen)),]

lct <- txlen[(txlen$utr5_len > 50) & (txlen$cds_len > 100) & (txlen$utr3_len > 50), ]
lct <- lct[lct$tx_name %in% matrix_of_everything[matrix_of_everything$shapes_cell256_TC > 0.3 & !is.na(matrix_of_everything$shapes_cell256_TC),]$n,]

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


save_path_tx_TC <- "/Volumes/USELESS/META/SHAPES/June2018/tx_profiles/tx_unit_length"

for (i in 1:length(libs)) {
  u <- get(load(file = c(file.path(libs_path, libs[i]))))
  u <- unlist(u)
  names(u) <- substring(names(u), 1, 18)
  gr <- GRanges(seqnames=Rle(names(u)),
                IRanges(start(u), width=1),
                strand='*')
  gr$TC.treated <- u$TC.treated
  gr <- split(gr, seqnames(gr))
  
  sev <- normalizeTxUL(gr, sam, norm="TC.treated")
  t <- paste("tx_profile, ", substr(libs[i], 1, (nchar(libs[i])-6)), ", TC.treated", sep="")
  p_tx <- ggplot() + geom_step(data = sev, aes(start, value), size=0.1) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) +
    xlab("leader, CDS, trailer") + ylab("sum(SHAPE)") + ggtitle(t)
  fp = c(file.path(save_path_tx_TC, paste("tx_", substr(libs[i], 1, (nchar(libs[i])-6)), "_TREATED.png", sep="")))
  ggsave(file = fp, p_tx)

  
}


##################################### utr5_cds_utr3 ###############################################

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


save_path_tx_TC <- "/Volumes/USELESS/META/SHAPES/June2018/tx_profiles/utr5_cds_utr3"

for (i in 1:length(libs)) {
  u <- get(load(file = c(file.path(libs_path, libs[i]))))
  u <- unlist(u)
  names(u) <- substring(names(u), 1, 18)
  gr <- GRanges(seqnames=Rle(names(u)),
                IRanges(start(u), width=1),
                strand='*')
  gr$TC.treated <- u$TC.treated
  gr <- split(gr, seqnames(gr))
  
  sev <- normalizeTx(gr, txlen, sam, norm="TC.treated")
  t <- paste("tx_profile, ", substr(libs[i], 1, (nchar(libs[i])-6)), ", TC.treated", sep="")
  p_tx <- ggplot() + geom_step(data = sev, aes(start, value), size=0.1) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) +
    xlab("leader, CDS, trailer") + ylab("sum(SHAPE)") + ggtitle(t)
  fp = c(file.path(save_path_tx_TC, paste("tx_", substr(libs[i], 1, (nchar(libs[i])-6)), "_TREATED.png", sep="")))
  ggsave(file = fp, p_tx)
  
}


