### periodicity in CDSs
library(GenomicFeatures)
library(ggplot2)
library(GeneCycle)
library(seqinr)
library(data.table)
txdb_can <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/danio_rerio/Danio_rerio.GRCz10.84.gtf", format = "gtf")
txLengths <- transcriptLengths(txdb_can, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
rownames(txLengths) <- txLengths$tx_name
txLengths <- txLengths[order(rownames(txLengths)),]

cds150 <- GRanges(Rle(rownames(txLengths), rep(1, nrow(txLengths))), IRanges(txLengths$utr5_len+1, width=rep(150, nrow(txLengths))))
cds150 <- split(cds150, seqnames(cds150))

# shape
# shape_cell24 <- get(load(file = "/Volumes/USELESS/test/normalized/2-4cell.Rsave"))
# shape_cell256 <- get(load(file = "/Volumes/USELESS/test/normalized/256cell.Rsave"))
# shape_oblong <- get(load(file = "/Volumes/USELESS/test/normalized/oblong_invivo.Rsave"))
# shape_oblong_CHX <- get(load(file = "/Volumes/USELESS/test/normalized/oblong_CHX_invivo.Rsave"))
# rm(shapeseq_norm)
# # shapes
# shapes_cell24 <- get(load(file = "/Volumes/USELESS/test/normalized/cell2_4_invivo_NAI_N3.Rsave"))
# shapes_cell1_invitro <- get(load(file = "/Volumes/USELESS/test/normalized/cell1_invitro_NAI_N3.Rsave"))
# shapes_cell1_invitro_nonsel <- get(load(file = "/Volumes/USELESS/test/normalized/cell1_invitro_NAI_N3_nonsel.Rsave"))
# rm(shapes_norm)

####### get 150nt fragments on CDSs
load(file = "/Volumes/USELESS/test/normalized/256cell.Rsave")
u <- unlist(shapeseq_norm)
u$log2ratio[u$log2ratio < 0] <- 0
#save(u, file = "/Volumes/USELESS/test/normalized/u_256cell.Rsave")

f150 <- subsetByOverlaps(u, cds150, ignore.strand=TRUE)
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
p4 <- ggplot(df, aes(x=periods, y=amp)) + geom_line() + ggtitle("FFT, shape_cell256, TC.control") + xlim(0, 8)
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/periodicity/p4.png", plot = p4)

amplitudes <- abs(fft(metat))
amp <- amplitudes[2:(length(amplitudes)/2+1)]
periods <- 1/periodogram(metat, method = "clone")$freq
df <- data.frame(periods, amp)
p5 <- ggplot(df, aes(x=periods, y=amp)) + geom_line() + ggtitle("FFT, shape_cell256, TC.treated") + xlim(0, 8)
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/periodicity/p5.png", plot = p5)

amplitudes <- abs(fft(metal))
amp <- amplitudes[2:(length(amplitudes)/2+1)]
periods <- 1/periodogram(metal, method = "clone")$freq
df <- data.frame(periods, amp)
p6 <- ggplot(df, aes(x=periods, y=amp)) + geom_line() + ggtitle("FFT, shape_cell256, log2ratio") + xlim(0, 8)
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/periodicity/p6.png", plot = p6)

########################################

###############
load(file = "/Volumes/USELESS/test/normalized/cell2_4_invivo_NAI_N3.Rsave")
u <- unlist(shapes_norm)
save(u, file = "/Volumes/USELESS/test/normalized/u_cell2_4_invivo_NAI_N3.Rsave")

f150 <- subsetByOverlaps(u, cds150, ignore.strand=TRUE)
f150 <- split(f150, seqnames(f150))
meta <- f150[sapply(f150, function(x){length(x) > 0})]

metat <- meta[sapply(meta, function(x){sum(x$TC) > 0})]
metaw <- meta[sapply(meta, function(x){sum(x$winsor) > 0})]

metat <- Reduce("+", lapply(metat, function(x){x$TC/sum(x$TC)}))
metaw <- Reduce("+", lapply(metaw, function(x){x$winsor/sum(x$winsor)}))

amplitudes <- abs(fft(metat))
amp <- amplitudes[2:(length(amplitudes)/2+1)]
periods <- 1/periodogram(metat, method = "clone")$freq
df <- data.frame(periods, amp)
p13 <- ggplot(df, aes(x=periods, y=amp)) + geom_line() + ggtitle("FFT, shapes_cell24, TC") + xlim(0, 8)
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/periodicity/p13.png", plot = p13)

amplitudes <- abs(fft(metaw))
amp <- amplitudes[2:(length(amplitudes)/2+1)]
periods <- 1/periodogram(metaw, method = "clone")$freq
df <- data.frame(periods, amp)
p14 <- ggplot(df, aes(x=periods, y=amp)) + geom_line() + ggtitle("FFT, shapes_cell24, winsor") + xlim(0, 8)
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/periodicity/p14.png", plot = p14)

########################################

# convert p1.png p2.png p3.png p4.png p5.png p6.png +append top.png
# convert p7.png p8.png p9.png p10.png p11.png p12.png +append mid.png
# convert p13.png p14.png p15.png p16.png p17.png p18.png +append bot.png
# convert top.png mid.png bot.png -append periodicity.png


#####

shape_libs <- c("shape_cell1_et1", "shape_cell1_et2", "shape_cell24", "shape_cell256", "shape_cell1K_et1", "shape_cell1K_et2", "shape_oblong", "shape_oblong_CHX")
shapes_libs <- c("shapes_cell1_invitro", "shapes_cell1_invitro_nonsel", "shapes_cell24")


for (i in 1:8) {
  
  u <- unlist(eval(parse(text = shape_libs[i])))
  f150 <- subsetByOverlaps(u, cds150, ignore.strand=TRUE)
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
  d <- data.frame(periods, amp)
  ggplot(d, aes(x=periods, y=amp)) + geom_line() + ggtitle(paste(shape_libs[i], "control")) + xlim(0, 8)
  fp <- paste0("/Volumes/USELESS/META/SHAPES_NEW/periodicity/p", i, "c.png", collapse = NULL)
  ggsave(fp)
  
  amplitudes <- abs(fft(metat))
  amp <- amplitudes[2:(length(amplitudes)/2+1)]
  periods <- 1/periodogram(metat, method = "clone")$freq
  d <- data.frame(periods, amp)
  ggplot(d, aes(x=periods, y=amp)) + geom_line() + ggtitle(paste(shape_libs[i], "treated")) + xlim(0, 8)
  fp <- paste0("/Volumes/USELESS/META/SHAPES_NEW/periodicity/p", i, "t.png", collapse = NULL)
  ggsave(fp)
  
  amplitudes <- abs(fft(metal))
  amp <- amplitudes[2:(length(amplitudes)/2+1)]
  periods <- 1/periodogram(metal, method = "clone")$freq
  d <- data.frame(periods, amp)
  ggplot(d, aes(x=periods, y=amp)) + geom_line() + ggtitle(paste(shape_libs[i], "log2ratio")) + xlim(0, 8)
  fp <- paste0("/Volumes/USELESS/META/SHAPES_NEW/periodicity/p", i, "l.png", collapse = NULL)
  ggsave(fp)
  
}


for (i in 1:3) {
  
  u <- unlist(eval(parse(text = shapes_libs[i])))
  f150 <- subsetByOverlaps(u, cds150, ignore.strand=TRUE)
  f150 <- split(f150, seqnames(f150))
  meta <- f150[sapply(f150, function(x){length(x) > 0})]
  
  metat <- meta[sapply(meta, function(x){sum(x$TC) > 0})]
  
  metat <- Reduce("+", lapply(metat, function(x){x$TC/sum(x$TC)}))
  
  amplitudes <- abs(fft(metat))
  amp <- amplitudes[2:(length(amplitudes)/2+1)]
  periods <- 1/periodogram(metat, method = "clone")$freq
  d <- data.frame(periods, amp)
  ggplot(d, aes(x=periods, y=amp)) + geom_line() + ggtitle(paste(shapes_libs[i], "TC")) + xlim(0, 8)
  fp <- paste0("/Volumes/USELESS/META/SHAPES_NEW/periodicity/s", i, ".png", collapse = NULL)
  ggsave(fp)
  
}

## convert *l.png +append log2rat.png
## convert *c.png +append control.png
## convert *t.png +append treated.png
## convert control.png treated.png log2rat.png -append periodicity_shape.png

## convert s* +append periodicity_shapes.png

##### periodicity on GC content
fasta_cds <- read.fasta("/Volumes/USELESS/DATA/fasta/danio_rerio/GRCz10/cds/Danio_rerio.GRCz10.cds.all.fa")
names(fasta_cds) <- sapply(names(fasta_cds), function(x){substr(x,1,18)})
fasta_cds <- fasta_cds[unlist(lapply(fasta_cds, function(x){length(x) %% 3 == 0}))]

### GC content
sam <- sample(names(fasta_cds), 1000)
samGC <- fasta_cds[names(fasta_cds) %in% sam]
gc <- data.table(a = 2, t = 2, g = 3, c = 3, n = 2)

gc_score <- sapply(samGC, function(s){as.vector(sapply(s, function(x){gc[[x]]}))})
gc_score <- gc_score[sapply(gc_score, function(x){length(x) > 150})]
meta <- Reduce("+", lapply(gc_score, function(x){x[1:150]}))

amplitudes <- abs(fft(meta))
amp <- amplitudes[2:(length(amplitudes)/2+1)]
periods <- 1/periodogram(meta, method = "clone")$freq
df <- data.frame(periods, amp)
ggplot(df, aes(x=periods, y=amp)) + geom_line() + ggtitle("FFT, GC content (G,C=3; A,T=2)") + xlim(0, 8)
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/periodicity/GC.png")


########### redo on 2-4 and 256cell
load(file = "/Volumes/USELESS/DATA/Shape-Seq/dTCR/PAIR_GC_cell24.Rsave")
GR <- split(GR, names(GR))

txlen <- txLengths[names(GR),]
txlen <- txlen[txlen$cds_len > 150,]
GR <- GR[names(GR) %in% rownames(txlen)]
utr5len <- txlen$utr5_len

dmso <- mapply(function(x,y){x$TC.control[(y+1):(y+151)]}, GR, utr5len)
nai <- mapply(function(x,y){x$TC.treated[(y+1):(y+151)]}, GR, utr5len)
dtcr <- mapply(function(x,y){x$dTCR[(y+1):(y+151)]}, GR, utr5len)

meta_dmso <- rowSums(dmso)
meta_nai <- rowSums(nai)
meta_dtcr <- rowSums(dtcr)


amplitudes <- abs(fft(meta_dmso))
amp <- amplitudes[2:(length(amplitudes)/2+1)]
periods <- 1/periodogram(meta_dmso, method = "clone")$freq
df <- data.frame(periods, amp)
ggplot(df, aes(x=periods, y=amp)) + geom_line() + ggtitle("FFT, shape_cell24, TC.control") + xlim(0, 8)
#ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/periodicity/p4.png", plot = p4)

amplitudes <- abs(fft(meta_nai))
amp <- amplitudes[2:(length(amplitudes)/2+1)]
periods <- 1/periodogram(meta_nai, method = "clone")$freq
df <- data.frame(periods, amp)
ggplot(df, aes(x=periods, y=amp)) + geom_line() + ggtitle("FFT, shape_cell24, TC.treated") + xlim(0, 8)
#ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/periodicity/p5.png", plot = p5)

amplitudes <- abs(fft(meta_dtcr))
amp <- amplitudes[2:(length(amplitudes)/2+1)]
periods <- 1/periodogram(meta_dtcr, method = "clone")$freq
df <- data.frame(periods, amp)
ggplot(df, aes(x=periods, y=amp)) + geom_line() + ggtitle("FFT, shape_cell24, dTRC") + xlim(0, 8)
#ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/periodicity/p5.png", plot = p5)



