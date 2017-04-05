### periodicity in CDSs
library(GenomicFeatures)
library(ggplot2)
library(GeneCycle)
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
