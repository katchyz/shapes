## periodicity (in SHAPE): similar SHAPE-Seq/RNA-Seq; high ribo vs low ribo
# 2-4cell, 256cell, 1Kcell

# cut off on RNA-seq
#2-4cell
d <- df[df$rna_cell24 > 1 & df$shape_cell24_log2ratio > 0,]
# get high ribo and low ribo
d <- d[!is.na(d$ribo_cell24),]
low_ribo <- d[{q<-rank(d$ribo_cell24[!is.na(d$ribo_cell24)])/length(d$ribo_cell24[!is.na(d$ribo_cell24)]);q<0.1},]$n
high_ribo <- d[{q<-rank(d$ribo_cell24[!is.na(d$ribo_cell24)])/length(d$ribo_cell24[!is.na(d$ribo_cell24)]);q>0.9},]$n
# x[{q<-rank(x)/length(x);q<0.1 | q>=0.9}]
low <- shape_cell24[names(shape_cell24) %in% low_ribo]
high <- shape_cell24[names(shape_cell24) %in% high_ribo]

# low
u <- unlist(low)
u$log2ratio[u$log2ratio < 0] <- 0
f150 <- subsetByOverlaps(u, cds150, ignore.strand=TRUE)
f150 <- split(f150, seqnames(f150))
meta <- f150[sapply(f150, function(x){length(x) == 150})]

metac <- meta[sapply(meta, function(x){sum(x$TC.control) > 0})]
metac <- Reduce("+", lapply(metac, function(x){x$TC.control/sum(x$TC.control)}))
amplitudes <- abs(fft(metac))
amp <- amplitudes[2:(length(amplitudes)/2+1)]
periods <- 1/periodogram(metac, method = "clone")$freq
d <- data.frame(periods, amp)
p1 <- ggplot(d, aes(x=periods, y=amp)) + geom_line() + ggtitle("cell24, low ribo, control") + xlim(0, 8)
ggsave(file="/Volumes/USELESS/META/SHAPES_NEW/periodicity/high_low_ribo/cell24_low_ctr.png", plot=p1)

metat <- meta[sapply(meta, function(x){sum(x$TC.treated) > 0})]
metat <- Reduce("+", lapply(metat, function(x){x$TC.treated/sum(x$TC.treated)}))
amplitudes <- abs(fft(metat))
amp <- amplitudes[2:(length(amplitudes)/2+1)]
periods <- 1/periodogram(metat, method = "clone")$freq
d <- data.frame(periods, amp)
p2 <- ggplot(d, aes(x=periods, y=amp)) + geom_line() + ggtitle("cell24, low ribo, treated") + xlim(0, 8)
ggsave(file="/Volumes/USELESS/META/SHAPES_NEW/periodicity/high_low_ribo/cell24_low_tre.png", plot=p2)


# high
u <- unlist(high)
u$log2ratio[u$log2ratio < 0] <- 0
f150 <- subsetByOverlaps(u, cds150, ignore.strand=TRUE)
f150 <- split(f150, seqnames(f150))
meta <- f150[sapply(f150, function(x){length(x) == 150})]

metac <- meta[sapply(meta, function(x){sum(x$TC.control) > 0})]
metac <- Reduce("+", lapply(metac, function(x){x$TC.control/sum(x$TC.control)}))
amplitudes <- abs(fft(metac))
amp <- amplitudes[2:(length(amplitudes)/2+1)]
periods <- 1/periodogram(metac, method = "clone")$freq
d <- data.frame(periods, amp)
p3 <- ggplot(d, aes(x=periods, y=amp)) + geom_line() + ggtitle("cell24, high ribo, control") + xlim(0, 8)
ggsave(file="/Volumes/USELESS/META/SHAPES_NEW/periodicity/high_low_ribo/cell24_high_ctr.png", plot=p3)

metat <- meta[sapply(meta, function(x){sum(x$TC.treated) > 0})]
metat <- Reduce("+", lapply(metat, function(x){x$TC.treated/sum(x$TC.treated)}))
amplitudes <- abs(fft(metat))
amp <- amplitudes[2:(length(amplitudes)/2+1)]
periods <- 1/periodogram(metat, method = "clone")$freq
d <- data.frame(periods, amp)
p4 <- ggplot(d, aes(x=periods, y=amp)) + geom_line() + ggtitle("cell24, high ribo, treated") + xlim(0, 8)
ggsave(file="/Volumes/USELESS/META/SHAPES_NEW/periodicity/high_low_ribo/cell24_high_tre.png", plot=p4)


#256cell
d <- df[df$rna_cell256 > 1 & df$shape_cell256_log2ratio > 0,]
# get high ribo and low ribo
d <- d[!is.na(d$ribo_cell256),]
low_ribo <- d[{q<-rank(d$ribo_cell256[!is.na(d$ribo_cell256)])/length(d$ribo_cell256[!is.na(d$ribo_cell256)]);q<0.1},]$n
high_ribo <- d[{q<-rank(d$ribo_cell256[!is.na(d$ribo_cell256)])/length(d$ribo_cell256[!is.na(d$ribo_cell256)]);q>0.9},]$n
# x[{q<-rank(x)/length(x);q<0.1 | q>=0.9}]
low <- shape_cell256[names(shape_cell256) %in% low_ribo]
high <- shape_cell256[names(shape_cell256) %in% high_ribo]

# low
u <- unlist(low)
u$log2ratio[u$log2ratio < 0] <- 0
f150 <- subsetByOverlaps(u, cds150, ignore.strand=TRUE)
f150 <- split(f150, seqnames(f150))
meta <- f150[sapply(f150, function(x){length(x) == 150})]

metac <- meta[sapply(meta, function(x){sum(x$TC.control) > 0})]
metac <- Reduce("+", lapply(metac, function(x){x$TC.control/sum(x$TC.control)}))
amplitudes <- abs(fft(metac))
amp <- amplitudes[2:(length(amplitudes)/2+1)]
periods <- 1/periodogram(metac, method = "clone")$freq
d <- data.frame(periods, amp)
p1 <- ggplot(d, aes(x=periods, y=amp)) + geom_line() + ggtitle("cell256, low ribo, control") + xlim(0, 8)
ggsave(file="/Volumes/USELESS/META/SHAPES_NEW/periodicity/high_low_ribo/cell256_low_ctr.png", plot=p1)

metat <- meta[sapply(meta, function(x){sum(x$TC.treated) > 0})]
metat <- Reduce("+", lapply(metat, function(x){x$TC.treated/sum(x$TC.treated)}))
amplitudes <- abs(fft(metat))
amp <- amplitudes[2:(length(amplitudes)/2+1)]
periods <- 1/periodogram(metat, method = "clone")$freq
d <- data.frame(periods, amp)
p2 <- ggplot(d, aes(x=periods, y=amp)) + geom_line() + ggtitle("cell256, low ribo, treated") + xlim(0, 8)
ggsave(file="/Volumes/USELESS/META/SHAPES_NEW/periodicity/high_low_ribo/cell256_low_tre.png", plot=p2)


# high
u <- unlist(high)
u$log2ratio[u$log2ratio < 0] <- 0
f150 <- subsetByOverlaps(u, cds150, ignore.strand=TRUE)
f150 <- split(f150, seqnames(f150))
meta <- f150[sapply(f150, function(x){length(x) == 150})]

metac <- meta[sapply(meta, function(x){sum(x$TC.control) > 0})]
metac <- Reduce("+", lapply(metac, function(x){x$TC.control/sum(x$TC.control)}))
amplitudes <- abs(fft(metac))
amp <- amplitudes[2:(length(amplitudes)/2+1)]
periods <- 1/periodogram(metac, method = "clone")$freq
d <- data.frame(periods, amp)
p3 <- ggplot(d, aes(x=periods, y=amp)) + geom_line() + ggtitle("cell256, high ribo, control") + xlim(0, 8)
ggsave(file="/Volumes/USELESS/META/SHAPES_NEW/periodicity/high_low_ribo/cell256_high_ctr.png", plot=p3)

metat <- meta[sapply(meta, function(x){sum(x$TC.treated) > 0})]
metat <- Reduce("+", lapply(metat, function(x){x$TC.treated/sum(x$TC.treated)}))
amplitudes <- abs(fft(metat))
amp <- amplitudes[2:(length(amplitudes)/2+1)]
periods <- 1/periodogram(metat, method = "clone")$freq
d <- data.frame(periods, amp)
p4 <- ggplot(d, aes(x=periods, y=amp)) + geom_line() + ggtitle("cell256, high ribo, treated") + xlim(0, 8)
ggsave(file="/Volumes/USELESS/META/SHAPES_NEW/periodicity/high_low_ribo/cell256_high_tre.png", plot=p4)



#1Kcell
d <- df[df$rna_cell1K > 1 & df$shape_cell1K_et1_log2ratio > 0,]
# get high ribo and low ribo
d <- d[!is.na(d$ribo_cell1K),]
low_ribo <- d[{q<-rank(d$ribo_cell1K[!is.na(d$ribo_cell1K)])/length(d$ribo_cell1K[!is.na(d$ribo_cell1K)]);q<0.1},]$n
high_ribo <- d[{q<-rank(d$ribo_cell1K[!is.na(d$ribo_cell1K)])/length(d$ribo_cell1K[!is.na(d$ribo_cell1K)]);q>0.9},]$n
# x[{q<-rank(x)/length(x);q<0.1 | q>=0.9}]
low <- shape_cell1K_et1[names(shape_cell1K_et1) %in% low_ribo]
high <- shape_cell1K_et1[names(shape_cell1K_et1) %in% high_ribo]

# low
u <- unlist(low)
u$log2ratio[u$log2ratio < 0] <- 0
f150 <- subsetByOverlaps(u, cds150, ignore.strand=TRUE)
f150 <- split(f150, seqnames(f150))
meta <- f150[sapply(f150, function(x){length(x) == 150})]

metac <- meta[sapply(meta, function(x){sum(x$TC.control) > 0})]
metac <- Reduce("+", lapply(metac, function(x){x$TC.control/sum(x$TC.control)}))
amplitudes <- abs(fft(metac))
amp <- amplitudes[2:(length(amplitudes)/2+1)]
periods <- 1/periodogram(metac, method = "clone")$freq
d <- data.frame(periods, amp)
p1 <- ggplot(d, aes(x=periods, y=amp)) + geom_line() + ggtitle("cell1K, low ribo, control") + xlim(0, 8)
ggsave(file="/Volumes/USELESS/META/SHAPES_NEW/periodicity/high_low_ribo/cell1K_low_ctr.png", plot=p1)

metat <- meta[sapply(meta, function(x){sum(x$TC.treated) > 0})]
metat <- Reduce("+", lapply(metat, function(x){x$TC.treated/sum(x$TC.treated)}))
amplitudes <- abs(fft(metat))
amp <- amplitudes[2:(length(amplitudes)/2+1)]
periods <- 1/periodogram(metat, method = "clone")$freq
d <- data.frame(periods, amp)
p2 <- ggplot(d, aes(x=periods, y=amp)) + geom_line() + ggtitle("cell1K, low ribo, treated") + xlim(0, 8)
ggsave(file="/Volumes/USELESS/META/SHAPES_NEW/periodicity/high_low_ribo/cell1K_low_tre.png", plot=p2)


# high
u <- unlist(high)
u$log2ratio[u$log2ratio < 0] <- 0
f150 <- subsetByOverlaps(u, cds150, ignore.strand=TRUE)
f150 <- split(f150, seqnames(f150))
meta <- f150[sapply(f150, function(x){length(x) == 150})]

metac <- meta[sapply(meta, function(x){sum(x$TC.control) > 0})]
metac <- Reduce("+", lapply(metac, function(x){x$TC.control/sum(x$TC.control)}))
amplitudes <- abs(fft(metac))
amp <- amplitudes[2:(length(amplitudes)/2+1)]
periods <- 1/periodogram(metac, method = "clone")$freq
d <- data.frame(periods, amp)
p3 <- ggplot(d, aes(x=periods, y=amp)) + geom_line() + ggtitle("cell1K, high ribo, control") + xlim(0, 8)
ggsave(file="/Volumes/USELESS/META/SHAPES_NEW/periodicity/high_low_ribo/cell1K_high_ctr.png", plot=p3)

metat <- meta[sapply(meta, function(x){sum(x$TC.treated) > 0})]
metat <- Reduce("+", lapply(metat, function(x){x$TC.treated/sum(x$TC.treated)}))
amplitudes <- abs(fft(metat))
amp <- amplitudes[2:(length(amplitudes)/2+1)]
periods <- 1/periodogram(metat, method = "clone")$freq
d <- data.frame(periods, amp)
p4 <- ggplot(d, aes(x=periods, y=amp)) + geom_line() + ggtitle("cell1K, high ribo, treated") + xlim(0, 8)
ggsave(file="/Volumes/USELESS/META/SHAPES_NEW/periodicity/high_low_ribo/cell1K_high_tre.png", plot=p4)



################################
#1Kcell
d <- df[df$rna_cell1K > 1 & df$shape_cell1K_et1_log2ratio > 0,]
# get high ribo and low ribo
d <- d[!is.na(d$ribo_cell1K),]
low_ribo <- d[{q<-rank(d$ribo_cell1K[!is.na(d$ribo_cell1K)])/length(d$ribo_cell1K[!is.na(d$ribo_cell1K)]);q<0.3},]$n
high_ribo <- d[{q<-rank(d$ribo_cell1K[!is.na(d$ribo_cell1K)])/length(d$ribo_cell1K[!is.na(d$ribo_cell1K)]);q>0.7},]$n
# x[{q<-rank(x)/length(x);q<0.1 | q>=0.9}]
low <- shape_cell1K_et1[names(shape_cell1K_et1) %in% low_ribo]
high <- shape_cell1K_et1[names(shape_cell1K_et1) %in% high_ribo]

# low
u <- unlist(low)
u$log2ratio[u$log2ratio < 0] <- 0
f150 <- subsetByOverlaps(u, cds150, ignore.strand=TRUE)
f150 <- split(f150, seqnames(f150))
meta <- f150[sapply(f150, function(x){length(x) == 150})]

metac <- meta[sapply(meta, function(x){sum(x$TC.control) > 0})]
metac <- Reduce("+", lapply(metac, function(x){x$TC.control/sum(x$TC.control)}))
amplitudes <- abs(fft(metac))
amp <- amplitudes[2:(length(amplitudes)/2+1)]
periods <- 1/periodogram(metac, method = "clone")$freq
d <- data.frame(periods, amp)
p1 <- ggplot(d, aes(x=periods, y=amp)) + geom_line() + ggtitle("cell1K, low ribo, control") + xlim(0, 8)

metat <- meta[sapply(meta, function(x){sum(x$TC.treated) > 0})]
metat <- Reduce("+", lapply(metat, function(x){x$TC.treated/sum(x$TC.treated)}))
amplitudes <- abs(fft(metat))
amp <- amplitudes[2:(length(amplitudes)/2+1)]
periods <- 1/periodogram(metat, method = "clone")$freq
d <- data.frame(periods, amp)
p2 <- ggplot(d, aes(x=periods, y=amp)) + geom_line() + ggtitle("cell1K, low ribo, treated") + xlim(0, 8)


# high
u <- unlist(high)
u$log2ratio[u$log2ratio < 0] <- 0
f150 <- subsetByOverlaps(u, cds150, ignore.strand=TRUE)
f150 <- split(f150, seqnames(f150))
meta <- f150[sapply(f150, function(x){length(x) == 150})]

metac <- meta[sapply(meta, function(x){sum(x$TC.control) > 0})]
metac <- Reduce("+", lapply(metac, function(x){x$TC.control/sum(x$TC.control)}))
amplitudes <- abs(fft(metac))
amp <- amplitudes[2:(length(amplitudes)/2+1)]
periods <- 1/periodogram(metac, method = "clone")$freq
d <- data.frame(periods, amp)
p3 <- ggplot(d, aes(x=periods, y=amp)) + geom_line() + ggtitle("cell1K, high ribo, control") + xlim(0, 8)

metat <- meta[sapply(meta, function(x){sum(x$TC.treated) > 0})]
metat <- Reduce("+", lapply(metat, function(x){x$TC.treated/sum(x$TC.treated)}))
amplitudes <- abs(fft(metat))
amp <- amplitudes[2:(length(amplitudes)/2+1)]
periods <- 1/periodogram(metat, method = "clone")$freq
d <- data.frame(periods, amp)
p4 <- ggplot(d, aes(x=periods, y=amp)) + geom_line() + ggtitle("cell1K, high ribo, treated") + xlim(0, 8)




