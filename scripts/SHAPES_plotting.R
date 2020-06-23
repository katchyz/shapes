library(GenomicFeatures)
library(ggplot2)
library(GeneCycle)
library(reshape2)

load("/Volumes/USELESS/OUT/SHAPES_out/rsave/gc_dtcr_2_4.Rsave")
load("/Volumes/USELESS/OUT/SHAPES_out/rsave/gc_dtcr_256.Rsave")

### get transcript dtabase
txdb_can <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/danio_rerio/zebrafish_GRCz10.gtf", format = "gtf")
cds <- cdsBy(txdb_can, by="tx", use.names=TRUE)
fiveUTR <- fiveUTRsByTranscript(txdb_can, use.names=TRUE)
threeUTR <- threeUTRsByTranscript(txdb_can, use.names=TRUE)
exons <- exonsBy(txdb_can, by="tx", use.names=TRUE)
##############

###### five_end ######
exon_2_4_list <- split(gc_dtcr_2_4, gc_dtcr_2_4$trnames)
exon_2_4_list <- exon_2_4_list[lengths(exon_2_4_list) > 50]
### get 5'end of transcripts
meta <- Reduce("+", lapply(exon_2_4_list, function(x){x[1:50]$dtcr}))
scale <- c(1:50)
df <- data.frame(scale,meta)
p1 <- ggplot(df, aes(x=scale, y=meta)) + geom_bar(stat = "identity") + xlab("position from 5'end of transcript") + ylab("sum(reactivities)") + ggtitle("2-4cell")
ggsave(file = "/Volumes/USELESS/META/SHAPES/five_end/p_2_4.png", plot = p1)

exon_256_list <- split(gc_dtcr_256, gc_dtcr_256$trnames)
exon_256_list <- exon_256_list[lengths(exon_256_list) > 50]
### get 5'end of transcripts
meta <- Reduce("+", lapply(exon_256_list, function(x){x[1:50]$dtcr}))
scale <- c(1:50)
df <- data.frame(scale,meta)
p2 <- ggplot(df, aes(x=scale, y=meta)) + geom_bar(stat = "identity") + xlab("position from 5'end of transcript") + ylab("sum(reactivities)") + ggtitle("256cell")
ggsave(file = "/Volumes/USELESS/META/SHAPES/five_end/p_256.png", plot = p2)

###### periodicity ######
cds_2_4 <- subsetByOverlaps(gc_dtcr_2_4, cds)
cds_2_4_list <- split(cds_2_4, cds_2_4$trnames)
cds_2_4_list <- cds_2_4_list[lengths(cds_2_4_list) > 150]
### get 5'end of CDSs
meta <- Reduce("+", lapply(cds_2_4_list, function(x){x[1:150]$dtcr}))
amplitudes <- abs(fft(meta))
amp <- amplitudes[2:(length(amplitudes)/2+1)]
periods <- 1/periodogram(meta, method = "clone")$freq
df <- data.frame(periods, amp)
p3 <- ggplot(df, aes(x=periods, y=amp)) + geom_line() + ggtitle("FFT, 2-4cell") + xlim(0, 10)
ggsave(file = "/Volumes/USELESS/META/SHAPES/periodicity/p_2_4.png", plot = p3)

cds_256 <- subsetByOverlaps(gc_dtcr_256, cds)
cds_256_list <- split(cds_256, cds_256$trnames)
cds_256_list <- cds_256_list[lengths(cds_256_list) > 150]
### get 5'end of CDSs
meta <- Reduce("+", lapply(cds_256_list, function(x){x[1:150]$dtcr}))
amplitudes <- abs(fft(meta))
amp <- amplitudes[2:(length(amplitudes)/2+1)]
periods <- 1/periodogram(meta, method = "clone")$freq
df <- data.frame(periods, amp)
p4 <- ggplot(df, aes(x=periods, y=amp)) + geom_line() + ggtitle("FFT, 256cell") + xlim(0, 10)
ggsave(file = "/Volumes/USELESS/META/SHAPES/periodicity/p_256.png", plot = p4)

###### densities ######
dens_cds <- sapply(cds_2_4_list, function(x){sum(x$dtcr)/length(x)})

five_2_4 <- subsetByOverlaps(gc_dtcr_2_4, fiveUTR)
five_2_4_list <- split(five_2_4, five_2_4$trnames)
dens_five <- sapply(five_2_4_list, function(x){sum(x$dtcr)/length(x)})

three_2_4 <- subsetByOverlaps(gc_dtcr_2_4, threeUTR)
three_2_4_list <- split(three_2_4, three_2_4$trnames)
dens_three <- sapply(three_2_4_list, function(x){sum(x$dtcr)/length(x)})

df <- data.frame(dens_five, dens_cds, dens_three)
colnames(df) <- c("leaders", "CDSs", "trailers")
df$scale <- c(1:length(dens_cds))
d <- melt(df, id.var="scale")
d <- d[complete.cases(d), ]
p5 <- ggplot(d, aes(x=value, fill=variable)) + geom_density(alpha=.3) + ggtitle("2-4cell") + xlim(0, 0.2)
p6 <- ggplot(d, aes(x=value, fill=variable)) + geom_density(alpha=.3) + ggtitle("2-4cell") + scale_x_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES/distributions/p_2_4.png", plot = p5)
ggsave(file = "/Volumes/USELESS/META/SHAPES/distributions/p_2_4log.png", plot = p6)


dens_cds <- sapply(cds_256_list, function(x){sum(x$dtcr)/length(x)})

five_256 <- subsetByOverlaps(gc_dtcr_256, fiveUTR)
five_256_list <- split(five_256, five_256$trnames)
dens_five <- sapply(five_256_list, function(x){sum(x$dtcr)/length(x)})

three_256 <- subsetByOverlaps(gc_dtcr_256, threeUTR)
three_256_list <- split(three_256, three_256$trnames)
dens_three <- sapply(three_256_list, function(x){sum(x$dtcr)/length(x)})

df <- data.frame(dens_five, dens_cds, dens_three)
colnames(df) <- c("leaders", "CDSs", "trailers")
df$scale <- c(1:length(dens_cds))
d <- melt(df, id.var="scale")
d <- d[complete.cases(d), ]
p7 <- ggplot(d, aes(x=value, fill=variable)) + geom_density(alpha=.3) + ggtitle("256cell") + xlim(0, 0.2)
p8 <- ggplot(d, aes(x=value, fill=variable)) + geom_density(alpha=.3) + ggtitle("256cell") + scale_x_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES/distributions/p_256.png", plot = p7)
ggsave(file = "/Volumes/USELESS/META/SHAPES/distributions/p_256log.png", plot = p8)

###### start ######
# get cdss where leaders are longer than 30
cds_2_4_list <- cds_2_4_list[lengths(five_2_4_list) > 30]
cds_2_4_list <- cds_2_4_list[lengths(cds_2_4_list) > 33]
five_2_4_list <- five_2_4_list[lengths(five_2_4_list) > 30]
five_2_4_list <- five_2_4_list[lengths(cds_2_4_list) > 33]
# get starts of cdss
meta_cds <- Reduce("+", lapply(cds_2_4_list, function(x){x[1:33]$dtcr}))
# get ends of leaders
meta_five <- Reduce("+", lapply(five_2_4_list, function(x){x[(length(x)-29):length(x)]$dtcr}))
# stitch together and plot
scale <- c(-30:33)[-31]
meta <- c(meta_five, meta_cds)
df <- data.frame(scale,meta)
p9 <- ggplot(df, aes(x=scale, y=meta)) + geom_bar(stat = "identity") + xlab("position from start codon") + ylab("sum(reactivities)") + ggtitle("2-4cell")
ggsave(file = "/Volumes/USELESS/META/SHAPES/start/p_2_4.png", plot = p9)

# get cdss where leaders are longer than 30
cds_256_list <- cds_256_list[lengths(five_256_list) > 30]
cds_256_list <- cds_256_list[lengths(cds_256_list) > 33]
five_256_list <- five_256_list[lengths(five_256_list) > 30]
five_256_list <- five_256_list[lengths(cds_256_list) > 33]
# get starts of cdss
meta_cds <- Reduce("+", lapply(cds_256_list, function(x){x[1:33]$dtcr}))
# get ends of leaders
meta_five <- Reduce("+", lapply(five_256_list, function(x){x[(length(x)-29):length(x)]$dtcr}))
# stitch together and plot
scale <- c(-30:33)[-31]
meta <- c(meta_five, meta_cds)
df <- data.frame(scale,meta)
p10 <- ggplot(df, aes(x=scale, y=meta)) + geom_bar(stat = "identity") + xlab("position from start codon") + ylab("sum(reactivities)") + ggtitle("256cell")
ggsave(file = "/Volumes/USELESS/META/SHAPES/start/p_256.png", plot = p10)


###### stop ######
# get cdss where trailers are longer than 30
cds_2_4_list <- cds_2_4_list[lengths(three_2_4_list) > 30]
cds_2_4_list <- cds_2_4_list[lengths(cds_2_4_list) > 33]
three_2_4_list <- three_2_4_list[lengths(three_2_4_list) > 30]
three_2_4_list <- three_2_4_list[lengths(cds_2_4_list) > 33]
# get ends of cdss
meta_cds <- Reduce("+", lapply(cds_2_4_list, function(x){x[(length(x)-32):length(x)]$dtcr}))
# get starts of trailers
meta_three <- Reduce("+", lapply(three_2_4_list, function(x){x[1:30]$dtcr}))
# stitch together and plot
scale <- c(-33:30)[-34]
meta <- c(meta_cds, meta_three)
df <- data.frame(scale,meta)
p11 <- ggplot(df, aes(x=scale, y=meta)) + geom_bar(stat = "identity") + xlab("position from stop codon") + ylab("sum(reactivities)") + ggtitle("2-4cell")
ggsave(file = "/Volumes/USELESS/META/SHAPES/stop/p_2_4.png", plot = p11)

# get cdss where trailers are longer than 30
cds_256_list <- cds_256_list[lengths(three_256_list) > 30]
cds_256_list <- cds_256_list[lengths(cds_256_list) > 33]
three_256_list <- three_256_list[lengths(three_256_list) > 30]
three_256_list <- three_256_list[lengths(cds_256_list) > 33]
# get ends of cdss
meta_cds <- Reduce("+", lapply(cds_256_list, function(x){x[(length(x)-32):length(x)]$dtcr}))
# get starts of trailers
meta_three <- Reduce("+", lapply(three_256_list, function(x){x[1:30]$dtcr}))
# stitch together and plot
scale <- c(-33:30)[-34]
meta <- c(meta_cds, meta_three)
df <- data.frame(scale,meta)
p12 <- ggplot(df, aes(x=scale, y=meta)) + geom_bar(stat = "identity") + xlab("position from stop codon") + ylab("sum(reactivities)") + ggtitle("256cell")
ggsave(file = "/Volumes/USELESS/META/SHAPES/stop/p_256.png", plot = p12)

###### TE ######
data = read.csv("/Users/kasia/Documents/PhD/TE.csv", header = TRUE)
rownames(data) <- data$X.transcript_id

names(dens_cds) <- sapply(names(dens_cds), function(x){substr(x,1,18)})
names(dens_five) <- sapply(names(dens_five), function(x){substr(x,1,18)})
data <- data[rownames(data) %in% names(dens_cds),]
data <- data[rownames(data) %in% names(dens_five),]
dens_cds <- dens_cds[names(dens_cds) %in% rownames(data)]
dens_five <- dens_five[names(dens_five) %in% rownames(data)]

#shape_te_cds <- merge(dens_cds, data.frame(rownames(data), data$X02_256cell_gene_te))
shape_te_cds <- data.frame(as.numeric(dens_cds), data$X01_2to4cell_gene_te)
colnames(shape_te_cds) <- c("cdss", "te_cds")
#shape_te_five <- merge(dens_five, data.frame(rownames(data), data$X02_256cell_lead_te))
shape_te_five <- data.frame(as.numeric(dens_five), data$X01_2to4cell_lead_te)
colnames(shape_te_five) <- c("leaders", "te_leader")

p15 <- ggplot(shape_te_five, aes(x=leaders, y=te_leader)) + geom_point() + xlab("sum fo reactivities/length") + ylab("TE") + ggtitle("leaders, 2_4cell") + scale_x_log10() + scale_y_log10()
p16 <- ggplot(shape_te_cds, aes(x=cdss, y=te_cds)) + geom_point() + xlab("sum fo reactivities/length") + ylab("TE") + ggtitle("CDSs, 2_4cell") + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES/te/p_lead_log_2_4.png", plot = p15)
ggsave(file = "/Volumes/USELESS/META/SHAPES/te/p_cds_log_2_4.png", plot = p16)

cds_256 <- subsetByOverlaps(gc_dtcr_256, cds)
cds_256_list <- split(cds_256, cds_256$trnames)
dens_cds <- sapply(cds_256_list, function(x){sum(x$dtcr)/length(x)})
five_256 <- subsetByOverlaps(gc_dtcr_256, fiveUTR)
five_256_list <- split(five_256, five_256$trnames)
dens_five <- sapply(five_256_list, function(x){sum(x$dtcr)/length(x)})
data = read.csv("/Users/kasia/Documents/PhD/TE.csv", header = TRUE)
rownames(data) <- data$X.transcript_id
names(dens_cds) <- sapply(names(dens_cds), function(x){substr(x,1,18)})
names(dens_five) <- sapply(names(dens_five), function(x){substr(x,1,18)})
data <- data[rownames(data) %in% names(dens_cds),]
data <- data[rownames(data) %in% names(dens_five),]
dens_cds <- dens_cds[names(dens_cds) %in% rownames(data)]
dens_five <- dens_five[names(dens_five) %in% rownames(data)]
shape_te_cds <- data.frame(as.numeric(dens_cds), data$X02_256cell_gene_te)
colnames(shape_te_cds) <- c("cdss", "te_cds")
shape_te_five <- data.frame(as.numeric(dens_five), data$X02_256cell_lead_te)
colnames(shape_te_five) <- c("leaders", "te_leader")

p17 <- ggplot(shape_te_five, aes(x=leaders, y=te_leader)) + geom_point() + xlab("sum fo reactivities/length") + ylab("TE") + ggtitle("leaders, 256cell") + scale_x_log10() + scale_y_log10()
p18 <- ggplot(shape_te_cds, aes(x=cdss, y=te_cds)) + geom_point() + xlab("sum fo reactivities/length") + ylab("TE") + ggtitle("CDSs, 256cell") + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES/te/p_lead_log_256.png", plot = p17)
ggsave(file = "/Volumes/USELESS/META/SHAPES/te/p_cds_log_256.png", plot = p18)
###### hindrance ######

# sum of leaders vs TE of CDS
sum_five <- sapply(five_2_4_list, function(x){sum(x$dtcr)})
names(sum_five) <- sapply(names(sum_five), function(x){substr(x,1,18)})
sum_five <- sum_five[names(sum_five) %in% rownames(data)]
data <- data[rownames(data) %in% names(sum_five),]
shape_five_tecds <- data.frame(as.numeric(sum_five), data$X01_2to4cell_gene_te)
colnames(shape_five_tecds) <- c("leader", "te_cds")
p19 <- ggplot(shape_five_tecds, aes(x=leader, y=te_cds, color=te_cds)) + geom_point() + xlab("sum of reactivities, leaders") + ylab("TE of CDS") + ggtitle("hindrance, 2_4cell") + scale_x_log10() + scale_y_log10() #+ scale_colour_gradient()
ggsave(file = "/Volumes/USELESS/META/SHAPES/te/p_lead_te_cds_2_4_log.png", plot = p19)

sum_five <- sapply(five_256_list, function(x){sum(x$dtcr)})
names(sum_five) <- sapply(names(sum_five), function(x){substr(x,1,18)})
sum_five <- sum_five[names(sum_five) %in% rownames(data)]
data <- data[rownames(data) %in% names(sum_five),]
shape_five_tecds <- data.frame(as.numeric(sum_five), data$X02_256cell_gene_te)
colnames(shape_five_tecds) <- c("leader", "te_cds")
p20 <- ggplot(shape_five_tecds, aes(x=leader, y=te_cds)) + geom_point() + xlab("sum fo reactivities, leaders") + ylab("TE of CDS") + ggtitle("hindrance, 256cell") #+ scale_x_log10() + scale_y_log10()# + scale_colour_gradient()
ggsave(file = "/Volumes/USELESS/META/SHAPES/te/p_lead_te_cds_256_log.png", plot = p20)

###### exons ######
exon_2_4_list <- exon_2_4_list[names(exon_2_4_list) %in% names(exons)]
exon_256_list <- exon_256_list[names(exon_256_list) %in% names(exons)]

splice_sites <- function(tr) {
  meta <- rep(0,61)
  ex <- 0
  l <- width(ranges(exons[[tr]]))
  if (length(l) > 1) {
    for (i in (1:(length(l)-1))) {
      ex <- ex + l[i]
      if ((l[i]>30)&(l[i+1]>30)) {
        meta <- meta + exon_2_4_list[[tr]]$dtcr[(ex-30):(ex+30)]
      }
    }
  }
  return(meta)
}

#meta_exon <- Reduce("+", lapply(names(exon_2_4_list), function(x){splice_sites(x)}))
meta_exon <- lapply(names(exon_2_4_list), function(x){splice_sites(x)})
load("/Volumes/USELESS/META/SHAPES/exons/meta_exon2_4.Rsave")
meta_exon <- sapply(meta_exon, function(x){x[!any(is.na(x))]})
meta_exon <- meta_exon[sapply(meta_exon, function(x){length(x) == 61})]
meta <- Reduce("+", meta_exon)

scale <- c(-30:31)[-31]
df <- data.frame(scale,meta)
p21 <- ggplot(df, aes(x=scale, y=meta)) + geom_bar(stat = "identity") + xlab("position from splice site") + ylab("sum(reactivities)") + ggtitle("2-4cell")
ggsave(file = "/Volumes/USELESS/META/SHAPES/exons/p_2_4.png", plot = p21)
p22 <- ggplot(df, aes(x=scale, y=meta)) + geom_bar(stat = "identity") + xlab("position from splice site") + ylab("sum(reactivities)") + ggtitle("256cell")
ggsave(file = "/Volumes/USELESS/META/SHAPES/exons/p_256.png", plot = p22)


#### normalized to 1
# norm_dtcr_256 <- lapply(exon_256_list, function(x) if (sum(x$dtcr) > 0) {x$dtcr/sum(x$dtcr)} else {x$dtcr})
# norm_256$dtcr <- unlist(norm_dtcr_256)
# norm_256_list <- split(norm_256, norm_256$trnames)
# save(norm_256_list, file="/Volumes/USELESS/META/SHAPES/norm_256_list.Rsave")