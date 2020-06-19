### stall sites
library(GenomicFeatures)
library(ggplot2)


# 2-4Cell
fc <- file("/Volumes/USELESS/META/SHAPES/gwp_2-4Cell.csv")
gwp_2_4 <- strsplit(readLines(fc), ",")
close(fc)

# 256Cell
fc <- file("/Volumes/USELESS/META/SHAPES/gwp_256Cell.csv")
gwp_256 <- strsplit(readLines(fc), ",")
close(fc)

# 1KCell
fc <- file("/Volumes/USELESS/META/SHAPES/gwp_1KCell.csv")
gwp_1K <- strsplit(readLines(fc), ",")
close(fc)

# transform the list of peaks into sth
peaks24 <- sapply(gwp_2_4, function(l) as.numeric(c(l[2:length(l)])))
peaks24 <- setNames(peaks24, sapply(gwp_2_4, function(l) l[[1]]))
names24 <- sapply(gwp_2_4, function(l) l[[1]])

peaks256 <- sapply(gwp_256, function(l) as.numeric(c(l[2:length(l)])))
peaks256 <- setNames(peaks256, sapply(gwp_256, function(l) l[[1]]))
names256 <- sapply(gwp_256, function(l) l[[1]])

peaks1K <- sapply(gwp_1K, function(l) as.numeric(c(l[2:length(l)])))
peaks1K <- setNames(peaks1K, sapply(gwp_1K, function(l) l[[1]]))
names1K <- sapply(gwp_1K, function(l) l[[1]])

names24 <- as.character(unlist(mapply(function(x,y){rep(x, length(y))}, names24, peaks24)))
df24 <- data.frame(n = names24, peak = as.vector(unlist(peaks24)))
df24$utr5_len <- txLengths[df24$n,]$utr5_len

names256 <- as.character(unlist(mapply(function(x,y){rep(x, length(y))}, names256, peaks256)))
df256 <- data.frame(n = names256, peak = as.vector(unlist(peaks256)))
df256$utr5_len <- txLengths[df256$n,]$utr5_len

names1K <- as.character(unlist(mapply(function(x,y){rep(x, length(y))}, names1K, peaks1K)))
df1K <- data.frame(n = names1K, peak = as.vector(unlist(peaks1K)))
df1K$utr5_len <- txLengths[df1K$n,]$utr5_len

### utr5_len + peak - 29 : utr5_len + peak + 30
## subset on control and treated separately
txdb_can <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/Danio_rerio.Zv9.79.gtf", format = "gtf")
txLengths <- transcriptLengths(txdb_can, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
rownames(txLengths) <- txLengths$tx_name
txLengths <- txLengths[order(rownames(txLengths)),]


stall <- GRanges(Rle(df24$n, rep(1, nrow(df24))), IRanges(df24$utr5_len+df24$peak-29, width=rep(60, nrow(df24))))
#stall <- split(stall, seqnames(stall))

scale <- c(-30:30)[-31]

####### SHAPE
load(file = "/Volumes/USELESS/DATA/Shape-Seq/shape_cell24.Rsave")
load(file = "/Volumes/USELESS/DATA/Shape-Seq/shape_cell256.Rsave")
load(file = "/Volumes/USELESS/DATA/Shape-Seq/shape_cell1K.Rsave")

########
u24 <- unlist(shape_cell24)
meta <- subsetByOverlaps(u24, stall, ignore.strand=TRUE)
meta <- split(meta, seqnames(meta))
meta <- meta[sapply(meta, function(x){length(x) == 60})]

metac <- meta[sapply(meta, function(x){sum(x$TC.control) > 0})]
metac <- Reduce("+", lapply(metac, function(x){x$TC.control/sum(x$TC.control)}))
metat <- meta[sapply(meta, function(x){sum(x$TC.treated) > 0})]
metat <- Reduce("+", lapply(metat, function(x){x$TC.treated/sum(x$TC.treated)}))

df <- data.frame(scale = scale, control = metac, treated = metat)
ggplot(df, aes(x = scale, y = control)) + geom_bar(stat = "identity") + ggtitle("STALL, cell24, control")
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/stall/control_cell24.png")
ggplot(df, aes(x = scale, y = treated)) + geom_bar(stat = "identity") + ggtitle("STALL, cell24, treated")
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/stall/treated_cell24.png")


########
u256 <- unlist(shape_cell256)
meta <- subsetByOverlaps(u256, stall, ignore.strand=TRUE)
meta <- split(meta, seqnames(meta))
meta <- meta[sapply(meta, function(x){length(x) == 60})]

metac <- meta[sapply(meta, function(x){sum(x$TC.control) > 0})]
metac <- Reduce("+", lapply(metac, function(x){x$TC.control/sum(x$TC.control)}))
metat <- meta[sapply(meta, function(x){sum(x$TC.treated) > 0})]
metat <- Reduce("+", lapply(metat, function(x){x$TC.treated/sum(x$TC.treated)}))

df <- data.frame(scale = scale, control = metac, treated = metat)
ggplot(df, aes(x = scale, y = control)) + geom_bar(stat = "identity") + ggtitle("STALL, cell256, control")
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/stall/control_cell256.png")
ggplot(df, aes(x = scale, y = treated)) + geom_bar(stat = "identity") + ggtitle("STALL, cell256, treated")
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/stall/treated_cell256.png")


########
u1K <- unlist(shape_cell1K)
meta <- subsetByOverlaps(u1K, stall, ignore.strand=TRUE)
meta <- split(meta, seqnames(meta))
meta <- meta[sapply(meta, function(x){length(x) == 60})]

metac <- meta[sapply(meta, function(x){sum(x$TC.control) > 0})]
metac <- Reduce("+", lapply(metac, function(x){x$TC.control/sum(x$TC.control)}))
metat <- meta[sapply(meta, function(x){sum(x$TC.treated) > 0})]
metat <- Reduce("+", lapply(metat, function(x){x$TC.treated/sum(x$TC.treated)}))

df <- data.frame(scale = scale, control = metac, treated = metat)
ggplot(df, aes(x = scale, y = control)) + geom_bar(stat = "identity") + ggtitle("STALL, cell1K, control")
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/stall/control_cell1K.png")
ggplot(df, aes(x = scale, y = treated)) + geom_bar(stat = "identity") + ggtitle("STALL, cell1K, treated")
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/stall/treated_cell1K.png")
