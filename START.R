### START
library(GenomicFeatures)
library(ggplot2)
library(plyr)

txdb_can <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/danio_rerio/Danio_rerio.GRCz10.84.gtf", format = "gtf")
txLengths <- transcriptLengths(txdb_can, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
rownames(txLengths) <- txLengths$tx_name
txLengths <- txLengths[order(rownames(txLengths)),]

# get longest transcript per gene
txlen <- arrange(txLengths, gene_id, desc(tx_len))
txlen <- txlen[!duplicated(txlen$gene_id),]

### tx with utr5 longer than 30
#txLengths <- txLengths[txLengths$utr5_len > 30,]
txlen <- txlen[txlen$utr5_len > 30,]
rownames(txlen) <- txlen$tx_name
txlen <- txlen[order(rownames(txlen)),]

#start <- GRanges(Rle(rownames(txLengths), rep(1, nrow(txLengths))), IRanges(txLengths$utr5_len-29, width=rep(60, nrow(txLengths))))
start <- GRanges(Rle(rownames(txlen), rep(1, nrow(txlen))), IRanges(txlen$utr5_len-29, width=rep(60, nrow(txlen))))
start <- split(start, seqnames(start))

scale <- c(-30:30)[-31]

###
load(file = "/Volumes/USELESS/test/normalized/u_2-4cell.Rsave")
meta <- subsetByOverlaps(u, start, ignore.strand=TRUE)
meta <- split(meta, seqnames(meta))
meta <- meta[sapply(meta, function(x){length(x) > 0})]
metal <- meta[sapply(meta, function(x){sum(x$log2ratio) > 0})]
cell24 <- metal
metal <- Reduce("+", lapply(metal, function(x){x$log2ratio/sum(x$log2ratio)}))
df <- data.frame(scale = scale, cell24 = metal)

load(file = "/Volumes/USELESS/test/normalized/u_256cell.Rsave")
meta <- subsetByOverlaps(u, start, ignore.strand=TRUE)
meta <- split(meta, seqnames(meta))
meta <- meta[sapply(meta, function(x){length(x) > 0})]
metal <- meta[sapply(meta, function(x){sum(x$log2ratio) > 0})]
cell256 <- metal
metal <- Reduce("+", lapply(metal, function(x){x$log2ratio/sum(x$log2ratio)}))
df$cell256 <- metal

load(file = "/Volumes/USELESS/test/normalized/u_oblong_invivo.Rsave")
meta <- subsetByOverlaps(u, start, ignore.strand=TRUE)
meta <- split(meta, seqnames(meta))
meta <- meta[sapply(meta, function(x){length(x) > 0})]
metal <- meta[sapply(meta, function(x){sum(x$log2ratio) > 0})]
oblong <- metal
metal <- Reduce("+", lapply(metal, function(x){x$log2ratio/sum(x$log2ratio)}))
df$oblong <- metal


### common tx (in cell24, cell256 & oblong)
cell24 <- cell24[names(cell24) %in% names(cell256)]
cell24 <- cell24[names(cell24) %in% names(oblong)]
cell256 <- cell256[names(cell256) %in% names(cell24)]
oblong <- oblong[names(oblong) %in% names(cell24)]

df_common <- data.frame(scale = scale,
                        cell24 = Reduce("+", lapply(cell24, function(x){x$log2ratio/sum(x$log2ratio)})),
                        cell256 = Reduce("+", lapply(cell256, function(x){x$log2ratio/sum(x$log2ratio)})),
                        oblong = Reduce("+", lapply(oblong, function(x){x$log2ratio/sum(x$log2ratio)})))

# ggplot(df, aes(x = scale, y = cell24)) + geom_bar(stat = "identity") + ggtitle("START, cell24")
# ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/start/cell24.png")
# ggplot(df, aes(x = scale, y = cell256)) + geom_bar(stat = "identity") + ggtitle("START, cell256")
# ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/start/cell256.png")
# ggplot(df, aes(x = scale, y = oblong)) + geom_bar(stat = "identity") + ggtitle("START, oblong")
# ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/start/oblong.png")

ggplot(df_common, aes(x = scale, y = cell24)) + geom_bar(stat = "identity") + ggtitle("START, cell24")
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/start/common_cell24.png")
ggplot(df_common, aes(x = scale, y = cell256)) + geom_bar(stat = "identity") + ggtitle("START, cell256")
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/start/common_cell256.png")
ggplot(df_common, aes(x = scale, y = oblong)) + geom_bar(stat = "identity") + ggtitle("START, oblong")
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/start/common_oblong.png")
