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
txlen <- txlen[txlen$utr5_len > 60,]
rownames(txlen) <- txlen$tx_name
txlen <- txlen[order(rownames(txlen)),]

#start <- GRanges(Rle(rownames(txLengths), rep(1, nrow(txLengths))), IRanges(txLengths$utr5_len-29, width=rep(60, nrow(txLengths))))
start <- GRanges(Rle(rownames(txlen), rep(1, nrow(txlen))), IRanges(txlen$utr5_len-59, width=rep(120, nrow(txlen))))
start <- split(start, seqnames(start))

scale <- c(-60:60)[-61]

###
load(file = "/Volumes/USELESS/DATA/Shape-Seq/shape_cell24.Rsave")
u <- unlist(shape_cell24)
#load(file = "/Volumes/USELESS/test/normalized/u_2-4cell.Rsave")
meta <- subsetByOverlaps(u, start, ignore.strand=TRUE)
meta <- split(meta, seqnames(meta))
meta <- meta[sapply(meta, function(x){length(x) > 0})]
metal <- meta[sapply(meta, function(x){sum(x$log2ratio) > 0})]
cell24 <- metal
metal <- Reduce("+", lapply(metal, function(x){x$log2ratio/sum(x$log2ratio)}))
df <- data.frame(scale = scale, cell24 = metal)

load(file = "/Volumes/USELESS/DATA/Shape-Seq/shape_cell256.Rsave")
u <- unlist(shape_cell256)
#load(file = "/Volumes/USELESS/test/normalized/u_256cell.Rsave")
meta <- subsetByOverlaps(u, start, ignore.strand=TRUE)
meta <- split(meta, seqnames(meta))
meta <- meta[sapply(meta, function(x){length(x) > 0})]
metal <- meta[sapply(meta, function(x){sum(x$log2ratio) > 0})]
cell256 <- metal
metal <- Reduce("+", lapply(metal, function(x){x$log2ratio/sum(x$log2ratio)}))
df$cell256 <- metal

load(file = "/Volumes/USELESS/DATA/Shape-Seq/shape_oblong.Rsave")
u <- unlist(shape_oblong)
#load(file = "/Volumes/USELESS/test/normalized/u_oblong_invivo.Rsave")
meta <- subsetByOverlaps(u, start, ignore.strand=TRUE)
meta <- split(meta, seqnames(meta))
meta <- meta[sapply(meta, function(x){length(x) > 0})]
metal <- meta[sapply(meta, function(x){sum(x$log2ratio) > 0})]
oblong <- metal
metal <- Reduce("+", lapply(metal, function(x){x$log2ratio/sum(x$log2ratio)}))
df$oblong <- metal

load(file = "/Volumes/USELESS/DATA/Shape-Seq/shape_cell1K.Rsave")
u <- unlist(shape_cell1K)
meta <- subsetByOverlaps(u, start, ignore.strand=TRUE)
meta <- split(meta, seqnames(meta))
meta <- meta[sapply(meta, function(x){length(x) > 0})]
metal <- meta[sapply(meta, function(x){sum(x$log2ratio) > 0})]
cell1K <- metal
metal <- Reduce("+", lapply(metal, function(x){x$log2ratio/sum(x$log2ratio)}))
df$cell1K <- metal

load(file = "/Volumes/USELESS/DATA/Shape-Seq/shape_cell1.Rsave")
u <- unlist(shape_cell1)
meta <- subsetByOverlaps(u, start, ignore.strand=TRUE)
meta <- split(meta, seqnames(meta))
meta <- meta[sapply(meta, function(x){length(x) > 0})]
metal <- meta[sapply(meta, function(x){sum(x$log2ratio) > 0})]
cell1 <- metal
metal <- Reduce("+", lapply(metal, function(x){x$log2ratio/sum(x$log2ratio)}))
df$cell1 <- metal



### common tx (in cell24, cell256 & oblong)
cell24 <- cell24[names(cell24) %in% names(cell256)]
cell24 <- cell24[names(cell24) %in% names(oblong)]
cell256 <- cell256[names(cell256) %in% names(cell24)]
oblong <- oblong[names(oblong) %in% names(cell24)]

cell1 <- cell1[names(cell1) %in% names(oblong)]
cell1K <- cell1K[names(cell1K) %in% names(oblong)]
cell1 <- cell1[names(cell1) %in% names(cell1K)]
cell1K <- cell1K[names(cell1K) %in% names(cell1)]
cell24 <- cell24[names(cell24) %in% names(cell1)]
cell256 <- cell256[names(cell256) %in% names(cell1)]
oblong <- oblong[names(oblong) %in% names(cell1)]



df_common <- data.frame(scale = scale,
                        cell1 = Reduce("+", lapply(cell1, function(x){x$log2ratio/sum(x$log2ratio)})),
                        cell24 = Reduce("+", lapply(cell24, function(x){x$log2ratio/sum(x$log2ratio)})),
                        cell256 = Reduce("+", lapply(cell256, function(x){x$log2ratio/sum(x$log2ratio)})),
                        cell1K = Reduce("+", lapply(cell1K, function(x){x$log2ratio/sum(x$log2ratio)})),
                        oblong = Reduce("+", lapply(oblong, function(x){x$log2ratio/sum(x$log2ratio)})))

# ggplot(df, aes(x = scale, y = cell24)) + geom_bar(stat = "identity") + ggtitle("START, cell24")
# ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/start/cell24.png")
# ggplot(df, aes(x = scale, y = cell256)) + geom_bar(stat = "identity") + ggtitle("START, cell256")
# ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/start/cell256.png")
# ggplot(df, aes(x = scale, y = oblong)) + geom_bar(stat = "identity") + ggtitle("START, oblong")
# ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/start/oblong.png")

ggplot(df_common, aes(x = scale, y = cell1)) + geom_bar(stat = "identity") + ggtitle("START, cell1")
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/start/long/common_cell1.png")
ggplot(df_common, aes(x = scale, y = cell24)) + geom_bar(stat = "identity") + ggtitle("START, cell24")
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/start/long/common_cell24.png")
ggplot(df_common, aes(x = scale, y = cell256)) + geom_bar(stat = "identity") + ggtitle("START, cell256")
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/start/long/common_cell256.png")
ggplot(df_common, aes(x = scale, y = cell1K)) + geom_bar(stat = "identity") + ggtitle("START, cell1K")
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/start/long/common_cell1K.png")
ggplot(df_common, aes(x = scale, y = oblong)) + geom_bar(stat = "identity") + ggtitle("START, oblong")
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/start/long/common_oblong.png")


############# new data
for (i in 1:8) {
  
  u <- unlist(eval(parse(text = shape_libs[i])))
  u$log2ratio[u$log2ratio < 0] <- 0
  meta <- subsetByOverlaps(u, start, ignore.strand=TRUE)
  meta <- split(meta, seqnames(meta))
  meta <- meta[sapply(meta, function(x){length(x) > 0})]
  meta <- meta[sapply(meta, function(x){sum(x$log2ratio) > 0})]
  
  assign(paste0("var",i), names(meta))
  meta <- Reduce("+", lapply(meta, function(x){x$log2ratio/sum(x$log2ratio)}))
  d <- data.frame(scale = scale, meta = meta)

  ggplot(d, aes(x = scale, y = meta)) + geom_bar(stat = "identity") + ggtitle(paste(shape_libs[i], "START"))
  fp <- paste0("/Volumes/USELESS/META/SHAPES_NEW/start/p", i, ".png", collapse = NULL)
  ggsave(fp)
  
}


for (i in 1:3) {
  
  u <- unlist(eval(parse(text = shapes_libs[i])))
  meta <- subsetByOverlaps(u, start, ignore.strand=TRUE)
  meta <- split(meta, seqnames(meta))
  meta <- meta[sapply(meta, function(x){length(x) > 0})]
  meta <- meta[sapply(meta, function(x){sum(x$TC) > 0})]
  
  assign(paste0("svar",i), names(meta))
  meta <- Reduce("+", lapply(meta, function(x){x$TC/sum(x$TC)}))
  d <- data.frame(scale = scale, meta = meta)
  
  ggplot(d, aes(x = scale, y = meta)) + geom_bar(stat = "identity") + ggtitle(paste(shapes_libs[i], "START"))
  fp <- paste0("/Volumes/USELESS/META/SHAPES_NEW/start/s", i, ".png", collapse = NULL)
  ggsave(fp)
  
}

# convert p* +append start_shape.png
# convert s* +append start_shapes.png

### get tx in common
common <- var1[var1 %in% var2]
common <- common[common %in% var3]
common <- common[common %in% var4]
common <- common[common %in% var5]
common <- common[common %in% var6]
common <- common[common %in% var7]
common <- common[common %in% var8]
### 605 genes 

commons <- svar1[svar1 %in% svar2]
commons <- commons[commons %in% svar3]
### 1433 genes

for (i in 1:8) {
  
  u <- unlist(eval(parse(text = shape_libs[i])))
  u$log2ratio[u$log2ratio < 0] <- 0
  meta <- subsetByOverlaps(u, start, ignore.strand=TRUE)
  meta <- split(meta, seqnames(meta))
  meta <- meta[names(meta) %in% common]
  meta <- meta[sapply(meta, function(x){sum(x$log2ratio) > 0})]
  
  meta <- Reduce("+", lapply(meta, function(x){x$log2ratio/sum(x$log2ratio)}))
  d <- data.frame(scale = scale, meta = meta)
  
  ggplot(d, aes(x = scale, y = meta)) + geom_bar(stat = "identity") + ggtitle(paste(shape_libs[i], "START, common"))
  fp <- paste0("/Volumes/USELESS/META/SHAPES_NEW/start/common_p", i, ".png", collapse = NULL)
  ggsave(fp)
  
}


for (i in 1:3) {
  
  u <- unlist(eval(parse(text = shapes_libs[i])))
  meta <- subsetByOverlaps(u, start, ignore.strand=TRUE)
  meta <- split(meta, seqnames(meta))
  meta <- meta[names(meta) %in% commons]
  
  meta <- Reduce("+", lapply(meta, function(x){x$TC/sum(x$TC)}))
  d <- data.frame(scale = scale, meta = meta)
  
  ggplot(d, aes(x = scale, y = meta)) + geom_bar(stat = "identity") + ggtitle(paste(shapes_libs[i], "START, common"))
  fp <- paste0("/Volumes/USELESS/META/SHAPES_NEW/start/common_s", i, ".png", collapse = NULL)
  ggsave(fp)
  
}



# convert common_p* +append start_shape_common.png
# convert common_s* +append start_shapes_common.png