### STOP

library(GenomicFeatures)
library(ggplot2)
library(reshape2)
library(plyr)
library(seqinr)

txdb_can <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/danio_rerio/Danio_rerio.GRCz10.84.gtf", format = "gtf")
txLengths <- transcriptLengths(txdb_can, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
rownames(txLengths) <- txLengths$tx_name
txLengths <- txLengths[order(rownames(txLengths)),]

# get longest transcript per gene
#txlen <- arrange(txLengths, gene_id, desc(tx_len))
#txlen <- txlen[!duplicated(txlen$gene_id),]

### tx with utr5 longer than 30
txLengths <- txLengths[txLengths$utr3_len > 30,]
#txlen <- txlen[txlen$utr3_len > 30,]
#rownames(txlen) <- txlen$tx_name
#txlen <- txlen[order(rownames(txlen)),]

stop <- GRanges(Rle(rownames(txLengths), rep(1, nrow(txLengths))), IRanges(txLengths$utr5_len+txLengths$cds_len-29, width=rep(60, nrow(txLengths))))
#stop <- GRanges(Rle(rownames(txlen), rep(1, nrow(txlen))), IRanges(txlen$utr5_len+txlen$cds_len-29, width=rep(60, nrow(txlen))))
stop <- split(stop, seqnames(stop))

scale <- c(-30:30)[-31]

### stop codon
fasta_cds <- read.fasta("/Volumes/USELESS/DATA/fasta/danio_rerio/GRCz10/cds/Danio_rerio.GRCz10.cds.all.fa")
names(fasta_cds) <- sapply(names(fasta_cds), function(x){substr(x,1,18)})
fasta_cds <- fasta_cds[unlist(lapply(fasta_cds, function(x){length(x) %% 3 == 0}))]

taa <- sort(names(fasta_cds[unlist(lapply(fasta_cds, function(x){paste(x[(length(x)-2):length(x)], collapse = "") == "taa"}))]))
tga <- sort(names(fasta_cds[unlist(lapply(fasta_cds, function(x){paste(x[(length(x)-2):length(x)], collapse = "") == "tga"}))]))
tag <- sort(names(fasta_cds[unlist(lapply(fasta_cds, function(x){paste(x[(length(x)-2):length(x)], collapse = "") == "tag"}))]))

# cell24
load(file = "/Volumes/USELESS/test/normalized/u_2-4cell.Rsave")
meta <- subsetByOverlaps(u, stop, ignore.strand=TRUE)
meta <- split(meta, seqnames(meta))
meta <- meta[sapply(meta, function(x){length(x) > 0})]
metal <- meta[sapply(meta, function(x){sum(x$log2ratio) > 0})]
cell24 <- metal

cell24_taa <- cell24[names(cell24) %in% taa]
cell24_tga <- cell24[names(cell24) %in% tga]
cell24_tag <- cell24[names(cell24) %in% tag]

df_cell24 <- data.frame(scale = scale,
                 cell24_taa = Reduce("+", lapply(cell24_taa, function(x){x$log2ratio/sum(x$log2ratio)})),
                 cell24_tga = Reduce("+", lapply(cell24_tga, function(x){x$log2ratio/sum(x$log2ratio)})),
                 cell24_tag = Reduce("+", lapply(cell24_tag, function(x){x$log2ratio/sum(x$log2ratio)})))

ggplot(df_cell24, aes(x=scale, y=cell24_taa)) + geom_bar(stat = "identity") + ggtitle("STOP, cell24, TAA")
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/stop/cell24_taa.png")
ggplot(df_cell24, aes(x=scale, y=cell24_tga)) + geom_bar(stat = "identity") + ggtitle("STOP, cell24, TGA")
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/stop/cell24_tga.png")
ggplot(df_cell24, aes(x=scale, y=cell24_tag)) + geom_bar(stat = "identity") + ggtitle("STOP, cell24, TAG")
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/stop/cell24_tag.png")

# cell256
load(file = "/Volumes/USELESS/test/normalized/u_256cell.Rsave")
meta <- subsetByOverlaps(u, stop, ignore.strand=TRUE)
meta <- split(meta, seqnames(meta))
meta <- meta[sapply(meta, function(x){length(x) > 0})]
metal <- meta[sapply(meta, function(x){sum(x$log2ratio) > 0})]
cell256 <- metal

cell256_taa <- cell256[names(cell256) %in% taa]
cell256_tga <- cell256[names(cell256) %in% tga]
cell256_tag <- cell256[names(cell256) %in% tag]

df_cell256 <- data.frame(scale = scale,
                        cell256_taa = Reduce("+", lapply(cell256_taa, function(x){x$log2ratio/sum(x$log2ratio)})),
                        cell256_tga = Reduce("+", lapply(cell256_tga, function(x){x$log2ratio/sum(x$log2ratio)})),
                        cell256_tag = Reduce("+", lapply(cell256_tag, function(x){x$log2ratio/sum(x$log2ratio)})))

ggplot(df_cell256, aes(x=scale, y=cell256_taa)) + geom_bar(stat = "identity") + ggtitle("STOP, cell256, TAA")
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/stop/cell256_taa.png")
ggplot(df_cell256, aes(x=scale, y=cell256_tga)) + geom_bar(stat = "identity") + ggtitle("STOP, cell256, TGA")
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/stop/cell256_tga.png")
ggplot(df_cell256, aes(x=scale, y=cell256_tag)) + geom_bar(stat = "identity") + ggtitle("STOP, cell256, TAG")
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/stop/cell256_tag.png")

# oblong
load(file = "/Volumes/USELESS/test/normalized/u_oblong_invivo.Rsave")
meta <- subsetByOverlaps(u, stop, ignore.strand=TRUE)
meta <- split(meta, seqnames(meta))
meta <- meta[sapply(meta, function(x){length(x) > 0})]
metal <- meta[sapply(meta, function(x){sum(x$log2ratio) > 0})]
oblong <- metal

oblong_taa <- oblong[names(oblong) %in% taa]
oblong_tga <- oblong[names(oblong) %in% tga]
oblong_tag <- oblong[names(oblong) %in% tag]

df_oblong <- data.frame(scale = scale,
                        oblong_taa = Reduce("+", lapply(oblong_taa, function(x){x$log2ratio/sum(x$log2ratio)})),
                        oblong_tga = Reduce("+", lapply(oblong_tga, function(x){x$log2ratio/sum(x$log2ratio)})),
                        oblong_tag = Reduce("+", lapply(oblong_tag, function(x){x$log2ratio/sum(x$log2ratio)})))

ggplot(df_oblong, aes(x=scale, y=oblong_taa)) + geom_bar(stat = "identity") + ggtitle("STOP, oblong, TAA")
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/stop/oblong_taa.png")
ggplot(df_oblong, aes(x=scale, y=oblong_tga)) + geom_bar(stat = "identity") + ggtitle("STOP, oblong, TGA")
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/stop/oblong_tga.png")
ggplot(df_oblong, aes(x=scale, y=oblong_tag)) + geom_bar(stat = "identity") + ggtitle("STOP, oblong, TAG")
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/stop/oblong_tag.png")


######################
#### common, OTPG ####
######################
cell24 <- cell24[names(cell24) %in% names(cell256)]
cell24 <- cell24[names(cell24) %in% names(oblong)]
cell256 <- cell256[names(cell256) %in% names(cell24)]
oblong <- oblong[names(oblong) %in% names(cell24)]

#########
ggplot(df_cell24, aes(x=scale, y=cell24_taa)) + geom_bar(stat = "identity") + ggtitle("STOP, cell24, TAA")
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/stop/common_cell24_taa.png")
ggplot(df_cell24, aes(x=scale, y=cell24_tga)) + geom_bar(stat = "identity") + ggtitle("STOP, cell24, TGA")
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/stop/common_cell24_tga.png")
ggplot(df_cell24, aes(x=scale, y=cell24_tag)) + geom_bar(stat = "identity") + ggtitle("STOP, cell24, TAG")
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/stop/common_cell24_tag.png")

ggplot(df_cell256, aes(x=scale, y=cell256_taa)) + geom_bar(stat = "identity") + ggtitle("STOP, cell256, TAA")
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/stop/common_cell256_taa.png")
ggplot(df_cell256, aes(x=scale, y=cell256_tga)) + geom_bar(stat = "identity") + ggtitle("STOP, cell256, TGA")
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/stop/common_cell256_tga.png")
ggplot(df_cell256, aes(x=scale, y=cell256_tag)) + geom_bar(stat = "identity") + ggtitle("STOP, cell256, TAG")
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/stop/common_cell256_tag.png")

ggplot(df_oblong, aes(x=scale, y=oblong_taa)) + geom_bar(stat = "identity") + ggtitle("STOP, oblong, TAA")
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/stop/common_oblong_taa.png")
ggplot(df_oblong, aes(x=scale, y=oblong_tga)) + geom_bar(stat = "identity") + ggtitle("STOP, oblong, TGA")
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/stop/common_oblong_tga.png")
ggplot(df_oblong, aes(x=scale, y=oblong_tag)) + geom_bar(stat = "identity") + ggtitle("STOP, oblong, TAG")
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/stop/common_oblong_tag.png")
