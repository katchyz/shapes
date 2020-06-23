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
scale <- c(-60:60)[-61]

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




###### libs
shape_libs <- c("shape_cell1", "shape_cell24", "shape_cell256", "shape_cell1K", "shape_oblong", "shape_oblong_CHX")
#shapes_libs <- c("shapes_cell1_invitro", "shapes_cell1_invitro_nonsel", "shapes_cell24")


###############
############# new data
for (i in 1:8) {
  
  u <- unlist(eval(parse(text = shape_libs[i])))
  u$log2ratio[u$log2ratio < 0] <- 0
  meta <- subsetByOverlaps(u, stop, ignore.strand=TRUE)
  meta <- split(meta, seqnames(meta))
  meta <- meta[sapply(meta, function(x){length(x) > 0})]
  meta <- meta[sapply(meta, function(x){sum(x$log2ratio) > 0})]
  
  assign(paste0("var",i), names(meta))
  sub_taa <- meta[names(meta) %in% taa]
  sub_tga <- meta[names(meta) %in% tga]
  sub_tag <- meta[names(meta) %in% tag]
  d <- data.frame(scale = scale,
                  taa = Reduce("+", lapply(sub_taa, function(x){x$log2ratio/sum(x$log2ratio)})),
                  tga = Reduce("+", lapply(sub_tga, function(x){x$log2ratio/sum(x$log2ratio)})),
                  tag = Reduce("+", lapply(sub_tag, function(x){x$log2ratio/sum(x$log2ratio)})))
  
  ggplot(d, aes(x = scale, y = taa)) + geom_bar(stat = "identity") + ggtitle(paste(shape_libs[i], "STOP, TAA"))
  fp <- paste0("/Volumes/USELESS/META/SHAPES_NEW/stop/p_taa", i, ".png", collapse = NULL)
  ggsave(fp)
  ggplot(d, aes(x = scale, y = tga)) + geom_bar(stat = "identity") + ggtitle(paste(shape_libs[i], "STOP, TGA"))
  fp <- paste0("/Volumes/USELESS/META/SHAPES_NEW/stop/p_tga", i, ".png", collapse = NULL)
  ggsave(fp)
  ggplot(d, aes(x = scale, y = tag)) + geom_bar(stat = "identity") + ggtitle(paste(shape_libs[i], "STOP, TAG"))
  fp <- paste0("/Volumes/USELESS/META/SHAPES_NEW/stop/p_tag", i, ".png", collapse = NULL)
  ggsave(fp)
  
}


for (i in 1:3) {
  
  u <- unlist(eval(parse(text = shapes_libs[i])))
  meta <- subsetByOverlaps(u, stop, ignore.strand=TRUE)
  meta <- split(meta, seqnames(meta))
  meta <- meta[sapply(meta, function(x){length(x) > 0})]
  meta <- meta[sapply(meta, function(x){sum(x$TC) > 0})]
  
  assign(paste0("svar",i), names(meta))
  sub_taa <- meta[names(meta) %in% taa]
  sub_tga <- meta[names(meta) %in% tga]
  sub_tag <- meta[names(meta) %in% tag]
  d <- data.frame(scale = scale,
                  taa = Reduce("+", lapply(sub_taa, function(x){x$TC/sum(x$TC)})),
                  tga = Reduce("+", lapply(sub_tga, function(x){x$TC/sum(x$TC)})),
                  tag = Reduce("+", lapply(sub_tag, function(x){x$TC/sum(x$TC)})))
  
  ggplot(d, aes(x = scale, y = taa)) + geom_bar(stat = "identity") + ggtitle(paste(shapes_libs[i], "STOP, TAA"))
  fp <- paste0("/Volumes/USELESS/META/SHAPES_NEW/stop/s_taa", i, ".png", collapse = NULL)
  ggsave(fp)
  ggplot(d, aes(x = scale, y = tga)) + geom_bar(stat = "identity") + ggtitle(paste(shapes_libs[i], "STOP, TGA"))
  fp <- paste0("/Volumes/USELESS/META/SHAPES_NEW/stop/s_tga", i, ".png", collapse = NULL)
  ggsave(fp)
  ggplot(d, aes(x = scale, y = tag)) + geom_bar(stat = "identity") + ggtitle(paste(shapes_libs[i], "STOP, TAG"))
  fp <- paste0("/Volumes/USELESS/META/SHAPES_NEW/stop/s_tag", i, ".png", collapse = NULL)
  ggsave(fp)
  
}


### convert p_taa* +append taa.png
### convert p_tga* +append tga.png
### convert p_tag* +append tag.png
### convert taa.png tga.png tag.png -append stop_shape.png

### convert s_taa* +append staa.png
### convert s_tga* +append stga.png
### convert s_tag* +append stag.png
### convert staa.png stga.png stag.png -append stop_shapes.png


### get tx in common
common <- var1[var1 %in% var2]
common <- common[common %in% var3]
common <- common[common %in% var4]
common <- common[common %in% var5]
common <- common[common %in% var6]
common <- common[common %in% var7]
common <- common[common %in% var8]
### 281 genes 

commons <- svar1[svar1 %in% svar2]
commons <- commons[commons %in% svar3]
### 1059 genes


for (i in 1:8) {
  
  u <- unlist(eval(parse(text = shape_libs[i])))
  u$log2ratio[u$log2ratio < 0] <- 0
  meta <- subsetByOverlaps(u, stop, ignore.strand=TRUE)
  meta <- split(meta, seqnames(meta))
  meta <- meta[names(meta) %in% common]
  meta <- meta[sapply(meta, function(x){sum(x$log2ratio) > 0})]
  
  sub_taa <- meta[names(meta) %in% taa]
  sub_tga <- meta[names(meta) %in% tga]
  sub_tag <- meta[names(meta) %in% tag]
  d <- data.frame(scale = scale,
                  taa = Reduce("+", lapply(sub_taa, function(x){x$log2ratio/sum(x$log2ratio)})),
                  tga = Reduce("+", lapply(sub_tga, function(x){x$log2ratio/sum(x$log2ratio)})),
                  tag = Reduce("+", lapply(sub_tag, function(x){x$log2ratio/sum(x$log2ratio)})))
  
  ggplot(d, aes(x = scale, y = taa)) + geom_bar(stat = "identity") + ggtitle(paste(shape_libs[i], "STOP, TAA, common"))
  fp <- paste0("/Volumes/USELESS/META/SHAPES_NEW/stop/common_p_taa", i, ".png", collapse = NULL)
  ggsave(fp)
  ggplot(d, aes(x = scale, y = tga)) + geom_bar(stat = "identity") + ggtitle(paste(shape_libs[i], "STOP, TGA, common"))
  fp <- paste0("/Volumes/USELESS/META/SHAPES_NEW/stop/common_p_tga", i, ".png", collapse = NULL)
  ggsave(fp)
  ggplot(d, aes(x = scale, y = tag)) + geom_bar(stat = "identity") + ggtitle(paste(shape_libs[i], "STOP, TAG, common"))
  fp <- paste0("/Volumes/USELESS/META/SHAPES_NEW/stop/common_p_tag", i, ".png", collapse = NULL)
  ggsave(fp)
  
}

for (i in 1:3) {
  
  u <- unlist(eval(parse(text = shapes_libs[i])))
  meta <- subsetByOverlaps(u, stop, ignore.strand=TRUE)
  meta <- split(meta, seqnames(meta))
  meta <- meta[names(meta) %in% commons]
  meta <- meta[sapply(meta, function(x){sum(x$TC) > 0})]
  
  sub_taa <- meta[names(meta) %in% taa]
  sub_tga <- meta[names(meta) %in% tga]
  sub_tag <- meta[names(meta) %in% tag]
  d <- data.frame(scale = scale,
                  taa = Reduce("+", lapply(sub_taa, function(x){x$TC/sum(x$TC)})),
                  tga = Reduce("+", lapply(sub_tga, function(x){x$TC/sum(x$TC)})),
                  tag = Reduce("+", lapply(sub_tag, function(x){x$TC/sum(x$TC)})))
  
  ggplot(d, aes(x = scale, y = taa)) + geom_bar(stat = "identity") + ggtitle(paste(shapes_libs[i], "STOP, TAA, common"))
  fp <- paste0("/Volumes/USELESS/META/SHAPES_NEW/stop/common_s_taa", i, ".png", collapse = NULL)
  ggsave(fp)
  ggplot(d, aes(x = scale, y = tga)) + geom_bar(stat = "identity") + ggtitle(paste(shapes_libs[i], "STOP, TGA, common"))
  fp <- paste0("/Volumes/USELESS/META/SHAPES_NEW/stop/common_s_tga", i, ".png", collapse = NULL)
  ggsave(fp)
  ggplot(d, aes(x = scale, y = tag)) + geom_bar(stat = "identity") + ggtitle(paste(shapes_libs[i], "STOP, TAG, common"))
  fp <- paste0("/Volumes/USELESS/META/SHAPES_NEW/stop/common_s_tag", i, ".png", collapse = NULL)
  ggsave(fp)
  
}


# convert common_p_taa* +append common_taa.png
# convert common_p_tga* +append common_tga.png
# convert common_p_tag* +append common_tag.png
# convert common_taa.png common_tga.png common_tag.png -append stop_shape_common.png

# convert common_s_taa* +append commons_taa.png
# convert common_s_tga* +append commons_tga.png
# convert common_s_tag* +append commons_tag.png
# convert commons_taa.png commons_tga.png commons_tag.png -append stop_shapes_common.png




###### control and treated

for (i in 1:6) {
  
  u <- unlist(eval(parse(text = shape_libs[i])))
  meta <- subsetByOverlaps(u, stop, ignore.strand=TRUE)
  meta <- split(meta, seqnames(meta))
  meta <- meta[sapply(meta, function(x){length(x) > 0})]
  meta <- meta[sapply(meta, function(x){sum(x$TC.control) > 0})]
  
  assign(paste0("var",i), names(meta))
  sub_taa <- meta[names(meta) %in% taa]
  sub_tga <- meta[names(meta) %in% tga]
  sub_tag <- meta[names(meta) %in% tag]
  d <- data.frame(scale = scale,
                  taa = Reduce("+", lapply(sub_taa, function(x){x$TC.control/sum(x$TC.control)})),
                  tga = Reduce("+", lapply(sub_tga, function(x){x$TC.control/sum(x$TC.control)})),
                  tag = Reduce("+", lapply(sub_tag, function(x){x$TC.control/sum(x$TC.control)})))
  
  ggplot(d, aes(x = scale, y = taa)) + geom_bar(stat = "identity") + ggtitle(paste(shape_libs[i], "STOP, TAA, ctr"))
  fp <- paste0("/Volumes/USELESS/META/SHAPES_NEW/stop/long/control/p_taa", i, ".png", collapse = NULL)
  ggsave(fp)
  ggplot(d, aes(x = scale, y = tga)) + geom_bar(stat = "identity") + ggtitle(paste(shape_libs[i], "STOP, TGA, ctr"))
  fp <- paste0("/Volumes/USELESS/META/SHAPES_NEW/stop/long/control/p_tga", i, ".png", collapse = NULL)
  ggsave(fp)
  ggplot(d, aes(x = scale, y = tag)) + geom_bar(stat = "identity") + ggtitle(paste(shape_libs[i], "STOP, TAG, ctr"))
  fp <- paste0("/Volumes/USELESS/META/SHAPES_NEW/stop/long/control/p_tag", i, ".png", collapse = NULL)
  ggsave(fp)
  
}




for (i in 1:6) {
  
  u <- unlist(eval(parse(text = shape_libs[i])))
  meta <- subsetByOverlaps(u, stop, ignore.strand=TRUE)
  meta <- split(meta, seqnames(meta))
  meta <- meta[sapply(meta, function(x){length(x) > 0})]
  meta <- meta[sapply(meta, function(x){sum(x$TC.treated) > 0})]
  
  assign(paste0("var",i), names(meta))
  sub_taa <- meta[names(meta) %in% taa]
  sub_tga <- meta[names(meta) %in% tga]
  sub_tag <- meta[names(meta) %in% tag]
  d <- data.frame(scale = scale,
                  taa = Reduce("+", lapply(sub_taa, function(x){x$TC.treated/sum(x$TC.treated)})),
                  tga = Reduce("+", lapply(sub_tga, function(x){x$TC.treated/sum(x$TC.treated)})),
                  tag = Reduce("+", lapply(sub_tag, function(x){x$TC.treated/sum(x$TC.treated)})))
  
  ggplot(d, aes(x = scale, y = taa)) + geom_bar(stat = "identity") + ggtitle(paste(shape_libs[i], "STOP, TAA, tre"))
  fp <- paste0("/Volumes/USELESS/META/SHAPES_NEW/stop/long/treated/p_taa", i, ".png", collapse = NULL)
  ggsave(fp)
  ggplot(d, aes(x = scale, y = tga)) + geom_bar(stat = "identity") + ggtitle(paste(shape_libs[i], "STOP, TGA, tre"))
  fp <- paste0("/Volumes/USELESS/META/SHAPES_NEW/stop/long/treated/p_tga", i, ".png", collapse = NULL)
  ggsave(fp)
  ggplot(d, aes(x = scale, y = tag)) + geom_bar(stat = "identity") + ggtitle(paste(shape_libs[i], "STOP, TAG, tre"))
  fp <- paste0("/Volumes/USELESS/META/SHAPES_NEW/stop/long/treated/p_tag", i, ".png", collapse = NULL)
  ggsave(fp)
  
}

