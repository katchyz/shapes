### distributions (utr5, CDS, utr3)

library(data.table)
library(ggplot2)
library(GenomicRanges)
library(GenomicFeatures)
library(plyr)

### get transcript database
txdb_can <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/danio_rerio/Danio_rerio.GRCz10.84.gtf", format = "gtf")
txLengths <- transcriptLengths(txdb_can, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
rownames(txLengths) <- txLengths$tx_name
txLengths <- txLengths[order(rownames(txLengths)),]

# get longest transcript per gene
txlen <- arrange(txLengths, gene_id, desc(tx_len))
txlen <- txlen[!duplicated(txlen$gene_id),]

txlen <- txlen[(txlen$utr5_len > 10) & (txlen$cds_len > 10) & (txlen$utr3_len > 10), ]

### load data
load(file="/Volumes/USELESS/META/SHAPES_NEW/general/clouds/cloud_df.Rsave")
# shape
load(file="/Volumes/USELESS/DATA/Shape-Seq/shape_cell1.Rsave")
load(file="/Volumes/USELESS/DATA/Shape-Seq/shape_cell24.Rsave")
load(file="/Volumes/USELESS/DATA/Shape-Seq/shape_cell256.Rsave")
load(file="/Volumes/USELESS/DATA/Shape-Seq/shape_cell1K.Rsave")
load(file="/Volumes/USELESS/DATA/Shape-Seq/shape_oblong.Rsave")
load(file="/Volumes/USELESS/DATA/Shape-Seq/shape_oblong_CHX.Rsave")

# shapes
shapes_cell1_invitro <- get(load(file = "/Volumes/USELESS/test/normalized/cell1_invitro_NAI_N3.Rsave"))
shapes_cell24 <- get(load(file = "/Volumes/USELESS/test/normalized/cell2_4_invivo_NAI_N3.Rsave"))
rm(shapes_norm)

#shape <- shape_cell1
shape <- shape_oblong_CHX
len <- txlen[txlen$tx_name %in% names(shape),]
len <- len[order(len$tx_name),]
shape <- shape[names(shape) %in% txlen$tx_name]
shape <- shape[order(names(shape))]
leader <- mapply(function(x,y){mean(x$log2ratio[2:y])}, shape, len$utr5_len)
shape <- mapply(function(x,y){x$log2ratio[(y+1):length(x)]}, shape, len$utr5_len)
cds <- mapply(function(x,y){mean(x[1:y])}, shape, len$cds_len)
trailer <- mapply(function(x,y){mean(x[(y+1):length(x)])}, shape, len$cds_len)


d <- data.frame(leader, cds, trailer)
d$scale <- c(1:nrow(d))
d <- melt(d, id.var="scale")
#ggplot(d, aes(x=value, fill=variable)) + geom_density(alpha=.3) + ggtitle("cell1") + scale_x_log10()
#ggsave(file="/Volumes/USELESS/META/SHAPES_NEW/distributions/cell1.png")
ggplot(d, aes(x=value, fill=variable)) + geom_density(alpha=.3) + ggtitle("oblong_CHX") + scale_x_log10()
ggsave(file="/Volumes/USELESS/META/SHAPES_NEW/distributions/oblong_CHX.png")

### NAI-N3
shape <- shapes_cell24
len <- txlen[txlen$tx_name %in% names(shape),]
len <- len[order(len$tx_name),]
shape <- shape[names(shape) %in% txlen$tx_name]
shape <- shape[order(names(shape))]
leader <- mapply(function(x,y){mean(x$TC[2:y])}, shape, len$utr5_len)
shape <- mapply(function(x,y){x$TC[(y+1):length(x)]}, shape, len$utr5_len)
cds <- mapply(function(x,y){mean(x[1:y])}, shape, len$cds_len)
trailer <- mapply(function(x,y){mean(x[(y+1):length(x)])}, shape, len$cds_len)


d <- data.frame(leader, cds, trailer)
d$scale <- c(1:nrow(d))
d <- melt(d, id.var="scale")
#ggplot(d, aes(x=value, fill=variable)) + geom_density(alpha=.3) + ggtitle("cell1") + scale_x_log10()
#ggsave(file="/Volumes/USELESS/META/SHAPES_NEW/distributions/cell1.png")
ggplot(d, aes(x=value, fill=variable)) + geom_density(alpha=.3) + ggtitle("cell24_TC") + scale_x_log10()
ggsave(file="/Volumes/USELESS/META/SHAPES_NEW/distributions/cell24_TC.png")

##### incl. 5'end
shape <- shape_cell256
len <- txlen[txlen$tx_name %in% names(shape),]
len <- len[order(len$tx_name),]
shape <- shape[names(shape) %in% txlen$tx_name]
shape <- shape[order(names(shape))]
leader <- mapply(function(x,y){mean(x$log2ratio[1:y])}, shape, len$utr5_len)
shape <- mapply(function(x,y){x$log2ratio[(y+1):length(x)]}, shape, len$utr5_len)
cds <- mapply(function(x,y){mean(x[1:y])}, shape, len$cds_len)
trailer <- mapply(function(x,y){mean(x[(y+1):length(x)])}, shape, len$cds_len)

d <- data.frame(leader, cds, trailer)
d$scale <- c(1:nrow(d))
d <- melt(d, id.var="scale")
ggplot(d, aes(x=value, fill=variable)) + geom_density(alpha=.3) + ggtitle("cell256 incl. 5'end") + scale_x_log10()


################################
#### Ribo-Seq
load(file="/Volumes/USELESS/DATA/Ribo-Seq/ribo24.Rsave")
load(file="/Volumes/USELESS/DATA/Ribo-Seq/ribo256.Rsave")
load(file="/Volumes/USELESS/DATA/Ribo-Seq/ribo1K.Rsave")


#########
ribo <- ribo1K
len <- txlen[txlen$tx_name %in% names(ribo),]
len <- len[order(len$tx_name),]
ribo <- ribo[names(ribo) %in% txlen$tx_name]
ribo <- ribo[order(names(ribo))]
leader <- mapply(function(x,y){mean(x[1:y])}, ribo, len$utr5_len)
ribo <- mapply(function(x,y){x[(y+1):length(x)]}, ribo, len$utr5_len)
cds <- mapply(function(x,y){mean(x[1:y])}, ribo, len$cds_len)
trailer <- mapply(function(x,y){mean(x[(y+1):length(x)])}, ribo, len$cds_len)


d <- data.frame(leader, cds, trailer)
d$scale <- c(1:nrow(d))
d <- melt(d, id.var="scale")
ggplot(d, aes(x=value, fill=variable)) + geom_density(alpha=.3) + ggtitle("ribo 1K") + scale_x_log10()
ggsave(file="/Volumes/USELESS/META/SHAPES_NEW/distributions/ribo/ribo1K.png")

