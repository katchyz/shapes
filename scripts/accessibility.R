### accessibility
library(GenomicFeatures)
library(ggplot2)

## Shape-seq
load(file = "/Volumes/USELESS/test/normalized/shape_cell24.Rsave")
load(file = "/Volumes/USELESS/test/normalized/shape_cell256.Rsave")
load(file = "/Volumes/USELESS/test/normalized/shape_oblong.Rsave")
load(file = "/Volumes/USELESS/test/normalized/shape_oblong_CHX.Rsave")

## ribo-seq & RNA-seq
load(file="/Volumes/USELESS/META/SHAPES/FPKM_24.Rdata")
load(file="/Volumes/USELESS/META/SHAPES/FPKM_256.Rdata")
ribo_fpkm <- data.frame(cell24_ribo = FPKM_24$exons_Ribo_fpkm, cell256_ribo = FPKM_256$exons_Ribo_fpkm, n = rownames(FPKM_24))
RNA_fpkm <- data.frame(cell24_rna = FPKM_24$exons_RNA_fpkm, cell256_rna = FPKM_256$exons_RNA_fpkm, n = rownames(FPKM_24))


rna <- RNA_fpkm[rownames(RNA_fpkm) %in% names(shape_cell24),]
shape <- shape_cell24[names(shape_cell24) %in% rownames(rna)]
rna <- rna[order(rownames(rna)),]
shape <- shape[order(names(shape))]
accessibility_cell24 <- mapply(function(x,y){x$accessibility <- x$log2ratio/y}, shape, rna$cell24_rna)
#accessibility_cell24 <- accessibility_cell24[!sapply(accessibility_cell24, function(x){any(is.na(x))})]
#accessibility_cell24 <- accessibility_cell24[sapply(accessibility_cell24, function(x){sum(x) > 0})]
#save(accessibility_cell24, file="/Users/kasia/Desktop/accessibility_cell24_dla_Kornelka.Rsave")

## add to GRanges
a24 <- unlist(accessibility_cell24)
s24 <- unlist(shape)
s24$accessibility <- a24
s24 <- split(s24, seqnames(s24))
s24 <- s24[sapply(s24, function(x){length(x) > 0})]

rna <- RNA_fpkm[rownames(RNA_fpkm) %in% names(shape_cell256),]
shape <- shape_cell256[names(shape_cell256) %in% rownames(rna)]
rna <- rna[order(rownames(rna)),]
shape <- shape[order(names(shape))]
accessibility_cell256 <- mapply(function(x,y){x$accessibility <- x$log2ratio/y}, shape, rna$cell256_rna)
#accessibility_cell256 <- accessibility_cell256[!sapply(accessibility_cell256, function(x){any(is.na(x))})]
#accessibility_cell256 <- accessibility_cell256[sapply(accessibility_cell256, function(x){sum(x) > 0})]
#save(accessibility_cell256, file="/Users/kasia/Desktop/chuje_muje_dzikie_weze.Rsave")

## add to GRanges
a256 <- unlist(accessibility_cell256)
s256 <- unlist(shape)
s256$accessibility <- a256
s256 <- split(s256, seqnames(s256))
s256 <- s256[sapply(s256, function(x){length(x) > 0})]

### mean accessibility vs TE
acc24 <- sapply(s24, function(x){mean(x$accessibility)})
acc256 <- sapply(s256, function(x){mean(x$accessibility)})

acc24 <- data.frame(n = names(acc24), acc24 = as.numeric(acc24))
acc256 <- data.frame(n = names(acc256), acc256 = as.numeric(acc256))

df <- merge(acc24, acc256, by="n", all=TRUE)
df <- merge(df, RNA_fpkm, by="n", all=TRUE)
df <- merge(df, ribo_fpkm, by="n", all=TRUE)
df$te24 <- df$cell24_ribo / df$cell24_rna
df$te256 <- df$cell256_ribo / df$cell256_rna

ggplot(df, aes(x = acc24, y = te24)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/accessibility/acc24_te24.png")
ggplot(df, aes(x = acc256, y = te256)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/accessibility/acc256_te256.png")

### accessibility over start vs TE
