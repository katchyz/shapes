### internal vs complete read-through
library(GenomicFeatures)
library(ggplot2)

## ribo-seq & RNA-seq
load(file="/Volumes/USELESS/META/SHAPES/FPKM_24.Rdata")
load(file="/Volumes/USELESS/META/SHAPES/FPKM_256.Rdata")
ribo_fpkm <- data.frame(cell24_ribo = FPKM_24$exons_Ribo_fpkm, cell256_ribo = FPKM_256$exons_Ribo_fpkm, n = rownames(FPKM_24))
RNA_fpkm <- data.frame(cell24_rna = FPKM_24$exons_RNA_fpkm, cell256_rna = FPKM_256$exons_RNA_fpkm, n = rownames(FPKM_24))

## Shape-seq
load(file = "/Volumes/USELESS/test/normalized/shape_cell24.Rsave")
load(file = "/Volumes/USELESS/test/normalized/shape_cell256.Rsave")
load(file = "/Volumes/USELESS/test/normalized/shape_oblong.Rsave")
load(file = "/Volumes/USELESS/test/normalized/shape_oblong_CHX.Rsave")


df_shape_cell24 <- data.frame(shape_cell24_TC.control = sapply(shape_cell24, function(x){mean(x$TC.control)}),
                              shape_cell24_TC.treated = sapply(shape_cell24, function(x){mean(x$TC.treated)}),
                              shape_cell24_log2ratio = sapply(shape_cell24, function(x){mean(x$log2ratio)}),
                              n = names(shape_cell24))

df_shape_cell256 <- data.frame(shape_cell256_TC.control = sapply(shape_cell256, function(x){mean(x$TC.control)}),
                               shape_cell256_TC.treated = sapply(shape_cell256, function(x){mean(x$TC.treated)}),
                               shape_cell256_log2ratio = sapply(shape_cell256, function(x){mean(x$log2ratio)}),
                               n = names(shape_cell256))

df_shape <- merge(df_shape_cell24, df_shape_cell256, by="n", all=TRUE)
df_shape <- merge(df_shape, ribo_fpkm, by="n", all=TRUE)
df_shape <- merge(df_shape, RNA_fpkm, by="n", all=TRUE)

p1 <- ggplot(df_shape, aes(x = shape_cell24_TC.control, y = cell24_rna)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
p2 <- ggplot(df_shape, aes(x = shape_cell24_TC.treated, y = cell24_rna)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
p3 <- ggplot(df_shape, aes(x = shape_cell24_log2ratio, y = cell24_rna)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
p4 <- ggplot(df_shape, aes(x = shape_cell256_TC.control, y = cell256_rna)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
p5 <- ggplot(df_shape, aes(x = shape_cell256_TC.treated, y = cell256_rna)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
p6 <- ggplot(df_shape, aes(x = shape_cell256_log2ratio, y = cell256_rna)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()

ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_RNA/p1.png", plot = p1)
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_RNA/p2.png", plot = p2)
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_RNA/p3.png", plot = p3)
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_RNA/p4.png", plot = p4)
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_RNA/p5.png", plot = p5)
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_RNA/p6.png", plot = p6)

#####
p1 <- ggplot(df_shape, aes(x = shape_cell24_TC.control, y = cell24_ribo)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_ribo/p1.png", plot = p1)
p2 <- ggplot(df_shape, aes(x = shape_cell24_TC.treated, y = cell24_ribo)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_ribo/p2.png", plot = p2)
p3 <- ggplot(df_shape, aes(x = shape_cell24_log2ratio, y = cell24_ribo)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_ribo/p3.png", plot = p3)

p4 <- ggplot(df_shape, aes(x = shape_cell256_TC.control, y = cell256_ribo)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_ribo/p4.png", plot = p4)
p5 <- ggplot(df_shape, aes(x = shape_cell256_TC.treated, y = cell256_ribo)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_ribo/p5.png", plot = p5)
p6 <- ggplot(df_shape, aes(x = shape_cell256_log2ratio, y = cell256_ribo)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_ribo/p6.png", plot = p6)


########
df_shape$te24 <- df_shape$cell24_ribo / df_shape$cell24_rna
df_shape$te256 <- df_shape$cell256_ribo / df_shape$cell256_rna

p1 <- ggplot(df_shape, aes(x = shape_cell24_TC.control, y = te24)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_TE/p1.png", plot = p1)
p2 <- ggplot(df_shape, aes(x = shape_cell24_TC.treated, y = te24)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_TE/p2.png", plot = p2)
p3 <- ggplot(df_shape, aes(x = shape_cell24_log2ratio, y = te24)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_TE/p3.png", plot = p3)

p4 <- ggplot(df_shape, aes(x = shape_cell256_TC.control, y = te256)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_TE/p4.png", plot = p4)
p5 <- ggplot(df_shape, aes(x = shape_cell256_TC.treated, y = te256)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_TE/p5.png", plot = p5)
p6 <- ggplot(df_shape, aes(x = shape_cell256_log2ratio, y = te256)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_TE/p6.png", plot = p6)


### 5'end peak (compare with Cage data?)
# Cage
fertilized <- read.csv(file = "/Volumes/USELESS/DATA/Cage_Nepal/leaders_zv10_zf_02_fertilized_egg.csv")
cells_512 <- read.csv(file = "/Volumes/USELESS/DATA/Cage_Nepal/leaders_zv10_zf_04_cells_512.csv")
# set threshold of 10 normalized counts
fertilized <- fertilized[fertilized$count_at_highest_peak > 10,]
cells_512 <- cells_512[cells_512$count_at_highest_peak > 10,]

txdb_can <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/danio_rerio/Danio_rerio.GRCz10.84.gtf", format = "gtf")
exons <- exonsBy(txdb_can, by="tx", use.names=TRUE)

cage_fertilized <- GRanges(Rle(gsub('chr', '', as.character(fertilized$chr)), rep(1, nrow(fertilized))), IRanges(fertilized$highest_peak, width=rep(1, nrow(fertilized))), strand = fertilized$dir)
cage_cells_512 <- GRanges(Rle(gsub('chr', '', as.character(cells_512$chr)), rep(1, nrow(cells_512))), IRanges(cells_512$highest_peak, width=rep(1, nrow(cells_512))), strand = cells_512$dir)

load(file = "/Volumes/USELESS/test/normalized/gc_2-4cell.Rsave")
load(file = "/Volumes/USELESS/test/normalized/gc_256cell.Rsave")

peaks24 <- subsetByOverlaps(gc_shape_cell24, cage_fertilized, ignore.strand=TRUE)
peaks256 <- subsetByOverlaps(gc_shape_cell256, cage_fertilized, ignore.strand=TRUE)

p24 <- data.frame(n = substr(peaks24@ranges@NAMES, 1, 18), TC.control = peaks24$TC.control, TC.treated = peaks24$TC.treated, log2ratio = peaks24$log2ratio)
p256 <- data.frame(n = substr(peaks256@ranges@NAMES, 1, 18), TC.control = peaks256$TC.control, TC.treated = peaks256$TC.treated, log2ratio = peaks256$log2ratio)

p24 <- merge(p24, RNA_fpkm, by="n", all=TRUE)
p256 <- merge(p256, RNA_fpkm, by="n", all=TRUE)

ggplot(p24, aes(x = TC.control, y = cell24_rna)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10() + ggtitle("cage, cell24")
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_RNA/cage_5end/p1.png")
ggplot(p24, aes(x = TC.treated, y = cell24_rna)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10() + ggtitle("cage, cell24")
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_RNA/cage_5end/p2.png")
ggplot(p24, aes(x = log2ratio, y = cell24_rna)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10() + ggtitle("cage, cell24")
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_RNA/cage_5end/p3.png")

ggplot(p256, aes(x = TC.control, y = cell256_rna)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10() + ggtitle("cage, cell256")
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_RNA/cage_5end/p4.png")
ggplot(p256, aes(x = TC.treated, y = cell256_rna)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10() + ggtitle("cage, cell256")
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_RNA/cage_5end/p5.png")
ggplot(p256, aes(x = log2ratio, y = cell256_rna)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10() + ggtitle("cage, cell256")
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_RNA/cage_5end/p6.png")


##########
df_shape_5end_cell24 <- data.frame(shape_cell24_TC.control = sapply(shape_cell24, function(x){x$TC.control[1]/sum(x$TC.control)}),
                                   shape_cell24_TC.treated = sapply(shape_cell24, function(x){x$TC.treated[1]/sum(x$TC.treated)}),
                                   shape_cell24_log2ratio = sapply(shape_cell24, function(x){x$log2ratio[1]/sum(x$log2ratio)}),
                                   n = names(shape_cell24))

df_shape_5end_cell256 <- data.frame(shape_cell256_TC.control = sapply(shape_cell256, function(x){x$TC.control[1]/sum(x$TC.control)}),
                                    shape_cell256_TC.treated = sapply(shape_cell256, function(x){x$TC.treated[1]/sum(x$TC.treated)}),
                                    shape_cell256_log2ratio = sapply(shape_cell256, function(x){x$log2ratio[1]/sum(x$log2ratio)}),
                                    n = names(shape_cell256))


df_5end <- merge(df_shape_5end_cell24, df_shape_5end_cell256, by="n", all=TRUE)
df_5end <- merge(df_5end, ribo_fpkm, by="n", all=TRUE)
df_5end <- merge(df_5end, RNA_fpkm, by="n", all=TRUE)

ggplot(df_5end, aes(x = shape_cell24_TC.control, y = cell24_rna)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_RNA/annotated_5end/p1.png")
ggplot(df_5end, aes(x = shape_cell24_TC.treated, y = cell24_rna)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_RNA/annotated_5end/p2.png")
ggplot(df_5end, aes(x = shape_cell24_log2ratio, y = cell24_rna)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_RNA/annotated_5end/p3.png")

ggplot(df_5end, aes(x = shape_cell256_TC.control, y = cell256_rna)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_RNA/annotated_5end/p4.png")
ggplot(df_5end, aes(x = shape_cell256_TC.treated, y = cell256_rna)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_RNA/annotated_5end/p5.png")
ggplot(df_5end, aes(x = shape_cell256_log2ratio, y = cell256_rna)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_RNA/annotated_5end/p6.png")
