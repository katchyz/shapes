### shape vs RNA
library(GenomicFeatures)
library(ggplot2)

## ribo-seq & RNA-seq
load(file="/Volumes/USELESS/META/SHAPES/FPKM_1K.Rdata")
load(file="/Volumes/USELESS/META/SHAPES/FPKM_24.Rdata")
ribo_fpkm <- data.frame(cell1K_ribo = FPKM_1K$exons_Ribo_fpkm, cell24_ribo = FPKM_24$exons_Ribo_fpkm, n = rownames(FPKM_1K))
RNA_fpkm <- data.frame(cell1K_rna = FPKM_1K$exons_RNA_fpkm, cell24_rna = FPKM_24$exons_RNA_fpkm, n = rownames(FPKM_1K))

## Shape-seq
# shape
shape_cell1K_et1 <- get(load(file = "/Volumes/USELESS/test2/normalized/cell1K_invivo_et1.Rsave"))
shape_cell1K_et2 <- get(load(file = "/Volumes/USELESS/test2/normalized/cell1K_invivo_et2.Rsave"))
shape_cell1_et1 <- get(load(file = "/Volumes/USELESS/test2/normalized/cell1_invivo_et1.Rsave"))
shape_cell1_et2 <- get(load(file = "/Volumes/USELESS/test2/normalized/cell1_invivo_et2.Rsave"))
rm(shapeseq_norm)


df_shape_cell1K_et1 <- data.frame(shape_cell1K_et1_TC.control = sapply(shape_cell1K_et1, function(x){mean(x$TC.control[2:length(x)])}),
                                  shape_cell1K_et1_TC.treated = sapply(shape_cell1K_et1, function(x){mean(x$TC.treated[2:length(x)])}),
                                  shape_cell1K_et1_log2ratio = sapply(shape_cell1K_et1, function(x){mean(x$log2ratio[2:length(x)])}),
                                  n = names(shape_cell1K_et1))

df_shape_cell1K_et2 <- data.frame(shape_cell1K_et2_TC.control = sapply(shape_cell1K_et2, function(x){mean(x$TC.control[2:length(x)])}),
                                  shape_cell1K_et2_TC.treated = sapply(shape_cell1K_et2, function(x){mean(x$TC.treated[2:length(x)])}),
                                  shape_cell1K_et2_log2ratio = sapply(shape_cell1K_et2, function(x){mean(x$log2ratio[2:length(x)])}),
                                  n = names(shape_cell1K_et2))

df_shape_cell1_et1 <- data.frame(shape_cell1_et1_TC.control = sapply(shape_cell1_et1, function(x){mean(x$TC.control[2:length(x)])}),
                                  shape_cell1_et1_TC.treated = sapply(shape_cell1_et1, function(x){mean(x$TC.treated[2:length(x)])}),
                                  shape_cell1_et1_log2ratio = sapply(shape_cell1_et1, function(x){mean(x$log2ratio[2:length(x)])}),
                                  n = names(shape_cell1_et1))

df_shape_cell1_et2 <- data.frame(shape_cell1_et2_TC.control = sapply(shape_cell1_et2, function(x){mean(x$TC.control[2:length(x)])}),
                                  shape_cell1_et2_TC.treated = sapply(shape_cell1_et2, function(x){mean(x$TC.treated[2:length(x)])}),
                                  shape_cell1_et2_log2ratio = sapply(shape_cell1_et2, function(x){mean(x$log2ratio[2:length(x)])}),
                                  n = names(shape_cell1_et2))



df_shape <- merge(df_shape_cell1K_et1, df_shape_cell1K_et2, by="n", all=TRUE)
df_shape <- merge(df_shape, df_shape_cell1_et1, by="n", all=TRUE)
df_shape <- merge(df_shape, df_shape_cell1_et2, by="n", all=TRUE)
df_shape <- merge(df_shape, ribo_fpkm, by="n", all=TRUE)
df_shape <- merge(df_shape, RNA_fpkm, by="n", all=TRUE)


p1 <- ggplot(df_shape, aes(x = shape_cell1_et1_TC.control, y = cell24_rna)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
p2 <- ggplot(df_shape, aes(x = shape_cell1_et1_TC.treated, y = cell24_rna)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
p3 <- ggplot(df_shape, aes(x = shape_cell1_et1_log2ratio, y = cell24_rna)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
p4 <- ggplot(df_shape, aes(x = shape_cell1_et2_TC.control, y = cell24_rna)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
p5 <- ggplot(df_shape, aes(x = shape_cell1_et2_TC.treated, y = cell24_rna)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
p6 <- ggplot(df_shape, aes(x = shape_cell1_et2_log2ratio, y = cell24_rna)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()

p7 <- ggplot(df_shape, aes(x = shape_cell1K_et1_TC.control, y = cell1K_rna)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
p8 <- ggplot(df_shape, aes(x = shape_cell1K_et1_TC.treated, y = cell1K_rna)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
p9 <- ggplot(df_shape, aes(x = shape_cell1K_et1_log2ratio, y = cell1K_rna)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
p10 <- ggplot(df_shape, aes(x = shape_cell1K_et2_TC.control, y = cell1K_rna)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
p11 <- ggplot(df_shape, aes(x = shape_cell1K_et2_TC.treated, y = cell1K_rna)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
p12 <- ggplot(df_shape, aes(x = shape_cell1K_et2_log2ratio, y = cell1K_rna)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()

ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_RNA/new/p1.png", plot = p1)
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_RNA/new/p2.png", plot = p2)
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_RNA/new/p3.png", plot = p3)
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_RNA/new/p4.png", plot = p4)
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_RNA/new/p5.png", plot = p5)
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_RNA/new/p6.png", plot = p6)

ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_RNA/new/p7.png", plot = p7)
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_RNA/new/p8.png", plot = p8)
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_RNA/new/p9.png", plot = p9)
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_RNA/new/p10.png", plot = p10)
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_RNA/new/p11.png", plot = p11)
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_RNA/new/p12.png", plot = p12)

#####
p1 <- ggplot(df_shape, aes(x = shape_cell1_et1_TC.control, y = cell24_ribo)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_ribo/new/p1.png", plot = p1)
p2 <- ggplot(df_shape, aes(x = shape_cell1_et1_TC.treated, y = cell24_ribo)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_ribo/new/p2.png", plot = p2)
p3 <- ggplot(df_shape, aes(x = shape_cell1_et1_log2ratio, y = cell24_ribo)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_ribo/new/p3.png", plot = p3)
p4 <- ggplot(df_shape, aes(x = shape_cell1_et2_TC.control, y = cell24_ribo)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_ribo/new/p4.png", plot = p4)
p5 <- ggplot(df_shape, aes(x = shape_cell1_et2_TC.treated, y = cell24_ribo)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_ribo/new/p5.png", plot = p5)
p6 <- ggplot(df_shape, aes(x = shape_cell1_et2_log2ratio, y = cell24_ribo)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_ribo/new/p6.png", plot = p6)

p7 <- ggplot(df_shape, aes(x = shape_cell1K_et1_TC.control, y = cell1K_ribo)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_ribo/new/p7.png", plot = p7)
p8 <- ggplot(df_shape, aes(x = shape_cell1K_et1_TC.treated, y = cell1K_ribo)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_ribo/new/p8.png", plot = p8)
p9 <- ggplot(df_shape, aes(x = shape_cell1K_et1_log2ratio, y = cell1K_ribo)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_ribo/new/p9.png", plot = p9)
p10 <- ggplot(df_shape, aes(x = shape_cell1K_et2_TC.control, y = cell1K_ribo)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_ribo/new/p10.png", plot = p10)
p11 <- ggplot(df_shape, aes(x = shape_cell1K_et2_TC.treated, y = cell1K_ribo)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_ribo/new/p11.png", plot = p11)
p12 <- ggplot(df_shape, aes(x = shape_cell1K_et2_log2ratio, y = cell1K_ribo)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_ribo/new/p12.png", plot = p12)

########
df_shape$te1 <- df_shape$cell24_ribo / df_shape$cell24_rna
df_shape$te1K <- df_shape$cell1K_ribo / df_shape$cell1K_rna

p1 <- ggplot(df_shape, aes(x = shape_cell1_et1_TC.control, y = te1)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_TE/new/p1.png", plot = p1)
p2 <- ggplot(df_shape, aes(x = shape_cell1_et1_TC.treated, y = te1)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_TE/new/p2.png", plot = p2)
p3 <- ggplot(df_shape, aes(x = shape_cell1_et1_log2ratio, y = te1)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_TE/new/p3.png", plot = p3)
p4 <- ggplot(df_shape, aes(x = shape_cell1_et2_TC.control, y = te1)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_TE/new/p4.png", plot = p4)
p5 <- ggplot(df_shape, aes(x = shape_cell1_et2_TC.treated, y = te1)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_TE/new/p5.png", plot = p5)
p6 <- ggplot(df_shape, aes(x = shape_cell1_et2_log2ratio, y = te1)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_TE/new/p6.png", plot = p6)

p7 <- ggplot(df_shape, aes(x = shape_cell1K_et1_TC.control, y = te1K)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_TE/new/p7.png", plot = p7)
p8 <- ggplot(df_shape, aes(x = shape_cell1K_et1_TC.treated, y = te1K)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_TE/new/p8.png", plot = p8)
p9 <- ggplot(df_shape, aes(x = shape_cell1K_et1_log2ratio, y = te1K)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_TE/new/p9.png", plot = p9)
p10 <- ggplot(df_shape, aes(x = shape_cell1K_et2_TC.control, y = te1K)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_TE/new/p10.png", plot = p10)
p11 <- ggplot(df_shape, aes(x = shape_cell1K_et2_TC.treated, y = te1K)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_TE/new/p11.png", plot = p11)
p12 <- ggplot(df_shape, aes(x = shape_cell1K_et2_log2ratio, y = te1K)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_TE/new/p12.png", plot = p12)

######################
#### in vitro ########

tc_cell1_invitro_et1 <- get(load(file = "/Volumes/USELESS/test2/normalized/RAW_cell1_invitro_NAI_et1.Rsave"))
tc_cell1_invitro_et1 <- split(tc_cell1_invitro_et1, seqnames(tc_cell1_invitro_et1))
tc_cell1_invitro_et1 <- tc_cell1_invitro_et1[sapply(tc_cell1_invitro_et1, function(x){length(x) > 0})]
tc_cell1_invitro_et2 <- get(load(file = "/Volumes/USELESS/test2/normalized/RAW_cell1_invitro_NAI_et2.Rsave"))
tc_cell1_invitro_et2 <- split(tc_cell1_invitro_et2, seqnames(tc_cell1_invitro_et2))
tc_cell1_invitro_et2 <- tc_cell1_invitro_et2[sapply(tc_cell1_invitro_et2, function(x){length(x) > 0})]
rm(treated_comp)

df_cell1_invitro_et1 <- data.frame(shape_cell1_invitro_et1_TC = sapply(tc_cell1_invitro_et1, function(x){mean(x$TC[2:length(x)])}),
                                  n = names(tc_cell1_invitro_et1))

df_cell1_invitro_et2 <- data.frame(shape_cell1_invitro_et2_TC = sapply(tc_cell1_invitro_et2, function(x){mean(x$TC[2:length(x)])}),
                                   n = names(tc_cell1_invitro_et2))

df_shape <- merge(df_shape, df_cell1_invitro_et1, by="n", all=TRUE)
df_shape <- merge(df_shape, df_cell1_invitro_et2, by="n", all=TRUE)

ggplot(df_shape, aes(x = shape_cell1_invitro_et1_TC, y = cell24_rna)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_RNA/new/invitro1.png")

ggplot(df_shape, aes(x = shape_cell1_invitro_et2_TC, y = cell24_rna)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_RNA/new/invitro2.png")

