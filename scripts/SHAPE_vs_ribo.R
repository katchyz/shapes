### mean shape vs Ribo-seq
library(ggplot2)

## ribo-seq
load(file="/Volumes/USELESS/META/SHAPES/FPKM_24.Rdata")
load(file="/Volumes/USELESS/META/SHAPES/FPKM_256.Rdata")
ribo_fpkm <- data.frame(cell24_ribo = FPKM_24$exons_Ribo_fpkm, cell256_ribo = FPKM_256$exons_Ribo_fpkm, n = rownames(FPKM_24))
RNA_fpkm <- data.frame(cell24_rna = FPKM_24$exons_RNA_fpkm, cell256_rna = FPKM_256$exons_RNA_fpkm, n = rownames(FPKM_24))

# shape
shape_cell24 <- get(load(file = "/Volumes/USELESS/test/normalized/2-4cell.Rsave"))
shape_cell256 <- get(load(file = "/Volumes/USELESS/test/normalized/256cell.Rsave"))
rm(shapeseq_norm)

df_shape_cell24 <- data.frame(shape_cell24_TC.control = sapply(shape_cell24, function(x){mean(x$TC.control)}),
                              shape_cell24_TC.treated = sapply(shape_cell24, function(x){mean(x$TC.treated)}),
                              shape_cell24_log2ratio = sapply(shape_cell24, function(x){mean(x$log2ratio)}),
                              n = names(shape_cell24))

df_shape_cell256 <- data.frame(shape_cell256_TC.control = sapply(shape_cell256, function(x){mean(x$TC.control)}),
                               shape_cell256_TC.treated = sapply(shape_cell256, function(x){mean(x$TC.treated)}),
                               shape_cell256_log2ratio = sapply(shape_cell256, function(x){mean(x$log2ratio)}),
                               n = names(shape_cell256))

df <- merge(df_shape_cell24, df_shape_cell256, by="n", all=TRUE)
df <- merge(df, ribo_fpkm, by="n", all=TRUE)
df <- merge(df, RNA_fpkm, by="n", all=TRUE)

#### plot
p1 <- ggplot(df, aes(x = shape_cell24_TC.control, y = cell24_ribo)) + geom_point() + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_ribo/p1.png", plot = p1)
p2 <- ggplot(df, aes(x = shape_cell24_TC.treated, y = cell24_ribo)) + geom_point() + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_ribo/p2.png", plot = p2)
p3 <- ggplot(df, aes(x = shape_cell24_log2ratio, y = cell24_ribo)) + geom_point() + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_ribo/p3.png", plot = p3)

p4 <- ggplot(df, aes(x = shape_cell256_TC.control, y = cell256_ribo)) + geom_point() + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_ribo/p4.png", plot = p4)
p5 <- ggplot(df, aes(x = shape_cell256_TC.treated, y = cell256_ribo)) + geom_point() + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_ribo/p5.png", plot = p5)
p6 <- ggplot(df, aes(x = shape_cell256_log2ratio, y = cell256_ribo)) + geom_point() + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_ribo/p6.png", plot = p6)

#### what is this 'piggyback' cloud? - shape_cell24_TC.treated ### RIBO ###
f <- function(x) {x + 13}
ggplot(df, aes(x = log2(shape_cell24_TC.treated), y = log2(cell24_ribo))) + geom_point() + stat_function(fun = f, colour = "red")

piggyback <-df$n[log2(df$cell24_ribo) > (log2(df$shape_cell24_TC.treated) + 13)]
piggyback <- piggyback[!is.na(piggyback)]
piggyback <- data.frame(n = piggyback, piggyback = rep("piggyback", length(piggyback)))
df <- merge(df, piggyback, by="n", all=TRUE)


ggplot(df, aes(x = log2(shape_cell24_TC.treated), y = log2(cell24_ribo))) + geom_point() + stat_function(fun = f, colour = "red") +
  geom_point(data = subset(df, piggyback == 'piggyback'), aes(x=log2(shape_cell24_TC.treated), y=log2(cell24_ribo), colour=piggyback)) + theme(legend.position="none")
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_ribo/piggyback_NAI.png")

piggyback_NAI_ribo <- as.character(df[!is.na(df$piggyback),]$n)

#### GO
go <- read.csv(file="/Users/kasia/Documents/PhD/matrix_go_fq.csv")
go <- data.frame(n = go$X.transcript_id, go_term = go$GO_Term_Name)
df <- merge(df, go, by="n", all=TRUE)
go_piggyback_NAI_ribo <- subset(df, piggyback == "piggyback")$go_term

########
#### what is this 'piggyback' cloud? - shape_cell24_TC.treated ### RNA ###
f <- function(x) {x + 10}
ggplot(df, aes(x = log2(shape_cell24_TC.treated), y = log2(cell24_rna))) + geom_point() + stat_function(fun = f, colour = "red")

piggyback <-df$n[log2(df$cell24_rna) > (log2(df$shape_cell24_TC.treated) + 10)]
piggyback <- piggyback[!is.na(piggyback)]
piggyback <- data.frame(n = piggyback, piggyback = rep("piggyback", length(piggyback)))
df <- merge(df, piggyback, by="n", all=TRUE)


ggplot(df, aes(x = log2(shape_cell24_TC.treated), y = log2(cell24_rna))) + geom_point() + stat_function(fun = f, colour = "red") +
  geom_point(data = subset(df, piggyback.y == 'piggyback'), aes(x=log2(shape_cell24_TC.treated), y=log2(cell24_rna), colour=piggyback.y)) + theme(legend.position="none")

piggyback_NAI_rna <- as.character(df[!is.na(df$piggyback.y),]$n)
go_piggyback_NAI_rna <- subset(df, piggyback.y == "piggyback")$go_term

# > sum(piggyback_NAI_ribo %in% piggyback_NAI_rna)
# [1] 2224

#########################################################
############# shape_vs_TE ###############################
#########################################################

df$te24 <- df$cell24_ribo / df$cell24_rna
df$te256 <- df$cell256_ribo / df$cell256_rna

p1 <- ggplot(df, aes(x = shape_cell24_TC.control, y = te24)) + geom_point() + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_TE/p1.png", plot = p1)
p2 <- ggplot(df, aes(x = shape_cell24_TC.treated, y = te24)) + geom_point() + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_TE/p2.png", plot = p2)
p3 <- ggplot(df, aes(x = shape_cell24_log2ratio, y = te24)) + geom_point() + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_TE/p3.png", plot = p3)

p4 <- ggplot(df, aes(x = shape_cell256_TC.control, y = te256)) + geom_point() + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_TE/p4.png", plot = p4)
p5 <- ggplot(df, aes(x = shape_cell256_TC.treated, y = te256)) + geom_point() + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_TE/p5.png", plot = p5)
p6 <- ggplot(df, aes(x = shape_cell256_log2ratio, y = te256)) + geom_point() + scale_x_log10() + scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/shape_vs_TE/p6.png", plot = p6)

