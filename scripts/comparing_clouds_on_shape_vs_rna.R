### plot cloud from every stage in different stages (shape vs RNA) (include ribo??)

## load data
# 1st compare cell1_et1 & cell1_et2; cell1K_et1 & cell1K_et2; oblong & oblong_CHX
shape_cell1_et1 <- get(load(file = "/Volumes/USELESS/test2/normalized/cell1_invivo_et1.Rsave"))
shape_cell1_et2 <- get(load(file = "/Volumes/USELESS/test2/normalized/cell1_invivo_et2.Rsave"))
shape_cell1K_et1 <- get(load(file = "/Volumes/USELESS/test2/normalized/cell1K_invivo_et1.Rsave"))
shape_cell1K_et2 <- get(load(file = "/Volumes/USELESS/test2/normalized/cell1K_invivo_et2.Rsave"))
shape_oblong <- get(load(file = "/Volumes/USELESS/test/normalized/shape_oblong.Rsave"))
shape_oblong_CHX <- get(load(file = "/Volumes/USELESS/test/normalized/shape_oblong_CHX.Rsave"))
rm(shapeseq_norm)

df_shape_cell1_et1 <- data.frame(shape_cell1_et1_log2ratio = sapply(shape_cell1_et1, function(x){mean(x$log2ratio[2:length(x)])}),
                                 n = names(shape_cell1_et1))
df_shape_cell1_et2 <- data.frame(shape_cell1_et2_log2ratio = sapply(shape_cell1_et2, function(x){mean(x$log2ratio[2:length(x)])}),
                                 n = names(shape_cell1_et2))
df_shape_cell1K_et1 <- data.frame(shape_cell1K_et1_log2ratio = sapply(shape_cell1K_et1, function(x){mean(x$log2ratio[2:length(x)])}),
                                 n = names(shape_cell1K_et1))
df_shape_cell1K_et2 <- data.frame(shape_cell1K_et2_log2ratio = sapply(shape_cell1K_et2, function(x){mean(x$log2ratio[2:length(x)])}),
                                 n = names(shape_cell1K_et2))
df_shape_oblong <- data.frame(shape_oblong_log2ratio = sapply(shape_oblong, function(x){mean(x$log2ratio[2:length(x)])}),
                                 n = names(shape_oblong))
df_shape_oblong_CHX <- data.frame(shape_oblong_CHX_log2ratio = sapply(shape_oblong_CHX, function(x){mean(x$log2ratio[2:length(x)])}),
                              n = names(shape_oblong_CHX))

# RNA: FPKM_24, FPKM_256, FPKM_1K, fpkm_lee (4h)
## ribo-seq & RNA-seq
load(file="/Volumes/USELESS/META/SHAPES/FPKM_24.Rdata")
load(file="/Volumes/USELESS/META/SHAPES/FPKM_256.Rdata")
load(file="/Volumes/USELESS/META/SHAPES/FPKM_1K.Rdata")
load(file="/Volumes/USELESS/META/SHAPES/FPKM_Lee.Rsave")
colnames(lee_RNA_FPKM) <- c("total_rna_2h", "total_rna_4h")
lee_RNA_FPKM$n <- rownames(lee_RNA_FPKM)

RNA_fpkm <- data.frame(rna_cell24 = FPKM_24$exons_RNA_fpkm, rna_cell256 = FPKM_256$exons_RNA_fpkm, rna_cell1K = FPKM_1K$exons_RNA_fpkm,
                       n = rownames(FPKM_1K))
RNA_fpkm <- merge(RNA_fpkm, lee_RNA_FPKM, by="n", all=TRUE)

ribo_fpkm <- data.frame(ribo_cell24 = FPKM_24$exons_Ribo_fpkm, ribo_cell256 = FPKM_256$exons_Ribo_fpkm,
                        ribo_cell1K = FPKM_1K$exons_Ribo_fpkm, n = rownames(FPKM_1K))

# SHAPE(S): cell1_invitro, cell1, cell24, cell24_N3, cell256, cell1K, oblong
# shapes
shapes_cell24 <- get(load(file = "/Volumes/USELESS/test/normalized/cell2_4_invivo_NAI_N3.Rsave"))
shapes_cell1_invitro <- get(load(file = "/Volumes/USELESS/test/normalized/cell1_invitro_NAI_N3.Rsave"))
shapes_cell1_invitro_nonsel <- get(load(file = "/Volumes/USELESS/test/normalized/cell1_invitro_NAI_N3_nonsel.Rsave"))
rm(shapes_norm)

shape_cell24 <- get(load(file = "/Volumes/USELESS/test/normalized/2-4cell.Rsave"))
u <- unlist(shape_cell24)
u$log2ratio[u$log2ratio < 0] <- 0
shape_cell24 <- split(u, seqnames(u))
shape_cell24 <- shape_cell24[sapply(shape_cell24, function(x){length(x) > 0})]
shape_cell256 <- get(load(file = "/Volumes/USELESS/test/normalized/256cell.Rsave"))
u <- unlist(shape_cell256)
u$log2ratio[u$log2ratio < 0] <- 0
shape_cell256 <- split(u, seqnames(u))
shape_cell256 <- shape_cell256[sapply(shape_cell256, function(x){length(x) > 0})]
rm(shapeseq_norm)

df_shapes_cell24 <- data.frame(shapes_cell24_TC = sapply(shapes_cell24, function(x){mean(x$TC[2:length(x)])}),
                               n = names(shapes_cell24))
df_shapes_cell1_invitro <- data.frame(shapes_cell1_invitro_TC = sapply(shapes_cell1_invitro, function(x){mean(x$TC[2:length(x)])}),
                                      n = names(shapes_cell1_invitro))
df_shapes_cell1_invitro_nonsel <- data.frame(shapes_cell1_invitro_nonsel_TC = sapply(shapes_cell1_invitro_nonsel,
                                            function(x){mean(x$TC[2:length(x)])}), n = names(shapes_cell1_invitro_nonsel))

df_shape_cell24 <- data.frame(shape_cell24_log2ratio = sapply(shape_cell24, function(x){mean(x$log2ratio[2:length(x)])}),
                              n = names(shape_cell24))
df_shape_cell256 <- data.frame(shape_cell256_log2ratio = sapply(shape_cell256, function(x){mean(x$log2ratio[2:length(x)])}),
                               n = names(shape_cell256))

## merge into data frame
df <- merge(df_shape_cell1_et1, df_shape_cell1_et2, by="n", all=TRUE)
df <- merge(df, df_shape_cell1K_et1, by="n", all=TRUE)
df <- merge(df, df_shape_cell1K_et2, by="n", all=TRUE)
df <- merge(df, df_shape_oblong, by="n", all=TRUE)
df <- merge(df, df_shape_oblong_CHX, by="n", all=TRUE)

df <- merge(df, df_shapes_cell24, by="n", all=TRUE)
df <- merge(df, df_shapes_cell1_invitro, by="n", all=TRUE)
df <- merge(df, df_shapes_cell1_invitro_nonsel, by="n", all=TRUE)
df <- merge(df, df_shape_cell24, by="n", all=TRUE)
df <- merge(df, df_shape_cell256, by="n", all=TRUE)

df <- merge(df, RNA_fpkm, by="n", all=TRUE)
df <- merge(df, ribo_fpkm, by="n", all=TRUE)

## plot log2(shape) vs log2(RNA)
## define what is 'cloud' and what is 'main' for each stage

# shape_cell1_et1_log2ratio
ggplot(df, aes(x = log2(shape_cell1_et1_log2ratio), y = log2(rna_cell24))) + geom_point(size=0.1)
f <- function(x) {1.2*x + 16}
ggplot(df, aes(x = log2(shape_cell1_et1_log2ratio), y = log2(rna_cell24))) + geom_point(size=0.1) + geom_point(size=0.1) + stat_function(fun = f, colour = "red")
piggyback <- df$n[log2(df$rna_cell24) > (1.2*log2(df$shape_cell1_et1_log2ratio) + 16) & log2(df$shape_cell1_et1_log2ratio) < -9 & log2(df$shape_cell1_et1_log2ratio) > -13]
piggyback <- piggyback[!is.na(piggyback)]
piggyback <- data.frame(n = piggyback, cloud_cell1_et1 = rep("cloud", length(piggyback)))
df <- merge(df, piggyback, by="n", all=TRUE)

# shape_cell1_et2_log2ratio
ggplot(df, aes(x = log2(shape_cell1_et2_log2ratio), y = log2(rna_cell24))) + geom_point(size=0.1)
f <- function(x) {1.2*x + 15}
ggplot(df, aes(x = log2(shape_cell1_et2_log2ratio), y = log2(rna_cell24))) + geom_point(size=0.1) + geom_point(size=0.1) + stat_function(fun = f, colour = "red")
piggyback <- df$n[log2(df$rna_cell24) > (1.2*log2(df$shape_cell1_et2_log2ratio) + 15) & log2(df$shape_cell1_et2_log2ratio) < -8 & log2(df$shape_cell1_et2_log2ratio) > -14]
piggyback <- piggyback[!is.na(piggyback)]
piggyback <- data.frame(n = piggyback, cloud_cell1_et2 = rep("cloud", length(piggyback)))
df <- merge(df, piggyback, by="n", all=TRUE)

# shape_cell1K_et1_log2ratio
ggplot(df, aes(x = log2(shape_cell1K_et1_log2ratio), y = log2(rna_cell1K))) + geom_point(size=0.1)
f <- function(x) {1.2*x + 16}
ggplot(df, aes(x = log2(shape_cell1K_et1_log2ratio), y = log2(rna_cell1K))) + geom_point(size=0.1) + geom_point(size=0.1) + stat_function(fun = f, colour = "red")
piggyback <- df$n[log2(df$rna_cell1K) > (1.2*log2(df$shape_cell1K_et1_log2ratio) + 16) & log2(df$shape_cell1K_et1_log2ratio) < -9 & log2(df$shape_cell1K_et1_log2ratio) > -13]
piggyback <- piggyback[!is.na(piggyback)]
piggyback <- data.frame(n = piggyback, cloud_cell1K_et1 = rep("cloud", length(piggyback)))
df <- merge(df, piggyback, by="n", all=TRUE)

# shape_cell1K_et2_log2ratio
ggplot(df, aes(x = log2(shape_cell1K_et2_log2ratio), y = log2(rna_cell1K))) + geom_point(size=0.1)
f <- function(x) {1.2*x + 16}
ggplot(df, aes(x = log2(shape_cell1K_et2_log2ratio), y = log2(rna_cell1K))) + geom_point(size=0.1) + geom_point(size=0.1) + stat_function(fun = f, colour = "red")
piggyback <- df$n[log2(df$rna_cell1K) > (1.2*log2(df$shape_cell1K_et2_log2ratio) + 16) & log2(df$shape_cell1K_et2_log2ratio) < -9 & log2(df$shape_cell1K_et2_log2ratio) > -13]
piggyback <- piggyback[!is.na(piggyback)]
piggyback <- data.frame(n = piggyback, cloud_cell1K_et2 = rep("cloud", length(piggyback)))
df <- merge(df, piggyback, by="n", all=TRUE)

# shape_oblong_log2ratio
ggplot(df, aes(x = log2(shape_oblong_log2ratio), y = log2(total_rna_4h))) + geom_point(size=0.1)
f <- function(x) {1.2*x + 7}
ggplot(df, aes(x = log2(shape_oblong_log2ratio), y = log2(total_rna_4h))) + geom_point(size=0.1) + geom_point(size=0.1) + stat_function(fun = f, colour = "red")
piggyback <- df$n[log2(df$total_rna_4h) > (1.2*log2(df$shape_oblong_log2ratio) + 7) & log2(df$shape_oblong_log2ratio) < -7.5]
piggyback <- piggyback[!is.na(piggyback)]
piggyback <- data.frame(n = piggyback, cloud_oblong = rep("cloud", length(piggyback)))
df <- merge(df, piggyback, by="n", all=TRUE)

# shape_oblong_CHX_log2ratio
ggplot(df, aes(x = log2(shape_oblong_CHX_log2ratio), y = log2(total_rna_4h))) + geom_point(size=0.1)
f <- function(x) {1.2*x + 3}
ggplot(df, aes(x = log2(shape_oblong_CHX_log2ratio), y = log2(total_rna_4h))) + geom_point(size=0.1) + geom_point(size=0.1) + stat_function(fun = f, colour = "red")
piggyback <- df$n[log2(df$total_rna_4h) > (1.2*log2(df$shape_oblong_CHX_log2ratio) + 3) & log2(df$shape_oblong_CHX_log2ratio) < -5]
piggyback <- piggyback[!is.na(piggyback)]
piggyback <- data.frame(n = piggyback, cloud_oblong_CHX = rep("cloud", length(piggyback)))
df <- merge(df, piggyback, by="n", all=TRUE)

# shapes_cell24_TC
ggplot(df, aes(x = log2(shapes_cell24_TC), y = log2(rna_cell24))) + geom_point(size=0.1)
f <- function(x) {1.2*x + 13}
ggplot(df, aes(x = log2(shapes_cell24_TC), y = log2(rna_cell24))) + geom_point(size=0.1) + geom_point(size=0.1) + stat_function(fun = f, colour = "red")
piggyback <- df$n[log2(df$rna_cell24) > (1.2*log2(df$shapes_cell24_TC) + 13) & log2(df$shapes_cell24_TC) < -5 & log2(df$shapes_cell24_TC) > -15]
piggyback <- piggyback[!is.na(piggyback)]
piggyback <- data.frame(n = piggyback, cloud_cell24_TC = rep("cloud", length(piggyback)))
df <- merge(df, piggyback, by="n", all=TRUE)

# shapes_cell1_invitro_TC
ggplot(df, aes(x = log2(shapes_cell1_invitro_TC), y = log2(rna_cell24))) + geom_point(size=0.1)
f <- function(x) {1.2*x + 16}
ggplot(df, aes(x = log2(shapes_cell1_invitro_TC), y = log2(rna_cell24))) + geom_point(size=0.1) + geom_point(size=0.1) + stat_function(fun = f, colour = "red")
piggyback <- df$n[log2(df$rna_cell24) > (1.2*log2(df$shapes_cell1_invitro_TC) + 16) & log2(df$shapes_cell1_invitro_TC) < -8 & log2(df$shapes_cell1_invitro_TC) > -15]
piggyback <- piggyback[!is.na(piggyback)]
piggyback <- data.frame(n = piggyback, cloud_cell1_invitro_TC = rep("cloud", length(piggyback)))
df <- merge(df, piggyback, by="n", all=TRUE)

# shapes_cell1_invitro_nonsel_TC
ggplot(df, aes(x = log2(shapes_cell1_invitro_nonsel_TC), y = log2(rna_cell24))) + geom_point(size=0.1)
f <- function(x) {1.2*x + 16}
ggplot(df, aes(x = log2(shapes_cell1_invitro_nonsel_TC), y = log2(rna_cell24))) + geom_point(size=0.1) + geom_point(size=0.1) + stat_function(fun = f, colour = "red")
piggyback <- df$n[log2(df$rna_cell24) > (1.2*log2(df$shapes_cell1_invitro_nonsel_TC) + 16) & log2(df$shapes_cell1_invitro_nonsel_TC) < -8 & log2(df$shapes_cell1_invitro_nonsel_TC) > -15]
piggyback <- piggyback[!is.na(piggyback)]
piggyback <- data.frame(n = piggyback, cloud_cell1_invitro_nonsel_TC = rep("cloud", length(piggyback)))
df <- merge(df, piggyback, by="n", all=TRUE)

# shape_cell24_log2ratio
ggplot(df, aes(x = log2(shape_cell24_log2ratio), y = log2(rna_cell24))) + geom_point(size=0.1)
f <- function(x) {1.2*x + 14}
ggplot(df, aes(x = log2(shape_cell24_log2ratio), y = log2(rna_cell24))) + geom_point(size=0.1) + geom_point(size=0.1) + stat_function(fun = f, colour = "red")
piggyback <- df$n[log2(df$rna_cell24) > (1.2*log2(df$shape_cell24_log2ratio) + 14) & log2(df$shape_cell24_log2ratio) < -8 & log2(df$shape_cell24_log2ratio) > -15]
piggyback <- piggyback[!is.na(piggyback)]
piggyback <- data.frame(n = piggyback, cloud_cell24 = rep("cloud", length(piggyback)))
df <- merge(df, piggyback, by="n", all=TRUE)

# shape_cell256_log2ratio
ggplot(df, aes(x = log2(shape_cell256_log2ratio), y = log2(rna_cell256))) + geom_point(size=0.1)
f <- function(x) {1.2*x + 12}
ggplot(df, aes(x = log2(shape_cell256_log2ratio), y = log2(rna_cell256))) + geom_point(size=0.1) + geom_point(size=0.1) + stat_function(fun = f, colour = "red")
piggyback <- df$n[log2(df$rna_cell256) > (1.2*log2(df$shape_cell256_log2ratio) + 12) & log2(df$shape_cell256_log2ratio) < -8 & log2(df$shape_cell256_log2ratio) > -15]
piggyback <- piggyback[!is.na(piggyback)]
piggyback <- data.frame(n = piggyback, cloud_cell256 = rep("cloud", length(piggyback)))
df <- merge(df, piggyback, by="n", all=TRUE)


#####
save(df, file="/Volumes/USELESS/META/SHAPES_NEW/general/piggyback/cloud_df.Rsave")
load(file="/Volumes/USELESS/META/SHAPES_NEW/general/piggyback/cloud_df.Rsave")

df <- merge(df, gn, by="n", all=TRUE)

# ggplot(df, aes(x = log2(shape_cell1_et2_log2ratio), y = log2(rna_cell24))) + geom_point(size=0.1) + stat_function(fun = f, colour = "red") + geom_point(size=0.1, data = subset(df, cloud_cell1_et2 == 'cloud'), aes(x=log2(shape_cell1_et2_log2ratio), y=log2(rna_cell24), colour=cloud_cell1_et2)) + theme(legend.position="none") + xlim(-15,0) + ylim(-10,15)


## matrix of plots

Xaes = c("shapes_cell1_invitro_TC", "shapes_cell1_invitro_nonsel_TC", "shape_cell1_et1_log2ratio", "shape_cell1_et2_log2ratio", "shape_cell24_log2ratio", "shapes_cell24_TC", "shape_cell256_log2ratio", "shape_cell1K_et1_log2ratio", "shape_cell1K_et2_log2ratio", "shape_oblong_log2ratio", "shape_oblong_CHX_log2ratio")

Yaes = c("rna_cell24", "rna_cell24", "rna_cell24", "rna_cell24", "rna_cell24", "rna_cell24", "rna_cell256", "rna_cell1K", "rna_cell1K", "total_rna_4h", "total_rna_4h")

clouds = c("cloud_cell1_invitro_TC", "cloud_cell1_invitro_nonsel_TC", "cloud_cell1_et1", "cloud_cell1_et2", "cloud_cell24", "cloud_cell24_TC", "cloud_cell256", "cloud_cell1K_et1", "cloud_cell1K_et2", "cloud_oblong", "cloud_oblong_CHX")


for (i in 1:11) {
  for (j in 1:11) {
    ggplot(df, aes(x = log2(df[Xaes[i]]), y = log2(df[Yaes[i]]), colour=df[,clouds[j]])) + geom_point(size=0.1) + theme(legend.position="none") + xlim(-15,0) + ylim(-10,15) + xlab(Xaes[i]) + ylab(Yaes[i]) + ggtitle(clouds[j])
    fp <- paste0("/Volumes/USELESS/META/SHAPES_NEW/general/piggyback/clouds/p", i, "_", j, ".png", collapse = NULL)
    ggsave(fp)
  }
}

## for i in {1..11}; do convert p${i}_1.png p${i}_2.png p${i}_3.png p${i}_4.png p${i}_5.png p${i}_6.png p${i}_7.png p${i}_8.png p${i}_9.png p${i}_10.png p${i}_11.png +append p${i}.png; done

## convert p1.png p2.png p3.png p4.png p5.png p6.png p7.png p8.png p9.png p10.png p11.png -append gene.png


## ensembl_ID, gene_name
ensembl_gn <- read.csv(file="/Volumes/USELESS/DATA/genomes/ensembl_gene_name.txt", header = T)
gn <- data.frame(n = ensembl_gn$Transcript.stable.ID, gene_name = ensembl_gn$Gene.name)

### common tx from the cloud
partdf <- data.frame(n = df$n, cloud_cell1K_et1 = df$cloud_cell1K_et1, cloud_cell256 = df$cloud_cell256, cloud_oblong = df$cloud_oblong)
partdf <- merge(partdf, gn, by="n", all=TRUE)
partdf <- partdf[complete.cases(partdf),]

candidates <- unique(partdf$gene_name)

for (i in 1:11) {
  for (j in 1:11) {
    ggplot(df, aes(x = log2(df[Xaes[i]]), y = log2(df[Yaes[i]]))) + geom_point(size=0.1) + theme(legend.position="none") + xlab(Xaes[i]) + ylab(Yaes[i]) + ggtitle(clouds[j]) + geom_point(data=subset(df, gene_name == "senp3a"), aes(x = log2(subset(df, gene_name == "senp3a")[Xaes[i]]), y = log2(subset(df, gene_name == "senp3a")[Yaes[i]]), colour="green"))
    fp <- paste0("/Volumes/USELESS/META/SHAPES_NEW/general/piggyback/temp/p", i, "_", j, ".png", collapse = NULL)
    ggsave(fp)
  }
}


ggplot(d, aes(x = a, y = a)) + geom_point() + geom_point(data=subset(d, a1 == "dupa"), aes(x=a, y=a), colour=a1)

for (i in 1:11) {
    ggplot(df, aes(x = log2(df[Xaes[i]]), y = log2(df[Yaes[i]]))) + geom_point(size=0.1) + theme(legend.position="none") + xlab(Xaes[i]) + ylab(Yaes[i]) + geom_point(data=subset(df, gene_name == "clec16a"), aes(x = log2(subset(df, gene_name == "clec16a")[Xaes[i]]), y = log2(subset(df, gene_name == "clec16a")[Yaes[i]]), colour="green"))
    fp <- paste0("/Volumes/USELESS/META/SHAPES_NEW/general/piggyback/temp/p", i, ".png", collapse = NULL)
    ggsave(fp)
}

### for ribo

Xaes = c("shapes_cell1_invitro_TC", "shapes_cell1_invitro_nonsel_TC", "shape_cell1_et1_log2ratio", "shape_cell1_et2_log2ratio", "shape_cell24_log2ratio", "shapes_cell24_TC", "shape_cell256_log2ratio", "shape_cell1K_et1_log2ratio", "shape_cell1K_et2_log2ratio")

Yaes = c("ribo_cell24", "ribo_cell24", "ribo_cell24", "ribo_cell24", "ribo_cell24", "ribo_cell24", "ribo_cell256", "ribo_cell1K", "ribo_cell1K")

gene = "ulk4"
for (i in 1:9) {
  ggplot(df, aes(x = log2(df[Xaes[i]]), y = log2(df[Yaes[i]]))) + geom_point(size=0.1) + theme(legend.position="none") + xlab(Xaes[i]) + ylab(Yaes[i]) + geom_point(data=subset(df, gene_name == gene), aes(x = log2(subset(df, gene_name == gene)[Xaes[i]]), y = log2(subset(df, gene_name == gene)[Yaes[i]]), colour="green"))
  fp <- paste0("/Volumes/USELESS/META/SHAPES_NEW/general/piggyback/temp/p", i, ".png", collapse = NULL)
  ggsave(fp)
}


#######
load(file = "/Volumes/USELESS/META/SHAPES/FPKM_1cell.Rdata")
cell1 <- data.frame(n = FPKM_1cell$n, rna_cell1 = FPKM_1cell$exons_RNA_fpkm)
df <- merge(df, cell1, by="n", all=TRUE)

ggplot(df, aes(x = log2(shape_cell1_et1_log2ratio), y = log2(rna_cell1))) + geom_point(size=0.1)

