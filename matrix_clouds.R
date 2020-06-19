### plot cloud from every stage in different stages (shape vs RNA)

## RNA: FPKM; poly(A): cell1, cell24, cell256, cell1K, 3.5h; total: cell1, 2h, 4h
## SHAPE: cell1, cell24, cell256, cell1K, oblong, oblong_CHX; SHAPES: cell1_invitro, cell24

## load data
load(file="/Volumes/USELESS/DATA/Shape-Seq/FPKM.Rdata")

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


# 'internal' structure
df_shape_cell1 <- data.frame(shape_cell1_log2ratio_internal = sapply(shape_cell1, function(x){mean(x$log2ratio[2:length(x)])}),
                                 n = names(shape_cell1))
df_shape_cell24 <- data.frame(shape_cell24_log2ratio_internal = sapply(shape_cell24, function(x){mean(x$log2ratio[2:length(x)])}),
                             n = names(shape_cell24))
df_shape_cell256 <- data.frame(shape_cell256_log2ratio_internal = sapply(shape_cell256, function(x){mean(x$log2ratio[2:length(x)])}),
                              n = names(shape_cell256))
df_shape_cell1K <- data.frame(shape_cell1K_log2ratio_internal = sapply(shape_cell1K, function(x){mean(x$log2ratio[2:length(x)])}),
                                  n = names(shape_cell1K))
df_shape_oblong <- data.frame(shape_oblong_log2ratio_internal = sapply(shape_oblong, function(x){mean(x$log2ratio[2:length(x)])}),
                              n = names(shape_oblong))
df_shape_oblong_CHX <- data.frame(shape_oblong_CHX_log2ratio_internal = sapply(shape_oblong_CHX, function(x){mean(x$log2ratio[2:length(x)])}),
                                  n = names(shape_oblong_CHX))

df_shapes_cell1_invitro <- data.frame(shapes_cell1_invitro_TC_internal = sapply(shapes_cell1_invitro, function(x){mean(x$TC[2:length(x)])}),
                                      n = names(shapes_cell1_invitro))
df_shapes_cell24 <- data.frame(shapes_cell24_TC_internal = sapply(shapes_cell24, function(x){mean(x$TC[2:length(x)])}),
                               n = names(shapes_cell24))

# 5'end coverage (on log2ratio / TC)
end_shape_cell1 <- data.frame(shape_cell1_log2ratio_5end = sapply(shape_cell1, function(x){x$log2ratio[1]}),
                             n = names(shape_cell1))
end_shape_cell24 <- data.frame(shape_cell24_log2ratio_5end = sapply(shape_cell24, function(x){x$log2ratio[1]}),
                              n = names(shape_cell24))
end_shape_cell256 <- data.frame(shape_cell256_log2ratio_5end = sapply(shape_cell256, function(x){x$log2ratio[1]}),
                               n = names(shape_cell256))
end_shape_cell1K <- data.frame(shape_cell1K_log2ratio_5end = sapply(shape_cell1K, function(x){x$log2ratio[1]}),
                              n = names(shape_cell1K))
end_shape_oblong <- data.frame(shape_oblong_log2ratio_5end = sapply(shape_oblong, function(x){x$log2ratio[1]}),
                              n = names(shape_oblong))
end_shape_oblong_CHX <- data.frame(shape_oblong_CHX_log2ratio_5end = sapply(shape_oblong_CHX, function(x){x$log2ratio[1]}),
                                  n = names(shape_oblong_CHX))

end_shapes_cell1_invitro <- data.frame(shapes_cell1_invitro_TC_5end = sapply(shapes_cell1_invitro, function(x){x$TC[1]}),
                                      n = names(shapes_cell1_invitro))
end_shapes_cell24 <- data.frame(shapes_cell24_TC_5end = sapply(shapes_cell24, function(x){x$TC[1]}),
                               n = names(shapes_cell24))

# 5'end coverage (on TC.treated)
endtr_shape_cell1 <- data.frame(shape_cell1_TC.treated_5end = sapply(shape_cell1, function(x){x$TC.treated[1]}),
                              n = names(shape_cell1))
endtr_shape_cell24 <- data.frame(shape_cell24_TC.treated_5end = sapply(shape_cell24, function(x){x$TC.treated[1]}),
                               n = names(shape_cell24))
endtr_shape_cell256 <- data.frame(shape_cell256_TC.treated_5end = sapply(shape_cell256, function(x){x$TC.treated[1]}),
                                n = names(shape_cell256))
endtr_shape_cell1K <- data.frame(shape_cell1K_TC.treated_5end = sapply(shape_cell1K, function(x){x$TC.treated[1]}),
                               n = names(shape_cell1K))
endtr_shape_oblong <- data.frame(shape_oblong_TC.treated_5end = sapply(shape_oblong, function(x){x$TC.treated[1]}),
                               n = names(shape_oblong))
endtr_shape_oblong_CHX <- data.frame(shape_oblong_CHX_TC.treated_5end = sapply(shape_oblong_CHX, function(x){x$TC.treated[1]}),
                                   n = names(shape_oblong_CHX))



## merge into data frame
df <- merge(df_shape_cell1, df_shape_cell24, by="n", all=TRUE)
df <- merge(df, df_shape_cell256, by="n", all=TRUE)
df <- merge(df, df_shape_cell1K, by="n", all=TRUE)
df <- merge(df, df_shape_oblong, by="n", all=TRUE)
df <- merge(df, df_shape_oblong_CHX, by="n", all=TRUE)
df <- merge(df, df_shapes_cell1_invitro, by="n", all=TRUE)
df <- merge(df, df_shapes_cell24, by="n", all=TRUE)

df <- merge(df, end_shape_cell1, by="n", all=TRUE)
df <- merge(df, end_shape_cell24, by="n", all=TRUE)
df <- merge(df, end_shape_cell256, by="n", all=TRUE)
df <- merge(df, end_shape_cell1K, by="n", all=TRUE)
df <- merge(df, end_shape_oblong, by="n", all=TRUE)
df <- merge(df, end_shape_oblong_CHX, by="n", all=TRUE)
df <- merge(df, end_shapes_cell1_invitro, by="n", all=TRUE)
df <- merge(df, end_shapes_cell24, by="n", all=TRUE)

df <- merge(df, endtr_shape_cell1, by="n", all=TRUE)
df <- merge(df, endtr_shape_cell24, by="n", all=TRUE)
df <- merge(df, endtr_shape_cell256, by="n", all=TRUE)
df <- merge(df, endtr_shape_cell1K, by="n", all=TRUE)
df <- merge(df, endtr_shape_oblong, by="n", all=TRUE)
df <- merge(df, endtr_shape_oblong_CHX, by="n", all=TRUE)

df <- merge(df, FPKM, by="n", all=TRUE)

## plot log2(shape) vs log2(RNA)
## compare internal vs 5end

Xaes = c("shape_cell1_log2ratio_internal", "shape_cell24_log2ratio_internal", "shape_cell256_log2ratio_internal", "shape_cell1K_log2ratio_internal", "shape_oblong_log2ratio_internal", "shape_oblong_CHX_log2ratio_internal", "shapes_cell1_invitro_TC_internal", "shapes_cell24_TC_internal")

Yaes = c("rna_cell1", "rna_cell24", "rna_cell256", "rna_cell1K", "rna_3_5h", "rna_3_5h", "rna_cell1", "rna_cell24")

for (i in 1:8) {
    ggplot(df, aes(x = log2(df[Xaes[i]]), y = log2(df[Yaes[i]]))) + geom_point(size=0.1) + xlab(Xaes[i]) + ylab(Yaes[i])
    fp <- paste0("/Volumes/USELESS/META/SHAPES_NEW/general/clouds/internal", i, ".png", collapse = NULL)
    ggsave(fp)
}

Xaes = c("shape_cell1_log2ratio_5end", "shape_cell24_log2ratio_5end", "shape_cell256_log2ratio_5end", "shape_cell1K_log2ratio_5end", "shape_oblong_log2ratio_5end", "shape_oblong_CHX_log2ratio_5end", "shapes_cell1_invitro_TC_5end", "shapes_cell24_TC_5end")

Yaes = c("rna_cell1", "rna_cell24", "rna_cell256", "rna_cell1K", "rna_3_5h", "rna_3_5h", "rna_cell1", "rna_cell24")

for (i in 1:8) {
  ggplot(df, aes(x = log2(df[Xaes[i]]), y = log2(df[Yaes[i]]))) + geom_point(size=0.1) + xlab(Xaes[i]) + ylab(Yaes[i])
  fp <- paste0("/Volumes/USELESS/META/SHAPES_NEW/general/clouds/5end", i, ".png", collapse = NULL)
  ggsave(fp)
}

Xaes = c("shape_cell1_TC.treated_5end", "shape_cell24_TC.treated_5end", "shape_cell256_TC.treated_5end", "shape_cell1K_TC.treated_5end", "shape_oblong_TC.treated_5end", "shape_oblong_CHX_TC.treated_5end", "shapes_cell1_invitro_TC_5end", "shapes_cell24_TC_5end")

Yaes = c("rna_cell1", "rna_cell24", "rna_cell256", "rna_cell1K", "rna_3_5h", "rna_3_5h", "rna_cell1", "rna_cell24")

for (i in 1:8) {
  ggplot(df, aes(x = log2(df[Xaes[i]]), y = log2(df[Yaes[i]]))) + geom_point(size=0.1) + xlab(Xaes[i]) + ylab(Yaes[i])
  fp <- paste0("/Volumes/USELESS/META/SHAPES_NEW/general/clouds/treated_5end", i, ".png", collapse = NULL)
  ggsave(fp)
}

## define what is 'cloud' and what is 'main' for each stage

# shape_cell1_log2ratio
ggplot(df, aes(x = log2(shape_cell1_log2ratio_internal), y = log2(rna_cell1))) + geom_point(size=0.1)
f <- function(x) {1.5*x + 18}
ggplot(df, aes(x = log2(shape_cell1_log2ratio_internal), y = log2(rna_cell1))) + geom_point(size=0.1) + stat_function(fun = f, colour = "red")
piggyback <- df$n[log2(df$rna_cell1) > (1.5*log2(df$shape_cell1_log2ratio_internal) + 18) & log2(df$shape_cell1_log2ratio_internal) > -15 & log2(df$shape_cell1_log2ratio_internal) < -8]
piggyback <- piggyback[!is.na(piggyback)]
piggyback <- data.frame(n = piggyback, cloud_cell1 = rep("cloud", length(piggyback)))
df <- merge(df, piggyback, by="n", all=TRUE)

# shape_cell24_log2ratio
ggplot(df, aes(x = log2(shape_cell24_log2ratio_internal), y = log2(rna_cell24))) + geom_point(size=0.1)
f <- function(x) {1.3*x + 13}
ggplot(df, aes(x = log2(shape_cell24_log2ratio_internal), y = log2(rna_cell24))) + geom_point(size=0.1) + stat_function(fun = f, colour = "red")
piggyback <- df$n[log2(df$rna_cell24) > (1.3*log2(df$shape_cell24_log2ratio_internal) + 13)]
piggyback <- piggyback[!is.na(piggyback)]
piggyback <- data.frame(n = piggyback, cloud_cell24 = rep("cloud", length(piggyback)))
df <- merge(df, piggyback, by="n", all=TRUE)

# shape_cell256_log2ratio
ggplot(df, aes(x = log2(shape_cell256_log2ratio_internal), y = log2(rna_cell256))) + geom_point(size=0.1)
f <- function(x) {1.2*x + 11}
ggplot(df, aes(x = log2(shape_cell256_log2ratio_internal), y = log2(rna_cell256))) + geom_point(size=0.1) + stat_function(fun = f, colour = "red")
piggyback <- df$n[log2(df$rna_cell256) > (1.2*log2(df$shape_cell256_log2ratio_internal) + 11)]
piggyback <- piggyback[!is.na(piggyback)]
piggyback <- data.frame(n = piggyback, cloud_cell256 = rep("cloud", length(piggyback)))
df <- merge(df, piggyback, by="n", all=TRUE)

# shape_cell1K_log2ratio
ggplot(df, aes(x = log2(shape_cell1K_log2ratio_internal), y = log2(rna_cell1K))) + geom_point(size=0.1)
f <- function(x) {1.5*x + 18}
ggplot(df, aes(x = log2(shape_cell1K_log2ratio_internal), y = log2(rna_cell1K))) + geom_point(size=0.1) + stat_function(fun = f, colour = "red")
piggyback <- df$n[log2(df$rna_cell1K) > (1.5*log2(df$shape_cell1K_log2ratio_internal) + 18)]
piggyback <- piggyback[!is.na(piggyback)]
piggyback <- data.frame(n = piggyback, cloud_cell1K = rep("cloud", length(piggyback)))
df <- merge(df, piggyback, by="n", all=TRUE)

# shape_oblong_log2ratio
ggplot(df, aes(x = log2(shape_oblong_log2ratio_internal), y = log2(rna_3_5h))) + geom_point(size=0.1)
f <- function(x) {1.3*x + 15}
ggplot(df, aes(x = log2(shape_oblong_log2ratio_internal), y = log2(rna_3_5h))) + geom_point(size=0.1) + stat_function(fun = f, colour = "red")
piggyback <- df$n[log2(df$rna_3_5h) > (1.3*log2(df$shape_oblong_log2ratio_internal) + 15)]
piggyback <- piggyback[!is.na(piggyback)]
piggyback <- data.frame(n = piggyback, cloud_oblong = rep("cloud", length(piggyback)))
df <- merge(df, piggyback, by="n", all=TRUE)

# shape_oblong_CHX_log2ratio
ggplot(df, aes(x = log2(shape_oblong_CHX_log2ratio_internal), y = log2(rna_3_5h))) + geom_point(size=0.1)
f <- function(x) {1.2*x + 10}
ggplot(df, aes(x = log2(shape_oblong_CHX_log2ratio_internal), y = log2(rna_3_5h))) + geom_point(size=0.1) + stat_function(fun = f, colour = "red")
piggyback <- df$n[log2(df$rna_3_5h) > (1.2*log2(df$shape_oblong_CHX_log2ratio_internal) + 10)]
piggyback <- piggyback[!is.na(piggyback)]
piggyback <- data.frame(n = piggyback, cloud_oblong_CHX = rep("cloud", length(piggyback)))
df <- merge(df, piggyback, by="n", all=TRUE)

# shapes_cell1_invitro_TC_internal
ggplot(df, aes(x = log2(shapes_cell1_invitro_TC_internal), y = log2(rna_cell1))) + geom_point(size=0.1)
f <- function(x) {1.5*x + 19}
ggplot(df, aes(x = log2(shapes_cell1_invitro_TC_internal), y = log2(rna_cell1))) + geom_point(size=0.1) + stat_function(fun = f, colour = "red")
piggyback <- df$n[log2(df$rna_cell1) > (1.5*log2(df$shapes_cell1_invitro_TC_internal) + 19) & log2(df$shapes_cell1_invitro_TC_internal) > -15 & log2(df$shapes_cell1_invitro_TC_internal) < -8]
piggyback <- piggyback[!is.na(piggyback)]
piggyback <- data.frame(n = piggyback, cloud_cell1_invitro = rep("cloud", length(piggyback)))
df <- merge(df, piggyback, by="n", all=TRUE)

# shapes_cell24_TC_internal
ggplot(df, aes(x = log2(shapes_cell24_TC_internal), y = log2(rna_cell24))) + geom_point(size=0.1)
f <- function(x) {1.4*x + 15}
ggplot(df, aes(x = log2(shapes_cell24_TC_internal), y = log2(rna_cell24))) + geom_point(size=0.1) + stat_function(fun = f, colour = "red")
piggyback <- df$n[log2(df$rna_cell24) > (1.4*log2(df$shapes_cell24_TC_internal) + 15)]
piggyback <- piggyback[!is.na(piggyback)]
piggyback <- data.frame(n = piggyback, cloud_cell24_TC = rep("cloud", length(piggyback)))
df <- merge(df, piggyback, by="n", all=TRUE)


#####
save(df, file="/Volumes/USELESS/META/SHAPES_NEW/general/clouds/cloud_df.Rsave")
load(file="/Volumes/USELESS/META/SHAPES_NEW/general/clouds/cloud_df.Rsave")


### matrix of plots
Xaes = c("shape_cell1_log2ratio_internal", "shape_cell24_log2ratio_internal", "shape_cell256_log2ratio_internal", "shape_cell1K_log2ratio_internal", "shape_oblong_log2ratio_internal", "shape_oblong_CHX_log2ratio_internal", "shapes_cell1_invitro_TC_internal", "shapes_cell24_TC_internal")

Yaes = c("rna_cell1", "rna_cell24", "rna_cell256", "rna_cell1K", "rna_3_5h", "rna_3_5h", "rna_cell1", "rna_cell24")

clouds = c("cloud_cell1", "cloud_cell24", "cloud_cell256", "cloud_cell1K", "cloud_oblong", "cloud_oblong_CHX", "cloud_cell1_invitro", "cloud_cell24_TC")

for (i in 1:8) {
  for (j in 1:8) {
    ggplot(df, aes(x = log2(df[Xaes[i]]), y = log2(df[Yaes[i]]), colour=df[,clouds[j]])) + geom_point(size=0.1) + theme(legend.position="none") + xlim(-15,0) + ylim(-10,15) + xlab(Xaes[i]) + ylab(Yaes[i]) + ggtitle(clouds[j])
    fp <- paste0("/Volumes/USELESS/META/SHAPES_NEW/general/clouds/matrix/p", i, "_", j, ".png", collapse = NULL)
    ggsave(fp)
  }
}


## for i in {1..8}; do convert p${i}_* +append p${i}.png; done

## convert p1.png p2.png p3.png p4.png p5.png p6.png p7.png p8.png -append matrix.png



go_cc <- read.csv(file="/Volumes/USELESS/DATA/genomes/GO/go_cc.tsv", sep="\t")
go_cc <- data.frame(n = go_cc$ID, GO_CC = go_cc$GO.Name)
df <- merge(df, go_cc, by="n", all=TRUE)

go_bp <- read.csv(file="/Volumes/USELESS/DATA/genomes/GO/go_bp.tsv", sep="\t")
go_bp <- data.frame(n = go_bp$ID, GO_BP = go_bp$GO.Name)
df <- merge(df, go_bp, by="n", all=TRUE)

go_mf <- read.csv(file="/Volumes/USELESS/DATA/genomes/GO/go_mf.tsv", sep="\t")
go_mf <- data.frame(n = go_mf$ID, GO_MF = go_mf$GO.Name)
df <- merge(df, go_mf, by="n", all=TRUE)

df <- merge(df, stop_codons, by="n", all=TRUE)



######### mini cloud on accessibility vs ribo
df$accessibility24 <- df$shape_cell24_log2ratio_internal / df$rna_cell24
df$accessibility256 <- df$shape_cell256_log2ratio_internal / df$rna_cell256
df$accessibility1K <- df$shape_cell1K_log2ratio_internal / df$rna_cell1K

ggplot(df, aes(x = log2(accessibility24), y = log2(ribo_cell24))) + geom_point(size=0.1)
a24 <- df[log2(df$accessibility24) > -2.5 & log2(df$ribo_cell24) > 2.5 & log2(df$accessibility24) < 3,]
a24 <- unique(a24[!is.na(a24$accessibility24) & !is.na(a24$ribo_cell24),]$gene)
## histones! and some from zgc (zebrafish gene collection) and si (sanger institute)

ggplot(df, aes(x = log2(accessibility256), y = log2(ribo_cell256))) + geom_point(size=0.1)
a256 <- df[log2(df$accessibility256) > -2.5 & log2(df$ribo_cell256) > 2.5 & log2(df$accessibility256) < 3,]
a256 <- unique(a256[!is.na(a256$accessibility256) & !is.na(a256$ribo_cell256),]$gene)

ggplot(df, aes(x = log2(accessibility1K), y = log2(ribo_cell1K))) + geom_point(size=0.1)
a1K <- df[log2(df$accessibility1K) > -5 & log2(df$ribo_cell1K) > 5 & log2(df$accessibility1K) < 0,]
a1K <- unique(a1K[!is.na(a1K$accessibility1K) & !is.na(a1K$ribo_cell1K),]$gene)
