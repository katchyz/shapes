### mean shape vs RNA-seq
library(ggplot2)
library(gridExtra)

# shape
shape_cell24 <- get(load(file = "/Volumes/USELESS/test/normalized/2-4cell.Rsave"))
shape_cell256 <- get(load(file = "/Volumes/USELESS/test/normalized/256cell.Rsave"))
shape_oblong <- get(load(file = "/Volumes/USELESS/test/normalized/oblong_invivo.Rsave"))
shape_oblong_CHX <- get(load(file = "/Volumes/USELESS/test/normalized/oblong_CHX_invivo.Rsave"))
rm(shapeseq_norm)
# shapes
shapes_cell24 <- get(load(file = "/Volumes/USELESS/test/normalized/cell2_4_invivo_NAI_N3.Rsave"))
shapes_cell1_invitro <- get(load(file = "/Volumes/USELESS/test/normalized/cell1_invitro_NAI_N3.Rsave"))
shapes_cell1_invitro_nonsel <- get(load(file = "/Volumes/USELESS/test/normalized/cell1_invitro_NAI_N3_nonsel.Rsave"))
rm(shapes_norm)

# RNA FPKM
load(file="/Volumes/USELESS/META/SHAPES/FPKM_Lee.Rsave")
lee_RNA_FPKM$n <- rownames(lee_RNA_FPKM)
load(file="/Volumes/USELESS/META/SHAPES/FPKM_24.Rdata")
FPKM_24 <- data.frame(cell24 = FPKM_24$exons_RNA_fpkm, n = rownames(FPKM_24))

######
df_shape_cell24 <- data.frame(shape_cell24_TC.control = sapply(shape_cell24, function(x){mean(x$TC.control)}),
                              shape_cell24_TC.treated = sapply(shape_cell24, function(x){mean(x$TC.treated)}),
                              shape_cell24_log2ratio = sapply(shape_cell24, function(x){mean(x$log2ratio)}),
                              n = names(shape_cell24))

df_shape_cell256 <- data.frame(shape_cell256_TC.control = sapply(shape_cell256, function(x){mean(x$TC.control)}),
                              shape_cell256_TC.treated = sapply(shape_cell256, function(x){mean(x$TC.treated)}),
                              shape_cell256_log2ratio = sapply(shape_cell256, function(x){mean(x$log2ratio)}),
                              n = names(shape_cell256))

df_shape_oblong <- data.frame(shape_oblong_TC.control = sapply(shape_oblong, function(x){mean(x$TC.control)}),
                              shape_oblong_TC.treated = sapply(shape_oblong, function(x){mean(x$TC.treated)}),
                              shape_oblong_log2ratio = sapply(shape_oblong, function(x){mean(x$log2ratio)}),
                              n = names(shape_oblong))

df_shape_oblong_CHX <- data.frame(shape_oblong_CHX_TC.control = sapply(shape_oblong_CHX, function(x){mean(x$TC.control)}),
                              shape_oblong_CHX_TC.treated = sapply(shape_oblong_CHX, function(x){mean(x$TC.treated)}),
                              shape_oblong_CHX_log2ratio = sapply(shape_oblong_CHX, function(x){mean(x$log2ratio)}),
                              n = names(shape_oblong_CHX))

df_shape <- merge(df_shape_cell24, df_shape_cell256, by="n", all=TRUE)
df_shape <- merge(df_shape, df_shape_oblong, by="n", all=TRUE)
df_shape <- merge(df_shape, df_shape_oblong_CHX, by="n", all=TRUE)


df_shapes_cell24 <- data.frame(shapes_cell24_TC = sapply(shapes_cell24, function(x){mean(x$TC)}),
                              shapes_cell24_winsor = sapply(shapes_cell24, function(x){mean(x$winsor)}),
                              n = names(shapes_cell24))

df_shapes_cell1_invitro <- data.frame(shapes_cell1_invitro_TC = sapply(shapes_cell1_invitro, function(x){mean(x$TC)}),
                               shapes_cell1_invitro_winsor = sapply(shapes_cell1_invitro, function(x){mean(x$winsor)}),
                               n = names(shapes_cell1_invitro))

df_shapes_cell1_invitro_nonsel <- data.frame(shapes_cell1_invitro_nonsel_TC = sapply(shapes_cell1_invitro_nonsel, function(x){mean(x$TC)}),
                              shapes_cell1_invitro_nonsel_winsor = sapply(shapes_cell1_invitro_nonsel, function(x){mean(x$winsor)}),
                              n = names(shapes_cell1_invitro_nonsel))

df_shapes <- merge(df_shapes_cell24, df_shapes_cell1_invitro, by="n", all=TRUE)
df_shapes <- merge(df_shapes, df_shapes_cell1_invitro_nonsel, by="n", all=TRUE)

df <- merge(df_shape, df_shapes, by="n", all=TRUE)
df <- merge(df, lee_RNA_FPKM, by="n", all=TRUE)
df <- merge(df, FPKM_24, by="n", all=TRUE)

## pick only protein-coding genes, check for correlation (why?? what happens to rRNA?)
## colour code by protein-coding /rRNA

# gtf_annot <- import("/Volumes/USELESS/DATA/genomes/GTF/danio_rerio/Danio_rerio.GRCz10.84.gtf", format = "gtf")
# gtf_annot <- gtf_annot[gtf_annot$type == "transcript"]
# transcript_df <- data.frame(n = gtf_annot$transcript_id, tx_biotype = gtf_annot$transcript_biotype)

# df <- merge(df, transcript_df, by="n", all=TRUE)

##############################################################################################################
# > colnames(df)
# [1] "n"                                  "shape_cell24_TC.control"            "shape_cell24_TC.treated"           
# [4] "shape_cell24_log2ratio"             "shape_cell256_TC.control"           "shape_cell256_TC.treated"          
# [7] "shape_cell256_log2ratio"            "shape_oblong_TC.control"            "shape_oblong_TC.treated"           
# [10] "shape_oblong_log2ratio"             "shape_oblong_CHX_TC.control"        "shape_oblong_CHX_TC.treated"       
# [13] "shape_oblong_CHX_log2ratio"         "shapes_cell24_TC"                   "shapes_cell24_winsor"              
# [16] "shapes_cell1_invitro_TC"            "shapes_cell1_invitro_winsor"        "shapes_cell1_invitro_nonsel_TC"    
# [19] "shapes_cell1_invitro_nonsel_winsor" "h2"                                 "h4"                                
# [22] "cell24"                             "tx_biotype"
##############################################################################################################

# ggplot(df, aes(x = shape_cell24_TC.control, y = cell24)) + geom_point() + scale_x_log10() + scale_y_log10() +
#   geom_point(data = subset(df, tx_biotype == 'rRNA'), aes(x=shape_cell24_TC.control, y=cell24, colour=tx_biotype)) +
#   theme(legend.position="none")

p1 <- ggplot(df, aes(x = shape_cell24_TC.control, y = cell24)) + geom_point() + scale_x_log10() + scale_y_log10()
p2 <- ggplot(df, aes(x = shape_cell24_TC.treated, y = cell24)) + geom_point() + scale_x_log10() + scale_y_log10()
p3 <- ggplot(df, aes(x = shape_cell24_log2ratio, y = cell24)) + geom_point() + scale_x_log10() + scale_y_log10()
p4 <- ggplot(df, aes(x = shape_cell256_TC.control, y = h2)) + geom_point() + scale_x_log10() + scale_y_log10()
p5 <- ggplot(df, aes(x = shape_cell256_TC.treated, y = h2)) + geom_point() + scale_x_log10() + scale_y_log10()
p6 <- ggplot(df, aes(x = shape_cell256_log2ratio, y = h2)) + geom_point() + scale_x_log10() + scale_y_log10()
p7 <- ggplot(df, aes(x = shape_oblong_TC.control, y = h4)) + geom_point() + scale_x_log10() + scale_y_log10()
p8 <- ggplot(df, aes(x = shape_oblong_TC.treated, y = h4)) + geom_point() + scale_x_log10() + scale_y_log10()
p9 <- ggplot(df, aes(x = shape_cell24_log2ratio, y = h4)) + geom_point() + scale_x_log10() + scale_y_log10()
p10 <- ggplot(df, aes(x = shape_oblong_CHX_TC.control, y = h4)) + geom_point() + scale_x_log10() + scale_y_log10()
p11 <- ggplot(df, aes(x = shape_oblong_CHX_TC.treated, y = h4)) + geom_point() + scale_x_log10() + scale_y_log10()
p12 <- ggplot(df, aes(x = shape_oblong_CHX_log2ratio, y = h4)) + geom_point() + scale_x_log10() + scale_y_log10()
p13 <- ggplot(df, aes(x = shapes_cell24_TC, y = cell24)) + geom_point() + scale_x_log10() + scale_y_log10()
p14 <- ggplot(df, aes(x = shapes_cell24_winsor, y = cell24)) + geom_point() + scale_x_log10() + scale_y_log10()
p15 <- ggplot(df, aes(x = shapes_cell1_invitro_TC, y = cell24)) + geom_point() + scale_x_log10() + scale_y_log10()
p16 <- ggplot(df, aes(x = shapes_cell1_invitro_winsor, y = cell24)) + geom_point() + scale_x_log10() + scale_y_log10()
p17 <- ggplot(df, aes(x = shapes_cell1_invitro_nonsel_TC, y = cell24)) + geom_point() + scale_x_log10() + scale_y_log10()
p18 <- ggplot(df, aes(x = shapes_cell1_invitro_nonsel_winsor, y = cell24)) + geom_point() + scale_x_log10() + scale_y_log10()

#p <- c(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18)
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/p1.png", plot = p1)
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/p2.png", plot = p2)
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/p3.png", plot = p3)
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/p4.png", plot = p4)
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/p5.png", plot = p5)
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/p6.png", plot = p6)
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/p7.png", plot = p7)
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/p8.png", plot = p8)
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/p9.png", plot = p9)
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/p10.png", plot = p10)
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/p11.png", plot = p11)
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/p12.png", plot = p12)
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/p13.png", plot = p13)
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/p14.png", plot = p14)
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/p15.png", plot = p15)
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/p16.png", plot = p16)
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/p17.png", plot = p17)
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/p18.png", plot = p18)


#grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, ncol=6, nrow =3)


#### what is this 'piggyback' cloud? - shape_cell24_TC.control
f <- function(x) {x + 9}
ggplot(df, aes(x = log2(shape_cell24_TC.control), y = log2(cell24))) + geom_point() + stat_function(fun = f, colour = "red")

piggyback <-df$n[log2(df$cell24) > (log2(df$shape_cell24_TC.control) + 9)]
piggyback <- piggyback[!is.na(piggyback)]
piggyback <- data.frame(n = piggyback, piggyback = rep("piggyback", length(piggyback)))
df <- merge(df, piggyback, by="n", all=TRUE)


ggplot(df, aes(x = log2(shape_cell24_TC.control), y = log2(cell24))) + geom_point() + stat_function(fun = f, colour = "red") +
  geom_point(data = subset(df, piggyback == 'piggyback'), aes(x=log2(shape_cell24_TC.control), y=log2(cell24), colour=piggyback)) +
  theme(legend.position="none")
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/piggyback.png")


##### add GO
go <- read.csv(file="/Users/kasia/Documents/PhD/matrix_go_fq.csv")
go <- data.frame(n = go$X.transcript_id, go_term = go$GO_Term_Name)

df <- merge(df, go, by="n", all=TRUE)
go_piggyback <- subset(df, piggyback == "piggyback")$go_term # mostly "binding" ??? 'membrane', 'nucleus', 'DNA binding' --- what is the cellular component?? ---- for CONTROL

## THEORY: are these nuclear binding proteins? SHAPE read coverage is lower relative to RNA-seq, can the reagents penetrate to nucleus?

#### what is this 'piggyback' cloud? - shape_cell24_TC.treated
f <- function(x) {x + 10}
ggplot(df, aes(x = log2(shape_cell24_TC.treated), y = log2(cell24))) + geom_point() + stat_function(fun = f, colour = "red")

piggyback <-df$n[log2(df$cell24) > (log2(df$shape_cell24_TC.treated) + 10)]
piggyback <- piggyback[!is.na(piggyback)]
piggyback <- data.frame(n = piggyback, piggyback = rep("piggyback", length(piggyback)))
df <- merge(df, piggyback, by="n", all=TRUE)


ggplot(df, aes(x = log2(shape_cell24_TC.treated), y = log2(cell24))) + geom_point() + stat_function(fun = f, colour = "red") +
  geom_point(data = subset(df, piggyback.y == 'piggyback'), aes(x=log2(shape_cell24_TC.treated), y=log2(cell24), colour=piggyback.y)) +
  theme(legend.position="none")
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/general/piggyback_NAI.png")

go_piggyback <- subset(df, piggyback.y == "piggyback")$go_term
go_regular <- subset(df, is.na(piggyback.y))$go_te

gop <- unlist(sapply(go_piggyback, function(x){strsplit(as.character(x), split=";")}))
gor <- unlist(sapply(go_regular, function(x){strsplit(as.character(x), split=";")}))

########
# summary(as.factor(gop))/length(gop)
# -                                             protein_binding 
# 0.064338362                                                 0.050079217 
# nucleus                                        transferase_activity 
# 0.022869739                                                 0.022594200 
# ATP_binding                                          cellular_component 
# 0.021492044                                                 0.021078735 
# membrane                                          nucleotide_binding 
# 0.019563271                                                 0.019563271 
# integral_component_of_membrane                                           metal_ion_binding 
# 0.016738996                                                 0.016601226 
# nucleic_acid_binding                                          biological_process 
# 0.015912379                                                 0.015705724 
# zinc_ion_binding                                          molecular_function 
# 0.015361301                                                 0.015223531 
# DNA_binding                                                   cytoplasm 
# 0.015085762                                                 0.013019219 

#####
# summary(as.factor(gor))/length(gor)
# -                                                 protein_binding 
# 0.132486686                                                     0.051887022 
# membrane                                  integral_component_of_membrane 
# 0.028895022                                                     0.025514919 
# metal_ion_binding                                            transferase_activity 
# 0.022001070                                                     0.018596649 
# cellular_component                                                     ATP_binding 
# 0.017271357                                                     0.015812319 
# nucleus                                              biological_process 
# 0.014766676                                                     0.012985434 
# zinc_ion_binding                                              molecular_function 
# 0.011927632                                                     0.011824283 
# DNA_binding                                              nucleotide_binding 
# 0.011258906                                                     0.011088685 
# nucleic_acid_binding                                             signal_transduction 
# 0.011033971                                                     0.009884979 