cai <- read.csv("/Volumes/USELESS/META/SHAPES/caiDf.csv")
cai$X <- NULL
rownames(cai) <- cai$X0
colnames(cai) <- c('transcript', 'CAI')

shape_2_4 <- shape_te_cds_2_4[rownames(shape_te_cds_2_4) %in% rownames(cai),]
shape_256 <- shape_te_cds_256[rownames(shape_te_cds_256) %in% rownames(cai),]

cai_2_4 <- cai[rownames(cai) %in% rownames(shape_2_4),]
cai_256 <- cai[rownames(cai) %in% rownames(shape_256),]

shape_2_4$cai <- cai_2_4$CAI
shape_256$cai <- cai_256$CAI

p_CAI_2_4 <- ggplot(shape_2_4, aes(x=shape_cds, y=cai)) + geom_point() + xlab("sum of reactivities/length") + ylab("CAI") + scale_x_log10() + scale_y_log10() + ggtitle("2_4cell")

p_CAI_256 <- ggplot(shape_256, aes(x=shape_cds, y=cai)) + geom_point() + xlab("sum of reactivities/length") + ylab("CAI") + scale_x_log10() + scale_y_log10() + ggtitle("256cell")

ggsave("/Volumes/USELESS/META/SHAPES/CAI/p_CAI_2_4.png", plot = p_CAI_2_4)
ggsave("/Volumes/USELESS/META/SHAPES/CAI/p_CAI_256.png", plot = p_CAI_256)

##########

p_CAI_TE_2_4 <- ggplot(shape_2_4, aes(x=te_cds, y=cai)) + geom_point() + xlab("TE") + ylab("CAI") + scale_x_log10() + scale_y_log10() + ggtitle("2_4cell")

p_CAI_TE_256 <- ggplot(shape_256, aes(x=te_cds, y=cai)) + geom_point() + xlab("TE") + ylab("CAI") + scale_x_log10() + scale_y_log10() + ggtitle("256cell")

ggsave("/Volumes/USELESS/META/SHAPES/CAI/p_CAI_TE_2_4.png", plot = p_CAI_TE_2_4)
ggsave("/Volumes/USELESS/META/SHAPES/CAI/p_CAI_TE_256.png", plot = p_CAI_TE_256)

shape_2_4 <- shape_2_4[shape_2_4$te_cds > 0,]
shape_256 <- shape_256[shape_256$te_cds > 0,]

##########
plot_shape_tb <- function(df, top_q=0.95, bot_q=0.05) {
  df = within(df, {
    bin = ifelse(rownames(df) %in% rownames(subset(df, df$te_cds > quantile(df$te_cds, probs = c(top_q)))), "H", "M")
  })
  df = within(df, {
    bin = ifelse(rownames(df) %in% rownames(subset(df, df$te_cds < quantile(df$te_cds, probs = c(bot_q)))), "L", df$bin)
  })
  p <- ggplot(df, aes(x = shape_cds, y = cai, fill = bin)) + geom_point() + facet_wrap(~ bin) + xlab("shape") + ylab("CAI") + guides(fill=FALSE) + scale_x_log10() + scale_y_log10()
  return(p)
}

p_2_4 <- plot_shape_tb(shape_2_4, top_q=0.95, bot_q=0.05) + ggtitle("2_4cell")
p_256 <- plot_shape_tb(shape_256, top_q=0.95, bot_q=0.05) + ggtitle("256cell")

ggsave("/Volumes/USELESS/META/SHAPES/CAI/p_str_CAI_by_TE_2_4.png", plot = p_2_4)
ggsave("/Volumes/USELESS/META/SHAPES/CAI/p_str_CAI_by_TE_256.png", plot = p_256)



#######
#######
### last 50

cai_last50 <- read.csv("/Volumes/USELESS/META/SHAPES/caiDf_last50.csv")
cai_last50$X <- NULL
rownames(cai_last50) <- cai_last50$X0
colnames(cai_last50) <- c('transcript', 'CAI')

cds2_4 <- cds_2_4_list[lengths(cds_2_4_list) > 150]
cds2_4 <- cds2_4[sapply(cds2_4, function(x){sum(x$dtcr) > 0})]
names(cds2_4) <- sapply(names(cds2_4), function(x){substr(x,1,18)})

cds256 <- cds_256_list[lengths(cds_256_list) > 150]
cds256 <- cds256[sapply(cds256, function(x){sum(x$dtcr) > 0})]
names(cds256) <- sapply(names(cds256), function(x){substr(x,1,18)})

shape_last50_2_4 <- sapply(cds2_4, function(x){mean(x[(length(x)-150):length(x)]$dtcr)})
shape_last50_256 <- sapply(cds256, function(x){mean(x[(length(x)-150):length(x)]$dtcr)})

cai2_4 <- cai_last50[rownames(cai_last50) %in% names(shape_last50_2_4),]
sh2_4 <- shape_last50_2_4[names(shape_last50_2_4) %in% rownames(cai2_4)]
sh_cai_2_4 <- data.frame(unname(sh2_4), cai2_4$CAI)
colnames(sh_cai_2_4) <- c("shape", "CAI")
rownames(sh_cai_2_4) <- rownames(cai2_4)
ggplot(sh_cai_2_4, mapping = aes(x = shape, y = CAI)) + geom_point() + scale_x_log10() + scale_y_log10() + ggtitle("2-4cell, last 50 codons of CDS")

cai256 <- cai_last50[rownames(cai_last50) %in% names(shape_last50_256),]
sh256 <- shape_last50_256[names(shape_last50_256) %in% rownames(cai256)]
sh_cai_256 <- data.frame(unname(sh256), cai256$CAI)
colnames(sh_cai_256) <- c("shape", "CAI")
rownames(sh_cai_256) <- rownames(cai256)
ggplot(sh_cai_256, mapping = aes(x = shape, y = CAI)) + geom_point() + scale_x_log10() + scale_y_log10() + ggtitle("256cell, last 50 codons of CDS")
