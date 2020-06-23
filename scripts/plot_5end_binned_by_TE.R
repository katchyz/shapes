#unsmoothed
load("/Volumes/USELESS/META/SHAPES/gc_unsmoothed_2_4.Rsave")
load("/Volumes/USELESS/META/SHAPES/gc_unsmoothed_256.Rsave")

gc_unsmoothed_2_4$dtcr[is.na(gc_unsmoothed_2_4$dtcr)] <- 0
gc_unsmoothed_256$dtcr[is.na(gc_unsmoothed_256$dtcr)] <- 0

exon_2_4_list <- split(gc_unsmoothed_2_4, gc_unsmoothed_2_4$trnames)
exon_2_4_list <- exon_2_4_list[lengths(exon_2_4_list) > 50]
names(exon_2_4_list) <- sapply(names(exon_2_4_list), function(x){substr(x,1,18)})

exon_256_list <- split(gc_unsmoothed_256, gc_unsmoothed_256$trnames)
exon_256_list <- exon_256_list[lengths(exon_256_list) > 50]
names(exon_256_list) <- sapply(names(exon_256_list), function(x){substr(x,1,18)})

load("/Volumes/USELESS/META/SHAPES/shape_te_cds_2_4.Rsave")
load("/Volumes/USELESS/META/SHAPES/shape_te_cds_256.Rsave")

# plot first 50nt stratified by TE
library(dplyr)
library(reshape2)
library(plyr)
shape_2_4 <- shape_te_cds_2_4[shape_te_cds_2_4$te_cds > 0,]
shape_256 <- shape_te_cds_256[shape_te_cds_256$te_cds > 0,]

exon_2_4_list <- exon_2_4_list[lengths(exon_2_4_list) > 50]
exon_256_list <- exon_256_list[lengths(exon_256_list) > 50]

exon_2_4_list <- exon_2_4_list[names(exon_2_4_list) %in% rownames(shape_2_4)]
exon_256_list <- exon_256_list[names(exon_256_list) %in% rownames(shape_256)]

exon_2_4_list <- exon_2_4_list[sapply(exon_2_4_list, function(x){sum(x[1:50]$dtcr) > 0})]
exon_256_list <- exon_256_list[sapply(exon_256_list, function(x){sum(x[1:50]$dtcr) > 0})]

shape_2_4 <- shape_2_4[rownames(shape_2_4) %in% names(exon_2_4_list),]
shape_256 <- shape_256[rownames(shape_256) %in% names(exon_256_list),]


### each NORMALIZED to 1
plot_shape <- function(df, exon_list, nbins = 5) {
  df$bin <- ntile(df$te_cds, nbins)
  names(exon_list) <- sapply(names(exon_list), function(x){substr(x,1,18)})
  exon_list <- exon_list[names(exon_list) %in% rownames(df)]
  ### divide based on bin
  # get fragments of 50, reduce based on bin
  df <- rbind(t(df), sapply(exon_list, function(x){as.numeric((x[1:50]$dtcr)/sum(x[1:50]$dtcr))}))
  df <- data.frame(t(df))
  df$bin <- as.factor(df$bin)
  meta <- by(df, df$bin, function(x){as.numeric(colSums(x[,6:55]))})
  print(meta)
  df <- ldply(meta, data.frame)
  colnames(df) <- c("bin", "reactivities")
  df$scale <- rep(c(1:50),nbins)
  p <- ggplot(df, aes(x = scale, y = reactivities, fill = bin)) + geom_bar(stat = "identity") + facet_wrap(~ bin) + xlab("position from 5'end of transcript") + guides(fill=FALSE)
  return(p)
}

### NOT NORMALIZED (to see how extreme it is)
plot_shape <- function(df, exon_list, nbins = 5) {
  df$bin <- ntile(df$te_cds, nbins)
  names(exon_list) <- sapply(names(exon_list), function(x){substr(x,1,18)})
  exon_list <- exon_list[names(exon_list) %in% rownames(df)]
  ### divide based on bin
  # get fragments of 50, reduce based on bin
  df <- rbind(t(df), sapply(exon_list, function(x){as.numeric((x[1:50]$dtcr))}))
  df <- data.frame(t(df))
  df$bin <- as.factor(df$bin)
  meta <- by(df, df$bin, function(x){as.numeric(colSums(x[,6:55]))})
  print(meta)
  df <- ldply(meta, data.frame)
  colnames(df) <- c("bin", "reactivities")
  df$scale <- rep(c(1:50),nbins)
  p <- ggplot(df, aes(x = scale, y = reactivities, fill = bin)) + geom_bar(stat = "identity") + facet_wrap(~ bin) + xlab("position from 5'end of transcript") + guides(fill=FALSE)
  return(p)
}

p1 <- plot_shape(shape_2_4, exon_2_4_list)
p2 <- plot_shape(shape_256, exon_256_list)
p1 <- p1 + xlab("position from 5'end of transcript") + ggtitle("2-4cell")
ggsave("/Volumes/USELESS/META/SHAPES/binned_TE/2_4cell_5end_unsm.png", plot = p1)
p2 <- p2 + xlab("position from 5'end of transcript") + ggtitle("256cell")
ggsave("/Volumes/USELESS/META/SHAPES/binned_TE/256cell_5end_unsm.png", plot = p2)

p1 <- p1 + xlab("position from 5'end of transcript") + ggtitle("2-4cell (not normalized)")
ggsave("/Volumes/USELESS/META/SHAPES/binned_TE/2_4cell_5end_unsm_unn.png", plot = p1)
p2 <- p2 + xlab("position from 5'end of transcript") + ggtitle("256cell (not normalized)")
ggsave("/Volumes/USELESS/META/SHAPES/binned_TE/256cell_5end_unsm_unn.png", plot = p2)


rRNA_2_4 <- shape_2_4[rownames(shape_2_4) %in% rownames(rDf),]
rRNA_256 <- shape_256[rownames(shape_256) %in% rownames(rDf),]
p_256_rRNA <- plot_shape(rRNA_256, exon_256_list)
p_2_4_rRNA <- plot_shape(rRNA_2_4, exon_2_4_list)

snoRNA_256 <- shape_256[rownames(shape_256) %in% rownames(snoDf),]
p_256_snoRNA <- plot_shape(snoRNA_256, exon_256_list)
##########
plot_shape_tb <- function(df, exon_list, top_q=0.95, bot_q=0.05) {
  df = within(df, {
    bin = ifelse(rownames(df) %in% rownames(subset(df, df$te_cds > quantile(df$te_cds, probs = c(top_q)))), "H", "M")
  })
  df = within(df, {
    bin = ifelse(rownames(df) %in% rownames(subset(df, df$te_cds < quantile(df$te_cds, probs = c(bot_q)))), "L", df$bin)
  })
  names(exon_list) <- sapply(names(exon_list), function(x){substr(x,1,18)})
  exon_list <- exon_list[names(exon_list) %in% rownames(df)]
  df <- rbind(t(df), sapply(exon_list, function(x){(x[1:50]$dtcr)/sum(x[1:50]$dtcr)}))
  df <- data.frame(t(df), stringsAsFactors = FALSE)
  #df$bin <- as.factor(df$bin)
  df <- cbind(df[,1:5], as.data.frame(sapply(df[,6:55], as.numeric)))
  meta <- by(df, df$bin, function(x){as.numeric(colSums(x[,6:55]))})
  print(meta)
  df <- ldply(meta, data.frame)
  colnames(df) <- c("bin", "reactivities")
  df$scale <- rep(c(1:50),3)
  p <- ggplot(df, aes(x = scale, y = reactivities, fill = bin)) + geom_bar(stat = "identity") + facet_wrap(~ bin) + xlab("position from 5'end of transcript") + guides(fill=FALSE)
  return(p)
}

p3 <- plot_shape_tb(shape_2_4, exon_2_4_list, top_q=0.95, bot_q=0.05)
p4 <- plot_shape_tb(shape_256, exon_256_list, top_q=0.95, bot_q=0.05)

p3 <- p3 + facet_wrap(~ bin, scales = "free_y") + ggtitle("2-4cell")
ggsave("/Volumes/USELESS/META/SHAPES/binned_TE/2_4cell_5end_5pr_unsm.png", plot = p3)
p4 <- p4 + facet_wrap(~ bin, scales = "free_y") + ggtitle("256cell")
ggsave("/Volumes/USELESS/META/SHAPES/binned_TE/256cell_5end_5pr_unsm.png", plot = p4)

###
p5 <- plot_shape_tb(shape_2_4, exon_2_4_list, top_q=0.95, bot_q=0.05)
p6 <- plot_shape_tb(shape_256, exon_256_list, top_q=0.95, bot_q=0.05)

p5 <- p5 + facet_wrap(~ bin, scales = "free_y") + ggtitle("2-4cell")
ggsave("/Volumes/USELESS/META/SHAPES/binned_TE/2_4cell_5end_5pr.png", plot = p5)
p6 <- p6 + facet_wrap(~ bin, scales = "free_y") + ggtitle("256cell")
ggsave("/Volumes/USELESS/META/SHAPES/binned_TE/256cell_5end_5pr.png", plot = p6)

