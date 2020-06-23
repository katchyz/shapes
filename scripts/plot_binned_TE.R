data = read.csv("/Users/kasia/Documents/PhD/TE.csv", header = TRUE)
rownames(data) <- data$X.transcript_id

data_2_4 <- data.frame(data$X01_2to4cell_gene_te)
rownames(data_2_4) <- rownames(data)
colnames(data_2_4) <- c("TE")
data_2_4 <- subset(data_2_4, data_2_4$TE > 0)

data_256 <- data.frame(data$X02_256cell_gene_te)
rownames(data_256) <- rownames(data)
colnames(data_256) <- c("TE")
data_256 <- subset(data_256, data_256$TE > 0)


## df <-df[order(df$TE), , drop = FALSE]
# quantiles
q_bottom <- quantile(data_2_4$TE, probs = c(0.33))
q_top <- quantile(data_2_4$TE, probs = c(0.67))
data_2_4_low <- subset(data_2_4, data_2_4$TE < q_bottom)
data_2_4_high <- subset(data_2_4, data_2_4$TE > q_top)

q_bottom <- quantile(data_256$TE, probs = c(0.33))
q_top <- quantile(data_256$TE, probs = c(0.67))
data_256_low <- subset(data_256, data_256$TE < q_bottom)
data_256_high <- subset(data_256, data_256$TE > q_top)

# shape_te_cds_2_4
# shape_te_cds_256

########
five_2_4 <- subsetByOverlaps(gc_dtcr_2_4, fiveUTR)
five_2_4_list <- split(five_2_4, five_2_4$trnames)
dens_five_2_4 <- sapply(five_2_4_list, function(x){sum(x$dtcr)/length(x)})

three_2_4 <- subsetByOverlaps(gc_dtcr_2_4, threeUTR)
three_2_4_list <- split(three_2_4, three_2_4$trnames)
dens_three_2_4 <- sapply(three_2_4_list, function(x){sum(x$dtcr)/length(x)})

five_256 <- subsetByOverlaps(gc_dtcr_256, fiveUTR)
five_256_list <- split(five_256, five_256$trnames)
dens_five_256 <- sapply(five_256_list, function(x){sum(x$dtcr)/length(x)})

three_256 <- subsetByOverlaps(gc_dtcr_256, threeUTR)
three_256_list <- split(three_256, three_256$trnames)
dens_three_256 <- sapply(three_256_list, function(x){sum(x$dtcr)/length(x)})

##############

# shape_te_cds_2_4
# shape_te_cds_256
names(dens_five_2_4) <- sapply(names(dens_five_2_4), function(x){substr(x,1,18)})
names(dens_three_2_4) <- sapply(names(dens_three_2_4), function(x){substr(x,1,18)})
names(dens_five_256) <- sapply(names(dens_five_256), function(x){substr(x,1,18)})
names(dens_three_256) <- sapply(names(dens_three_256), function(x){substr(x,1,18)})

# ...add that stuff to shape_te_cds_2_4
shape_te_cds_2_4$shape_five <- as.numeric(dens_five_2_4[names(dens_five_2_4) %in% rownames(shape_te_cds_2_4)])
shape_te_cds_2_4$shape_three <- as.numeric(dens_three_2_4[names(dens_three_2_4) %in% rownames(shape_te_cds_2_4)])

shape_te_cds_256$shape_five <- as.numeric(dens_five_256[names(dens_five_256) %in% rownames(shape_te_cds_256)])
shape_te_cds_256$shape_three <- as.numeric(dens_three_256[names(dens_three_256) %in% rownames(shape_te_cds_256)])

save(shape_te_cds_2_4, file="/Volumes/USELESS/META/SHAPES/shape_te_cds_2_4.Rsave")
save(shape_te_cds_256, file="/Volumes/USELESS/META/SHAPES/shape_te_cds_256.Rsave")

####### plot
#ggplot(shape_te_cds_2_4, aes(x=te, y=shape_three, fill=te)) + geom_boxplot() + guides(fill=FALSE) + scale_y_log10()

shape_te_cds_2_4 <- subset(shape_te_cds_2_4, te_cds > 0)
shape_te_cds_256 <- subset(shape_te_cds_256, te_cds > 0)


plot_shape <- function(df, top_q=0.95, bot_q=0.05, shape_feature) {
  df = within(df, {
    te = ifelse(rownames(df) %in% rownames(subset(df, df$te_cds > quantile(df$te_cds, probs = c(top_q)))), "H", "M")
  })
  df = within(df, {
    te = ifelse(rownames(df) %in% rownames(subset(df, df$te_cds < quantile(df$te_cds, probs = c(bot_q)))), "L", df$te)
  })
  p <- ggplot(temp, aes(x=te, y=eval(parse(text = shape_feature)), fill=te)) + geom_boxplot() + guides(fill=FALSE) + scale_y_log10() + ylim(0.001, 0.010) + ylab(shape_feature)
  return(p)
}


ggplot(d, aes(x = variable, y = value, fill = variable)) + geom_boxplot() + facet_wrap(~ variable)


















plot_shape <- function(df, top_q=0.95, bot_q=0.05) {
  df = within(df, {
    te = ifelse(rownames(df) %in% rownames(subset(df, df$te_cds > quantile(df$te_cds, probs = c(top_q)))), "H", "M")
  })
  df = within(df, {
    te = ifelse(rownames(df) %in% rownames(subset(df, df$te_cds < quantile(df$te_cds, probs = c(bot_q)))), "L", df$te)
  })
  d <- melt(df[0:4], id.var="te_cds")
  d$te <- rep(df$te, 3)
  p <- ggplot(d, aes(x = te, y = value, fill = te)) + geom_boxplot(outlier.size = NA) + facet_wrap(~ variable, scales = "free_y") + scale_y_log10() + ylab("avg shape reactivities") + guides(fill=FALSE) #+ ylim(0.001, 0.010)
  p_cds <- ggplot(df, aes(x = te, y = shape_cds, fill = te)) + geom_boxplot(outlier.size = NA) + scale_y_log10() + ylab("avg shape reactivities") + guides(fill=FALSE) #+ ylim(0.001, 0.010)
  return(p)
}

p1a <- plot_shape(shape_te_cds_2_4, top_q = 0.99, bot_q = 0.01) + ggtitle("2-4cell, 1%")
ggsave("/Volumes/USELESS/META/SHAPES/binned_TE/2_4cell_1pr.png", plot = p1a)

p1b <- plot_shape(shape_te_cds_256, top_q = 0.99, bot_q = 0.01) + ggtitle("256cell, 1%")
ggsave("/Volumes/USELESS/META/SHAPES/binned_TE/256cell_1pr.png", plot = p1b)

p25a <- plot_shape(shape_te_cds_2_4, top_q = 0.75, bot_q = 0.25) + ggtitle("2-4cell, 25%")
ggsave("/Volumes/USELESS/META/SHAPES/binned_TE/2_4cell_25pr.png", plot = p25a)

p25b <- plot_shape(shape_te_cds_256, top_q = 0.75, bot_q = 0.25) + ggtitle("256cell, 25%")
ggsave("/Volumes/USELESS/META/SHAPES/binned_TE/256cell_25pr.png", plot = p25b)


