load("/Volumes/USELESS/META/SHAPES/gc_control_comp_2_4.Rsave")
load("/Volumes/USELESS/META/SHAPES/gc_control_comp_256.Rsave")

load("/Volumes/USELESS/META/SHAPES/gc_treated_comp_2_4.Rsave")
load("/Volumes/USELESS/META/SHAPES/gc_treated_comp_256.Rsave")

gc_control_comp_2_4$TC[is.na(gc_control_comp_2_4$TC)] <- 0
gc_control_comp_256$TC[is.na(gc_control_comp_256$TC)] <- 0
gc_treated_comp_2_4$TC[is.na(gc_treated_comp_2_4$TC)] <- 0
gc_treated_comp_256$TC[is.na(gc_treated_comp_256$TC)] <- 0

gc_control_comp_2_4$TC[is.nan(gc_control_comp_2_4$TC)] <- 0
gc_control_comp_256$TC[is.nan(gc_control_comp_256$TC)] <- 0
gc_treated_comp_2_4$TC[is.nan(gc_treated_comp_2_4$TC)] <- 0
gc_treated_comp_256$TC[is.nan(gc_treated_comp_256$TC)] <- 0

# ctrl_2_4
ctrl_2_4 <- split(gc_control_comp_2_4, gc_control_comp_2_4$trnames)
names(ctrl_2_4) <- sapply(names(ctrl_2_4), function(x){substr(x,1,18)})
ctrl_2_4 <- ctrl_2_4[lengths(ctrl_2_4) > 50]

# ctrl_256
ctrl_256 <- split(gc_control_comp_256, gc_control_comp_256$trnames)
names(ctrl_256) <- sapply(names(ctrl_256), function(x){substr(x,1,18)})
ctrl_256 <- ctrl_256[lengths(ctrl_256) > 50]

# trea_2_4
trea_2_4 <- split(gc_treated_comp_2_4, gc_treated_comp_2_4$trnames)
names(trea_2_4) <- sapply(names(trea_2_4), function(x){substr(x,1,18)})
trea_2_4 <- trea_2_4[lengths(trea_2_4) > 50]

# trea_256
trea_256 <- split(gc_treated_comp_256, gc_treated_comp_256$trnames)
names(trea_256) <- sapply(names(trea_256), function(x){substr(x,1,18)})
trea_256 <- trea_256[lengths(trea_256) > 50]


########### PLOT first 50
### each NORMALIZED to 1
plot_shape <- function(df, exon_list, nbins = 5) {
  df$bin <- ntile(df$te_cds, nbins)
  names(exon_list) <- sapply(names(exon_list), function(x){substr(x,1,18)})
  exon_list <- exon_list[names(exon_list) %in% rownames(df)]
  df <- df[rownames(df) %in% names(exon_list),]
  first50 <- sapply(exon_list, function(x){as.numeric((x[1:50]$TC)/sum(x[1:50]$TC))})
  first50[is.nan(first50)] <- 0
  df <- rbind(t(df), first50)
  df <- data.frame(t(df))
  df$bin <- as.factor(df$bin)
  meta <- by(df, df$bin, function(x){as.numeric(colSums(x[,7:56]))})
  print(meta)
  df <- ldply(meta, data.frame)
  colnames(df) <- c("bin", "reactivities")
  df$scale <- rep(c(1:50),nbins)
  p <- ggplot(df, aes(x = scale, y = reactivities, fill = bin)) + geom_bar(stat = "identity") + facet_wrap(~ bin) + xlab("position from 5'end of transcript") + guides(fill=FALSE)
  return(p)
}


pctrl_2_4 <- plot_shape(shape_2_4, ctrl_2_4)
pctrl_2_4 <- pctrl_2_4 + ggtitle("control, 2-4cell")
ggsave("/Volumes/USELESS/META/SHAPES/binned_TE/pctrl_2_4.png", plot = pctrl_2_4)

pctrl_256 <- plot_shape(shape_256, ctrl_256)
pctrl_256 <- pctrl_256 + ggtitle("control, 256cell")
ggsave("/Volumes/USELESS/META/SHAPES/binned_TE/pctrl_256.png", plot = pctrl_256)

ptrea_2_4 <- plot_shape(shape_2_4, trea_2_4)
ptrea_2_4 <- ptrea_2_4 + ggtitle("treated, 2-4cell")
ggsave("/Volumes/USELESS/META/SHAPES/binned_TE/ptrea_2_4.png", plot = ptrea_2_4)

ptrea_256 <- plot_shape(shape_256, trea_256)
ptrea_256 <- ptrea_256 + ggtitle("treated, 256cell")
ggsave("/Volumes/USELESS/META/SHAPES/binned_TE/ptrea_256.png", plot = ptrea_256)

##############
avg_ctr24 <- sapply(ctrl_2_4, function(x){mean(x$TC)})
avg_ctr256 <- sapply(ctrl_256, function(x){mean(x$TC)})
avg_tre24 <- sapply(trea_2_4, function(x){mean(x$TC)})
avg_tre256 <- sapply(trea_256, function(x){mean(x$TC)})

## cut-off on RPKM ?
fpkm24 <- FPKM_24[FPKM_24$exons_RNA_fpkm > 1,]
fpkm256 <- FPKM_256[FPKM_256$exons_RNA_fpkm > 1,]

## plot vs RNA fpkm
shape <- avg_ctr24[names(avg_ctr24) %in% rownames(fpkm24)]
rna <- fpkm24[rownames(fpkm24) %in% names(shape),]$exons_RNA_fpkm
df <- data.frame(shape = as.vector(shape), rna = rna)
ggplot(df, aes(x = shape, y = rna)) + geom_point() + scale_x_log10() + scale_y_log10() + ggtitle("2-4cell, control")
ggsave(file = "/Volumes/USELESS/META/SHAPES/accessibility/control_vs_RNA_24.png")

shape <- avg_ctr256[names(avg_ctr256) %in% rownames(fpkm256)]
rna <- fpkm256[rownames(fpkm256) %in% names(shape),]$exons_RNA_fpkm
df <- data.frame(shape = as.vector(shape), rna = rna)
ggplot(df, aes(x = shape, y = rna)) + geom_point() + scale_x_log10() + scale_y_log10() + ggtitle("256cell, control")
ggsave(file = "/Volumes/USELESS/META/SHAPES/accessibility/control_vs_RNA_256.png")

#
shape <- avg_tre24[names(avg_tre24) %in% rownames(fpkm24)]
rna <- fpkm24[rownames(fpkm24) %in% names(shape),]$exons_RNA_fpkm
df <- data.frame(shape = as.vector(shape), rna = rna)
ggplot(df, aes(x = shape, y = rna)) + geom_point() + scale_x_log10() + scale_y_log10() + ggtitle("2-4cell, treated")
ggsave(file = "/Volumes/USELESS/META/SHAPES/accessibility/treated_vs_RNA_24.png")

shape <- avg_tre256[names(avg_tre256) %in% rownames(fpkm256)]
rna <- fpkm256[rownames(fpkm256) %in% names(shape),]$exons_RNA_fpkm
df <- data.frame(shape = as.vector(shape), rna = rna)
ggplot(df, aes(x = shape, y = rna)) + geom_point() + scale_x_log10() + scale_y_log10() + ggtitle("256cell, treated")
ggsave(file = "/Volumes/USELESS/META/SHAPES/accessibility/treated_vs_RNA_256.png")


c <- ctrl_256[names(ctrl_256) %in% names(trea_256)]
t <- trea_256[names(trea_256) %in% names(c)]



