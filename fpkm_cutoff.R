library(GenomicRanges)
library(GenomicFeatures)

fpkm_2_4 <- read.csv("/Volumes/USELESS/META/SHAPES/fpkm/Pauli_rna/0_2to4Cell_fpkm.csv")
fpkm_256 <- read.csv("/Volumes/USELESS/META/SHAPES/fpkm/Pauli_rna/1_256Cell_fpkm.csv")

#ggplot(df, aes(x = rpkm)) + geom_density(alpha=.3)

#  loose cut-off at .5 RPKM, or 1, to be more stringent
# most of the papers arbitrarily define expression threshold i.e, >1 FPKM/RPKM to identify an expressed transcript

df <- data.frame(fpkm_2_4$rkpm)
colnames(df) <- c("rpkm")
rownames(df) <- fpkm_2_4$tx_id
rpkm_2_4 <- subset(df, rpkm > 1) ### rownames(rpkm_2_4)

df <- data.frame(fpkm_256$rkpm)
colnames(df) <- c("rpkm")
rownames(df) <- fpkm_256$tx_id
rpkm_256 <- subset(df, rpkm > 1) ### rownames(rpkm_256)


### TE
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

###
txdb_can <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/danio_rerio/zebrafish_GRCz10.gtf", format = "gtf")
cds <- cdsBy(txdb_can, by="tx", use.names=TRUE)
fiveUTR <- fiveUTRsByTranscript(txdb_can, use.names=TRUE)
threeUTR <- threeUTRsByTranscript(txdb_can, use.names=TRUE)

load("/Volumes/USELESS/META/SHAPES/gc_unsmoothed_2_4.Rsave")
load("/Volumes/USELESS/META/SHAPES/gc_unsmoothed_256.Rsave")

gc_unsmoothed_2_4$dtcr[is.na(gc_unsmoothed_2_4$dtcr)] <- 0
gc_unsmoothed_256$dtcr[is.na(gc_unsmoothed_256$dtcr)] <- 0

### CDS
# 2-4cell
cds_2_4 <- subsetByOverlaps(gc_unsmoothed_2_4, cds)
cds_2_4[is.na(cds_2_4$dtcr)]$dtcr <- 0
cds_2_4 <- split(cds_2_4, cds_2_4$trnames)
# 256cell
cds_256 <- subsetByOverlaps(gc_unsmoothed_256, cds)
cds_256[is.na(cds_256$dtcr)]$dtcr <- 0
cds_256 <- split(cds_256, cds_256$trnames)
# names
names(cds_2_4) <- sapply(names(cds_2_4), function(x){substr(x,1,18)})
names(cds_256) <- sapply(names(cds_256), function(x){substr(x,1,18)})
### FILTER ON RPKM
cds_2_4 <- cds_2_4[names(cds_2_4) %in% rownames(rpkm_2_4)]
cds_256 <- cds_256[names(cds_256) %in% rownames(rpkm_256)]
# sum(dtcr) > 0
cds_2_4 <- cds_2_4[sapply(cds_2_4, function(x){sum(x$dtcr) > 0})]
cds_256 <- cds_256[sapply(cds_256, function(x){sum(x$dtcr) > 0})]
cds_2_4 <- cds_2_4[sapply(cds_2_4, function(x){length(x) > 0})]
cds_256 <- cds_256[sapply(cds_256, function(x){length(x) > 0})]

##############
### 5UTR
# 2-4cell
five_2_4 <- subsetByOverlaps(gc_unsmoothed_2_4, fiveUTR)
five_2_4[is.na(five_2_4$dtcr)]$dtcr <- 0
five_2_4 <- split(five_2_4, five_2_4$trnames)
# 256cell
five_256 <- subsetByOverlaps(gc_unsmoothed_256, fiveUTR)
five_256[is.na(five_256$dtcr)]$dtcr <- 0
five_256 <- split(five_256, five_256$trnames)
# names
names(five_2_4) <- sapply(names(five_2_4), function(x){substr(x,1,18)})
names(five_256) <- sapply(names(five_256), function(x){substr(x,1,18)})
### FILTER ON RPKM
five_2_4 <- five_2_4[names(five_2_4) %in% rownames(rpkm_2_4)]
five_256 <- five_256[names(five_256) %in% rownames(rpkm_256)]
# sum(dtcr) > 0
#five_2_4 <- five_2_4[sapply(five_2_4, function(x){sum(x$dtcr) > 0})]
#five_256 <- five_256[sapply(five_256, function(x){sum(x$dtcr) > 0})]
five_2_4 <- five_2_4[sapply(five_2_4, function(x){length(x) > 0})]
five_256 <- five_256[sapply(five_256, function(x){length(x) > 0})]

##############
### 3UTR
# 2-4cell
three_2_4 <- subsetByOverlaps(gc_unsmoothed_2_4, threeUTR)
three_2_4[is.na(three_2_4$dtcr)]$dtcr <- 0
three_2_4 <- split(three_2_4, three_2_4$trnames)
# 256cell
three_256 <- subsetByOverlaps(gc_unsmoothed_256, threeUTR)
three_256[is.na(three_256$dtcr)]$dtcr <- 0
three_256 <- split(three_256, three_256$trnames)
# names
names(three_2_4) <- sapply(names(three_2_4), function(x){substr(x,1,18)})
names(three_256) <- sapply(names(three_256), function(x){substr(x,1,18)})
### FILTER ON RPKM
three_2_4 <- three_2_4[names(three_2_4) %in% rownames(rpkm_2_4)]
three_256 <- three_256[names(three_256) %in% rownames(rpkm_256)]
# sum(dtcr) > 0
#three_2_4 <- three_2_4[sapply(three_2_4, function(x){sum(x$dtcr) > 0})]
#three_256 <- three_256[sapply(three_256, function(x){sum(x$dtcr) > 0})]
three_2_4 <- three_2_4[sapply(three_2_4, function(x){length(x) > 0})]
three_256 <- three_256[sapply(three_256, function(x){length(x) > 0})]

### get average dtcr values (leader, CDS, trailer); combine with TE
# five_2_4, cds_2_4, three_2_4; data_2_4

# 2-4cell
dens_five_2_4 <- sapply(five_2_4, function(x){sum(x$dtcr)/length(x)})
dens_cds_2_4 <- sapply(cds_2_4, function(x){sum(x$dtcr)/length(x)})
dens_three_2_4 <- sapply(three_2_4, function(x){sum(x$dtcr)/length(x)})

# 256cell
dens_five_256 <- sapply(five_256, function(x){sum(x$dtcr)/length(x)})
dens_cds_256 <- sapply(cds_256, function(x){sum(x$dtcr)/length(x)})
dens_three_256 <- sapply(three_256, function(x){sum(x$dtcr)/length(x)})

##################
########### SAVED
load("/Volumes/USELESS/META/SHAPES/shape_te_cds_2_4.Rsave")
load("/Volumes/USELESS/META/SHAPES/shape_te_cds_256.Rsave")

shape_2_4 <- shape_te_cds_2_4[rownames(shape_te_cds_2_4) %in% rownames(rpkm_2_4),]
shape_256 <- shape_te_cds_256[rownames(shape_te_cds_256) %in% rownames(rpkm_256),]

## get TE higher than 0
shape_2_4 <- subset(shape_2_4, te_cds > 0)
shape_256 <- subset(shape_256, te_cds > 0)

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
  return(p)
}

p1a <- plot_shape(shape_2_4, top_q = 0.9, bot_q = 0.1) + ggtitle("2-4cell, 10%")
ggsave("/Volumes/USELESS/META/SHAPES/binned_TE/2_4cell_1pr.png", plot = p1a)

p1b <- plot_shape(shape_256, top_q = 0.9, bot_q = 0.1) + ggtitle("256cell, 10%")
ggsave("/Volumes/USELESS/META/SHAPES/binned_TE/256cell_1pr.png", plot = p1b)
