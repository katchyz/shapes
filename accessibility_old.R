## ACCESSIBILITY

library(data.table)
library(ggplot2)
library(GenomicRanges)
library(GenomicFeatures)
library(plyr)

### get transcript database
#txdb_can <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/danio_rerio/Danio_rerio.GRCz10.84.gtf", format = "gtf")
txdb_can <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/danio_rerio/zebrafish_GRCz10.gtf", format = "gtf")
cds <- cdsBy(txdb_can, by="tx", use.names=TRUE)
fiveUTR <- fiveUTRsByTranscript(txdb_can, use.names=TRUE)
threeUTR <- threeUTRsByTranscript(txdb_can, use.names=TRUE)
exons <- exonsBy(txdb_can, by="tx", use.names=TRUE)

names(cds) <- sapply(names(cds), function(x){substr(x,1,18)})
names(fiveUTR) <- sapply(names(fiveUTR), function(x){substr(x,1,18)})
names(threeUTR) <- sapply(names(threeUTR), function(x){substr(x,1,18)})
names(exons) <- sapply(names(exons), function(x){substr(x,1,18)})

txLengths <- transcriptLengths(txdb_can, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
rownames(txLengths) <- txLengths$tx_name
rownames(txLengths) <- sapply(rownames(txLengths), function(x){substr(x,1,18)})

### features
leader24 <- subsetByOverlaps(gc_unsmoothed_2_4, fiveUTR)
leader24 <- split(leader24, leader24$trnames)
names(leader24) <- sapply(names(leader24), function(x){substr(x,1,18)})

cds24 <- subsetByOverlaps(gc_unsmoothed_2_4, cds)
cds24 <- split(cds24, cds24$trnames)
names(cds24) <- sapply(names(cds24), function(x){substr(x,1,18)})

trailer24 <- subsetByOverlaps(gc_unsmoothed_2_4, threeUTR)
trailer24 <- split(trailer24, trailer24$trnames)
names(trailer24) <- sapply(names(trailer24), function(x){substr(x,1,18)})

exons24 <- subsetByOverlaps(gc_unsmoothed_2_4, exons)
exons24 <- split(exons24, exons24$trnames)
names(exons24) <- sapply(names(exons24), function(x){substr(x,1,18)})

leader256 <- subsetByOverlaps(gc_unsmoothed_256, fiveUTR)
leader256 <- split(leader256, leader256$trnames)
names(leader256) <- sapply(names(leader256), function(x){substr(x,1,18)})

cds256 <- subsetByOverlaps(gc_unsmoothed_256, cds)
cds256 <- split(cds256, cds256$trnames)
names(cds256) <- sapply(names(cds256), function(x){substr(x,1,18)})

trailer256 <- subsetByOverlaps(gc_unsmoothed_256, threeUTR)
trailer256 <- split(trailer256, trailer256$trnames)
names(trailer256) <- sapply(names(trailer256), function(x){substr(x,1,18)})

exons256 <- subsetByOverlaps(gc_unsmoothed_256, exons)
exons256 <- split(exons256, exons256$trnames)
names(exons256) <- sapply(names(exons256), function(x){substr(x,1,18)})

# FPKM_24
load(file = "/Volumes/USELESS/META/SHAPES/FPKM_24.Rdata")
# FPKM_256
load(file = "/Volumes/USELESS/META/SHAPES/FPKM_256.Rdata")

# leaders, CDSs, trailers (all present)
lct <- rownames(txLengths[(txLengths$utr5_len > 30) & (txLengths$cds_len > 100) & (txLengths$utr3_len > 30), ])
##
# cut-off on exons_RNA_fpkm
fpkm <- FPKM_24[rownames(FPKM_24) %in% lct,]
fpkm24 <- rownames(fpkm[fpkm$exons_RNA_fpkm > 1,])
fpkm <- FPKM_256[rownames(FPKM_256) %in% lct,]
fpkm256 <- rownames(fpkm[fpkm$exons_RNA_fpkm > 1,])

####
txdb_can <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/danio_rerio/Danio_rerio.GRCz10.84.gtf", format = "gtf")
# get longest transcript per gene ################### OTPG
txLengths <- transcriptLengths(txdb_can, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
rownames(txLengths) <- txLengths$tx_name
rownames(txLengths) <- sapply(rownames(txLengths), function(x){substr(x,1,18)})

txlen <- arrange(txLengths, gene_id, desc(tx_len))
#txlen$tx_name <- sapply(txlen$tx_name, function(x){substr(x,1,18)})
#txlen$gene_id <- sapply(txlen$gene_id, function(x){substr(x,1,18)})
txlen <- txlen[!duplicated(txlen$gene_id),]

fpkm24 <- fpkm24[fpkm24 %in% txlen$tx_name] ### OTPG
fpkm256 <- fpkm256[fpkm256 %in% txlen$tx_name] ### OTPG

## ???????????????
## SHAPE for fpkm24 & fpkm256
## whole tx, 5'utr, cds, 3'utr, ~start codon
## divide by exons_RNA_fpkm (optionally by feature's fpkm)

## exons by exons
ex24 <- exons24[names(exons24) %in% fpkm24]
exons_by_exons_24 <- sapply(names(ex24), function(x){sum(ex24[[x]]$dtcr) / FPKM_24[x,]$exons_RNA_fpkm})
ex256 <- exons256[names(exons256) %in% fpkm256]
exons_by_exons_256 <- sapply(names(ex256), function(x){sum(ex256[[x]]$dtcr) / FPKM_256[x,]$exons_RNA_fpkm})

## utr5 by utr5
l24 <- leader24[names(leader24) %in% fpkm24]
utr5_by_utr5_24 <- sapply(names(l24), function(x){sum(l24[[x]]$dtcr) / FPKM_24[x,]$utr5_RNA_fpkm})
l256 <- leader256[names(leader256) %in% fpkm256]
utr5_by_utr5_256 <- sapply(names(l256), function(x){sum(l256[[x]]$dtcr) / FPKM_256[x,]$utr5_RNA_fpkm})
## utr5 by exons
utr5_by_exons_24 <- sapply(names(l24), function(x){sum(l24[[x]]$dtcr) / FPKM_24[x,]$exons_RNA_fpkm})
utr5_by_exons_256 <- sapply(names(l256), function(x){sum(l256[[x]]$dtcr) / FPKM_256[x,]$exons_RNA_fpkm})

## cds by cds
c24 <- cds24[names(cds24) %in% fpkm24]
cds_by_cds_24 <- sapply(names(c24), function(x){sum(c24[[x]]$dtcr) / FPKM_24[x,]$cds_RNA_fpkm})
c256 <- cds256[names(cds256) %in% fpkm256]
cds_by_cds_256 <- sapply(names(c256), function(x){sum(c256[[x]]$dtcr) / FPKM_256[x,]$cds_RNA_fpkm})
## cds by exons
cds_by_exons_24 <- sapply(names(c24), function(x){sum(c24[[x]]$dtcr) / FPKM_24[x,]$exons_RNA_fpkm})
cds_by_exons_256 <- sapply(names(c256), function(x){sum(c256[[x]]$dtcr) / FPKM_256[x,]$exons_RNA_fpkm})

## utr3 by utr3
t24 <- trailer24[names(trailer24) %in% fpkm24]
utr3_by_utr3_24 <- sapply(names(t24), function(x){sum(t24[[x]]$dtcr) / FPKM_24[x,]$utr3_RNA_fpkm})
t256 <- trailer256[names(trailer256) %in% fpkm256]
utr3_by_utr3_256 <- sapply(names(t256), function(x){sum(t256[[x]]$dtcr) / FPKM_256[x,]$utr3_RNA_fpkm})
## utr3 by exons
utr3_by_exons_24 <- sapply(names(t24), function(x){sum(t24[[x]]$dtcr) / FPKM_24[x,]$exons_RNA_fpkm})
utr3_by_exons_256 <- sapply(names(t256), function(x){sum(t256[[x]]$dtcr) / FPKM_256[x,]$exons_RNA_fpkm})

## ~start by exons
s24 <- exons24[names(exons24) %in% fpkm24]
start_by_exons_24 <- sapply(names(s24), function(x){sum(s24[[x]]$dtcr[(txLengths[x,]$utr5_len-20):(txLengths[x,]$utr5_len)]) / FPKM_24[x,]$exons_RNA_fpkm})
s256 <- exons256[names(exons256) %in% fpkm256]
start_by_exons_256 <- sapply(names(s256), function(x){sum(s256[[x]]$dtcr[(txLengths[x,]$utr5_len-20):(txLengths[x,]$utr5_len)]) / FPKM_256[x,]$exons_RNA_fpkm})


########### vs TE
TE_24 <- data.frame(exons = FPKM_24$exons_Ribo_fpkm/FPKM_24$exons_RNA_fpkm, utr5 = FPKM_24$utr5_Ribo_fpkm/FPKM_24$utr5_RNA_fpkm, cds = FPKM_24$cds_Ribo_fpkm/FPKM_24$cds_RNA_fpkm, utr3 = FPKM_24$utr3_Ribo_fpkm/FPKM_24$utr3_RNA_fpkm, row.names = rownames(FPKM_24))

TE_256 <- data.frame(exons = FPKM_256$exons_Ribo_fpkm/FPKM_256$exons_RNA_fpkm, utr5 = FPKM_256$utr5_Ribo_fpkm/FPKM_256$utr5_RNA_fpkm, cds = FPKM_256$cds_Ribo_fpkm/FPKM_256$cds_RNA_fpkm, utr3 = FPKM_256$utr3_Ribo_fpkm/FPKM_256$utr3_RNA_fpkm, row.names = rownames(FPKM_256))

# exons_by_exons_24 exons_by_exons_256
df <- data.frame(accessibility = as.vector(exons_by_exons_24[order(names(exons_by_exons_24))]), TE = TE_24[rownames(TE_24) %in% names(exons_by_exons_24),]$cds)
ggplot(df, aes(x = accessibility, y = TE)) + geom_point() + scale_x_log10() + scale_y_log10() + ggtitle("2-4cell, whole TX") + ylab("TE, cds")
ggsave(file = "/Volumes/USELESS/META/SHAPES/accessibility/tx_24.png")

df <- data.frame(accessibility = as.vector(exons_by_exons_256[order(names(exons_by_exons_256))]), TE = TE_256[rownames(TE_256) %in% names(exons_by_exons_256),]$cds)
ggplot(df, aes(x = accessibility, y = TE)) + geom_point() + scale_x_log10() + scale_y_log10() + ggtitle("256cell, whole TX") + ylab("TE, cds")
ggsave(file = "/Volumes/USELESS/META/SHAPES/accessibility/tx_256.png")

# utr5_by_utr5_24, utr5_by_utr5_256, utr5_by_exons_24, utr5_by_exons_256
df <- data.frame(accessibility = as.vector(utr5_by_utr5_256[order(names(utr5_by_utr5_256))]), TE = TE_256[rownames(TE_256) %in% names(utr5_by_utr5_256),]$cds)
ggplot(df, aes(x = accessibility, y = TE)) + geom_point() + scale_x_log10() + scale_y_log10() + ggtitle("256cell, utr5") + ylab("TE, cds")
ggsave(file = "/Volumes/USELESS/META/SHAPES/accessibility/utr5_256.png")

# cds_by_cds_24, cds_by_cds_256, cds_by_exons_24, cds_by_exons_256
df <- data.frame(accessibility = as.vector(cds_by_cds_24[order(names(cds_by_cds_24))]), TE = TE_24[rownames(TE_24) %in% names(cds_by_cds_24),]$cds)
ggplot(df, aes(x = accessibility, y = TE)) + geom_point() + scale_x_log10() + scale_y_log10() + ggtitle("24cell, cds") + ylab("TE, cds")
ggsave(file = "/Volumes/USELESS/META/SHAPES/accessibility/cds_24.png")

# utr3_by_utr3_24, utr3_by_utr3_256, utr3_by_exons_24, utr3_by_exons_256
df <- data.frame(accessibility = as.vector(utr3_by_utr3_256[order(names(utr3_by_utr3_256))]), TE = TE_256[rownames(TE_256) %in% names(utr3_by_utr3_256),]$cds)
ggplot(df, aes(x = accessibility, y = TE)) + geom_point() + scale_x_log10() + scale_y_log10() + ggtitle("256cell, utr3") + ylab("TE, cds")
ggsave(file = "/Volumes/USELESS/META/SHAPES/accessibility/utr3_256.png")

# start_by_exons_24, start_by_exons_256
df <- data.frame(accessibility = as.vector(start_by_exons_24[order(names(start_by_exons_24))]), TE = TE_24[rownames(TE_24) %in% names(start_by_exons_24),]$cds)
ggplot(df, aes(x = accessibility, y = TE)) + geom_point() + scale_x_log10() + scale_y_log10() + ggtitle("2-4cell, start (-20:-1)") + ylab("TE, cds")
ggsave(file = "/Volumes/USELESS/META/SHAPES/accessibility/start_24.png")

df <- data.frame(accessibility = as.vector(start_by_exons_256[order(names(start_by_exons_256))]), TE = TE_256[rownames(TE_256) %in% names(start_by_exons_256),]$cds)
ggplot(df, aes(x = accessibility, y = TE)) + geom_point() + scale_x_log10() + scale_y_log10() + ggtitle("256cell, start (-20:-1)") + ylab("TE, cds")
ggsave(file = "/Volumes/USELESS/META/SHAPES/accessibility/start_256.png")


#######
# shape vs RNA-seq
cds_avg_24 <- sapply(c24, function(x){mean(x$dtcr)})
cds_avg_256 <- sapply(c256, function(x){mean(x$dtcr)})

df <- data.frame(avg_shape = as.vector(cds_avg_24), rna = FPKM_24[rownames(FPKM_24) %in% names(cds_avg_24),]$cds_RNA_fpkm)
ggplot(df, aes(x = log2(rna), y = avg_shape)) + geom_point() + ggtitle("2-4cell, CDS")
ggsave(file = "/Volumes/USELESS/META/SHAPES/accessibility/avg_shape_vs_RNA_CDS_24.png")

df <- data.frame(avg_shape = as.vector(cds_avg_256), rna = FPKM_256[rownames(FPKM_256) %in% names(cds_avg_256),]$cds_RNA_fpkm)
ggplot(df, aes(x = log2(rna), y = avg_shape)) + geom_point() + ggtitle("256cell, CDS")
ggsave(file = "/Volumes/USELESS/META/SHAPES/accessibility/avg_shape_vs_RNA_CDS_256.png")

# shape vs Ribo-seq
df <- data.frame(avg_shape = as.vector(cds_avg_24), ribo = FPKM_24[rownames(FPKM_24) %in% names(cds_avg_24),]$cds_Ribo_fpkm)
ggplot(df, aes(x = log2(ribo), y = avg_shape)) + geom_point() + ggtitle("2-4cell, CDS")
ggsave(file = "/Volumes/USELESS/META/SHAPES/accessibility/ribo_vs_shape_24.png")

df <- data.frame(avg_shape = as.vector(cds_avg_256), ribo = FPKM_256[rownames(FPKM_256) %in% names(cds_avg_256),]$cds_Ribo_fpkm)
ggplot(df, aes(x = log2(ribo), y = avg_shape)) + geom_point() + ggtitle("256cell, CDS")
ggsave(file = "/Volumes/USELESS/META/SHAPES/accessibility/ribo_vs_shape_256.png")

## shape 20nt before start vs RNA
start_sum_24 <- sapply(names(s24), function(x){sum(s24[[x]]$dtcr[(txLengths[x,]$utr5_len-20):(txLengths[x,]$utr5_len)])})

df <- data.frame(sum_shape_20nt_upstream_of_start = as.vector(start_sum_24), rna = FPKM_24[rownames(FPKM_24) %in% names(exons_sum_24),]$exons_RNA_fpkm)
ggplot(df, aes(x = rna, y = sum_shape_20nt_upstream_of_start)) + geom_point() + ggtitle("2-4cell") + scale_x_log10() +  scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES/accessibility/start_vs_RNA_CDS_24.png")

df <- data.frame(sum_shape_20nt_upstream_of_start = as.vector(start_sum_24), TE = TE_24[rownames(TE_24) %in% names(start_sum_24),]$cds)
ggplot(df, aes(x = TE, y = sum_shape_20nt_upstream_of_start)) + geom_point() + ggtitle("2-4cell") + scale_x_log10() +  scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES/accessibility/start_vs_TE_CDS_24.png")

## shape 20nt before start vs RNA
start_sum_256 <- sapply(names(s256), function(x){sum(s256[[x]]$dtcr[(txLengths[x,]$utr5_len-20):(txLengths[x,]$utr5_len)])})

df <- data.frame(sum_shape_20nt_upstream_of_start = as.vector(start_sum_256), rna = FPKM_256[rownames(FPKM_256) %in% names(exons_sum_256),]$exons_RNA_fpkm)
ggplot(df, aes(x = rna, y = sum_shape_20nt_upstream_of_start)) + geom_point() + ggtitle("256cell") + scale_x_log10() +  scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES/accessibility/start_vs_RNA_CDS_256.png")

df <- data.frame(sum_shape_20nt_upstream_of_start = as.vector(start_sum_256), TE_cds = TE_256[rownames(TE_256) %in% names(start_sum_256),]$cds)
ggplot(df, aes(x = TE_cds, y = sum_shape_20nt_upstream_of_start)) + geom_point() + ggtitle("256cell") + scale_x_log10() +  scale_y_log10()
ggsave(file = "/Volumes/USELESS/META/SHAPES/accessibility/start_vs_TE_CDS_256.png")

##
# shape vs RNA-seq
exons_avg_24 <- sapply(ex24, function(x){mean(x$dtcr)})
exons_avg_256 <- sapply(ex256, function(x){mean(x$dtcr)})

df <- data.frame(avg_shape = as.vector(exons_avg_24), rna = FPKM_24[rownames(FPKM_24) %in% names(exons_avg_24),]$exons_RNA_fpkm)
ggplot(df, aes(x = log2(rna), y = avg_shape)) + geom_point() + ggtitle("2-4cell, tx")
ggsave(file = "/Volumes/USELESS/META/SHAPES/accessibility/avg_shape_vs_RNA_tx_24.png")

df <- data.frame(avg_shape = as.vector(exons_avg_256), rna = FPKM_256[rownames(FPKM_256) %in% names(exons_avg_256),]$exons_RNA_fpkm)
ggplot(df, aes(x = log2(rna), y = avg_shape)) + geom_point() + ggtitle("256cell, tx")
ggsave(file = "/Volumes/USELESS/META/SHAPES/accessibility/avg_shape_vs_RNA_tx_256.png")

###############################################
# start_by_exons_24, start_by_exons_256, Ribo
df <- data.frame(accessibility = as.vector(start_by_exons_24[order(names(start_by_exons_24))]), ribo = FPKM_24[rownames(FPKM_24) %in% names(start_by_exons_24),]$cds_Ribo_fpkm)
ggplot(df, aes(x = accessibility, y = ribo)) + geom_point() + scale_x_log10() + scale_y_log10() + ggtitle("2-4cell, start (-20:-1)")
ggsave(file = "/Volumes/USELESS/META/SHAPES/accessibility/start_Ribo_24.png")

df <- data.frame(accessibility = as.vector(start_by_exons_256[order(names(start_by_exons_256))]), ribo = FPKM_256[rownames(FPKM_256) %in% names(start_by_exons_256),]$cds_Ribo_fpkm)
ggplot(df, aes(x = accessibility, y = ribo)) + geom_point() + scale_x_log10() + scale_y_log10() + ggtitle("256cell, start (-20:-1)")
ggsave(file = "/Volumes/USELESS/META/SHAPES/accessibility/start_Ribo_256.png")

# exons by exons, Ribo
df <- data.frame(accessibility = as.vector(exons_by_exons_24[order(names(exons_by_exons_24))]), ribo = FPKM_24[rownames(FPKM_24) %in% names(exons_by_exons_24),]$cds_Ribo_fpkm)
ggplot(df, aes(x = accessibility, y = ribo)) + geom_point() + scale_x_log10() + scale_y_log10() + ggtitle("2-4cell, whole tx")
ggsave(file = "/Volumes/USELESS/META/SHAPES/accessibility/tx_Ribo_24.png")

df <- data.frame(accessibility = as.vector(exons_by_exons_256[order(names(exons_by_exons_256))]), ribo = FPKM_256[rownames(FPKM_256) %in% names(exons_by_exons_256),]$cds_Ribo_fpkm)
ggplot(df, aes(x = accessibility, y = ribo)) + geom_point() + scale_x_log10() + scale_y_log10() + ggtitle("256cell, whole tx")
ggsave(file = "/Volumes/USELESS/META/SHAPES/accessibility/tx_Ribo_256.png")

# utr5 by exons, Ribo
df <- data.frame(accessibility = as.vector(utr5_by_exons_24[order(names(utr5_by_exons_24))]), ribo = FPKM_24[rownames(FPKM_24) %in% names(utr5_by_exons_24),]$cds_Ribo_fpkm)
ggplot(df, aes(x = accessibility, y = ribo)) + geom_point() + scale_x_log10() + scale_y_log10() + ggtitle("2-4cell, utr5")
ggsave(file = "/Volumes/USELESS/META/SHAPES/accessibility/utr5_Ribo_24.png")

df <- data.frame(accessibility = as.vector(utr5_by_exons_256[order(names(utr5_by_exons_256))]), ribo = FPKM_256[rownames(FPKM_256) %in% names(utr5_by_exons_256),]$cds_Ribo_fpkm)
ggplot(df, aes(x = accessibility, y = ribo)) + geom_point() + scale_x_log10() + scale_y_log10() + ggtitle("256cell, utr5")
ggsave(file = "/Volumes/USELESS/META/SHAPES/accessibility/utr5_Ribo_256.png")

# cds by exons, Ribo
df <- data.frame(accessibility = as.vector(cds_by_exons_24[order(names(cds_by_exons_24))]), ribo = FPKM_24[rownames(FPKM_24) %in% names(cds_by_exons_24),]$cds_Ribo_fpkm)
ggplot(df, aes(x = accessibility, y = ribo)) + geom_point() + scale_x_log10() + scale_y_log10() + ggtitle("2-4cell, cds")
ggsave(file = "/Volumes/USELESS/META/SHAPES/accessibility/cds_Ribo_24.png")

df <- data.frame(accessibility = as.vector(cds_by_exons_256[order(names(cds_by_exons_256))]), ribo = FPKM_256[rownames(FPKM_256) %in% names(cds_by_exons_256),]$cds_Ribo_fpkm)
ggplot(df, aes(x = accessibility, y = ribo)) + geom_point() + scale_x_log10() + scale_y_log10() + ggtitle("256cell, cds")
ggsave(file = "/Volumes/USELESS/META/SHAPES/accessibility/cds_Ribo_256.png")

# utr3 by exons, Ribo
df <- data.frame(accessibility = as.vector(utr3_by_exons_24[order(names(utr3_by_exons_24))]), ribo = FPKM_24[rownames(FPKM_24) %in% names(utr3_by_exons_24),]$cds_Ribo_fpkm)
ggplot(df, aes(x = accessibility, y = ribo)) + geom_point() + scale_x_log10() + scale_y_log10() + ggtitle("2-4cell, utr3")
ggsave(file = "/Volumes/USELESS/META/SHAPES/accessibility/utr3_Ribo_24.png")

df <- data.frame(accessibility = as.vector(utr3_by_exons_256[order(names(utr3_by_exons_256))]), ribo = FPKM_256[rownames(FPKM_256) %in% names(utr3_by_exons_256),]$cds_Ribo_fpkm)
ggplot(df, aes(x = accessibility, y = ribo)) + geom_point() + scale_x_log10() + scale_y_log10() + ggtitle("256cell, utr3")
ggsave(file = "/Volumes/USELESS/META/SHAPES/accessibility/utr3_Ribo_256.png")
