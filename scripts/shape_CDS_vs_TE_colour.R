library(GenomicFeatures)
library(ggplot2)

data = read.csv("/Users/kasia/Documents/PhD/TE.csv", header = TRUE)
rownames(data) <- data$X.transcript_id
data_orig <- data # save just in case

load("/Volumes/USELESS/OUT/SHAPES_out/rsave/gc_dtcr_2_4.Rsave")
load("/Volumes/USELESS/OUT/SHAPES_out/rsave/gc_dtcr_256.Rsave")
txdb_can <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/danio_rerio/zebrafish_GRCz10.gtf", format = "gtf")
cds <- cdsBy(txdb_can, by="tx", use.names=TRUE)

cds_2_4 <- subsetByOverlaps(gc_dtcr_2_4, cds)
cds_2_4_list <- split(cds_2_4, cds_2_4$trnames)

cds_256 <- subsetByOverlaps(gc_dtcr_256, cds)
cds_256_list <- split(cds_256, cds_256$trnames)

# X01_2to4cell_gene_te
# X02_256cell_gene_te

data <- data_orig
sum_cds <- sapply(cds_2_4_list, function(x){sum(x$dtcr)/length(x$dtcr)})
names(sum_cds) <- sapply(names(sum_cds), function(x){substr(x,1,18)})
sum_cds <- sum_cds[names(sum_cds) %in% rownames(data)]
data_2_4 <- data[rownames(data) %in% names(sum_cds),]
shape_te_cds <- data.frame(as.numeric(sum_cds), data_2_4$X01_2to4cell_gene_te)
colnames(shape_te_cds) <- c("shape_cds", "te_cds")
shape_te_cds_2_4 <- shape_te_cds

data <- data_orig
sum_cds <- sapply(cds_256_list, function(x){sum(x$dtcr)/length(x$dtcr)})
names(sum_cds) <- sapply(names(sum_cds), function(x){substr(x,1,18)})
sum_cds <- sum_cds[names(sum_cds) %in% rownames(data)]
data_256 <- data[rownames(data) %in% names(sum_cds),]
shape_te_cds <- data.frame(as.numeric(sum_cds), data_256$X02_256cell_gene_te)
colnames(shape_te_cds) <- c("shape_cds", "te_cds")
shape_te_cds_256 <- shape_te_cds

######

te_log_2_4 <- data.frame(log10(shape_te_cds_2_4$shape_cds), log10(shape_te_cds_2_4$te_cds))
te_log_256 <- data.frame(log10(shape_te_cds_256$shape_cds), log10(shape_te_cds_256$te_cds))
rownames(te_log_2_4) <- rownames(data_2_4)
rownames(te_log_256) <- rownames(data_256)

colnames(te_log_2_4) <- c("shape_cds", "te_cds")
colnames(te_log_256) <- c("shape_cds", "te_cds")

high_2_4 <- subset(te_log_2_4, shape_cds > -4)
high_2_4 <- subset(high_2_4, te_cds > 2)
high_256 <- subset(te_log_256, shape_cds > -4)
high_256 <- subset(high_256, te_cds > 2)

#data <- data_orig
#data <- data[rownames(data) %in% names(sum_five_2_4),]
#high_GO_2_4_trnames <- names(sum_five_2_4[as.numeric(rownames(high_2_4))])
#high_GO_2_4 <- data[as.numeric(rownames(high_2_4)),]$GO_Term_Name

#data <- data_orig
#data <- data[rownames(data) %in% names(sum_five_256),]
#high_GO_256_trnames <- names(sum_five_256[as.numeric(rownames(high_256))])
#high_GO_256 <- data[as.numeric(rownames(high_256)),]$GO_Term_Name


### high reactivities (less structure) - highly expressed in 2-4; does the structure increase (drop in reactivities) in 256?
# sum_five_2_4[high_GO_2_4_trnames]
# sum_five_256[high_GO_2_4_trnames]

### low reactivities (more structure) - lower expression in 2-4; does the structure decrease (increase in reactivities) in 256?
# sum_five_2_4[high_GO_256_trnames]
# sum_five_256[high_GO_256_trnames]

te_log_2_4 = within(te_log_2_4, {
  high_te = ifelse(rownames(te_log_2_4) %in% rownames(high_256), "H", "N")
})

p_2_4 <- ggplot(te_log_2_4, aes(x=shape_cds, y=te_cds, color=high_te)) + geom_point() + xlab("sum of reactivities, CDS") + ylab("TE of CDS") + ggtitle("2_4cell") + geom_point(data = subset(te_log_2_4, high_te == 'H'), aes(x=shape_cds, y=te_cds, color=high_te)) + theme(legend.position="none")

ggsave(file = "/Volumes/USELESS/META/SHAPES/shape_cds_2_4.png", plot = p_2_4)


te_log_256 = within(te_log_256, {
  high_te = ifelse(rownames(te_log_256) %in% rownames(high_2_4), "H", "N")
})

p_256 <- ggplot(te_log_256, aes(x=shape_cds, y=te_cds, color=high_te)) + geom_point() + xlab("sum of reactivities, CDS") + ylab("TE of CDS") + ggtitle("256cell") + geom_point(data = subset(te_log_256, high_te == 'H'), aes(x=shape_cds, y=te_cds, color=high_te)) + theme(legend.position="none")

ggsave(file = "/Volumes/USELESS/META/SHAPES/shape_cds_256.png", plot = p_256)



write(rownames(high_2_4), file = "/Volumes/USELESS/META/SHAPES/tower_2_4cell.txt", ncolumns = 1)
###############
data_2_4[rownames(high_2_4),]$gene_name
data_2_4[rownames(high_2_4),]$GO_Term_Name

df_2_4 <- data.frame(data_2_4[rownames(high_2_4),]$gene_name, data_2_4[rownames(high_2_4),]$X01_2to4cell_gene_te, data_2_4[rownames(high_2_4),]$GO_Term_Name)
rownames(df_2_4) <- rownames(high_2_4)
colnames(df_2_4) <- c("gene_name", "TE", "GO_Term_Name")
write.table(df_2_4, file = "/Volumes/USELESS/META/SHAPES/tower_2_4cell_GO.txt", sep = "\t")

df_256 <- data.frame(data_256[rownames(high_256),]$gene_name, data_256[rownames(high_256),]$X02_256cell_gene_te, data_256[rownames(high_256),]$GO_Term_Name)
rownames(df_256) <- rownames(high_256)
colnames(df_256) <- c("gene_name", "TE", "GO_Term_Name")
write.table(df_256, file = "/Volumes/USELESS/META/SHAPES/tower_256cell_GO.txt", sep = "\t")
