library(GenomicFeatures)
library(ggplot2)

data = read.csv("/Users/kasia/Documents/PhD/TE.csv", header = TRUE)
rownames(data) <- data$X.transcript_id
data_orig <- data # save just in case

load("/Volumes/USELESS/OUT/SHAPES_out/rsave/gc_dtcr_2_4.Rsave")
load("/Volumes/USELESS/OUT/SHAPES_out/rsave/gc_dtcr_256.Rsave")
txdb_can <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/danio_rerio/zebrafish_GRCz10.gtf", format = "gtf")
fiveUTR <- fiveUTRsByTranscript(txdb_can, use.names=TRUE)

five_2_4 <- subsetByOverlaps(gc_dtcr_2_4, fiveUTR)
five_2_4_list <- split(five_2_4, five_2_4$trnames)

five_256 <- subsetByOverlaps(gc_dtcr_256, fiveUTR)
five_256_list <- split(five_256, five_256$trnames)

# X01_2to4cell_gene_te
# X02_256cell_gene_te

sum_five <- sapply(five_2_4_list, function(x){sum(x$dtcr)})
names(sum_five) <- sapply(names(sum_five), function(x){substr(x,1,18)})
sum_five <- sum_five[names(sum_five) %in% rownames(data)]
data <- data[rownames(data) %in% names(sum_five),]
shape_five_tecds <- data.frame(as.numeric(sum_five), data$X01_2to4cell_gene_te)
colnames(shape_five_tecds) <- c("leader", "te_cds")
shape_five_tecds_2_4 <- shape_five_tecds

sum_five <- sapply(five_256_list, function(x){sum(x$dtcr)})
names(sum_five) <- sapply(names(sum_five), function(x){substr(x,1,18)})
sum_five <- sum_five[names(sum_five) %in% rownames(data)]
data <- data[rownames(data) %in% names(sum_five),]
shape_five_tecds <- data.frame(as.numeric(sum_five), data$X02_256cell_gene_te)
colnames(shape_five_tecds) <- c("leader", "te_cds")
shape_five_tecds_256 <- shape_five_tecds

######
######## remember to update data!

te_log_2_4 <- data.frame(log10(shape_five_tecds_2_4$leader), log10(shape_five_tecds_2_4$te_cds))
te_log_256 <- data.frame(log10(shape_five_tecds_256$leader), log10(shape_five_tecds_256$te_cds))

colnames(te_log_2_4) <- c("leader", "te_cds")
colnames(te_log_256) <- c("leader", "te_cds")

high_2_4 <- subset(te_log_2_4, leader > 0)
high_2_4 <- subset(high_2_4, te_cds > 2)
high_256 <- subset(te_log_256, leader > 0)
high_256 <- subset(high_256, te_cds > 2)

data <- data_orig
data <- data[rownames(data) %in% names(sum_five_2_4),]
high_GO_2_4_trnames <- names(sum_five_2_4[as.numeric(rownames(high_2_4))])
#high_GO_2_4 <- data[as.numeric(rownames(high_2_4)),]$GO_Term_Name

data <- data_orig
data <- data[rownames(data) %in% names(sum_five_256),]
high_GO_256_trnames <- names(sum_five_256[as.numeric(rownames(high_256))])
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

p_2_4 <- ggplot(te_log_2_4, aes(x=leader, y=te_cds, color=high_te)) + geom_point() + xlab("sum of reactivities, leaders") + ylab("TE of CDS") + ggtitle("hindrance, 2_4cell") + geom_point(data = subset(te_log_2_4, high_te == 'H'), aes(x=leader, y=te_cds, color=high_te)) + theme(legend.position="none")

ggsave(file = "/Volumes/USELESS/META/SHAPES/hindrance_2_4.png", plot = p_2_4)


te_log_256 = within(te_log_256, {
  high_te = ifelse(rownames(te_log_256) %in% rownames(high_2_4), "H", "N")
})

p_256 <- ggplot(te_log_256, aes(x=leader, y=te_cds, color=high_te)) + geom_point() + xlab("sum of reactivities, leaders") + ylab("TE of CDS") + ggtitle("hindrance, 256cell") + geom_point(data = subset(te_log_256, high_te == 'H'), aes(x=leader, y=te_cds, color=high_te)) + theme(legend.position="none")

ggsave(file = "/Volumes/USELESS/META/SHAPES/hindrance_256.png", plot = p_256)
