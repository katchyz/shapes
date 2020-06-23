####################################
### make big shapeSeq data frame ###
####################################
library(GenomicRanges)
library(GenomicFeatures)

### data
data = read.csv("/Users/kasia/Documents/PhD/TE.csv", header = TRUE)
rownames(data) <- data$X.transcript_id

### get transcript dtabase
txdb_can <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/danio_rerio/zebrafish_GRCz10.gtf", format = "gtf")
cds <- cdsBy(txdb_can, by="tx", use.names=TRUE)
fiveUTR <- fiveUTRsByTranscript(txdb_can, use.names=TRUE)
threeUTR <- threeUTRsByTranscript(txdb_can, use.names=TRUE)

txLengths <- transcriptLengths(txdb_can, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
rownames(txLengths) <- txLengths$tx_name
rownames(txLengths) <- sapply(rownames(txLengths), function(x){substr(x,1,18)})

#unsmoothed
load("/Volumes/USELESS/META/SHAPES/gc_unsmoothed_2_4.Rsave")
load("/Volumes/USELESS/META/SHAPES/gc_unsmoothed_256.Rsave")

# 2-4cell
cds_2_4 <- subsetByOverlaps(gc_unsmoothed_2_4, cds)
cds_2_4 <- split(cds_2_4, cds_2_4$trnames)
five_2_4 <- subsetByOverlaps(gc_unsmoothed_2_4, fiveUTR)
five_2_4 <- split(five_2_4, five_2_4$trnames)
three_2_4 <- subsetByOverlaps(gc_unsmoothed_2_4, threeUTR)
three_2_4 <- split(three_2_4, three_2_4$trnames)

# 256cell
cds_256 <- subsetByOverlaps(gc_unsmoothed_256, cds)
cds_256 <- split(cds_256, cds_256$trnames)
five_256 <- subsetByOverlaps(gc_unsmoothed_256, fiveUTR)
five_256 <- split(five_256, five_256$trnames)
three_256 <- subsetByOverlaps(gc_unsmoothed_256, threeUTR)
three_256 <- split(three_256, three_256$trnames)

features <- c(cds_2_4, five_2_4, three_2_4, cds_256, five_256, three_256)

names(cds_2_4) <- sapply(names(cds_2_4), function(x){substr(x,1,18)})
names(five_2_4) <- sapply(names(five_2_4), function(x){substr(x,1,18)})
names(three_2_4) <- sapply(names(three_2_4), function(x){substr(x,1,18)})
names(cds_256) <- sapply(names(cds_256), function(x){substr(x,1,18)})
names(five_256) <- sapply(names(five_256), function(x){substr(x,1,18)})
names(three_256) <- sapply(names(three_256), function(x){substr(x,1,18)})

shapeSeq <- data.frame(row.names = union(names(cds_2_4), names(cds_256)))
#shapeSeq <- cbind(shapeSeq, data[, "gene_name"][match(rownames(shapeSeq), rownames(data))])

shapeSeq <- cbind(shapeSeq,data[, "gene_name"][match(rownames(shapeSeq), rownames(data))],stringsAsFactors = FALSE)
shapeSeq <- cbind(shapeSeq,data[, "GO_Term_Name"][match(rownames(shapeSeq), rownames(data))],stringsAsFactors = FALSE)

shapeSeq <- cbind(shapeSeq,data[, "X01_2to4cell_gene_te"][match(rownames(shapeSeq), rownames(data))],stringsAsFactors = FALSE)
shapeSeq <- cbind(shapeSeq,data[, "X02_256cell_gene_te"][match(rownames(shapeSeq), rownames(data))],stringsAsFactors = FALSE)

shapeSeq <- cbind(shapeSeq,data[, "X01_2to4cell_lead_te"][match(rownames(shapeSeq), rownames(data))],stringsAsFactors = FALSE)
shapeSeq <- cbind(shapeSeq,data[, "X02_256cell_lead_te"][match(rownames(shapeSeq), rownames(data))],stringsAsFactors = FALSE)

#colnames(shapeSeq) <- c("geneName", "GO", "TE_cds_2_4", "TE_cds_256", "TE_five_2_4", "TE_five_256")

# lengths
shapeSeq <- cbind(shapeSeq, txLengths[, "cds_len"][match(rownames(shapeSeq), rownames(txLengths))],stringsAsFactors = FALSE)
shapeSeq <- cbind(shapeSeq, txLengths[, "utr5_len"][match(rownames(shapeSeq), rownames(txLengths))],stringsAsFactors = FALSE)
shapeSeq <- cbind(shapeSeq, txLengths[, "utr3_len"][match(rownames(shapeSeq), rownames(txLengths))],stringsAsFactors = FALSE)

colnames(shapeSeq) <- c("geneName", "GO", "TE_cds_2_4", "TE_cds_256", "TE_five_2_4", "TE_five_256", "cds_len", "utr5_len", "utr3_len")

# avg reactivities: cds, five, three

######## frequencies of paired nucleotides
# get ones with sum(dtcr) > 0
# iterate over 1st, 2nd and 3rd nts, if it's > 0, then add +1; divide over sum

# remove NAs
cds_2_4 <- unlist(cds_2_4)
cds_2_4[is.na(cds_2_4$dtcr)]$dtcr <- 0
cds_2_4 <- split(cds_2_4, cds_2_4$trnames)

cds_256 <- unlist(cds_256)
cds_256[is.na(cds_256$dtcr)]$dtcr <- 0
cds_256 <- split(cds_256, cds_256$trnames)

# sum(dtcr) > 0
cds_2_4 <- cds_2_4[sapply(cds_2_4, function(x){sum(x$dtcr) > 0})]
cds_256 <- cds_256[sapply(cds_256, function(x){sum(x$dtcr) > 0})]

cds_2_4 <- cds_2_4[sapply(cds_2_4, function(x){length(x) > 0})]
cds_256 <- cds_256[sapply(cds_256, function(x){length(x) > 0})]

## 2-4Cell
# orf1
cds_2_4_orf1 <- sapply(cds_2_4, function(x){x[seq(1, length(x), 3)]$dtcr > 0})
cds_2_4_orf1_count <- sum(sapply(cds_2_4_orf1, function(x){length(which(x))}))
cds_2_4_orf1_total <- sum(sapply(cds_2_4_orf1, function(x){length(x)}))
cds_2_4_orf1_freq <- cds_2_4_orf1_count / cds_2_4_orf1_total * 100

# orf2
cds_2_4_orf2 <- sapply(cds_2_4, function(x){x[seq(2, length(x), 3)]$dtcr > 0})
cds_2_4_orf2_count <- sum(sapply(cds_2_4_orf2, function(x){length(which(x))}))
cds_2_4_orf2_total <- sum(sapply(cds_2_4_orf2, function(x){length(x)}))
cds_2_4_orf2_freq <- cds_2_4_orf2_count / cds_2_4_orf2_total * 100

# orf3
cds_2_4_orf3 <- sapply(cds_2_4, function(x){x[seq(3, length(x), 3)]$dtcr > 0})
cds_2_4_orf3_count <- sum(sapply(cds_2_4_orf3, function(x){length(which(x))}))
cds_2_4_orf3_total <- sum(sapply(cds_2_4_orf3, function(x){length(x)}))
cds_2_4_orf3_freq <- cds_2_4_orf3_count / cds_2_4_orf3_total * 100


## 2-4Cell
# orf1
cds_256_orf1 <- sapply(cds_256, function(x){x[seq(1, length(x), 3)]$dtcr > 0})
cds_256_orf1_count <- sum(sapply(cds_256_orf1, function(x){length(which(x))}))
cds_256_orf1_total <- sum(sapply(cds_256_orf1, function(x){length(x)}))
cds_256_orf1_freq <- cds_256_orf1_count / cds_256_orf1_total * 100

# orf2
cds_256_orf2 <- sapply(cds_256, function(x){x[seq(2, length(x), 3)]$dtcr > 0})
cds_256_orf2_count <- sum(sapply(cds_256_orf2, function(x){length(which(x))}))
cds_256_orf2_total <- sum(sapply(cds_256_orf2, function(x){length(x)}))
cds_256_orf2_freq <- cds_256_orf2_count / cds_256_orf2_total * 100

# orf3
cds_256_orf3 <- sapply(cds_2_4, function(x){x[seq(3, length(x), 3)]$dtcr > 0})
cds_256_orf3_count <- sum(sapply(cds_256_orf3, function(x){length(which(x))}))
cds_256_orf3_total <- sum(sapply(cds_256_orf3, function(x){length(x)}))
cds_256_orf3_freq <- cds_256_orf3_count / cds_256_orf3_total * 100


