# ...running on kjempetuja 11670.pts-14.kjempetuja

### make big shapeSeq data frame
library(GenomicRanges)
library(GenomicFeatures)


### get transcript dtabase
txdb_can <- makeTxDbFromGFF("/Home/ii/katchyz/DATA/genomes/GTF/zebrafish_GRCz10.gtf", format = "gtf")
cds <- cdsBy(txdb_can, by="tx", use.names=TRUE)

#unsmoothed
load("/Home/ii/katchyz/META/SHAPES/gc_unsmoothed_2_4.Rsave")
load("/Home/ii/katchyz/META/SHAPES/gc_unsmoothed_256.Rsave")
load("/Volumes/USELESS/META/SHAPES/gc_unsmoothed_256.Rsave")

# 2-4cell
cds_2_4 <- subsetByOverlaps(gc_unsmoothed_2_4, cds)
cds_2_4 <- split(cds_2_4, cds_2_4$trnames)


# 256cell
cds_256 <- subsetByOverlaps(gc_unsmoothed_256, cds)
cds_256 <- split(cds_256, cds_256$trnames)


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

###########
# 2-4cell
orf <- c("1", "2", "3")
frequency_paired_nts <- c(0.05002587, 0.04904469, 0.04951098)
df <- data.frame(orf, frequency_paired_nts)

ggplot(df, aes(x=orf, y=frequency_paired_nts)) + geom_bar(stat = "identity")


## 2-4Cell
# orf1
cds_2_4_orf1 <- sum(sapply(cds_2_4, function(x){sum(x[seq(1, length(x), 3)]$dtcr)}))
# orf2
cds_2_4_orf2 <- sum(sapply(cds_2_4, function(x){sum(x[seq(2, length(x), 3)]$dtcr)}))
# orf3
cds_2_4_orf3 <- sum(sapply(cds_2_4, function(x){sum(x[seq(3, length(x), 3)]$dtcr)}))


## 2-4Cell
# orf1
cds_256_orf1 <- sum(sapply(cds_256, function(x){sum(x[seq(1, length(x), 3)]$dtcr)}))
# orf2
cds_256_orf2 <- sum(sapply(cds_256, function(x){sum(x[seq(2, length(x), 3)]$dtcr)}))
# orf3
cds_256_orf3 <- sum(sapply(cds_2_4, function(x){sum(x[seq(3, length(x), 3)]$dtcr)}))



######
######

# 2-4cell
orf <- c("1", "2", "3")
sum_reactivities <- c(86801.97, 82753.9, 89249.82)
df <- data.frame(orf, sum_reactivities)
p_2_4 <- ggplot(df, aes(x=orf, y=sum_reactivities)) + geom_bar(stat = "identity")
ggsave("/Volumes/USELESS/META/SHAPES/periodicity/orf_2_4.png", plot = p_2_4)

# 256cell
orf <- c("1", "2", "3")
sum_reactivities <- c(111118.7, 102985.3, 89249.82)
df <- data.frame(orf, sum_reactivities)
p_256 <- ggplot(df, aes(x=orf, y=sum_reactivities)) + geom_bar(stat = "identity")
ggsave("/Volumes/USELESS/META/SHAPES/periodicity/orf_256.png", plot = p_256)


