### CAI vs ORF score (avg?)

cai <- read.csv("/Volumes/USELESS/META/SHAPES/caiDf.csv")

### get unsmoothed CDSs !!!!!!!!!!!!!!!!!!!

names(cds_2_4) <- sapply(names(cds_2_4), function(x){substr(x,1,18)})
names(cds_256) <- sapply(names(cds_256), function(x){substr(x,1,18)})

sam <- cds_2_4[1:10]
orf1_sam <- sapply(sam, function(x){sum(x[seq(1, length(x), 3)]$dtcr)/floor(length(x)/3)})

## 2-4Cell
# orf1
orf1_2_4 <- sapply(cds_2_4, function(x){sum(x[seq(1, length(x), 3)]$dtcr)/floor(length(x)/3)})
# orf2
orf2_2_4 <- sapply(cds_2_4, function(x){sum(x[seq(2, length(x), 3)]$dtcr)/floor(length(x)/3)})
# orf3
orf3_2_4 <- sapply(cds_2_4, function(x){sum(x[seq(3, length(x), 3)]$dtcr)/floor(length(x)/3)})


## 256Cell
# orf1
orf1_256 <- sapply(cds_256, function(x){sum(x[seq(1, length(x), 3)]$dtcr)/floor(length(x)/3)})
# orf2
orf2_256 <- sapply(cds_256, function(x){sum(x[seq(2, length(x), 3)]$dtcr)/floor(length(x)/3)})
# orf3
orf3_256 <- sapply(cds_256, function(x){sum(x[seq(3, length(x), 3)]$dtcr)/floor(length(x)/3)})

### put orf values in data frame with CAI, plot
cai_2_4 <- cai[rownames(cai) %in% names(cds_2_4),]
cai_256 <- cai[rownames(cai) %in% names(cds_256),]

orf1_2_4 <- orf1_2_4[names(orf1_2_4) %in% rownames(cai_2_4)]
orf2_2_4 <- orf2_2_4[names(orf2_2_4) %in% rownames(cai_2_4)]
orf3_2_4 <- orf3_2_4[names(orf3_2_4) %in% rownames(cai_2_4)]

CAI_ORF_2_4 <- data.frame(cai_2_4$CAI, orf1_2_4, orf2_2_4, orf3_2_4) 

orf1_256 <- orf1_256[names(orf1_256) %in% rownames(cai_256)]
orf2_256 <- orf2_256[names(orf2_256) %in% rownames(cai_256)]
orf3_256 <- orf3_256[names(orf3_256) %in% rownames(cai_256)]

CAI_ORF_256 <- data.frame(cai_256$CAI, orf1_256, orf2_256, orf3_256)

save(CAI_ORF_2_4, file = "/Volumes/USELESS/META/SHAPES/CAI_ORF_2_4.Rsave")
save(CAI_ORF_256, file = "/Volumes/USELESS/META/SHAPES/CAI_ORF_256.Rsave")

ggplot(CAI_ORF_256, aes(x = orf1_256, y = cai_256.CAI)) + geom_point() + scale_x_log10() + scale_y_log10() + xlab("orf1") + ylab("CAI")
ggsave("/Volumes/USELESS/META/SHAPES/CAI/orf1_256.png")

ggplot(CAI_ORF_2_4, aes(x = orf3_2_4, y = cai_2_4.CAI)) + geom_point() + scale_x_log10() + scale_y_log10() + xlab("orf3") + ylab("CAI")
ggsave("/Volumes/USELESS/META/SHAPES/CAI/orf3_2_4.png")

