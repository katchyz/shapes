exon_2_4_list

dtcr5_2_4 <- sapply(exon_2_4_list, function(x){as.numeric(x[1:5]$dtcr)})
dtcr5_256 <- sapply(exon_256_list, function(x){as.numeric(x[1:5]$dtcr)})

# get common names, put in one matrix (or data frame)

dtcr5_2_4 <- t(dtcr5_2_4)
dtcr5_256 <- t(dtcr5_256)

dtcr5_2_4 <- dtcr5_2_4[rownames(dtcr5_2_4) %in% rownames(dtcr5_256),]
dtcr5_256 <- dtcr5_256[rownames(dtcr5_256) %in% rownames(dtcr5_2_4),]

dtcr5 <- cbind(dtcr5_2_4, dtcr5_256)
colnames(dtcr5) <- c('2-4cell_1', '2-4cell_2', '2-4cell_3', '2-4cell_4', '2-4cell_5', '256cell_1', '256cell_2', '256cell_3', '256cell_4', '256cell_5')

save(dtcr5, file = "/Volumes/USELESS/META/SHAPES/dtcr5.Rsave")
write.csv(dtcr5, file = "/Volumes/USELESS/META/SHAPES/dtcr5.csv")
