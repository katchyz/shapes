data = read.csv("/Users/kasia/Documents/PhD/TE.csv", header = TRUE)
rownames(data) <- data$X.transcript_id

# 2-4 cell, slograt
leaders <- c()
cdss <- c()
#trailers <- c()
te_leader <- c()
te_cds <- c()

for (tr in seqnames(seqinfo(slograt_2_4))) {
  slr <- elementMetadata(slograt_2_4)$slograt[as.vector(seqnames(slograt_2_4)) == tr]
  tr_name <- substr(tr,1,18)
  if (length(fiveUTR[[tr_name]]) > 0) {
    if (length(threeUTR[[tr_name]]) > 0) { # if fiveUTR AND threeUTR exist
      if (data[tr_name, "X01_2to4cell_lead_te"] < 10) {
        if (data[tr_name, "X01_2to4cell_gene_te"] < 10) {
          len_f <- sum(width(ranges(fiveUTR[[tr_name]])))
          len_c <- sum(width(ranges(cds[[tr_name]])))
          #len_t <- sum(width(ranges(threeUTR[[tr_name]])))
          leader <- slr[1:len_f]
          coding <- slr[(len_f+1):(len_f+len_c)]
          #trailer <- slr[(len_f+len_c+1):length(slr)]
          leaders <- c(leaders, sum(leader)/len_f)
          cdss <- c(cdss, sum(coding)/len_c)
          #trailers <- c(trailers, sum(trailer)/len_t)
          te_leader <- c(te_leader, data[tr_name, "X01_2to4cell_lead_te"])
          te_cds <- c(te_cds, data[tr_name, "X01_2to4cell_gene_te"])
        }
      }
    }
  }
}

df_l <- data.frame(leaders,te_leader)
df_c <- data.frame(cdss,te_cds)

p1a <- ggplot(df_l, aes(x=leaders, y=te_leader)) + geom_point() + xlab("sum fo reactivities/length") + ylab("TE") + ggtitle("leaders, 2-4cell, slograt")

p1b <- ggplot(df_l, aes(x=cdss, y=te_cds)) + geom_point() + xlab("sum fo reactivities/length") + ylab("TE") + ggtitle("CDSs, 2-4cell, slograt")

# 2-4 cell, dtcr
leaders <- c()
cdss <- c()
#trailers <- c()
te_leader <- c()
te_cds <- c()

for (tr in seqnames(seqinfo(dtcr_2_4))) {
  slr <- elementMetadata(dtcr_2_4)$dtcr[as.vector(seqnames(dtcr_2_4)) == tr]
  tr_name <- substr(tr,1,18)
  if (length(fiveUTR[[tr_name]]) > 0) {
    if (length(threeUTR[[tr_name]]) > 0) { # if fiveUTR AND threeUTR exist
      if (data[tr_name, "X01_2to4cell_lead_te"] < 10) {
        if (data[tr_name, "X01_2to4cell_gene_te"] < 10) {
          len_f <- sum(width(ranges(fiveUTR[[tr_name]])))
          len_c <- sum(width(ranges(cds[[tr_name]])))
          #len_t <- sum(width(ranges(threeUTR[[tr_name]])))
          leader <- slr[1:len_f]
          coding <- slr[(len_f+1):(len_f+len_c)]
          #trailer <- slr[(len_f+len_c+1):length(slr)]
          leaders <- c(leaders, sum(leader)/len_f)
          cdss <- c(cdss, sum(coding)/len_c)
          #trailers <- c(trailers, sum(trailer)/len_t)
          te_leader <- c(te_leader, data[tr_name, "X01_2to4cell_lead_te"])
          te_cds <- c(te_cds, data[tr_name, "X01_2to4cell_gene_te"])
        }
      }
    }
  }
}

df_l <- data.frame(leaders,te_leader)
df_c <- data.frame(cdss,te_cds)

p2a <- ggplot(df_l, aes(x=leaders, y=te_leader)) + geom_point() + xlab("sum fo reactivities/length") + ylab("TE") + ggtitle("leaders, 2-4cell, dtcr") + scale_x_log10() + scale_y_log10()

p2b <- ggplot(df_l, aes(x=cdss, y=te_cds)) + geom_point() + xlab("sum fo reactivities/length") + ylab("TE") + ggtitle("CDSs, 2-4cell, dtcr") + scale_x_log10() + scale_y_log10()


#####

# 256 cell, slograt
leaders <- c()
cdss <- c()
#trailers <- c()
te_leader <- c()
te_cds <- c()

for (tr in seqnames(seqinfo(slograt_256))) {
  slr <- elementMetadata(slograt_256)$slograt[as.vector(seqnames(slograt_256)) == tr]
  tr_name <- substr(tr,1,18)
  if (length(fiveUTR[[tr_name]]) > 0) {
    if (length(threeUTR[[tr_name]]) > 0) { # if fiveUTR AND threeUTR exist
      if (data[tr_name, "X02_256cell_lead_te"] < 10) {
        if (data[tr_name, "X02_256cell_gene_te"] < 10) {
          len_f <- sum(width(ranges(fiveUTR[[tr_name]])))
          len_c <- sum(width(ranges(cds[[tr_name]])))
          #len_t <- sum(width(ranges(threeUTR[[tr_name]])))
          leader <- slr[1:len_f]
          coding <- slr[(len_f+1):(len_f+len_c)]
          #trailer <- slr[(len_f+len_c+1):length(slr)]
          leaders <- c(leaders, sum(leader)/len_f)
          cdss <- c(cdss, sum(coding)/len_c)
          #trailers <- c(trailers, sum(trailer)/len_t)
          te_leader <- c(te_leader, data[tr_name, "X02_256cell_lead_te"])
          te_cds <- c(te_cds, data[tr_name, "X02_256cell_gene_te"])
        }
      }
    }
  }
}

df_l <- data.frame(leaders,te_leader)
df_c <- data.frame(cdss,te_cds)

p3a <- ggplot(df_l, aes(x=leaders, y=te_leader)) + geom_point() + xlab("sum fo reactivities/length") + ylab("TE") + ggtitle("leaders, 256cell, slograt")

p3b <- ggplot(df_l, aes(x=cdss, y=te_cds)) + geom_point() + xlab("sum fo reactivities/length") + ylab("TE") + ggtitle("CDSs, 256cell, slograt")


# 256 cell, dtcr
leaders <- c()
cdss <- c()
#trailers <- c()
te_leader <- c()
te_cds <- c()

for (tr in seqnames(seqinfo(dtcr_256))) {
  slr <- elementMetadata(dtcr_256)$dtcr[as.vector(seqnames(dtcr_256)) == tr]
  tr_name <- substr(tr,1,18)
  if (length(fiveUTR[[tr_name]]) > 0) {
    if (length(threeUTR[[tr_name]]) > 0) { # if fiveUTR AND threeUTR exist
      if (data[tr_name, "X02_256cell_lead_te"] < 10) {
        if (data[tr_name, "X02_256cell_gene_te"] < 10) {
          len_f <- sum(width(ranges(fiveUTR[[tr_name]])))
          len_c <- sum(width(ranges(cds[[tr_name]])))
          #len_t <- sum(width(ranges(threeUTR[[tr_name]])))
          leader <- slr[1:len_f]
          coding <- slr[(len_f+1):(len_f+len_c)]
          #trailer <- slr[(len_f+len_c+1):length(slr)]
          leaders <- c(leaders, sum(leader)/len_f)
          cdss <- c(cdss, sum(coding)/len_c)
          #trailers <- c(trailers, sum(trailer)/len_t)
          te_leader <- c(te_leader, data[tr_name, "X02_256cell_lead_te"])
          te_cds <- c(te_cds, data[tr_name, "X02_256cell_gene_te"])
        }
      }
    }
  }
}

df_l <- data.frame(leaders,te_leader)
df_c <- data.frame(cdss,te_cds)

p4a <- ggplot(df_l, aes(x=leaders, y=te_leader)) + geom_point() + xlab("sum fo reactivities/length") + ylab("TE") + ggtitle("leaders, 256cell, dtcr") + scale_x_log10() + scale_y_log10()

p4b <- ggplot(df_l, aes(x=cdss, y=te_cds)) + geom_point() + xlab("sum fo reactivities/length") + ylab("TE") + ggtitle("CDSs, 256cell, dtcr") + scale_x_log10() + scale_y_log10()



ggsave(file = "/Volumes/USELESS/META/SHAPES/te/p1a.png", plot = p1a)
ggsave(file = "/Volumes/USELESS/META/SHAPES/te/p2a.png", plot = p2a)
ggsave(file = "/Volumes/USELESS/META/SHAPES/te/p3a.png", plot = p3a)
ggsave(file = "/Volumes/USELESS/META/SHAPES/te/p4a.png", plot = p4a)
ggsave(file = "/Volumes/USELESS/META/SHAPES/te/p1b.png", plot = p1b)
ggsave(file = "/Volumes/USELESS/META/SHAPES/te/p2b.png", plot = p2b)
ggsave(file = "/Volumes/USELESS/META/SHAPES/te/p3b.png", plot = p3b)
ggsave(file = "/Volumes/USELESS/META/SHAPES/te/p4b.png", plot = p4b)

#convert p1a.png p1b.png p2a.png p2b.png +append 2-4cell.png
#convert p3a.png p3b.png p4a.png p4b.png +append 256cell.png
#convert 2-4cell.png 256cell.png -append shape_vs_te.png


### log
ggsave(file = "/Volumes/USELESS/META/SHAPES/te/log/p2a.png", plot = p2a)
ggsave(file = "/Volumes/USELESS/META/SHAPES/te/log/p4a.png", plot = p4a)
ggsave(file = "/Volumes/USELESS/META/SHAPES/te/log/p2b.png", plot = p2b)
ggsave(file = "/Volumes/USELESS/META/SHAPES/te/log/p4b.png", plot = p4b)

#convert p2a.png p2b.png +append 2-4cell.png
#convert p4a.png p4b.png +append 256cell.png
#convert 2-4cell.png 256cell.png -append shape_vs_te_log.png