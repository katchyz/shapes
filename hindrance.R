data = read.csv("/Users/kasia/Documents/PhD/TE.csv", header = TRUE)
rownames(data) <- data$X.transcript_id

####### not normalized by length
### hindrance
library(reshape2)

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
          leaders <- c(leaders, sum(leader))
          cdss <- c(cdss, sum(coding))
          #trailers <- c(trailers, sum(trailer)/len_t)
          te_leader <- c(te_leader, data[tr_name, "X01_2to4cell_lead_te"])
          te_cds <- c(te_cds, data[tr_name, "X01_2to4cell_gene_te"])
        }
      }
    }
  }
}


df <- data.frame(leaders, cdss)
df_te <- data.frame(te_leader, te_cds)

df$scale <- c(1:length(leaders))
df_te$scale <- c(1:length(leaders))

d <- melt(df, id.var="scale")
d_te <- melt(df_te, id.var="scale")

dat <- data.frame(d, d_te)
colnames(dat) <- c("scale", "variable", "value", "scale_te", "variable_te", "value_te")

p2 <- ggplot(dat, aes(x=value, y=value_te, color=variable)) + geom_point() + xlab("sum fo reactivities") + ylab("TE") + ggtitle("2-4cell")  + scale_x_log10() + scale_y_log10()

#####
#####
df <- data.frame(leaders, cdss)
df_te <- data.frame(te_cds, te_cds)
df$scale <- c(1:length(leaders))
df_te$scale <- c(1:length(leaders))
d <- melt(df, id.var="scale")
d_te <- melt(df_te, id.var="scale")
dat <- data.frame(d, d_te)
colnames(dat) <- c("scale", "variable", "value", "scale_te", "variable_te", "value_te")

p1 <- ggplot(dat, aes(x=value, y=value_te, color=variable)) + geom_point() + xlab("sum fo reactivities") + ylab("TE of CDS") + ggtitle("2-4cell") + scale_x_log10() + scale_y_log10()

# leaders only
ldf <- data.frame(leaders, te_cds)
p1 <- ggplot(ldf, aes(x=leaders, y=te_cds)) + geom_point() + xlab("sum fo reactivities for leaders") + ylab("TE of CDS") + ggtitle("2-4cell")
p1_log <- ggplot(ldf, aes(x=leaders, y=te_cds)) + geom_point() + xlab("sum fo reactivities for leaders") + ylab("TE of CDS") + ggtitle("2-4cell") + scale_x_log10() + scale_y_log10()

###

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
          leaders <- c(leaders, sum(leader))
          cdss <- c(cdss, sum(coding))
          #trailers <- c(trailers, sum(trailer)/len_t)
          te_leader <- c(te_leader, data[tr_name, "X02_256cell_lead_te"])
          te_cds <- c(te_cds, data[tr_name, "X02_256cell_gene_te"])
        }
      }
    }
  }
}

df <- data.frame(leaders, cdss)
df_te <- data.frame(te_leader, te_cds)

df$scale <- c(1:length(leaders))
df_te$scale <- c(1:length(leaders))

d <- melt(df, id.var="scale")
d_te <- melt(df_te, id.var="scale")

dat <- data.frame(d, d_te)
colnames(dat) <- c("scale", "variable", "value", "scale_te", "variable_te", "value_te")

p4 <- ggplot(dat, aes(x=value, y=value_te, color=variable)) + geom_point() + xlab("sum fo reactivities") + ylab("TE") + ggtitle("256cell")  + scale_x_log10() + scale_y_log10()

#####
#####
df <- data.frame(leaders, cdss)
df_te <- data.frame(te_cds, te_cds)
df$scale <- c(1:length(leaders))
df_te$scale <- c(1:length(leaders))
d <- melt(df, id.var="scale")
d_te <- melt(df_te, id.var="scale")
dat <- data.frame(d, d_te)
colnames(dat) <- c("scale", "variable", "value", "scale_te", "variable_te", "value_te")

p3 <- ggplot(dat, aes(x=value, y=value_te, color=variable)) + geom_point() + xlab("sum fo reactivities") + ylab("TE of CDS") + ggtitle("256cell") + scale_x_log10() + scale_y_log10()

# leaders only
ldf <- data.frame(leaders, te_cds)
p3 <- ggplot(ldf, aes(x=leaders, y=te_cds)) + geom_point() + xlab("sum fo reactivities for leaders") + ylab("TE of CDS") + ggtitle("256cell")
p3_log <- ggplot(ldf, aes(x=leaders, y=te_cds)) + geom_point() + xlab("sum fo reactivities for leaders") + ylab("TE of CDS") + ggtitle("256cell") + scale_x_log10() + scale_y_log10()


### log
ggsave(file = "/Volumes/USELESS/META/SHAPES/te/hindrance/log/p2.png", plot = p2)
ggsave(file = "/Volumes/USELESS/META/SHAPES/te/hindrance/log/p4.png", plot = p4)

#convert p2.png p4.png -append hindrance_log.png

### leaders vs CDS
ggsave(file = "/Volumes/USELESS/META/SHAPES/te/hindrance/leader_cds/p1.png", plot = p1)
ggsave(file = "/Volumes/USELESS/META/SHAPES/te/hindrance/leader_cds/p1_log.png", plot = p1_log)
ggsave(file = "/Volumes/USELESS/META/SHAPES/te/hindrance/leader_cds/p3.png", plot = p3)
ggsave(file = "/Volumes/USELESS/META/SHAPES/te/hindrance/leader_cds/p3_log.png", plot = p3_log)

#convert p1.png p1_log.png +append 2-4cell.png
#convert p3.png p3_log.png +append 256cell.png
#convert 2-4cell.png 256cell.png -append hindrance.png