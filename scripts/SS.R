### stall sites - avg. reactivities before and after
library(ggplot2)
library(reshape2)

load("/Volumes/USELESS/META/SHAPES/riboseq_2_4.Rdata")
load("/Volumes/USELESS/META/SHAPES/riboseq_256.Rdata")

# cut-off, well expressed genes
we_2_4 <- riboseq_2_4[sapply(riboseq_2_4, function(x){median(unname(tapply(x, (seq_along(x)-1) %/% 3, sum))) > 0})]
we_256 <- riboseq_256[sapply(riboseq_256, function(x){median(unname(tapply(x, (seq_along(x)-1) %/% 3, sum))) > 0})]

cds_we_2_4 <- cds_2_4[names(cds_2_4) %in% names(we_2_4)]
ribo_we_2_4 <- we_2_4[names(we_2_4) %in% names(cds_we_2_4)]
ribo_we_2_4 <- ribo_we_2_4[order(names(ribo_we_2_4))]

cds_we_256 <- cds_256[names(cds_256) %in% names(we_256)]
ribo_we_256 <- we_256[names(we_256) %in% names(cds_we_256)]
ribo_we_256 <- ribo_we_256[order(names(ribo_we_256))]

## dividable by 3
ribo_we_2_4 <- ribo_we_2_4[sapply(ribo_we_2_4, function(x){length(x) %% 3 == 0})]
ribo_we_256 <- ribo_we_256[sapply(ribo_we_256, function(x){length(x) %% 3 == 0})]

cds_we_2_4 <- cds_2_4[names(cds_2_4) %in% names(ribo_we_2_4)]
cds_we_256 <- cds_256[names(cds_256) %in% names(ribo_we_256)]

# z-scores
zscore <- function(x){
  z <- (x-mean(x))/sd(x)
  return(z)
}

# calculate z-scores on CODON coverage, WITHOUT start and stop codons
z_2_4 <- sapply(ribo_we_2_4, function(x){c(0, zscore(unname(tapply(x[4:(length(x)-3)], (seq_along(x[4:(length(x)-3)])-1) %/% 3, sum))), 0)})
z_256 <- sapply(ribo_we_256, function(x){c(0, zscore(unname(tapply(x[4:(length(x)-3)], (seq_along(x[4:(length(x)-3)])-1) %/% 3, sum))), 0)})

# get genes with peaks - and peak position (on CODON) --- to convert to nt: codon*3-2
gwp_2_4 <- sapply(z_2_4, function(x){which(x > 8.0)})
gwp_2_4 <- gwp_2_4[sapply(gwp_2_4, function(x){length(x) > 0})]

gwp_256 <- sapply(z_256, function(x){which(x > 8.0)})
gwp_256 <- gwp_256[sapply(gwp_256, function(x){length(x) > 0})]

one_ss_2_4 <- gwp_2_4[sapply(gwp_2_4, function(x){length(x) == 1})]
before_2_4 <- list()
after_2_4 <- list()
for (name in names(one_ss_2_4)) {
  nt <- one_ss_2_4[[name]]*3-2
  if (nt < (length(cds_2_4[[name]]) - 15)) {
    ### average shape before and after
    before_2_4[name] <- sum(cds_2_4[[name]]$dtcr[1:nt]) / nt
    after_2_4[name] <- sum(cds_2_4[[name]]$dtcr[nt:(length(cds_2_4[[name]]) - nt)]) / (length(cds_2_4[[name]]) - nt)
  }
}

mean(sapply(before_2_4, function(x){x}))
mean(sapply(after_2_4, function(x){x}))

one_ss_256 <- gwp_256[sapply(gwp_256, function(x){length(x) == 1})]
before_256 <- list()
after_256 <- list()
for (name in names(one_ss_256)) {
  nt <- one_ss_256[[name]]*3-2
  if (nt < (length(cds_256[[name]]) - 15)) {
    ### average shape before and after
    before_256[name] <- sum(cds_256[[name]]$dtcr[1:nt]) / nt
    after_256[name] <- sum(cds_256[[name]]$dtcr[nt:(length(cds_256[[name]]) - nt)]) / (length(cds_256[[name]]) - nt)
  }
}

mean(sapply(before_256, function(x){x}))
mean(sapply(after_256, function(x){x}))

lb <- list()
la <- list()
for (name in names(one_ss_256)) {
  nt <- one_ss_256[[name]]*3-2
  if (nt < (length(cds_256[[name]]) - 15)) {
    ### average shape before and after
    lb[name] <- nt
    la[name] <- length(cds_256[[name]]) - nt
  }
}

# data frames for plotting
befaft_2_4 <- data.frame(before = unname(sapply(before_2_4, function(x){x})), after = unname(sapply(after_2_4, function(x){x})))
befaft_256 <- data.frame(before = unname(sapply(before_256, function(x){x})), after = unname(sapply(after_256, function(x){x})))

df <- melt(befaft_2_4)
ggplot(df, aes(x=value, colour=variable)) + geom_density() + xlim(0,0.04) + ggtitle("2-4cell, avg.shape before/after stall sites") + xlab("avg. SHAPE reactivities")
ggsave(file = "/Volumes/USELESS/META/SHAPES/stall_sites/shape_before_after_SS_2_4.png")

df <- melt(befaft_256)
ggplot(df, aes(x=value, colour=variable)) + geom_density() + xlim(0,0.04) + ggtitle("256cell, avg.shape before/after stall sites") + xlab("avg. SHAPE reactivities")
ggsave(file = "/Volumes/USELESS/META/SHAPES/stall_sites/shape_before_after_SS_256.png")

###### control - well expressed, but no stall site (z-scores < 5)
ctr_2_4 <- z_2_4[sapply(z_2_4, function(x){all(x < 6.0)})]
ctr_256 <- z_256[sapply(z_256, function(x){all(x < 6.0)})]

before_ctr_2_4 <- list()
after_ctr_2_4 <- list()
for (name in names(ctr_2_4)) {
  pos <- sample(seq(15:(length(cds_2_4[[name]]) - 15)), 1)
  before_ctr_2_4[name] <- sum(cds_2_4[[name]]$dtcr[1:pos]) / pos
  after_ctr_2_4[name] <- sum(cds_2_4[[name]]$dtcr[pos:(length(cds_2_4[[name]]) - pos)]) / (length(cds_2_4[[name]]) - pos)
}

befaft_ctr_2_4 <- data.frame(before = unname(sapply(before_ctr_2_4, function(x){x})), after = unname(sapply(after_ctr_2_4, function(x){x})))

df <- melt(befaft_ctr_2_4)
ggplot(df, aes(x=value, colour=variable)) + geom_density() + xlim(0,0.04) + ggtitle("2-4cell, avg.shape before/after random position") + xlab("avg. SHAPE reactivities")
ggsave(file = "/Volumes/USELESS/META/SHAPES/stall_sites/shape_before_after_random_pos_2_4.png")

#
before_ctr_256 <- list()
after_ctr_256 <- list()
for (name in names(ctr_256)) {
  pos <- sample(seq(15:(length(cds_256[[name]]) - 15)), 1)
  before_ctr_256[name] <- sum(cds_256[[name]]$dtcr[1:pos]) / pos
  after_ctr_256[name] <- sum(cds_256[[name]]$dtcr[pos:(length(cds_256[[name]]) - pos)]) / (length(cds_256[[name]]) - pos)
}

befaft_ctr_256 <- data.frame(before = unname(sapply(before_ctr_256, function(x){x})), after = unname(sapply(after_ctr_256, function(x){x})))

df <- melt(befaft_ctr_256)
ggplot(df, aes(x=value, colour=variable)) + geom_density() + xlim(0,0.04) + ggtitle("256cell, avg.shape before/after random position") + xlab("avg. SHAPE reactivities")
ggsave(file = "/Volumes/USELESS/META/SHAPES/stall_sites/shape_before_after_random_pos_256.png")

#####
kasior <- t(data.frame(c(123.1, 145.0, 20.5), 
                       c(155.2, 170.1, 12.1), 
                       c(121.3, 135.4, 10.1), 
                       c(155.1, 175.8, 28.5),
                       c(130.2, 140.3, 20.7),
                       c(150.1, 170.9, 10.9)))
colnames(kasior) <- c("start", "end", "value")
