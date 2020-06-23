library(GeneCycle)

cds_256 <- subsetByOverlaps(gc_unsmoothed_256, cds)
cds_256_list <- split(cds_256, cds_256$trnames)
cds_256_list <- cds_256_list[lengths(cds_256_list) > 150]
names(cds_256_list) <- sapply(names(cds_256_list), function(x){substr(x,1,18)})

cds_256_list <- cds_256_list[names(cds_256_list) %in% rownames(shape_256)]
te_256 <- shape_256[rownames(shape_256) %in% names(cds_256_list),]

###
top_q <- 0.90
te_256 = within(te_256, {
  bin = ifelse(rownames(te_256) %in% rownames(subset(te_256, te_256$te_cds > quantile(te_256$te_cds, probs = c(top_q)))), "H", "M")
})
bot_q <- 0.10
te_256 = within(te_256, {
  bin = ifelse(rownames(te_256) %in% rownames(subset(te_256, te_256$te_cds < quantile(te_256$te_cds, probs = c(bot_q)))), "L", "M")
})

high <- subset(te_256, bin == "H")
rest <- subset(te_256, bin == "M")
low <- subset(te_256, bin == "L")

high_cds_list <- cds_256_list[names(cds_256_list) %in% rownames(high)]
rest_cds_list <- cds_256_list[names(cds_256_list) %in% rownames(rest)]
low_cds_list <- cds_256_list[names(cds_256_list) %in% rownames(low)]

### HIGH -- get 5'end of CDSs
meta <- lapply(high_cds_list, function(x){(x[1:150]$dtcr)/sum(x[1:150]$dtcr)})
meta <- rapply(meta, f=function(x) ifelse(is.na(x),0,x), how="replace")
meta <- Reduce("+", meta)
amplitudes <- abs(fft(meta))
amp <- amplitudes[2:(length(amplitudes)/2+1)]
periods <- 1/periodogram(meta, method = "clone")$freq
df <- data.frame(periods, amp)
p_high <- ggplot(df, aes(x=periods, y=amp)) + geom_line() + ggtitle("FFT, 256cell, high") + xlim(0, 10)
ggsave(file = "/Volumes/USELESS/META/SHAPES/periodicity/TE_high_256.png", plot = p_high)

### REST -- get 5'end of CDSs
meta <- lapply(rest_cds_list, function(x){(x[1:150]$dtcr)/sum(x[1:150]$dtcr)})
meta <- rapply(meta, f=function(x) ifelse(is.na(x),0,x), how="replace")
meta <- Reduce("+", meta)
amplitudes <- abs(fft(meta))
amp <- amplitudes[2:(length(amplitudes)/2+1)]
periods <- 1/periodogram(meta, method = "clone")$freq
df <- data.frame(periods, amp)
p_rest <- ggplot(df, aes(x=periods, y=amp)) + geom_line() + ggtitle("FFT, 256cell, rest") + xlim(0, 10)
ggsave(file = "/Volumes/USELESS/META/SHAPES/periodicity/TE_rest_256.png", plot = p_rest)


###########
###########
###########
###########
###########

cds_2_4 <- subsetByOverlaps(gc_unsmoothed_2_4, cds)
cds_2_4_list <- split(cds_2_4, cds_2_4$trnames)
cds_2_4_list <- cds_2_4_list[lengths(cds_2_4_list) > 150]
names(cds_2_4_list) <- sapply(names(cds_2_4_list), function(x){substr(x,1,18)})

cds_2_4_list <- cds_2_4_list[names(cds_2_4_list) %in% rownames(shape_2_4)]
te_2_4 <- shape_2_4[rownames(shape_2_4) %in% names(cds_2_4_list),]

###
top_q <- 0.90
te_2_4 = within(te_2_4, {
  bin = ifelse(rownames(te_2_4) %in% rownames(subset(te_2_4, te_2_4$te_cds > quantile(te_2_4$te_cds, probs = c(top_q)))), "H", "M")
})
bot_q <- 0.10
te_2_4 = within(te_2_4, {
  bin = ifelse(rownames(te_2_4) %in% rownames(subset(te_2_4, te_2_4$te_cds < quantile(te_2_4$te_cds, probs = c(bot_q)))), "L", "M")
})

high <- subset(te_2_4, bin == "H")
rest <- subset(te_2_4, bin == "M")
low <- subset(te_2_4, bin == "L")

high_cds_list <- cds_2_4_list[names(cds_2_4_list) %in% rownames(high)]
rest_cds_list <- cds_2_4_list[names(cds_2_4_list) %in% rownames(rest)]
low_cds_list <- cds_2_4_list[names(cds_2_4_list) %in% rownames(low)]

### HIGH -- get 5'end of CDSs
meta <- lapply(high_cds_list, function(x){(x[1:150]$dtcr)/sum(x[1:150]$dtcr)})
meta <- rapply(meta, f=function(x) ifelse(is.na(x),0,x), how="replace")
meta <- Reduce("+", meta)
amplitudes <- abs(fft(meta))
amp <- amplitudes[2:(length(amplitudes)/2+1)]
periods <- 1/periodogram(meta, method = "clone")$freq
df <- data.frame(periods, amp)
p_high <- ggplot(df, aes(x=periods, y=amp)) + geom_line() + ggtitle("FFT, 2_4cell, high") + xlim(0, 10)
ggsave(file = "/Volumes/USELESS/META/SHAPES/periodicity/TE_high_2_4.png", plot = p_high)

### REST -- get 5'end of CDSs
meta <- lapply(rest_cds_list, function(x){(x[1:150]$dtcr)/sum(x[1:150]$dtcr)})
meta <- rapply(meta, f=function(x) ifelse(is.na(x),0,x), how="replace")
meta <- Reduce("+", meta)
amplitudes <- abs(fft(meta))
amp <- amplitudes[2:(length(amplitudes)/2+1)]
periods <- 1/periodogram(meta, method = "clone")$freq
df <- data.frame(periods, amp)
p_rest <- ggplot(df, aes(x=periods, y=amp)) + geom_line() + ggtitle("FFT, 2_4cell, rest") + xlim(0, 10)
ggsave(file = "/Volumes/USELESS/META/SHAPES/periodicity/TE_rest_2_4.png", plot = p_rest)

###### LOW
### LOW -- get 5'end of CDSs
meta <- lapply(low_cds_list, function(x){(x[1:150]$dtcr)/sum(x[1:150]$dtcr)})
meta <- rapply(meta, f=function(x) ifelse(is.na(x),0,x), how="replace")
meta <- Reduce("+", meta)
amplitudes <- abs(fft(meta))
amp <- amplitudes[2:(length(amplitudes)/2+1)]
periods <- 1/periodogram(meta, method = "clone")$freq
df <- data.frame(periods, amp)
p_low <- ggplot(df, aes(x=periods, y=amp)) + geom_line() + ggtitle("FFT, 256cell, low") + xlim(0, 10)
ggsave(file = "/Volumes/USELESS/META/SHAPES/periodicity/TE_low_256.png", plot = p_low)

### LOW -- get 5'end of CDSs
meta <- lapply(low_cds_list, function(x){(x[1:150]$dtcr)/sum(x[1:150]$dtcr)})
meta <- rapply(meta, f=function(x) ifelse(is.na(x),0,x), how="replace")
meta <- Reduce("+", meta)
amplitudes <- abs(fft(meta))
amp <- amplitudes[2:(length(amplitudes)/2+1)]
periods <- 1/periodogram(meta, method = "clone")$freq
df <- data.frame(periods, amp)
p_low <- ggplot(df, aes(x=periods, y=amp)) + geom_line() + ggtitle("FFT, 2-4cell, low") + xlim(0, 10)
ggsave(file = "/Volumes/USELESS/META/SHAPES/periodicity/TE_low_2_4.png", plot = p_low)


#########################
####### RPKM ############
#########################
fpkm_2_4 <- read.csv("/Volumes/USELESS/META/SHAPES/fpkm/Pauli_rna/0_2to4Cell_fpkm.csv")
fpkm_256 <- read.csv("/Volumes/USELESS/META/SHAPES/fpkm/Pauli_rna/1_256Cell_fpkm.csv")

df <- data.frame(fpkm_2_4$rkpm)
colnames(df) <- c("rpkm")
rownames(df) <- fpkm_2_4$tx_id
rpkm_2_4 <- subset(df, rpkm > 0)

df <- data.frame(fpkm_256$rkpm)
colnames(df) <- c("rpkm")
rownames(df) <- fpkm_256$tx_id
rpkm_256 <- subset(df, rpkm > 0)

plot_shape_tb <- function(df_rpkm, cds_list, top_q=0.90, bot_q=0.10) {
  df_rpkm = within(df_rpkm, {
    bin = ifelse(rownames(df_rpkm) %in% rownames(subset(df_rpkm, df_rpkm$rpkm > quantile(df_rpkm$rpkm, probs = c(top_q)))), "H", "M")
  })
  df_rpkm = within(df_rpkm, {
    bin = ifelse(rownames(df_rpkm) %in% rownames(subset(df_rpkm, df_rpkm$rpkm < quantile(df_rpkm$rpkm, probs = c(bot_q)))), "L", df_rpkm$bin)
  })
  cds_list <- cds_list[names(cds_list) %in% rownames(df_rpkm)]
  meta <- lapply(cds_list, function(x){(x[1:150]$dtcr)/sum(x[1:150]$dtcr)})
  meta <- rapply(meta, f=function(x) ifelse(is.na(x),0,x), how="replace")
  ### add up meta for H, M, L
  meta_H <- meta[names(meta) %in% rownames(subset(df_rpkm, bin == "H"))]
  meta_L <- meta[names(meta) %in% rownames(subset(df_rpkm, bin == "L"))]
  meta_M <- meta[names(meta) %in% rownames(subset(df_rpkm, bin == "M"))]
  meta_H <- Reduce("+", meta_H)
  meta_L <- Reduce("+", meta_L)
  meta_M <- Reduce("+", meta_M)
  
  amplitudes_H <- abs(fft(meta_H))
  amp <- amplitudes_H[2:(length(amplitudes_H)/2+1)]
  periods_H <- 1/periodogram(meta_H, method = "clone")$freq
  bin <- rep("H", 75)
  df <- data.frame(periods, amp, bin)
  
  amplitudes_L <- abs(fft(meta_L))
  amp <- amplitudes_L[2:(length(amplitudes_L)/2+1)]
  periods_L <- 1/periodogram(meta_L, method = "clone")$freq
  bin <- rep("L", 75)
  df <- rbind(df, data.frame(periods, amp, bin))
  
  amplitudes_M <- abs(fft(meta_M))
  amp <- amplitudes_M[2:(length(amplitudes_M)/2+1)]
  periods_M <- 1/periodogram(meta_M, method = "clone")$freq
  bin <- rep("M", 75)
  df <- rbind(df, data.frame(periods, amp, bin))
  
  p <- ggplot(df, aes(x=periods, y=amp)) + geom_line() + xlim(1, 8) + facet_wrap(~ bin)
  return(p)
}


p_rpkm_2_4 <- plot_shape_tb(rpkm_2_4, cds_2_4_list, top_q=0.90, bot_q=0.10)
p_rpkm_2_4 <- p_rpkm_2_4 + ggtitle("FFT, RPKM, 2_4cell")
ggsave("/Volumes/USELESS/META/SHAPES/periodicity/rpkm_2_4.png", plot = p_rpkm_2_4)


p_rpkm_256 <- plot_shape_tb(rpkm_256, cds_256_list, top_q=0.90, bot_q=0.10)
p_rpkm_256 <- p_rpkm_256 + ggtitle("FFT, RPKM, 256cell")
ggsave("/Volumes/USELESS/META/SHAPES/periodicity/rpkm_256.png", plot = p_rpkm_256)


###############

data = read.csv("/Users/kasia/Documents/PhD/TE.csv", header = TRUE)
rownames(data) <- data$X.transcript_id

ribo_rna <- data.frame(data$X01_2to4cell_gene_ribo, data$X02_256cell_gene_ribo, data$X01_2to4cell_gene_rna, data$X02_256cell_gene_rna)
rownames(ribo_rna) <- data$X.transcript_id
colnames(ribo_rna) <- c("ribo_2_4", "ribo_256", "rna_2_4", "rna_256")


plot_shape_tb <- function(df_rpkm, cds_list, top_q=0.90, bot_q=0.10, col="ribo_2_4") {
  df_rpkm <- subset(df_rpkm, eval(parse(text=col)) > 0)
  df_rpkm = within(df_rpkm, {
    bin = ifelse(rownames(df_rpkm) %in% rownames(subset(df_rpkm, df_rpkm[[eval(col)]] > quantile(df_rpkm[[eval(col)]], probs = c(top_q)))), "H", "M")
  })
  df_rpkm = within(df_rpkm, {
    bin = ifelse(rownames(df_rpkm) %in% rownames(subset(df_rpkm, df_rpkm[[eval(col)]] < quantile(df_rpkm[[eval(col)]], probs = c(bot_q)))), "L", df_rpkm$bin)
  })
  cds_list <- cds_list[names(cds_list) %in% rownames(df_rpkm)]
  meta <- lapply(cds_list, function(x){(x[1:150]$dtcr)/sum(x[1:150]$dtcr)})
  meta <- rapply(meta, f=function(x) ifelse(is.na(x),0,x), how="replace")
  ### add up meta for H, M, L
  meta_H <- meta[names(meta) %in% rownames(subset(df_rpkm, bin == "H"))]
  meta_L <- meta[names(meta) %in% rownames(subset(df_rpkm, bin == "L"))]
  meta_M <- meta[names(meta) %in% rownames(subset(df_rpkm, bin == "M"))]
  meta_H <- Reduce("+", meta_H)
  meta_L <- Reduce("+", meta_L)
  meta_M <- Reduce("+", meta_M)
  meta_M <- meta_M/8
  
  amplitudes_H <- abs(fft(meta_H))
  amp <- amplitudes_H[2:(length(amplitudes_H)/2+1)]
  periods_H <- 1/periodogram(meta_H, method = "clone")$freq
  bin <- rep("H", 75)
  df <- data.frame(periods, amp, bin)
  
  amplitudes_L <- abs(fft(meta_L))
  amp <- amplitudes_L[2:(length(amplitudes_L)/2+1)]
  periods_L <- 1/periodogram(meta_L, method = "clone")$freq
  bin <- rep("L", 75)
  df <- rbind(df, data.frame(periods, amp, bin))
  
  amplitudes_M <- abs(fft(meta_M))
  amp <- amplitudes_M[2:(length(amplitudes_M)/2+1)]
  periods_M <- 1/periodogram(meta_M, method = "clone")$freq
  bin <- rep("M", 75)
  df <- rbind(df, data.frame(periods, amp, bin))
  
  p <- ggplot(df, aes(x=periods, y=amp)) + geom_line() + xlim(1, 8) + facet_wrap(~ bin)
  return(p)
}


p_ribo_2_4 <- plot_shape_tb(ribo_rna, cds_2_4_list, top_q=0.90, bot_q=0.10, col="ribo_2_4")
p_ribo_2_4 <- p_ribo_2_4 + ggtitle("FFT, 2_4cell, ribo")
ggsave("/Volumes/USELESS/META/SHAPES/periodicity/p_ribo_2_4.png", plot = p_ribo_2_4)

p_ribo_256 <- plot_shape_tb(ribo_rna, cds_256_list, top_q=0.90, bot_q=0.10, col="ribo_256")
p_ribo_256 <- p_ribo_256  + ggtitle("FFT, 256cell, ribo")
ggsave("/Volumes/USELESS/META/SHAPES/periodicity/p_ribo_256.png", plot = p_ribo_256)


p_rna_2_4 <- plot_shape_tb(ribo_rna, cds_2_4_list, top_q=0.90, bot_q=0.10, col="rna_2_4")
p_rna_2_4 <- p_rna_2_4 + ggtitle("FFT, 2_4cell, rna")
ggsave("/Volumes/USELESS/META/SHAPES/periodicity/p_rna_2_4.png", plot = p_rna_2_4)

p_rna_256 <- plot_shape_tb(ribo_rna, cds_256_list, top_q=0.90, bot_q=0.10, col="rna_256")
p_rna_256 <- p_rna_256 + ggtitle("FFT, 256cell, rna")
ggsave("/Volumes/USELESS/META/SHAPES/periodicity/p_rna_256.png", plot = p_rna_256)



#### ZERO RIBO COVERAGE
# 2-4cell
df_rpkm <- ribo_rna
df_rpkm <- subset(df_rpkm, ribo_2_4 == 0)
cds_list <- cds_2_4_list
cds_list <- cds_list[names(cds_list) %in% rownames(df_rpkm)]
meta <- lapply(cds_list, function(x){(x[1:150]$dtcr)/sum(x[1:150]$dtcr)})
meta <- rapply(meta, f=function(x) ifelse(is.na(x),0,x), how="replace")
meta <- Reduce("+", meta)
amplitudes <- abs(fft(meta))
amp <- amplitudes[2:(length(amplitudes)/2+1)]
periods <- 1/periodogram(meta, method = "clone")$freq
df <- data.frame(periods, amp)
p_NOribo_2_4 <- ggplot(df, aes(x=periods, y=amp)) + geom_line() + xlim(1, 8) + ggtitle("FFT, 2_4cell, NO_ribo")
ggsave("/Volumes/USELESS/META/SHAPES/periodicity/p_NOribo_2_4.png", plot = p_NOribo_2_4)

# 256cell
df_rpkm <- ribo_rna
df_rpkm <- subset(df_rpkm, ribo_256 == 0)
cds_list <- cds_256_list
cds_list <- cds_list[names(cds_list) %in% rownames(df_rpkm)]
meta <- lapply(cds_list, function(x){(x[1:150]$dtcr)/sum(x[1:150]$dtcr)})
meta <- rapply(meta, f=function(x) ifelse(is.na(x),0,x), how="replace")
meta <- Reduce("+", meta)
amplitudes <- abs(fft(meta))
amp <- amplitudes[2:(length(amplitudes)/2+1)]
periods <- 1/periodogram(meta, method = "clone")$freq
df <- data.frame(periods, amp)
p_NOribo_256 <- ggplot(df, aes(x=periods, y=amp)) + geom_line() + xlim(1, 8) + ggtitle("FFT, 256cell, NO_ribo")
ggsave("/Volumes/USELESS/META/SHAPES/periodicity/p_NOribo_256.png", plot = p_NOribo_256)

#####

df_rpkm <- ribo_rna
df_rpkm <- subset(df_rpkm, ribo_2_4 > 0)
cds_list <- cds_2_4_list
col = "ribo_2_4"
df_rpkm = within(df_rpkm, {
  bin = ifelse(rownames(df_rpkm) %in% rownames(subset(df_rpkm, df_rpkm[[eval(col)]] > quantile(df_rpkm[[eval(col)]], probs = c(top_q)))), "H", "M")
})
df_rpkm = within(df_rpkm, {
  bin = ifelse(rownames(df_rpkm) %in% rownames(subset(df_rpkm, df_rpkm[[eval(col)]] < quantile(df_rpkm[[eval(col)]], probs = c(bot_q)))), "L", df_rpkm$bin)
})

H <- data.frame(rownames(subset(df_rpkm, bin == "H")))
colnames(H) <- "transcript"
write.csv(H, file = "~/Desktop/H.csv")

L <- data.frame(rownames(subset(df_rpkm, bin == "L")))
colnames(L) <- "transcript"
write.csv(L, file = "~/Desktop/L.csv")

M <- data.frame(rownames(subset(df_rpkm, bin == "M")))
colnames(M) <- "transcript"
write.csv(M, file = "~/Desktop/M.csv")


#######
metaDF <- read.csv("~/Desktop/metaDF.csv")
metaDF$X <- NULL
meta_H <- metaDF$X0
meta_L <- metaDF$X1
meta_M <- metaDF$X2

p <- p + ggtitle('ribo coverage, 256cell')
ggsave("/Volumes/USELESS/META/SHAPES/periodicity/p_riboCov_256.png", plot = p)

p0 <- p0 + ggtitle('ribo coverage, 2-4cell')
ggsave("/Volumes/USELESS/META/SHAPES/periodicity/p_riboCov_2_4.png", plot = p0)



