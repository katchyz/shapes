# CAI vs ribo-seq FPKM

TE_fpkm <- read.csv("/Volumes/USELESS/OUT/fpkm/Chew_TE_matrix.csv")
rownames(TE_fpkm) <- TE_fpkm$X.transcript_id

ribo_fpkm <- data.frame(TE_fpkm$ribo_0_2.4Cell.wig, TE_fpkm$ribo_1_256Cell.wig)
colnames(ribo_fpkm) <- c("ribo_2_4", "ribo_256")
rownames(ribo_fpkm) <- rownames(TE_fpkm)

caiS <- cai[rownames(cai) %in% rownames(ribo_fpkm),]
ribo_fpkm <- ribo_fpkm[rownames(ribo_fpkm) %in% rownames(caiS),]

ribo_fpkm$cai <- caiS$CAI

p_2_4 <- ggplot(ribo_fpkm, aes(x = ribo_2_4, y = cai)) + geom_point() + xlab("ribo_fpkm") + ylab("CAI") + ggtitle("2-4cell") + scale_y_log10() + scale_x_log10()

ggsave("/Volumes/USELESS/META/SHAPES/CAI/ribo_fpkm_CAI_2_4.png", plot = p_2_4)

p_256 <- ggplot(ribo_fpkm, aes(x = ribo_256, y = cai)) + geom_point() + xlab("ribo_fpkm") + ylab("CAI") + ggtitle("256cell") + scale_y_log10() + scale_x_log10()

ggsave("/Volumes/USELESS/META/SHAPES/CAI/ribo_fpkm_CAI_256.png", plot = p_256)


################
## phasing, binned by ribo

###
ribo_fpkm <- subset(ribo_fpkm, ribo_256 > 0)

top_q <- 0.90
ribo_fpkm = within(ribo_fpkm, {
  bin = ifelse(rownames(ribo_fpkm) %in% rownames(subset(ribo_fpkm, ribo_fpkm$ribo_256 > quantile(ribo_fpkm$ribo_256, probs = c(top_q)))), "H", "M")
})
bot_q <- 0.10
ribo_fpkm = within(ribo_fpkm, {
  bin = ifelse(rownames(ribo_fpkm) %in% rownames(subset(ribo_fpkm, ribo_fpkm$ribo_256 < quantile(ribo_fpkm$ribo_256, probs = c(bot_q)))), "L", ribo_fpkm$bin)
})

high <- subset(ribo_fpkm, bin == "H")
rest <- subset(ribo_fpkm, bin == "M")
low <- subset(ribo_fpkm, bin == "L")

#########
## get subsets of shape-seq based on high/rest/low, calculate ORF score

# cds_256_list
cds_256_list <- cds_256
names(cds_256_list) <- sapply(names(cds_256_list), function(x){substr(x,1,18)})

# ORF score
# H
h <- cds_256_list[names(cds_256_list) %in% rownames(high)]
orf1_h <- sum(sapply(h, function(x){sum(x[seq(1, length(x), 3)]$dtcr)}))
orf2_h <- sum(sapply(h, function(x){sum(x[seq(2, length(x), 3)]$dtcr)}))
orf3_h <- sum(sapply(h, function(x){sum(x[seq(3, length(x), 3)]$dtcr)}))

# H - normalized
h <- cds_256_list[names(cds_256_list) %in% rownames(high)]
orf1_h <- sum(sapply(h, function(x){sum(x[seq(1, length(x), 3)]$dtcr)}))
orf2_h <- sum(sapply(h, function(x){sum(x[seq(2, length(x), 3)]$dtcr)}))
orf3_h <- sum(sapply(h, function(x){sum(x[seq(3, length(x), 3)]$dtcr)}))

# M
m <- cds_256_list[names(cds_256_list) %in% rownames(rest)]
orf1_m <- sum(sapply(m, function(x){sum(x[seq(1, length(x), 3)]$dtcr)}))
orf2_m <- sum(sapply(m, function(x){sum(x[seq(2, length(x), 3)]$dtcr)}))
orf3_m <- sum(sapply(m, function(x){sum(x[seq(3, length(x), 3)]$dtcr)}))

# L
l <- cds_256_list[names(cds_256_list) %in% rownames(low)]
orf1_l <- sum(sapply(l, function(x){sum(x[seq(1, length(x), 3)]$dtcr)}))
orf2_l <- sum(sapply(l, function(x){sum(x[seq(2, length(x), 3)]$dtcr)}))
orf3_l <- sum(sapply(l, function(x){sum(x[seq(3, length(x), 3)]$dtcr)}))

# ALL
orf1 <- sum(sapply(cds_256_list, function(x){sum(x[seq(1, length(x), 3)]$dtcr)}))
orf2 <- sum(sapply(cds_256_list, function(x){sum(x[seq(2, length(x), 3)]$dtcr)}))
orf3 <- sum(sapply(cds_256_list, function(x){sum(x[seq(3, length(x), 3)]$dtcr)}))

### PLOT
orf <- c("1", "2", "3")
sum_reactivities_h <- c(orf1_h, orf2_h, orf3_h)
df_h <- data.frame(orf, sum_reactivities_h)
ph <- ggplot(df_h, aes(x=orf, y=sum_reactivities_h)) + geom_bar(stat = "identity")
ggsave("/Volumes/USELESS/META/SHAPES/periodicity/phasing_Hribo.png", plot = ph)

sum_reactivities_m <- c(orf1_m, orf2_m, orf3_m)
df_m <- data.frame(orf, sum_reactivities_m)
pm <- ggplot(df_m, aes(x=orf, y=sum_reactivities_m)) + geom_bar(stat = "identity")
ggsave("/Volumes/USELESS/META/SHAPES/periodicity/phasing_Mribo.png", plot = pm)

sum_reactivities_l <- c(orf1_l, orf2_l, orf3_l)
df_l <- data.frame(orf, sum_reactivities_l)
pl <- ggplot(df_l, aes(x=orf, y=sum_reactivities_l)) + geom_bar(stat = "identity")
ggsave("/Volumes/USELESS/META/SHAPES/periodicity/phasing_Lribo.png", plot = pl)

sum_reactivities <- c(orf1, orf2, orf3)
df <- data.frame(orf, sum_reactivities)
p <- ggplot(df, aes(x=orf, y=sum_reactivities)) + geom_bar(stat = "identity")
ggsave("/Volumes/USELESS/META/SHAPES/periodicity/phasing_ribo.png", plot = p)

######## CONTROL and TREATED #########
orf1 <- sum(sapply(ctrl_2_4, function(x){sum(x[seq(1, length(x), 3)]$TC)}))
orf2 <- sum(sapply(ctrl_2_4, function(x){sum(x[seq(2, length(x), 3)]$TC)}))
orf3 <- sum(sapply(ctrl_2_4, function(x){sum(x[seq(3, length(x), 3)]$TC)}))
orf <- c("1", "2", "3")
sum_reactivities <- c(orf1, orf2, orf3)
df <- data.frame(orf, sum_reactivities)
p_ctrl_2_4 <- ggplot(df, aes(x=orf, y=sum_reactivities)) + geom_bar(stat = "identity")
ggsave("/Volumes/USELESS/META/SHAPES/periodicity/orf_CTRL_2_4.png", plot = p_ctrl_2_4)

orf1 <- sum(sapply(ctrl_256, function(x){sum(x[seq(1, length(x), 3)]$TC)}))
orf2 <- sum(sapply(ctrl_256, function(x){sum(x[seq(2, length(x), 3)]$TC)}))
orf3 <- sum(sapply(ctrl_256, function(x){sum(x[seq(3, length(x), 3)]$TC)}))
orf <- c("1", "2", "3")
sum_reactivities <- c(orf1, orf2, orf3)
df <- data.frame(orf, sum_reactivities)
p_ctrl_256 <- ggplot(df, aes(x=orf, y=sum_reactivities)) + geom_bar(stat = "identity")
ggsave("/Volumes/USELESS/META/SHAPES/periodicity/orf_CTRL_256.png", plot = p_ctrl_256)

orf1 <- sum(sapply(trea_2_4, function(x){sum(x[seq(1, length(x), 3)]$TC)}))
orf2 <- sum(sapply(trea_2_4, function(x){sum(x[seq(2, length(x), 3)]$TC)}))
orf3 <- sum(sapply(trea_2_4, function(x){sum(x[seq(3, length(x), 3)]$TC)}))
orf <- c("1", "2", "3")
sum_reactivities <- c(orf1, orf2, orf3)
df <- data.frame(orf, sum_reactivities)
p_trea_2_4 <- ggplot(df, aes(x=orf, y=sum_reactivities)) + geom_bar(stat = "identity")
ggsave("/Volumes/USELESS/META/SHAPES/periodicity/orf_TREA_2_4.png", plot = p_trea_2_4)

orf1 <- sum(sapply(trea_256, function(x){sum(x[seq(1, length(x), 3)]$TC)}))
orf2 <- sum(sapply(trea_256, function(x){sum(x[seq(2, length(x), 3)]$TC)}))
orf3 <- sum(sapply(trea_256, function(x){sum(x[seq(3, length(x), 3)]$TC)}))
orf <- c("1", "2", "3")
sum_reactivities <- c(orf1, orf2, orf3)
df <- data.frame(orf, sum_reactivities)
p_trea_256 <- ggplot(df, aes(x=orf, y=sum_reactivities)) + geom_bar(stat = "identity")
ggsave("/Volumes/USELESS/META/SHAPES/periodicity/orf_TREA_256.png", plot = p_trea_256)

