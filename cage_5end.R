### Shape-seq analysis based on Cage peaks

# Cage
fertilized <- read.csv(file = "/Volumes/USELESS/DATA/Cage_Nepal/leaders_zv10_zf_02_fertilized_egg.csv")
cells_512 <- read.csv(file = "/Volumes/USELESS/DATA/Cage_Nepal/leaders_zv10_zf_04_cells_512.csv")

# set threshold of 10 normalized counts
fertilized <- fertilized[fertilized$count_at_highest_peak > 10,]
cells_512 <- cells_512[cells_512$count_at_highest_peak > 10,]

# Shape-seq
load("/Volumes/USELESS/META/SHAPES/gc_unsmoothed_2_4.Rsave")
load("/Volumes/USELESS/META/SHAPES/gc_unsmoothed_256.Rsave")


## make GRanges object with a (50nt?) fragment from highest Cage peak (GRangesList for each tx)
leaders <- GRangesList()
for (i in 1:nrow(fertilized)) {
  row <- fertilized[i,]
  if (row$dir == "+") {
    leaders[[as.character(row$X.gene_id)]] <- GRanges(seqnames = as.character(row$chr), ranges = IRanges(row$highest_peak, row$highest_peak+50), strand = "+")
  } else if (row$dir == "-") {
    leaders[[as.character(row$X.gene_id)]] <- GRanges(seqnames = as.character(row$chr), ranges = IRanges(row$highest_peak-50, row$highest_peak), strand = "-")
  }
}

leaders512 <- GRangesList()
for (i in 1:nrow(cells_512)) {
  row <- cells_512[i,]
  if (row$dir == "+") {
    leaders512[[as.character(row$X.gene_id)]] <- GRanges(seqnames = as.character(row$chr), ranges = IRanges(row$highest_peak, row$highest_peak+50), strand = "+")
  } else if (row$dir == "-") {
    leaders512[[as.character(row$X.gene_id)]] <- GRanges(seqnames = as.character(row$chr), ranges = IRanges(row$highest_peak-50, row$highest_peak), strand = "-")
  }
}

# makeLeaderGRanges <- function(row) {
#   if (row$dir == "+") {
#     leaders[[as.character(row$X.gene_id)]] <- GRanges(seqnames = as.character(row$chr), ranges = IRanges(row$highest_peak, row$highest_peak+50), strand = "+")
#   } else if (row$dir == "-") {
#     leaders[[as.character(row$X.gene_id)]] <- GRanges(seqnames = as.character(row$chr), ranges = IRanges(row$highest_peak-50, row$highest_peak), strand = "-")
#   }
# }
# 
# apply(fertilized, 1, makeLeaderGRanges)

## subset gc_unsmoothed ??????????????? or some other data file ????????????????
l24 <- subsetByOverlaps(gc_unsmoothed_2_4, leaders)
l24 <- split(l24, l24$trnames)
names(l24) <- sapply(names(l24), function(x){substr(x,1,18)})
l24 <- l24[sapply(l24, function(x){length(x) == 51})]

l256 <- subsetByOverlaps(gc_unsmoothed_256, leaders512)
l256 <- split(l256, l256$trnames)
names(l256) <- sapply(names(l256), function(x){substr(x,1,18)})
l256 <- l256[sapply(l256, function(x){length(x) == 51})]

## plot
scale <- seq(1,51)
meta24 <- rep(0,51)
for (i in seq(1, length(l24))) {
  if (sum(is.na(l24[[i]]$dtcr)) == 0) {
    meta24 <- meta24 + l24[[i]]$dtcr
  }
}
meta256 <- rep(0,51)
for (i in seq(1, length(l256))) {
  if (sum(is.na(l256[[i]]$dtcr)) == 0) {
    meta256 <- meta256 + l256[[i]]$dtcr
  }
}

df24 <- data.frame(scale,meta24)
ggplot(df24, aes(x=scale, y=meta24)) + geom_bar(stat = "identity") + xlab("position from 5'end of transcript") + ylab("sum(reactivities)") + ggtitle("2-4cell, dtcr")
ggsave(file = "/Volumes/USELESS/META/SHAPES/cage_5end/cage_5end_24.png")

df256 <- data.frame(scale,meta256)
ggplot(df256, aes(x=scale, y=meta256)) + geom_bar(stat = "identity") + xlab("position from 5'end of transcript") + ylab("sum(reactivities)") + ggtitle("256cell, dtcr")
ggsave(file = "/Volumes/USELESS/META/SHAPES/cage_5end/cage_5end_256.png")

### treated and control
load("/Volumes/USELESS/META/SHAPES/gc_control_comp_2_4.Rsave")
load("/Volumes/USELESS/META/SHAPES/gc_control_comp_256.Rsave")

load("/Volumes/USELESS/META/SHAPES/gc_treated_comp_2_4.Rsave")
load("/Volumes/USELESS/META/SHAPES/gc_treated_comp_256.Rsave")

gc_control_comp_2_4$TC[is.na(gc_control_comp_2_4$TC)] <- 0
gc_control_comp_256$TC[is.na(gc_control_comp_256$TC)] <- 0
gc_treated_comp_2_4$TC[is.na(gc_treated_comp_2_4$TC)] <- 0
gc_treated_comp_256$TC[is.na(gc_treated_comp_256$TC)] <- 0

gc_control_comp_2_4$TC[is.nan(gc_control_comp_2_4$TC)] <- 0
gc_control_comp_256$TC[is.nan(gc_control_comp_256$TC)] <- 0
gc_treated_comp_2_4$TC[is.nan(gc_treated_comp_2_4$TC)] <- 0
gc_treated_comp_256$TC[is.nan(gc_treated_comp_256$TC)] <- 0

###
t24 <- subsetByOverlaps(gc_treated_comp_2_4, leaders)
t24 <- split(t24, t24$trnames)
names(t24) <- sapply(names(t24), function(x){substr(x,1,18)})
t24 <- t24[sapply(t24, function(x){length(x) == 51})]

t256 <- subsetByOverlaps(gc_treated_comp_256, leaders512)
t256 <- split(t256, t256$trnames)
names(t256) <- sapply(names(t256), function(x){substr(x,1,18)})
t256 <- t256[sapply(t256, function(x){length(x) == 51})]

c24 <- subsetByOverlaps(gc_control_comp_2_4, leaders)
c24 <- split(c24, c24$trnames)
names(c24) <- sapply(names(c24), function(x){substr(x,1,18)})
c24 <- c24[sapply(c24, function(x){length(x) == 51})]

c256 <- subsetByOverlaps(gc_control_comp_256, leaders512)
c256 <- split(c256, c256$trnames)
names(c256) <- sapply(names(c256), function(x){substr(x,1,18)})
c256 <- c256[sapply(c256, function(x){length(x) == 51})]


## plot
# treated
scale <- seq(1,51)
meta24tr <- rep(0,51)
for (i in seq(1, length(t24))) {
  if (sum(is.na(t24[[i]]$TC)) == 0) {
    if (sum(t24[[i]]$TC) > 0) {
      meta24tr <- meta24tr + (t24[[i]]$TC / sum(t24[[i]]$TC))
    }
  }
}
meta256tr <- rep(0,51)
for (i in seq(1, length(t256))) {
  if (sum(is.na(t256[[i]]$TC)) == 0) {
    if (sum(t256[[i]]$TC) > 0) {
      meta256tr <- meta256tr + (t256[[i]]$TC / sum(t256[[i]]$TC))
    }
  }
}

df24tr <- data.frame(scale,meta24tr)
ggplot(df24tr, aes(x=scale, y=meta24tr)) + geom_bar(stat = "identity") + xlab("position from 5'end of transcript") + ylab("sum(reactivities)") + ggtitle("2-4cell, treated")
ggsave(file = "/Volumes/USELESS/META/SHAPES/cage_5end/cage_5end_24tr.png")

df256tr <- data.frame(scale,meta256tr)
ggplot(df256tr, aes(x=scale, y=meta256tr)) + geom_bar(stat = "identity") + xlab("position from 5'end of transcript") + ylab("sum(reactivities)") + ggtitle("256cell, treated")
ggsave(file = "/Volumes/USELESS/META/SHAPES/cage_5end/cage_5end_256tr.png")


# control
scale <- seq(1,51)
meta24ct <- rep(0,51)
for (i in seq(1, length(c24))) {
  if (sum(is.na(c24[[i]]$TC)) == 0) {
    if (sum(c24[[i]]$TC) > 0) {
      meta24ct <- meta24ct + (c24[[i]]$TC / sum(c24[[i]]$TC))
    }
  }
}
meta256ct <- rep(0,51)
for (i in seq(1, length(c256))) {
  if (sum(is.na(c256[[i]]$TC)) == 0) {
    if (sum(c256[[i]]$TC) > 0) {
      meta256ct <- meta256ct + (c256[[i]]$TC / sum(c256[[i]]$TC))
    }
  }
}

df24ct <- data.frame(scale,meta24ct)
ggplot(df24ct, aes(x=scale, y=meta24ct)) + geom_bar(stat = "identity") + xlab("position from 5'end of transcript") + ylab("sum(reactivities)") + ggtitle("2-4cell, control")
ggsave(file = "/Volumes/USELESS/META/SHAPES/cage_5end/cage_5end_24ct.png")

df256ct <- data.frame(scale,meta256ct)
ggplot(df256ct, aes(x=scale, y=meta256ct)) + geom_bar(stat = "identity") + xlab("position from 5'end of transcript") + ylab("sum(reactivities)") + ggtitle("256cell, control")
ggsave(file = "/Volumes/USELESS/META/SHAPES/cage_5end/cage_5end_256ct.png")


###
cc <- subsetByOverlaps(gc_control_comp_256, lead)
cc <- split(cc, cc$trnames)
names(cc) <- sapply(names(cc), function(x){substr(x,1,18)})
cc <- cc[sapply(cc, function(x){length(x) == 51})]

met <- rep(0,51)
for (i in seq(1, length(cc))) {
  if (sum(is.na(cc[[i]]$TC)) == 0) {
    if (sum(cc[[i]]$TC) > 0) {
      met <- met + (cc[[i]]$TC / sum(cc[[i]]$TC))
    }
  }
}

sca <- seq(0,50)
df <- data.frame(sca,met)
ggplot(df, aes(x=sca, y=met)) + geom_bar(stat = "identity") + xlab("position from 5'end of transcript") + ylab("sum(reactivities)") + ggtitle("256cell, control")
ggsave(file = "/Volumes/USELESS/META/SHAPES/cage_5end/peak0.png")
