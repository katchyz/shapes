### stall sites

# 2-4Cell
fc <- file("/Volumes/USELESS/META/SHAPES/gwp_2-4Cell.csv")
gwp_2_4 <- strsplit(readLines(fc), ",")
close(fc)

# 256Cell
fc <- file("/Volumes/USELESS/META/SHAPES/gwp_256Cell.csv")
gwp_256 <- strsplit(readLines(fc), ",")
close(fc)

# transform the list of peaks into sth
names_2_4 <- sapply(gwp_2_4, function(l) l[[1]])
values_2_4 <- sapply(gwp_2_4, function(l) c(l[2:length(l)]))

names_256 <- sapply(gwp_256, function(l) l[[1]])
values_256 <- sapply(gwp_256, function(l) c(l[2:length(l)]))

# subset granges based on transcript name
names(cds_2_4) <- sapply(names(cds_2_4), function(x){substr(x,1,18)})
names(cds_256) <- sapply(names(cds_256), function(x){substr(x,1,18)})

ss_2_4 <- cds_2_4[names(cds_2_4) %in% names_2_4]
ss_256 <- cds_256[names(cds_256) %in% names_256]

# subset names & values based on indices (the ones in ss)

# add up dtcr around stall sites

ss_2_4_sum <- rep(0,73)
for (i in 1:length(names_2_4)) {
  if (names_2_4[i] %in% names(ss_2_4)) {
    for (p in values_2_4[[i]]) {
      # get fragment around peak from ss_2_4
      if (as.numeric(p)*3 > 35) {
        if (as.numeric(p)*3 < (length(ss_2_4[[names_2_4[[i]]]]$dtcr)-38)) {
          start <- as.numeric(p)*3 - 35
          stop <- as.numeric(p)*3 + 37
          ss_2_4_sum <- ss_2_4_sum + ss_2_4[[names_2_4[[i]]]]$dtcr[start:stop]
          #print(as.numeric(p))
        }
      }
    }
  }
}

ss_256_sum <- rep(0,73)
for (i in 1:length(names_256)) {
  if (names_256[i] %in% names(ss_256)) {
    for (p in values_256[[i]]) {
      # get fragment around peak from ss_2_4
      if (as.numeric(p)*3 > 35) {
        if (as.numeric(p)*3 < (length(ss_256[[names_256[[i]]]]$dtcr)-38)) {
          start <- as.numeric(p)*3 - 35 + 1
          stop <- as.numeric(p)*3 + 37 + 1
          if (sum(ss_256[[names_256[[i]]]]$dtcr[start:stop]) > 0) {
            ss_256_sum <- ss_256_sum + (ss_256[[names_256[[i]]]]$dtcr[start:stop] / sum(ss_256[[names_256[[i]]]]$dtcr[start:stop]))
            #print(as.numeric(p))
          }
        }
      }
    }
  }
}

###
scale <- c(-35:38)[-36]

df <- data.frame(ss_2_4_sum, scale)
p_2_4_ss <- ggplot(df, aes(x=scale, y=ss_2_4_sum)) + geom_bar(stat = "identity")
ggsave("/Volumes/USELESS/META/SHAPES/stall_sites/p_2_4_ss2.png", plot = p_2_4_ss)

df <- data.frame(ss_256_sum, scale)
p_256_ss <- ggplot(df, aes(x=scale, y=ss_256_sum)) + geom_bar(stat = "identity")
ggsave("/Volumes/USELESS/META/SHAPES/stall_sites/p_256_ss2.png", plot = p_256_ss)


#### CONTROL

ctrl_256_sum <- rep(0,33)
for (i in 1:length(names_256)) {
  if (names_256[i] %in% names(ss_256)) {
    for (p in values_256[[i]]) {
      # get fragment around peak from ss_2_4
      if (as.numeric(p)*3 > 45) {
        if (as.numeric(p)*3 < (length(ss_256[[names_256[[i]]]]$dtcr)-18)) {
          start <- as.numeric(p)*3 - 45 + 1
          stop <- as.numeric(p)*3 - 13 + 1
          if (sum(ss_256[[names_256[[i]]]]$dtcr[start:stop]) > 0) {
            ctrl_256_sum <- ctrl_256_sum + (ss_256[[names_256[[i]]]]$dtcr[start:stop] / sum(ss_256[[names_256[[i]]]]$dtcr[start:stop]))
            #print(as.numeric(p))
          }
        }
      }
    }
  }
}

dfc <- data.frame(ctrl_256_sum, scale)
p_256_ssc <- ggplot(dfc, aes(x=scale, y=ctrl_256_sum)) + geom_bar(stat = "identity")

ggsave("/Volumes/USELESS/META/SHAPES/stall_sites/p_256_ss.png", plot = p_256_ss)
ggsave("/Volumes/USELESS/META/SHAPES/stall_sites/p_256_ssc.png", plot = p_256_ssc)

p_256_ss <- p_256_ss + ylab("sum of SHAPE reactivities") + ggtitle("STALL SITES") + xlab(NULL)
p_256_ssc <- p_256_ssc + ylab("sum of SHAPE reactivities") + ggtitle("CONTROL") + xlab(NULL)
