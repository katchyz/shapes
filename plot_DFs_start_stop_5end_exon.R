path_norm <- "/Volumes/USELESS/META/SHAPES/normalized"

dfs <- c("df_256_5end.Rsave", "df_256_exon.Rsave", "df_256_start.Rsave", "df_256_stop.Rsave", "df_2_4_5end.Rsave",  "df_2_4_exon.Rsave",  "df_2_4_start.Rsave", "df_2_4_stop.Rsave")

for (d in dfs) {
  load(file.path(path_norm, d))
}

########
# 5'end
ggplot(df_2_4_5end, aes(x=scale, y=meta)) + geom_bar(stat = "identity") + xlab("position from 5'end of transcript") + ylab("sum(reactivities), norm") + ggtitle("2-4cell")
ggsave(file = "/Volumes/USELESS/META/SHAPES/normalized/p_2_4_5end.png")

ggplot(df_256_5end, aes(x=scale, y=meta)) + geom_bar(stat = "identity") + xlab("position from 5'end of transcript") + ylab("sum(reactivities), norm") + ggtitle("256cell")
ggsave(file = "/Volumes/USELESS/META/SHAPES/normalized/p_256_5end.png")

########
# start
ggplot(df_2_4_start, aes(x=scale, y=meta)) + geom_bar(stat = "identity") + xlab("position from start codon") + ylab("sum(reactivities), norm") + ggtitle("2-4cell")
ggsave(file = "/Volumes/USELESS/META/SHAPES/normalized/p_2_4_start.png")

ggplot(df_256_start, aes(x=scale, y=meta)) + geom_bar(stat = "identity") + xlab("position from start codon") + ylab("sum(reactivities), norm") + ggtitle("256cell")
ggsave(file = "/Volumes/USELESS/META/SHAPES/normalized/p_256_start.png")

########
# stop
ggplot(df_2_4_stop, aes(x=scale, y=meta)) + geom_bar(stat = "identity") + xlab("position from stop codon") + ylab("sum(reactivities), norm") + ggtitle("2-4cell")
ggsave(file = "/Volumes/USELESS/META/SHAPES/normalized/p_2_4_stop.png")

ggplot(df_256_stop, aes(x=scale, y=meta)) + geom_bar(stat = "identity") + xlab("position from stop codon") + ylab("sum(reactivities), norm") + ggtitle("256cell")
ggsave(file = "/Volumes/USELESS/META/SHAPES/normalized/p_256_stop.png")

########
# exon
ggplot(df_2_4_exon, aes(x=scale, y=meta)) + geom_bar(stat = "identity") + xlab("position from splice site") + ylab("sum(reactivities), norm") + ggtitle("2-4cell")
ggsave(file = "/Volumes/USELESS/META/SHAPES/normalized/p_2_4_exon.png")

ggplot(df_256_exon, aes(x=scale, y=meta)) + geom_bar(stat = "identity") + xlab("position from splice site") + ylab("sum(reactivities), norm") + ggtitle("256cell")
ggsave(file = "/Volumes/USELESS/META/SHAPES/normalized/p_256_exon.png")
