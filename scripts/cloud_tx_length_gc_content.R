### length and GC content of the cloud -> any different?
library(seqinr)
library(ggplot2)

txl <- data.frame(n = txLengths$tx_name, tx_len = txLengths$tx_len,
                  cds_len = txLengths$cds_len, utr5_len = txLengths$utr5_len, utr3_len = txLengths$utr3_len)

df <- merge(df, txl, by="n", all=TRUE)
###save(df, file="/Volumes/USELESS/META/SHAPES_NEW/general/clouds/cloud_df.Rsave")

hybr_cloud <- df[df$cloud_cell24 == "cloud" & !is.na(df$cloud_cell24),]$n
hybr_main <- df[(df$shape_cell24_log2ratio_internal > 0) & (df$rna_cell24 > 0),]$n
hybr_main <- hybr_main[!(hybr_main %in% hybr_cloud)]
hybr_main <- hybr_main[complete.cases(hybr_main)]

cld <- data.frame(values = df[df$n %in% hybr_cloud,]$tx_len,
                  type = rep("cloud", length(df[df$n %in% hybr_cloud,]$tx_len)))
rst <- data.frame(values = df[df$n %in% hybr_main,]$tx_len,
                  type = rep("main", length(df[df$n %in% hybr_main,]$tx_len)))
#d <- rbind(cld, sample_n(rst, nrow(cld)))
d <- rbind(cld, rst)

ggplot(d, aes(values, fill=type)) + geom_density(alpha=.3) + ggtitle("2-4cell") + scale_x_log10() + xlab("tx_lengths")
ggsave(file="/Volumes/USELESS/META/SHAPES_NEW/general/lengths/txlen_cell24.png")

### plot per utr5, cds, utr3
hybr_cloud <- df[df$cloud_oblong == "cloud" & !is.na(df$cloud_oblong),]$n
hybr_main <- df[(df$shape_oblong_log2ratio_internal > 0) & (df$rna_3_5h > 0),]$n
hybr_main <- hybr_main[!(hybr_main %in% hybr_cloud)]
hybr_main <- hybr_main[complete.cases(hybr_main)]

cld <- data.frame(values = df[df$n %in% hybr_cloud,]$utr5_len,
                  type = rep("cloud", length(df[df$n %in% hybr_cloud,]$utr5_len)))
rst <- data.frame(values = df[df$n %in% hybr_main,]$utr5_len,
                  type = rep("main", length(df[df$n %in% hybr_main,]$utr5_len)))
#d <- rbind(cld, sample_n(rst, nrow(cld)))
d <- rbind(cld, rst)

ggplot(d, aes(values, fill=type)) + geom_density(alpha=.3) + ggtitle("oblong") + scale_x_log10() + xlab("utr5_lengths")
ggsave(file="/Volumes/USELESS/META/SHAPES_NEW/general/lengths/utr5len_oblong.png")

###########
####GC content
fasta_cdna <- read.fasta("/Volumes/USELESS/DATA/fasta/danio_rerio/GRCz10/cdna/Danio_rerio.GRCz10.cdna.all.fa")
names(fasta_cdna) <- sapply(names(fasta_cdna), function(x){substr(x,1,18)})

gc <- sapply(fasta_cdna, function(x){GC(x)})
gc <- data.frame(n = names(gc), gc = gc)
df <- merge(df, gc, by="n", all=TRUE)

##
hybr_cloud <- df[df$cloud_cell1_invitro == "cloud" & !is.na(df$cloud_cell1_invitro),]$n
hybr_main <- df[(df$shapes_cell1_invitro_TC_internal > 0) & (df$rna_cell1 > 0),]$n
hybr_main <- hybr_main[!(hybr_main %in% hybr_cloud)]
hybr_main <- hybr_main[complete.cases(hybr_main)]

cld <- data.frame(values = df[df$n %in% hybr_cloud,]$gc,
                  type = rep("cloud", length(df[df$n %in% hybr_cloud,]$gc)))
rst <- data.frame(values = df[df$n %in% hybr_main,]$gc,
                  type = rep("main", length(df[df$n %in% hybr_main,]$gc)))
#d <- rbind(cld, sample_n(rst, nrow(cld)))
d <- rbind(cld, rst)

ggplot(d, aes(values, fill=type)) + geom_density(alpha=.3) + ggtitle("1cell invitro") + scale_x_log10() + xlab("GC content")
ggsave(file="/Volumes/USELESS/META/SHAPES_NEW/general/gc/gc_cell1_invitro.png")


########
## Vienna folding
rnafold <- read.csv(file="/Volumes/USELESS/DATA/fasta/danio_rerio/GRCz10/cdna/min_free_energies.txt", header = T, sep = "\t")
rnafold <- rnafold[order(rnafold$n),]
rownames(rnafold) <- NULL

hybr_cloud <- df[df$cloud_oblong_CHX == "cloud" & !is.na(df$cloud_oblong_CHX),]$n
hybr_main <- df[(df$shape_oblong_CHX_log2ratio_internal > 0) & (df$rna_3_5h > 0),]$n
hybr_main <- hybr_main[!(hybr_main %in% hybr_cloud)]
hybr_main <- hybr_main[complete.cases(hybr_main)]
cld <- data.frame(values = df[df$n %in% hybr_cloud,]$minfe,
                  type = rep("cloud", length(df[df$n %in% hybr_cloud,]$minfe)))
rst <- data.frame(values = df[df$n %in% hybr_main,]$minfe,
                  type = rep("main", length(df[df$n %in% hybr_main,]$minfe)))
d <- rbind(cld, rst)

ggplot(d, aes(values, fill=type)) + geom_density(alpha=.3) + ggtitle("oblong CHX") + xlab("min folding E")
ggsave(file="/Volumes/USELESS/META/SHAPES_NEW/general/minfe/minfe_oblong_CHX.png")

