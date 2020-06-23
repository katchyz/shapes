#### check for binding sites of RT primer
#### (are they in the transcripts in the 'cloud')

library(Biostrings)
library(seqinr)
library(dplyr)

### RT primer: AGACGTGTGCTCTTCCGATCTNNNNNNNNNNNNNNN

RTprimer <- DNAString("AGACGTGTGCTCTTCCGATCT")
RT <- tolower("AGACGTGTGCTCTTCCGATCT")
#RTprimer_comp <- complement(RTprimer)
#RTprimer_revComp <- reverseComplement(RTprimer)
#RTc <- as.character(RTprimer_comp)
#RTrc <- as.character(RTprimer_revComp)

RTc <- tolower(RTc)
RTrc <- tolower(RTrc)

#as.numeric(stringDist(c(RTc, "ATGACCGTAA")))

### get fasta sequences, subset cloud, run distances
fasta_cdna <- read.fasta("/Volumes/USELESS/DATA/fasta/danio_rerio/GRCz10/cdna/Danio_rerio.GRCz10.cdna.all.fa")
names(fasta_cdna) <- sapply(names(fasta_cdna), function(x){substr(x,1,18)})
fasta_cdna <- sapply(fasta_cdna, function(x){paste(as.character(x), sep="", collapse="")})

RTprimer_aln_width <- sapply(fasta_cdna, function(x){width(pattern(pairwiseAlignment(RT, x, type = "local")))})

r <- data.frame(n = names(RTprimer_aln_width), RTprimer_aln_width = RTprimer_aln_width)

df <- merge(df, r, by="n", all=TRUE)
###save(df, file="/Volumes/USELESS/META/SHAPES_NEW/general/clouds/cloud_df.Rsave")

RTc_aln_width <- sapply(fasta_cdna, function(x){width(pattern(pairwiseAlignment(RTc, x, type = "local")))})
r <- data.frame(n = names(RTc_aln_width), RTc_aln_width = RTc_aln_width)
df <- merge(df, r, by="n", all=TRUE)

RTrc_aln_width <- sapply(fasta_cdna, function(x){width(pattern(pairwiseAlignment(RTrc, x, type = "local")))})
r <- data.frame(n = names(RTrc_aln_width), RTrc_aln_width = RTrc_aln_width)
df <- merge(df, r, by="n", all=TRUE)

##### plot cloud and sample of main on histogram !!!!!!! (either RTc or RTrc > 7)

### plot distributions of cloud and rest
hybr_cloud <- df[df$cloud_cell24 == "cloud" & !is.na(df$cloud_cell24),]$n
hybr_main <- df[(df$shape_cell24_log2ratio_internal > 0) & (df$rna_cell24 > 0),]$n
hybr_main <- hybr_main[!(hybr_main %in% hybr_cloud)]
hybr_main <- hybr_main[complete.cases(hybr_main)]


cld <- data.frame(values = df[df$n %in% hybr_cloud,]$RTc_aln_width,
                  type = rep("cloud", length(df[df$n %in% hybr_cloud,]$RTc_aln_width)))
rst <- data.frame(values = df[df$n %in% hybr_main,]$RTc_aln_width,
                  type = rep("main", length(df[df$n %in% hybr_main,]$RTc_aln_width)))
d <- rbind(cld, sample_n(rst, nrow(cld)))

ggplot(d, aes(values, fill = type)) + geom_histogram(position = "dodge") + ggtitle("2-4Cell")
ggsave(file="/Volumes/USELESS/META/SHAPES_NEW/general/RTprimer/RTc_cell24.png")

#################################

###
cloud <- df[df$cloud_cell24 == "cloud" & !is.na(df$cloud_cell24),]$n
seqs <- fasta_cdna[names(fasta_cdna) %in% cloud]
seqs <- sapply(seqs, function(x){paste(as.character(x), sep="", collapse="")})
main <- df[is.na(df$cloud_cell24) & df$shape_cell24_log2ratio_internal > 0 & df$rna_cell24 > 0,]$n
rest <- fasta_cdna[names(fasta_cdna) %in% main]
rest <- sapply(rest, function(x){paste(as.character(x), sep="", collapse="")})

alns <- sapply(seqs, function(x){width(pattern(pairwiseAlignment(RT, x, type = "local")))})
alns_rest <- sapply(rest, function(x){width(pattern(pairwiseAlignment(RT, x, type = "local")))})

### plot distributions of cloud and rest
cld <- data.frame(values = alns, type = rep("cloud", length(alns)))
rst <- data.frame(values = alns_rest, type = rep("main", length(alns_rest)))
d <- rbind(cld, sample_n(rst, nrow(cld)))

ggplot(d, aes(values, fill = type)) + geom_histogram(position = "dodge") + ggtitle("2-4Cell")
ggsave(file="/Volumes/USELESS/META/SHAPES_NEW/general/RTprimer/RTprimer_binding_cell24.png")

a_cell24 <- alns
ar_cell24 <- alns_rest
####

cloud <- df[df$cloud_cell256 == "cloud" & !is.na(df$cloud_cell256),]$n
seqs <- fasta_cdna[names(fasta_cdna) %in% cloud]
seqs <- sapply(seqs, function(x){paste(as.character(x), sep="", collapse="")})
main <- df[is.na(df$cloud_cell256) & df$shape_cell256_log2ratio_internal > 0 & df$rna_cell256 > 0,]$n
rest <- fasta_cdna[names(fasta_cdna) %in% main]
rest <- sapply(rest, function(x){paste(as.character(x), sep="", collapse="")})

alns <- sapply(seqs, function(x){width(pattern(pairwiseAlignment(RT, x, type = "local")))})
alns_rest <- sapply(rest, function(x){width(pattern(pairwiseAlignment(RT, x, type = "local")))})

### plot distributions of cloud and rest
cld <- data.frame(values = alns, type = rep("cloud", length(alns)))
rst <- data.frame(values = alns_rest, type = rep("main", length(alns_rest)))
d <- rbind(cld, sample_n(rst, nrow(cld)))

ggplot(d, aes(values, fill = type)) + geom_histogram(position = "dodge") + ggtitle("256Cell")
ggsave(file="/Volumes/USELESS/META/SHAPES_NEW/general/RTprimer/RTprimer_binding_cell256.png")

a_cell256 <- alns
ar_cell256 <- alns_rest


#################
#################
### filter out all with alignment width > 7 from the cloud; make GO analysis on the rest ???

ggplot(df, aes(x = log2(shape_cell24_log2ratio_internal), y = log2(rna_cell24), colour=RTprimer_aln_width)) + geom_point(size=0.1) + scale_colour_gradient2(low = "white", mid = "white", high = "blue", midpoint = 7)
ggsave(file="/Volumes/USELESS/META/SHAPES_NEW/general/RTprimer/RTprimer_binding_cell24_cloud.png")

ggplot(df, aes(x = log2(shape_cell24_log2ratio_internal), y = log2(rna_cell24))) + 
  geom_point(aes(colour = cut(RTprimer_aln_width, c(0, 7, Inf))), size = 0.1) + 
  scale_colour_manual(values = c("[0,7]" = "white", "(7,Inf]" = "blue"))


##########

### plot distributions of cloud and rest
hybr_cloud <- df[df$cloud_oblong == "cloud" & !is.na(df$cloud_oblong),]$n
hybr_main <- df[(df$shape_oblong_log2ratio_internal > 0) & (df$rna_3_5h > 0),]$n
hybr_main <- hybr_main[!(hybr_main %in% hybr_cloud)]
hybr_main <- hybr_main[complete.cases(hybr_main)]


cld <- data.frame(values = df[df$n %in% hybr_cloud,]$RTrc_aln_width,
                  type = rep("cloud", length(df[df$n %in% hybr_cloud,]$RTrc_aln_width)))
rst <- data.frame(values = df[df$n %in% hybr_main,]$RTrc_aln_width,
                  type = rep("main", length(df[df$n %in% hybr_main,]$RTrc_aln_width)))
d <- rbind(cld, sample_n(rst, nrow(cld)))

ggplot(d, aes(values, fill = type)) + geom_histogram(position = "dodge") + ggtitle("oblong")
ggsave(file="/Volumes/USELESS/META/SHAPES_NEW/general/RTprimer/RTrc_oblong.png")
