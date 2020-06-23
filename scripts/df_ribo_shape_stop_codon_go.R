### table: tx, stop_codon, shape24, shape256, ribo24, ribo256 (all sum over stop codon, relative to the whole tx), GOterm(?)

library(seqinr)
library(GenomicFeatures)

fasta_cds <- read.fasta("/Volumes/USELESS/DATA/fasta/danio_rerio/GRCz10/cds/Danio_rerio.GRCz10.cds.all.fa")
names(fasta_cds) <- sapply(names(fasta_cds), function(x){substr(x,1,18)})

fasta_cds <- fasta_cds[unlist(lapply(fasta_cds, function(x){length(x) %% 3 == 0}))]

# names of tx
taa <- names(fasta_cds[unlist(lapply(fasta_cds, function(x){paste(x[(length(x)-2):length(x)], collapse = "") == "taa"}))])
tga <- names(fasta_cds[unlist(lapply(fasta_cds, function(x){paste(x[(length(x)-2):length(x)], collapse = "") == "tga"}))])
tag <- names(fasta_cds[unlist(lapply(fasta_cds, function(x){paste(x[(length(x)-2):length(x)], collapse = "") == "tag"}))])

# load shape-seq and ribo-seq, subset CDS
# save stop codon value (sum of last three) in the table

## ribo-seq
load(file = "/Volumes/USELESS/META/SHAPES/riboseq_2_4.Rdata")
load(file = "/Volumes/USELESS/META/SHAPES/riboseq_256.Rdata")

txdb_can <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/danio_rerio/zebrafish_GRCz10.gtf", format = "gtf")
cds <- cdsBy(txdb_can, by="tx", use.names=TRUE)

## shape-seq
cds24 <- subsetByOverlaps(gc_unsmoothed_2_4, cds)
cds24 <- split(cds24, cds24$trnames)
names(cds24) <- sapply(names(cds24), function(x){substr(x,1,18)})

cds256 <- subsetByOverlaps(gc_unsmoothed_256, cds)
cds256 <- split(cds256, cds256$trnames)
names(cds256) <- sapply(names(cds256), function(x){substr(x,1,18)})

## GO
go <- read.csv("/Users/kasia/Documents/PhD/matrix_go_fq.csv", header = TRUE)
rownames(go) <- go$X.transcript_id

# taa, tag, tga, cds24, cds256, riboseq_2_4, riboseq_256, go

######################### TAA ###########
ribo24_taa <- riboseq_2_4[names(riboseq_2_4) %in% taa]
ribo256_taa <- riboseq_256[names(riboseq_256) %in% taa]
shape24_taa <- cds24[names(cds24) %in% taa]
shape256_taa <- cds256[names(cds256) %in% taa]
go_taa <- go[rownames(go) %in% taa,]

## merge(a, b, by = "n", all = TRUE)
r24_taa <- data.frame(ribo24_stop_codon = sapply(ribo24_taa, function(x){sum(x[(length(x)-2):length(x)])/sum(x)}))
r24_taa$n <- rownames(r24_taa)
r256_taa <- data.frame(ribo256_stop_codon = sapply(ribo256_taa, function(x){sum(x[(length(x)-2):length(x)])/sum(x)}))
r256_taa$n <- rownames(r256_taa)

df_taa <- merge(r24_taa, r256_taa, by = "n", all = TRUE)

shape24_taa <- shape24_taa[sapply(shape24_taa, function(x){length(x) > 10})]
shape24_taa <- shape24_taa[sapply(shape24_taa, function(x){sum(x$dtcr) > 0})]
s24_taa <- data.frame(shape24_stop_codon = sapply(shape24_taa, function(x){sum(x[(length(x)-2):length(x)]$dtcr)/sum(x$dtcr)}))

shape256_taa <- shape256_taa[sapply(shape256_taa, function(x){length(x) > 10})]
shape256_taa <- shape256_taa[sapply(shape256_taa, function(x){sum(x$dtcr) > 0})]
s256_taa <- data.frame(shape256_stop_codon = sapply(shape256_taa, function(x){sum(x[(length(x)-2):length(x)]$dtcr)/sum(x$dtcr)}))

s24_taa$n <- rownames(s24_taa)
s256_taa$n <- rownames(s256_taa)

df_taa <- merge(df_taa, s24_taa, by = "n", all = TRUE)
df_taa <- merge(df_taa, s256_taa, by = "n", all = TRUE)

df_taa$stop_codon <- rep("taa", nrow(df_taa))

go_taa <- go[rownames(go) %in% df_taa$n,]
gotaa <- data.frame(go = go_taa$GO_Term_Name, n = rownames(go_taa))

df_taa <- merge(df_taa, gotaa, by = "n", all = TRUE)

######################### TAG ###########
ribo24_tag <- riboseq_2_4[names(riboseq_2_4) %in% tag]
ribo256_tag <- riboseq_256[names(riboseq_256) %in% tag]
shape24_tag <- cds24[names(cds24) %in% tag]
shape256_tag <- cds256[names(cds256) %in% tag]
go_tag <- go[rownames(go) %in% tag,]

## merge(a, b, by = "n", all = TRUE)
r24_tag <- data.frame(ribo24_stop_codon = sapply(ribo24_tag, function(x){sum(x[(length(x)-2):length(x)])/sum(x)}))
r24_tag$n <- rownames(r24_tag)
r256_tag <- data.frame(ribo256_stop_codon = sapply(ribo256_tag, function(x){sum(x[(length(x)-2):length(x)])/sum(x)}))
r256_tag$n <- rownames(r256_tag)

df_tag <- merge(r24_tag, r256_tag, by = "n", all = TRUE)

shape24_tag <- shape24_tag[sapply(shape24_tag, function(x){length(x) > 10})]
shape24_tag <- shape24_tag[sapply(shape24_tag, function(x){sum(x$dtcr) > 0})]
s24_tag <- data.frame(shape24_stop_codon = sapply(shape24_tag, function(x){sum(x[(length(x)-2):length(x)]$dtcr)/sum(x$dtcr)}))

shape256_tag <- shape256_tag[sapply(shape256_tag, function(x){length(x) > 10})]
shape256_tag <- shape256_tag[sapply(shape256_tag, function(x){sum(x$dtcr) > 0})]
s256_tag <- data.frame(shape256_stop_codon = sapply(shape256_tag, function(x){sum(x[(length(x)-2):length(x)]$dtcr)/sum(x$dtcr)}))

s24_tag$n <- rownames(s24_tag)
s256_tag$n <- rownames(s256_tag)

df_tag <- merge(df_tag, s24_tag, by = "n", all = TRUE)
df_tag <- merge(df_tag, s256_tag, by = "n", all = TRUE)

df_tag$stop_codon <- rep("tag", nrow(df_tag))

go_tag <- go[rownames(go) %in% df_tag$n,]
gotag <- data.frame(go = go_tag$GO_Term_Name, n = rownames(go_tag))

df_tag <- merge(df_tag, gotag, by = "n", all = TRUE)

######################### TGA ###########
ribo24_tga <- riboseq_2_4[names(riboseq_2_4) %in% tga]
ribo256_tga <- riboseq_256[names(riboseq_256) %in% tga]
shape24_tga <- cds24[names(cds24) %in% tga]
shape256_tga <- cds256[names(cds256) %in% tga]
go_tga <- go[rownames(go) %in% tga,]

## merge(a, b, by = "n", all = TRUE)
r24_tga <- data.frame(ribo24_stop_codon = sapply(ribo24_tga, function(x){sum(x[(length(x)-2):length(x)])/sum(x)}))
r24_tga$n <- rownames(r24_tga)
r256_tga <- data.frame(ribo256_stop_codon = sapply(ribo256_tga, function(x){sum(x[(length(x)-2):length(x)])/sum(x)}))
r256_tga$n <- rownames(r256_tga)

df_tga <- merge(r24_tga, r256_tga, by = "n", all = TRUE)

shape24_tga <- shape24_tga[sapply(shape24_tga, function(x){length(x) > 10})]
shape24_tga <- shape24_tga[sapply(shape24_tga, function(x){sum(x$dtcr) > 0})]
s24_tga <- data.frame(shape24_stop_codon = sapply(shape24_tga, function(x){sum(x[(length(x)-2):length(x)]$dtcr)/sum(x$dtcr)}))

shape256_tga <- shape256_tga[sapply(shape256_tga, function(x){length(x) > 10})]
shape256_tga <- shape256_tga[sapply(shape256_tga, function(x){sum(x$dtcr) > 0})]
s256_tga <- data.frame(shape256_stop_codon = sapply(shape256_tga, function(x){sum(x[(length(x)-2):length(x)]$dtcr)/sum(x$dtcr)}))

s24_tga$n <- rownames(s24_tga)
s256_tga$n <- rownames(s256_tga)

df_tga <- merge(df_tga, s24_tga, by = "n", all = TRUE)
df_tga <- merge(df_tga, s256_tga, by = "n", all = TRUE)

df_tga$stop_codon <- rep("tga", nrow(df_tga))

go_tga <- go[rownames(go) %in% df_tga$n,]
gotga <- data.frame(go = go_tga$GO_Term_Name, n = rownames(go_tga))

df_tga <- merge(df_tga, gotga, by = "n", all = TRUE)

####################
rownames(df_taa) <- df_taa$n
rownames(df_tag) <- df_tag$n
rownames(df_tga) <- df_tga$n

ribo_shape_stop_codon_go <- rbind(df_taa, df_tag, df_tga)

ribo_shape_stop_codon_go$n <- NULL

write.csv(ribo_shape_stop_codon_go, file = "/Volumes/USELESS/META/SHAPES/ribo_shape_stop_codon_go.csv")
