### plot START for 2-4 cell and 256 cell
### unsmoothed, per nt
### same scale! (facet? or set it manually)

# load unsmoothed, subset exons, txlen
# gc_unsmoothed_2_4; gc_unsmoothed_256

load(file = "/Volumes/USELESS/META/SHAPES/SHAPE_FPKM.Rdata")

### get transcript dtabase
txdb_can <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/danio_rerio/zebrafish_GRCz10.gtf", format = "gtf")
exons <- exonsBy(txdb_can, by="tx", use.names=TRUE)
exons <- exons[order(names(exons))]

cds <- cdsBy(txdb_can, by="tx", use.names=TRUE)
cds <- cds[order(names(cds))]

txLengths <- transcriptLengths(txdb_can, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
rownames(txLengths) <- txLengths$tx_name
txLengths <- txLengths[order(rownames(txLengths)),]
rownames(txLengths) <- sapply(rownames(txLengths), function(x){substr(x,1,18)})

exons24 <- subsetByOverlaps(gc_unsmoothed_2_4, exons)
exons24 <- split(exons24, exons24$trnames)
names(exons24) <- sapply(names(exons24), function(x){substr(x,1,18)})

exons256 <- subsetByOverlaps(gc_unsmoothed_256, exons)
exons256 <- split(exons256, exons256$trnames)
names(exons256) <- sapply(names(exons256), function(x){substr(x,1,18)})

# subset on shape fpkm
wes <- rownames(SHAPE_FPKM[SHAPE_FPKM$stage_256cell > 1 & SHAPE_FPKM$stage_2_4cell > 1,])
# subset on utr5 and cds length
txl <- txLengths[rownames(txLengths) %in% wes,]
txl <- rownames(txl[txl$utr5_len > 30 & txl$cds_len > 33,])

# pick the same transcripts in 2-4cell and 256cell
t24 <- exons24[names(exons24) %in% txl]
t256 <- exons256[names(exons256) %in% txl]

t24 <- t24[sapply(t24, function(x){sum(x$dtcr) > 0})]
t256 <- t256[sapply(t256, function(x){sum(x$dtcr) > 0})]

### OTPG
txdb_can <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/danio_rerio/Danio_rerio.GRCz10.84.gtf", format = "gtf")
txLengths <- transcriptLengths(txdb_can, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
rownames(txLengths) <- txLengths$tx_name
txLengths <- txLengths[order(rownames(txLengths)),]
rownames(txLengths) <- sapply(rownames(txLengths), function(x){substr(x,1,18)})
# get longest transcript per gene
txlen <- arrange(txLengths, gene_id, desc(tx_len))
txlen <- txlen[!duplicated(txlen$gene_id),]
# get transcripts with at least 33 on CDS and 30 on utr3
otpg <- sort(txlen[(txlen$cds_len > 33) & (txlen$utr3_len > 30), ]$tx_name)

###
t24 <- t24[names(t24) %in% otpg]
t256 <- t256[names(t256) %in% otpg]

t24 <- t24[names(t24) %in% names(t256)]
t256 <- t256[names(t256) %in% names(t24)]

# plot
meta24 <- Reduce("+", lapply(names(t24), function(x){t24[[x]][(txLengths[x,]$utr5_len-29):(txLengths[x,]$utr5_len+33)]$dtcr / sum(t24[[x]]$dtcr)}))

meta256 <- Reduce("+", lapply(names(t256), function(x){t256[[x]][(txLengths[x,]$utr5_len-29):(txLengths[x,]$utr5_len+33)]$dtcr / sum(t256[[x]]$dtcr)}))

for (name in names(t256)) {
  if (txLengths[name,]$utr5_len+32 > length(t256[[name]])) {
    print(name)
  }
}

t256[["ENSDART00000166917"]] <- NULL
t256[["ENSDART00000167973"]] <- NULL
t256[["ENSDART00000054202"]] <- NULL

#
scale <- c(-30:33)[-31]
df24 <- data.frame(scale,meta24)
p24 <- ggplot(df24, aes(x=scale, y=meta24)) + geom_bar(stat = "identity") + xlab("position from start codon") + ylab("SHAPE") + ggtitle("2-4cell") + ylim(0,13)
ggsave(file = "/Volumes/USELESS/META/SHAPES/start/p24_unsm.png", plot = p24)

df256 <- data.frame(scale,meta256)
p256 <- ggplot(df256, aes(x=scale, y=meta256)) + geom_bar(stat = "identity") + xlab("position from start codon") + ylab("SHAPE") + ggtitle("256cell") + ylim(0,13)
ggsave(file = "/Volumes/USELESS/META/SHAPES/start/p256_unsm.png", plot = p256)

####### smoothed
exons24 <- subsetByOverlaps(gc_dtcr_2_4, exons)
exons24 <- split(exons24, exons24$trnames)
names(exons24) <- sapply(names(exons24), function(x){substr(x,1,18)})

exons256 <- subsetByOverlaps(gc_dtcr_256, exons)
exons256 <- split(exons256, exons256$trnames)
names(exons256) <- sapply(names(exons256), function(x){substr(x,1,18)})

s24 <- exons24[names(exons24) %in% names(t24)]
s256 <- exons256[names(exons256) %in% names(t256)]

# plot
meta24 <- Reduce("+", lapply(names(s24), function(x){s24[[x]][(txLengths[x,]$utr5_len-29):(txLengths[x,]$utr5_len+33)]$dtcr / sum(s24[[x]]$dtcr)}))

meta256 <- Reduce("+", lapply(names(s256), function(x){s256[[x]][(txLengths[x,]$utr5_len-29):(txLengths[x,]$utr5_len+33)]$dtcr / sum(s256[[x]]$dtcr)}))

#
scale <- c(-30:33)[-31]
df24 <- data.frame(scale,meta24)
p24 <- ggplot(df24, aes(x=scale, y=meta24)) + geom_bar(stat = "identity") + xlab("position from start codon") + ylab("SHAPE") + ggtitle("2-4cell")
ggsave(file = "/Volumes/USELESS/META/SHAPES/start/p24_sm.png", plot = p24)

df256 <- data.frame(scale,meta256)
p256 <- ggplot(df256, aes(x=scale, y=meta256)) + geom_bar(stat = "identity") + xlab("position from start codon") + ylab("SHAPE") + ggtitle("256cell")
ggsave(file = "/Volumes/USELESS/META/SHAPES/start/p256_sm.png", plot = p256)


##
load(file = "/Volumes/USELESS/OUT/SHAPES_out/shapeseq_norm_dtcr_2_4.Rsave")
old24 <- shapeseq_norm
load(file = "/Volumes/USELESS/OUT/SHAPES_out/shapeseq_norm_dtcr_256.Rsave")
old256 <- shapeseq_norm

old24 <- split(old24, seqnames(old24))
old256 <- split(old256, seqnames(old256))

names(old24) <- sapply(names(old24), function(x){substr(x,1,18)})
names(old256) <- sapply(names(old256), function(x){substr(x,1,18)})

s24 <- old24[names(old24) %in% names(t24)]
s256 <- old256[names(old256) %in% names(t256)]

u24 <- unlist(s24)
u24[is.na(u24$dtcr)]$dtcr <- 0
u24 <- split(u24, names(u24))

u256 <- unlist(s256)
u256[is.na(u256$dtcr)]$dtcr <- 0
u256 <- split(u256, names(u256))

# plot
meta24 <- Reduce("+", lapply(names(u24), function(x){u24[[x]][(txLengths[x,]$utr5_len-29):(txLengths[x,]$utr5_len+33)]$dtcr / sum(u24[[x]]$dtcr)}))

meta256 <- Reduce("+", lapply(names(u256), function(x){u256[[x]][(txLengths[x,]$utr5_len-29):(txLengths[x,]$utr5_len+33)]$dtcr / sum(u256[[x]]$dtcr)}))


##
load(file = "/Volumes/USELESS/OUT/SHAPES_out/norm_2_4_list.Rsave")
load(file = "/Volumes/USELESS/OUT/SHAPES_out/norm_256_list.Rsave")

names(norm_2_4_list) <- sapply(names(norm_2_4_list), function(x){substr(x,1,18)})
names(norm_256_list) <- sapply(names(norm_256_list), function(x){substr(x,1,18)})

s24 <- norm_2_4_list[names(norm_2_4_list) %in% names(t24)]
s256 <- norm_256_list[names(norm_256_list) %in% names(t256)]


load(file = "/Volumes/USELESS/OUT/SHAPES_out/df_2_4_start.Rsave")
load(file = "/Volumes/USELESS/OUT/SHAPES_out/df_256_start.Rsave")

p24 <- ggplot(df_2_4_start, aes(x=scale, y=meta)) + geom_bar(stat = "identity") + xlab("position from start codon") + ylab("SHAPE") + ggtitle("2-4cell") + ylim(0,55)
ggsave(file = "/Volumes/USELESS/META/SHAPES/start/p24_bump.png", plot = p24)

p256 <- ggplot(df_256_start, aes(x=scale, y=meta)) + geom_bar(stat = "identity") + xlab("position from start codon") + ylab("SHAPE") + ggtitle("256cell") + ylim(0,55)
ggsave(file = "/Volumes/USELESS/META/SHAPES/start/p256_bump.png", plot = p256)


######### get FASTA
library(seqinr)

fasta_cdna <- read.fasta("/Volumes/USELESS/DATA/fasta/danio_rerio/GRCz10/cdna/Danio_rerio.GRCz10.cdna.all.fa")
names(fasta_cdna) <- sapply(names(fasta_cdna), function(x){substr(x,1,18)})

fasta <- fasta_cdna[names(fasta_cdna) %in% names(t24)]

logo <- lapply(names(fasta), function(x){fasta[[x]][(txLengths[x,]$utr5_len-29):(txLengths[x,]$utr5_len+33)]})

cds <- cdsBy(txdb_can, by="tx", use.names=TRUE)
cds <- cds[order(names(cds))]

cds24 <- subsetByOverlaps(gc_unsmoothed_2_4, cds)
cds24 <- split(cds24, cds24$trnames)
names(cds24) <- sapply(names(cds24), function(x){substr(x,1,18)})

cds256 <- subsetByOverlaps(gc_unsmoothed_256, cds)
cds256 <- split(cds256, cds256$trnames)
names(cds256) <- sapply(names(cds256), function(x){substr(x,1,18)})

c24 <- cds24[names(cds24) %in% names(t24)]
c256 <- cds256[names(cds256) %in% names(t256)]

## phasing - normalize by sum_tx

n <- names(c24)
c24 <- lapply(names(c24), function(x){c24[[x]]$dtcr / sum(t24[[x]]$dtcr)})
names(c24) <- n

n <- names(c256)
c256 <- lapply(names(c256), function(x){c256[[x]]$dtcr / sum(t256[[x]]$dtcr)})
names(c256) <- n

phasing24 <- lapply(c24, function(x){c(sum(x[seq(1, length(x), 3)]), sum(x[seq(2, length(x), 3)]), sum(x[seq(3, length(x), 3)]))})
phasing256 <- lapply(c256, function(x){c(sum(x[seq(1, length(x), 3)]), sum(x[seq(2, length(x), 3)]), sum(x[seq(3, length(x), 3)]))})

all24 <- c(sum(sapply(phasing24, function(x){x[1]})), sum(sapply(phasing24, function(x){x[2]})), sum(sapply(phasing24, function(x){x[3]})))
all256 <- c(sum(sapply(phasing256, function(x){x[1]})), sum(sapply(phasing256, function(x){x[2]})), sum(sapply(phasing256, function(x){x[3]})))

orf <- c(1,2,3)
df <- data.frame(shape = all24, orf = orf)
ggplot(df, aes(x=orf, y=shape)) + geom_bar(stat = "identity") + ggtitle("2-4cell")
ggsave(file = "/Volumes/USELESS/META/SHAPES/start/orf24.png")
df <- data.frame(shape = all256, orf = orf)
ggplot(df, aes(x=orf, y=shape)) + geom_bar(stat = "identity") + ggtitle("256cell")
ggsave(file = "/Volumes/USELESS/META/SHAPES/start/orf256.png")

h24max <- data.frame(max = sapply(phasing24, function(x){which.max(x)}))
h24min <- data.frame(min = sapply(phasing24, function(x){which.min(x)}))

ggplot(h24max, aes(max)) + geom_histogram()
ggplot(h24min, aes(min)) + geom_histogram()

####
#### FASTA per codon
fasta_cds <- read.fasta("/Volumes/USELESS/DATA/fasta/danio_rerio/GRCz10/cds/Danio_rerio.GRCz10.cds.all.fa")
names(fasta_cds) <- sapply(names(fasta_cds), function(x){substr(x,1,18)})
fasta <- fasta_cds[names(fasta_cds) %in% names(t24)]

## collapse 3nt
fasta <- fasta[unlist(lapply(fasta, function(x){length(x) %% 3 == 0}))]
fasta <- sapply(fasta, function(x){as.character(tapply(x, rep(1:(length(x)/3), each = 3), paste, collapse = ""))})
## sum 3nt shape

sum24 <- sapply(c24, function(x){as.vector(tapply(x[1:(as.integer(length(x)/3)*3)], rep(1:(length(x)/3), each = 3), sum))})
sum256 <- sapply(c256, function(x){as.vector(tapply(x[1:(as.integer(length(x)/3)*3)], rep(1:(length(x)/3), each = 3), sum))})

sum24 <- sapply(sum24, function(x){x/sum(x)})
sum256 <- sapply(sum256, function(x){x/sum(x)})

## add codon and corresponding value to data frame
i = 0
for (name in names(fasta)) {
  if (length(fasta[[name]]) != length(sum24[[name]])) {
    fasta[[name]] <- NULL
    sum24[[name]] <- NULL
    sum256[[name]] <- NULL
  }
}

sum24 <- sum24[names(sum24) %in% names(fasta)]
sum256 <- sum256[names(sum256) %in% names(fasta)]

#
d <- data.frame(val = c(1,2,3,4), cod = c('A', 'B', 'C', 'A'))
d[] <- lapply(d, function(x) type.convert(as.character(x)))
aggregate(. ~ cod, d, sum)

df <- data.frame(fasta = as.character(unlist(fasta)), sum24 = as.numeric(unlist(sum24)), sum256 = as.numeric(unlist(sum256)))
agg <- aggregate(. ~ fasta, df, sum)

agg <- agg[!grepl('n', agg$fasta),]
save(agg, file = "/Volumes/USELESS/META/SHAPES/codon/codon_24_256.Rsave")

agg <- agg[order(-agg$sum24),]
ggplot(agg, aes(x=fasta, y=sum24)) + geom_bar(stat = "identity")

## codon usage
cu <- read.csv(file = "/Volumes/USELESS/DATA/fasta/danio_rerio/GRCz10/codon_usage.csv")
cu$codon <- gsub("u", "t", tolower(cu$codon))

agg <- agg[order(agg$fasta),]
cu <- cu[order(cu$codon),]

df <- data.frame(codon = cu$codon, freq = cu$freq, shape = agg$sum24)
ggplot(df, aes(x=shape, y=freq)) + geom_point()


####### MEAN !!!!!!!!!!!!!!!!!!!!!
df <- data.frame(fasta = as.character(unlist(fasta)), sum24 = as.numeric(unlist(sum24)), sum256 = as.numeric(unlist(sum256)))
agg <- aggregate(. ~ fasta, df, mean)

agg <- agg[!grepl('n', agg$fasta),]
save(agg, file = "/Volumes/USELESS/META/SHAPES/codon/codon_24_256_MEAN.Rsave")

agg <- agg[order(-agg$sum24),]
ggplot(agg, aes(x=fasta, y=sum24)) + geom_bar(stat = "identity")

agg <- agg[order(agg$fasta),]
cu <- cu[order(cu$codon),]

df <- data.frame(codon = cu$codon, freq = cu$freq, shape = agg$sum24)
ggplot(df, aes(x=shape, y=freq, label=codon)) + geom_point() + ggtitle("mean shape per codon VS frequency in the genome") + geom_text(aes(label=codon),hjust=0, vjust=0, size=3)
ggsave(file = "/Volumes/USELESS/META/SHAPES/codon/mean_shape_per_codon_VS_freq_in_the_genome.png")

ggplot(agg, aes(x=sum24, y=sum256, label=fasta, colour=aa)) + geom_point() + ggtitle("mean shape per codon in 2-4cell vs 256cell") + geom_text(aes(label=fasta),hjust=0, vjust=0, size=3)
ggsave(file = "/Volumes/USELESS/META/SHAPES/codon/shape_per_codon_24cell_vs_256cell.png")

#######
### shape per codon vs relative adaptiveness
w <- read.csv("/Volumes/USELESS/META/SHAPES/w.csv", header = TRUE, row.names = 1)
rownames(w) <- tolower(rownames(w))


df <- data.frame(codon = cu$codon, freq = cu$freq, shape256 = agg$sum256, w = w$w, aa = cu$aa)

ggplot(df, aes(x=shape256, y=freq, label=codon, colour=aa)) + geom_point() + ggtitle("mean shape per codon VS frequency in the genome") + geom_text(aes(label=codon),hjust=0, vjust=0, size=3)
ggsave(file = "/Volumes/USELESS/META/SHAPES/codon/mean_shape_per_codon_256_VS_freq_in_the_genome.png")

ggplot(df, aes(x=shape256, y=w, label=codon, colour=aa)) + geom_point() + ggtitle("mean shape per codon VS relative codon adaptiveness") + geom_text(aes(label=codon),hjust=0, vjust=0, size=3)
ggsave(file = "/Volumes/USELESS/META/SHAPES/codon/mean_shape_per_codon_256_vs_relative_codon_adaptiveness.png")

### TGA??

table(as.character(sapply(fasta, function(x){tail(x,1)})))
tga <- fasta[sapply(fasta, function(x){tail(x,1) == "tga"})]

table(as.character(sapply(tga, function(x){tail(x,2)[1]})))

n <- names(fasta)
save(n, file = "/Users/kasia/Desktop/START_names.Rsave")
# load(file = "/Users/kasia/Desktop/START_names.Rsave")


sum24 <- sapply(cds24, function(x){as.vector(tapply(x$dtcr[1:(as.integer(length(x)/3)*3)], rep(1:(length(x)/3), each = 3), sum))})
sum256 <- sapply(cds256, function(x){as.vector(tapply(x$dtcr[1:(as.integer(length(x)/3)*3)], rep(1:(length(x)/3), each = 3), sum))})

norm24 <- sapply(sum24, function(x){x/sum(x)})
norm256 <- sapply(sum256, function(x){x/sum(x)})

## STOP codon vs SHAPE_FPKM
fasta <- fasta[order(names(fasta))]

names <- n[order(n)]
stop_codon <- as.character(sapply(fasta, function(x){tail(x,1)}))
stop_codon_shape24 <- as.numeric(sapply(norm24, function(x){tail(x,1)}))
stop_codon_shape256 <- as.numeric(sapply(norm256, function(x){tail(x,1)}))

s <- SHAPE_FPKM[rownames(SHAPE_FPKM) %in% names,]

df <- data.frame(row.names = names, stop_codon = stop_codon, stop_codon_shape24 = stop_codon_shape24, stop_codon_shape256 = stop_codon_shape256, shape_fpkm24 = s$stage_2_4cell, shape_fpkm256 = s$stage_256cell)

df <- df[df$stop_codon %in% c('taa', 'tga', 'tag'),]

ggplot(df, aes(x = stop_codon_shape24, y = shape_fpkm24, colour = stop_codon)) + geom_point() + scale_x_log10() + scale_y_log10() + ggtitle("2-4cell") + facet_wrap(~ stop_codon)
ggsave(file = "/Volumes/USELESS/META/SHAPES/codon/stop_codon/stop_codon_shape_vs_SHAPE_FPKM_24.png")

# START codon vs SHAPE_FPKM
start_codon_shape24 <- as.numeric(sapply(norm24, function(x){head(x,1)}))
start_codon_shape256 <- as.numeric(sapply(norm256, function(x){head(x,1)}))

df <- data.frame(row.names = names, stop_codon = stop_codon, stop_codon_shape24 = stop_codon_shape24, stop_codon_shape256 = stop_codon_shape256, shape_fpkm24 = s$stage_2_4cell, shape_fpkm256 = s$stage_256cell, start_codon_shape24 = start_codon_shape24, start_codon_shape256 = start_codon_shape256)

ggplot(df, aes(x = start_codon_shape256, y = shape_fpkm256)) + geom_point() + scale_x_log10() + scale_y_log10() + ggtitle("256cell")
ggsave(file = "/Volumes/USELESS/META/SHAPES/codon/stop_codon/start_codon_shape_vs_SHAPE_FPKM_256.png")

# 2nd last codon vs SHAPE_FPKM (split by stop codon)

second_last_codon <- as.character(sapply(fasta, function(x){head(tail(x,2),1)}))
second_last_codon_shape24 <- as.numeric(sapply(norm24, function(x){head(tail(x,1),1)}))
second_last_codon_shape256 <- as.numeric(sapply(norm256, function(x){head(tail(x,1),1)}))

# slc <- data.frame(second_last_codon = names(sort(table(second_last_codon), decreasing = TRUE)), occurence = as.numeric(sort(table(second_last_codon), decreasing = TRUE)))

df <- data.frame(row.names = names, stop_codon = stop_codon, stop_codon_shape24 = stop_codon_shape24, stop_codon_shape256 = stop_codon_shape256, shape_fpkm24 = s$stage_2_4cell, shape_fpkm256 = s$stage_256cell, start_codon_shape24 = start_codon_shape24, start_codon_shape256 = start_codon_shape256, second_last_codon = second_last_codon, second_last_codon_shape24 = second_last_codon_shape24, second_last_codon_shape256 = second_last_codon_shape256)
df <- df[df$stop_codon %in% c('taa', 'tga', 'tag'),]

ggplot(df, aes(x = second_last_codon_shape256, y = shape_fpkm256, colour = stop_codon)) + geom_point() + scale_x_log10() + scale_y_log10() + ggtitle("256cell") #+ facet_wrap(~ stop_codon)
ggsave(file = "/Volumes/USELESS/META/SHAPES/codon/stop_codon/second_last_codon_shape_vs_SHAPE_FPKM_256.png")

# Nth codons vs SHAPE_FPKM (control)
nth_codon_shape24 <- as.numeric(sapply(norm24, function(x){head(tail(x,5),1)}))
nth_codon_shape256 <- as.numeric(sapply(norm256, function(x){head(tail(x,5),1)}))

df <- data.frame(row.names = names, stop_codon = stop_codon, stop_codon_shape24 = stop_codon_shape24, stop_codon_shape256 = stop_codon_shape256, shape_fpkm24 = s$stage_2_4cell, shape_fpkm256 = s$stage_256cell, start_codon_shape24 = start_codon_shape24, start_codon_shape256 = start_codon_shape256, second_last_codon = second_last_codon, second_last_codon_shape24 = second_last_codon_shape24, second_last_codon_shape256 = second_last_codon_shape256, nth_codon_shape24 = nth_codon_shape24, nth_codon_shape256 = nth_codon_shape256)

ggplot(df, aes(x = nth_codon_shape24, y = shape_fpkm24)) + geom_point() + scale_x_log10() + scale_y_log10() + ggtitle("2-4cell")
ggplot(df, aes(x = nth_codon_shape256, y = shape_fpkm256)) + geom_point() + scale_x_log10() + scale_y_log10() + ggtitle("256cell")

ggsave(file = "/Volumes/USELESS/META/SHAPES/codon/stop_codon/nth_codon_shape_vs_SHAPE_FPKM_24.png")

## mean shape
mean_shape24 <- as.numeric(sapply(norm24, function(x){mean(x)}))
mean_shape256 <- as.numeric(sapply(norm256, function(x){mean(x)}))

df <- data.frame(row.names = names, stop_codon = stop_codon, stop_codon_shape24 = stop_codon_shape24, stop_codon_shape256 = stop_codon_shape256, shape_fpkm24 = s$stage_2_4cell, shape_fpkm256 = s$stage_256cell, start_codon_shape24 = start_codon_shape24, start_codon_shape256 = start_codon_shape256, second_last_codon = second_last_codon, second_last_codon_shape24 = second_last_codon_shape24, second_last_codon_shape256 = second_last_codon_shape256, nth_codon_shape24 = nth_codon_shape24, nth_codon_shape256 = nth_codon_shape256, mean_shape24 = mean_shape24, mean_shape256 = mean_shape256)

ggplot(df, aes(x = mean_shape24, y = shape_fpkm24)) + geom_point() + scale_x_log10() + scale_y_log10() + ggtitle("2-4cell")
ggsave(file = "/Volumes/USELESS/META/SHAPES/codon/stop_codon/mean_shape_CDS_vs_SHAPE_FPKM_24.png")

### SHAPE vs Ribo-Seq on codon basis

# Ribo-seq? ----> Load from whole_tx_ribo.R !!!!!!!!!
# cov_cds_24
# cov_cds_256

ribo24 <- sapply(cov_cds_24, function(x){as.vector(tapply(x[1:(as.integer(length(x)/3)*3)], rep(1:(length(x)/3), each = 3), sum))})
ribo256 <- sapply(cov_cds_256, function(x){as.vector(tapply(x[1:(as.integer(length(x)/3)*3)], rep(1:(length(x)/3), each = 3), sum))})

ribo24 <- sapply(ribo24, function(x){x/sum(x)})
ribo256 <- sapply(ribo256, function(x){x/sum(x)})

#
ribo256 <- ribo256[names(ribo256) %in% names(ribo24)]
f <- fasta[names(fasta) %in% names(ribo256)]
ribo24 <- ribo24[names(ribo24) %in% names(ribo256)]

for (name in names(f)) {
  if (length(f[[name]]) != length(ribo24[[name]])) {
    f[[name]] <- NULL
    ribo24[[name]] <- NULL
    ribo256[[name]] <- NULL
  }
}

i = 0
for (name in names(f)) {
  if (length(f[[name]]) == length(ribo256[[name]])) {
    i <- i + 1
  }
}


df <- data.frame(fasta = as.character(unlist(f)), ribo24 = as.numeric(unlist(ribo24)), ribo256 = as.numeric(unlist(ribo256)))
agg <- aggregate(. ~ fasta, df, sum)


####
## median shape
median_shape24 <- as.numeric(sapply(norm24, function(x){median(x)}))
median_shape256 <- as.numeric(sapply(norm256, function(x){median(x)}))

df <- data.frame(row.names = names, stop_codon = stop_codon, stop_codon_shape24 = stop_codon_shape24, stop_codon_shape256 = stop_codon_shape256, shape_fpkm24 = s$stage_2_4cell, shape_fpkm256 = s$stage_256cell, start_codon_shape24 = start_codon_shape24, start_codon_shape256 = start_codon_shape256, second_last_codon = second_last_codon, second_last_codon_shape24 = second_last_codon_shape24, second_last_codon_shape256 = second_last_codon_shape256, nth_codon_shape24 = nth_codon_shape24, nth_codon_shape256 = nth_codon_shape256, median_shape24 = median_shape24, median_shape256 = median_shape256)

ggplot(df, aes(x = median_shape256, y = shape_fpkm256)) + geom_point() + scale_x_log10() + scale_y_log10() + ggtitle("256cell")
ggsave(file = "/Volumes/USELESS/META/SHAPES/codon/stop_codon/median_shape_CDS_vs_SHAPE_FPKM_256.png")
