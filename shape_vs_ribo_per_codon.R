library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)
library(ggplot2)

txdb_can <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/danio_rerio/Danio_rerio.GRCz10.84.gtf", format = "gtf")
cds <- cdsBy(txdb_can, by="tx", use.names=TRUE)
cds <- cds[order(names(cds))]
names(cds) <- sapply(names(cds), function(x){substr(x,1,18)})

txLengths <- transcriptLengths(txdb_can, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
rownames(txLengths) <- txLengths$tx_name
rownames(txLengths) <- sapply(rownames(txLengths), function(x){substr(x,1,18)})
txLengths <- txLengths[order(rownames(txLengths)),]

### Ribo-Seq

############### 2-4cell #################
# read BigWigs
# substitute strand, merge
bw24fw <- import.bw("/Volumes/USELESS/OUT/zebrafish_our_out/GRCz10/2-4Cell_fw.bw")
strand(bw24fw) <- rep("+", length(bw24fw))
bw24rv <- import.bw("/Volumes/USELESS/OUT/zebrafish_our_out/GRCz10/2-4Cell_rv.bw")
strand(bw24rv) <- rep("-", length(bw24rv))
bw24 <- c(bw24fw, bw24rv)

# map to transcripts, subset CDSs
ribo_cds_24 <- mapToTranscripts(bw24, cds, ignore.strand = FALSE)
mcols(ribo_cds_24) <- cbind(mcols(ribo_cds_24), DataFrame(bw24[ribo_cds_24$xHits]))
cds_reads <- sum(score(ribo_cds_24))
ribo_cds_24 <- split(ribo_cds_24, seqnames(ribo_cds_24))
ribo_cds_24 <- ribo_cds_24[order(names(ribo_cds_24))]

# make granges with 0 score of transcript length, one by one
# substitute values
seqlengths(ribo_cds_24)[order(names(seqlengths(ribo_cds_24)))] <- txLengths[rownames(txLengths) %in% names(ribo_cds_24),]$cds_len

# cds
cov <- coverage(ribo_cds_24)
cov <- cov[order(names(cov))]
plus <- ribo_cds_24[strand(ribo_cds_24) == "+"]
plus <- plus[sapply(plus, function(x){length(x) > 0})]
minus <- ribo_cds_24[strand(ribo_cds_24) == "-"]
minus <- minus[sapply(minus, function(x){length(x) > 0})]
# whole vector
# for - strand:rev()
#as.vector(cov[["ENSDART00000162928"]])
cov_plus <- sapply(cov[names(cov) %in% names(plus)], function(x){as.vector(x)})
cov_minus <- sapply(cov[names(cov) %in% names(minus)], function(x){rev(as.vector(x))})
for (name in names(cov_plus)) {
  cov_plus[[name]][cov_plus[[name]] == 1] <- ribo_cds_24[[name]]$score
}
for (name in names(cov_minus)) {
  cov_minus[[name]][cov_minus[[name]] == 1] <- rev(ribo_cds_24[[name]]$score)
}
cov <- c(cov_plus, cov_minus)
cov <- cov[order(names(cov))]
cov_cds_24 <- cov

##
ribo24 <- sapply(cov_cds_24, function(x){as.vector(tapply(x[1:(as.integer(length(x)/3)*3)], rep(1:(length(x)/3), each = 3), sum))})
ribo24 <- sapply(ribo24, function(x){x/sum(x)})

############### 256cell #################
# read BigWigs
# substitute strand, merge
bw256fw <- import.bw("/Volumes/USELESS/OUT/zebrafish_our_out/GRCz10/256Cell_fw.bw")
strand(bw256fw) <- rep("+", length(bw256fw))
bw256rv <- import.bw("/Volumes/USELESS/OUT/zebrafish_our_out/GRCz10/256Cell_rv.bw")
strand(bw256rv) <- rep("-", length(bw256rv))
bw256 <- c(bw256fw, bw256rv)

# map to transcripts, subset CDSs
ribo_cds_256 <- mapToTranscripts(bw256, cds, ignore.strand = FALSE)
mcols(ribo_cds_256) <- cbind(mcols(ribo_cds_256), DataFrame(bw256[ribo_cds_256$xHits]))
cds_reads <- sum(score(ribo_cds_256))
ribo_cds_256 <- split(ribo_cds_256, seqnames(ribo_cds_256))
ribo_cds_256 <- ribo_cds_256[order(names(ribo_cds_256))]

# make granges with 0 score of transcript length, one by one
# substitute values
seqlengths(ribo_cds_256)[order(names(seqlengths(ribo_cds_256)))] <- txLengths[rownames(txLengths) %in% names(ribo_cds_256),]$cds_len

# cds
cov <- coverage(ribo_cds_256)
cov <- cov[order(names(cov))]
plus <- ribo_cds_256[strand(ribo_cds_256) == "+"]
plus <- plus[sapply(plus, function(x){length(x) > 0})]
minus <- ribo_cds_256[strand(ribo_cds_256) == "-"]
minus <- minus[sapply(minus, function(x){length(x) > 0})]
# whole vector
# for - strand:rev()
#as.vector(cov[["ENSDART00000162928"]])
cov_plus <- sapply(cov[names(cov) %in% names(plus)], function(x){as.vector(x)})
cov_minus <- sapply(cov[names(cov) %in% names(minus)], function(x){rev(as.vector(x))})
for (name in names(cov_plus)) {
  cov_plus[[name]][cov_plus[[name]] == 1] <- ribo_cds_256[[name]]$score
}
for (name in names(cov_minus)) {
  cov_minus[[name]][cov_minus[[name]] == 1] <- rev(ribo_cds_256[[name]]$score)
}
cov <- c(cov_plus, cov_minus)
cov <- cov[order(names(cov))]
cov_cds_256 <- cov

##
ribo256 <- sapply(cov_cds_256, function(x){as.vector(tapply(x[1:(as.integer(length(x)/3)*3)], rep(1:(length(x)/3), each = 3), sum))})
ribo256 <- sapply(ribo256, function(x){x/sum(x)})

##############
ribo256 <- ribo256[names(ribo256) %in% names(ribo24)]
f <- fasta[names(fasta) %in% names(ribo256)]
ribo24 <- ribo24[names(ribo24) %in% names(ribo256)]

##### PLOT ########
df <- data.frame(fasta = as.character(unlist(f)), ribo24 = as.numeric(unlist(ribo24)), ribo256 = as.numeric(unlist(ribo256)))
agg_ribo <- aggregate(. ~ fasta, df, mean)
agg_ribo <- agg_ribo[!grepl('n', agg_ribo$fasta),]

cu <- read.csv(file = "/Volumes/USELESS/DATA/fasta/danio_rerio/GRCz10/codon_usage.csv")
cu$codon <- gsub("u", "t", tolower(cu$codon))
cu <- cu[order(cu$codon),]

agg_ribo$aa <- cu$aa

## shape
load(file = "/Volumes/USELESS/META/SHAPES/codon/codon_24_256_MEAN.Rsave") #agg

agg_ribo$shape24 <- agg$sum24
agg_ribo$shape256 <- agg$sum256

###
ggplot(agg_ribo, aes(x=shape24, y=ribo24, label=fasta, colour=aa)) + geom_point() + ggtitle("mean shape VS mean ribo") + geom_text(aes(label=fasta),hjust=0, vjust=0, size=3)
ggsave(file = "/Volumes/USELESS/META/SHAPES/codon/ribo/mean_shape_VS_mean_ribo24.png")

ggplot(agg_ribo, aes(x=shape256, y=ribo256, label=fasta, colour=aa)) + geom_point() + ggtitle("mean shape VS mean ribo") + geom_text(aes(label=fasta),hjust=0, vjust=0, size=3)
ggsave(file = "/Volumes/USELESS/META/SHAPES/codon/ribo/mean_shape_VS_mean_ribo256.png")
