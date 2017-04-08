### shape vs TE per codon
# (scatter plot, per codon??? divide ribo coverage per RNA, add up codons)
library(GenomicFeatures)
library(seqinr)
library(ggplot2)

## cds
load(file = "/Volumes/USELESS/META/SHAPES/riboseq_2_4.Rdata")
load(file = "/Volumes/USELESS/META/SHAPES/riboseq_256.Rdata")

load(file = "/Volumes/USELESS/META/SHAPES/FPKM_24.Rdata")
load(file = "/Volumes/USELESS/META/SHAPES/FPKM_256.Rdata")

# pick those with RNA_FPKM > 1
fpkm24 <- FPKM_24[FPKM_24$exons_RNA_fpkm > 1,]
fpkm256 <- FPKM_256[FPKM_256$exons_RNA_fpkm > 1,]

ribo24 <- riboseq_2_4[names(riboseq_2_4) %in% rownames(fpkm24)]
fpkm24 <- fpkm24[rownames(fpkm24) %in% names(ribo24),]
ribo256 <- riboseq_256[names(riboseq_256) %in% rownames(fpkm256)]
fpkm256 <- fpkm256[rownames(fpkm256) %in% names(ribo256),]

ribo24 <- ribo24[order(names(ribo24))]
ribo256 <- ribo256[order(names(ribo256))]

# divide by RNA_fpkm
te24 <- mapply(function(x,y){x/y}, ribo24, fpkm24$exons_RNA_fpkm)
te256 <- mapply(function(x,y){x/y}, ribo256, fpkm256$exons_RNA_fpkm)

# get fasta
#### FASTA per codon
fasta_cds <- read.fasta("/Volumes/USELESS/DATA/fasta/danio_rerio/GRCz10/cds/Danio_rerio.GRCz10.cds.all.fa")
names(fasta_cds) <- sapply(names(fasta_cds), function(x){substr(x,1,18)})

u <- unique(c(names(te24), names(te256)))
fasta <- fasta_cds[names(fasta_cds) %in% u]
## collapse 3nt
fasta <- fasta[unlist(lapply(fasta, function(x){length(x) %% 3 == 0}))]
fasta <- sapply(fasta, function(x){as.character(tapply(x, rep(1:(length(x)/3), each = 3), paste, collapse = ""))})

# add te by codons
# fasta24 & te24
fasta <- fasta[names(fasta) %in% names(te24)]
fasta <- fasta[names(fasta) %in% names(te256)]
te24 <- te24[names(te24) %in% names(fasta)]
te256 <- te256[names(te256) %in% names(fasta)]
## order!!!
fasta <- fasta[order(names(fasta))]
te24 <- te24[order(names(te24))]
te256 <- te256[order(names(te256))]

sum24 <- sapply(te24, function(x){as.vector(tapply(x[1:(as.integer(length(x)/3)*3)], rep(1:(length(x)/3), each = 3), sum))})
sum256 <- sapply(te256, function(x){as.vector(tapply(x[1:(as.integer(length(x)/3)*3)], rep(1:(length(x)/3), each = 3), sum))})
sum24 <- sapply(sum24, function(x){x/sum(x)})
sum256 <- sapply(sum256, function(x){x/sum(x)})

i = 0
for (name in names(fasta)) {
  if (length(fasta[[name]]) != length(sum24[[name]]) | length(fasta[[name]]) != length(sum256[[name]])) {
    print(name)
    fasta[[name]] <- NULL
    sum24[[name]] <- NULL
    sum256[[name]] <- NULL
  }
}

df <- data.frame(fasta = as.character(unlist(fasta)), te24 = as.numeric(unlist(sum24)), te256 = as.numeric(unlist(sum256)))
agg <- aggregate(. ~ fasta, df, mean)

agg <- agg[!grepl('n', agg$fasta),]

agg <- agg[order(-agg$te24),]
ggplot(agg, aes(x=fasta, y=te24)) + geom_bar(stat = "identity")

# load shape CDSs, aggregate
# add shape by codons
# scatter plot


