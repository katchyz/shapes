library(seqinr)
library(GenomicFeatures)
library(ggplot2)
library(GeneCycle)
library(car)
txdb_can <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/danio_rerio/Danio_rerio.GRCz10.84.gtf", format = "gtf")
txLengths <- transcriptLengths(txdb_can, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
rownames(txLengths) <- txLengths$tx_name
txLengths <- txLengths[order(rownames(txLengths)),]
txLengths <- txLengths[txLengths$cds_len > 150,]

fasta_vienna <- read.fasta("/Volumes/USELESS/DATA/fasta/danio_rerio/GRCz10/cdna/structures_only.txt")
fasta_vienna <- sapply(fasta_vienna, function(x){recode(getSequence(x),"'.'=1;'('=0;')'=0")})
fasta_vienna <- fasta_vienna[names(fasta_vienna) %in% rownames(txLengths)]

txLengths <- txLengths[rownames(txLengths) %in% names(fasta_vienna),]
utr5_len <- setNames(txLengths$utr5_len, rownames(txLengths))

subset150 <- mapply(function(x,y){x[(y+1):(y+150)]}, fasta_vienna, utr5_len)
meta <- rowSums(subset150, na.rm=TRUE)

meta <- f150[sapply(f150, function(x){length(x) == 150})]
meta <- meta[sapply(meta, function(x){sum(x$TC) > 0})]
m <- Reduce("+", lapply(meta, function(x){x$TC/sum(x$TC)}))

amplitudes <- abs(fft(m))
amp <- amplitudes[2:(length(amplitudes)/2+1)]
periods <- 1/periodogram(m, method = "clone")$freq
df <- data.frame(periods, amp)
ggplot(df, aes(x=periods, y=amp)) + geom_line() + ggtitle("cell1 in vitro no control") + xlim(0, 8)


ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/periodicity/tc_cell1_invitro_et2.png")

