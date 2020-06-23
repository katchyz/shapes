library(seqinr)
library(data.table)
library(GeneCycle)
library(ggplot2)

fasta_cds <- read.fasta("/Volumes/USELESS/DATA/fasta/danio_rerio/GRCz10/cds/Danio_rerio.GRCz10.cds.all.fa")
names(fasta_cds) <- sapply(names(fasta_cds), function(x){substr(x,1,18)})

# pick sample
gc <- data.table(a = 2, t = 2, g = 3, c = 3)
gc_cds <- sapply(fasta_cds, function(s){as.vector(sapply(s, function(x){gc[[x]]}))})
gc_cds <- gc_cds[sapply(gc_cds, function(x){length(x) > 150})]
gc_cds <- gc_cds[sapply(gc_cds, function(x){typeof(x) == "double"})]

meta <- Reduce("+", lapply(gc_cds, function(x){x[1:150]/sum(x[1:150])}))

amplitudes <- abs(fft(meta))
amp <- amplitudes[2:(length(amplitudes)/2+1)]
periods <- 1/periodogram(meta, method = "clone")$freq
d <- data.frame(periods, amp)

ggplot(d, aes(x=periods, y=amp)) + geom_line() + ggtitle("GC") + xlim(0, 8)
ggsave(file = "/Volumes/USELESS/META/SHAPES_NEW/periodicity/gc.png")

