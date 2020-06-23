library(rtracklayer)
library(GenomicFeatures)

load(file="/Users/kasia/Downloads/gc_cell1.Rsave")
txdb_can <- makeTxDbFromGFF("/Home/ii/katchyz/DATA/genomes/GTF/Danio_rerio.GRCz10.84.gtf", format = "gtf")
exons <- exonsBy(txdb_can, by="tx", use.names=TRUE)
names(exons) <- sapply(names(exons), function(x){substr(x,1,18)})

u <- unlist(shape_cell1)
gc_cell1 <- mapFromTranscripts(u, exons, ignore.strand=FALSE)
gc_cell1$TC.control <- u$TC.control
gc_cell1$TC.treated <- u$TC.treated
gc_cell1$log2ratio <- u$log2ratio
save(gc_cell1, file="gc_cell1.Rsave")

## get chromosome sizes, assign seqlengths
#chr_sizes <- read.csv(file="/Volumes/USELESS/DATA/genomes/chrom_sizes_GRCz10_84.txt", header=T, sep="\t")
#chr366 <- chr_sizes[chr_sizes$chr %in% names(seqlengths(gc_cell1)),]

#seqlen <- chr366$length
#names(seqlen) <- as.character(chr366$chr)
#seqlen <- seqlen[order(match(names(seqlen), levels(seqnames(gc_cell1))))]

#seqlengths(gc_cell1) <- seqlen


gc_cell1$score <- gc_cell1$log2ratio
export.bedGraph(gc_cell1, "/Users/kasia/Downloads/test.bedGraph")

######
gc_cell1$score <- gc_cell1$log2ratio
export.bedGraph(gc_cell1, "cell1.bedGraph")

#########
gc_cell1$score <- gc_cell1$TC.treated
plus <- gc_cell1[strand(gc_cell1) == "+"]
minus <- gc_cell1[strand(gc_cell1) == "-"]
export.bedGraph(plus, "cell1_fw.bedGraph")
export.bedGraph(minus, "cell1_rv.bedGraph")

gc_cell24$score <- gc_cell24$TC.treated
plus <- gc_cell24[strand(gc_cell24) == "+"]
minus <- gc_cell24[strand(gc_cell24) == "-"]
export.bedGraph(plus, "cell24_fw.bedGraph")
export.bedGraph(minus, "cell24_rv.bedGraph")

gc_cell256$score <- gc_cell256$TC.treated
plus <- gc_cell256[strand(gc_cell256) == "+"]
minus <- gc_cell256[strand(gc_cell256) == "-"]
export.bedGraph(plus, "cell256_fw.bedGraph")
export.bedGraph(minus, "cell256_rv.bedGraph")

gc_cell1K$score <- gc_cell1K$TC.treated
plus <- gc_cell1K[strand(gc_cell1K) == "+"]
minus <- gc_cell1K[strand(gc_cell1K) == "-"]
export.bedGraph(plus, "cell1K_fw.bedGraph")
export.bedGraph(minus, "cell1K_rv.bedGraph")

gc_oblong$score <- gc_oblong$TC.treated
plus <- gc_oblong[strand(gc_oblong) == "+"]
minus <- gc_oblong[strand(gc_oblong) == "-"]
export.bedGraph(plus, "oblong_fw.bedGraph")
export.bedGraph(minus, "oblong_rv.bedGraph")

gc_oblong_CHX$score <- gc_oblong_CHX$TC.treated
plus <- gc_oblong_CHX[strand(gc_oblong_CHX) == "+"]
minus <- gc_oblong_CHX[strand(gc_oblong_CHX) == "-"]
export.bedGraph(plus, "oblong_CHX_fw.bedGraph")
export.bedGraph(minus, "oblong_CHX_rv.bedGraph")

gc_cell1_invitro$TC <- u$TC
gc_cell1_invitro$score <- gc_cell1_invitro$TC
plus <- gc_cell1_invitro[strand(gc_cell1_invitro) == "+"]
minus <- gc_cell1_invitro[strand(gc_cell1_invitro) == "-"]
export.bedGraph(plus, "cell1_invitro_fw.bedGraph")
export.bedGraph(minus, "cell1_invitro_rv.bedGraph")

gc_cell24_TC$TC <- u$TC
gc_cell24_TC$score <- gc_cell24_TC$TC
plus <- gc_cell24_TC[strand(gc_cell24_TC) == "+"]
minus <- gc_cell24_TC[strand(gc_cell24_TC) == "-"]
export.bedGraph(plus, "cell24_TC_fw.bedGraph")
export.bedGraph(minus, "cell24_TC_rv.bedGraph")


