### process single-end data, after trimming untemplated nucleotides
### TC counts log2ratio
library(GenomicFeatures)
#library(GenomicRanges)

rtgu_control <- get(load(file = "/export/valenfs/data/raw_data/Shape-Seq_Aug2017/save_rt_granges/1Kcell_polyA_vitro_DMSO.Rsave"))
rtgu_treated <- get(load(file = "/export/valenfs/data/raw_data/Shape-Seq_Aug2017/save_rt_granges/1Kcell_polyA_vitro_NAI.Rsave"))
rm(rt_granges_unique)
save_file <- "/export/valenfs/data/raw_data/Shape-Seq_Aug2017/save_rt_granges/GC_1Kcell_polyA_vitro"

# RT position
#### CONTROL
plus <- rtgu_control[strand(rtgu_control) == "+"]
minus <- rtgu_control[strand(rtgu_control) == "-"]
plus@ranges <- IRanges(start = start(plus), end = start(plus), width = rep(1, length(plus)))
minus@ranges <- IRanges(start = end(minus), end = end(minus), width = rep(1, length(minus)))
control_euc <- c(plus, minus)
#### TREATED
plus <- rtgu_treated[strand(rtgu_treated) == "+"]
minus <- rtgu_treated[strand(rtgu_treated) == "-"]
plus@ranges <- IRanges(start = start(plus), end = start(plus), width = rep(1, length(plus)))
minus@ranges <- IRanges(start = end(minus), end = end(minus), width = rep(1, length(minus)))
treated_euc <- c(plus, minus)

txdb_can <- makeTxDbFromGFF("/Home/ii/katchyz/DATA/genomes/GTF/Danio_rerio.GRCz10.84.gtf", format = "gtf")
exons <- exonsBy(txdb_can, by="tx", use.names=TRUE)

#### CONTROL
control_tx <- mapToTranscripts(control_euc, exons, ignore.strand = FALSE)
control_tx$EUC <- control_euc[control_tx$xHits]$unique_barcodes
control_TC <- coverage(control_tx, weight=control_tx$EUC) ### RleList - TC for each tx
control_TC <- control_TC[order(names(control_TC))]
exons_control <- exons[names(exons) %in% names(control_TC)]
exons_control <- exons_control[order(names(exons_control))]
#### TREATED
treated_tx <- mapToTranscripts(treated_euc, exons, ignore.strand = FALSE)
treated_tx$EUC <- treated_euc[treated_tx$xHits]$unique_barcodes
treated_TC <- coverage(treated_tx, weight=treated_tx$EUC) ### RleList - TC for each tx
treated_TC <- treated_TC[order(names(treated_TC))]
exons_treated <- exons[names(exons) %in% names(treated_TC)]
exons_treated <- exons_treated[order(names(exons_treated))]

f <- function(x,y){Rle(c(x, rep(0,(sum(y@ranges@width) - length(x)))))} ## padding with zeros to tx length
control_TC <- mapply(f, control_TC, exons_control)
treated_TC <- mapply(f, treated_TC, exons_treated)


####### get both control and treated in one granges

g_control <- function(x,y){GRanges(seqnames = rep(as.character(seqnames(y))[1], length(x)),
                           IRanges(start = 1:length(x), width = 1),
                           strand = rep(as.character(strand(y))[1], length(x)),
                           TC.control = as.vector(x))} # converting Rle into GRanges

g_treated <- function(x,y){GRanges(seqnames = rep(as.character(seqnames(y))[1], length(x)),
                                   IRanges(start = 1:length(x), width = 1),
                                   strand = rep(as.character(strand(y))[1], length(x)),
                                   TC.treated = as.vector(x))} # converting Rle into GRanges


GR_control <- GRangesList(mapply(g_control, control_TC, exons_control))
GR_treated <- GRangesList(mapply(g_treated, treated_TC, exons_treated))

## put the control and treated together
GR_control <- GR_control[names(GR_control) %in% names(GR_treated)]
GR_treated <- GR_treated[names(GR_treated) %in% names(GR_control)]

GR <- unlist(GR_control)
GR$TC.treated <- unlist(GR_treated)$TC.treated


## normalize (log2ratio)
P <- 1
GR$log2ratio <- log2(GR$TC.treated+P) - log2(GR$TC.control+P)


save(GR, file = save_file)


