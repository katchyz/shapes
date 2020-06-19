### process single-end data, after trimming untemplated nucleotides, DMSO + NAI (SHAPE)
### TC counts, log2ratio
library(GenomicFeatures)

rtgu_treated <- get(load(file = "/export/valenfs/data/raw_data/SHAPE/final_results/save_rt_granges/norm1_oblong_CHX_invivo_NAI.Rsave"))
rtgu_control <- get(load(file = "/export/valenfs/data/raw_data/SHAPE/final_results/save_rt_granges/norm1_oblong_CHX_invivo_DMSO.Rsave"))
rm(rt_granges_unique)
save_file <- "/export/valenfs/data/raw_data/SHAPE/final_results/save_rt_granges/norm2_oblong_CHX_invivo.Rsave"

rtgu_treated <- rtgu_treated[start(rtgu_treated) != 0]
rtgu_control <- rtgu_control[start(rtgu_control) != 0]
# RT position and coverage 

#### TREATED
plus <- rtgu_treated[strand(rtgu_treated) == "+"]
minus <- rtgu_treated[strand(rtgu_treated) == "-"]

plus@ranges <- IRanges(start = start(plus), end = start(plus), width = rep(1, length(plus)))
minus@ranges <- IRanges(start = end(minus), end = end(minus), width = rep(1, length(minus)))

### TC
plus_TC <- coverage(plus, weight=plus$unique_barcodes)
plus_TC <- plus_TC[sapply(plus_TC, function(x){length(x) > 0})]
minus_TC <- coverage(minus, weight=minus$unique_barcodes)
minus_TC <- minus_TC[sapply(minus_TC, function(x){length(x) > 0})]

###### change plus into unique, add plus_TC > 0 #############################
plus <- unique(plus)
pTC <- as.numeric(unlist(plus_TC))
pTC <- pTC[pTC > 0]
plus$TC <- pTC

minus <- unique(minus)
mTC <- as.numeric(unlist(minus_TC))
mTC <- mTC[mTC > 0]
minus$TC <- mTC

treated_euc <- c(plus, minus)

#### CONTROL
plus <- rtgu_control[strand(rtgu_control) == "+"]
minus <- rtgu_control[strand(rtgu_control) == "-"]

plus@ranges <- IRanges(start = start(plus), end = start(plus), width = rep(1, length(plus)))
minus@ranges <- IRanges(start = end(minus), end = end(minus), width = rep(1, length(minus)))

### TC
plus_TC <- coverage(plus, weight=plus$unique_barcodes)
plus_TC <- plus_TC[sapply(plus_TC, function(x){length(x) > 0})]
minus_TC <- coverage(minus, weight=minus$unique_barcodes)
minus_TC <- minus_TC[sapply(minus_TC, function(x){length(x) > 0})]

###### change plus into unique, add plus_TC > 0 #############################
plus <- unique(plus)
pTC <- as.numeric(unlist(plus_TC))
pTC <- pTC[pTC > 0]
plus$TC <- pTC

minus <- unique(minus)
mTC <- as.numeric(unlist(minus_TC))
mTC <- mTC[mTC > 0]
minus$TC <- mTC

control_euc <- c(plus, minus)


#####################
###### map to transcripts
txdb_can <- makeTxDbFromGFF("/Home/ii/katchyz/DATA/genomes/GTF/Danio_rerio.GRCz10.84.gtf", format = "gtf")
exons <- exonsBy(txdb_can, by="tx", use.names=TRUE)

#### TREATED
treated_tx <- mapToTranscripts(treated_euc, exons, ignore.strand = FALSE)
treated_tx$TC <- treated_euc[treated_tx$xHits]$TC
##
treated_TC <- coverage(treated_tx, weight=treated_tx$TC) ### RleList - TC for each tx
treated_TC <- treated_TC[order(names(treated_TC))]

exons_treated <- exons[names(exons) %in% names(treated_TC)]
exons_treated <- exons_treated[order(names(exons_treated))]

f <- function(x,y){Rle(c(x, rep(0,(sum(y@ranges@width) - length(x)))))} ## padding with zeros to tx length
treated_TC <- mapply(f, treated_TC, exons_treated)


####### get both control and treated in one granges

g_treated <- function(x,y){GRanges(seqnames = rep(as.character(seqnames(y))[1], length(x)),
                                   IRanges(start = 1:length(x), width = 1),
                                   strand = rep(as.character(strand(y))[1], length(x)),
                                   TC.treated = as.vector(x))} # converting Rle into GRanges


GR_treated <- GRangesList(mapply(g_treated, treated_TC, exons_treated))

################
#### CONTROL
control_tx <- mapToTranscripts(control_euc, exons, ignore.strand = FALSE)
control_tx$TC <- control_euc[control_tx$xHits]$TC
##
control_TC <- coverage(control_tx, weight=control_tx$TC) ### RleList - TC for each tx
control_TC <- control_TC[order(names(control_TC))]

exons_control <- exons[names(exons) %in% names(control_TC)]
exons_control <- exons_control[order(names(exons_control))]

f <- function(x,y){Rle(c(x, rep(0,(sum(y@ranges@width) - length(x)))))} ## padding with zeros to tx length
control_TC <- mapply(f, control_TC, exons_control)


####### get both control and treated in one granges

g_control <- function(x,y){GRanges(seqnames = rep(as.character(seqnames(y))[1], length(x)),
                                   IRanges(start = 1:length(x), width = 1),
                                   strand = rep(as.character(strand(y))[1], length(x)),
                                   TC.control = as.vector(x))} # converting Rle into GRanges


GR_control <- GRangesList(mapply(g_control, control_TC, exons_control))



############ join treated and control
GR_treated <- GR_treated[names(GR_treated) %in% names(GR_control)]
GR_control <- GR_control[names(GR_control) %in% names(GR_treated)]

GR <- unlist(GR_control)
GR$TC.treated <- unlist(GR_treated)$TC.treated

## log2ratio
P <- 1
GR$log2ratio <- log2(GR$TC.treated+P) - log2(GR$TC.control+P)
# ReLU
GR$log2ratio[GR$log2ratio < 0] <- 0

GR <- split(GR, names(GR))

save(GR, file = save_file)


