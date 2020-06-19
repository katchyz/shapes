### process paired-end data, after trimming untemplated nucleotides
##### FIX !!!!!!! ######## - counts... as in norm2_shape
##########################
### TC counts dTCR
library(GenomicFeatures)
#library(GenomicRanges)

rtgu_control <- get(load(file = "/export/valenfs/data/raw_data/Shape-Seq_Aug2017/save_rt_granges/PAIRED_1Kcell_RZ_DMSO.Rsave"))
rtgu_treated <- get(load(file = "/export/valenfs/data/raw_data/Shape-Seq_Aug2017/save_rt_granges/PAIRED_1Kcell_RZ_NAI.Rsave"))
rm(rt_granges_unique)
save_file <- "/export/valenfs/data/raw_data/Shape-Seq_Aug2017/save_rt_granges/PAIR_GC_1Kcell_RZ.Rsave"

rtgu_control <- rtgu_control[start(rtgu_control) != 0]
rtgu_treated <- rtgu_treated[start(rtgu_treated) != 0]
# RT position and coverage 
#### CONTROL
plus <- rtgu_control[strand(rtgu_control) == "+"]
minus <- rtgu_control[strand(rtgu_control) == "-"]

# split for getting coverage
plus <- split(plus, seqnames(plus))
plus <- plus[sapply(plus, function(x){length(x) > 0})]
minus <- split(minus, seqnames(minus))
minus <- minus[sapply(minus, function(x){length(x) > 0})]
### coverage
cov_plus <- lapply(plus, function(x){coverage(x, weight = x$unique_barcodes)})
cov_plus <- unlist(lapply(cov_plus, function(x){as.list(x)}))
cov_plus <- cov_plus[sapply(cov_plus, function(x){length(x) > 0})]
cov_minus <- lapply(minus, function(x){coverage(x, weight = x$unique_barcodes)})
cov_minus <- unlist(lapply(cov_minus, function(x){as.list(x)}))
cov_minus <- cov_minus[sapply(cov_minus, function(x){length(x) > 0})]
# as.numeric(cov_plus[[1]][start(plus[[1]])])
get_cov_plus <- function(x,y){as.numeric(x[start(y)])}
coverage_plus <- as.numeric(unlist(mapply(get_cov_plus, cov_plus, plus)))
get_cov_minus <- function(x,y){as.numeric(x[end(y)])}
coverage_minus <- as.numeric(unlist(mapply(get_cov_minus, cov_minus, minus)))

plus <- unlist(plus)
minus <- unlist(minus)

plus@ranges <- IRanges(start = start(plus), end = start(plus), width = rep(1, length(plus)))
minus@ranges <- IRanges(start = end(minus), end = end(minus), width = rep(1, length(minus)))
plus$coverage <- coverage_plus
minus$coverage <- coverage_minus

### TC
plus_TC <- coverage(plus, weight=plus$unique_barcodes)
plus_TC <- plus_TC[sapply(plus_TC, function(x){length(x) > 0})]
minus_TC <- coverage(minus, weight=minus$unique_barcodes)
minus_TC <- minus_TC[sapply(minus_TC, function(x){length(x) > 0})]

# split for getting coverage (TC)
plus <- split(plus, seqnames(plus))
plus <- plus[sapply(plus, function(x){length(x) > 0})]
minus <- split(minus, seqnames(minus))
minus <- minus[sapply(minus, function(x){length(x) > 0})]

get_TC_plus <- function(x,y){as.numeric(x[start(y)])}
TC_plus <- as.numeric(unlist(mapply(get_TC_plus, plus_TC, plus)))
get_TC_minus <- function(x,y){as.numeric(x[end(y)])}
TC_minus <- as.numeric(unlist(mapply(get_TC_minus, minus_TC, minus)))

plus <- unlist(plus)
minus <- unlist(minus)

plus$TC <- TC_plus
minus$TC <- TC_minus

control_euc <- c(plus, minus)
control_euc$TCR.control <- control_euc$TC / control_euc$coverage
control_euc$TCR.control[is.na(control_euc$TCR.control)] <- 0

#### TREATED
plus <- rtgu_treated[strand(rtgu_treated) == "+"]
minus <- rtgu_treated[strand(rtgu_treated) == "-"]

# split for getting coverage
plus <- split(plus, seqnames(plus))
plus <- plus[sapply(plus, function(x){length(x) > 0})]
minus <- split(minus, seqnames(minus))
minus <- minus[sapply(minus, function(x){length(x) > 0})]
### coverage
cov_plus <- lapply(plus, function(x){coverage(x, weight = x$unique_barcodes)})
cov_plus <- unlist(lapply(cov_plus, function(x){as.list(x)}))
cov_plus <- cov_plus[sapply(cov_plus, function(x){length(x) > 0})]
cov_minus <- lapply(minus, function(x){coverage(x, weight = x$unique_barcodes)})
cov_minus <- unlist(lapply(cov_minus, function(x){as.list(x)}))
cov_minus <- cov_minus[sapply(cov_minus, function(x){length(x) > 0})]
# as.numeric(cov_plus[[1]][start(plus[[1]])])
get_cov_plus <- function(x,y){as.numeric(x[start(y)])}
coverage_plus <- as.numeric(unlist(mapply(get_cov_plus, cov_plus, plus)))
get_cov_minus <- function(x,y){as.numeric(x[end(y)])}
coverage_minus <- as.numeric(unlist(mapply(get_cov_minus, cov_minus, minus)))

plus <- unlist(plus)
minus <- unlist(minus)

plus@ranges <- IRanges(start = start(plus), end = start(plus), width = rep(1, length(plus)))
minus@ranges <- IRanges(start = end(minus), end = end(minus), width = rep(1, length(minus)))
plus$coverage <- coverage_plus
minus$coverage <- coverage_minus

### TC
plus_TC <- coverage(plus, weight=plus$unique_barcodes)
plus_TC <- plus_TC[sapply(plus_TC, function(x){length(x) > 0})]
minus_TC <- coverage(minus, weight=minus$unique_barcodes)
minus_TC <- minus_TC[sapply(minus_TC, function(x){length(x) > 0})]

###### change plus into unique, add plus_TC > 0 #############################
plus <- unique(plus)
pTC <- as.numeric(unlist(plus_TC))
pTC <- pTC[pTC > 0]
plus$TC.treated <- pTC

minus <- unique(minus)
mTC <- as.numeric(unlist(minus_TC))
mTC <- mTC[mTC > 0]
minus$TC.treated <- mTC

# split for getting coverage (TC)
plus <- split(plus, seqnames(plus))
plus <- plus[sapply(plus, function(x){length(x) > 0})]
minus <- split(minus, seqnames(minus))
minus <- minus[sapply(minus, function(x){length(x) > 0})]

get_TC_plus <- function(x,y){as.numeric(x[start(y)])}
TC_plus <- as.numeric(unlist(mapply(get_TC_plus, plus_TC, plus)))
get_TC_minus <- function(x,y){as.numeric(x[end(y)])}
TC_minus <- as.numeric(unlist(mapply(get_TC_minus, minus_TC, minus)))

plus <- unlist(plus)
minus <- unlist(minus)

plus$TC <- TC_plus
minus$TC <- TC_minus

treated_euc <- c(plus, minus)
treated_euc$TCR.treated <- treated_euc$TC / treated_euc$coverage
treated_euc$TCR.treated[is.na(treated_euc$TCR.treated)] <- 0


###### map to transcripts
txdb_can <- makeTxDbFromGFF("/Home/ii/katchyz/DATA/genomes/GTF/Danio_rerio.GRCz10.84.gtf", format = "gtf")
exons <- exonsBy(txdb_can, by="tx", use.names=TRUE)

#### CONTROL
control_tx <- mapToTranscripts(control_euc, exons, ignore.strand = FALSE)
control_tx$TCR.control <- control_euc[control_tx$xHits]$TCR.control
control_tx$TC <- control_euc[control_tx$xHits]$TC
control_tx$coverage <- control_euc[control_tx$xHits]$coverage

#control_TC <- coverage(control_tx, weight = control_tx$TCR.control) ### RleList - TC for each tx
#control_TC <- control_TC[order(names(control_TC))]

##
control_TCR <- coverage(control_tx, weight = control_tx$TCR.control) ### RleList - TC for each tx
control_TCR <- control_TCR[order(names(control_TCR))]
control_TC <- coverage(control_tx, weight = control_tx$TC) ### RleList - TC for each tx
control_TC <- control_TC[order(names(control_TC))]
control_cov <- coverage(control_tx, weight = control_tx$coverage) ### RleList - TC for each tx
control_cov <- control_cov[order(names(control_cov))]
#
exons_control <- exons[names(exons) %in% names(control_TC)]
exons_control <- exons_control[order(names(exons_control))]

#### TREATED
treated_tx <- mapToTranscripts(treated_euc, exons, ignore.strand = FALSE)
treated_tx$TCR.treated <- treated_euc[treated_tx$xHits]$TCR.treated
treated_tx$TC <- treated_euc[treated_tx$xHits]$TC
treated_tx$coverage <- treated_euc[treated_tx$xHits]$coverage
##
treated_TCR <- coverage(treated_tx, weight=treated_tx$TCR.treated) ### RleList - TC for each tx
treated_TCR <- treated_TCR[order(names(treated_TCR))]
treated_TC <- coverage(treated_tx, weight=treated_tx$TC) ### RleList - TC for each tx
treated_TC <- treated_TC[order(names(treated_TC))]
treated_cov <- coverage(treated_tx, weight=treated_tx$coverage) ### RleList - TC for each tx
treated_cov <- treated_cov[order(names(treated_cov))]

exons_treated <- exons[names(exons) %in% names(treated_TC)]
exons_treated <- exons_treated[order(names(exons_treated))]

f <- function(x,y){Rle(c(x, rep(0,(sum(y@ranges@width) - length(x)))))} ## padding with zeros to tx length
control_TC <- mapply(f, control_TC, exons_control)
control_TCR <- mapply(f, control_TCR, exons_control)
control_cov <- mapply(f, control_cov, exons_control)

treated_TC <- mapply(f, treated_TC, exons_treated)
treated_TCR <- mapply(f, treated_TCR, exons_treated)
treated_cov <- mapply(f, treated_cov, exons_treated)


####### get both control and treated in one granges

g_control <- function(x,y){GRanges(seqnames = rep(as.character(seqnames(y))[1], length(x)),
                                   IRanges(start = 1:length(x), width = 1),
                                   strand = rep(as.character(strand(y))[1], length(x)),
                                   TCR.control = as.vector(x))} # converting Rle into GRanges

g_treated <- function(x,y){GRanges(seqnames = rep(as.character(seqnames(y))[1], length(x)),
                                   IRanges(start = 1:length(x), width = 1),
                                   strand = rep(as.character(strand(y))[1], length(x)),
                                   TCR.treated = as.vector(x))} # converting Rle into GRanges


GR_control <- GRangesList(mapply(g_control, control_TCR, exons_control))
GR_control <- unlist(GR_control)
GR_control$TC.control <- as.numeric(unlist(sapply(control_TC, as.vector)))
GR_control$cov.control <- as.numeric(unlist(sapply(control_cov, as.vector)))
GR_control <- split(GR_control, names(GR_control))
#GR_control <- GR_control[sapply(GR_control, function(x){length(x) > 0})]

GR_treated <- GRangesList(mapply(g_treated, treated_TCR, exons_treated))
GR_treated <- unlist(GR_treated)
GR_treated$TC.treated <- as.numeric(unlist(sapply(treated_TC, as.vector)))
GR_treated$cov.treated <- as.numeric(unlist(sapply(treated_cov, as.vector)))
GR_treated <- split(GR_treated, names(GR_treated))
#GR_treated <- GR_treated[sapply(GR_treated, function(x){length(x) > 0})]


## put the control and treated together
GR_control <- GR_control[names(GR_control) %in% names(GR_treated)]
GR_treated <- GR_treated[names(GR_treated) %in% names(GR_control)]

GR <- unlist(GR_control)
GRt <- unlist(GR_treated)
GR$TCR.treated <- GRt$TCR.treated
GR$TC.treated <- GRt$TC.treated
GR$cov.treated <- GRt$cov.treated


## normalize (dTCR)
dtcr <- (GR$TCR.treated - GR$TCR.control) / (1 - GR$TCR.control)
dtcr[dtcr < 0] <- 0
dtcr[is.na(dtcr)] <- 0
dtcr[dtcr == Inf] <- 0
GR$dTCR <- dtcr


save(GR, file = save_file)


