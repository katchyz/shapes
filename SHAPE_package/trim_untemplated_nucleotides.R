### trim untemplated nucleotides ###

library(GenomicAlignments)
library(Rsamtools)
library(purrr)

# read in BAM file (single end? paired end?)

bam_file <- "/Volumes/USELESS/OUT/SHAPES_out/256cell_DMSO_control/lane1/accepted_hits.bam"
#bam_reads <- readGAlignmentPairs(bam_file)
# get 'MD:Z:' tag

#myBf <- BamFile(bam_file, index=paste0(bam_file, ".bai"), asMates=TRUE)
#md <- ScanBamParam(what=scanBamWhat(), tag="MD", flag=scanBamFlag(isUnmappedQuery=FALSE))

###
param <- ScanBamParam(tag="MD")
#bam <- readGAlignments(bam_file, use.names=TRUE, param=param)
bam_pairs <- readGAlignmentPairs(bam_file, use.names=TRUE, param=param)
#bam_list <- readGAlignmentsList(bam_file, use.names=TRUE, param=param)
#bam_gap <- readGappedReads(bam_file, use.names=TRUE, param=param)
#mcols(bam) # access MD tag

# get all mated reads
##bam_list[mcols(bam_list)$mate_status == "mated"]

### bam_pairs
# get all concordant alignments
bam_pairs <- bam_pairs[!is.na(seqnames(bam_pairs))]

################
# get 1st read
first <- first(bam_pairs)

# look at MD tag
mcols(first(bam_pairs))$MD

## sam <- c("0A13", "0T20", "0C14", "15G10", "1C0T31", "0C5T12")
## grepl("^0[ATCG]", sam)

trim1_fwd <- c("^0[ATCG][2-9][0-9]*")
trim2_fwd <- c("^1[ATCG][1-9][0-9]*", "^0[ATCG]0[ATCG][1-9][0-9]*")
trim3_fwd <- c("^2[ATCG]", "^1[ATCG]0[ATCG]", "^0[ATCG]1[ATCG]", "^0[ATCG]0[ATCG]0[ATCG]")

## sam <- c("13A0", "20T0", "14C0", "15G10", "310T1C", "12T5C0")
## grepl("[ATCG]0$", sam)

trim1_rev <- c("[2-9][0-9]*[ATCG]0$")
trim2_rev <- c("[1-9][0-9]*[ATCG]1$", "[1-9][0-9]*[ATCG]0[ATCG]0$")
trim3_rev <- c("[ATCG]2$", "[ATCG]0[ATCG]1$", "[ATCG]1[ATCG]0$", "[ATCG]0[ATCG]0[ATCG]0$")

#filter(d, grepl(paste(pat, collapse="|"), b))
plus <- first[strand(first) == "+"]
minus <- first[strand(first) == "-"]

#filter(data.frame(mcols(plus)), grepl(paste(trim2_fwd, collapse="|")), MD)
#grep(paste(c("[2-9][0-9]*[ATCG]0$", "[1-9][0-9]*[ATCG]1$"),collapse="|"), d$b)

#grep(paste(trim2_fwd, collapse="|"), mcols(plus)$MD) ### indices of reads in 'plus' subject to trimming 2nt

# move the start/end depending on offset (MD) and strand
# plus strand - offset start, minus strand - offset end (look at cigar too, to avoid introns?)
# cigar - if at the beginning/end there is [12]M - offset adding/subtracting number of Ns


## i <- extractAlignmentRangesOnReference("3M100N17M", pos = 555)
## unlist(tile(i, width=1))[5]

## cs <- cumsum(width(i))
## idx <- which(cs > 5)
## start(i[idx]) + (5 - cs[idx-1]) - 1

# MD:Z:0


trim_untemplated_nt <- function(bam_pairs) {
  
  ### trimming patterns
  trim1_fwd <- c("^0[ATCG][2-9][0-9]*")
  trim2_fwd <- c("^1[ATCG][1-9][0-9]*", "^0[ATCG]0[ATCG][1-9][0-9]*")
  trim3_fwd <- c("^2[ATCG]", "^1[ATCG]0[ATCG]", "^0[ATCG]1[ATCG]", "^0[ATCG]0[ATCG]0[ATCG]")
  
  trim1_rev <- c("[2-9][0-9]*[ATCG]0$")
  trim2_rev <- c("[1-9][0-9]*[ATCG]1$", "[1-9][0-9]*[ATCG]0[ATCG]0$")
  trim3_rev <- c("[ATCG]2$", "[ATCG]0[ATCG]1$", "[ATCG]1[ATCG]0$", "[ATCG]0[ATCG]0[ATCG]0$")
  
  ### trim first read, split by plus and minus strand
  first <- first(bam_pairs)
  plus <- first[strand(first) == "+"]
  minus <- first[strand(first) == "-"]
  
  ### get indices of reads to be trimmed
  t1_fwd <- grep(paste(trim1_fwd, collapse="|"), mcols(plus)$MD)
  t2_fwd <- grep(paste(trim2_fwd, collapse="|"), mcols(plus)$MD)
  t3_fwd <- grep(paste(trim3_fwd, collapse="|"), mcols(plus)$MD)
  
  t1_rev <- grep(paste(trim1_rev, collapse="|"), mcols(minus)$MD)
  t2_rev <- grep(paste(trim2_rev, collapse="|"), mcols(minus)$MD)
  t3_rev <- grep(paste(trim3_rev, collapse="|"), mcols(minus)$MD)
  
  ### add mcol with offsets?
  mcols(plus)$offset <- rep(0, length(plus))
  mcols(minus)$offset <- rep(0, length(minus))
  
  mcols(plus)$offset[t1_fwd] <- 1
  mcols(plus)$offset[t2_fwd] <- 2
  mcols(plus)$offset[t3_fwd] <- 3
  
  mcols(minus)$offset[t1_rev] <- 1
  mcols(minus)$offset[t2_rev] <- 2
  mcols(minus)$offset[t3_rev] <- 3
  
  mcols(plus)$newstart <- rep(0, length(plus))
  mcols(minus)$newend <- rep(0, length(minus))
  
  ### offset - plus ###
  i <- extractAlignmentRangesOnReference(cigar(plus), start(plus))
  cs <- cumsum(width(i))
  #idx <- map(as.list(which(cs > mcols(plus)$offset)), 1)
  idx <- which(cs > mcols(plus)$offset)
  idx <- lapply(idx, function(x){x[1]})
  
  mcols(plus)$newstart[which(idx == 1)] <- sapply(start(i[which(idx == 1)]), function(x){x[1]}) + 
    mcols(plus)$offset[which(idx == 1)]
  
  mcols(plus)$newstart[which(idx > 1)] <- 
    mapply(function(x,y){x[y]}, as.list(start(i[which(idx > 1)])), idx[which(idx > 1)]) + 
    (mcols(plus)$offset[which(idx > 1)] - mapply(function(x,y){x[y-1]}, cs[which(idx > 1)], idx[which(idx > 1)]))
  
  ### offset - minus ###
  i <- extractAlignmentRangesOnReference(cigar(minus), start(minus))
  i[width(i@partitioning) > 1] <- IRangesList(lapply(i[width(i@partitioning) > 1], function(x){reverse(x)}))
  #ri <- lapply(i, function(x){reverse(x)})
  #ri <- lapply(extractAlignmentRangesOnReference(cigar(minus), start(minus)), function(x){rev(x)})
  cs <- cumsum(width(i))
  #idx <- map(as.list(which(cs > mcols(minus)$offset)), 1)
  idx <- which(cs > mcols(minus)$offset)
  idx <- lapply(idx, function(x){x[1]})
  
  mcols(minus)$newend[which(idx == 1)] <- sapply(end(i[which(idx == 1)]), function(x){x[1]}) - 
    mcols(minus)$offset[which(idx == 1)]
  
  mcols(minus)$newend[which(idx > 1)] <- 
    mapply(function(x,y){x[y]}, as.list(end(i[which(idx > 1)])), idx[which(idx > 1)]) - 
    (mcols(minus)$offset[which(idx > 1)] - mapply(function(x,y){x[y-1]}, cs[which(idx > 1)], idx[which(idx > 1)]))
  
  
  ###############################################
  ## sub newstart for start, newend for end
  plus@start <- as.integer(mcols(plus)$newstart)
  minus@end <- as.integer(mcols(minus)$newend)
  
  ## remove extra columns
  mcols(plus)$newstart <- NULL
  mcols(plus)$offset <- NULL
  mcols(minus)$newend <- NULL
  mcols(minus)$offset <- NULL
  
  ## merge plus and minus, sort
  newfirst <- c(plus, minus)
  newfirst <- newfirst[order(match(names(newfirst), names(first)))]
  
  ## sub first into bam_pairs
  bam_pairs@first <- newfirst
  
}


offset_plus <- function(read, offset=1) {
  i <- extractAlignmentRangesOnReference(cigar(read), start(read))
  cs <- cumsum(width(i))
  idx <- as.numeric(which(cs > offset)[1])
  if (idx > 1) {
    start <- as.numeric(start(i[idx]) + (offset - cs[idx-1]))
  } else {
    start <- as.numeric(start(i[1]) + offset)
  }
  return(start)
}

offset_minus <- function(read, offset=1) {
  ri <- rev(extractAlignmentRangesOnReference(cigar(read), start(read)))
  cs <- cumsum(width(ri))
  idx <- as.numeric(which(cs > offset)[1])
  if (idx > 1) {
    end <- as.numeric(end(ri[idx]) - (offset - cs[idx-1]))
  } else {
    end <- as.numeric(end(ri[1]) - offset)
  }
  return(end)
}


################
################

barcode_file <- "/Volumes/USELESS/DATA/Shape-Seq/OUR_zebrafish/preproc_256cell_ctrl/bar.txt"
barcodes <- read.csv(file=barcode_file, sep="\t", header = F)
colnames(barcodes) <- c("read", "barcode")

bar <- barcodes[barcodes$read %in% names(bam_pairs),]
