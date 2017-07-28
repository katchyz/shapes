### trim untemplated nucleotides ###

library(GenomicAlignments)
library(Rsamtools)

# read in BAM file (single end? paired end?)

bam_file <- "/Volumes/USELESS/DATA/Shape-Seq/demo_SHAPE_package/256cell_ctrl.bam"
# get 'MD:Z:' tag
param <- ScanBamParam(tag="MD")
bam_pairs <- readGAlignmentPairs(bam_file, use.names=TRUE, param=param)

### bam_pairs
# get all concordant alignments
bam_pairs <- bam_pairs[!is.na(seqnames(bam_pairs))]

################
# get 1st read
#first <- first(bam_pairs)

# look at MD tag
#mcols(first(bam_pairs))$MD

## sam <- c("0A13", "0T20", "0C14", "15G10", "1C0T31", "0C5T12")
## grepl("^0[ATCG]", sam)

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



##############
##############
bp <- c(bam_pairs[11:15], bam_pairs[194:195], bam_pairs[585:587])
# convert into GRanges (merges read pair)
m <- minus[11:15] 
p <- plus[1:5]
gbp <- GRanges(bp)
# sub newstart and newend for plus and minus strand
# offset by 1
off <- 1
mgr <- GRanges(seqnames(m), IRanges(start(m) - off, width = (mcols(m)$newend - start(m) + 1)), strand = strand(m))
names(mgr) <- names(m)
pgr <- GRanges(seqnames(p), IRanges(mcols(p)$newstart + off, width = (end(p) - mcols(p)$newstart + 1)), strand = strand(p))
names(pgr) <- names(p)

gr <- c(mgr, pgr)

# add barcodes
b <- bar[bar$read %in% names(gr),]
b <- b[order(match(b$read, names(gr))),]

gr$barcode <- b$barcode

# count unique barcodes
gru <- unique(gr)
gru$unique_barcodes <- countOverlaps(gru, gr, type="equal")

# get TC position ??? (start for pgr, end for mgr)



###############################################
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

barcode_file <- "/Volumes/USELESS/DATA/Shape-Seq/demo_SHAPE_package/bar.txt"
barcodes <- read.csv(file=barcode_file, sep="\t", header = F)
colnames(barcodes) <- c("read", "barcode")

bar <- barcodes[barcodes$read %in% names(bam_pairs),]
