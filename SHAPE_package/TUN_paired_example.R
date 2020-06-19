### trim untemplated nucleotides ###

### SERVER VERSION ###

### PAIRED END ###


library(GenomicAlignments)
library(Rsamtools)
#library(readr)
#library(data.table)

lib <- "1Kcell_RZ_vitro_NAI"
shape_path <- "/export/valenfs/data/raw_data/Shape-Seq_Aug2017"

barcode_file <- c(file.path(shape_path, paste("temp_", lib, sep=""), "barcodes.txt"))
save_file <- c(file.path(shape_path, "save_rt_granges", paste("PAIRED_", lib, ".Rsave", sep="")))
bam_file <- c(file.path(shape_path, paste("PAIRED_", lib, sep=""), "accepted_hits.bam"))

# read in BAM file (single end? paired end?)
# get 'MD:Z:' tag
param <- ScanBamParam(tag="MD")
bam_pairs <- readGAlignmentPairs(bam_file, use.names=TRUE, param=param)

### bam_pairs
# get all concordant alignments
bam_pairs <- bam_pairs[!is.na(seqnames(bam_pairs))]
# get uniquely mapping reads (discard multimappers)
bam_pairs <- bam_pairs[unique(names(bam_pairs))]

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


# convert into GRanges (merges read pair)
granges_bam_pairs <- GRanges(bam_pairs)

# split into "+" and "-"
gr_bp_plus <- granges_bam_pairs[strand(granges_bam_pairs) == "+"]
gr_bp_minus <- granges_bam_pairs[strand(granges_bam_pairs) == "-"]

plus <- plus[names(plus) %in% names(gr_bp_plus)]
minus <- minus[names(minus) %in% names(gr_bp_minus)]

# sub newstart and newend -/+ offset of 1
off <- 1
gr_bp_plus <- GRanges(seqnames(gr_bp_plus), IRanges(mcols(plus)$newstart - off, width = (end(gr_bp_plus) - start(gr_bp_plus) - 1)), strand = strand(gr_bp_plus))
names(gr_bp_plus) <- names(plus)

gr_bp_minus <- GRanges(seqnames(gr_bp_minus), IRanges(start(gr_bp_minus), width = (mcols(minus)$newend - start(gr_bp_minus) + off + 1)), strand = strand(gr_bp_minus))
names(gr_bp_minus) <- names(minus)

# merge into granges
rt_granges <- c(gr_bp_plus, gr_bp_minus) 

# add barcodes
barcodes <- data.table::fread(file=barcode_file, sep="\t")
data.table::setDF(barcodes)
colnames(barcodes) <- c("read", "barcode")

bar <- barcodes[barcodes$read %in% names(rt_granges),]
rt_granges <- rt_granges[names(rt_granges) %in% bar$read]
bar <- bar[order(match(bar$read, names(rt_granges))),]

rt_granges$barcode <- bar$barcode

############
unique_granges <- function(sites, sum.counts = FALSE, counts.col = NULL){
  require(dplyr)
  # Checks and balance
  if(!class(sites) == "GRanges"){
    stop("Sites object is not a GRanges class.")}
  if(sum.counts & is.null(counts.col)){
    stop("Please specify the names of the column with count information.")}
  if(!is.null(counts.col)){
    if(!counts.col %in% names(mcols(sites))){
      stop("Could not find counts column name in sites object.")}}
  
  # Convert sites to a data.frame and remove duplicates
  if(!length(names(sites)) == length(unique(names(sites)))){
    message("Dropping rownames for data.frame conversion.")
    df <- GenomicRanges::as.data.frame(sites, row.names = NULL)
  }else{
    df <- GenomicRanges::as.data.frame(sites)
  }
  cols <- names(df)
  
  if(sum.counts){
    counts.pos <- match(counts.col, cols)}
  
  # Sum counts if needed
  if(!sum.counts){
    df <- dplyr::distinct(df)
  }else{
    df$counts <- df[,cols[counts.pos]]
    groups <- lapply(cols[-counts.pos], as.symbol)
    df <- dplyr::group_by_(df, .dots = groups) %>%
      dplyr::summarise(counts = sum(counts))
    names(df) <- c(cols[-counts.pos], cols[counts.pos])
  }
  
  # Rebuild GRanges object
  gr <- GRanges(
    seqnames = df$seqnames,
    ranges = IRanges(start = df$start, end = df$end),
    strand = df$strand,
    seqinfo = seqinfo(sites)
  )
  
  mcols(gr) <- df[,6:length(df)]
  gr
}
############

# count unique barcodes
rt_granges_unique <- unique(rt_granges)
rt_granges_unique$unique_barcodes <- countOverlaps(rt_granges_unique, unique_granges(rt_granges), type="equal")

#### save
save(rt_granges_unique, file=save_file)
#### pass to dtcr

