#source("https://bioconductor.org/biocLite.R")
library(GenomicAlignments)
library(rtracklayer)

#bam_files = c("2-4cell_NAI", "2-4cell_ctrl", "256cell_NAI", "256cell_ctrl")
sample = "2-4cell_50_mM_NAI"
bw_out = "2-4cell_NAI"
lanes = c("lane1", "lane2", "lane3", "lane4")
lanes = c("lane3", "lane4")
path_io = "/Volumes/USELESS/OUT/SHAPES_out"

for (lane in lanes) {
  bam <- file.path(path_io, sample, lane, "accepted_hits.bam")
  bw_fwd_file <- file.path(path_io, "tracks", paste0(bw_out, "_", lane, "_fwd.bw"))
  bw_rev_file <- file.path(path_io, "tracks", paste0(bw_out, "_", lane, "_rev.bw"))
  #
  RiboReadsRangesResized <- resize(granges(readGAlignmentPairs(bam), use.mcols=T), 1)
  #
  red <- reduce(RiboReadsRangesResized)
  uni <- unique(RiboReadsRangesResized)
  count <- countOverlaps(red, RiboReadsRangesResized)
  count <- countOverlaps(uni, RiboReadsRangesResized)
  red$score <- count
  #
  plus <- red[as.vector(strand(red) == "+")]
  minus <- red[as.vector(strand(red) == "-")]
  #
  export.bw(plus, bw_fwd_file)
  export.bw(minus, bw_rev_file)
  print(lane)
}



bam_file <- "/Volumes/USELESS/OUT/SHAPES_out/256cell_DMSO_control/lane1/accepted_hits.bam"
bw_fwd_file <- "/Users/kasia/Desktop/sam_fwd.bw"
bw_rev_file <- "/Users/kasia/Desktop/sam_rev.bw"

RiboReads <- readGAlignmentPairs(bam_file)
RiboReadsRanges <- granges(RiboReads, use.mcols=T)
RiboReadsRangesResized <- resize(RiboReadsRanges, 1)

red <- reduce(RiboReadsRangesResized)
count <- countOverlaps(red, RiboReadsRangesResized)
red$score <- count

plus <- red[as.vector(strand(red) == "+")]
minus <- red[as.vector(strand(red) == "-")]

export.bw(plus, bw_fwd_file)
export.bw(minus, bw_rev_file)
