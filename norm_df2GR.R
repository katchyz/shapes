norm_df2GR <- function(norm_df){
  GRanges(seqnames = norm_df$RNAid, IRanges(start = norm_df$Pos, width = 1), strand = norm_df$strand, norm_df[names(norm_df) != "RNAid" & names(norm_df) != "Pos" & names(norm_df) != "strand"])
}
