GR2norm_df <- function(norm_GR){
  norm_methods <- names(mcols(norm_GR))
  RNAid <- levels(seqnames(norm_GR))
  norm_GR <- norm_GR[seqnames(norm_GR) %in% RNAid, norm_methods]
  data.frame(RNAid = as.character(seqnames(norm_GR)), Pos = as.integer(start(norm_GR)),
             strand = as.character(strand(norm_GR)), mcols(norm_GR))
}

