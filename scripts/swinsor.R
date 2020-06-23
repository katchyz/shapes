swinsor <- function(Comp_GR, winsor_level=0.9, only_top= TRUE, nt_offset=1){
  
  
  ###Function body:
  grl <- split(Comp_GR, seqnames(Comp_GR))
  
  #Do processing (smooth, offset, p-values)
  #normalized <- GRangesList(lapply(grl, FUN=.winsor_oneRNA, winsor_level, only_top, nt_offset))
  
  normalized <- unlist(suppressWarnings(endoapply(grl, FUN=.winsor_oneRNA, winsor_level, only_top, nt_offset)))
  
  #Comp_GR <- unlist(suppressWarnings(endoapply(euc_by_RNA, FUN=.process_oneRNA_euc, TC_all)))
  
  #ugrl <- unlist(normalized)
  #ugrl$winsor[is.na(ugrl$winsor)] <- 0
  normalized$winsor[is.na(normalized$winsor)] <- 0
  normalized <- split(normalized, seqnames(normalized))
  
  normalized
}

###Auxiliary functions:

#Process single RNA:
.winsor_oneRNA <- function(oneRNA_comp_df, winsor_level, only_top, nt_offset){
  w <- winsor(oneRNA_comp_df@elementMetadata$TC,
                                    winsor_level=winsor_level,
                                    only_top=only_top)
  oneRNA_comp_df@elementMetadata$winsor <- w[(1+nt_offset):(length(oneRNA_comp_df)+nt_offset)]
  #oneRNA_comp_df[1:(nrow(oneRNA_comp_df)-nt_offset),]
  oneRNA_comp_df$TC <- oneRNA_comp_df$TC[(1+nt_offset):(length(oneRNA_comp_df)+nt_offset)]
  oneRNA_comp_df
}

