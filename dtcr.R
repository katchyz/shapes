dtcr <- function(control_GR, treated_GR, exons){
  
  ###Main Body:
  control <- GR2norm_df(control_GR)
  treated <- GR2norm_df(treated_GR)
  
  control$nt <- NULL
  control$TC <- NULL
  control$Cover <- NULL
  control$PC <- NULL
  treated$nt <- NULL
  treated$TC <- NULL
  treated$Cover <- NULL
  treated$PC <- NULL
  
  #Merge control and treated
  comp_merg <- merge(control, treated, by=c("RNAid", "Pos"), all=TRUE, suffixes=c(".control",".treated"))
  
  #Repair ordering after merging
  comp_merg <- comp_merg[order(comp_merg$RNAid, comp_merg$Pos),]
  #Changes NA to 0 - treat no events as 0 events
  #comp_merg[is.na(comp_merg)] <- 0
  comp_merg$TCR.treated[is.na(comp_merg$TCR.treated)] <- 0
  comp_merg$TCR.control[is.na(comp_merg$TCR.control)] <- 0
  
  #save(comp_merg, file="/Users/kasia/Desktop/comp_merg.Rsave")
  
  ## fix strand
  comp_merg$strand.control[is.na(comp_merg$strand.control)] <- comp_merg$strand.treated[is.na(comp_merg$strand.control)]
  comp_merg$strand.treated[is.na(comp_merg$strand.treated)] <- comp_merg$strand.control[is.na(comp_merg$strand.treated)]
  
  
  ###Calculate deltaTCR as described in HRF-Seq paper:
  #Calculate probabilistic difference
  dtcr <- (comp_merg$TCR.treated - comp_merg$TCR.control)/(1 - comp_merg$TCR.control)
  #If bring_to_zero=TRUE, all negative deltaTCRs bring to 0, else change -Inf
  #to NA [possible when treated=0 and control=1]
  dtcr[dtcr < 0] <- 0
  
  #Import dtcr to merged data frame
  comp_merg$dtcr <- dtcr
  ###
  
  gr <- norm_df2GR(data.frame(RNAid = comp_merg$RNAid, Pos = comp_merg$Pos, strand = comp_merg$strand.treated, dtcr = comp_merg$dtcr))
  grl <- split(gr, seqnames(gr))
  grl <- grl[grl@unlistData@seqnames@lengths > 2]
  ###
  ugrl <- unlist(grl)
  ugrl$dtcr[is.na(ugrl$dtcr)] <- 0
  grl <- split(ugrl, seqnames(ugrl))
  ######################### offset tx by 1
  #grl <- GRangesList(lapply(grl, function(x){GRanges(seqnames(x), ranges(x), dtcr = c(x$dtcr[2:length(x)], 0))}))
  ## USE with OR by
  
  # # rev - reverse elements (for minus strand...)

  # ## SUBSTITUTE strand!!!
  # for (name in names(strand(grl))) {
  #   runValue(strand(grl)[[name]]) <- runValue(strand(exons)[[name]])
  # }
  # 
  # as.character(runValue(strand(exons[[1]])))
  # e <- exons[names(exons) %in% names(grl)]
  # e <- e[order(names(e))]
  # as.character(runValue(strand(e)))
  
  grl
}


