comp <- function(euc_GR, cutoff=1, exons){
  
  ###Function body:
  #Remove inserts shorter than cutoff (keep removed in removed_GR:
  good_length <- (width(euc_GR) >= cutoff)
  removed_GR <- euc_GR[!good_length]
  euc_GR_good <- euc_GR[good_length]
  
  
  #Calculate TC using coverage function. If element length is set to 1 at the
  #stop site then it corresponds to termination count.
  euc_forTC <- euc_GR_good
  #end(euc_forTC) <- start(euc_forTC)
  TC_all <- coverage(euc_forTC, weight=euc_forTC$EUC)
  n <- names(TC_all)
  TC_all <- lapply(names(TC_all), function(x){Rle(c(as.numeric(TC_all[[x]]), rep(0,(sum(exons[[x]]@ranges@width)-length(as.numeric(TC_all[[x]]))))))})
  names(TC_all) <- n
  
  
  ###Run the single RNA processing for all the RNAs:
  x <- BiocGenerics::as.vector(seqnames(euc_GR_good))
  euc_by_RNA <- split(euc_GR_good, f=x, drop=TRUE)
  Comp_GR <- unlist(suppressWarnings(endoapply(euc_by_RNA,
                                               FUN=.process_oneRNA_euc, TC_all)))
  
  
  #Print info on fraction of removed EUC's:
  percent_removed <- sum(removed_GR$EUC)/(sum(removed_GR$EUC) +
                                            sum(euc_GR_good$EUC))*100
  message(paste(round(percent_removed,2), "% of EUCs removed due to cutoff"))
  
  Comp_GR
}

###Auxiliary functions

.process_oneRNA_euc <- function(oneRNA_euc, TC_all){
  
  #Name of analysed RNA
  RNAid <- as.character(seqnames(oneRNA_euc[1]))
  
  #Which coverage vector in the coverage list corresponds to our gene. It is
  #calculated once and assumed to be the same for all three coverage vectors:
  #TC (termination counts), PC (priming counts) and coverage.
  RNA_order <- which(names(TC_all)==RNAid)
  
  if(length(RNA_order) > 1)
    stop(paste("ERROR: More than one gene with the same ID provided. Check",
               RNAid))
  
  #Extracting coverage vectors from coverage list for TC, PC and coverage:
  gene_TC <- as.numeric(TC_all[[RNA_order]]) #coverage to vector
  
  #And constructing GRanges object for single RNA:
  GRanges(seqnames=RNAid,
          IRanges(start=1:length(gene_TC), width=1),
          strand=as.character(runValue(strand(oneRNA_euc))),
          TC=gene_TC)
}