readsamples <- function(samples, euc="counts", m="", k2n_files=""){
  
  raw_data <- lapply(samples, read.table)
  
  ncols <- sapply(raw_data, ncol)
  if (is.element("TRUE", ncols != 4)){
    stop("All input files should contain exactly 4 columns")
  }
  
  if(euc=="HRF-Seq")
    k2n_values <- lapply(k2n_files, scan, quiet=TRUE)
  
  ###Function body:
  
  #Run proper function depending on euc setting:
  processed_data <- switch(which(euc==c("counts","Fu","HRF-Seq")),
                           .ubar(raw_data), .Fu(raw_data, m),
                           .HRF_EUC(raw_data, k2n_values))
  
  colnames(processed_data) <- c("RNAid", "Start", "End",
                                "Count")
  
  #Lines without end position info - make it equal to start position:
  no_end_info <- is.na(processed_data$End)
  processed_data$End[no_end_info] <- processed_data$Start[no_end_info]
  
  #Modify into GRanges:
  processed_data <- GRanges(seqnames=processed_data$RNAid,
                            IRanges(start=processed_data$Start,
                                    end=processed_data$End),
                            strand="*", EUC=processed_data$Count)
  
  if(is.element(Inf, processed_data$EUC)) {
    Message <- "Barcodes oversaturated. Inf returned.
    Running correct_oversaturation() strongly recommended."
    message(strwrap(Message))
  }
  
  sort(processed_data)
}

###Auxiliary functions

##EUC functions:
#if euc=="counts" - function merging and returning merged data frames,
#if no EUC calculation requested:
.ubar <- function(rdf_list){
  rdf <- do.call("rbind", rdf_list)
  message("Reporting unique barcodes count, no EUC calculation")
  
  rdf
}

#if euc=="Fu" - Function calculating EUC based on number of observed barcodes
#following the formula: n=log((m-k)/m)/log((m-1)/m) derived from
#k=m*(1-(1-1/m)**n), where k is number of observed barcodes, m - number of all
#possible barcodes, n - estimated unique count (number of underlying, target
#molecules).
#Formula from Fu GK et al. PNAS 2011 (binomial distribution calculation).
#Results rounded to nearest integer.
.Fu <- function(rdf_list, m){
  rdf <- do.call("rbind", rdf_list)
  m <- as.integer(m)
  
  #Stop if any record has more observed barcodes than possible (m):
  if(max(rdf[,4]) > m) {
    Message <- "provided 'm' is smaller than the highest observed
    unique barcode count. Revise 'm'"
    stop(strwrap(Message))
  }
  
  rdf[,4] <- round(log((m - rdf[,4])/m)/log((m - 1)/m))
  message("Reporting estimated unique counts according to Fu et al.")
  
  rdf
}

#if euc=="HRF-Seq" - Function calculating EUC based on number of observed
#barcodes following method described in Kielpinski and Vinther, NAR 2014
#(similar to Fu et al. but allows for different barcodes to have different
#attachment probability)
.HRF_EUC <- function(rdf_list, k2n_values){
  for(input_count in 1:length(rdf_list)){
    rdf_list[[input_count]][,4] <-
      k2n_values[[((input_count-1)%%length(k2n_values)+1)]][rdf_list[[
        input_count]][,4]]
  }
  rdf <- do.call("rbind", rdf_list)
  message("Reporting estimated unique counts according to HRF-Seq method")
  
  rdf
}