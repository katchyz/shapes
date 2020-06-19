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