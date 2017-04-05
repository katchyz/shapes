cc <- split(control_comp, seqnames(control_comp))
tc <- split(treated_comp, seqnames(treated_comp))

## add to control ones that exist in treated (with zeros)
toadd <- unlist(tc[!(names(tc) %in% names(cc))])
toadd$TC <- rep(0, length(toadd$TC))
toadd <- split(toadd, seqnames(toadd))
toadd <- toadd[sapply(toadd, function(x){length(x) > 0})]

cc <- c(cc, toadd)
cc <- cc[names(cc) %in% names(tc)]
cc <- cc[order(names(cc))]

ucc <- unlist(cc)
ucc$TC.control <- ucc$TC
ucc$TC.treated <- unlist(tc)$TC
ucc$TC <- NULL
P <- 1
ucc$log2ratio <- log2(ucc$TC.treated+P) - log2(ucc$TC.control+P)

shapeseq_norm <- split(ucc, seqnames(ucc))
shapeseq_norm <- shapeseq_norm[sapply(shapeseq_norm, function(x){length(x) > 0})]

###########################
library(GenomicRanges)

.process_oneRNA_off <- function(rna, df){
  rna@elementMetadata <- rbind(rna@elementMetadata[2:nrow(rna@elementMetadata),], df)
  rna
}

load(file="oblong_invivo.Rsave")
df <- data.frame(TC.control = 0, TC.treated = 0, log2ratio = 0)
shapeseq_norm <- endoapply(shapeseq_norm, FUN=.process_oneRNA_off, df)


save(shapeseq_norm, file="oblong_invivo_OFF.Rsave")


#save(shapeseq_norm, file="/Volumes/USELESS/test/normalized/oblong_CHX_invivo.Rsave")
