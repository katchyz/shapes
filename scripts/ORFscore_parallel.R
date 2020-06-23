library(GenomicRanges)
library(GenomicFeatures)
library(parallel)
txdb_can <- makeTxDbFromGFF("/Home/ii/katchyz/DATA/genomes/GTF/zebrafish_GRCz10.gtf", format = "gtf")
cds <- cdsBy(txdb_can, by="tx", use.names=TRUE)

#path <- "/Volumes/USELESS/META/SHAPES"
path <- "/Home/ii/katchyz/META/SHAPES"

files <- c("gc_control_comp_256.Rsave", "gc_control_comp_2_4.Rsave", "gc_treated_comp_256.Rsave", "gc_treated_comp_2_4.Rsave", "gc_unsmoothed_256.Rsave", "gc_unsmoothed_2_4.Rsave")

for (file in files) {
  load(file.path(path, file))
}


## Number of workers (R processes) to use:
numWorkers <- 8


delete.NULLs  <-  function(x.list){   # delele null/empty entries in a list
  x.list[unlist(mclapply(x.list, length, mc.cores = numWorkers) != 0)]
}

#control
ctrl_256 <- subsetByOverlaps(gc_control_comp_256, cds)
ctrl_256$TC[is.na(ctrl_256$TC)] <- 0
ctrl_256 <- split(ctrl_256, ctrl_256$trnames)
ctrl_256 <- delete.NULLs(ctrl_256)

ctrl_2_4 <- subsetByOverlaps(gc_control_comp_2_4, cds)
ctrl_2_4$TC[is.na(ctrl_2_4$TC)] <- 0
ctrl_2_4 <- split(ctrl_2_4, ctrl_2_4$trnames)
ctrl_2_4 <- delete.NULLs(ctrl_2_4)

# treated
trea_256 <- subsetByOverlaps(gc_treated_comp_256, cds)
trea_256$TC[is.na(trea_256$TC)] <- 0
trea_256 <- split(trea_256, trea_256$trnames)
trea_256 <- delete.NULLs(trea_256)

trea_2_4 <- subsetByOverlaps(gc_treated_comp_2_4, cds)
trea_2_4$TC[is.na(trea_2_4$TC)] <- 0
trea_2_4 <- split(trea_2_4, trea_2_4$trnames)
trea_2_4 <- delete.NULLs(trea_2_4)

#unsmoothed
unsm_256 <- subsetByOverlaps(gc_unsmoothed_256, cds)
unsm_256$dtcr[is.na(unsm_256$dtcr)] <- 0
unsm_256 <- split(unsm_256, unsm_256$trnames)
unsm_256 <- delete.NULLs(unsm_256)

unsm_2_4 <- subsetByOverlaps(gc_unsmoothed_2_4, cds)
unsm_2_4$dtcr[is.na(unsm_2_4$dtcr)] <- 0
unsm_2_4 <- split(unsm_2_4, unsm_2_4$trnames)
unsm_2_4 <- delete.NULLs(unsm_2_4)

#### ORFs
# control
ctrl_256_orf1 <- sum(unlist(mclapply(ctrl_256, function(x){sum(x[seq(1, length(x), 3)]$TC)}, mc.cores = numWorkers)))
ctrl_256_orf2 <- sum(unlist(mclapply(ctrl_256, function(x){sum(x[seq(2, length(x), 3)]$TC)}, mc.cores = numWorkers)))
ctrl_256_orf3 <- sum(unlist(mclapply(ctrl_256, function(x){sum(x[seq(3, length(x), 3)]$TC)}, mc.cores = numWorkers)))

ctrl_2_4_orf1 <- sum(unlist(mclapply(ctrl_2_4, function(x){sum(x[seq(1, length(x), 3)]$TC)}, mc.cores = numWorkers)))
ctrl_2_4_orf2 <- sum(unlist(mclapply(ctrl_2_4, function(x){sum(x[seq(2, length(x), 3)]$TC)}, mc.cores = numWorkers)))
ctrl_2_4_orf3 <- sum(unlist(mclapply(ctrl_2_4, function(x){sum(x[seq(3, length(x), 3)]$TC)}, mc.cores = numWorkers)))

# treated
trea_256_orf1 <- sum(unlist(mclapply(trea_256, function(x){sum(x[seq(1, length(x), 3)]$TC)}, mc.cores = numWorkers)))
trea_256_orf2 <- sum(unlist(mclapply(trea_256, function(x){sum(x[seq(2, length(x), 3)]$TC)}, mc.cores = numWorkers)))
trea_256_orf3 <- sum(unlist(mclapply(trea_256, function(x){sum(x[seq(3, length(x), 3)]$TC)}, mc.cores = numWorkers)))

trea_2_4_orf1 <- sum(unlist(mclapply(trea_2_4, function(x){sum(x[seq(1, length(x), 3)]$TC)}, mc.cores = numWorkers)))
trea_2_4_orf2 <- sum(unlist(mclapply(trea_2_4, function(x){sum(x[seq(2, length(x), 3)]$TC)}, mc.cores = numWorkers)))
trea_2_4_orf3 <- sum(unlist(mclapply(trea_2_4, function(x){sum(x[seq(3, length(x), 3)]$TC)}, mc.cores = numWorkers)))

# unsmoothed
#unsm_256_orf1 <- sum(unlist(mclapply(unsm_256, function(x){sum(x[seq(1, length(x), 3)]$dtcr)}, mc.cores = numWorkers)))
unsm_256_orf2 <- sum(unlist(mclapply(unsm_256, function(x){sum(as.numeric(x[seq(2, length(x), 3)]$dtcr))}, mc.cores = numWorkers)))
unsm_256_orf3 <- sum(unlist(mclapply(unsm_256, function(x){sum(as.numeric(x[seq(3, length(x), 3)]$dtcr))}, mc.cores = numWorkers)))

#unsm_2_4_orf1 <- sum(unlist(mclapply(unsm_2_4, function(x){sum(x[seq(1, length(x), 3)]$dtcr)}, mc.cores = numWorkers)))
unsm_2_4_orf2 <- sum(unlist(mclapply(unsm_2_4, function(x){sum(as.numeric(x[seq(2, length(x), 3)]$dtcr))}, mc.cores = numWorkers)))
unsm_2_4_orf3 <- sum(unlist(mclapply(unsm_2_4, function(x){sum(as.numeric(x[seq(3, length(x), 3)]$dtcr))}, mc.cores = numWorkers)))


###########
unsm_2_4_orf2 <- sum(sapply(unsm_2_4, function(x){sum(x[seq(2, length(x), 3)]$dtcr)}))



for (name in names(unsm_2_4)) {
  print(name)
  print(sum(unsm_2_4[[name]][seq(2, length(unsm_2_4[[name]]), 3)]$dtcr))
}

