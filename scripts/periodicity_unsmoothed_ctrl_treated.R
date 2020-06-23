load("/Volumes/USELESS/META/SHAPES/normalized/shapeseq_norm_dtcr_2_4.Rsave")
unsmoothed_2_4 <- shapeseq_norm
load("/Volumes/USELESS/META/SHAPES/normalized/shapeseq_norm_dtcr_256.Rsave")
unsmoothed_256 <- shapeseq_norm
rm(shapeseq_norm)

#######
# exon_minus <- exons[strand(exons) == "-"]
# strand(shapeseq_norm[seqnames(shapeseq_norm) %in% names(exon_minus)]) <- rep("-", length(shapeseq_norm[seqnames(shapeseq_norm) %in% names(exon_minus)]))
# gc_unsmoothed_2_4 <- mapFromTranscripts(shapeseq_norm, exons, ignore.strand = FALSE)
# save(gc_unsmoothed_2_4, file="/Home/ii/katchyz/META/SHAPES/gc_unsmoothed_2_4.Rsave")
# 
# exon_minus <- exons[strand(exons) == "-"]
# strand(shapeseq_norm[seqnames(shapeseq_norm) %in% names(exon_minus)]) <- rep("-", length(shapeseq_norm[seqnames(shapeseq_norm) %in% names(exon_minus)]))
# gc_unsmoothed_256 <- mapFromTranscripts(shapeseq_norm, exons, ignore.strand = FALSE)
# save(gc_unsmoothed_256, file="/Home/ii/katchyz/META/SHAPES/gc_unsmoothed_256.Rsave")
#########
# exon_minus <- names(unlist(exons)[strand(unlist(exons)) == "-"])
# exon_plus <- names(unlist(exons)[strand(unlist(exons)) == "+"])

# strand(shapeseq_norm[seqnames(shapeseq_norm) %in% exon_minus]) <- rep("-", length(shapeseq_norm[seqnames(shapeseq_norm) %in% exon_minus]))
# strand(shapeseq_norm[seqnames(shapeseq_norm) %in% exon_plus]) <- rep("+", length(shapeseq_norm[seqnames(shapeseq_norm) %in% exon_plus]))

shapeseq_norm <- shapeseq_norm[seqnames(shapeseq_norm) %in% names(exons)]
gc_unsmoothed_2_4 <- mapFromTranscripts(shapeseq_norm, exons, ignore.strand = FALSE)
gc_unsmoothed_2_4$dtcr <- elementMetadata(shapeseq_norm)$dtcr
gc_unsmoothed_2_4$trnames <- seqnames(shapeseq_norm)
save(gc_unsmoothed_2_4, file="/Home/ii/katchyz/META/SHAPES/gc_unsmoothed_2_4.Rsave")


shapeseq_norm <- shapeseq_norm[seqnames(shapeseq_norm) %in% names(exons)]
gc_unsmoothed_256 <- mapFromTranscripts(shapeseq_norm, exons, ignore.strand = FALSE)
gc_unsmoothed_256$dtcr <- elementMetadata(shapeseq_norm)$dtcr
gc_unsmoothed_256$trnames <- seqnames(shapeseq_norm)
save(gc_unsmoothed_256, file="/Home/ii/katchyz/META/SHAPES/gc_unsmoothed_256.Rsave")
########
load("/Volumes/USELESS/META/SHAPES/gc_unsmoothed_2_4.Rsave")

###### periodicity ######
cds_2_4 <- subsetByOverlaps(gc_unsmoothed_2_4, cds)
cds_2_4_list <- split(cds_2_4, cds_2_4$trnames)
cds_2_4_list <- cds_2_4_list[lengths(cds_2_4_list) > 150]
### get 5'end of CDSs

meta <- lapply(cds_2_4_list, function(x){x[1:150]$dtcr})
meta <- rapply(meta, f=function(x) ifelse(is.na(x),0,x), how="replace")
## REDUCE!!
meta <- Reduce("+", meta)
amplitudes <- abs(fft(meta))
amp <- amplitudes[2:(length(amplitudes)/2+1)]
periods <- 1/periodogram(meta, method = "clone")$freq
df <- data.frame(periods, amp)
p3 <- ggplot(df, aes(x=periods, y=amp)) + geom_line() + ggtitle("FFT, 2-4cell") + xlim(0, 10)
ggsave(file = "/Volumes/USELESS/META/SHAPES/periodicity/p_2_4_unsmoothed.png", plot = p3)

##
load("/Volumes/USELESS/META/SHAPES/gc_unsmoothed_256.Rsave")

cds_256 <- subsetByOverlaps(gc_unsmoothed_256, cds)
cds_256_list <- split(cds_256, cds_256$trnames)
cds_256_list <- cds_256_list[lengths(cds_256_list) > 150]
### get 5'end of CDSs
meta <- lapply(cds_256_list, function(x){x[1:150]$dtcr})
meta <- rapply(meta, f=function(x) ifelse(is.na(x),0,x), how="replace")
## REDUCE!!
meta <- Reduce("+", meta)

amplitudes <- abs(fft(meta))
amp <- amplitudes[2:(length(amplitudes)/2+1)]
periods <- 1/periodogram(meta, method = "clone")$freq
df <- data.frame(periods, amp)
p4 <- ggplot(df, aes(x=periods, y=amp)) + geom_line() + ggtitle("FFT, 256cell") + xlim(0, 10)
ggsave(file = "/Volumes/USELESS/META/SHAPES/periodicity/p_256_unsmoothed.png", plot = p4)

###################
###################
exon_minus <- names(unlist(exons)[strand(unlist(exons)) == "-"])
exon_plus <- names(unlist(exons)[strand(unlist(exons)) == "+"])

# control
strand(control_comp[seqnames(control_comp) %in% exon_minus]) <- rep("-", length(control_comp[seqnames(control_comp) %in% exon_minus]))
strand(control_comp[seqnames(control_comp) %in% exon_plus]) <- rep("+", length(control_comp[seqnames(control_comp) %in% exon_plus]))
control_comp <- control_comp[seqnames(control_comp) %in% names(exons)]
gc_control_comp_256 <- mapFromTranscripts(control_comp, exons, ignore.strand = FALSE)
gc_control_comp_256$Cover <- elementMetadata(control_comp)$Cover

gc_control_comp_2_4$Cover <- NULL
gc_control_comp_2_4$TC <- elementMetadata(control_comp)$TC
#gc_control_comp_256$trnames <- seqnames(control_comp)
save(gc_control_comp_2_4, file="/Home/ii/katchyz/META/SHAPES/gc_control_comp_2_4.Rsave")
# treated
#strand(treated_comp[seqnames(treated_comp) %in% exon_minus]) <- rep("-", length(treated_comp[seqnames(treated_comp) %in% exon_minus]))
#strand(treated_comp[seqnames(treated_comp) %in% exon_plus]) <- rep("+", length(treated_comp[seqnames(treated_comp) %in% exon_plus]))
#treated_comp <- treated_comp[seqnames(treated_comp) %in% names(exons)]
#gc_treated_comp_256 <- mapFromTranscripts(treated_comp, exons, ignore.strand = FALSE)
#gc_treated_comp_256$Cover <- elementMetadata(treated_comp)$Cover
gc_treated_comp_2_4$Cover <- NULL
gc_treated_comp_2_4$TC <- elementMetadata(treated_comp)$TC
#gc_treated_comp_256$trnames <- seqnames(treated_comp)
save(gc_treated_comp_2_4, file="/Home/ii/katchyz/META/SHAPES/gc_treated_comp_2_4.Rsave")

##############
## periodicity on treated and control
path <- "/Volumes/USELESS/META/SHAPES"

samples <- c("gc_control_comp_2_4.Rsave", "gc_control_comp_256.Rsave", "gc_treated_comp_2_4.Rsave", "gc_treated_comp_256.Rsave")

for (s in samples) {
  load(file.path(path, s))
}

# control_comp_2_4
cds_2_4 <- subsetByOverlaps(gc_control_comp_2_4, cds)
cds_2_4_list <- split(cds_2_4, cds_2_4$trnames)
cds_2_4_list <- cds_2_4_list[lengths(cds_2_4_list) > 150]
meta <- lapply(cds_2_4_list, function(x){x[1:150]$TC})
meta <- rapply(meta, f=function(x) ifelse(is.na(x),0,x), how="replace")
meta <- Reduce("+", meta)
amplitudes <- abs(fft(meta))
amp <- amplitudes[2:(length(amplitudes)/2+1)]
periods <- 1/periodogram(meta, method = "clone")$freq
df <- data.frame(periods, amp)
p1 <- ggplot(df, aes(x=periods, y=amp)) + geom_line() + ggtitle("FFT, 2-4cell, control") + xlim(0, 10)
ggsave(file = "/Volumes/USELESS/META/SHAPES/periodicity/p_2_4_control.png", plot = p1)

# treated_comp_2_4
cds_2_4 <- subsetByOverlaps(gc_treated_comp_2_4, cds)
cds_2_4_list <- split(cds_2_4, cds_2_4$trnames)
cds_2_4_list <- cds_2_4_list[lengths(cds_2_4_list) > 150]
meta <- lapply(cds_2_4_list, function(x){x[1:150]$TC})
meta <- rapply(meta, f=function(x) ifelse(is.na(x),0,x), how="replace")
meta <- Reduce("+", meta)
amplitudes <- abs(fft(meta))
amp <- amplitudes[2:(length(amplitudes)/2+1)]
periods <- 1/periodogram(meta, method = "clone")$freq
df <- data.frame(periods, amp)
p2 <- ggplot(df, aes(x=periods, y=amp)) + geom_line() + ggtitle("FFT, 2-4cell, treated") + xlim(0, 10)
ggsave(file = "/Volumes/USELESS/META/SHAPES/periodicity/p_2_4_treated.png", plot = p2)

# control_comp_256
cds_256 <- subsetByOverlaps(gc_control_comp_256, cds)
cds_256_list <- split(cds_256, cds_256$trnames)
cds_256_list <- cds_256_list[lengths(cds_256_list) > 150]
meta <- lapply(cds_256_list, function(x){x[1:150]$TC})
meta <- rapply(meta, f=function(x) ifelse(is.na(x),0,x), how="replace")
meta <- Reduce("+", meta)
amplitudes <- abs(fft(meta))
amp <- amplitudes[2:(length(amplitudes)/2+1)]
periods <- 1/periodogram(meta, method = "clone")$freq
df <- data.frame(periods, amp)
p3 <- ggplot(df, aes(x=periods, y=amp)) + geom_line() + ggtitle("FFT, 256cell, control") + xlim(0, 10)
ggsave(file = "/Volumes/USELESS/META/SHAPES/periodicity/p_256_control.png", plot = p3)

# treated_comp_256
cds_256 <- subsetByOverlaps(gc_treated_comp_256, cds)
cds_256_list <- split(cds_256, cds_256$trnames)
cds_256_list <- cds_256_list[lengths(cds_256_list) > 150]
meta <- lapply(cds_256_list, function(x){x[1:150]$TC})
meta <- rapply(meta, f=function(x) ifelse(is.na(x),0,x), how="replace")
meta <- Reduce("+", meta)
amplitudes <- abs(fft(meta))
amp <- amplitudes[2:(length(amplitudes)/2+1)]
periods <- 1/periodogram(meta, method = "clone")$freq
df <- data.frame(periods, amp)
p4 <- ggplot(df, aes(x=periods, y=amp)) + geom_line() + ggtitle("FFT, 256cell, treated") + xlim(0, 10)
ggsave(file = "/Volumes/USELESS/META/SHAPES/periodicity/p_256_treated.png", plot = p4)

