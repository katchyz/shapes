library(RNAprobR)

# reading in datasets
#samples_path = "/export/valenfs/data/raw_data/SHAPES/"
samples_path = "/Volumes/USELESS/test"
control = c(file.path(samples_path, "summarize_oblong_invivo_DMSO", "unique_barcodes.txt"))
treated = c(file.path(samples_path, "summarize_oblong_invivo_NAI", "unique_barcodes.txt"))

control_counts = readsamples(control, euc="counts")
treated_counts = readsamples(treated, euc="counts")

txdb_can <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/danio_rerio/Danio_rerio.GRCz10.84.gtf", format = "gtf")
exons <- exonsBy(txdb_can, by="tx", use.names=TRUE)

control_tx <- mapToTranscripts(control_counts, exons, ignore.strand = FALSE)
control_tx$EUC <- control_counts[control_tx$xHits]$EUC
#control_tx <- split(control_tx, seqnames(control_tx))

treated_tx <- mapToTranscripts(treated_counts, exons)
treated_tx$EUC <- treated_counts[treated_tx$xHits]$EUC

# compiling positional data
control_comp = comp(control_tx, cutoff=1, exons)
treated_comp = comp(treated_tx, cutoff=1, exons)

# save...
#save(control_comp, file="/export/valenfs/data/raw_data/SHAPES/normalized/RAW_oblong_CHX_invivo_DMSO.Rsave")
save(control_comp, file="/Volumes/USELESS/test/normalized/RAW_oblong_invivo_DMSO.Rsave")
#save(treated_comp, file="/export/valenfs/data/raw_data/SHAPES/normalized/RAW_oblong_CHX_invivo_NAI.Rsave")
save(treated_comp, file="/Volumes/USELESS/test/normalized/RAW_oblong_invivo_NAI.Rsave")




# normalization
shapeseq_norm <- dtcr(control_GR=control_comp, treated_GR=treated_comp)

# save
save(shapeseq_norm, file="/Volumes/USELESS/test/normalized/oblong_CHX_invivo.Rsave")





# export
bed_zebrafish <- "/Home/ii/katchyz/DATA/genomes/BED/zebrafish_GRCz10.bed"
bedgraph_out <- "/export/valenfs/data/raw_data/SHAPES/normalized/bedgraph/swinsor_oblong_CHX_invivo"

norm2bedgraph(shapeseq_norm, bed_file = bed_zebrafish, genome_build = "danRer10", bedgraph_out_file = bedgraph_out, track_name = "swinsor_oblong_CHX_invivo", track_description = "swinsor normalization")


