library(RNAprobR)

# reading in datasets
#samples_path = "/export/valenfs/data/raw_data/SHAPES/"
samples_path = "/Volumes/USELESS/test"
treated = c(file.path(samples_path, "summarize_cell1_invitro_NAI_N3_nonsel", "unique_barcodes.txt"))

treated_counts = readsamples(treated, euc="counts")

txdb_can <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/danio_rerio/Danio_rerio.GRCz10.84.gtf", format = "gtf")
exons <- exonsBy(txdb_can, by="tx", use.names=TRUE)

treated_tx <- mapToTranscripts(treated_counts, exons)
treated_tx$EUC <- treated_counts[treated_tx$xHits]$EUC

# compiling positional data
treated_comp = comp(treated_tx, cutoff=1, exons)

## extend to the full length of TX

# save...
#save(treated_comp, file="/export/valenfs/data/raw_data/SHAPES/normalized/RAW_oblong_CHX_invivo_NAI.Rsave")
save(treated_comp, file="/Volumes/USELESS/test/normalized/RAW_cell1_invitro_NAI_N3_nonsel.Rsave")



load(file="/Volumes/USELESS/test/normalized/RAW_cell2_4_invivo_NAI_N3.Rsave")

# normalization
shapes_norm <- swinsor(Comp_GR=treated_comp, only_top = TRUE)

# save
save(shapes_norm, file="/Volumes/USELESS/test/normalized/cell2_4_invivo_NAI_N3.Rsave")





# export
bed_zebrafish <- "/Home/ii/katchyz/DATA/genomes/BED/zebrafish_GRCz10.bed"
bedgraph_out <- "/export/valenfs/data/raw_data/SHAPES/normalized/bedgraph/swinsor_oblong_CHX_invivo"

norm2bedgraph(shapeseq_norm, bed_file = bed_zebrafish, genome_build = "danRer10", bedgraph_out_file = bedgraph_out, track_name = "swinsor_oblong_CHX_invivo", track_description = "swinsor normalization")


