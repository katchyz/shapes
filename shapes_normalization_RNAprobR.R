library(RNAprobR)

# reading in datasets
samples_path = "/export/valenfs/data/raw_data/SHAPES/"
treated = c(file.path(samples_path, "summarize_cell2_4_invivo_NAI_N3", "unique_barcodes.txt"))

treated_counts = readsamples(treated, euc="counts")

# compiling positional data
treated_comp = comp(treated_counts)
save(treated_comp, file="/export/valenfs/data/raw_data/SHAPES/normalized/RAW_cell2_4_invivo_NAI_N3.Rsave")


# normalization
shapes_norm <- swinsor(treated_comp)

# save
save(shapes_norm, file="/export/valenfs/data/raw_data/SHAPES/normalized/cell2_4_invivo_NAI_N3.Rsave")

# export
bed_zebrafish <- "/Home/ii/katchyz/DATA/genomes/BED/zebrafish_GRCz10.bed"
bedgraph_out <- "/export/valenfs/data/raw_data/SHAPES/normalized/bedgraph/swinsor_cell2_4_invivo_NAI_N3"

norm2bedgraph(shapeseq_norm, bed_file = bed_zebrafish, genome_build = "danRer10", bedgraph_out_file = bedgraph_out, track_name = "swinsor_cell2_4_invivo_NAI_N3", track_description = "swinsor normalization")


###############################################################################
######## save raw !!!!
save(treated_comp, file="/export/valenfs/data/raw_data/SHAPES/normalized/RAW_cell2_4_invivo_NAI_N3.Rsave")
