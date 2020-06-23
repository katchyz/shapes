library(RNAprobR)

# reading in datasets
samples_path = "/Home/ii/katchyz/OUT/SHAPES_out"
control = c(file.path(samples_path, "summarize_2-4cell_ctrl", "unique_barcodes.txt"))
treated = c(file.path(samples_path, "summarize_2-4cell_NAI", "unique_barcodes.txt"))

control_counts = readsamples(control, euc="counts")
treated_counts = readsamples(treated, euc="counts")

# compiling positional data
control_comp = comp(control_counts)
treated_comp = comp(treated_counts)

# save...
save(control_comp, file="/Home/ii/katchyz/OUT/SHAPES_out/control_comp_2-4cell.Rsave")
save(treated_comp, file="/Home/ii/katchyz/OUT/SHAPES_out/treated_comp_2-4cell.Rsave")

# normalization
shapeseq_norm <- slograt(control_GR=control_comp, treated_GR=treated_comp, depth_correction = "RNA")
#shapeseq_norm <- dtcr(control_GR=control_comp, treated_GR=treated_comp)

# export
bed_zebrafish <- "/Home/ii/katchyz/DATA/genomes/BED/zebrafish_GRCz10.bed"
bedgraph_out <- "/Home/ii/katchyz/OUT/SHAPES_out/dtcr_2-4cell"

norm2bedgraph(shapeseq_norm, bed_file = bed_zebrafish, genome_build = "danRer10", bedgraph_out_file = bedgraph_out, track_name = "dtcr", track_description = "dtcr normalization")
