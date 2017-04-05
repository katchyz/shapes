library(rtracklayer)

gtf_annot <- import("/Volumes/USELESS/DATA/genomes/GTF/danio_rerio/Danio_rerio.GRCz10.84.gtf", format = "gtf")
gtf_annot <- gtf_annot[gtf_annot$type == "transcript"]
transcript_df <- data.frame(n = gtf_annot$transcript_id, tx_biotype = gtf_annot$transcript_biotype)

protein_coding <- as.character(subset(transcript_df, tx_biotype == "protein_coding")$n)
protein_coding <- protein_coding[order(protein_coding)]

pt_shape_cell24 <- shape_cell24[names(shape_cell24) %in% protein_coding]
pt_shape_cell256 <- shape_cell256[names(shape_cell256) %in% protein_coding]
pt_shape_oblong <- shape_oblong[names(shape_oblong) %in% protein_coding]
pt_shape_oblong_CHX <- shape_oblong_CHX[names(shape_oblong_CHX) %in% protein_coding]
pt_shapes_cell24 <- shapes_cell24[names(shapes_cell24) %in% protein_coding]
pt_shapes_cell1_invitro <- shapes_cell1_invitro[names(shapes_cell1_invitro) %in% protein_coding]
pt_shapes_cell1_invitro_nonsel <- shapes_cell1_invitro_nonsel[names(shapes_cell1_invitro_nonsel) %in% protein_coding]

Reduce("+", sapply(pt_shape_cell24, function(x){sum(x$TC.control)}))        # 9037691
Reduce("+", sapply(pt_shape_cell24, function(x){sum(x$TC.treated)}))        # 6288159
Reduce("+", sapply(pt_shape_cell256, function(x){sum(x$TC.control)}))       # 5750132
Reduce("+", sapply(pt_shape_cell256, function(x){sum(x$TC.treated)}))       # 8569845
Reduce("+", sapply(pt_shape_oblong, function(x){sum(x$TC.control)}))        # 2648256
Reduce("+", sapply(pt_shape_oblong, function(x){sum(x$TC.treated)}))        # 1467983
Reduce("+", sapply(pt_shape_oblong_CHX, function(x){sum(x$TC.control)}))    # 2418685
Reduce("+", sapply(pt_shape_oblong_CHX, function(x){sum(x$TC.treated)}))    # 21600594

Reduce("+", sapply(pt_shapes_cell24, function(x){sum(x$TC)}))               # 2566584
Reduce("+", sapply(pt_shapes_cell1_invitro, function(x){sum(x$TC)}))        # 408662
Reduce("+", sapply(pt_shapes_cell1_invitro_nonsel, function(x){sum(x$TC)})) # 553527

