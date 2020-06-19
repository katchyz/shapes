### matrix of everything
library(GenomicFeatures)

# tx / gene id
# tx / utr5 / cds / utr3 length
# RNA
# Ribo
# GC content
# stop codon
# max of RTprimer_aln_width, RTc_aln_widt, RTrc_aln_width
# GO_cc
# polyA tail length
# CAI
# min free E (rnafold)
# Shape

# cloud / not cloud
# m6A ? (position?)
# stall site ? (in cds? in utr3?)

# -------------------------------------------------------
# combine features... (TE, accessibility), bin into high/low

txdb_can <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/danio_rerio/Danio_rerio.GRCz10.84.gtf", format = "gtf")
txLengths <- transcriptLengths(txdb_can, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
rownames(txLengths) <- txLengths$tx_name
txLengths <- txLengths[order(rownames(txLengths)),]

matrix_of_everything <- data.frame(n = txLengths$tx_name, gene_id = txLengths$gene_id, tx_len = txLengths$tx_len,
                                   cds_len = txLengths$cds_len, utr5_len = txLengths$utr5_len, utr3_len = txLengths$utr3_len)


gene_names <- read.csv(file = "/Volumes/USELESS/DATA/genomes/ensembl_gene_name.txt", header = TRUE)
gene_names <- data.frame(n = gene_names$Transcript.stable.ID, gene_name = gene_names$Gene.name)
gene_names <- gene_names[gene_names$n %in% matrix_of_everything$n,]
matrix_of_everything <- merge(matrix_of_everything, gene_names, by="n", all=TRUE)


load(file="/Volumes/USELESS/META/SHAPES_NEW/general/clouds/cloud_df.Rsave")
df <- df[df$n %in% matrix_of_everything$n,]
rownames(df) <- df$n
df <- df[order(rownames(df)),]
df$n <- matrix_of_everything$n

matrix_of_everything$rna_cell1 <- df$rna_cell1
matrix_of_everything$rna_cell1_total <- df$rna_cell1_total
matrix_of_everything$rna_cell24 <- df$rna_cell24
matrix_of_everything$rna_2h_total <- df$rna_2h_total
matrix_of_everything$rna_cell256 <- df$rna_cell256
matrix_of_everything$rna_cell1K <- df$rna_cell1K
matrix_of_everything$rna_3_5h <- df$rna_3_5h
matrix_of_everything$rna_4h_total <- df$rna_4h_total

matrix_of_everything$ribo_cell24 <- df$ribo_cell24
matrix_of_everything$ribo_cell256 <- df$ribo_cell256
matrix_of_everything$ribo_cell1K <- df$ribo_cell1K
# load(file = "/Volumes/USELESS/META/SHAPES_OLD/FPKM_4h_ribo.Rdata")
FPKM_4h <- FPKM_4h[FPKM_4h$n %in% matrix_of_everything$n,]
matrix_of_everything <- merge(matrix_of_everything, FPKM_4h, by="n", all=TRUE)

### save(matrix_of_everything, file = "/Volumes/USELESS/META/SHAPES/matrix_of_everything.Rdata")

matrix_of_everything$gc <- df$gc
matrix_of_everything$stop_codon <- df$stop_codon
# df$RT_primer_aln_width <- pmax(df$RTprimer_aln_width, df$RTc_aln_width, df$RTrc_aln_width)
matrix_of_everything$RT_primer_aln_width <- df$RT_primer_aln_width

go_cc <- read.csv(file = "/Volumes/USELESS/DATA/genomes/GO/go_cc.tsv", sep = "\t", header = TRUE)
go_cc <- arrange(go_cc, ID)
go_cc <- go_cc[!duplicated(go_cc$ID),]
go_cc <- go_cc[go_cc$ID %in% matrix_of_everything$n,]
go_cc <- data.frame(n = go_cc$ID, go_cc = go_cc$GO.Name)
matrix_of_everything <- merge(matrix_of_everything, go_cc, by="n", all=TRUE)
# pal from giant_matrix
matrix_of_everything <- merge(matrix_of_everything, pal, by="n", all=TRUE)
# cai from giant_matrix
matrix_of_everything <- merge(matrix_of_everything, cai, by="n", all=TRUE)
# rnafold from giant_matrix
matrix_of_everything <- merge(matrix_of_everything, rnafold, by="n", all=TRUE)

gtf_annot <- import("/Volumes/USELESS/DATA/genomes/GTF/danio_rerio/Danio_rerio.GRCz10.84.gtf", format = "gtf")
gtf_annot <- gtf_annot[gtf_annot$type == "transcript"]
transcript_df <- data.frame(n = gtf_annot$transcript_id, tx_biotype = gtf_annot$transcript_biotype)
matrix_of_everything <- merge(matrix_of_everything, transcript_df, by="n", all=TRUE)

###### SHAPE
libs_path <- "/Volumes/USELESS/DATA/SHAPES/SE"
libs <- c("norm2_cell1_invitro_NAIN3.Rsave", "norm2_sphere_invitro_A.Rsave", "norm2_sphere_invitro_C.Rsave",
          "norm2_cell24_invivo.Rsave", "norm2_cell256_invivo.Rsave", "norm2_oblong_CHX_invivo.Rsave", "norm2_oblong_invivo.Rsave")

shape <- get(load(file = c(file.path(libs_path, libs[1]))))
TCtr <- sapply(shape, function(x){mean(x$TC.treated)})
sdf <- data.frame(n = names(TCtr), shapes_cell1_invitro_TC = as.numeric(TCtr))

shape_df <- sdf

shape <- get(load(file = c(file.path(libs_path, libs[2]))))
TCtr <- sapply(shape, function(x){mean(x$TC.treated)})
sdf <- data.frame(n = names(TCtr), shapes_sphere_invitro_A_TC = as.numeric(TCtr))

shape_df <- merge(shape_df, sdf, by="n", all=TRUE)

shape <- get(load(file = c(file.path(libs_path, libs[3]))))
TCtr <- sapply(shape, function(x){mean(x$TC.treated)})
sdf <- data.frame(n = names(TCtr), shapes_sphere_invitro_C_TC = as.numeric(TCtr))

shape_df <- merge(shape_df, sdf, by="n", all=TRUE)

shape <- get(load(file = c(file.path(libs_path, libs[4]))))
TCtr <- sapply(shape, function(x){mean(x$TC.treated)})
sdf <- data.frame(n = names(TCtr), shape_cell24_invivo_log2ratio = as.numeric(TCtr))

shape_df <- merge(shape_df, sdf, by="n", all=TRUE)

shape <- get(load(file = c(file.path(libs_path, libs[5]))))
TCtr <- sapply(shape, function(x){mean(x$TC.treated)})
sdf <- data.frame(n = names(TCtr), shape_cell256_invivo_log2ratio = as.numeric(TCtr))

shape_df <- merge(shape_df, sdf, by="n", all=TRUE)

shape <- get(load(file = c(file.path(libs_path, libs[6]))))
TCtr <- sapply(shape, function(x){mean(x$TC.treated)})
sdf <- data.frame(n = names(TCtr), shape_oblong_CHX_invivo_log2ratio = as.numeric(TCtr))

shape_df <- merge(shape_df, sdf, by="n", all=TRUE)

shape <- get(load(file = c(file.path(libs_path, libs[7]))))
TCtr <- sapply(shape, function(x){mean(x$TC.treated)})
sdf <- data.frame(n = names(TCtr), shape_oblong_invivo_log2ratio = as.numeric(TCtr))

shape_df <- merge(shape_df, sdf, by="n", all=TRUE)

matrix_of_everything <- merge(matrix_of_everything, shape_df, by="n", all=TRUE)


#####################
###### SHAPE_internal
libs_path <- "/Volumes/USELESS/DATA/SHAPES/SE"
libs <- c("norm2_cell1_invitro_NAIN3.Rsave", "norm2_sphere_invitro_A.Rsave", "norm2_sphere_invitro_C.Rsave",
          "norm2_cell24_invivo.Rsave", "norm2_cell256_invivo.Rsave", "norm2_oblong_CHX_invivo.Rsave", "norm2_oblong_invivo.Rsave")

shape <- get(load(file = c(file.path(libs_path, libs[1]))))
TCtr <- sapply(shape, function(x){mean(x[2:length(x)]$TC.treated)})
sdf <- data.frame(n = names(TCtr), shapes_cell1_invitro_TC_internal = as.numeric(TCtr))

shape_df <- sdf

shape <- get(load(file = c(file.path(libs_path, libs[2]))))
TCtr <- sapply(shape, function(x){mean(x[2:length(x)]$TC.treated)})
sdf <- data.frame(n = names(TCtr), shapes_sphere_invitro_A_TC_internal = as.numeric(TCtr))

shape_df <- merge(shape_df, sdf, by="n", all=TRUE)

shape <- get(load(file = c(file.path(libs_path, libs[3]))))
TCtr <- sapply(shape, function(x){mean(x[2:length(x)]$TC.treated)})
sdf <- data.frame(n = names(TCtr), shapes_sphere_invitro_C_TC_internal = as.numeric(TCtr))

shape_df <- merge(shape_df, sdf, by="n", all=TRUE)

shape <- get(load(file = c(file.path(libs_path, libs[4]))))
TCtr <- sapply(shape, function(x){mean(x[2:length(x)]$TC.treated)})
sdf <- data.frame(n = names(TCtr), shape_cell24_invivo_log2ratio_internal = as.numeric(TCtr))

shape_df <- merge(shape_df, sdf, by="n", all=TRUE)

shape <- get(load(file = c(file.path(libs_path, libs[5]))))
TCtr <- sapply(shape, function(x){mean(x[2:length(x)]$TC.treated)})
sdf <- data.frame(n = names(TCtr), shape_cell256_invivo_log2ratio_internal = as.numeric(TCtr))

shape_df <- merge(shape_df, sdf, by="n", all=TRUE)

shape <- get(load(file = c(file.path(libs_path, libs[6]))))
TCtr <- sapply(shape, function(x){mean(x[2:length(x)]$TC.treated)})
sdf <- data.frame(n = names(TCtr), shape_oblong_CHX_invivo_log2ratio_internal = as.numeric(TCtr))

shape_df <- merge(shape_df, sdf, by="n", all=TRUE)

shape <- get(load(file = c(file.path(libs_path, libs[7]))))
TCtr <- sapply(shape, function(x){mean(x[2:length(x)]$TC.treated)})
sdf <- data.frame(n = names(TCtr), shape_oblong_invivo_log2ratio_internal = as.numeric(TCtr))

shape_df <- merge(shape_df, sdf, by="n", all=TRUE)

#matrix_of_everything <- merge(matrix_of_everything, shape_df, by="n", all=TRUE)

