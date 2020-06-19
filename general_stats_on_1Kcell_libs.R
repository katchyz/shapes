##### general stats on 1Kcell libraries #####
### load shape-seq
library(rtracklayer)

GR_RZ <- get(load(file = "/Volumes/USELESS/DATA/Shape-Seq/1Kcell/GC_1Kcell_RZ.Rsave"))
GR_RZ_vitro <- get(load(file = "/Volumes/USELESS/DATA/Shape-Seq/1Kcell/GC_1Kcell_RZ_vitro.Rsave"))
GR_polyA <- get(load(file = "/Volumes/USELESS/DATA/Shape-Seq/1Kcell/GC_1Kcell_polyA.Rsave"))
GR_polyA_vitro <- get(load(file = "/Volumes/USELESS/DATA/Shape-Seq/1Kcell/GC_1Kcell_polyA_vitro.Rsave"))
rm(GR)

## get number of 5'ends (sum TC.control and TC.treated)

sum(GR_RZ$TC.control)
sum(GR_RZ$TC.treated)
sum(GR_RZ_vitro$TC.control)
sum(GR_RZ_vitro$TC.treated)
sum(GR_polyA$TC.control)
sum(GR_polyA$TC.treated)
sum(GR_polyA_vitro$TC.control)
sum(GR_polyA_vitro$TC.treated)

### mapping to protein-coding
gtf_annot <- import("/Volumes/USELESS/DATA/genomes/GTF/danio_rerio/Danio_rerio.GRCz10.84.gtf", format = "gtf")
gtf_annot <- gtf_annot[gtf_annot$type == "transcript"]
transcript_df <- data.frame(n = gtf_annot$transcript_id, tx_biotype = gtf_annot$transcript_biotype)

protein_coding <- as.character(subset(transcript_df, tx_biotype == "protein_coding")$n)
protein_coding <- protein_coding[order(protein_coding)]


###### get just the longest transcripts - this way I'm not counting TC multiple times
txdb_can <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/danio_rerio/Danio_rerio.GRCz10.84.gtf", format = "gtf")
txLengths <- transcriptLengths(txdb_can, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
rownames(txLengths) <- txLengths$tx_name
txLengths <- txLengths[order(rownames(txLengths)),]
rownames(txLengths) <- sapply(rownames(txLengths), function(x){substr(x,1,18)})
# get longest transcript per gene
txlen <- arrange(txLengths, gene_id, desc(tx_len))
txlen <- txlen[!duplicated(txlen$gene_id),]
## txlen$tx_name

grl_RZ <- split(GR_RZ, names(GR_RZ))
grl_RZ <- grl_RZ[names(grl_RZ) %in% txlen$tx_name]
pt_RZ <- grl_RZ[(names(grl_RZ) %in% protein_coding)]
pt_RZ <- unlist(pt_RZ)
sum(pt_RZ$TC.control)
sum(pt_RZ$TC.treated)

grl_RZ_vitro <- split(GR_RZ_vitro, names(GR_RZ_vitro))
grl_RZ_vitro <- grl_RZ_vitro[names(grl_RZ_vitro) %in% txlen$tx_name]
pt_RZ_vitro <- grl_RZ_vitro[(names(grl_RZ_vitro) %in% protein_coding)]
pt_RZ_vitro <- unlist(pt_RZ_vitro)
sum(pt_RZ_vitro$TC.control)
sum(pt_RZ_vitro$TC.treated)

grl_polyA <- split(GR_polyA, names(GR_polyA))
grl_polyA <- grl_polyA[names(grl_polyA) %in% txlen$tx_name]
pt_polyA <- grl_polyA[(names(grl_polyA) %in% protein_coding)]
pt_polyA <- unlist(pt_polyA)
sum(pt_polyA$TC.control)
sum(pt_polyA$TC.treated)

grl_polyA_vitro <- split(GR_polyA_vitro, names(GR_polyA_vitro))
grl_polyA_vitro <- grl_polyA_vitro[names(grl_polyA_vitro) %in% txlen$tx_name]
pt_polyA_vitro <- grl_polyA_vitro[(names(grl_polyA_vitro) %in% protein_coding)]
pt_polyA_vitro <- unlist(pt_polyA_vitro)
sum(pt_polyA_vitro$TC.control)
sum(pt_polyA_vitro$TC.treated)




### rRNA
rRNA <- as.character(subset(transcript_df, tx_biotype == "rRNA")$n)
rRNA <- rRNA[order(rRNA)]

notpt <- as.character(subset(transcript_df, tx_biotype != "protein_coding")$n)
notpt <- notpt[order(notpt)]

grl_RZ <- split(GR_RZ, names(GR_RZ))
rrna_RZ <- grl_RZ[(names(grl_RZ) %in% rRNA)]
rrna_RZ <- unlist(rrna_RZ)
sum(rrna_RZ$TC.control)
sum(rrna_RZ$TC.treated)

grl_RZ_vitro <- split(GR_RZ_vitro, names(GR_RZ_vitro))
rrna_RZ_vitro <- grl_RZ_vitro[(names(grl_RZ_vitro) %in% rRNA)]
rrna_RZ_vitro <- unlist(rrna_RZ_vitro)
sum(rrna_RZ_vitro$TC.control)
sum(rrna_RZ_vitro$TC.treated)

grl_polyA <- split(GR_polyA, names(GR_polyA))
rrna_polyA <- grl_polyA[(names(grl_polyA) %in% rRNA)]
rrna_polyA <- unlist(rrna_polyA)
sum(rrna_polyA$TC.control)
sum(rrna_polyA$TC.treated)

grl_polyA_vitro <- split(GR_polyA_vitro, names(GR_polyA_vitro))
rrna_polyA_vitro <- grl_polyA_vitro[(names(grl_polyA_vitro) %in% rRNA)]
rrna_polyA_vitro <- unlist(rrna_polyA_vitro)
sum(rrna_polyA_vitro$TC.control)
sum(rrna_polyA_vitro$TC.treated)
