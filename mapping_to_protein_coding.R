library(rtracklayer)
library(GenomicFeatures)
library(dplyr)

gtf_annot <- import("/Volumes/USELESS/DATA/genomes/GTF/danio_rerio/Danio_rerio.GRCz10.84.gtf", format = "gtf")
gtf_annot <- gtf_annot[gtf_annot$type == "transcript"]
transcript_df <- data.frame(n = gtf_annot$transcript_id, tx_biotype = gtf_annot$transcript_biotype)

protein_coding <- as.character(subset(transcript_df, tx_biotype == "protein_coding")$n)
protein_coding <- protein_coding[order(protein_coding)]

rRNA <- as.character(subset(transcript_df, tx_biotype == "rRNA")$n)
rRNA <- rRNA[order(rRNA)]

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

### data
load(file = "/Volumes/USELESS/DATA/Shape-Seq/shape_cell24.Rsave")
load(file = "/Volumes/USELESS/DATA/Shape-Seq/shape_cell256.Rsave")
load(file = "/Volumes/USELESS/DATA/Shape-Seq/shape_oblong.Rsave")
load(file = "/Volumes/USELESS/DATA/Shape-Seq/shape_oblong_CHX.Rsave")

shapes_cell24 <- get(load(file = "/Volumes/USELESS/test/normalized/cell2_4_invivo_NAI_N3.Rsave"))
shapes_cell1_invitro <- get(load(file = "/Volumes/USELESS/test/normalized/cell1_invitro_NAI_N3.Rsave"))
shapes_cell1_invitro_nonsel <- get(load(file = "/Volumes/USELESS/test/normalized/cell1_invitro_NAI_N3_nonsel.Rsave"))

## shape_cell24
shape_cell24 <- shape_cell24[names(shape_cell24) %in% txlen$tx_name]
u_shape_cell24 <- unlist(shape_cell24)
sum(u_shape_cell24$TC.control)
sum(u_shape_cell24$TC.treated)
pt_shape_cell24 <- shape_cell24[(names(shape_cell24) %in% protein_coding)]
pt_shape_cell24 <- unlist(pt_shape_cell24)
sum(pt_shape_cell24$TC.control)
sum(pt_shape_cell24$TC.treated)
rrna_shape_cell24 <- shape_cell24[(names(shape_cell24) %in% rRNA)]
rrna_shape_cell24 <- unlist(rrna_shape_cell24)
sum(rrna_shape_cell24$TC.control)
sum(rrna_shape_cell24$TC.treated)

## shape_cell256
shape_cell256 <- shape_cell256[names(shape_cell256) %in% txlen$tx_name]
u_shape_cell256 <- unlist(shape_cell256)
sum(u_shape_cell256$TC.control)
sum(u_shape_cell256$TC.treated)
pt_shape_cell256 <- shape_cell256[(names(shape_cell256) %in% protein_coding)]
pt_shape_cell256 <- unlist(pt_shape_cell256)
sum(pt_shape_cell256$TC.control)
sum(pt_shape_cell256$TC.treated)
rrna_shape_cell256 <- shape_cell256[(names(shape_cell256) %in% rRNA)]
rrna_shape_cell256 <- unlist(rrna_shape_cell256)
sum(rrna_shape_cell256$TC.control)
sum(rrna_shape_cell256$TC.treated)

## shape_oblong
shape_oblong <- shape_oblong[names(shape_oblong) %in% txlen$tx_name]
u_shape_oblong <- unlist(shape_oblong)
sum(u_shape_oblong$TC.control)
sum(u_shape_oblong$TC.treated)
pt_shape_oblong <- shape_oblong[(names(shape_oblong) %in% protein_coding)]
pt_shape_oblong <- unlist(pt_shape_oblong)
sum(pt_shape_oblong$TC.control)
sum(pt_shape_oblong$TC.treated)
rrna_shape_oblong <- shape_oblong[(names(shape_oblong) %in% rRNA)]
rrna_shape_oblong <- unlist(rrna_shape_oblong)
sum(rrna_shape_oblong$TC.control)
sum(rrna_shape_oblong$TC.treated)

## shape_oblong_CHX
shape_oblong_CHX <- shape_oblong_CHX[names(shape_oblong_CHX) %in% txlen$tx_name]
u_shape_oblong_CHX <- unlist(shape_oblong_CHX)
sum(u_shape_oblong_CHX$TC.control)
sum(u_shape_oblong_CHX$TC.treated)
pt_shape_oblong_CHX <- shape_oblong_CHX[(names(shape_oblong_CHX) %in% protein_coding)]
pt_shape_oblong_CHX <- unlist(pt_shape_oblong_CHX)
sum(pt_shape_oblong_CHX$TC.control)
sum(pt_shape_oblong_CHX$TC.treated)
rrna_shape_oblong_CHX <- shape_oblong_CHX[(names(shape_oblong_CHX) %in% rRNA)]
rrna_shape_oblong_CHX <- unlist(rrna_shape_oblong_CHX)
sum(rrna_shape_oblong_CHX$TC.control)
sum(rrna_shape_oblong_CHX$TC.treated)

## shapes_cell24
shapes_cell24 <- shapes_cell24[names(shapes_cell24) %in% txlen$tx_name]
u_shapes_cell24 <- unlist(shapes_cell24)
sum(u_shapes_cell24$TC)
pt_shapes_cell24 <- shapes_cell24[(names(shapes_cell24) %in% protein_coding)]
pt_shapes_cell24 <- unlist(pt_shapes_cell24)
sum(pt_shapes_cell24$TC)
rrna_shapes_cell24 <- shapes_cell24[(names(shapes_cell24) %in% rRNA)]
rrna_shapes_cell24 <- unlist(rrna_shapes_cell24)
sum(rrna_shapes_cell24$TC)

## shapes_cell1_invitro
shapes_cell1_invitro <- shapes_cell1_invitro[names(shapes_cell1_invitro) %in% txlen$tx_name]
u_shapes_cell1_invitro <- unlist(shapes_cell1_invitro)
sum(u_shapes_cell1_invitro$TC)
pt_shapes_cell1_invitro <- shapes_cell1_invitro[(names(shapes_cell1_invitro) %in% protein_coding)]
pt_shapes_cell1_invitro <- unlist(pt_shapes_cell1_invitro)
sum(pt_shapes_cell1_invitro$TC)
rrna_shapes_cell1_invitro <- shapes_cell1_invitro[(names(shapes_cell1_invitro) %in% rRNA)]
rrna_shapes_cell1_invitro <- unlist(rrna_shapes_cell1_invitro)
sum(rrna_shapes_cell1_invitro$TC)

## shapes_cell1_invitro_nonsel
shapes_cell1_invitro_nonsel <- shapes_cell1_invitro_nonsel[names(shapes_cell1_invitro_nonsel) %in% txlen$tx_name]
u_shapes_cell1_invitro_nonsel <- unlist(shapes_cell1_invitro_nonsel)
sum(u_shapes_cell1_invitro_nonsel$TC)
pt_shapes_cell1_invitro_nonsel <- shapes_cell1_invitro_nonsel[(names(shapes_cell1_invitro_nonsel) %in% protein_coding)]
pt_shapes_cell1_invitro_nonsel <- unlist(pt_shapes_cell1_invitro_nonsel)
sum(pt_shapes_cell1_invitro_nonsel$TC)
rrna_shapes_cell1_invitro_nonsel <- shapes_cell1_invitro_nonsel[(names(shapes_cell1_invitro_nonsel) %in% rRNA)]
rrna_shapes_cell1_invitro_nonsel <- unlist(rrna_shapes_cell1_invitro_nonsel)
sum(rrna_shapes_cell1_invitro_nonsel$TC)

#load(file = "/Volumes/USELESS/DATA/Shape-Seq/shape_cell1K.Rsave")
#load(file = "/Volumes/USELESS/DATA/Shape-Seq/shape_cell1.Rsave")

load(file = "/Volumes/USELESS/test2/normalized/cell1_invivo_et1.Rsave")
## shapeseq_norm
shapeseq_norm <- shapeseq_norm[names(shapeseq_norm) %in% txlen$tx_name]
u_shapeseq_norm <- unlist(shapeseq_norm)
sum(u_shapeseq_norm$TC.control)
sum(u_shapeseq_norm$TC.treated)
pt_shapeseq_norm <- shapeseq_norm[(names(shapeseq_norm) %in% protein_coding)]
pt_shapeseq_norm <- unlist(pt_shapeseq_norm)
sum(pt_shapeseq_norm$TC.control)
sum(pt_shapeseq_norm$TC.treated)
rrna_shapeseq_norm <- shapeseq_norm[(names(shapeseq_norm) %in% rRNA)]
rrna_shapeseq_norm <- unlist(rrna_shapeseq_norm)
sum(rrna_shapeseq_norm$TC.control)
sum(rrna_shapeseq_norm$TC.treated)

load(file = "/Volumes/USELESS/test2/normalized/cell1_invivo_et2.Rsave")
load(file = "/Volumes/USELESS/test2/normalized/cell1K_invivo_et1.Rsave")
load(file = "/Volumes/USELESS/test2/normalized/cell1K_invivo_et2.Rsave")



############
load(file = "/Volumes/USELESS/test2/normalized/RAW_cell1_invitro_NAI_et1.Rsave")
treated_comp <- split(treated_comp, names(treated_comp))
treated_comp <- treated_comp[names(treated_comp) %in% txlen$tx_name]
u_treated_comp <- unlist(treated_comp)
sum(u_treated_comp$TC)
pt_treated_comp <- treated_comp[(names(treated_comp) %in% protein_coding)]
pt_treated_comp <- unlist(pt_treated_comp)
sum(pt_treated_comp$TC)
rrna_treated_comp <- treated_comp[(names(treated_comp) %in% rRNA)]
rrna_treated_comp <- unlist(rrna_treated_comp)
sum(rrna_treated_comp$TC)

load(file = "/Volumes/USELESS/test2/normalized/RAW_cell1_invitro_NAI_et2.Rsave")

##################################
###### OLD and slow ##############
##################################

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



############
Reduce("+", sapply(shape_cell24, function(x){sum(x$TC.control)}))        
Reduce("+", sapply(shape_cell24, function(x){sum(x$TC.treated)}))        
Reduce("+", sapply(shape_cell256, function(x){sum(x$TC.control)}))       
Reduce("+", sapply(shape_cell256, function(x){sum(x$TC.treated)}))       
Reduce("+", sapply(shape_oblong, function(x){sum(x$TC.control)}))        
Reduce("+", sapply(shape_oblong, function(x){sum(x$TC.treated)}))        
Reduce("+", sapply(shape_oblong_CHX, function(x){sum(x$TC.control)}))    
Reduce("+", sapply(shape_oblong_CHX, function(x){sum(x$TC.treated)}))    


Reduce("+", sapply(rrna_shape_cell24, function(x){sum(x$TC.control)}))        
Reduce("+", sapply(rrna_shape_cell24, function(x){sum(x$TC.treated)}))        
Reduce("+", sapply(rrna_shape_cell256, function(x){sum(x$TC.control)}))       
Reduce("+", sapply(rrna_shape_cell256, function(x){sum(x$TC.treated)}))       
Reduce("+", sapply(rrna_shape_oblong, function(x){sum(x$TC.control)}))        
Reduce("+", sapply(rrna_shape_oblong, function(x){sum(x$TC.treated)}))        
Reduce("+", sapply(rrna_shape_oblong_CHX, function(x){sum(x$TC.control)}))    
Reduce("+", sapply(rrna_shape_oblong_CHX, function(x){sum(x$TC.treated)}))


Reduce("+", sapply(shapes_cell24, function(x){sum(x$TC)}))               
Reduce("+", sapply(shapes_cell1_invitro, function(x){sum(x$TC)}))        ??
Reduce("+", sapply(shapes_cell1_invitro_nonsel, function(x){sum(x$TC)})) 

Reduce("+", sapply(rrna_shapes_cell24, function(x){sum(x$TC)}))               
Reduce("+", sapply(rrna_shapes_cell1_invitro, function(x){sum(x$TC)}))        
Reduce("+", sapply(rrna_shapes_cell1_invitro_nonsel, function(x){sum(x$TC)})) 

### get one tx per gene, sum TC again

################
library(rtracklayer)
library(GenomicFeatures)

# shape
shape_cell1K_et1 <- get(load(file = "/Volumes/USELESS/test2/normalized/cell1K_invivo_et1.Rsave"))
shape_cell1K_et2 <- get(load(file = "/Volumes/USELESS/test2/normalized/cell1K_invivo_et2.Rsave"))
shape_cell1_et1 <- get(load(file = "/Volumes/USELESS/test2/normalized/cell1_invivo_et1.Rsave"))
shape_cell1_et2 <- get(load(file = "/Volumes/USELESS/test2/normalized/cell1_invivo_et2.Rsave"))
rm(shapeseq_norm)

tc_cell1_invitro_et1 <- get(load(file = "/Volumes/USELESS/test2/normalized/RAW_cell1_invitro_NAI_et1.Rsave"))
tc_cell1_invitro_et1 <- split(tc_cell1_invitro_et1, seqnames(tc_cell1_invitro_et1))
tc_cell1_invitro_et1 <- tc_cell1_invitro_et1[sapply(tc_cell1_invitro_et1, function(x){length(x) > 0})]
tc_cell1_invitro_et2 <- get(load(file = "/Volumes/USELESS/test2/normalized/RAW_cell1_invitro_NAI_et2.Rsave"))
tc_cell1_invitro_et2 <- split(tc_cell1_invitro_et2, seqnames(tc_cell1_invitro_et2))
tc_cell1_invitro_et2 <- tc_cell1_invitro_et2[sapply(tc_cell1_invitro_et2, function(x){length(x) > 0})]
rm(treated_comp)

gtf_annot <- import("/Volumes/USELESS/DATA/genomes/GTF/danio_rerio/Danio_rerio.GRCz10.84.gtf", format = "gtf")
gtf_annot <- gtf_annot[gtf_annot$type == "transcript"]
transcript_df <- data.frame(n = gtf_annot$transcript_id, tx_biotype = gtf_annot$transcript_biotype)

protein_coding <- as.character(subset(transcript_df, tx_biotype == "protein_coding")$n)
protein_coding <- protein_coding[order(protein_coding)]


pt_shape_cell1K_et1 <- shape_cell1K_et1[names(shape_cell1K_et1) %in% protein_coding]
pt_shape_cell1K_et2 <- shape_cell1K_et2[names(shape_cell1K_et2) %in% protein_coding]
pt_shape_cell1_et1 <- shape_cell1_et1[names(shape_cell1_et1) %in% protein_coding]
pt_shape_cell1_et2 <- shape_cell1_et2[names(shape_cell1_et2) %in% protein_coding]

pt_tc_cell1_invitro_et1 <- tc_cell1_invitro_et1[names(tc_cell1_invitro_et1) %in% protein_coding]
pt_tc_cell1_invitro_et2 <- tc_cell1_invitro_et2[names(tc_cell1_invitro_et2) %in% protein_coding]



Reduce("+", sapply(pt_shape_cell1K_et1, function(x){sum(x$TC.control)}))        # 571987
Reduce("+", sapply(pt_shape_cell1K_et1, function(x){sum(x$TC.treated)}))        # 559194
Reduce("+", sapply(pt_shape_cell1K_et2, function(x){sum(x$TC.control)}))        # 569169
Reduce("+", sapply(pt_shape_cell1K_et2, function(x){sum(x$TC.treated)}))        # 294462
Reduce("+", sapply(pt_shape_cell1_et1, function(x){sum(x$TC.control)}))         # 632694
Reduce("+", sapply(pt_shape_cell1_et1, function(x){sum(x$TC.treated)}))         # 436589
Reduce("+", sapply(pt_shape_cell1_et2, function(x){sum(x$TC.control)}))         # 680845
Reduce("+", sapply(pt_shape_cell1_et2, function(x){sum(x$TC.treated)}))         # 1158475

Reduce("+", sapply(pt_tc_cell1_invitro_et1, function(x){sum(x$TC)}))            # 418732
Reduce("+", sapply(pt_tc_cell1_invitro_et2, function(x){sum(x$TC)}))            # 616902


#######
rRNA <- as.character(subset(transcript_df, tx_biotype == "rRNA")$n)
rRNA <- rRNA[order(rRNA)]
rrna_shape_cell1K_et1 <- shape_cell1K_et1[names(shape_cell1K_et1) %in% rRNA]
rrna_shape_cell1_et1 <- shape_cell1_et1[names(shape_cell1_et1) %in% rRNA]
Reduce("+", sapply(rrna_shape_cell1K_et1, function(x){sum(x$TC.treated)})) # 2229179
Reduce("+", sapply(rrna_shape_cell1_et1, function(x){sum(x$TC.treated)})) # 867669


notpt <- as.character(subset(transcript_df, tx_biotype != "protein_coding")$n)
notpt <- notpt[order(notpt)]
notpt_shape_cell1K_et1 <- shape_cell1K_et1[names(shape_cell1K_et1) %in% notpt]
notpt_shape_cell1_et1 <- shape_cell1_et1[names(shape_cell1_et1) %in% notpt]
Reduce("+", sapply(notpt_shape_cell1K_et1, function(x){sum(x$TC.treated)}))
Reduce("+", sapply(notpt_shape_cell1_et1, function(x){sum(x$TC.treated)}))


########
rrna_tc_cell1_invitro_et1 <- tc_cell1_invitro_et1[names(tc_cell1_invitro_et1) %in% rRNA]
rrna_tc_cell1_invitro_et2 <- tc_cell1_invitro_et2[names(tc_cell1_invitro_et2) %in% rRNA]
Reduce("+", sapply(tc_cell1_invitro_et1, function(x){sum(x$TC)})) 
Reduce("+", sapply(tc_cell1_invitro_et2, function(x){sum(x$TC)}))
Reduce("+", sapply(rrna_tc_cell1_invitro_et1, function(x){sum(x$TC)})) 
Reduce("+", sapply(rrna_tc_cell1_invitro_et2, function(x){sum(x$TC)}))


#######
libs_path <- "/Volumes/USELESS/DATA/Shape-Seq/shapes_Dec2017"
selA <- get(load(file = c(file.path(libs_path, "PAIR_GC_SelA.Rsave"))))
selC <- get(load(file = c(file.path(libs_path, "PAIR_GC_SelC.Rsave"))))
rm(GR_treated)

## selA/C
shape <- selC[names(selC) %in% txlen$tx_name]
u_shape <- unlist(shape)
sum(u_shape$TC.treated)
pt_shape <- shape[(names(shape) %in% protein_coding)]
pt_shape <- unlist(pt_shape)
sum(pt_shape$TC.treated)
rrna_shape <- shape[(names(shape) %in% rRNA)]
rrna_shape <- unlist(rrna_shape)
sum(rrna_shape$TC.treated)


############## remapped, all single-end

libs_path <- "/Volumes/USELESS/DATA/SHAPES/SE"
libs <- c("norm2_cell1_invitro_NAIN3.Rsave", "norm2_cell24_invivo.Rsave", "norm2_cell256_invivo.Rsave",
          "norm2_oblong_CHX_invivo.Rsave", "norm2_oblong_invivo.Rsave", "norm2_sphere_invitro_A.Rsave", "norm2_sphere_invitro_C.Rsave")


for (i in 1:length(libs)) {
  print(libs[i])
  shape <- get(load(file = c(file.path(libs_path, libs[i]))))
  shape <- shape[names(shape) %in% txlen$tx_name]
  u_shape <- unlist(shape)
  print(sum(u_shape$TC.treated))
  print(sum(u_shape$TC.control))
  pt_shape <- shape[(names(shape) %in% protein_coding)]
  pt_shape <- unlist(pt_shape)
  print(sum(pt_shape$TC.treated))
  print(sum(pt_shape$TC.control))
}

############## June 2018

libs_path <- "/Volumes/USELESS/DATA/SHAPES/SE/June2018"
libs <- list.files(path = libs_path)
libs <- libs[10:18]

for (i in 1:length(libs)) {
  print(libs[i])
  shape <- get(load(file = c(file.path(libs_path, libs[i]))))
  shape <- shape[names(shape) %in% txlen$tx_name]
  u_shape <- unlist(shape)
  print(sum(u_shape$TC.treated))
  pt_shape <- shape[(names(shape) %in% protein_coding)]
  pt_shape <- unlist(pt_shape)
  print(sum(pt_shape$TC.treated))
}


