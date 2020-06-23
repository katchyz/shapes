library(GenomicFeatures)

txdb_can <- makeTxDbFromGFF("/Home/ii/katchyz/DATA/genomes/GTF/Danio_rerio.GRCz10.84.gtf", format = "gtf")
exons <- exonsBy(txdb_can, by="tx", use.names=TRUE)

load(file = "/export/valenfs/data/raw_data/SHAPES/normalized/u_2-4cell.Rsave")
gc_shape_cell24 <- mapFromTranscripts(u, exons, ignore.strand = FALSE)

save(gc_shape_cell24, file = "/export/valenfs/data/raw_data/SHAPES/normalized/gc_2-4cell.Rsave")

######

library(GenomicFeatures)

txdb_can <- makeTxDbFromGFF("/Home/ii/katchyz/DATA/genomes/GTF/Danio_rerio.GRCz10.84.gtf", format = "gtf")
exons <- exonsBy(txdb_can, by="tx", use.names=TRUE)

load(file = "/export/valenfs/data/raw_data/SHAPES/normalized/u_oblong_CHX_invivo.Rsave")
gc_shape_oblong_CHX <- mapFromTranscripts(u, exons, ignore.strand = FALSE)

save(gc_shape_oblong_CHX, file = "/export/valenfs/data/raw_data/SHAPES/normalized/gc_oblong_CHX.Rsave")

