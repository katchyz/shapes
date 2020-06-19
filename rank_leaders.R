### a ranked list of GRCz10 transcripts with the most to least 5' leader structure
library(GenomicFeatures)
library(plyr)
library(ggplot2)

## norm2*.Rsave
libs_path <- "/Volumes/USELESS/DATA/SHAPES/normalized"
#libs <- list.files(path = libs_path, pattern = "^norm2")

cell256 <- c("norm2_cell256_HP_NAIN3.Rsave", "norm2_cell256_NAIN3_old.Rsave", "norm2_cell256_P_NAIN3.Rsave",
          "norm2_cell256_invivo.Rsave")
oblong <- c("norm2_oblong_CHX_invivo.Rsave", "norm2_oblong_NAIN3.Rsave", "norm2_oblong_NAIN3_old.Rsave",
            "norm2_oblong_invivo.Rsave")


## get utr5 length for subsetting
txdb_can <- makeTxDbFromGFF("/Volumes/USELESS/DATA/genomes/GTF/danio_rerio/Danio_rerio.GRCz10.84.gtf", format = "gtf")
txLengths <- transcriptLengths(txdb_can, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE)
rownames(txLengths) <- txLengths$tx_name
txLengths <- txLengths[order(rownames(txLengths)),]

# get longest transcript per gene
txlen <- arrange(txLengths, gene_id, desc(tx_len))
txlen <- txlen[!duplicated(txlen$gene_id),]

### tx with utr5 longer than 30
txlen <- txlen[txlen$utr5_len > 30,]
rownames(txlen) <- txlen$tx_name
txlen <- txlen[order(rownames(txlen)),]

leaders <- GRanges(Rle(rownames(txlen), rep(1, nrow(txlen))), IRanges(rep(1, nrow(txlen)), txlen$utr5_len))
leaders <- split(leaders, seqnames(leaders))

## RNA-Seq
load(file="/Volumes/USELESS/META/SHAPES/matrix_of_everything.Rdata")
rna_cell256 <- data.frame(n = matrix_of_everything$n, rna = matrix_of_everything$rna_cell256)
rna_oblong <- data.frame(n = matrix_of_everything$n, rna = matrix_of_everything$rna_3_5h)

for (i in 1:length(cell256)) {
  u <- get(load(file = c(file.path(libs_path, cell256[i]))))
  u <- unlist(u)
  names(u) <- substring(names(u), 1, 18)
  gr <- GRanges(seqnames=Rle(names(u)),
                IRanges(start(u), width=1),
                strand='*')
  gr$TC.treated <- u$TC.treated
  l <- subsetByOverlaps(gr, leaders, ignore.strand=TRUE)
  l <- split(l, seqnames(l))
  l <- l[sapply(l, function(x){length(x) > 0})]
  l <- l[sapply(l, function(x){sum(x$TC.treated) > 0})]
  ## normalize per RNA-seq
  rna <- rna_cell256[rna_cell256$n %in% names(l),]
  if (all(names(l) == rna$n)) {
    ## mean structure
    mean_structure <- sapply(l, function(x){mean(x$TC.treated)})
    mean_structure_norm <- mean_structure / rna$rna
    mean_structure_norm <- sort(mean_structure_norm, decreasing = TRUE)
    mean_structure_norm <- mean_structure_norm[mean_structure_norm != Inf]
    file <- c(file.path("/Volumes/USELESS/META/SHAPES/ranked_leaders", paste("ranked_leaders_", substr(cell256[i], 7, (nchar(cell256[i]))), sep="")))
    save(mean_structure_norm, file = file)
  }
}

for (i in 1:length(oblong)) {
  u <- get(load(file = c(file.path(libs_path, oblong[i]))))
  u <- unlist(u)
  names(u) <- substring(names(u), 1, 18)
  gr <- GRanges(seqnames=Rle(names(u)),
                IRanges(start(u), width=1),
                strand='*')
  gr$TC.treated <- u$TC.treated
  l <- subsetByOverlaps(gr, leaders, ignore.strand=TRUE)
  l <- split(l, seqnames(l))
  l <- l[sapply(l, function(x){length(x) > 0})]
  l <- l[sapply(l, function(x){sum(x$TC.treated) > 0})]
  ## normalize per RNA-seq
  rna <- rna_oblong[rna_oblong$n %in% names(l),]
  if (all(names(l) == rna$n)) {
    ## mean structure
    mean_structure <- sapply(l, function(x){mean(x$TC.treated)})
    mean_structure_norm <- mean_structure / rna$rna
    mean_structure_norm <- sort(mean_structure_norm, decreasing = TRUE)
    mean_structure_norm <- mean_structure_norm[mean_structure_norm != Inf]
    file <- c(file.path("/Volumes/USELESS/META/SHAPES/ranked_leaders", paste("ranked_leaders_", substr(oblong[i], 7, (nchar(oblong[i]))), sep="")))
    save(mean_structure_norm, file = file)
  }
}


### check if they correlate
fpath <- "/Volumes/USELESS/META/SHAPES/ranked_leaders"
libs <- list.files(path = fpath)


cell1 <- get(load(file = c(file.path(fpath, libs[1]))))
cell1 <- data.frame(n = names(cell1), shape_cell1 = cell1)
cell2 <- get(load(file = c(file.path(fpath, libs[2]))))
cell2 <- data.frame(n = names(cell2), shape_cell2 = cell2)
cell3 <- get(load(file = c(file.path(fpath, libs[3]))))
cell3 <- data.frame(n = names(cell3), shape_cell3 = cell3)
cell4 <- get(load(file = c(file.path(fpath, libs[4]))))
cell4 <- data.frame(n = names(cell4), shape_cell4 = cell4)

oblong1 <- get(load(file = c(file.path(fpath, libs[5])))))
oblong1 <- data.frame(n = names(oblong1), shape_oblong1 = oblong1)
oblong2 <- get(load(file = c(file.path(fpath, libs[6]))))
oblong2 <- data.frame(n = names(oblong2), shape_oblong2 = oblong2)
oblong3 <- get(load(file = c(file.path(fpath, libs[7]))))
oblong3 <- data.frame(n = names(oblong3), shape_oblong3 = oblong3)
oblong4 <- get(load(file = c(file.path(fpath, libs[8]))))
oblong4 <- data.frame(n = names(oblong4), shape_oblong4 = oblong4)

df <- cell1
df <- merge(df, cell2, by="n", all=TRUE)
df <- merge(df, cell3, by="n", all=TRUE)
df <- merge(df, cell4, by="n", all=TRUE)
df <- merge(df, oblong1, by="n", all=TRUE)
df <- merge(df, oblong2, by="n", all=TRUE)
df <- merge(df, oblong3, by="n", all=TRUE)
df <- merge(df, oblong4, by="n", all=TRUE)

ggplot(df, aes(shape_cell2, shape_cell4)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()
ggplot(df, aes(shape_oblong2, shape_oblong3)) + geom_point(size=0.1) + scale_x_log10() + scale_y_log10()

