### few tx with high coverage (in vitro / in vivo; across stages)
load(file="/Volumes/USELESS/META/SHAPES_NEW/general/clouds/cloud_df.Rsave")

# shape
load(file="/Volumes/USELESS/DATA/Shape-Seq/shape_cell1.Rsave")
load(file="/Volumes/USELESS/DATA/Shape-Seq/shape_cell24.Rsave")
load(file="/Volumes/USELESS/DATA/Shape-Seq/shape_cell256.Rsave")
load(file="/Volumes/USELESS/DATA/Shape-Seq/shape_cell1K.Rsave")
load(file="/Volumes/USELESS/DATA/Shape-Seq/shape_oblong.Rsave")
load(file="/Volumes/USELESS/DATA/Shape-Seq/shape_oblong_CHX.Rsave")

# shapes
shapes_cell1_invitro <- get(load(file = "/Volumes/USELESS/test/normalized/cell1_invitro_NAI_N3.Rsave"))
shapes_cell24 <- get(load(file = "/Volumes/USELESS/test/normalized/cell2_4_invivo_NAI_N3.Rsave"))
rm(shapes_norm)

# get max coverage (or highest number non-zero elements)

sums_cell1 <- sapply(shape_cell1, function(x){sum(x$log2ratio)})
sums_cell24 <- sapply(shape_cell24, function(x){sum(x$log2ratio)})
sums_cell256 <- sapply(shape_cell256, function(x){sum(x$log2ratio)})
sums_cell1K <- sapply(shape_cell1K, function(x){sum(x$log2ratio)})
sums_oblong <- sapply(shape_oblong, function(x){sum(x$log2ratio)})
sums_oblong_CHX <- sapply(shape_oblong_CHX, function(x){sum(x$log2ratio)})

sums_cell1_invitro <- sapply(shapes_cell1_invitro, function(x){sum(x$TC)})
sums_cell24_TC <- sapply(shapes_cell24, function(x){sum(x$TC)})

# gene names
gn <- read.csv(file="/Volumes/USELESS/DATA/genomes/ensembl_gene_name.txt", header=T)

top <- names(tail(sort(sums_cell1_invitro),100))[names(tail(sort(sums_cell1_invitro),100)) %in% names(tail(sort(sums_cell1),100))]
gn[gn$Transcript.stable.ID %in% top,]

head(df[with(df, order(-shape_cell1_log2ratio_internal,shapes_cell1_invitro_TC_internal)), ])

##
ndf <- df[(log2(df$ribo_cell24) - log2(df$ribo_cell256)) < -2,]
ndf <- ndf[!is.na(ndf$gene),]
ndf <- ndf[ndf$ribo_cell24 > 100 & ndf$ribo_cell256 > 100,]


top_cell1 <- unique(ndf[with(ndf, order(-shape_cell1_log2ratio_internal)),]$gene[1:100])
top_cell1_invitro <- unique(ndf[with(ndf, order(-shapes_cell1_invitro_TC_internal)),]$gene[1:100])
#top <- top_cell1[top_cell1 %in% top_cell1_invitro]

top_cell24 <- unique(ndf[with(ndf, order(-shape_cell24_log2ratio_internal)),]$gene[1:100])
top_cell256 <- unique(ndf[with(ndf, order(-shape_cell256_log2ratio_internal)),]$gene[1:100])
top_cell1K <- unique(ndf[with(ndf, order(-shape_cell1K_log2ratio_internal)),]$gene[1:100])
top_oblong <- unique(ndf[with(ndf, order(-shape_oblong_log2ratio_internal)),]$gene[1:100])
top_oblong_CHX <- unique(ndf[with(ndf, order(-shape_oblong_CHX_log2ratio_internal)),]$gene[1:100])

top_cell1[top_cell1 %in% top_cell24]

