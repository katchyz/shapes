diNT <- read.csv("/Volumes/USELESS/META/SHAPES/diNTrelAdf.csv")
rownames(diNT) <- diNT$X

save(diNT, file = "/Volumes/USELESS/META/SHAPES/diNT.Rsave")

dn <- diNT[rownames(diNT) %in% rownames(shape_256),]

df <- melt(dn, index = shape)
df$shape <- rep(shape_256$shape_cds, 48)

p <- ggplot(df, aes(x = value, y = shape)) + geom_point() + facet_wrap(~ variable) + scale_x_log10() + scale_y_log10()
ggsave("/Volumes/USELESS/META/SHAPES/diNT/matrix.png", plot = p)

##############
# diNT vs ORFscore

dn <- diNT[rownames(diNT) %in% names(orf1_256),]
df <- melt(dn)
df$orf1 <- rep(orf1_256, 48)
ggplot(df, aes(x = value, y = orf1)) + geom_point() + facet_wrap(~ variable) + scale_x_log10() + scale_y_log10()

############
# split by high/low TE (or just colour on the plot)
dn <- diNT[rownames(diNT) %in% names(orf3_256),]
dn <- dn[rownames(dn) %in% rownames(shape_256),]

sh <- shape_256[rownames(shape_256) %in% rownames(dn),]
orf <- orf3_256[names(orf3_256) %in% rownames(dn)]
df <- melt(dn)
df$orf3 <- rep(orf, 48)
df$te <- rep(sh$te_cds, 48)

high_te <- subset(df, te > 5.56)
low_te <- subset(df, te < 0.16)

ggplot(high_te, aes(x = value, y = orf3)) + geom_point() + facet_wrap(~ variable) + scale_x_log10() + scale_y_log10()
ggplot(low_te, aes(x = value, y = orf3)) + geom_point() + facet_wrap(~ variable) + scale_x_log10() + scale_y_log10()

######
# boxplots: high TE / low TE
# 5utr and 3utr

### GC content
GCcontent <- read.csv("/Volumes/USELESS/META/SHAPES/GCdf.csv")
GCcontent$X <- NULL
rownames(GCcontent) <- GCcontent$X0

sh <- shape_256[rownames(shape_256) %in% rownames(gc),]
gc <- GCcontent[rownames(GCcontent) %in% rownames(sh),]
gc$shape <- sh$shape_cds

ggplot(gc, aes(x = X1, y = shape)) + geom_point() + scale_x_log10() + scale_y_log10()
