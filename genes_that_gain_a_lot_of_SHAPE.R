# What are the genes that gain a lot of shape? > 10 on the y-axis?

# get shape_fpkm and GO_term

load(file = "/Volumes/USELESS/META/SHAPES/SHAPE_FPKM.Rdata")

## GO
go <- read.csv("/Users/kasia/Documents/PhD/matrix_go_fq.csv", header = TRUE)
rownames(go) <- go$X.transcript_id
go$n <- rownames(go)
go <- data.frame(n = go$n, GO = go$GO_Term_Name)

#SHAPE_FPKM <- merge(SHAPE_FPKM, go, by = "n", all = TRUE)
#s <- SHAPE_FPKM[complete.cases(SHAPE_FPKM),]

shape24 <- split(gc_unsmoothed_2_4, gc_unsmoothed_2_4$trnames)
shape256 <- split(gc_unsmoothed_256, gc_unsmoothed_256$trnames)

s24 <- sapply(shape24, function(x){mean(x$dtcr)})
names(s24) <- sapply(names(s24), function(x){substr(x,1,18)})
s24 <- data.frame(s24)
s24$n <- rownames(s24)

s256 <- sapply(shape256, function(x){mean(x$dtcr)})
names(s256) <- sapply(names(s256), function(x){substr(x,1,18)})
s256 <- data.frame(s256)
s256$n <- rownames(s256)

shape <- merge(s24, s256, by = "n", all = TRUE)

shape$shape256by24 <- shape$s24 / shape$s256
shape <- shape[complete.cases(shape$shape256by24),]
shape <- shape[is.finite(shape$shape256by24),]
shape <- merge(shape, go, by = "n", all = TRUE)


