### SLOGRAT

load(file = "/Volumes/USELESS/OUT/SHAPES_out/slograt_unsmoothed24.Rsave")
load(file = "/Volumes/USELESS/OUT/SHAPES_out/slograt_unsmoothed256.Rsave")

slog24 <- split(slograt_unsmoothed24, seqnames(slograt_unsmoothed24))
slog256 <- split(slograt_unsmoothed256, seqnames(slograt_unsmoothed256))

names(slog24) <- sapply(names(slog24), function(x){substr(x,1,18)})
names(slog256) <- sapply(names(slog256), function(x){substr(x,1,18)})

s24 <- sapply(slog24, function(x){mean(x$slograt)})
s256 <- sapply(slog256, function(x){mean(x$slograt)})

### RNA-Seq & Ribo-Seq
# FPKM_24
load(file = "/Volumes/USELESS/META/SHAPES/FPKM_24.Rdata")
# FPKM_256
load(file = "/Volumes/USELESS/META/SHAPES/FPKM_256.Rdata")

s <- s256[names(s256) %in% rownames(FPKM_256)]
f <- FPKM_256[rownames(FPKM_256) %in% names(s),]
f$slograt <- s

ggplot(f, aes(x = exons_RNA_fpkm, y = slograt)) + geom_point() + scale_x_log10() + scale_y_log10() + ggtitle("256cell")
ggsave(file = "/Volumes/USELESS/META/SHAPES/general/slograt/txFPKM_vs_mean_slograt256.png")
