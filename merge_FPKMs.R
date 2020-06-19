### make one FPKM data frame

## RNA
# 1cell
# 1cell (total)
# 2-4cell
# 256cell
# 1Kcell
# 2h (total)
# 4h (total)

## Ribo
# 2-4cell
# 256cell
# 1Kcell

load(file = "/Volumes/USELESS/META/SHAPES/FPKM_24.Rdata")
load(file = "/Volumes/USELESS/META/SHAPES/FPKM_256.Rdata")
load(file = "/Volumes/USELESS/META/SHAPES/FPKM_1K.Rdata")
load(file = "/Volumes/USELESS/META/SHAPES/FPKM_Lee.Rsave")

cell24 <- data.frame(n = rownames(FPKM_24), rna_cell24 = FPKM_24$exons_RNA_fpkm, ribo_cell24 = FPKM_24$exons_Ribo_fpkm)
cell256 <- data.frame(n = rownames(FPKM_256), rna_cell256 = FPKM_256$exons_RNA_fpkm, ribo_cell256 = FPKM_256$exons_Ribo_fpkm)
cell1K <- data.frame(n = rownames(FPKM_1K), rna_cell1K = FPKM_1K$exons_RNA_fpkm, ribo_cell1K = FPKM_1K$exons_Ribo_fpkm)
lee <- data.frame(n = rownames(lee_RNA_FPKM), rna_2h_total = lee_RNA_FPKM$h2, rna_4h_total = lee_RNA_FPKM$h4)

FPKM <- merge(cell24, cell256, by="n", all=TRUE)
FPKM <- merge(FPKM, cell1K, by="n", all=TRUE)
FPKM <- merge(FPKM, lee, by="n", all=TRUE)

save(FPKM, file="/Volumes/USELESS/DATA/Shape-Seq/FPKM.Rdata")

load(file = "/Volumes/USELESS/META/SHAPES/FPKM_1cell_aanes.Rdata")
cell1_fpkm <- data.frame(n = FPKM_1cell_aanes$n, rna_cell1 = FPKM_1cell_aanes$exons_RNA_fpkm)
FPKM <- merge(FPKM, cell1_fpkm, by="n", all=TRUE)
save(FPKM, file="/Volumes/USELESS/DATA/Shape-Seq/FPKM.Rdata")

load(file = "/Volumes/USELESS/META/SHAPES/FPKM_3_5h_aanes.Rdata")
h35 <- data.frame(n = FPKM_3_5h_aanes$n, rna_3_5h = FPKM_3_5h_aanes$exons_RNA_fpkm)
FPKM <- merge(FPKM, h35, by="n", all=TRUE)
save(FPKM, file="/Volumes/USELESS/DATA/Shape-Seq/FPKM.Rdata")

FPKM <- merge(FPKM, cell1, by="n", all=TRUE)
