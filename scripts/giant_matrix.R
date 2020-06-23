### giant matrix with everything!
library(GenomicFeatures)
library(seqinr)

## libraries:
## "/Volumes/USELESS/DATA/Shape-Seq/*.Rsave"
# 1cell invitro
# 1cell invivo
# 2-4cell
# 256cell
# 1Kcell
# oblong
# 2h, 4h (PAL-Seq, totalRNA)

# FPKM (RNA, Ribo)
load(file="/Volumes/USELESS/DATA/Shape-Seq/FPKM.Rdata")

## data types:
# SHAPE-Seq
# RNA-Seq
# Ribo-Seq
# PAL-Seq
pal2h <- read.csv(file="/Volumes/USELESS/DATA/PAL-Seq/subtelny_2014/Dre_mock_2hpf_absolute_no_intron_es1.csv", header = T, sep = "\t")
colnames(pal2h) <- c("n", "tail_lengths_2h", "tag_count_2h", "mean_length_2h")
pal4h <- read.csv(file="/Volumes/USELESS/DATA/PAL-Seq/subtelny_2014/Dre_mock_4hpf_absolute_no_intron_es1.csv", header = T, sep = "\t")
colnames(pal4h) <- c("n", "tail_lengths_4h", "tag_count_4h", "mean_length_4h")
pal <- merge(pal2h, pal4h, by="n", all=TRUE)
pal <- data.frame(n = pal$n, pal_2h = pal$mean_length_2h, pal_4h = pal$mean_length_4h)
pal <- pal[pal$n %in% matrix_of_everything$n,]

## ensembl_ID, gene_name
ensembl_gn <- read.csv(file="/Volumes/USELESS/DATA/genomes/ensembl_gene_name.txt", header = T)
## tx_features: lengths, (UTRs/CDS), CAI, stop_codon, GO_cc, _bp, _mf
## CAI
cai <- read.csv("/Volumes/USELESS/META/SHAPES_OLD/caiDf.csv")
cai$X <- NULL
colnames(cai) <- c('n', 'CAI')
cai <- cai[cai$n %in% matrix_of_everything$n,]

df <- merge(df, cai, by="n", all=TRUE)
for (i in 1:11) {
  ggplot(df, aes(x = log2(df[Xaes[i]]), y = CAI)) + geom_point(size=0.1) + theme(legend.position="none") + xlab(Xaes[i]) + ylim(0.65,0.95)
  fp <- paste0("/Volumes/USELESS/META/SHAPES_NEW/CAI/p", i, ".png", collapse = NULL)
  ggsave(fp)
}

## stop codon
fasta_cds <- read.fasta("/Volumes/USELESS/DATA/fasta/danio_rerio/GRCz10/cds/Danio_rerio.GRCz10.cds.all.fa")
names(fasta_cds) <- sapply(names(fasta_cds), function(x){substr(x,1,18)})
fasta_cds <- fasta_cds[unlist(lapply(fasta_cds, function(x){length(x) %% 3 == 0}))]

stop_codons <- lapply(fasta_cds, function(x){paste(x[(length(x)-2):length(x)], collapse = "")})
stop_codons <- data.frame(n = names(stop_codons), stop_codon = as.vector(unlist(stop_codons)))
stop_codons <- stop_codons[stop_codons$stop_codon == "taa" | stop_codons$stop_codon == "tga" | stop_codons$stop_codon == "tag",]
stop_codons <- stop_codons[order(stop_codons$n),]

## tx_biotype
## Vienna folding (minE)
rnafold <- read.csv(file="/Volumes/USELESS/DATA/fasta/danio_rerio/GRCz10/cdna/min_free_energies.txt", header = T, sep = "\t")
rnafold <- rnafold[order(rnafold$n),]
rownames(rnafold) <- NULL
rnafold <- rnafold[rnafold$n %in% matrix_of_everything$n,]
