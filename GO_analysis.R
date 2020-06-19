### GO analysis

go_cc <- read.csv(file="/Volumes/USELESS/DATA/genomes/GO/go_cc.tsv", sep="\t")
go_cc <- data.frame(n = go_cc$ID, GO_CC = go_cc$GO.Name)
go_cc <- unique(go_cc)
go_cc <- go_cc[!duplicated(go_cc$n),]

df <- merge(df, go_cc, by="n", all=TRUE)


hybr_cloud <- df[df$cloud_oblong_CHX == "cloud" & !is.na(df$cloud_oblong_CHX),]$n
hybr_main <- df[(df$shape_oblong_CHX_log2ratio_internal > 0) & (df$rna_3_5h > 0),]$n
hybr_main <- hybr_main[!(hybr_main %in% hybr_cloud)]
hybr_main <- hybr_main[complete.cases(hybr_main)]


cld <- df[df$n %in% hybr_cloud,]$GO_CC
cld <- cld[complete.cases(cld)]
rst <- df[df$n %in% hybr_main,]$GO_CC
rst <- rst[complete.cases(rst)]

head(summary(cld)/length(cld))
head(summary(rst)/length(rst))


#########
## cell1
> head(summary(cld)/length(cld))
membrane            nucleus cellular_component          cytoplasm      intracellular      mitochondrion 
0.20274390         0.17225610         0.12195122         0.07926829         0.07926829         0.01981707 
> head(summary(rst)/length(rst))
membrane            nucleus cellular_component          cytoplasm      intracellular    plasma membrane 
0.20540815         0.15767424         0.13609242         0.07744065         0.06918878         0.03199188
## cell1 invitro
> head(summary(cld)/length(cld))
nucleus           membrane cellular_component      intracellular          cytoplasm    plasma membrane 
0.24600639         0.15654952         0.13418530         0.08945687         0.07348243         0.02555911 
> head(summary(rst)/length(rst))
membrane            nucleus cellular_component          cytoplasm      intracellular    plasma membrane 
0.19810896         0.15878733         0.14092751         0.07909350         0.06618640         0.02761519
## cell24
> head(summary(cld)/length(cld))
nucleus           membrane cellular_component          cytoplasm      intracellular    plasma membrane 
0.20282120         0.18604651         0.11170416         0.07815478         0.07548608         0.02706824 
> head(summary(rst)/length(rst))
membrane            nucleus cellular_component          cytoplasm      intracellular    plasma membrane 
0.23218413         0.14474191         0.12596302         0.07155239         0.06490755         0.04140986 
## cell256
> head(summary(cld)/length(cld))
membrane            nucleus cellular_component          cytoplasm      intracellular    plasma membrane 
0.20526794         0.18952467         0.10868907         0.07992734         0.07538601         0.02633969 
> head(summary(rst)/length(rst))
membrane            nucleus cellular_component          cytoplasm      intracellular    plasma membrane 
0.23873874         0.14000735         0.12419562         0.07096893         0.06499356         0.04780290
## cell1K
> head(summary(cld)/length(cld))
nucleus           membrane cellular_component          cytoplasm      intracellular    plasma membrane 
0.18523878         0.18379161         0.12301013         0.08393632         0.06222865         0.02894356 
> head(summary(rst)/length(rst))
membrane            nucleus cellular_component          cytoplasm      intracellular    plasma membrane 
0.20322403         0.15954808         0.14025902         0.07784514         0.06764949         0.03196473 
## oblong
> head(summary(cld)/length(cld))
nucleus           membrane cellular_component          cytoplasm      intracellular      mitochondrion 
0.20222991         0.17377932         0.11534025         0.07958478         0.07920031         0.02268358 
> head(summary(rst)/length(rst))
membrane            nucleus cellular_component          cytoplasm      intracellular    plasma membrane 
0.22666667         0.14729167         0.13541667         0.07531250         0.06531250         0.03979167
## oblong_CHX
> head(summary(cld)/length(cld))
membrane            nucleus cellular_component      intracellular          cytoplasm    plasma membrane 
0.20139635         0.18152524         0.11922664         0.07769424         0.07482993         0.02649481 
> head(summary(rst)/length(rst))
membrane            nucleus cellular_component          cytoplasm      intracellular    plasma membrane 
0.23512158         0.14249219         0.12905668         0.07351689         0.06566373         0.04144195


########### for GOrilla
hybr_cloud <- df[df$cloud_oblong == "cloud" & !is.na(df$cloud_oblong),]$gene
hybr_main <- df[(df$shape_oblong_log2ratio_internal > 0) & (df$rna_3_5h > 0),]$gene
hybr_main <- hybr_main[!(hybr_main %in% hybr_cloud)]
hybr_main <- hybr_main[complete.cases(hybr_main)]
cloud <- as.character(unique(hybr_cloud))
main <- as.character(unique(hybr_main))
write(cloud, file="/Users/kasia/Downloads/cloud_oblong.txt", sep="\n")
write(main, file="/Users/kasia/Downloads/main_oblong.txt", sep="\n")

## http://cbl-gorilla.cs.technion.ac.il/GOrilla/rd5olrq4/GOResults.html