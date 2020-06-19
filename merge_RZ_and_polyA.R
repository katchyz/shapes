### merge RZ and polyA libs
### start/stop on raw
### start/stop smoothed


libs_path <- "/Volumes/USELESS/DATA/Shape-Seq/dTCR"

## dTCR

####### VIVO
paired_polyA <- get(load(file = c(file.path(libs_path, "PAIR_GC_1Kcell_polyA.Rsave"))))
paired_RZ <- get(load(file = c(file.path(libs_path, "PAIR_GC_1Kcell_RZ.Rsave"))))


rm(GR)

#names(paired_RZ) <- sapply(names(paired_RZ), function(x){substr(x,1,18)})

## MERGE
polyA <- split(paired_polyA, names(paired_polyA))
RZ <- split(paired_RZ, names(paired_RZ))

polyA <- polyA[names(polyA) %in% names(RZ)]
RZ <- RZ[names(RZ) %in% names(polyA)]

polyA <- unlist(polyA)
RZ <- unlist(RZ)
vivo <- polyA
mcols(vivo) <- NULL
vivo$TC.control <- polyA$TC.control + RZ$TC.control
vivo$TC.treated <- polyA$TC.treated + RZ$TC.treated
vivo$cov.control <- polyA$cov.control + RZ$cov.control
vivo$cov.treated <- polyA$cov.treated + RZ$cov.treated
## calculate TCR
vivo$TCR.control <- vivo$TC.control / vivo$cov.control
vivo$TCR.treated <- vivo$TC.treated / vivo$cov.treated
vivo$TCR.control[is.na(vivo$TCR.control)] <- 0
vivo$TCR.treated[is.na(vivo$TCR.treated)] <- 0

dtcr <- (vivo$TCR.treated - vivo$TCR.control) / (1 - vivo$TCR.control)
dtcr[dtcr < 0] <- 0
dtcr[is.na(dtcr)] <- 0
vivo$dTCR <- dtcr

save(vivo, file=c(file.path(libs_path, "1Kcell_vivo_dTCR.Rsave")))


####### VITRO
paired_polyA_vitro <- get(load(file = c(file.path(libs_path, "PAIR_GC_1Kcell_polyA_vitro.Rsave"))))
paired_RZ_vitro <- get(load(file = c(file.path(libs_path, "PAIR_GC_1Kcell_RZ_vitro.Rsave"))))

rm(GR)

## MERGE
polyA_vitro <- split(paired_polyA_vitro, names(paired_polyA_vitro))
RZ_vitro <- split(paired_RZ_vitro, names(paired_RZ_vitro))

polyA_vitro <- polyA_vitro[names(polyA_vitro) %in% names(RZ_vitro)]
RZ_vitro <- RZ_vitro[names(RZ_vitro) %in% names(polyA_vitro)]

polyA_vitro <- unlist(polyA_vitro)
RZ_vitro <- unlist(RZ_vitro)
vitro <- polyA_vitro
mcols(vitro) <- NULL
vitro$TC.control <- polyA_vitro$TC.control + RZ_vitro$TC.control
vitro$TC.treated <- polyA_vitro$TC.treated + RZ_vitro$TC.treated
vitro$cov.control <- polyA_vitro$cov.control + RZ_vitro$cov.control
vitro$cov.treated <- polyA_vitro$cov.treated + RZ_vitro$cov.treated
## calculate TCR
vitro$TCR.control <- vitro$TC.control / vitro$cov.control
vitro$TCR.treated <- vitro$TC.treated / vitro$cov.treated
vitro$TCR.control[is.na(vitro$TCR.control)] <- 0
vitro$TCR.treated[is.na(vitro$TCR.treated)] <- 0

dtcr <- (vitro$TCR.treated - vitro$TCR.control) / (1 - vitro$TCR.control)
dtcr[dtcr < 0] <- 0
dtcr[is.na(dtcr)] <- 0
vitro$dTCR <- dtcr

save(vitro, file=c(file.path(libs_path, "1Kcell_vitro_dTCR.Rsave")))

