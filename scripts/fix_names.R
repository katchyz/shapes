### FIX NAMES

##
libs_path <- "/Volumes/USELESS/DATA/Shape-Seq/dTCR"

libs <- c("PAIR_GC_cell24.Rsave", "PAIR_GC_cell256.Rsave", "1Kcell_vivo_dTCR.Rsave",
          "1Kcell_vitro_dTCR.Rsave", "PAIR_GC_1Kcell_RZ.Rsave",
          "PAIR_GC_1Kcell_RZ_vitro.Rsave", "PAIR_GC_1Kcell_polyA.Rsave",
          "PAIR_GC_1Kcell_polyA_vitro.Rsave")


for (i in 1:length(libs)) {
  GR <- get(load(file = c(file.path(libs_path, libs[i]))))
  
  GR <- split(GR, names(GR))
  
  names_GR <- substring(names(GR), 1, 18)
  names(GR) <- names_GR
  
  names_GR <- rep(names_GR, width(GR@partitioning))
  
  GR <- unlist(GR)
  names(GR) <- names_GR
  
  save(GR, file = c(file.path(libs_path, libs[i])))
}


###### single
libs_path <- "/Volumes/USELESS/DATA/Shape-Seq/"
libs  <- c("shapes_cell1_invitro_nonsel.Rsave")

libs <- c("shape_cell1.Rsave", "shape_cell1K.Rsave", "shape_oblong.Rsave", "shape_oblong_CHX.Rsave")

for (i in 1:length(libs)) {
  GR <- get(load(file = c(file.path(libs_path, libs[i]))))
  
  names_GR <- rep(names(GR), width(GR@partitioning))
  
  GR <- unlist(GR)
  names(GR) <- names_GR
  
  save(GR, file = c(file.path(libs_path, libs[i])))
}
