############### cell1 ################

shape_cell1_et1 <- get(load(file = "/Volumes/USELESS/test2/normalized/cell1_invivo_et1.Rsave"))
shape_cell1_et2 <- get(load(file = "/Volumes/USELESS/test2/normalized/cell1_invivo_et2.Rsave"))

common1 <- shape_cell1_et1[names(shape_cell1_et1) %in% names(shape_cell1_et2)]
common2 <- shape_cell1_et2[names(shape_cell1_et2) %in% names(shape_cell1_et1)]
unique <- c(shape_cell1_et1[!names(shape_cell1_et1) %in% names(shape_cell1_et2)], shape_cell1_et2[!names(shape_cell1_et2) %in% names(shape_cell1_et1)])

u1 <- unlist(common1)
u2 <- unlist(common2)
u <- unlist(unique)

u1$TC.control <- u1$TC.control + u2$TC.control
u1$TC.treated <- u1$TC.treated + u2$TC.treated
u1$log2ratio <- log2(u1$TC.treated+1) - log2(u1$TC.control+1)

uu <- c(u1, u)
uu$log2ratio[uu$log2ratio < 0] <- 0

shape_cell1 <- split(uu, seqnames(uu))
shape_cell1 <- shape_cell1[sapply(shape_cell1, function(x){length(x) > 0})]
shape_cell1 <- shape_cell1[order(names(shape_cell1))]

save(shape_cell1, file = "/Volumes/USELESS/DATA/Shape-Seq/shape_cell1.Rsave")

############## cell24 ###############
shape_cell24 <- get(load(file = "/Volumes/USELESS/test/normalized/2-4cell.Rsave"))
u <- unlist(shape_cell24)
u$log2ratio[u$log2ratio < 0] <- 0
shape_cell24 <- split(u, seqnames(u))
shape_cell24 <- shape_cell24[sapply(shape_cell24, function(x){length(x) > 0})]
shape_cell24 <- shape_cell24[order(names(shape_cell24))]
save(shape_cell24, file = "/Volumes/USELESS/DATA/Shape-Seq/shape_cell24.Rsave")

############## cell256 ###############
shape_cell256 <- get(load(file = "/Volumes/USELESS/test/normalized/256cell.Rsave"))
u <- unlist(shape_cell256)
u$log2ratio[u$log2ratio < 0] <- 0
shape_cell256 <- split(u, seqnames(u))
shape_cell256 <- shape_cell256[sapply(shape_cell256, function(x){length(x) > 0})]
shape_cell256 <- shape_cell256[order(names(shape_cell256))]
save(shape_cell256, file = "/Volumes/USELESS/DATA/Shape-Seq/shape_cell256.Rsave")


############## cell1K ##############
shape_cell1K_et1 <- get(load(file = "/Volumes/USELESS/test2/normalized/cell1K_invivo_et1.Rsave"))
shape_cell1K_et2 <- get(load(file = "/Volumes/USELESS/test2/normalized/cell1K_invivo_et2.Rsave"))

common1 <- shape_cell1K_et1[names(shape_cell1K_et1) %in% names(shape_cell1K_et2)]
common2 <- shape_cell1K_et2[names(shape_cell1K_et2) %in% names(shape_cell1K_et1)]
unique <- c(shape_cell1K_et1[!names(shape_cell1K_et1) %in% names(shape_cell1K_et2)], shape_cell1K_et2[!names(shape_cell1K_et2) %in% names(shape_cell1K_et1)])

u1 <- unlist(common1)
u2 <- unlist(common2)
u <- unlist(unique)

u1$TC.control <- u1$TC.control + u2$TC.control
u1$TC.treated <- u1$TC.treated + u2$TC.treated
u1$log2ratio <- log2(u1$TC.treated+1) - log2(u1$TC.control+1)

uu <- c(u1, u)
uu$log2ratio[uu$log2ratio < 0] <- 0

shape_cell1K <- split(uu, seqnames(uu))
shape_cell1K <- shape_cell1K[sapply(shape_cell1K, function(x){length(x) > 0})]
shape_cell1K <- shape_cell1K[order(names(shape_cell1K))]

save(shape_cell1K, file = "/Volumes/USELESS/DATA/Shape-Seq/shape_cell1K.Rsave")


############ oblong ##############
shape_oblong <- get(load(file = "/Volumes/USELESS/test/normalized/shape_oblong.Rsave"))
u <- unlist(shape_oblong)
u$log2ratio[u$log2ratio < 0] <- 0
shape_oblong <- split(u, seqnames(u))
shape_oblong <- shape_oblong[sapply(shape_oblong, function(x){length(x) > 0})]
shape_oblong <- shape_oblong[order(names(shape_oblong))]
save(shape_oblong, file = "/Volumes/USELESS/DATA/Shape-Seq/shape_oblong.Rsave")

############ oblong CHX ##############
shape_oblong_CHX <- get(load(file = "/Volumes/USELESS/test/normalized/shape_oblong_CHX.Rsave"))
u <- unlist(shape_oblong_CHX)
u$log2ratio[u$log2ratio < 0] <- 0
shape_oblong_CHX <- split(u, seqnames(u))
shape_oblong_CHX <- shape_oblong_CHX[sapply(shape_oblong_CHX, function(x){length(x) > 0})]
shape_oblong_CHX <- shape_oblong_CHX[order(names(shape_oblong_CHX))]
save(shape_oblong_CHX, file = "/Volumes/USELESS/DATA/Shape-Seq/shape_oblong_CHX.Rsave")

######
shapes_cell24 <- get(load(file = "/Volumes/USELESS/test/normalized/cell2_4_invivo_NAI_N3.Rsave"))
shapes_cell1_invitro <- get(load(file = "/Volumes/USELESS/test/normalized/cell1_invitro_NAI_N3.Rsave"))
shapes_cell1_invitro_nonsel <- get(load(file = "/Volumes/USELESS/test/normalized/cell1_invitro_NAI_N3_nonsel.Rsave"))
rm(shapes_norm)

