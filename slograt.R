library(RNAprobR)

# save...
load(file="/Home/ii/katchyz/OUT/SHAPES_out/control_comp_2-4cell.Rsave")
load(file="/Home/ii/katchyz/OUT/SHAPES_out/treated_comp_2-4cell.Rsave")

# normalization
slograt_unsmoothed24 <- slograt(control_GR = control_comp, treated_GR = treated_comp, window_size = 1)

# save
save(slograt_unsmoothed24, file="/Home/ii/katchyz/OUT/SHAPES_out/slograt_unsmoothed24.Rsave")



