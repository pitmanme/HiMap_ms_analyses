###################################################################
## This file computes the determinant and trace of a pregenerated design
###################################################################

###################################################################
## Load packages, set ggplot theme, and import data

# Load required packages
require( igraph)
require( tidyverse)

# VARS
k_numb = 39


# pfkfb3 radial
ligand.list <- c("m0", "m1", "m2", "m3", "m4", "m5", "m6", "m7", "m8", "m9", "m10", "m11", "m12", "m13", "m14", "m15", "m16", "m17", "m18", "m19", "m20", "m21", "m22", "m23", "m24", "m25", "m26", "m27", "m28", "m29", "m30", "m31", "m32", "m33", "m34", "m35", "m36", "m37", "m38", "m39")
reference.ligand <- c("m10")
LIGAND_1 <- c("m10", "m10", "m10", "m10", "m10", "m10", "m10", "m10", "m10", "m10", "m10", "m10", "m10", "m10", "m10", "m10", "m10", "m10", "m10", "m10", "m10", "m10", "m10", "m10", "m10", "m10", "m10", "m10", "m10", "m10", "m10", "m10", "m10", "m10", "m10", "m10", "m10", "m10", "m10")
LIGAND_2 <- c("m0", "m1", "m2", "m3", "m4", "m5", "m6", "m7", "m8", "m9", "m11", "m12", "m13", "m14", "m15", "m16", "m17", "m18", "m19", "m20", "m21", "m22", "m23", "m24", "m25", "m26", "m27", "m28", "m29", "m30", "m31", "m32", "m33", "m34", "m35", "m36", "m37", "m38", "m39")
SIM_WEIGHT <- c(0.44932896, 0.33287108, 0.4965853, 0.40656966, 0.36787944, 1.0, 0.40656966, 0.4965853, 0.44932896, 0.40656966, 0.90483742, 0.16529889, 0.082085, 0.40656966, 0.18268352, 0.33287108, 0.40656966, 0.12245643, 0.06720551, 0.10025884, 0.10025884, 0.12245643, 0.90483742, 0.74081822, 0.74081822, 0.90483742, 0.90483742, 0.74081822, 0.30119421, 0.0022313, 0.30119421, 0.74081822, 0.74081822, 0.74081822, 0.81873075, 0.67032005, 0.60653066, 0.90483742, 0.90483742)

#WEIGHT <- c((SIM_WEIGHT^2 - min(SIM_WEIGHT)^2)/(max(SIM_WEIGHT)^2 - min(SIM_WEIGHT)^2))
#print(WEIGHT)

analysis.dat<- data.frame(LIGAND_1, LIGAND_2, SIM_WEIGHT)

#analysis.dat <- dataframe %>%
  #as.matrix() %>%
  #reshape2::melt( as.is = TRUE) %>%
  #as_tibble() %>%
  #setNames( c( 'LIGAND_1', 'LIGAND_2', 'SIM'))


print( "analysis.dat")
print( analysis.dat)

# Compute the optimality criteria for each design.
crit.dat <- analysis.dat %>%
  tidyr::crossing( as_tibble( t( setNames( rep( 0, length( ligand.list)), ligand.list)))) %>%
  rowwise() %>%
  do( mutate_if( as_tibble(.), grepl( .$LIGAND_1, names(.)), function(x){-1})) %>%
  do( mutate_if( as_tibble(.), grepl( .$LIGAND_2, names(.)), function(x){1})) %>%
  ungroup() %>%
  do({
    design.mat <- {.} %>%
      select_at( ligand.list) %>%
      as.matrix() %>%
      rbind( diag( length( ligand.list))[ which( ligand.list == reference.ligand),])
    
    sim.weight <- c( .$SIM_WEIGHT, 2)
    weighted.inner =  t( design.mat) %*% diag( sim.weight) %*% design.mat
    
    print(" design.mat")
    print( design.mat)
    print(" sim.weight")
    print( sim.weight)

    print("weighted.inner")
    print(weighted.inner)
    print("Determinant:")
    print(det( solve( t( design.mat) %*% diag( sim.weight) %*% design.mat))^(1/length( ligand.list)))
    print("Trace:")
    print(sum( diag( solve( t( design.mat) %*% diag( sim.weight) %*% design.mat))))
    # Removed toggle for ap vs lomap weight.
    data.frame(
      A.ap = sum( diag( solve( t( design.mat) %*% diag( sim.weight) %*% design.mat))),
      D.ap = det( solve( t( design.mat) %*% diag( sim.weight) %*% design.mat))^(1/length( ligand.list)))
  })


