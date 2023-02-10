###################################################################
## This file computes the determinant and trace of a pregenerated design
###################################################################

###################################################################
## Load packages, set ggplot theme, and import data

# Load required packages
require( igraph)
require( tidyverse)

# VARS
k_numb = 22



#ligand.list <- c("m0", "m1", "m2", "m3", "m4", "m5", "m6", "m7", "m8", "m9")
#reference.ligand <- c("m5")
#LIGAND_1 <- c("m0", "m0", "m1", "m1", "m2", "m3", "m4", "m4", "m5", "m6", "m6", "m7")
#LIGAND_2 <- c("m4", "m2", "m8", "m9", "m3", "m4", "m8", "m5", "m8", "m8", "m9", "m8")
#SIM_WEIGHT <- c(0.8187307530779818, 0.9048374180359595, 0.6703200460356393, 0.6703200460356393, 0.8187307530779818, 0.7408182206817179, 0.6703200460356393, 0.951229424500714, 0.6376281516217733, 0.8187307530779818, 0.951229424500714, 0.40656965974059905)


# ptp1b
ligand.list <- c("m0", "m1", "m2", "m3", "m4", "m5", "m6", "m7", "m8", "m9", "m10", "m11", "m12", "m13", "m14", "m15", "m16", "m17", "m18", "m19", "m20", "m21", "m22")
reference.ligand <- c("m4")
LIGAND_1 <- c("m4", "m4", "m4", "m4", "m4", "m4", "m4", "m4", "m4", "m4", "m4", "m4", "m4", "m4", "m4", "m4", "m4", "m4", "m4", "m4", "m4", "m4")
LIGAND_2 <- c("m0", "m1", "m2", "m3", "m5", "m6", "m7", "m8", "m9", "m10", "m11", "m12", "m13", "m14", "m15", "m16", "m17", "m18", "m19", "m20", "m21", "m22")
SIM_WEIGHT <- c(0.01887343, 0.25410696, 0.31036694, 0.13945686, 0.76337949, 0.46301307, 0.20804518, 0.31036694, 0.46301307, 0.13945686, 0.22992549, 0.25410696, 0.37908304, 0.25410696, 0.31036694, 0.17033299, 0.11417762, 0.062662, 0.02305206, 0.01545226, 0.01707739, 0.01831564)

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


