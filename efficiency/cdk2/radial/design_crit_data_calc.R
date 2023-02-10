###################################################################
## This file computes the determinant and trace of a pregenerated design
###################################################################

###################################################################
## Load packages, set ggplot theme, and import data

# Load required packages
require( igraph)
require( tidyverse)

# VARS
k_numb = 15


# cdk2
ligand.list <- c("m0", "m1", "m2", "m3", "m4", "m5", "m6", "m7", "m8", "m9", "m10", "m11", "m12", "m13", "m14", "m15")
reference.ligand <- c("m3")
LIGAND_1 <- c("m3", "m3", "m3", "m3", "m3", "m3", "m3", "m3", "m3", "m3", "m3", "m3", "m3", "m3", "m3")
LIGAND_2 <- c("m0", "m1", "m2", "m4", "m5", "m6", "m7", "m8", "m9", "m10", "m11", "m12", "m13", "m14", "m15")
SIM_WEIGHT <- c(0.01657268, 0.02024191, 0.01657268, 0.40656966, 0.10025884, 0.60653066, 0.14956862, 0.14956862, 0.14956862, 0.33287108, 0.90483742, 0.74081822, 1.0, 0.44932896, 0.74081822)

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


