###################################################################
## This file computes the determinant and trace of a pregenerated design
###################################################################

###################################################################
## Load packages, set ggplot theme, and import data

# Load required packages
require( igraph)
require( tidyverse)

# VARS
k_numb = 56


# pfkfb3
ligand.list <- c("m0", "m1", "m2", "m3", "m4", "m5", "m6", "m7", "m8", "m9", "m10", "m11", "m12", "m13", "m14", "m15", "m16", "m17", "m18", "m19", "m20", "m21", "m22", "m23", "m24", "m25", "m26", "m27", "m28", "m29", "m30", "m31", "m32", "m33", "m34", "m35", "m36", "m37", "m38", "m39")
reference.ligand <- c("m10")
LIGAND_1 <- c("m0", "m0", "m1", "m1", "m2", "m2", "m2", "m2", "m2", "m2", "m3", "m4", "m4", "m5", "m5", "m5", "m5", "m6", "m9", "m10", "m10", "m10", "m10", "m10", "m10", "m13", "m13", "m13", "m14", "m15", "m17", "m17", "m17", "m18", "m18", "m19", "m19", "m20", "m23", "m23", "m23", "m24", "m25", "m26", "m29", "m29", "m31", "m32", "m32", "m33", "m33", "m34", "m35", "m35", "m35", "m36", "m38")
LIGAND_2 <- c("m2", "m3", "m15", "m17", "m3", "m6", "m7", "m8", "m9", "m10", "m8", "m9", "m12", "m9", "m11", "m23", "m27", "m7", "m12", "m11", "m23", "m26", "m27", "m35", "m39", "m14", "m16", "m30", "m27", "m16", "m21", "m22", "m27", "m21", "m22", "m20", "m21", "m21", "m24", "m25", "m28", "m25", "m28", "m27", "m31", "m39", "m39", "m33", "m35", "m34", "m35", "m35", "m36", "m37", "m38", "m37", "m39")
WEIGHT <- c(0.951229424500714, 0.8187307530779818, 0.7408182206817179, 0.8913661439068313, 0.7788007830714049, 0.8187307530779818, 0.8187307530779818, 0.8187307530779818, 0.7408182206817179, 0.6703200460356393, 0.951229424500714, 0.951229424500714, 0.522045776761016, 0.6703200460356392, 0.8607079764250578, 0.8607079764250578, 0.8607079764250578, 0.9048374180359595, 0.5488116360940264, 0.9048374180359595, 0.9048374180359595, 0.9048374180359595, 0.9048374180359595, 0.8607079764250578, 0.9048374180359595, 0.8187307530779818, 0.8607079764250578, 0.49658530379140947, 0.6376281516217733, 0.7046880897187134, 0.7408182206817179, 0.8187307530779818, 0.6376281516217733, 0.8913661439068313, 0.951229424500714, 0.7788007830714049, 0.7408182206817178, 0.951229424500714, 0.9048374180359595, 0.9048374180359595, 0.8607079764250578, 0.8187307530779818, 0.9048374180359595, 0.951229424500714, 0.8187307530779818, 0.5769498103804865, 0.5769498103804865, 0.951229424500714, 0.9048374180359595, 0.9048374180359595, 0.9048374180359595, 0.9048374180359595, 0.9048374180359595, 0.8607079764250578, 0.8607079764250578, 0.951229424500714, 0.951229424500714)

# Try to augment the SIM_WEIGHTs here.
maxim <- max(WEIGHT)
minim <- min(WEIGHT)

SIM_WEIGHT <- c((WEIGHT^2 - minim^2+0.00000001)/(maxim^2 - minim^2))
print(SIM_WEIGHT)

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
    
    # Compute the diagonal values of the hat matrix, the leverages
    print("H diags")
    print(diag(((design.mat %*% solve(t(design.mat) %*% diag( sim.weight) %*% design.mat) %*% t(design.mat) %*% diag( sim.weight)))))
    
    
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


