# This script it adapted by M Pitman from  https://doi.org/10.1002/jcc.26095
# The original publication is:
# Yang Q, Burchett W, Steeno GS, Liu S, Yang M, Mobley DL, Hou X.
# J. Comput. Chem. 2020 jan; 41(3):247–257.

###################################################################
## This file simulates data from FEP experiments, estimates the dG values with
## different designs, and summarizes various summary statistics about these estimates
###################################################################

###################################################################
## Load packages

# Load required packages
require(igraph)
require(ggplot2)
require(scales)
require(reshape2)
require(dplyr)
source( 'exchange_algorithm.R')

# Set ggplot2 theme
theme_set(theme_gray(base_size = 18))
theme_update( plot.title = element_text(hjust = 0.5))

# Define concordance correlation coefficient
pc <- function(x, y){
  return((2*var(x, y))/(var(x) + var(y) + (mean(x)-mean(y))^2))
}

######################################################
## Define simulation parameters

# Define the number of replicates for each simulation scenario
M <- 5

# Define simulation scenarios
results <- expand.grid(
  N = 15,
  N.pair = c( 33, 35),
  # 35, 40, 45
  N.known = c( 1),
  optimality = c( 'D', 'A', 'random'),
  ddG.MSE = 1,
  dG.weight = 2,
  rep = 1:M,
  stringsAsFactors = FALSE
) %>% filter( N.pair >= N) %>% filter( N.pair < choose( N, 2))
results$Label <- paste0( results$N, ' ligands, ', results$N.pair, ' pairs (', results$optimality, ')') # Include all parameters that vary in the label

# Function to calculate number of cycles.
FindCycles = function(g) {
    Cycles = NULL
    for(v1 in V(g)) {
        if(degree(g, v1, mode="in") == 0) { next }
        GoodNeighbors = neighbors(g, v1, mode="out")
        GoodNeighbors = GoodNeighbors[GoodNeighbors > v1]
        for(v2 in GoodNeighbors) {
            TempCyc = lapply(all_simple_paths(g, v2,v1, mode="out"), function(p) c(v1,p))
            TempCyc = TempCyc[which(sapply(TempCyc, length) > 3)]
          TempCyc = TempCyc[sapply(TempCyc, min) == sapply(TempCyc, `[`, 1)]
          Cycles  = c(Cycles, TempCyc)
        }
    }
    Cycles
}

######################################################
## Create the designs

design.dat <- results %>% select( N, N.pair, N.known, optimality, dG.weight, Label) %>%
  unique() %>% group_by( Label) %>% do( mat = {
  
  # Generate full design
  A.full <- matrix( 0, ncol=.$N, nrow = choose( .$N, 2))
  for( i in 1:nrow( A.full)){
    A.full[ i, t( combn( 1:(.$N), 2))[i,]] <- c( 1, -1)
  }
  
  # Define the A and W matrices for a full design
  known.ind <- 1:(.$N.known)
  A.use <- rbind( A.full, diag( ncol( A.full))[ known.ind,])
  W.use <- diag( c( rep( 1, nrow( A.full)), rep( .$dG.weight, length( known.ind))))
  
  # Compute the optimal design
  starting.rows <- Federov.exchange( A.use, W.use, c( 1:(.$N.pair), nrow( A.use)), 'random')$rows
  current.design <- Federov.exchange( A.use, W.use, starting.rows, .$optimality)$A
  
  current.design
})

######################################################
## Run the simulations

# Run analysis for each different scenario and store results
for( current.rep in 1:nrow( results)){
  
  message( paste0( 'Compute Results ', current.rep, ' of ', nrow( results), '.'))
  
  # Generate true dG values
  sim.dat <- data.frame( Ligand = paste0( 'Ligand ', 1:results$N[current.rep]), dG.true = rnorm( results$N[current.rep], 0, 1), stringsAsFactors=FALSE)
  # Generate experimental dG values (the precision relative to the ddG variability defined by the weight parameter)
  sim.dat$dG.exp <- rnorm( nrow( sim.dat), sim.dat$dG.true, sqrt( results$ddG.MSE[current.rep]/results$dG.weight[current.rep]))
  
  # Generate full design
  A.full <- matrix( 0, ncol = results$N[current.rep], nrow = choose( results$N[current.rep], 2))
  for( i in 1:nrow( A.full)){
    A.full[ i, t( combn( 1:results$N[current.rep], 2))[i,]] <- c( 1, -1)
  }
  
  # Select optimal design
  if( results$optimality[current.rep] == 'random'){ # Generate a random design when random is specified
    known.ind <- 1:results$N.known[current.rep]
    A.use <- rbind( A.full, diag( ncol( A.full))[ known.ind,])
    W.use <- diag( c( rep( 1, nrow( A.full)), rep( results$dG.weight[current.rep], length( known.ind))))
    A <- Federov.exchange( A.use, W.use, c( 1:results$N.pair[current.rep], nrow( A.use)), 'random')$A
  }else{
    A <- design.dat[ which( design.dat$Label == results$Label[current.rep]), 'mat'][[1]][[1]]
  }
  A <- A[which(rowSums(abs(A)) > 1),]
  known.ind <- 1:results$N.known[current.rep]
  #print("A")
  #print(A)

  # Compute ddG values
  sim.dat.pairs <- data.frame( Ligand1 = sim.dat$Ligand[ apply( A, 1, function(x){ which( x == 1)})], Ligand2 = sim.dat$Ligand[ apply( A, 1, function(x){ which( x == -1)})], ddG.true = as.numeric( A %*% sim.dat$dG.true))
  #print("sim.dat.pairs")
  #print(sim.dat.pairs)
  
  # Generate a graph from simulation data
  c <- graph_from_data_frame(sim.dat.pairs, directed = FALSE)
  #print(c)
  #print("cycles")
  # Calculate the number of cycles. This is np complete, comment out for large designs.
  cycle.number <- length(FindCycles(c))
  #print(cycle.number)
  
  sim.dat.pairs$ddG.pred <- rnorm( nrow( sim.dat.pairs), sim.dat.pairs$ddG.true, sqrt( results$ddG.MSE[current.rep]))
  
  # Modify design matrix to include known experimental dG information
  tmp.A <- rbind(A, diag( ncol(A))[ known.ind,])
  # Compute weights for weighted least squares based on relative variability between experimental dGs and computed ddGs
  w <- c(rep(1, nrow( sim.dat.pairs)), rep( results$dG.weight[current.rep], length( known.ind)))
  # Estimate dG values
  dG.pred <- as.data.frame( coef( summary( lm( c( sim.dat.pairs$ddG.pred, sim.dat$dG.exp[ known.ind]) ~ tmp.A - 1, weights = w)))[, c( 'Estimate', 'Std. Error')])
  rownames( dG.pred) <- sim.dat$Ligand
  
  # Compute analytic variability
  inv.inf <- solve( t(tmp.A) %*% diag(w) %*% tmp.A)
  # Compute variability for all pairs
  all.pairs.var <- NULL
  for( i in 1:(nrow( inv.inf) - 1)){
    for( j in (i+1):nrow( inv.inf)){
      all.pairs.var <- c( all.pairs.var, results$ddG.MSE[current.rep]*(inv.inf[i,i] + inv.inf[j,j] - 2*inv.inf[i,j]))
    }
  }
  # Compute variability for FEP pairs
  fep.pairs.var <- NULL
  pair.ind <- cbind( apply(A, 1, function(x){which(x==1)}), apply(A, 1, function(x){which(x==-1)}))
  for( i in 1:nrow(pair.ind)){
    fep.pairs.var <- c( fep.pairs.var, results$ddG.MSE[current.rep]*(inv.inf[pair.ind[i,1],pair.ind[i,1]] + inv.inf[pair.ind[i,2],pair.ind[i,2]] - 2*inv.inf[pair.ind[i,1],pair.ind[i,2]]))
  }
  # Compute variability for all dGs
  dG.var <- results$ddG.MSE[current.rep]*diag(inv.inf)
  
  # Store results
  results[ current.rep, 'DOptim'] <- log( det( t( tmp.A) %*% diag(w) %*% tmp.A))
  results[ current.rep, 'MSE.dG'] <- mean((sim.dat$dG.true - dG.pred$Estimate)^2)
  results[ current.rep, 'rank.corr'] <- cor( sim.dat$dG.true, dG.pred$Estimate, method = 'spearman')
  results[ current.rep, 'correct.min'] <- as.numeric( which.min( sim.dat$dG.true) == which.min( dG.pred$Estimate))
  results[ current.rep, 'MSE.ddG.fep'] <- mean((sim.dat.pairs$ddG.true - sim.dat.pairs$ddG.pred)^2)
  results[ current.rep, 'MSE.ddG.back'] <- mean((sim.dat.pairs$ddG.true - as.numeric( A %*% dG.pred$Estimate))^2)
  results[ current.rep, 'MSE.ddG.back.all'] <- mean((as.numeric( A.full %*% sim.dat$dG.true) - as.numeric( A.full %*% dG.pred$Estimate))^2)
  results[ current.rep, 'MSE.dG.theory'] <- mean( dG.var)
  results[ current.rep, 'MSE.ddG.back.theory'] <- mean( fep.pairs.var)
  results[ current.rep, 'MSE.ddG.back.all.theory'] <- mean( all.pairs.var)
  results[ current.rep, 'N.cycles'] <- cycle.number
}
#print(results)

######################################################
## Summarize the results (these plots and tables may need to change depending on which simulation parameters are varying)

# Format labels to look nicer
results$N <- paste0( results$N, ' ligands')
results$N.pair <- paste0( results$N.pair, ' pairs')
results$optimality[ which( results$optimality != 'random')] <- paste0( results$optimality[ which( results$optimality != 'random')], '-Optimal')
results$optimality[ which( results$optimality == 'random')] <- 'Random'

# Summarize results in a table
summary.dat <- as.data.frame(
  results %>%
    group_by( N, N.pair, optimality, Label) %>%
    summarize(
      MSE.dG.theory = median( MSE.dG.theory),
      N.cycles = mean(N.cycles),
      MSE.ddG.back.theory = median( MSE.ddG.back.theory),
      MSE.ddG.back.all.theory = median( MSE.ddG.back.all.theory),
      rank.corr.mean = mean( rank.corr),
      rank.corr.median = median( rank.corr),
      correct.min = paste0( round( 100*mean( correct.min), 1), '%')
    )
)
print( summary.dat[4:6])

# Plot distribution of rank correlation between delta G estimates and true values
ggplot( results, aes( x = 1, y = rank.corr)) +
  geom_violin( width = 0.7, aes( fill = optimality)) + facet_grid( N.pair ~ optimality) + theme( legend.position = 'top', axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  geom_point( data = summary.dat, aes( y = rank.corr.mean), x = 1) + guides( fill = "none") +
  scale_y_continuous( breaks = seq( -0.2, 1, by = 0.2)) +
  labs( x = '', y = expression('Rank Correlation between '*Delta*'G Estimates and True '*Delta*"Gs")) + scale_x_continuous( limits = c(0, 2), minor_breaks = NULL, breaks = 1)

# Plot distribution of theoretical MSEs for delta G
p1 <- ggplot( subset( results, optimality == 'Random'), aes( x = MSE.dG.theory)) +
  geom_density( aes( fill = 'Random Design'), alpha = 0.6, size = 1) + facet_wrap( ~ N.pair, ncol = 1) +
  theme( legend.position = 'top', legend.box = "horizontal") +
  geom_vline( size = 1, aes( xintercept = MSE.dG.theory, color = 'D-Optimal Design'), data = subset( results, optimality == 'D-Optimal')) + 
  geom_vline( size = 1, aes( xintercept = MSE.dG.theory, color = 'A-Optimal Design'), data = subset( results, optimality == 'A-Optimal')) + 
  labs( x = expression( 'Theoretical Mean MSE of '*Delta*'G Estimates'), y = '', fill = '', color = "") +
  scale_fill_manual( values = hue_pal()(3)[3]) +
  scale_color_manual( values = hue_pal()(3)[1:2]) +
  theme_bw( base_size = 18) + theme( legend.position = 'top', legend.box = 'horizontal') +
  scale_x_continuous( trans = log2_trans(), breaks = c( 0.5, 1, 2, 4, 8), minor_breaks = 1:16)

p2 <- ggplot( results, aes( x = MSE.dG)) +
  geom_density( aes( fill = optimality), alpha = 0.6, size = 1) + facet_wrap( ~ N.pair, ncol = 1) +
  theme( legend.position = 'top', legend.box = "horizontal") +
  labs( x = expression( 'Distribution of MSEs of '*Delta*'G Estimates from Simulation'), y = '', fill = '', color = "") +
  scale_fill_manual( values = hue_pal()(3)[1:3]) +
  theme_bw( base_size = 18) +
  scale_x_continuous( trans = log2_trans(), breaks = c( 0.25, 0.5, 1, 2, 4, 8, 16, 32), minor_breaks = 1:64, labels = as.numeric)

ggpubr::ggarrange( p1, p2, common.legend = TRUE, legend = 'top')

# Analytical MSE by pairs plot
plot.dat <- results %>%
  filter( optimality != 'Random') %>%
  mutate( N.pair = as.numeric( gsub( ' pairs', '', N.pair)))
ggplot( plot.dat, aes( x = N.pair, y = MSE.dG.theory, color = optimality)) +
  geom_point( size = 2) + geom_line( size = 1) + facet_wrap( ~ N, scales = 'free_x') +
  labs( x = 'Number of Pairs', y = expression( 'Theoretical Average MSE of '*Delta*'G'), color = '') +
  theme_bw( base_size = 18) +
  theme( legend.position = 'top', legend.box = 'vertical')

