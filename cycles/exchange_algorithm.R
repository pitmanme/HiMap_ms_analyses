# This script it adapted by M Pitman from  https://doi.org/10.1002/jcc.26095
# The original publication is:
# Yang Q, Burchett W, Steeno GS, Liu S, Yang M, Mobley DL, Hou X.
# J. Comput. Chem. 2020 jan; 41(3):247–257.

###################################################################
# This function uses the Federov exchange algorithm to find the design which maximizes an optimality criterion
# The inputs are the full design matrix (A), the weight matrix (W), the starting rows, and the criterion to maximize
###################################################################

Federov.exchange <- function( A.use, W.use, starting.rows, criterion, maxiter = 1e6, messages = TRUE){
  
  # Define the functions which compute the optimality criteria (return -Inf in cases where the information matrix is not invertible)
  if( criterion == 'D'){ # D-optimal (minimizes determinant of the information matrix)
    crit.func <- function( A, W){
      return( det( t(A) %*% W %*% A))
    }
  }else if( criterion == 'A'){ # A-optimal (minimizes the trace of the inverse information matrix, equiavalent to minimizing the average variability of the estimates)
    crit.func <- function( A, W){
      inf.mat <- t(A) %*% W %*% A
      if( class( try( solve( inf.mat), silent = TRUE))[1] == "try-error"){
        return( -Inf)
      }else{
        return( -sum( diag( solve( inf.mat))))
      }
    }
  }else if( criterion == 'P'){ # Minimizes the average variability of all possible pairwise differences
    crit.func <- function( A, W){
      inf.mat <- t(A) %*% W %*% A
      if( class( try( solve( inf.mat), silent = TRUE))[1] == "try-error"){
        return( -Inf)
      }else{
        inv.inf <- solve( inf.mat)
        all.pairs.var <- NULL
        for( i in 1:(nrow( inv.inf) - 1)){
          for( j in (i + 1):nrow( inv.inf)){
            all.pairs.var <- c( all.pairs.var, ( inv.inf[i,i] + inv.inf[j,j] - 2*inv.inf[i,j]))
          }
        }
        return( -sum( all.pairs.var))
      }
    }
  }else if( criterion == 'mA'){ # Minimizes the maximum variability of the estimates
    crit.func <- function( A, W){
      inf.mat <- t(A) %*% W %*% A
      if( class( try( solve( inf.mat), silent = TRUE))[1] == "try-error"){
        return( -Inf)
      }else{
        return( -max( diag( solve( inf.mat))))
      }
    }
  }else if( criterion == 'mP'){ # Minimizes the maximum variability of the pairwise differences
    crit.func <- function( A, W){
      inf.mat <- t(A) %*% W %*% A
      if( class( try( solve( inf.mat), silent = TRUE))[1] == "try-error"){
        return( -Inf)
      }else{
        inv.inf <- solve( inf.mat)
        all.pairs.var <- NULL
        for( i in 1:(nrow( inv.inf) - 1)){
          for( j in (i + 1):nrow( inv.inf)){
            all.pairs.var <- c( all.pairs.var, (inv.inf[i,i] + inv.inf[j,j] - 2*inv.inf[i,j]))
          }
        }
        return( -max( all.pairs.var))
      }
    }
  }else if( criterion == 'negD'){ # Finds the worst design for the D-optimal criteria
    crit.func <- function( A, W){
      inf.mat <- t(A) %*% W %*% A
      if( class( try( solve( inf.mat), silent = TRUE))[1] == "try-error"){
        return( -Inf)
      }else{
        return( -det( t(A) %*% W %*% A))
      }
    }
  }else if( criterion == 'negA'){ # Finds the worst design for the A-optimal criteria
    crit.func <- function( A, W){
      inf.mat <- t(A) %*% W %*% A
      if( class( try( solve( inf.mat), silent = TRUE))[1] == "try-error"){
        return( -Inf)
      }else{
        return( sum( diag( solve( inf.mat))))
      }
    }
  }else if( criterion == 'random'){ # Finds a random design that's also identifiable
    n.pairs <- length( starting.rows[ -which( starting.rows %in% which( rowSums( abs( A.use)) == 1))])
    repeat{
      current.rows <- c( sample( (1:nrow( A.use))[-which( rowSums( abs( A.use)) == 1)], n.pairs), which( rowSums( abs( A.use)) == 1))
      inf.mat <- t( A.use[ current.rows,]) %*% W.use[ current.rows, current.rows] %*% A.use[ current.rows,]
      if( !isTRUE( all.equal( det( inf.mat), 0))){
        if( class( try( solve( inf.mat), silent = TRUE))[1] != "try-error"){
          break
        }
      }
    }
    return(
      list(
        A = A.use[ current.rows,],
        W = W.use[ current.rows, current.rows],
        criterion = NA,
        rows = current.rows
      )
    )
  }else{
    warning( 'Unknown Criterion')
    return(NULL)
  }
  
  # Initialize the design matrix
  current.rows <- starting.rows
  
  # Keep track of iterations
  current.iter <- 1
  
  # Run the algorithm
  repeat{
    
    # Compute the criterion for the current design
    current.D <- crit.func( A.use[ current.rows,], W.use[ current.rows, current.rows])
    
    if( messages){
      message( paste0( 'Iteration ', current.iter, ' (criterion = ', signif( current.D, 4), ')'))
    }
    
    # Find all possible couples of current rows and potential rows
    candidate.rows <- setdiff( 1:nrow( A.use), current.rows)
    couples <- expand.grid( current = current.rows, candidate = candidate.rows, D = NA)
    
    # Compute the difference in criteria between all possible exchanges
    for( current.couple in 1:nrow( couples)){
      current.design <- A.use[ c( current.rows[ -which( current.rows == couples$current[ current.couple])], couples$candidate[ current.couple]),]
      current.weights <- W.use[ c( current.rows[ -which( current.rows == couples$current[ current.couple])], couples$candidate[ current.couple]), c(current.rows[-which(current.rows == couples$current[current.couple])], couples$candidate[current.couple])]
      couples$D[ current.couple] <- crit.func( current.design, current.weights)
    }
    couples$D.diff <- couples$D - current.D
    
    # If no exchanges improve the design, stop.  Otherwise, make the best swap and continue the algorithm
    if( sum( couples$D.diff > 0) == 0){
      break
    }else{
      current.rows <- c( current.rows[ -which( current.rows == couples$current[ which.max( couples$D.diff)[1]])], couples$candidate[ which.max( couples$D.diff)[1]])
      current.iter <- current.iter + 1
    }
  }
  
  # Return the optimal design, weight matrix, criteria, and rows
  return(
    list(
      A = A.use[ current.rows,],
      W = W.use[ current.rows,current.rows],
      criterion = current.D,
      rows = current.rows
    )
  )
}
