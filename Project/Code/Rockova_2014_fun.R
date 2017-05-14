


EMVS <- function(X, y, nu0, nu1, a, b, nu, lambda, Niter, tol = 1E-10){
  # EMVS Algorithm for variable selection
  # -------------------------------------
  # Args: 
  #   - X: design matrix
  #   - y: response vector
  #   - nu0, nu1: spike and slab parameters
  #   - a, b: hyperparameters for the prior probability of inclusion
  #   - nu, lambda: precision gamma hyperparameters
  #   - Niter: number of iterations
  #   - tol: convergence tolerance
  # Returns: 
  #   - list of parameters at each iteration
  # -------------------------------
  n <- nrow(X)
  p <- ncol(X)
  
  # Chain initialization
  sigma.chain <- array(NA, dim = Niter) 
  betas.chain <- array(NA, dim = c(Niter, p)) 
  theta.chain <- array(NA, dim = Niter)
  sigma.chain[1] <- 1
  betas.chain[1,] <- rep(1, p)
  theta.chain[1] <- 0.5
  
  # Precache operations
  XtX <- crossprod(X)
  Xty <- crossprod(X, y)
  for (iter in 2:Niter){

    # E-STEP:
    ai <- dnorm(betas.chain[iter-1,], mean = 0, sd = sigma.chain[iter-1] * sqrt(nu1)) * theta.chain[iter-1]
    bi <- dnorm(betas.chain[iter-1,], mean = 0, sd = sigma.chain[iter-1] * sqrt(nu0)) * (1 - theta.chain[iter-1])
    p.star <- ai/(ai + bi)
    d.star <- (1 - p.star)/nu0 + p.star/nu1
      
    # M-STEP:
    # Beta update:
    betas.chain[iter,] <- solve(XtX + diag(d.star), Xty)
    
    # Sigma update
    sigma.chain[iter] <- sqrt((crossprod(y - X %*% betas.chain[iter,]) + 
                                 crossprod(sqrt(d.star) * betas.chain[iter,]) + 
                                 nu * lambda)  / (n + p + nu))
    
    # Theta update:
    theta.chain[iter] <- (sum(p.star) + a - 1) / (a + b + p - 2)
    
    # Convergence check:
    if (sum(abs((betas.chain[iter,] - betas.chain[iter-1,])/betas.chain[iter-1,])) < tol){
      break;
    }
  }
  
  return(list('betas.chain' = betas.chain[1:iter,], 'sigma.chain' = sigma.chain[1:iter], 'theta.chain' = theta.chain[1:iter]))
}



