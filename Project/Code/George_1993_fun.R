

as.binary <- function(n, p, base=2){
  # Conversion from integer to binary number
  # ----------------------------------------
  # Args: 
  #   - n: integer number to convert
  #   - p: total number of digits in the binary representation (length of the output)
  #   - base: base, 2 for binary
  # Returns: 
  #   - out: binary number
  # ----------------------
  out <- NULL
  while(n > 0){
    out <- c(n %% base, out)
    n <- n %/% base
  }
  
  out <- c(rep(0, p - length(out)), out)
  return(out)
}


gibbs.SVSS <- function(X, y, c, tau, nu0, lambda0, Niter, burnin, thin){
  # Stochastic Search variable selection as described in George
  # ----------------------------------------
  # Args: 
  #   - X: design matrix
  #   - y: response vector
  #   - c, tau: spike and slab parameters
  #   - nu0, lambda0: precision gamma hyperparameters
  #   - Niter, burnin, thin: Gibbs parameters
  # Returns: 
  #   - MCMC of the joint posterior
  # -------------------------------
  n <- nrow(X)
  p <- ncol(X) - 1
  
  # Chain initialization
  lambda.chain <- array(NA, dim = Niter) 
  betas.chain <- array(NA, dim = c(Niter, p+1)) 
  gammas.chain <- array(NA, dim = c(Niter, p+1))
  wi.chain <- array(NA, dim = c(Niter, p+1))
  lambda.chain[1] <- 0.1
  #betas.chain[1,] <- solve(crossprod(X), crossprod(X, y))
  betas.chain[1,] <- rep(0, p+1)
  gammas.chain[1,] <- rep(1, p+1)
  wi.chain[1,] <- rep(0.5, p+1)
  
  for (iter in 2:Niter){
    # D Matrix computation (diagonal matrix if R = I, with spike and slab variances on 
    # the diagonal)
    Dinv <- matrix(rep(0, (p+1)*(p+1)), p+1, p+1)
    diag(Dinv)[gammas.chain[iter-1,] == 0] <- 1/tau^2
    diag(Dinv)[gammas.chain[iter-1,] == 1] <- 1/(c * tau)^2
    
    # Beta update:
    Sigma.star <- solve(crossprod(X) * lambda.chain[iter-1] + Dinv)
    mu.star <- lambda.chain[iter-1] * Sigma.star %*% crossprod(X, y)
    betas.chain[iter,] <- rmvnorm(1, mu.star, Sigma.star)
    
    # Lambda update
    nu.star <- (n + nu0)/2
    lambda.star <- (nu0 * lambda0 + crossprod(y - X %*% betas.chain[iter,]))/2
    lambda.chain[iter] <- rgamma(1, nu.star, lambda.star)
    
    # Gamma_i's update in random order to speed up the convergence
    perm <- sample(1:(p+1), size = p+1, replace = FALSE)
    for (j in perm){
      a <- dnorm(betas.chain[iter,j], 0, tau * c) * wi.chain[iter-1,j]
      b <- dnorm(betas.chain[iter,j], 0, tau) * (1 - wi.chain[iter-1,j])
      gammas.chain[iter,j] <- sample(0:1, 1, prob = c(b, a))
    }
    gammas.chain[iter,1] <- 1 # force the intercept in the model
    
    # w_i's update
    wi.chain[iter,] <- rbeta(p+1, gammas.chain[iter,] + 1, 1 - gammas.chain[iter,] + 1)
  }
  
  gammas.chain <- gammas.chain[seq(burnin+1, Niter, by = thin),]
  betas.chain <- betas.chain[seq(burnin+1, Niter, by = thin),]
  lambda.chain <- lambda.chain[seq(burnin+1, Niter, by = thin)]
  wi.chain <- wi.chain[seq(burnin+1, Niter, by = thin),]
  
  return(list('gammas.chain' = gammas.chain, 'betas.chain' = betas.chain, 'lambda.chain' = lambda.chain, 'wi.chain' = wi.chain))
}

