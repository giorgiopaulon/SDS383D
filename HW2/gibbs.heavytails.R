
gibbs.heavytails <- function(X, y, m, K, d, eta, h, niter = 6000, burnin = 1000, thin = 1){
  # Gibbs sampling algorithm to sample from the posterior of the heavy tailed model
  # -------------------------------------------------------------------------------
  # Args: 
  #   - X: design matrix
  #   - y: response vector
  #   - m, K, d, eta, h: hyperparameters
  #   - niter: number of iterations
  #   - burnin: initial samples to discard
  #   - thin: thinning factor of the chain
  # Returns: 
  #   - the posterior chains
  # ------------------------
  p <- dim(X)[2] - 1
  n <- dim(X)[1]

  # Set MCMC objects
  betas <- array(NA, dim = c(niter, p+1))
  lambdas <- array(NA, dim = c(niter, n))
  omegas <- array(NA, dim = niter)
  
  # Set initial values
  betas[1,] <- rep(0, p+1)
  lambdas[1,] <- rep(2, n)
  omegas[1] <- 2

  for (iter in 2:niter){
    # Compute auxiliary quantities for the Normal distribution
    Lambda <- diag(lambdas[iter-1,])
    XtLambdaX <- crossprod(X, Lambda) %*% X
    beta.hat <- solve(XtLambdaX) %*% crossprod(X, Lambda) %*% y
    K.new <- K + XtLambdaX
    m.new <- solve(K.new) %*% (K %*% m + XtLambdaX %*% beta.hat)
    
    # Update beta step
    betas[iter,] <- rmvnorm(1, mean = m.new, sigma = (1/omegas[iter-1])*solve(K.new))
    
    # Compute auxiliary quantities for the Gamma distribution
    S <- t(y - X %*% beta.hat) %*% Lambda %*% (y - X %*% beta.hat)
    r <- t(m - beta.hat) %*% solve(solve(K) + solve(XtLambdaX)) %*% (m - beta.hat)
    eta.new <- as.numeric(eta + S + r)
    
    # Update omega step
    omegas[iter] <- rgamma(1, (d + n)/2, eta.new/2)
    
    # Compute auxiliary quantities for the Gamma distribution
    new.rate <- (1/2) * (omegas[iter] * (y - X %*% betas[iter,])^2 + h)
    
    # Update lambda step
    lambdas[iter,] <- rgamma(n, rep((h+1)/2, n), as.numeric(new.rate))
  }

  # Thin the chains
  betas <- betas[seq(burnin + 1, niter, by = thin),]
  lambdas <- lambdas[seq(burnin + 1, niter, by = thin),]
  omegas <- omegas[seq(burnin + 1, niter, by = thin)]
  
  return(list('betas' = betas, 'lambdas' = lambdas, 'omegas' = omegas))
}
