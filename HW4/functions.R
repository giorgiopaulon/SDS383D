
rtrnorm <- function(n, mean, sd, a, b){
  # Sample from a truncated normal distribution on [a, b]
  # -----------------------------------------------------
  # Args:
  #   - n: number of samples needed
  #   - mean: mean of the normal distribution
  #   - sd: sd of the normal distribution
  #   - a: left bound of the truncation domain
  #   - b: right bound of the truncation domain
  # Returns:
  #   - x: sample from the truncated normal distribution
  # ----------------------------------------------------
  Fa <- pnorm(a, mean, sd)
  Fb <- pnorm(b, mean, sd)
  
  u <- runif(n, 0, 1)
  x <- qnorm(u * (Fb - Fa) + Fa, mean, sd)
  
  return(x)
}




gibbs.HLM <- function(yi, Xi, Zi, ni, store, Niter, burnin, thin){
  # Gibbs Sampler for the Hierarchical linear model
  # -----------------------------------------------
  # Args:
  #   - yi: list of response variables grouped by store
  #   - Xi: list of fixed effects covariate matrices grouped by store
  #   - Zi: list of random effects covariate matrices grouped by store
  #   - ni: group sample sizes
  #   - store: indicator for the groups
  #   - Niter, burnin, thin: MCMC parameters
  # Returns:
  #   - fit: MCMC output
  # --------------------

  n <- sum(ni)
  p <- dim(Xi[[1]])[2]
  q <- dim(Zi[[1]])[2]
  I <- length(ni)
  
  # Initialize the chain
  gammas.chain <- array(NA, dim = c(Niter, q, I))
  betas.chain <- array(NA, dim = c(Niter, p))
  Dinv.chain <- array(NA, dim = c(Niter, q, q))
  lambdas.chain <- array(NA, dim = Niter)
  gammas.chain[1,,] <- matrix(rep(0, I*q))
  betas.chain[1,] <- rep(0, p)
  Dinv.chain[1,,] <- diag(rep(0.1, q))
  lambdas.chain[1] <- 0.01
  
  # Hyperparameters
  nu <- q + 1
  psi <- diag(rep(1, q))
  lambda0 <- 0.01
  
  # Precaching operations
  SS <- 0
  for (i in 1:I){
    SS <- SS + crossprod(Xi[[i]])
  }
  SSinv <- solve(SS)
  
  # Gibbs sampler
  pb <- txtProgressBar(style = 3)
  for (iter in 2:Niter){
    
    # Update beta
    mean.post <- rep(0, p)
    for (i in 1:I){
      wi <- yi[[i]] - as.numeric(Zi[[i]] %*% gammas.chain[iter-1,,i])
      mean.post <- mean.post + crossprod(Xi[[i]], wi)
    }
    mean.post <- (lambdas.chain[iter-1] / (lambdas.chain[iter-1] + lambda0)) * SSinv %*% mean.post
    var.post <- SSinv/(lambdas.chain[iter-1] + lambda0)
    betas.chain[iter,] <- rmvnorm(1, mean.post, var.post)
    
    # Update gamma
    Dinv <- Dinv.chain[iter-1,,]
    for (i in 1:I){
      mi <- yi[[i]] - as.numeric(Xi[[i]] %*% betas.chain[iter,])
      var.post <- solve(lambdas.chain[iter-1] * crossprod(Zi[[i]]) + Dinv)
      mean.post <- as.numeric(lambdas.chain[iter-1] * var.post %*% crossprod(Zi[[i]], mi))
      gammas.chain[iter,,i] <- rmvnorm(1, mean.post, var.post)
    }
    
    # Update lambda
    S.gamma.beta <- 0
    for (i in 1:I){
      y.tilde <- yi[[i]] - as.numeric(Zi[[i]] %*% gammas.chain[iter,,i]) - 
        as.numeric(Xi[[i]] %*% betas.chain[iter,])
      S.gamma.beta <- S.gamma.beta + crossprod(y.tilde)
    }
    lambdas.chain[iter] <- rgamma(1, n/2, S.gamma.beta/2)
    
    # Update D
    S.gamma <- tcrossprod(gammas.chain[iter,,])
    Dinv.chain[iter,,] <- rwish(nu + I, psi + S.gamma)
    
    setTxtProgressBar(pb, iter/Niter)
  }
  # Thin the chains
  gammas.chain <- gammas.chain[seq(burnin + 1, Niter, by = thin),,]
  betas.chain <- betas.chain[seq(burnin + 1, Niter, by = thin),]
  Dinv.chain <- Dinv.chain[seq(burnin + 1, Niter, by = thin),,]
  lambdas.chain <- lambdas.chain[seq(burnin + 1, Niter, by = thin)]
  
  return(list('betas.chain' = betas.chain, 'gammas.chain' = gammas.chain, 'Dinv.chain' = Dinv.chain, 'lambdas.chain' = lambdas.chain))
}



gibbs.probit <- function(yi, Xi, Wi, ni, state, Niter, burnin, thin){
  # Gibbs Sampler for the Probit model
  # ----------------------------------
  # Args:
  #   - yi: list of response variables grouped by state
  #   - Xi: list of fixed effects covariate matrices grouped by state
  #   - Wi: list of random effects covariate matrices grouped by state
  #   - ni: group sample sizes
  #   - state: indicator for the groups
  #   - Niter, burnin, thin: MCMC parameters
  # Returns:
  #   - fit: MCMC output
  # --------------------
  n <- sum(ni)
  p <- dim(Xi[[1]])[2]
  q <- dim(Wi[[1]])[2]
  I <- length(ni)
  
  # Initialize the chain
  betas.chain <- array(NA, dim = c(Niter, p))
  gammas.chain <- array(NA, dim = c(Niter, q, I))
  Sigmainv.chain <- array(NA, dim = c(Niter, q, q))
  z.chain <- array(NA, dim = n)
  y.pred <- array(NA, dim = c(Niter, n))
  betas.chain[1,] <- matrix(rep(0, p))
  gammas.chain[1,,] <- matrix(rep(0, I*q))
  Sigmainv.chain[1,,] <- diag(rep(0.1, q))
  z.chain[y == 1] <- 1
  z.chain[y == 0] <- -1
  
  # Hyperparameters
  nu <- q + 1
  psi <- diag(rep(1, q))

  # Precaching operations
  SS <- 0
  for (i in 1:I){
    SS <- SS + crossprod(Xi[[i]])
  }
  SSinv <- solve(SS)
  
  pb <- txtProgressBar(style = 3)
  for (iter in 2:Niter){
    # Update beta
    var.post <- SSinv
    mean.post <- 0
    for (i in 1:I){
      bi <- z.chain[state == i] - as.numeric(Wi[[i]] %*% gammas.chain[iter-1,,i])
      mean.post <- mean.post + crossprod(Xi[[i]], bi)
    }
    mean.post <- SSinv %*% mean.post
    betas.chain[iter,] <- rmvnorm(1, mean.post, var.post)
    
    # Update the gamma_i's
    Sigmainv <- Sigmainv.chain[iter-1,,]
    for (i in 1:I){
      mi <- z.chain[state == i] - as.numeric(Xi[[i]] %*% betas.chain[iter,])
      var.post <- solve(crossprod(Wi[[i]]) + Sigmainv)
      mean.post <- as.numeric(var.post %*% crossprod(Wi[[i]], mi))
      gammas.chain[iter,,i] <- rmvnorm(1, mean.post, var.post)
    }
    
    # Update Sigma
    S.gamma <- tcrossprod(gammas.chain[iter,,])
    Sigmainv.chain[iter,,] <- rwish(nu + I, psi + S.gamma)
    
    # We impute the missing values
    #y[idx.na] <- as.numeric(z.chain[idx.na] <= 0)
    
    # Update the z_i's
    for (i in 1:I){
      low <- rep(0, ni[i]); low[yi[[i]] == 1] <- 0; low[yi[[i]] == 0] <- -Inf
      upp <- rep(0, ni[i]); upp[yi[[i]] == 1] <- +Inf; upp[yi[[i]] == 0] <- 0
      
      y.pred[iter, state == i] <- rnorm(ni[i], as.numeric(Xi[[i]] %*% betas.chain[iter,] - Wi[[i]] %*% gammas.chain[iter,,i]))
      z.chain[state == i] <- rtrnorm(ni[i], as.numeric(Xi[[i]] %*% betas.chain[iter,] - Wi[[i]] %*% gammas.chain[iter,,i]), rep(1, ni[i]), low, upp)
    }
    y.pred[iter, y.pred[iter,] > 0] <- 1
    y.pred[iter, y.pred[iter,] < 0] <- 0
    z.chain[z.chain == +Inf] <- 10
    z.chain[z.chain == -Inf] <- -10
    
    setTxtProgressBar(pb, iter/Niter)
  }
  # Thin the chains
  betas.chain <- betas.chain[seq(burnin + 1, Niter, by = thin),]
  gammas.chain <- gammas.chain[seq(burnin + 1, Niter, by = thin),,]
  Sigmainv.chain <- Sigmainv.chain[seq(burnin + 1, Niter, by = thin),,]
  y.pred <- y.pred[seq(burnin + 1, Niter, by = thin),]
  
  return(list('betas.chain' = betas.chain, 'gammas.chain' = gammas.chain, 'Sigmainv.chain' = Sigmainv.chain, 'y.pred' = y.pred))
}


