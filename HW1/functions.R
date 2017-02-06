

rmultinorm <- function(n, mean, sigma){
  p <- length(mean)
  res <- array(NA, dim = c(n, p))
  
  for (j in 1:p){
    res[,j] <- rnorm(n)
  }
  L <- chol(sigma)
  res <- res %*% L
  res <- t(mean + t(res))
}


norm.loglik <- function(params, X){
  n <- dim(X)[1]
  p <- dim(X)[2]
  mean <- params[1:p]
  lag <- 1
  sigma <- matrix(0, nrow = p, ncol = p)
  for (i in 1:p){
    for (j in i:p){
      sigma[i,j] <- sigma[j,i] <- params[p + lag]
      lag <- lag + 1
    }
  }
  
  if (det(sigma) < 0){
    loglik <- 1000
  }
  else{
    sigmainv <- solve(sigma)
    loglik <- array(NA, dim = n)
    for (i in 1:n){
      loglik[i] <- (X[i,] - mean) %*% sigmainv %*% (X[i,] - mean)
    }
    loglik <- (n/2) * log(det(sigma)) + (1/2) * sum(loglik)
  }
  
  return (loglik)
}

my.MLE <- function(sample){
  
  mu <- colMeans(sample)
  Sigma <- cov(sample)
  
  # res <- optim(par = c(0,0,1,0,1), fn = norm.loglik, X = sample, method = 'Nelder-Mead')
  # res$par
  # res$value
  return (list('mu' = mu, 'Sigma' = Sigma))
}

