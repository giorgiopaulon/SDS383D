

boot.res <- function(X, y, B = 1e4){
  # Computes the bootstrapped estimated covariance of beta.hat, sampling the residuals
  # ----------------------------------------------------------------------------------
  # Args: 
  #   - X: matrix of the features (n*p)
  #   - y: vector of responses
  #   - B: number of bootstrap iterations
  # Returns: 
  #   - cov.boot: estimated covariance matrix of the MLE
  # ----------------------------------------------------
  n <- nrow(X)
  p <- ncol(X) - 1
  
  # Calculate beta.hat
  XtXinv <- solve(crossprod(X))
  beta.hat <- XtXinv %*% crossprod(X, y)
  
  # Calculate predicted values and residuals
  y.hat <- X %*% beta.hat
  res <- y - y.hat
  
  beta.boot <- matrix(nrow = B, ncol = p + 1)
  for(i in 1:B) {
    idx <- sample(1:n, n, replace = T)
    res.boot <- res[idx]
    
    y.boot <- y.hat + res.boot
    beta.boot[i, ] <- XtXinv %*% crossprod(X, y.boot)
  }
  cov.boot <- cov(beta.boot)
  
  return(cov.boot)
}

boot.xy <- function(X, y, B = 1e4){
  # Computes the bootstrapped estimated covariance of beta.hat, sampling (X, y)
  # ---------------------------------------------------------------------------
  # Args: 
  #   - X: matrix of the features (n*p)
  #   - y: vector of responses
  #   - B: number of bootstrap iterations
  # Returns: 
  #   - cov.boot: estimated covariance matrix of the MLE
  # ----------------------------------------------------
  n <- nrow(X)
  p <- ncol(X) - 1
  
  beta.boot <- matrix(nrow = B, ncol = p + 1)
  for(i in 1:B) {
    idx <- sample(1:n, n, replace = T)
    
    X.boot <- X[idx, ]
    y.boot <- y[idx]
    
    XtX.boot <- crossprod(X.boot)
    
    beta.boot[i, ] <- as.numeric(solve(XtX.boot) %*% crossprod(X.boot, y.boot))
  }
  cov.boot <- cov(beta.boot)
  
  return(cov.boot)
}


rmultinorm <- function(n, mean, sigma){
  # Randomly generates a sample from a multivariate normal distribution via Choleski
  # --------------------------------------------------------------------------------
  # Args: 
  #   - n: number of desired samples
  #   - mean: vector of means
  #   - sigma: covariance matrix
  # Returns: 
  #   - res: random sample N(mean, sigma)
  # -------------------------------------
  p <- length(mean)
  res <- array(NA, dim = c(n, p))
  
  for (j in 1:p){
    res[,j] <- rnorm(n)
  }
  L <- chol(sigma)
  res <- res %*% L
  res <- t(mean + t(res))
  
  return (res)
}



