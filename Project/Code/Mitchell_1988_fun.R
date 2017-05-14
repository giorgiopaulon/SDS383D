
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


post.prob.model <- function(X, y, gamma){
  # ------------------------------------------------------
  # Args: 
  #   - X: design matrix (no intercept has to be included)
  #   - y: vector of responses
  #   - gamma: ratio between the height of the spike and the height of the slab 
  #     (can be interpreted as a penalty)
  # Returns: 
  #   - weights for each model and selected models
  # ----------------------------------------------
  n <- nrow(X)
  p <- ncol(X) # p does not include the intercept
  
  # In the wi vector, the null model (model 0, which is in position 1 in the vector)
  # represents the model with only the intercept. Model 1 (in position 2) represents
  # the model with intercept and beta_1, and so on...
  log.wi <- array(NA, dim = 2^p)
  select <- array(NA, dim = c(2^p, p))
  for (i in 1:(2^p)){
    # Subselect the predictors according to the current model
    select[i,] <- as.binary(i - 1, p)[p:1]
    Xm <- as.matrix(cbind(rep(1, n), as.matrix(X[,select[i,] == 1])))
    km <- sum(select[i,]) + 1
    
    # Compute the auxiliary quantities beta^hat and S_m^2
    beta.hat <- solve(crossprod(Xm), crossprod(Xm, y))
    Sm2 <- crossprod(y - Xm %*% beta.hat)
    
    # Compute the log-probability for the current model
    log.wi[i] <- km * (- log(gamma) + 0.5 * log(pi)) + log(gamma(0.5 * (n - km))) - 
      0.5 * log(det(crossprod(Xm))) - 0.5 * (n - km) * log(Sm2)
  }
  
  # Normalization trick via log-sum-exp
  c.norm <- max(log.wi)
  denom <- sum(exp(log.wi - c.norm))
  wi <- exp(log.wi - c.norm)/denom
  
  return(list('wi' = wi, 'select' = select))
}
