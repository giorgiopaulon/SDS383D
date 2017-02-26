

gauss.kernel <- function(x){
  # Computes the Gaussian kernel function
  # -------------------------------------------------------------------------------
  # Args: 
  #   - x: distance between the two points, scaled by h
  # Returns: 
  #   - K: value of the Gaussian kernel
  # -----------------------------------
  K <- (1/sqrt(2*pi)) * exp(-(1/2) * x^2)
  
  return (K)
}

unif.kernel <- function(x){
  # Computes the Uniform kernel function 
  # -------------------------------------------------------------------------------
  # Args: 
  #   - distance: distance between the two points, scaled by h
  # Returns: 
  #   - K: value of the Gaussian kernel
  # -----------------------------------
  K <- (1/2) * as.numeric(abs(x) <= 1)
  
  return (K)
}

weights <- function(x, xnew, h, type = 'gaussian'){
  # Computes the weigths for the new observation xnew given the previous ones x
  # ---------------------------------------------------------------------------
  # Args: 
  #   - x: vector containing the points
  #   - xnew: new point for which we compute the weights
  #   - h: bandwidth
  # Returns: 
  #   - w: vector of weights for the new observation xnew
  # -----------------------------------------------------
  n <- length(x)
  w <- array(NA, dim = n)
  
  if (type == 'gaussian'){
    for (i in 1:n){
      w[i] <- (1/h) * gauss.kernel((x[i] - xnew)/h)
    }
  }
  else if (type == 'uniform'){
    for (i in 1:n){
      w[i] <- (1/h) * unif.kernel((x[i] - xnew)/h)
    }
  }
  
  return (w)
}

linsmooth <- function(x, y, xnew, h, type = 'gaussian'){
  # Computes the predicted values of the linear smoother for the values xnew
  # ------------------------------------------------------------------------
  # Args: 
  #   - x: vector of the data to fit the linear smoother
  #   - y: response vector to fit the linear smoother
  #   - xnew: grid on which we wish to evaluate the linear smoother
  #   - h: bandwidth
  # Returns: 
  #   - ynew: values of the linear smoother on the grid xnew
  # --------------------------------------------------------
  ynew <- array(NA, dim = length(xnew))
  
  for (i in 1:length(ynew)){
    w <- weights(x, xnew[i], h, type)
    if (sum(w) > 0)
      w <- w/sum(w)
    ynew[i] <- sum(w * y)
  }
  
  return(ynew)
}

cv.linsmooth <- function(x.train, y.train, x.test, y.test, hgrid){
  # Computes the linear smoother functions for several values of the bandwidth h
  # and the MSE for each case
  # ----------------------------------------------------------------------------
  # Args: 
  #   - x.train: vector of the data to fit the linear smoother
  #   - y.train: response vector to fit the linear smoother
  #   - x.test: vector of the test data to compute the MSE
  #   - y.test: response vector of the test data to compute the MSE
  #   - hgrid: grid of values for the bandwidth
  # Returns: 
  #   - y.pred: values of the linear smoother on the grid x.test for each of the models
  #   - MSE: MSE values for each one of the models
  # ----------------------------------------------
  k <- length(hgrid)
  
  y.new <- array(NA, dim = c(k, length(y.test)))
  MSE <- array(NA, dim = k)
  for (i in 1:k){
    y.new[i,] <- linsmooth(x.train, y.train, x.test, hgrid[i])
    MSE[i] <- sum((y.new[i,] - y.test)^2)
  }
  
  return (list('y.pred' = y.new, 'MSE' = MSE))
}



polregr <- function(x, y, xnew, D, h, type = 'gaussian'){
  # Computes the predicted values of the linear polynomial regression for the
  # values xnew
  # -------------------------------------------------------------------------
  # Args: 
  #   - x: vector of the data to fit the linear smoother
  #   - y: response vector to fit the linear smoother
  #   - xnew: grid on which we wish to evaluate the linear smoother
  #   - D: maximum degree of the polynomial
  #   - h: bandwidth
  # Returns: 
  #   - ynew: values of the linear smoother on the grid xnew
  #   - H: the projection matrix, i.e. yhat = H %*% y
  # --------------------------------------------------------
  n <- length(x)
  nstar <- length(xnew)
  
  # Initialization of the data structures
  ynew <- array(NA, dim = nstar)
  H <- array(NA, dim = c(nstar, n))
  Rx <- array(NA, dim = c(n, D+1))
  Rx[,1] <- rep(1, n) # Initialization of the R matrix
  for (i in 1:length(ynew)){ # loop through the observations to predict
    # Compute the weights
    W <- diag(weights(x, xnew[i], h, type))
    
    # We compute the Rx matrix for the target point xnew[i]
    if (D >= 1){
      Rx[,2:(D+1)] <- rep(x - xnew[i], D)
      for (d in 1:(D+1))
        Rx[,d] <- (Rx[,d])^(d-1)
    }
    
    # We solve the WLS problem
    RtW <- crossprod(Rx, W)        # precaching
    b <- solve(RtW %*% Rx) %*% RtW # solution of the WLS problem
    
    # Compute the projection matrix for the ith data point
    H[i,] <- b[1,]
    
    # Remember: at the target point, the estimated function is a_0^hat
    ynew[i] <- (b %*% y)[1]
  }
  return(list('ynew' = ynew, 'H' = H))
}


cv.polregr <- function(x, y, D, hgrid, type = 'gaussian'){
  # Computes the LOOCV error for the linear polynomial regression for several values 
  # of the bandwidth h
  # ------------------
  # Args: 
  #   - x: vector of the data to fit the linear smoother
  #   - y: response vector to fit the linear smoother
  #   - D: maximum degree of the polynomial
  #   - hgrid: grid of values for the bandwidth
  # Returns: 
  #   - LOOCV: leave-one-out cross validation error
  # -----------------------------------------------
  k <- length(hgrid)
  n <- length(x)
  H <- array(NA, dim = c(n, n))
  
  LOOCV <- array(NA, dim = k)
  for (h in 1:k){
    # We compute the H matrix for each possible value of the bandwidth
    H <- polregr(x, y, x, D, hgrid[h])$H
    
    # We use the LOOCV lemma
    LOOCV[h] <- sum((y - H %*% y)^2 / (1 - diag(H))^2)
  }
  
  return (LOOCV)
}



