
my.dist <- function(x, y, b1, b2){
  # Function that computes the distance matrix D, where D_{ij} is x_i - y_j
  # -----------------------------------------------------------------------
  # Args: 
  #   - x: first vector of points 
  #   - y: second vector of points 
  # Returns: 
  #   - D: matrix of pairwise distances
  # -----------------------------------
  n <- dim(x)[1]
  m <- dim(y)[1]
  
  D <- array(NA, dim = c(n, m))
  
  X1 <- matrix(rep(x[,1], m), nrow = n, ncol = m, byrow = F)
  Y1 <- matrix(rep(y[,1], n), nrow = n, ncol = m, byrow = T)
  X2 <- matrix(rep(x[,2], m), nrow = n, ncol = m, byrow = F)
  Y2 <- matrix(rep(y[,2], n), nrow = n, ncol = m, byrow = T)
  
  D <- sqrt( ((X1 - Y1) / b1)^2 + ((X2 - Y2) / b2)^2 )
  
  return (D)
}


C.squaredexp <- function(x, y, b, tau1sq, tau2sq, b1, b2){
  # Function that computes the squared exponential covariance function
  # ------------------------------------------------------------------
  # Args: 
  #   - x: first vector of points
  #   - y: second vector of points
  #   - b, tau1sq, tau2sq: parameters of the SE covariance function
  # Returns: 
  #   - CSE: squared exponential covariance matrix
  # ----------------------------------------------
  n <- dim(x)[1]
  m <- dim(y)[1]
  CSE <- array(NA, dim = c(n, m))
  
  d <- my.dist(x, y, b1, b2)
  CSE <- tau1sq * exp(-(1/2)*(d/b)^2)
  diag(CSE) <- diag(CSE) + tau2sq
  
  return(CSE)
}


C.matern52 <- function(x, y, b, tau1sq, tau2sq, b1, b2){
  # Function that computes the Matern 5/2 covariance function
  # ---------------------------------------------------------
  # Args: 
  #   - x: first vector of points
  #   - y: second vector of points
  #   - b, tau1sq, tau2sq: parameters of the M52 covariance function
  # Returns: 
  #   - CM52: Matern 5/2 covariance matrix
  n <- length(x)
  m <- length(y)
  CM52 <- array(NA, dim = c(n, m))
  
  d <- my.dist(x, y, b1, b2)
  CM52 <- tau1sq * (1 + sqrt(5)*d/b + (5*d^2)/(3*b^2)) * exp(- sqrt(5) * d/b)
  diag(CM52) <- diag(CM52) + tau2sq
  
  return(CM52)
}


rgaussproc <- function(x, m, C){
  # Samples one trajectory from a Gaussian Process
  # ----------------------------------------------
  # Args: 
  #   - x: vector where evaluating the process
  #   - m: mean function for the process evaluated at x
  #   - C: covariance function of the process evaluated at x
  # Returns: 
  #   - values of the process at the target points x
  # ------------------------------------------------
  return (as.numeric(rmvnorm(1, mean = m, sigma = C)))
}


pred.GP <- function(x, y, xstar, b, tau1sq, tau2sq, sigma2, b1, b2, type = 'squaredexp'){
  # Computes the log-marginal likelihood of the Gaussian Process for a given set of
  # parameters
  # ----------
  # Args: 
  #   - x: vector where we observe the process
  #   - y: observed values for y
  #   - xstar: values where to predict mean and variance of f(xstar)
  #   - b, tau1sq, tau2sq: parameters for the Gaussian Process
  #   - sigma2: noise (unknown) of the data
  # Returns: 
  #   - pred.mean: predicted mean of f(xstar)
  #   - pred.var: predicted variance of f(xstar)
  # --------------------------------------------
  n <- dim(X)[1]
  
  if (type == 'squaredexp'){
    Cxx <- C.squaredexp(x, x, b, tau1sq, tau2sq, b1, b2)
    Czx <- C.squaredexp(xstar, x, b, tau1sq, tau2sq, b1, b2)
    Czz <- C.squaredexp(xstar, xstar, b, tau1sq, tau2sq, b1, b2)
  }
  else if (type == 'matern'){
    Cxx <- C.matern52(x, x, b, tau1sq, tau2sq, b1, b2)
    Czx <- C.matern52(xstar, x, b, tau1sq, tau2sq, b1, b2)
    Czz <- C.matern52(xstar, xstar, b, tau1sq, tau2sq, b1, b2)
  }
  else{
    stop ('Select a valid Covariance function family')
  }
  
  Vinv = solve(Cxx + diag(sigma2, n))
  
  pred.mean <- Czx %*% Vinv %*% y
  pred.var <- Czz - Czx %*% Vinv %*% t(Czx)
  
  return (pred.mean)
}


marginal.likelihood <- function(x, y, b, tau1sq, tau2sq, sigma2, b1, b2, type = 'squaredexp'){
  # Computes the log-marginal likelihood of the Gaussian Process for a given set of
  # parameters
  # ----------
  # Args: 
  #   - x: vector where we observe the process
  #   - y: observed values for y
  #   - b, tau1sq, tau2sq: parameters for the Gaussian Process
  #   - sigma2: noise (unknown) of the data
  # Returns: 
  #   - the value of the log-marginal likelihood
  # --------------------------------------------
  n <- dim(x)[1]
  
  if (type == 'squaredexp')
    C <- C.squaredexp(x, x, b, tau1sq, tau2sq, b1, b2)
  else if (type == 'matern')
    C <- C.matern52(x, x, b, tau1sq, tau2sq, b1, b2)
  else
    stop ('Select a valid Covariance function family')
  
  V = C + diag(sigma2, n)
  
  return (dmvnorm(y,  mean = rep(0, n), sigma = V, log=TRUE))
}

