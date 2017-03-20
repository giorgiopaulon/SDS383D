
my.dist <- function(x, y){
  # Function that computes the distance matrix D, where D_{ij} is x_i - y_j
  # -----------------------------------------------------------------------
  # Args: 
  #   - x: first vector of points 
  #   - y: second vector of points 
  # Returns: 
  #   - D: matrix of pairwise distances
  # -----------------------------------
  n <- length(x)
  m <- length(y)

  D <- array(NA, dim = c(n, m))

  # In order to speed up the code, we vectorize the distance computations: we create
  # two matrices, X and Y, by replicating the vectors x and y. The only operation we do
  # is the difference between the two matrices.
  X <- matrix(rep(x, m), nrow = n, ncol = m, byrow = F)
  Y <- matrix(rep(y, n), nrow = n, ncol = m, byrow = T)

  D <- abs(X - Y)

  return (D)
}


C.squaredexp <- function(x, y, b, tau1sq, tau2sq){
  # Function that computes the squared exponential covariance function
  # ------------------------------------------------------------------
  # Args: 
  #   - x: first vector of points
  #   - y: second vector of points
  #   - b, tau1sq, tau2sq: parameters of the SE covariance function
  # Returns: 
  #   - CSE: squared exponential covariance matrix
  # ----------------------------------------------
  n <- length(x)
  m <- length(y)
  CSE <- array(NA, dim = c(n, m))

  d <- my.dist(x, y)
  CSE <- tau1sq * exp(-(1/2)*(d/b)^2)
  diag(CSE) <- diag(CSE) + tau2sq

  return(CSE)
}


C.matern52 <- function(x, y, b, tau1sq, tau2sq){
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

  d <- my.dist(x, y)
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


pred.GP <- function(x, y, xstar, b, tau1sq, tau2sq, sigma2, type = 'squaredexp'){
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
  n <- length(x)
  
  if (type == 'squaredexp'){
    Cxx <- C.squaredexp(x, x, b, tau1sq, tau2sq)
    Czx <- C.squaredexp(xstar, x, b, tau1sq, tau2sq)
    Czz <- C.squaredexp(xstar, xstar, b, tau1sq, tau2sq)
  }
  else if (type == 'matern'){
    Cxx <- C.matern52(x, x, b, tau1sq, tau2sq)
    Czx <- C.matern52(xstar, x, b, tau1sq, tau2sq)
    Czz <- C.matern52(xstar, xstar, b, tau1sq, tau2sq)
  }
  else{
    stop ('Select a valid Covariance function family')
  }
  
  Vinv = solve(Cxx + diag(sigma2, n))
  
  pred.mean <- Czx %*% Vinv %*% y
  pred.var <- Czz - Czx %*% Vinv %*% t(Czx)
  
  return (list('pred.mean' = pred.mean, 'pred.var' = pred.var))
}


marginal.likelihood <- function(x, y, b, tau1sq, tau2sq, sigma2, type = 'squaredexp'){
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
  n <- length(x)
  
  if (type == 'squaredexp')
    C <- C.squaredexp(x, x, b, tau1sq, tau2sq)
  else if (type == 'matern')
    C <- C.matern52(x, x, b, tau1sq, tau2sq)
  else
    stop ('Select a valid Covariance function family')
    
  V = C + diag(sigma2, n)
  
  return (dmvnorm(y,  mean = rep(0, n), sigma = V, log=TRUE))
}
