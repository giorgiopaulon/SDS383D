
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



C.squaredexp <- function(x, y, A, ls){
  # Function that computes the squared exponential covariance function
  # ------------------------------------------------------------------
  # Args: 
  #   - x: first vector of points
  #   - y: second vector of points
  #   - b, tau1sq: parameters of the SE covariance function
  # Returns: 
  #   - CSE: squared exponential covariance matrix
  # ----------------------------------------------
  n <- length(x)
  m <- length(y)
  CSE <- array(NA, dim = c(n, m))

  d <- my.dist(x, y)
  CSE <- A * exp(-(1/2) * (d/ls)^2)

  return(CSE)
}



data.cov <- function(y, hyperpar){

  hyperpar.h <- hyperpar[1:2]
  hyperpar.n <- hyperpar[3:4]
  hyperpar.r <- hyperpar[5:6]
  sigma <- hyperpar[7]
  
  # Choose one gene to create the upper left block of sigma. We can pick one gene 
  # because all of the genes are observed at the same times and therefore the blocks will
  # be the same
  ni <- length(unique(y$gene))
  yn <- y[y$gene == y$gene[1],]
  
  repl <- unique(y$replicate)
  nr <- length(repl)
  
  # We build the full covariance matrix starting from the building blocks: first Kr, 
  # the within-replicate covariance matrix (12*12), then the Kn blocks, the within-gene
  # covariance matrix(36*36) and finally the full covariance matrix
  # ===========
  # Building Kr
  # ===========
  Kr.list <- list()
  for (k in 1:nr){
    times <- yn$time[yn$replicate == repl[k]]
    Kr.list[[k]] <- C.squaredexp(times, times, 
                                 hyperpar.r[1], hyperpar.r[2])
  }
  
  tn <- yn$time   # vector of times for all replicates in one gene
  D <- length(tn) # number of observations across replicates for one gene
  
  # ===========
  # Building Kn
  # ===========
  Kn <- C.squaredexp(tn, tn, hyperpar.n[1], hyperpar.n[2])
  
  # ===========
  # Building Kh
  # ===========
  # Within-cluster covariance
  t <- y$time   # vector of times for all replicates in all genes
  Kh <- C.squaredexp(t, t, hyperpar.h[1], hyperpar.h[2])
  
  # Covariance matrix for one gene (it is identical for the others)
  Sigma.n <- Kn + bdiag(Kr.list) + sigma^2 * Diagonal(D)
  
  # Covariance matrix for all genes in the cluster
  Sigma <- as.matrix(Kh + bdiag(rep(list(Sigma.n), ni)))

  return(Sigma)
}



marginal.negloglikelihood <- function(hyperpar, y){
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
  n <- dim(y)[1]
  
  V = data.cov(y, exp(hyperpar))
  
  return (- dmvnorm(y$log2exp,  mean = rep(0, n), sigma = V, log=TRUE))
}




