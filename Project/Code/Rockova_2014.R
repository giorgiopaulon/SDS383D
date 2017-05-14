
setwd('~/Desktop/Semester 2/Statistical Modeling II/Project/')

library(RColorBrewer)
library(mvtnorm)
library(plotrix)

rm(list=ls())
source('Code/Rockova_2014_fun.R')


# X <- read.table('~/Desktop/communities.txt', sep = ',')[,-(1:5)]
# y <- X[,123]
# X <- X[,-c(97:113,117:120,122,123)]
# X <- X[-131,]
# y <- y[-131]
# X$V31 <- as.numeric(as.character(X$V31))


# # Load the data
# X <- read.table('./Data/prostate.txt')[,1:8]
# y <- read.table('./Data/prostate.txt')[,9]
# 
# n <- nrow(X)
# p <- ncol(X)

# # We need to center and standardize the variables that are not categorical
# X[,-5] <- as.matrix(scale(X[,-5], center = F))
# X <- as.matrix(X)
# X <- as.matrix(scale(X, center = F))
# X <- as.matrix(X)

# # De-mean the response
# y <- y - mean(y)

n <- 100
p <- 1000

Sigma <- array(NA, dim = c(p, p))
for (i in 1:p){
  for (j in 1:p){
    Sigma[i,j] <- 0.6^(abs(i - j))
  }
}

beta <- rep(0, p)
beta[1] <- 3
beta[2] <- 2
beta[3] <- 1
sigma <- sqrt(3)

X <- rmvnorm(n, mean = rep(0, p), sigma = Sigma)
y <- rnorm(n, mean = X %*% beta, sd = sigma)


a <- b <- 1
nu <- lambda <- 1
nu0 <- 0.5
nu1 <- 100

maxiter <- 30
  
fit <- EMVS(X, y, nu0, nu1, a, b, nu, lambda, maxiter)
Niter <- length(fit$sigma.chain)

plot(fit$betas.chain[Niter,], pch = 16)
abline(0, 1)

ci <- dnorm(fit$betas.chain[Niter,], mean = 0, sd = fit$sigma.chain[Niter] * sqrt(nu1)) * fit$theta.chain[Niter]
di <- dnorm(fit$betas.chain[Niter,], mean = 0, sd = fit$sigma.chain[Niter] * sqrt(nu0)) * (1 - fit$theta.chain[Niter])
prob.inclusion <- ci/(ci + di)


nu0.grid <- seq(0.05, 0.5, by = 0.01)
betas <- array(NA, dim = c(length(nu0.grid), p))
thresh <- array(NA, dim = length(nu0.grid))
for (i in 1:length(nu0.grid)){
  fit <- EMVS(X, y, nu0.grid[i], nu1, a, b, nu, lambda, maxiter)
  Niter <- length(fit$sigma.chain)
  betas[i,] <- fit$betas.chain[Niter,]
  c2 <- nu1/nu0.grid[i]
  thresh[i] <- fit$sigma.chain[Niter] * sqrt(2 * nu0.grid[i] * log(sqrt(c2) * (1 - fit$theta.chain[Niter]) / fit$theta.chain[Niter]) * c2/(c2 - 1))
  cat("Iteration ", i, " out of ", length(nu0.grid), "\n")
}

save(nu0.grid, betas, thresh, file = 'EMVSplot.Rdata')
matplot(nu0.grid, betas, type = 'b', lwd = 2, pch = 16, cex = 0.8, ylim = c(-1,3.5))
lines(nu0.grid, thresh, lwd = 2, col = 'seagreen2')
lines(nu0.grid, - thresh, lwd = 2, col = 'seagreen2')

