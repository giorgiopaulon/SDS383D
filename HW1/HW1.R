
setwd('~/Desktop/Semester 2/Statistical Modeling II/')
rm(list=ls())

library(mlbench)
library(mvtnorm)

#########################
### LINEAR REGRESSION ###
#########################

# Load the data
ozone = data(Ozone, package='mlbench')

# Scrub the missing values
# Extract the relevant columns 
ozone = na.omit(Ozone)[,4:13]

y = ozone[,1]
x = as.matrix(ozone[,2:10])

# add an intercept
x = cbind(1,x)

# compute the estimator
betahat = solve(crossprod(x)) %*% t(x) %*% y

# Compute the standard deviations of each beta
sigmahat <- sum((y - x %*% betahat)^2)/(dim(x)[1] - dim(x)[2])
betacov = sqrt(sigmahat * diag(solve(crossprod(x))))

# Now compare to lm
# the 'minus 1' notation says not to fit an intercept (we've already hard-coded it as 
# an extra column)
lm1 = lm(y~x-1)

summary(lm1)
betacovlm = vcov(lm1)
sqrt(diag(betacovlm))

# Comparison between R's function and our estimate
comp <- cbind(betacov, sqrt(diag(betacovlm)))
colnames(comp) <- c("Proposed estimate", "lm package")
rownames(comp)[1] <- "Intercept"
# toLatex(xtable(comp, align=c('|c|c|c|'),digits=6))


#################
### BOOTSTRAP ###
#################

rm(list=ls())

# Load the data
ozone = data(Ozone, package='mlbench')

# Scrub the missing values
# Extract the relevant columns 
ozone = na.omit(Ozone)[,4:13]

y = ozone[,1]
x = as.matrix(ozone[,2:10])

# add an intercept
x = cbind(1,x)
n <- dim(x)[1]
p <- dim(x)[2] - 1

# FIRST METHOD:
# Nonparametric Bootstrap, by considering the vector (X, Y) as random

betahat = solve(crossprod(x)) %*% t(x) %*% y
sigmahat <- sum((y - x %*% betahat)^2)/(n - p + 1)
betacov = sigmahat * solve(crossprod(x))


B <- 10000
boot.betacov <- matrix(rep(0, (p+1)*(p+1)), p+1, p+1)
for (b in 1:B){
  idx <- sample(1:n, n, replace = T)
  betahat = solve(t(x[idx,]) %*% x[idx,]) %*% t(x[idx,]) %*% y[idx]
  sigmahat <- sum((y[idx] - x[idx,] %*% betahat)^2)/(n - p + 1)
  boot.betacov = boot.betacov + sigmahat * solve(crossprod(x[idx,]))
}
boot.betacov <- boot.betacov/B
MSE <- mean((betacov - boot.betacov)^2)


# SECOND METHOD:
# Nonparametric Bootstrap, by considering X fixed, and resampling only the residuals
# from their empirical cdf

betahat = solve(crossprod(x)) %*% t(x) %*% y
sigmahat <- sum((y - x %*% betahat)^2)/(n - p + 1)
betacov = sigmahat * solve(crossprod(x))

res <- y - x %*% betahat

B <- 10000
boot.betacov <- matrix(rep(0, (p+1)*(p+1)), p+1, p+1)
for (b in 1:B){
  res.new <- sample(res, n, replace = T)
  y.new <- x %*% betahat + res.new
  betahat.boot = solve(crossprod(x)) %*% t(x) %*% y.new
  sigmahat <- sum((y.new - x %*% betahat.boot)^2)/(n - p + 1)
  boot.betacov = boot.betacov + sigmahat * solve(crossprod(x))
}
boot.betacov <- boot.betacov/B
MSE <- mean((betacov - boot.betacov)^2)
toLatex(xtable(boot.betacov, align = '|c|c|c|c|c|c|c|c|c|c|c|', digits = 4))


# BOOTSTRAPPING THE MLE:

rm(list=ls())
source('./SDS383D/HW1/functions.R')

N <- 200
p <- 2
mu = c(-5, 5)
sigma = matrix(c(0.75, 0.2, 0.2, 0.75), 2, 2, byrow = T)
X <- rmvnorm(N, mu, sigma)

B <- 10000
boot.mu <- rep(0, p)
boot.sigma <- matrix(rep(0, p*p), p, p)
for (b in 1:B){
  boot.X <- X[sample(1:N, size = N, replace = T),]
  mle <- my.MLE(boot.X)
  boot.mu <- boot.mu + mle$mu
  boot.sigma <- boot.sigma + mle$Sigma
}
boot.mu <- boot.mu/B
boot.sigma <- boot.sigma/B





