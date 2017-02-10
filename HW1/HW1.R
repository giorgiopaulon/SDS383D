
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
source('SDS383D/HW1/functions.R')

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


betacov.boot <- boot.xy(x, y, 10000)

XtXinv <- solve(crossprod(x))
betahat = XtXinv %*% t(x) %*% y
sigmahat <- sum((y - x %*% betahat)^2)/(n - p + 1)
betacov = sigmahat * solve(crossprod(x))


# SECOND METHOD:
# Nonparametric Bootstrap, by considering X fixed, and resampling only the residuals
# from their empirical cdf

betacov.boot <- boot.res(x, y, 10000)

XtXinv <- solve(crossprod(x))
betahat = XtXinv %*% t(x) %*% y
sigmahat <- sum((y - x %*% betahat)^2)/(n - p + 1)
betacov = sigmahat * XtXinv


# BOOTSTRAPPING THE MLE:

rm(list=ls())
source('./SDS383D/HW1/functions.R')

N <- 200
p <- 2
mu = c(-5, 5)
s = matrix(c(0.75, 0.2, 0.2, 0.75), 2, 2, byrow = T)
X <- rmultinorm(N, mu, sigma)

B <- 10000
mu.hat <- array(NA, dim = c(B, p))
sigma.hat <- array(NA, dim = c(B, p * p))
for (b in 1:B){
  boot.X <- X[sample(1:N, size = N, replace = T),]
  mu.hat[b,] <- as.numeric(colMeans(boot.X))
  sigma.hat[b,] <- as.numeric(cov(boot.X))
}

colMeans(mu.hat)

par(mar=c(2,4,2,2), mfrow = c(2,1))
hist(mu.hat[,1], border = 'white', col = 'gray', xlab = '', nclass = 20, freq = F, main = bquote('Histogram of'~hat(mu)[1]))
abline(v = colMeans(mu.hat)[1], lwd = 2, col = 'firebrick3')
hist(mu.hat[,2], border = 'white', col = 'gray', xlab = '', nclass = 20, freq = F, main = bquote('Histogram of'~hat(mu)[2]))
abline(v = colMeans(mu.hat)[2], lwd = 2, col = 'firebrick3')

matrix(colMeans(sigma.hat), 2, 2, byrow = T)



