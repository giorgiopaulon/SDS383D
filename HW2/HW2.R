
setwd('~/Desktop/Semester 2/Statistical Modeling II/')
rm(list=ls())

library(mvtnorm)

# pdf('./SDS383D/HW2/Notes/Img/name.pdf', width = 8.5, height = 6)

###########################################
### THE CONJUGATE GAUSSIAN LINEAR MODEL ###
###########################################
gdp <- read.csv(file = 'SDS383D-master/data/gdpgrowth.csv')

n <- dim(gdp)[1]
X <- cbind(rep(1, n), gdp$DEF60)
y <- gdp$GR6096
p <- dim(X)[2] - 1

# FREQUENTIST LINEAR REGRESSION
lm1 = lm(y ~ X-1)
summary(lm1)

par(mar=c(4,4,2,2))
plot(X[,2], y, pch = 16, cex = 0.8, asp = 1, xlab = 'Defense Spending', ylab = 'GDP Growth Rate')
abline(a = lm1$coefficients[1], b = lm1$coefficients[2], lwd = 2, col = 'firebrick3')

# BAYESIAN LINEAR REGRESSION
# Y | beta, omega \sim N(X %*% beta; (omega * Lambda^{-1}))
# beta | omega \sim N(m; (omega * K^{-1}))
# omega \sim Gamma(d/2; eta/2)

# Hyperparameters
Lambda <- diag(rep(1, n))
m <- rep(0, p + 1)
K <- diag(c(0.01, 0.01))
d <- 0.02
eta <- 0.02


# Compute the MLE
XtLambdaX <- t(X) %*% Lambda %*% X
beta.hat <- solve(XtLambdaX) %*% crossprod(X, Lambda) %*% y

# Update the hyperparameters for the posterior distribution
K.new <- K + XtLambdaX
m.new <- solve(K.new) %*% (K %*% m + XtLambdaX %*% beta.hat)
d.new <- d + n
S <- t(y - X %*% beta.hat) %*% Lambda %*% (y - X %*% beta.hat)
r <- t(m - beta.hat) %*% solve(solve(K) + solve(XtLambdaX)) %*% (m - beta.hat)
# r <- t(beta.hat) %*% XtLambdaX %*% beta.hat + t(m) %*% K %*% m - t(m.new) %*% K.new %*% m.new
eta.new <- as.numeric(eta + S + r)

# Plot the regression line: the model is conjugate so we already have the MAP estimates
par(mar=c(4,4,2,2))
plot(X[,2], y, pch = 16, cex = 0.8, asp = 1, xlab = 'Defense Spending', ylab = 'GDP Growth Rate')
abline(a = lm1$coefficients[1], b = lm1$coefficients[2], lwd = 2, col = 'firebrick3')
abline(a = m.new[1], b = m.new[2], lwd = 2, col = 'dodgerblue3')
legend('topright', legend = c('Frequentist LM','Bayesian LM'), lwd = 2, lty = 1, col = c('firebrick3', 'dodgerblue3'))


# If we did not have a closed form solution for the posterior, we should have integrated
# the MCMC sample.
# Grid on which we evaluate the bayesian regression line, that is the line
# y = E[beta[1] | data] + E[beta[2] | data] * x
samp.size <- 5000
x.gr <- seq(min(X[,p+1]) - 0.1*diff(range(X[,p+1])), max(X[,p+1]) + 0.1*diff(range(X[,p+1])), length.out = 200)

pred <- matrix(ncol = 200, nrow = samp.size)
beta.post <- rmvt(samp.size, m.new, sigma = (eta.new/d.new) * solve(K.new), df = d.new)
sig.post <- 1/rgamma(samp.size, d.new/2, eta.new/2)

for (i in 1:samp.size){
  pred[i,] <- beta.post[i,1] + beta.post[i,2] * x.gr + rnorm(1, mean = 0, sd = sqrt(sig.post))
}

# Draw the average regression line
predit <- apply(pred, 2, mean)
lines(x.gr, predit, type="l", col = 'darkorchid2', lwd=2)
predit.quant <- apply(pred, 2, quantile, prob=c(0.05,0.95)) 
lines(x.gr,predit.quant[1,], col = 'darkorchid2', lwd = 2,lty = 2)
lines(x.gr,predit.quant[2,], col = 'darkorchid2', lwd = 2,lty = 2)

# Contour plot for the marginals of the beta parameters
xgrid <- seq(-0.4, 0.4, length.out = 200)
ygrid <- seq(-0.4, 0.4, length.out = 200)
z.prior <- matrix(nrow = length(xgrid), ncol = length(ygrid))
z.post <- matrix(nrow = length(xgrid), ncol = length(ygrid))
for (i in 1:length(xgrid)){
  for (j in 1:length(ygrid)){
    z.prior[i,j] <- dmvt(c(xgrid[i], ygrid[j]), m, sigma = (eta/d) * solve(K), df = d)
    z.post[i,j] <- dmvt(c(xgrid[i], ygrid[j]), m.new, sigma = (eta.new/d.new) * solve(K.new), df = d.new)
  }
}

par(mfrow=c(1,2), mar=c(4,4,2,2), cex = 1.1)
grays <- gray((200:0)/200)
image(xgrid, ygrid, exp(z.prior), col = grays, xlab = bquote(beta[0]), ylab = bquote(beta[1]), main = bquote("Prior density of"~beta)) 
points(m[1], m[2], pch = 4, col = 'firebrick3', lwd = 2)
image(xgrid, ygrid, exp(z.post), col = grays, xlab = bquote(beta[0]), ylab = bquote(beta[1]), main = bquote("Posterior density of"~beta)) 
points(m.new[1], m.new[2], pch = 4, col = 'firebrick3', lwd = 2)

##################################
### A HEAVY TAILED ERROR MODEL ###
##################################

rm(list=ls())
source('SDS383D/HW2/gibbs.heavytails.R')

gdp <- read.csv(file = 'SDS383D-master/data/gdpgrowth.csv')

n <- dim(gdp)[1]
X <- cbind(rep(1, n), gdp$DEF60)
y <- gdp$GR6096
p <- dim(X)[2] - 1

# Decide the number of iterations for the Gibbs Sampler, the thinning and the burnin
niter <- 12000
burnin <- 2000
thin <- 2

# Set hyperparameters
m <- rep(0, p + 1)
K <- diag(c(0.01, 0.01))
d <- 0.02
eta <- 0.02
h <- 0.02

posterior <- gibbs.heavytails(X, y, m, K, d, eta, h)
betas <- posterior$betas
lambdas <- posterior$lambdas
omegas <- posterior$omegas

# Select the first 5 observations whose precision is lower
idx <- NULL
for (i in 1:5){
  idx <- c(idx, which(colMeans(lambdas) == sort(colMeans(lambdas))[i]))
}

# Display those countries
par(mar=c(4,4,2,2), family = 'Palatino', cex = 1.1)
plot(X[,2], y, pch = 16, cex = 0.8, asp = 1, xlab = 'Defense Spending', ylab = 'GDP Growth Rate')
abline(a = colMeans(betas)[1], b = colMeans(betas)[2], lwd = 2, col = 'dodgerblue3')
points(X[idx,2], y[idx], pch = 16, cex = 0.9, col = 'red')
text(X[idx,2], y[idx], as.character(gdp$CODE[idx]), cex=.75, adj=1.25)
# This problem can be reinterpreted as an outlier detection. The regression line is
# severely affected by the outliers, that cause a big uncertainty. 


# Plot the traceplots for beta
par(mar=c(4,4,1,1), mfrow = c(2,1), cex = 1.1, family = 'Palatino')
plot(betas[,1], type = 'l', xlab = 'Iterations', ylab = bquote(beta[0]))
plot(betas[,2], type = 'l', xlab = 'Iterations', ylab = bquote(beta[1]))

# Plot the marginals for beta
par(mfrow=c(2,1),mar=c(2,2,2,2), cex=1, family = 'Palatino')
for (i in 1:2){
  plot(density(betas[,i]), xlab='values',ylab='density',lwd=2, cex.main=1.2, main=bquote("Density of "~beta[.(i-1)]))
  xapp <- seq(quantile(betas[,i], prob=c(0.025)), quantile(betas[,i],prob=c(0.975)), length=100)
  polygon(c(quantile(betas[,i], prob=c(0.025)), xapp, quantile(betas[,i],prob=c(0.975))), c(0, approxfun(density(betas[,i]))(xapp), 0), col='lightskyblue2', lty=0)
  segments(x0 = quantile(betas[,i], prob=c(0.025)), y0 = 0, y1 = approxfun(density(betas[,i]))(quantile(betas[,i], prob=c(0.025))), lwd=2)
  segments(x0 = quantile(betas[,i], prob=c(0.5)), y0 = 0, y1 = approxfun(density(betas[,i]))(quantile(betas[,i], prob=c(0.5))), lwd=2)
  segments(x0 = quantile(betas[,i], prob=c(0.975)), y0 = 0, y1 = approxfun(density(betas[,i]))(quantile(betas[,i], prob=c(0.975))), lwd=2)
  lines(density(betas[,i]), lwd=2)
}


