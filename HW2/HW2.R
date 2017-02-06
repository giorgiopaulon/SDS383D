
setwd('~/Desktop/Semester 2/Statistical Modeling II/')
rm(list=ls())

library(mvtnorm)

###########################################
### THE CONJUGATE GAUSSIAN LINEAR MODEL ###
###########################################
gdp <- read.csv(file = 'SDS383D-master/data/gdpgrowth.csv')

n <- dim(gdp)[1]
X <- cbind(rep(1, n), gdp$DEF60)
y <- gdp$GR6096
p <- dim(X)[2] - 1

# Frequentist model
lm1 = lm(y ~ X-1)
summary(lm1)

par(mar=c(4,4,2,2))
plot(X[,2], y, pch = 16, cex = 0.8, asp = 1, xlab = 'X', ylab = 'y')
abline(a = lm1$coefficients[1], b = lm1$coefficients[2], lwd = 2, col = 'firebrick3')

# Bayesian linear regression
# Y | beta, omega \sim N(X %*% beta; (omega * Lambda^{-1}))
# beta | omega \sim N(m; (omega * K^{-1}))
# omega \sim Gamma(d/2; eta/2)

Lambda <- diag(rep(1, n))
m <- rep(0, p + 1)
K <- diag(c(0.01, 0.01))
d <- 0.01
eta <- 0.01

XtLambdaX <- t(X) %*% Lambda %*% X

# Compute the MLE
beta.hat <- solve(XtLambdaX) %*% crossprod(X, Lambda) %*% y

# Update the hyperparameters for the posterior distribution
K.new <- K + XtLambdaX
m.new <- solve(K.new) %*% (K %*% m + XtLambdaX %*% beta.hat)
d.new <- d + n
S <- (1/2) * t(y - X %*% beta.hat) %*% Lambda %*% (y - X %*% beta.hat)
r <- t(m - beta.hat) %*% solve(solve(K) + solve(XtLambdaX)) %*% (m - beta.hat)
eta.new <- as.numeric(eta + S + r)


par(mar=c(4,4,2,2))
plot(X[,2], y, pch = 16, cex = 0.8, asp = 1, xlab = 'X', ylab = 'y')
abline(a = lm1$coefficients[1], b = lm1$coefficients[2], lwd = 2, col = 'firebrick3')
abline(a = m.new[1], b = m.new[2], lwd = 2, col = 'dodgerblue3')
legend('topright', legend = c('Frequentist LM','Bayesian LM'), lwd = 2, lty = 1, col = c('firebrick3', 'dodgerblue3'))



# Contour plot for the marginals of the beta parameters
xgrid <- seq(-2, 2, length.out = 200)
ygrid <- seq(-2, 2, length.out = 200)
z.prior <- matrix(nrow = length(xgrid), ncol = length(ygrid))
z.post <- matrix(nrow = length(xgrid), ncol = length(ygrid))
for (i in 1:length(xgrid)){
  for (j in 1:length(ygrid)){
    z.prior[i,j] <- dmvt(c(xgrid[i], ygrid[j]), m, (eta/d) * solve(K), df = d)
    z.post[i,j] <- dmvt(c(xgrid[i], ygrid[j]), m.new, (eta.new/d.new) * solve(K.new), df = d.new)
  }
}

# par(mfrow=c(1,2), mar = c(4,4,2,2))
# contour(xgrid, ygrid, nlevels = 5, z.prior, xlim = range(xgrid),
#         ylim = range(ygrid), main = bquote("Log-prior density of"~beta),
#         xlab = bquote(beta[0]), ylab = bquote(beta[1]), lwd = 2, col = 'dodgerblue3')
# points(m[1], m[2], pch = 16, col = 'firebrick3')
# contour(xgrid, ygrid, levels = seq(-6.5, -8.5, by = -.5), z.post, xlim = range(xgrid),
#         ylim = range(ygrid), main = bquote("Log-posterior density of"~beta),
#         xlab = bquote(beta[0]), ylab = bquote(beta[1]), lwd = 2, col = 'dodgerblue3')
# points(m.new[1], m.new[2], pch = 16, col = 'firebrick3')

par(mfrow=c(1,2), mar=c(4,4,2,2))
grays <- gray((100:0)/100)
image(xgrid, ygrid, exp(z.prior), col = grays, xlab = bquote(beta[0]), ylab = bquote(beta[1]), main = bquote("Prior density of"~beta)) 
points(m[1], m[2], pch = 16, col = 'firebrick3')
image(xgrid, ygrid, exp(z.post), col = grays, xlab = bquote(beta[0]), ylab = bquote(beta[1]), main = bquote("Posterior density of"~beta)) 
points(m.new[1], m.new[2], pch = 16, col = 'firebrick3')


##################################
### A HEAVY TAILED ERROR MODEL ###
##################################





