

setwd('~/Desktop/Semester 2/Statistical Modeling II/')

rm(list=ls())
library(lme4)
library(lattice)
# ======================================
# ==== Price elasticity of demand:  ====
# ======================================

cheese <- read.csv(file = 'SDS383D-master/data/cheese.csv')

cheese$logvol <- log(cheese$vol)
cheese$logprice <- log(cheese$price)

# Let us see the distribution of the scores for each school
par(mar=c(2,2,1,1), family = 'Palatino')
boxplot(cheese$vol ~ cheese$store, col = 'gray', pch = 16, cex = 0.8, lwd = 1.2, names = F)
abline(h = mean(cheese$vol), col = 'indianred3', lwd = 2)

# Probably the suggested log-log transformation is useful
par(mar=c(2,2,1,1), family = 'Palatino')
boxplot(cheese$logvol ~ cheese$store, col = 'gray', pch = 16, cex = 0.8, lwd = 1.2, names = F)
abline(h = mean(cheese$logvol), col = 'indianred3', lwd = 2)

# Is there an effect of the display? Here all the stores are grouped
par(mar=c(4,4,2,2), family = 'Palatino')
plot(cheese$logprice, cheese$logvol, pch = 21, bg = c('dodgerblue2','firebrick2')[cheese$disp+1], cex = 0.8, xlab = 'Price', ylab = 'Quantity')
legend('topright', pch = 16, legend = c('No','Yes'), col = c('dodgerblue2','firebrick2'))
# When we have displays we sell more cheap cheese. Blue points seem to be shifted left: 
# display happens when the cheese is cheaper: confounding. Price is a confounder: it's
# correlated with both the predictor (display) and the response (quantity).
boxplot(log(price) ~ disp, data = cheese, col = 'gray', pch = 16)


# Use lattice to reproduce small multiples plots
cols <- c('dodgerblue2','firebrick2')
xyplot(logvol ~ logprice | store, data = cheese, group = disp, col = cols, col.line = cols, pch = 16, cex = 0.3, type = c("p", "r"), par.strip.text=list(cex=.35))
# We could order the plots on the basis of the steepness of the lines


# Fit a HLM with random effects
hlm = lmer(logvol ~ (1 + logprice + disp + logprice:disp | store), data = cheese)
summary(hlm)
# For each single store, we are fitting 4 parameters \beta_i. The prior is 
# \beta_i \sim N(\mu, \Sigma). The matrix in the output is Sigma.
# The intercepts (2.24) and the logprice (2.16) vary a lot. 
# The correlations are significant as well

betas <- coef(hlm)$store


# Example: plot 4 stores
xgrid <- seq(min(cheese$logprice), max(cheese$logprice), length.out = 100)
par(mar=c(2,2,2,2), mfrow = c(2,2))
for (i in c(1,7,11,61)){
  plot(logvol ~ logprice, data=subset(cheese, as.numeric(store) == i), xlim = c(min(cheese$logprice), max(cheese$logprice)), ylim = c(min(cheese$logvol), max(cheese$logvol)))
  points(logvol ~ logprice, pch=19, col='dodgerblue2', data=subset(cheese, as.numeric(store) == i & disp == 0))
  points(logvol ~ logprice, pch=19, col='firebrick2', data=subset(cheese, as.numeric(store) == i & disp == 1))
  lines(xgrid, betas[i,4] + xgrid*betas[i,1], col='blue', lwd = 2)
  lines(xgrid, betas[i,4] + betas[i,2] + xgrid*(betas[i,1] + betas[i,3]), col='red', lwd = 2)
}

# How to determine an optimal price for each store?
# We will see an example in the Bayesian setting


rm(list=ls())
library(MCMCpack)
library(mvtnorm)
library(Matrix)
# ======================================================
# ==== Price elasticity of demand: Bayesian setting ====
# ======================================================

# Bayesian Hierarchical Model
source('SDS383D/HW4/functions.R')

cheese <- read.csv(file = 'SDS383D-master/data/cheese.csv')

y <- log(cheese$vol)                             #log-transformed data
ybar <- aggregate(y, list(cheese$store), mean)$x # aggregate means for groups
ni <- as.numeric(table(cheese$store))            # sample sizes for each group
n <- sum(ni)                                     # total number of observations
I <- length(unique(cheese$store))                # number of groups 
store <- as.numeric(cheese$store)                # indicators for the groups


# Matrix of the fixed effects
X <- cbind(rep(1, n), log(cheese$price), cheese$disp, log(cheese$price) * cheese$disp)

# Matrix of the random effects
Z <- cbind(rep(1, n), log(cheese$price), cheese$disp, log(cheese$price) * cheese$disp)

p <- dim(X)[2]
q <- dim(Z)[2]

# Initialize the store data and design matrices in a list
yi = list()
Xi = list()
Zi = list()
for (i in 1:I){	
  yi[[i]] = y[store==i]
  Xi[[i]] = X[store==i,]
  Zi[[i]] = Z[store==i,]
}

# Alternative parametrization for the Linear model with random effects:
# Y_i | beta_i, gamma_i, lambda \sim N (Z_i gamma_i, lambda^{-1} I)
#       gamma_i | D \iid N (0, D)
#       D \sim IW (nu, Psi)
#       beta \sim N (0, lambda_0^{-1}, I)
#       lambda \sim 1/lambda

# Run the Gibbs Sampler
Niter <- 6000
burnin <- 1000
thin <- 1

fit <- gibbs.HLM(yi, Xi, Zi, ni, store, Niter, burnin, thin)
  

# Analyse the convergence via traceplots
par(mar=c(4,4,2,2), mfrow = c(1,1), family = 'Palatino')
plot(fit$lambdas.chain, type = 'l', xlab = 'Iterations', ylab = bquote(lambda))
par(mar=c(4,4,2,2), mfrow = c(2,2), family = 'Palatino')
for (i in 1:4)
  plot(fit$betas.chain[,i], type = 'l', xlab = 'Iterations', ylab = bquote(beta[.(i)]))

par(mar=c(2,2,1,1), mfrow = c(5,5), family = 'Palatino')
for (i in 1:25)
  plot(fit$gammas.chain[,1,i], type = 'l', xlab = 'Iterations', ylab = bquote(beta[.(i)]))


# Compute the posterior point estimates
betas.post <- colMeans(fit$betas.chain)
gammas.post <- apply(fit$gammas.chain, c(2,3), mean)

# Draw the demand curves
cols <- c('dodgerblue2','firebrick2')
par(mar=c(2,2,1,1), mfrow = c(5,5), family = 'Palatino')
for (i in 1:25){
  xgrid <- seq(min(X[store==i,2]), max(X[store==i,2]), length.out = 100)
  plot(exp(X[store==i,2]), exp(y[store == i]), type = 'p', col = cols[X[store==i,3]+1], pch = 16, cex = 0.8)
  lines(exp(xgrid), exp((betas.post[1] + gammas.post[1,i]) + xgrid * (betas.post[2] + gammas.post[2,i])), col = cols[1], lwd = 2)
  lines(exp(xgrid), exp((betas.post[1] + betas.post[3] + gammas.post[1,i]+ gammas.post[3,i]) + xgrid * (betas.post[2] + betas.post[4] + gammas.post[2,i] + gammas.post[4,i])), col = cols[2], lwd = 2)              
}


# Fit four stores that are peculiar (few data points for one of the categories)
cols <- c('dodgerblue2','firebrick2')
xgrid <- seq(min(X[,2]), max(X[,2]), length.out = 100)
par(mar=c(2,2,2,2), mfrow = c(2,2), cex = 0.8, family = 'Palatino')
for (i in c(1,9,18,34)){
  plot(X[store==i,2], y[store == i], type = 'p', main = as.character(cheese$store[i]), col = cols[X[store==i,3]+1], pch = 16, cex = 0.8, xlim = c(min(X[,2]), max(X[,2])), ylim = c(min(y), max(y)))
  lines(xgrid, (betas.post[1] + gammas.post[1,i]) + xgrid * (betas.post[2] + gammas.post[2,i]), col = cols[1], lwd = 2)
  lines(xgrid, (betas.post[1] + betas.post[3] + gammas.post[1,i]+ gammas.post[3,i]) + xgrid * (betas.post[2] + betas.post[4] + gammas.post[2,i] + gammas.post[4,i]), col = cols[2], lwd = 2)              
}

# Same four stores via OLS
cols <- c('dodgerblue2','firebrick2')
xgrid <- seq(min(X[,2]), max(X[,2]), length.out = 100)
par(mar=c(2,2,2,2), mfrow = c(2,2), cex = 0.8, family = 'Palatino')
for (i in c(1,9,18,34)){
  plot(X[store==i,2], y[store == i], type = 'p', main = as.character(cheese$store[i]), col = cols[X[store==i,3]+1], pch = 16, cex = 0.8, xlim = c(min(X[,2]), max(X[,2])), ylim = c(min(y), max(y)))
  fit = lm(y[store == i] ~ -1 + X[store==i,])
  betas.lm <- fit$coefficients
  lines(xgrid, (betas.lm[1]) + xgrid * (betas.lm[2]), col = cols[1], lwd = 2)
  lines(xgrid, (betas.lm[1]+ betas.lm[3]) + xgrid * (betas.lm[2] + betas.lm[4]), col = cols[2], lwd = 2)              
}
# xtabs(~store + disp, data=cheese)

# The shrinkage effect is evident for the stores whose OLS estimates were distorted.
