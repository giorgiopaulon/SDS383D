
setwd('~/Desktop/Semester 2/Statistical Modeling II/Project/')

library(RColorBrewer)
library(mvtnorm)
library(plotrix)

rm(list=ls())
source('Code/George_1993_fun.R')

# Load the data
X <- read.table('./Data/prostate.txt')[,1:8]
y <- read.table('./Data/prostate.txt')[,9]

n <- nrow(X)
p <- ncol(X)

# We need to center and standardize the variables that are not categorical
X[,-5] <- as.matrix(scale(X[,-5]))
# We add an intercept
X <- as.matrix(cbind(rep(1, n), X))
# We do not need to center and scale the response variable


# Delta is the threshold to declare a beta_i significant
delta <- 0.1
c <- 10

# We compute consequently tau
eps <- sqrt(2*log(c)*c^2/(c^2-1))
tau <- delta/eps

# # Plot the chosen spike & slab configuration
# par(mar=c(2,4,2,2), cex = 1.1, family = 'Palatino')
# xgrid <- seq(-1, 1, length.out = 200)
# plot(xgrid, dnorm(xgrid, 0, tau), xlab = '', col='firebrick2', type = 'l', lwd=2, main='Spike & Slab SVSS', ylab='Density')
# lines(xgrid, dnorm(xgrid, 0, c*tau), col='dodgerblue3', lwd=2)
# lines(xgrid, 0.5 * dnorm(xgrid, 0, tau) + 0.5 * dnorm(xgrid, 0, c*tau), col= 'darkorchid2', type='l', lwd=2)
# legend('topright', c("Spike", "Slab", "Mixture"), col = c('firebrick2', 'dodgerblue3', 'darkorchid2'), lty = rep(1, 3), lwd = 2)
# text(-1, 11, labels = bquote(tau^2 == .(tau^2)), col='firebrick2', pos=4)
# text(-1, 10, labels = bquote(c^2~tau^2 == .(c^2 * tau^2)), col='dodgerblue3', pos=4)
# abline(v=c(-delta, delta), col='turquoise',lwd=2, lty=2)


# MCMC setting
Niter <- 30000
burnin <- 5000
thin <- 5

# Hyperparameters for the precision lambda
nu0 <- 1
lambda0 <- 1

SVSS <- gibbs.SVSS(X, y, c, tau, nu0, lambda0, Niter, burnin, thin)


par(mar=c(4,4,2,2), mfrow = c(3,3), family = 'Palatino')
for (i in 1:(p+1))
  plot(SVSS$betas.chain[,i], type = 'l', xlab = 'Iterations', ylab = bquote(beta[.(i)]))
dev.off()

par(mar=c(4,4,2,2), mfrow = c(3,3), family = 'Palatino')
for (i in 1:(p+1))
  plot(SVSS$wi.chain[,i], type = 'l', xlab = 'Iterations', ylab = bquote(w[.(i)]))
dev.off()

gammas.post <- apply(SVSS$gammas.chain, 2, mean)

# ========================
# MEDIAN PROBABILITY MODEL
# ========================
# Retain the variables whose posterior probability of being
# nonzero is larger than 0.5: just analyze the posterior mean of the gamma_i's
par(mar=c(4,4,2,2), family = 'Palatino')
barplot(gammas.post[2:(p+1)], names.arg=1:p, xlab = "Variable", ylab = "Probability")
abline(h = 0.5, lwd = 2, lty = 4, col = 'firebrick2')



# =========================
# HIGHEST POSTERIOR DENSITY
# =========================
# Take the model with highest posterior density
# Visualize the posterior probabilities for each model
# These three lines convert the 01's sequences in the model numbers
pot <- c(0, 2^(1:p-1))
dec <- SVSS$gammas.chain%*%pot
# The model 0 is the model with the only intercept, the model 1, the one 
# with intercept and beta_1, ...

length(table(dec)) # number of visited models: it will change according to c and tau

# Plot the probability for each model
dec <- factor(dec, levels = 1:2^p-1)
plot(1:2^p-1, as.numeric(table(dec)/length(dec)), type = 'h', lwd = 2, xlab = 'Model', ylab = 'Posterior Probability')

# Find the top10 models
dec.top10 <- as.numeric(names(sort(table(dec), decreasing=T)))[1:10]
dec.top10

# Find the optimal model and convert it to binary
opt.model <- dec.top10[1]
bin.max <- as.binary(opt.model, p)[p:1]
bin.max



# ==============
# HARD SHRINKAGE
# ==============
# Just look at the marginal posterior densities for each beta and see
# if the 0 is in the CI 95%

betas.quant <- apply(SVSS$betas.chain, 2, quantile, prob = c(0.025, 0.5, 0.975))

par(mar = c(4,4,2,2), family = 'Palatino', cex = 1.1)
plotCI(1:p, betas.quant[2,2:(p+1)], ui = betas.quant[3,2:(p+1)], li = betas.quant[1,2:(p+1)], pch=c(rep(21,p)), scol='dodgerblue2', col = 'dodgerblue2', lwd=2, 
       xaxt = 'n', ylab='', xlab = '')
axis(1, at = 1:p, labels = as.expression(c(bquote(beta[1]), bquote(beta[2]), bquote(beta[3]), bquote(beta[4]), 
                                           bquote(beta[5]), bquote(beta[6]), bquote(beta[7]), bquote(beta[8]))))
abline(h=0, lwd = 2, lty = 2, col='firebrick2')



