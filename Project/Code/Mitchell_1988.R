
setwd('~/Desktop/Semester 2/Statistical Modeling II/Project/')

library(RColorBrewer)

rm(list=ls())
source('Code/Mitchell_1988_fun.R')

# Load the data
X <- read.table('./Data/prostate.txt')[,1:8]
y <- read.table('./Data/prostate.txt')[,9]

n <- nrow(X)
p <- ncol(X)

# We need to center and standardize the variables that are not categorical
X[,-5] <- as.matrix(scale(X[,-5]))
# We do not need to center and scale the response variable


# X does not include the intercept: it will be included in the function that computes 
# the posterior distributions for each model. We here assume that the intercept is 
# exempt from deletion

# # Plot the Spike and Slab setting as described by Mitchell and Beauchamp (1988)
# h0 <- 10
# h1 <- 1
# par(mar=c(4,4,2,2), family = 'Palatino')
# plot(0, h0, type = 'p', pch = 16, xlim = c(-20, 20), ylim = c(0, h0), col = 'firebrick2', ylab = 'Density', xlab = bquote(beta[j]), xaxt = 'n')
# polygon(x = c(-10, -10, 10, 10), y = c(0, h1, h1, 0), lwd = 2, density = 2, col = 'darkorange2', lty = 1, angle = 45)
# abline(v=0, h=0)
# axis(1, at = c(-10, 10), labels = as.expression(c(bquote(-f[j]), bquote(f[j]))))
# text(x = 2, y = h0 - 0.5, labels = bquote(h[0~j]))
# text(x = 5, y = h1 + 0.5, labels = bquote(h[1~j]))
# text(x = -2.5, y = h1 - 0.5, labels = bquote(2~f[j]~h[1~j]))

# Choose a value for the penalty parameter (ratio of the heights of the spike and of the 
# slab)
gamma <- 10
wi <- post.prob.model(X, y, gamma)$wi
# The model 0 (wi[1]) is the model with the only intercept, the model 1 (wi[2]), the one 
# with intercept and beta_1, ...


# ========================
# MEDIAN PROBABILITY MODEL
# ========================
# Retain the variables whose posterior probability of being
# nonzero is larger than 0.5
post.prob <- array(NA, dim = p)
for (j in 1:p){
  post.model <- post.prob.model(X, y, gamma)
  wi <- post.model$wi
  select <- post.model$select
  post.prob[j] <- sum(wi[which(select[,j] == 1)])
}
par(mar=c(4,4,2,2), family = 'Palatino')
barplot(post.prob, names.arg=1:p, xlab = "Variable", ylab = "Probability")
abline(h = 0.5, lwd = 2, lty = 4, col = 'firebrick2')
# Beta_1, Beta_2, Beta_5 are retained



# =========================
# HIGHEST POSTERIOR DENSITY
# =========================
# Take the model with highest posterior density
# Visualize the posterior probabilities for each model
par(mar=c(4,4,2,2), family = 'Palatino')
plot(1:2^p-1, wi, type = 'h', lwd = 2, xlab = 'Model', ylab = 'Posterior Probability')

# Find the optimal model (we need to convert the model number to binary number in order 
# to understand which variables have been selected)
opt.model <- which.max(wi)
as.binary(opt.model - 1, p)[p:1]
# Beta_1, Beta_2, Beta_5 are retained



# Let us now plot how the variable selection behaves according to the penalty parameter
# gamma: 
# first compute the posterior probabilities for each beta according to the gamma value
gamma.grid <- exp(seq(-3, 3, length.out = 50))
post.prob <- array(NA, dim = c(length(gamma.grid), p))
for (i in 1:length(gamma.grid)){
  for (j in 1:p){
    post.model <- post.prob.model(X, y, gamma.grid[i])
    wi <- post.model$wi
    select <- post.model$select
    post.prob[i,j] <- sum(wi[which(select[,j] == 1)])
  }
}

# Then plot how the posterior probability of each beta changes according to the value
# of gamma
par(mar=c(4,4.2,2,2), family = 'Palatino', cex = 1.1)
cols <- colorRampPalette(brewer.pal(9, "Set1"))(p)
plot(gamma.grid, post.prob[,1], type = 'b', pch = 16, cex = 0.6, ylim = c(0,1), col = cols[1], xlab = bquote(gamma), ylab = bquote(P(beta[j] != 0)))
lines(gamma.grid, post.prob[,1], lwd = 2, col = cols[1])
text(x = gamma.grid[length(gamma.grid)] + 0.5, y = post.prob[length(gamma.grid),1], labels = 1, col = cols[1], cex = 0.7, lwd = 3)
for (j in 2:p){
  lines(gamma.grid, post.prob[,j], type = 'b', pch = 16, cex = 0.7, col = cols[j])
  lines(gamma.grid, post.prob[,j], lwd = 2, col = cols[j])
  text(x = gamma.grid[length(gamma.grid)] + 0.5, y = post.prob[length(gamma.grid),j], labels = bquote(.(j)), col = cols[j], cex = 0.7, lwd = 3)
}

