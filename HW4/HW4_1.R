
setwd('~/Desktop/Semester 2/Statistical Modeling II/')

rm(list=ls())
# ====================
# ==== Math Tests ====
# ====================

math <- read.csv(file = 'SDS383D-master/data/mathtest.csv')

y <- math$mathscore
ybar <- aggregate(y, list(math$school), mean)$x
ni <- as.numeric(table(math$school))
n <- sum(ni)
I <- length(unique(math$school))

# Let us see the distribution of the scores for each school
par(mar=c(2,2,1,1))
boxplot(y ~ math$school, col = 'gray', pch = 16, cex = 0.8, lwd = 1.2)
abline(h = mean(y), col = 'indianred3', lwd = 2)

# Let us plot the average scores for each school vs the sample size of that school
par(mar=c(4,4,2,2), family = 'Palatino', cex = 1.1)
plot(ni, ybar, pch = 16, xlab = 'Sample size', ylab = 'Group mean')
# We notice that extreme average values of the scores are obtained for school with few 
# students sampled. This happens because the distribution of ybar for each school has a 
# variance of sigma^2/ni: the smaller ni, the larger the variability of ybar around the 
# grand mean.


# Run the Gibbs Sampler
Niter <- 11000
burnin <- 1000
thin <- 2

# Initialize the chain
thetas.chain <- array(NA, dim = c(Niter, I))
mu.chain <- array(NA, dim = Niter)
sigma2.chain <- array(NA, dim = Niter)
tau2.chain <- array(NA, dim = Niter)
thetas.chain[1,] <- rep(0, I)
mu.chain[1] <- 0
sigma2.chain[1] <- 1
tau2.chain[1] <- 1

for (i in 2:Niter){
  # Update thetas
  var.post <- tau2.chain[i-1] * sigma2.chain[i-1] / (ni * tau2.chain[i-1] + 1)
  mean.post <- (mu.chain[i-1] + tau2.chain[i-1] * ni * ybar) / (ni * tau2.chain[i-1] + 1)
  thetas.chain[i,] <- rnorm(I, mean.post, sqrt(var.post))
  
  # Update mu
  theta.bar <- mean(thetas.chain[i,])
  mu.chain[i] <- rnorm(1, theta.bar, sqrt(sigma2.chain[i-1] * tau2.chain[i-1] / I))
  
  # Update sigma2
  S.theta <- sum((thetas.chain[i,] - mu.chain[i])^2)
  S.y <- sum((y - rep(thetas.chain[i,], times = ni))^2)
  rate.new <- (1/2) * (S.y + S.theta / tau2.chain[i-1])
  sigma2.chain[i] <- 1/rgamma(1, (n + I)/2, rate.new)
  
  # Update tau2
  rate.new <- S.theta / (2 * sigma2.chain[i])
  tau2.chain[i] <- 1/rgamma(1, I/2 - 1, rate.new)
}
# Thin the chains
thetas.chain <- thetas.chain[seq(burnin + 1, Niter, by = thin),]
mu.chain <- mu.chain[seq(burnin + 1, Niter, by = thin)]
sigma2.chain <- sigma2.chain[seq(burnin + 1, Niter, by = thin)]
tau2.chain <- tau2.chain[seq(burnin + 1, Niter, by = thin)]


# Let us see how the Bayesian estimates differ from the sample means
par(mar=c(4,4,2,2), family = 'Palatino', cex = 1.1)
plot(ybar, colMeans(thetas.chain), xlab=bquote(bar(y)), ylab=bquote(hat(theta)), pch = 16)
abline(0, 1, col = 'indianred3', lwd = 2)
# The slope of this line is smaller than 1, that is, high values of ybar_i correspond to
# slightly less high values of the Bayesian estimates of theta_i; low values 
# of ybar_i correspond to slightly less low values of the Bayesian estimates of 
# theta_i. This is the shrinkage effect towards the grand mean (partial pooling).
# The farther away from the mean you are, the more you are shrunk.

#pdf('./SDS383D/HW4/Notes/Img/shrinkage.pdf', width = 8.5, height = 6)
par(mar=c(4,4,2,2), family = 'Palatino', cex = 1.1)
kappa <- (ybar - colMeans(thetas.chain)) / ybar
plot(ni, abs(kappa), ylab=bquote(kappa), xlab="Sample size", pch = 16)
abline(h=0, lwd = 2, col = 'indianred3')
# Groups with low sample size get shrunk the most, whereas groups with large 
# sample size hardly get shrunk at all. The larger the sample size for a group, the more
# information we have for that group and the less information we need to borrow from the
# rest of the population.


# Let us see how the shrinkage effect affects the posterior estimates
par(mar=c(2,2,1,1))
boxplot(y ~ math$school, col = 'gray', pch = 16, cex = 0.8, lwd = 1.2)
points(1:I, colMeans(thetas.chain), col = 'indianred3', pch = 16, cex = 0.8)
abline(h = mean(mu.chain), col = 'dodgerblue3', lwd = 2)