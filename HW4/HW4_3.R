
setwd('~/Desktop/Semester 2/Statistical Modeling II/')
rm(list=ls())
library(mvtnorm)
library(MCMCpack)
library(lattice)
library(ggplot2)
library(lme4)
library(plotrix)
# ===========================================================
# ==== A hierarchical probit model via data augmentation ====
# ===========================================================


polls <- read.csv(file = 'SDS383D-master/data/polls.csv')
polls <- polls[,-c(1:3,10)]
idx.na <- is.na(polls$bush)
polls <- polls[-which(idx.na),]

levels(polls$edu) <- c('NoHS','HS','SomeColl','Bacc')

n <- dim(polls)[1]
ni <- as.numeric(table(polls$state))
I <- length(unique(polls$state))
X <- model.matrix( ~ edu + age + female + black, polls)
y <- polls$bush
ybar <- aggregate(bush ~ state, data = polls, mean)$bush

# Let us plot the average voting for each state vs the sample size of that state
par(mar=c(4,4,2,2), family = 'Palatino', cex = 1.1)
plot(ni, ybar, pch = 16, xlab = 'Sample size', ylab = 'Group mean')


# Let us see the distribution of the votes for each state
cols <- ifelse(ybar > 0.5, 'firebrick1', 'dodgerblue1')
cols[ybar == 0.5] <- 'gray'
par(mar=c(2,2,1,1), family = 'Palatino')
barchart(ybar ~ levels(state), data = polls, scales=list(x=list(rot=90,cex=0.8)), col = cols)


# Fit a Generalized HLM with random effects
hlm <- glmer(bush ~ edu + age + female + black + (1 | state), data = polls, family = binomial(link = 'probit'))
summary(hlm)
# For each single state, we are fitting a random intercept. 
betas <- coef(hlm)$state

randeff = ranef(hlm, condVar = TRUE)
plrand = dotplot(randeff)
plot(plrand$state)

y.pred <- array(NA, dim = n)
for (i in 1:I){
  y.pred[as.numeric(polls$state) == i] <- X[as.numeric(polls$state) == i,] %*% as.numeric(betas[i,])
}
y.pred <- ifelse(y.pred < 0, 0, 1)
sum(y == y.pred)/n

# What kinds of descriptives can we do for these categorical data?
# Does the marginal effect of being female voting for republican change from state to 
# state? Relative Risk of voting R for male vs female in NY and in Ca (log-odds ratio).
# Do contingency table!

rm(list=ls())
library(mvtnorm)
library(MCMCpack)
library(ggplot2)
# ===========================================================
# ==== A hierarchical probit model via data augmentation ====
# ===========================================================

source('./SDS383D/HW4/functions.R')

polls <- read.csv(file = 'SDS383D-master/data/polls.csv')
polls <- polls[,-c(1:3, 10)]

idx.na <- is.na(polls$bush)
polls <- polls[-which(idx.na),]

y <- polls$bush
ni <- as.numeric(table(polls$state))
n <- sum(ni)
I <- length(unique(polls$state))
state <- as.numeric(polls$state)

xtabs(~ state + bush, data = polls)


# Matrix of the random effects
X <- model.matrix( ~ edu + age + female + black, polls)

# Matrix of the random effects
W <- model.matrix( ~ edu + age + female + black, polls)


p <- ifelse(length(dim(X)[2]) > 0, dim(X)[2], 1)
q <- ifelse(length(dim(W)[2]) > 0, dim(W)[2], 1)

# Initialize the store data and design matrices in a list
yi = list()
Xi = list()
Wi = list()
for (i in 1:I){	
  yi[[i]] = y[state==i]
  Xi[[i]] = X[state==i,]
  Wi[[i]] = W[state==i,]
}

# Run the Gibbs Sampler
Niter <- 6000
burnin <- 1000
thin <- 1

fit <- gibbs.HLM(yi, Xi, Wi, ni, state, Niter, burnin, thin)

# save(fit, file = './SDS383D/HW4/pollsMCMC.Rdata')
# load('SDS383D/HW4/pollsMCMC.Rdata')

par(mar=c(4,4,2,2), mfrow = c(2,2), family = 'Palatino')
for (i in 1:p)
  plot(fit$betas.chain[,i], type = 'l', xlab = 'Iterations', ylab = bquote(beta[.(i)]))

par(mar=c(2,2,1,1), mfrow = c(5,5), family = 'Palatino')
for (i in 1:49)
  plot(fit$gammas.chain[,1,i], type = 'l', xlab = 'Iterations', ylab = bquote(gamma[.(i)]))

par(mar=c(2,2,1,1), mfrow = c(4,4), family = 'Palatino')
for (i in 1:p)
  plot(fit$Sigmainv.chain[,1,i], type = 'l', xlab = 'Iterations', ylab = bquote(gamma[.(i)]))


gammas.quantile <- apply(fit$gammas.chain, c(2,3), quantile, prob = c(0.025, 0.5, 0.975))
betas.quantile <- apply(fit$betas.chain, 2, quantile, prob = c(0.025, 0.5, 0.975))

par(mar = c(4,2,2,2), family = 'Palatino', cex = 0.8)
plotCI(x = gammas.quantile[2,1,], ui = gammas.quantile[3,1,], li = gammas.quantile[1,1,], pch=c(rep(21,I)), scol = rep('dodgerblue3',I), col = rep('dodgerblue3',I),
       ylab='', xlab = '', lwd=2, xaxt = 'n', cex = 1.2)
abline(h = 0, lwd = 2, lty = 2)
axis(1, at = 1:I, labels = levels(polls$state), las = 2)



# Compute accuracy?
gammas.post <- apply(fit$gammas.chain, c(2,3), mean)
betas.post <- colMeans(fit$betas.chain)

prob.pred <- colMeans(fit$y.pred)

# Majority vote
y.pred <- ifelse(prob.pred < 0.5, 0, 1)
sum(y == y.pred)/n

# ROC curve
threshold <- seq(0, 1, length.out = 100)
FP <- array(NA, dim = 100)
TP <- array(NA, dim = 100)
for (i in 1:length(threshold)){
  y.pred <- ifelse(prob.pred < threshold[i], 0, 1)
  FP[i] <- sum(y.pred == 1 & y == 0)/n
  TP[i] <- sum(y.pred == 1 & y == 1)/n
}


