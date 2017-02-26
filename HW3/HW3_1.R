
setwd('~/Desktop/Semester 2/Statistical Modeling II/')
rm(list=ls())


#########################################
### CURVE FITTING BY LINEAR SMOOTHING ###
#########################################

source('SDS383D/HW3/smoothing.R')

# pdf('./SDS383D/HW3/Notes/Img/h.pdf', width = 8.5, height = 6)

# Data generation
n <- 50
x <- runif(n, 0, 2*pi)
eps <- rnorm(n, 0, 0.25)
y <- sin(x) + eps

# Centering step
x <- x - mean(x)
y <- y - mean(y)

par(mar=c(4,4,2,2), family = 'Palatino', cex = 1.1)
plot(x, y, pch = 16, cex = 0.8, xlab = 'x', ylab = 'y')

h <- 1

xnew <- seq(min(x), max(x), length.out = 100)
ynew <- linsmooth(x, y, xnew, h)

lines(xnew, ynew, col = 'firebrick3', lwd = 2)

# pdf('./SDS383D/HW3/Notes/Img/112.pdf', width = 8.5, height = 6)
# par(mar=c(4,4,2,2), family = 'Palatino', cex = 1.1)
# plot(x, y, pch = 16, cex = 0.8, xlab = 'x', ylab = 'y')
# 
# h <- c(.5, 1, 2, 5)
# colors <- c('dodgerblue3','indianred2','goldenrod2','seagreen3','darkorchid2')
# for (i in 1:length(h)){
#   ynew <- linsmooth(x, y, xnew, h[i], type = 'uniform')
#   lines(xnew, ynew, col = colors[i], lwd = 2)
# }
# legend('topright', legend = c('h = 0.5','h = 1','h = 2','h = 5'), lwd = 2, col = colors[1:4])
# dev.off()

########################
### CROSS VALIDATION ###
########################

rm(list=ls())

source('SDS383D/HW3/smoothing.R')


# Case 1: wiggly function with highly noisy observations
n <- 500
x <- runif(n, 0, 1)
eps <- rnorm(n, 0, 0.75)
y <- sin(10*pi*x) + eps

# Centering step
mx <- mean(x); x <- x - mx
y <- y - mean(y)

x.train <- x[1:(n/2)]
y.train <- y[1:(n/2)]
x.test <- x[(n/2 + 1):n]
y.test <- y[(n/2 + 1):n]

par(mar=c(4,4,2,2), family = 'Palatino', cex = 1.1)
plot(x.train, y.train, pch = 16, cex = 0.8, xlab = 'x', ylab = 'y', main = 'Wiggly function, noisy obs.')
lines(sort(x.train), sin(10*pi*(sort(x.train) + mx)), lwd = 2.5, col = 'seagreen3', lty = 2)

hgrid <- exp(seq(-10, 1, length.out = 50))
cv.lin <- cv.linsmooth(x.train, y.train, x.test, y.test, hgrid)
h.opt <- hgrid[which.min(cv.lin$MSE)]

xgrid <- seq(min(x.train), max(x.train), length.out = 100)
y.pred <- linsmooth(x.train, y.train, xgrid, h.opt)
lines(xgrid, y.pred, lwd = 2.5, col = 'firebrick2')


# Case 2: wiggly function with not so noisy observations
n <- 500
x <- runif(n, 0, 1)
eps <- rnorm(n, 0, 0.2)
y <- sin(10*pi*x) + eps

# Centering step
mx <- mean(x); x <- x - mx
y <- y - mean(y)

x.train <- x[1:(n/2)]
y.train <- y[1:(n/2)]
x.test <- x[(n/2 + 1):n]
y.test <- y[(n/2 + 1):n]

par(mar=c(4,4,2,2), family = 'Palatino', cex = 1.1)
plot(x.train, y.train, pch = 16, cex = 0.8, xlab = 'x', ylab = 'y', main = 'Wiggly function, not noisy obs.')
lines(sort(x.train), sin(10*pi*(sort(x.train) + mx)), lwd = 2.5, col = 'seagreen3', lty = 2)

hgrid <- exp(seq(-10, 1, length.out = 50))
cv.lin <- cv.linsmooth(x.train, y.train, x.test, y.test, hgrid)
h.opt <- hgrid[which.min(cv.lin$MSE)]

xgrid <- seq(min(x.train), max(x.train), length.out = 100)
y.pred <- linsmooth(x.train, y.train, xgrid, h.opt)
lines(xgrid, y.pred, lwd = 2.5, col = 'firebrick2')


# Case 3: smooth function with highly noisy observations
n <- 500
x <- runif(n, 0, 1)
eps <- rnorm(n, 0, 0.75)
y <- sin(2*pi*x) + eps

# Centering step
mx <- mean(x); x <- x - mx
y <- y - mean(y)

x.train <- x[1:(n/2)]
y.train <- y[1:(n/2)]
x.test <- x[(n/2 + 1):n]
y.test <- y[(n/2 + 1):n]

par(mar=c(4,4,2,2), family = 'Palatino', cex = 1.1)
plot(x.train, y.train, pch = 16, cex = 0.8, xlab = 'x', ylab = 'y', main = 'Smooth function, noisy obs.')
lines(sort(x.train), sin(2*pi*(sort(x.train) + mx)), lwd = 2.5, col = 'seagreen3', lty = 2)

hgrid <- exp(seq(-10, 1, length.out = 50))
cv.lin <- cv.linsmooth(x.train, y.train, x.test, y.test, hgrid)
h.opt <- hgrid[which.min(cv.lin$MSE)]

xgrid <- seq(min(x.train), max(x.train), length.out = 100)
y.pred <- linsmooth(x.train, y.train, xgrid, h.opt)
lines(xgrid, y.pred, lwd = 2.5, col = 'firebrick2')


# Case 4: smooth function with not so noisy observations
n <- 500
x <- runif(n, 0, 1)
eps <- rnorm(n, 0, 0.1)
y <- sin(2*pi*x) + eps

# Centering step
mx <- mean(x); x <- x - mx
y <- y - mean(y)

x.train <- x[1:(n/2)]
y.train <- y[1:(n/2)]
x.test <- x[(n/2 + 1):n]
y.test <- y[(n/2 + 1):n]

par(mar=c(4,4,2,2), family = 'Palatino', cex = 1.1)
plot(x.train, y.train, pch = 16, cex = 0.8, xlab = 'x', ylab = 'y', main = 'Smooth function, not noisy obs.')
lines(sort(x.train), sin(2*pi*(sort(x.train) + mx)), lwd = 2.5, col = 'seagreen3', lty = 2)

hgrid <- exp(seq(-10, 1, length.out = 50))
cv.lin <- cv.linsmooth(x.train, y.train, x.test, y.test, hgrid)
h.opt <- hgrid[which.min(cv.lin$MSE)]

xgrid <- seq(min(x.train), max(x.train), length.out = 100)
y.pred <- linsmooth(x.train, y.train, xgrid, h.opt)
lines(xgrid, y.pred, lwd = 2.5, col = 'firebrick2')


###################################
### LOCAL POLYNOMIAL REGRESSION ###
###################################

rm(list=ls())

source('SDS383D/HW3/smoothing.R')


utilities <- read.csv(file = 'SDS383D-master/data/utilities.csv')

y <- log(utilities[,2]/utilities[,3])
x <- utilities[,1]

n <- length(y)

# # Centering step
# x <- x - mean(x)
# y <- y - mean(y)


par(mar=c(4,4,2,2), family = 'Palatino', cex = 1.1)
plot(x, y, pch = 16, cex = 0.8, xlab = 'Average Temperature', ylab = 'Average daily gas bill')


# Try with different values of h and with different bandwidths
h <- 1
D <- 1
xnew <- seq(min(x), max(x), length.out = 100)
ynew <- polregr(x, y, xnew, D, h)$ynew
lines(xnew, ynew, col = 'seagreen2', lwd = 2)
ynew <- linsmooth(x, y, xnew, h)
lines(xnew, ynew, col = 'darkorchid2', lwd = 2)

# We now pick the optimal h via LOOCV
# pdf('./SDS383D/HW3/Notes/Img/prova2.pdf', width = 8.5, height = 6)
par(mar=c(4,4,2,2), family = 'Palatino', cex = 1.1)
plot(x, y, pch = 16, cex = 0.8, xlab = 'Average Temperature', ylab = 'Average daily gas bill')

D <- 1
hgrid <- seq(1, 10, length.out = 100)
LOOCV <- cv.polregr(x, y, D, hgrid)
h.opt <- hgrid[which.min(LOOCV)]

xgrid <- seq(min(x), max(x), length.out = 100)
fit <- polregr(x, y, xgrid, D, h.opt)
y.pred <- fit$ynew
H <- fit$H
lines(xgrid, y.pred, lwd = 2, col = 'seagreen3')


# Inspect the residuals of this optimal model
yhat <- polregr(x, y, x, D, h.opt)$ynew
par(mar=c(4,4,2,2), family = 'Palatino', cex = 1.1)
plot(x, y - yhat, pch = 16, cex = 0.8, xlab = 'Average Temperature', ylab = 'Residuals')


# Let us now construct confidence intervals for the regression line
fit <- polregr(x, y, x, D, h.opt)
y.pred <- fit$ynew
H <- fit$H
sigma2hat <- sum((y - as.numeric(H %*% y))^2) / (n - 2 * sum(diag(H)) + sum(diag(crossprod(H))))

par(mar=c(4,4,2,2), family = 'Palatino', cex = 1.1)
plot(x, y, pch = 16, cex = 0.8, xlab = 'Average Temperature', ylab = 'Average daily gas bill')
lines(sort(x), y.pred[order(x)], lwd = 2, col = 'firebrick2')
polygon(c(rev(sort(x)), sort(x)), c(rev(y.pred[order(x)] - sqrt(sigma2hat) * qnorm(0.975)), y.pred[order(x)] + sqrt(sigma2hat) * qnorm(0.975)), col = rgb(0.83, 0.83, 0.83, 0.5), border = rgb(0.83, 0.83, 0.83, 0.5))    

