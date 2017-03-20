
setwd('~/Desktop/Semester 2/Statistical Modeling II/')
rm(list=ls())

library(mvtnorm)


# ================================================
# ==== Basics on Gaussian Processes Smoothing ====
# ================================================

source('SDS383D/HW3/gaussianprocess.R')

b <- 1
tau1sq <- 1
tau2sq <- 1E-6
# Role of the parameters:
#   - b: controls how much distant points are correlated. Small values of b produce a non
#     correlated process (more noisy), large values of b produce a process with high 
#     correlation and, subsequently, smoother functions.
#   - tau1sq: scaling factor that controls the amplitude of the covariances.
#   - tau2sq: scaling factor that controls the amplitude of the variances, that is, 
#             how noisy each point is.


colors <- c('indianred2','dodgerblue2','seagreen3','goldenrod2','darkorchid2')
#pdf('./SDS383D/HW3/Notes/Img/gaussianexp6.pdf', width = 8.5, height = 6)
# Squared Exponential covariance function
par(mar=c(4,4,1,2), family = 'Palatino', cex = 1.1)
xgrid <- seq(0, 1, length.out = 100)
plot(xgrid, rgaussproc(xgrid, m = rep(0, length(xgrid)), C = C.squaredexp(xgrid, xgrid, b, tau1sq, tau2sq)), type ='l', col = colors[1], lwd = 2, ylim = c(-3, 3), xlab = 'x', ylab = 'f(x)')
for (i in 2:5){
  lines(xgrid, rgaussproc(xgrid, m = rep(0, length(xgrid)), C = C.squaredexp(xgrid, xgrid, b, tau1sq, tau2sq)), col = colors[i], lwd = 2)
}

Npred <- 250
pred <- array(NA, dim = c(Npred, length(xgrid)))
for (i in 1:Npred){
  pred[i,] <- rgaussproc(xgrid, m = rep(0, length(xgrid)), C = C.squaredexp(xgrid, xgrid, b, tau1sq, tau2sq))
}
pred.quant <- as.matrix(apply(pred, 2, quantile, prob = c(0.025, 0.5, 0.975)))
lines(xgrid, pred.quant[2,], lty = 2, lwd = 2)
polygon(c(xgrid, rev(xgrid)), c(pred.quant[3,], rev(pred.quant[1,])), col = rgb(0.83,0.83,0.83,0.5), border = rgb(0.83, 0.83, 0.83, 0.5))
dev.off()

#pdf('./SDS383D/HW3/Notes/Img/gaussianmatern6.pdf', width = 8.5, height = 6)
# Matern covariance function
par(mar=c(4,4,1,2), family = 'Palatino', cex = 1.1)
xgrid <- seq(0, 1, length.out = 100)
plot(xgrid, rgaussproc(xgrid, m = rep(0, length(xgrid)), C = C.matern52(xgrid, xgrid, b, tau1sq, tau2sq)), type ='l', lwd = 2, col = colors[1], ylim = c(-3, 3), xlab = 'x', ylab = 'f(x)')
for (i in 2:5){
  lines(xgrid, rgaussproc(xgrid, m = rep(0, length(xgrid)), C = C.matern52(xgrid, xgrid, b, tau1sq, tau2sq)), col = colors[i], lwd = 2)
}
# We notice that, in general, using the same parameters as before the Matern covariance 
# function yields to a more noisy Gaussian process.

Npred <- 250
pred <- array(NA, dim = c(Npred, length(xgrid)))
for (i in 1:Npred){
  pred[i,] <- rgaussproc(xgrid, m = rep(0, length(xgrid)), C = C.matern52(xgrid, xgrid, b, tau1sq, tau2sq))
}
pred.quant <- as.matrix(apply(pred, 2, quantile, prob = c(0.025, 0.5, 0.975)))
lines(xgrid, pred.quant[2,], lty = 2, lwd = 2)
polygon(c(xgrid, rev(xgrid)), c(pred.quant[3,], rev(pred.quant[1,])), col = rgb(0.83,0.83,0.83,0.5), border = rgb(0.83, 0.83, 0.83, 0.5))
dev.off()



# ==================================
# ==== Nonparametric Regression ====
# ==================================

rm(list = ls())

source('SDS383D/HW3/gaussianprocess.R')

utilities <- read.csv(file = 'SDS383D-master/data/utilities.csv')

y <- utilities[,2]/utilities[,3]
x <- utilities[,1]

# Centering Step               
y <- y - mean(y)

n <- length(y)
sigma2 <- 1

#pdf('./SDS383D/HW3/Notes/Img/gaussianfit.pdf', width = 8.5, height = 6)
par(mar=c(4,4,2,2), family = 'Palatino', cex = 1.1)
plot(x, y, pch = 16, cex = 0.8, xlab = 'Average Temperature', ylab = 'Average daily gas bill')

b <- 10
tau1sq <- 10
tau2sq <- 1E-6

xgrid <- seq(min(x), max(x), length.out = 100)

pred <- pred.GP(x, y, xgrid, b, tau1sq, tau2sq, sigma2, type = 'matern')
lines(xgrid, pred$pred.mean, pch = 16, lwd = 2, col = 'firebrick3')

Npred <- 500
pred.grid <- array(NA, dim = c(Npred, length(xgrid)))
for (i in 1:Npred){
  pred.grid[i,] <- rgaussproc(xgrid, m = pred$pred.mean, C = pred$pred.var)
}
pred.quant <- as.matrix(apply(pred.grid, 2, quantile, prob = c(0.025, 0.5, 0.975)))
polygon(c(xgrid, rev(xgrid)), c(pred.quant[1,], rev(pred.quant[3,])), col = rgb(0.83,0.83,0.83,0.5), lwd = 0.1)


# Pick the hyperparameters by marginal likelihood
pred <- pred.GP(x, y, x, b, tau1sq, tau2sq, sigma2, type = 'matern')
sigma2 <- as.numeric(var(y - pred$pred.mean))

tau2sq <- 1E-6

D <- 20
b.grid <- seq(90, 110, length.out = D)
tau1sq.grid <- seq(30, 40, length.out = D)
sigma2.grid <- seq(0.01, 5, length.out = D)  

marg.loglik = array(NA, dim = c(D, D, D))
for(i in 1:D){
  for(j in 1:D){
    for(k in 1:D){
      marg.loglik[i,j,k] <- marginal.likelihood(x, y, b.grid[i], tau1sq.grid[j], tau2sq, sigma2.grid[k], type = 'matern')
    }
  }
  cat('Iteration: ', i, "\n")
}
b.opt <- b.grid[which(marg.loglik == max(marg.loglik), arr.ind = T)[1]]
tau1sq.opt <- tau1sq.grid[which(marg.loglik == max(marg.loglik), arr.ind = T)[2]]
sigma2.opt <- sigma2.grid[which(marg.loglik == max(marg.loglik), arr.ind = T)[3]]

# Check if we are looking in the correct region
par(mar=c(4,4.5,2,2), family = 'Palatino')
contour(b.grid, tau1sq.grid, marg.loglik[,,which(marg.loglik == max(marg.loglik), arr.ind = T)[3]], nlevels=20, xlab = 'b', ylab = bquote(tau[1]^2), col = 'dodgerblue4', lwd = 2)

points(b.opt, tau1sq.opt, pch = 16, col = 'indianred3')

par(mar=c(4,4,2,2), family = 'Palatino', cex = 1.1)
plot(x, y, pch = 16, cex = 0.8, xlab = 'Average Temperature', ylab = 'Average daily gas bill')

pred <- pred.GP(x, y, xgrid, b.opt, tau1sq.opt, tau2sq, sigma2.opt, type = 'matern')
lines(xgrid, pred$pred.mean, pch = 16, lwd = 2, col = 'firebrick3')


Npred <- 250
pred.grid <- array(NA, dim = c(Npred, length(xgrid)))
for (i in 1:Npred){
  pred.grid[i,] <- rgaussproc(xgrid, m = pred$pred.mean, C = pred$pred.var)
}
pred.quant <- as.matrix(apply(pred.grid, 2, quantile, prob = c(0.025, 0.5, 0.975)))
polygon(c(xgrid, rev(xgrid)), c(pred.quant[1,], rev(pred.quant[3,])), col = rgb(0.83,0.83,0.83,0.5), lwd = 0.1)


# ===========================
# ==== Spatial Smoothing ====
# ===========================

rm(list = ls())

source('SDS383D/HW3/spatialsmoothing.R')
library(fields)
library(ggmap)
library(RColorBrewer)
cols <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(100))


# Load the data
weather <- read.csv(file = 'SDS383D-master/data/weather.csv')
n <- dim(weather)[1]

# Download the map and plot the data
limits <- c(min(weather$lon), min(weather$lat), max(weather$lon), max(weather$lat))
bw.map <- get_map(location = limits, maptype = 'terrain-background')
ggmap(bw.map) + 
  scale_fill_distiller(palette = 'RdYlBu') + 
  geom_point(data = weather, aes(x = lon, y = lat, fill = temperature), colour = 'gray', size = 2.5, pch = 21)

# Choice of the anisotropy parameters
center.lat <- min(weather$lat) + diff(range(weather$lat))/2
b1 <- 1/sin(center.lat) # relative to the longitude
b2 <- 1                 # relative to the latitude

# #pdf('./SDS383D/HW3/Notes/Img/pressure.pdf', width = 8.5, height = 6)
# par(mar=c(4,4,2,2), family = 'Palatino')
# quilt.plot(weather$lon, weather$lat, weather$pressure, xlab = 'Longitude', ylab = 'Latitude', main = 'Heatmap of Pressure')
# dev.off()
# #pdf('./SDS383D/HW3/Notes/Img/temperature.pdf', width = 8.5, height = 6)
# par(mar=c(4,4,2,2), family = 'Palatino')
# quilt.plot(weather$lon, weather$lat, weather$temperature, xlab = 'Longitude', ylab = 'Latitude', main = 'Heatmap of Temperature')
# dev.off()

# Rename the data
X <- cbind(weather$lon, weather$lat)
z <- weather$temperature - mean(weather$temperature)

# First guess for sigma:
# - temperature data range from -9 to 9 (use sigma2 = 1)
# - pressure data range from -600 to 700 (use sigma2 = 5000)
sigma2 <- 1

# Plot for a guess of the parameters on a coarse grid
b <- 1
tau1sq <- sigma2
tau2sq <- 1E-6

D <- 64
longrid <- seq(min(weather$lon), max(weather$lon), length.out = D)
latgrid <- seq(min(weather$lat), max(weather$lat), length.out = D)
grid <- expand.grid(longrid, latgrid)

pred <- pred.GP(X, z, grid, b, tau1sq, tau2sq, sigma2, b1, b2, type = 'matern')
par(mar=c(4,4,2,2), family = 'Palatino')
quilt.plot(grid[,1], grid[,2], pred, nx = D, ny = D)

# Estimate the variance of the model
pred <- pred.GP(X, z, X, b, tau1sq, tau2sq, sigma2, b1, b2, type = 'matern')
sigma2 <- as.numeric(var(z - pred)) # estimate the variance of the model

# Pick the hyperparameters by marginal likelihood
tau2sq <- 1E-6 # this is fixed
D <- 20
b.grid <- seq(0.5, 1.5, length.out = D)
#b.grid <- seq(1, 1.6, length.out = D)
tau1sq.grid <- seq(1, 15, length.out = D)
#tau1sq.grid <- seq(50000, 80000, length.out = D)

marg.loglik = array(NA, dim = c(D, D))
for(i in 1:D){
  for(j in 1:D){
      marg.loglik[i,j] <- marginal.likelihood(X, z, b.grid[i], tau1sq.grid[j], tau2sq, sigma2, b1, b2, type = 'matern')
  }
  cat('Iteration: ', i, "\n")
}

par(mar=c(4,4.5,2,2), family = 'Palatino')
contour(b.grid, tau1sq.grid, marg.loglik, nlevels=20, xlab = 'b', ylab = bquote(tau[1]^2), col = 'dodgerblue4', lwd = 2)

b.opt <- b.grid[which(marg.loglik == max(marg.loglik), arr.ind = T)[1]]
tau1sq.opt <- tau1sq.grid[which(marg.loglik == max(marg.loglik), arr.ind = T)[2]]
points(b.opt, tau1sq.opt, pch = 16, col = 'indianred3')

# Fit the final model (fine grid, optimal hyperparameters)
D <- 64
longrid <- seq(min(weather$lon), max(weather$lon), length.out = D)
latgrid <- seq(min(weather$lat), max(weather$lat), length.out = D)
grid <- expand.grid(longrid, latgrid)

pred <- pred.GP(X, z, grid, b.opt, tau1sq.opt, tau2sq, sigma2, b1, b2, type = 'matern')

#load(file = './SDS383D/HW3/pred_press.Rdata')
#load(file = './SDS383D/HW3/pred_temp.Rdata')

#pdf('./SDS383D/HW3/Notes/Img/tempsmooth.pdf', width = 8.5, height = 6)
par(mar=c(4,4,2,2), family = 'Palatino')
quilt.plot(grid[,1], grid[,2], pred, xlab = 'Longitude', ylab = 'Latitude', main = 'Heatmap of Temperature', nx = D, ny = D)
contour(x = longrid, y = latgrid, matrix(pred, D, D, byrow = F), nlevels = 10, add = T, lwd = 1.5, labcex = 0.7)
dev.off()

pred.values <- data.frame(cbind(grid, pred))
names(pred.values) <- c('lon','lat','pred')
#pdf('./SDS383D/HW3/Notes/Img/presssmoothmap.pdf', width = 8.5, height = 6)
ggmap(bw.map) +
  geom_tile(data = pred.values, aes(x = lon, y = lat, fill = pred), alpha = 0.7) + 
  scale_fill_distiller(palette = 'RdYlBu') +
  geom_contour(data = pred.values, aes(x = lon, y = lat, z = pred, colour = ..level..), size = 1) + 
  scale_colour_gradientn(colors = cols) + 
  geom_point(data = weather, aes(x = lon, y = lat, fill = temperature), colour = 'gray', size = 2, pch = 21)
dev.off()


# ===========================
# ==== Spatial Smoothing ====
# ===========================

rm(list = ls())

source('SDS383D/HW3/spatialsmoothing.R')
library(fields)

weather <- read.csv(file = 'SDS383D-master/data/weather.csv')


#pdf('./SDS383D/HW3/Notes/Img/pressure.pdf', width = 8.5, height = 6)
par(mar=c(4,4,2,2), family = 'Palatino')
quilt.plot(weather$lon, weather$lat, weather$pressure, xlab = 'Longitude', ylab = 'Latitude', main = 'Heatmap of Pressure')
dev.off()
#pdf('./SDS383D/HW3/Notes/Img/temperature.pdf', width = 8.5, height = 6)
par(mar=c(4,4,2,2), family = 'Palatino')
quilt.plot(weather$lon, weather$lat, weather$temperature, xlab = 'Longitude', ylab = 'Latitude', main = 'Heatmap of Temperature')
dev.off()

# Rename the data
X <- cbind(weather$lon, weather$lat)
z <- weather$pressure
w <- weather$temperature

plot(w, z, pch = 16, cex = 0.8)
fit <- lm(z ~ w)
abline(a = fit$coefficients[1], b = fit$coefficients[2], lwd = 2, col = 'indianred3')
res <- z - fit$coefficients[1] - fit$coefficients[2] * w

hist(res, col = 'gray', border = 'white', nclass = 20)

# First guess for sigma:
# - temperature data range from -9 to 9 (use sigma2 = 1)
# - pressure data range from -600 to 700 (use sigma2 = 10000)
sigma2 <- 1000

# Plot for a guess of the parameters on a coarse grid
b <- 2
tau1sq <- 10000
tau2sq <- 1E-6

D <- 64
longrid <- seq(min(weather$lon), max(weather$lon), length.out = D)
latgrid <- seq(min(weather$lat), max(weather$lat), length.out = D)
grid <- expand.grid(longrid, latgrid)

pred <- pred.GP(X, res, grid, b, tau1sq, tau2sq, sigma2)
quilt.plot(grid[,1], grid[,2], pred, nx = D, ny = D)
quilt.plot(weather$lon, weather$lat, res, nx = D, ny = D)


# Estimate the variance of the model
pred <- pred.GP(X, res, X, b, tau1sq, tau2sq, sigma2)
sigma2 <- as.numeric(var(z - pred)) # estimate the variance of the model

# Pick the hyperparameters by marginal likelihood
tau2sq <- 1E-6 # this is fixed
D <- 40
#b.grid <- seq(0.5, 1.5, length.out = D)
b.grid <- seq(1.3, 1.4, length.out = D)
#tau1sq.grid <- seq(1, 15, length.out = D)
tau1sq.grid <- seq(40000, 55000, length.out = D)

marg.loglik = array(NA, dim = c(D, D))
for(i in 1:D){
  for(j in 1:D){
    marg.loglik[i,j] <- marginal.likelihood(X, z, b.grid[i], tau1sq.grid[j], tau2sq, sigma2)
  }
  cat('Iteration: ', i, "\n")
}

par(mar=c(4,4.5,2,2), family = 'Palatino')
contour(b.grid, tau1sq.grid, marg.loglik, nlevels=20, xlab = 'b', ylab = bquote(tau[1]^2), col = 'dodgerblue4', lwd = 2)

b.opt <- b.grid[which(marg.loglik == max(marg.loglik), arr.ind = T)[1]]
tau1sq.opt <- tau1sq.grid[which(marg.loglik == max(marg.loglik), arr.ind = T)[2]]
points(b.opt, tau1sq.opt, pch = 16, col = 'indianred3')

# Fit the final model (fine grid, optimal hyperparameters)
D <- 64
longrid <- seq(min(weather$lon), max(weather$lon), length.out = D)
latgrid <- seq(min(weather$lat), max(weather$lat), length.out = D)
grid <- expand.grid(longrid, latgrid)

pred <- pred.GP(X, res, grid, b.opt, tau1sq.opt, tau2sq, sigma2)

#load(file = './SDS383D/HW3/pred_press.Rdata')
#load(file = './SDS383D/HW3/pred_temp.Rdata')

#pdf('./SDS383D/HW3/Notes/Img/tempsmooth.pdf', width = 8.5, height = 6)
par(mar=c(4,4,2,2), family = 'Palatino')
quilt.plot(grid[,1], grid[,2], pred, xlab = 'Longitude', ylab = 'Latitude', main = 'Heatmap of Temperature', nx = D, ny = D)
contour(x = longrid, y = latgrid, matrix(pred, D, D, byrow = F), nlevels = 10, add = T, lwd = 1.5, labcex = 0.7)
dev.off()



pred.mat <- matrix(pred, D, D, byrow = F)
# Define colors
nbcol = 100
color = rev(rainbow(nbcol, start = 0/6, end = 4.5/6))
zcol  = cut(pred.mat, nbcol)
# Plot 3D
open3d()
persp3d(longrid, latgrid, pred.mat, theta=50, phi=25, expand=0.75, col=color[zcol],
        ticktype="detailed", xlab="Longitude", ylab="Latitude", zlab="", axes=TRUE, alpha = 0.7)
zcol  = cut(z, nbcol)
points3d(X[,1], X[,2], z, col = color[zcol], cex = 1.5)


