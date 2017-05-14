
setwd('~/Desktop/Semester 2/Statistical Modeling II/')
rm(list=ls())

library(lattice)
library(RColorBrewer)
library(Matrix)
library(mvtnorm)
source('./SDS383D/HW4/gaussianprocess.R')
# ========================================================
# ==== Hierarchical models: gene expression over time ====
# ========================================================

y <- read.csv(file = 'SDS383D-master/data/droslong.csv')[,-3]

y <- y[with(y, order(group, gene, replicate, time)),]

# xyplot(log2exp~time | group, data=droslong)
# 
# xyplot(log2exp~time | gene, data=droslong)
# 
# # Replicates of the same gene are almost identical: they differ by noise
# xyplot(log2exp~time | replicate, data=droslong)

cols <- rev(colorRampPalette(brewer.pal(8, "Accent"))(14))
par(mar=c(4,4,2,2), mfrow = c(1,3), family = 'Palatino')
for (j in unique(y$group))
  plot(y$time[y$group == j], y$log2exp[y$group == j], type = 'p', col = cols[as.numeric(y$gene[y$group == j])], cex = 1.2, pch = 16, ylim = c(min(y$log2exp), max(y$log2exp)))

H <- length(unique(y$group)) # number of groups



# ================================================
# Optimize hyperparameters via Marginal Likelihood
# ================================================
opt.hyperpar <- array(NA, dim = c(2*H+1, H))
colnames(opt.hyperpar) = as.character(unique(y$group))
rownames(opt.hyperpar) = c('A.h','ls.h','A.n','ls.n','A.r','ls.r','sigma')
  
# Initial values
A.h <- 2
ls.h <- 2
A.n <- 1
ls.n <- 1
A.r <- -2
ls.r <- -2
sigma <- -1
init.par <- c(A.h, ls.h, A.n, ls.n, A.r, ls.r, sigma)
init.par <- rep(0, 7)
for (h in 1:H){
  y.h = y[y$group==unique(y$group)[h],]
  opt <- optim(par = init.par, fn = marginal.negloglikelihood, y = y.h, 
               method = 'Nelder-Mead', control= list(maxit = 2000, abstol = 1e-50))
  if (opt$convergence == 0)
    opt.hyperpar[,h] = opt$par
}

hyperpar <- exp(opt.hyperpar)

# ======================================
# Predict values on a grid using the HGP
# ======================================
# Choose a grid where to predict the values
D <- 100
t.star <- seq(0.5, 12.5, length.out = D)

N <- aggregate(data=y, gene ~ group, function(x) length(unique(x)))$gene
R <- aggregate(data=y, replicate ~ gene + group, function(x) length(unique(x)))$replicate

g.h <- array(NA, dim = c(H, D))
g.h.lo <- array(NA, dim = c(H, D))
g.h.hi <- array(NA, dim = c(H, D))

g.hn <- array(NA, dim = c(sum(N), D))
g.hn.lo <- array(NA, dim = c(sum(N), D))
g.hn.hi <- array(NA, dim = c(sum(N), D))

g.hnr <- array(NA, dim = c(sum(R), D))
g.hnr.lo <- array(NA, dim = c(sum(R), D))
g.hnr.hi <- array(NA, dim = c(sum(R), D))

N <- aggregate(data=y, gene ~ group, function(x) length(unique(x)))$gene
R <- aggregate(data=y, replicate ~ gene + group, function(x) length(unique(x)))$replicate

iter <- 1
iter2 <- 1

for (h in 1:H){
  y.group <- y[y$group == unique(y$group)[h],]
  
  Sigma <- data.cov(y.group, hyperpar[,h])
  Sigmainv <- solve(Sigma)

  # Parameters for cluster-level function
  params.h <- hyperpar[1:2,h]
  
  # Off-diagonal block of covariance matrix of data and g.star
  Kh.star <- C.squaredexp(t.star, y.group$time, params.h[1], params.h[2])
  # Marginal covariance of g.star
  Kh.starstar <- C.squaredexp(t.star, t.star, params.h[1], params.h[2])
  
  # Posterior mean of h function
  g.h[h,] <- Kh.star %*% Sigmainv %*% y.group$log2exp
  
  # Posterior covariance matrix of hi.est
  hi.covmat <- Kh.starstar - Kh.star %*% Sigmainv %*% t(Kh.star)
  hi.var <- diag(hi.covmat)
  
  # 95% confidence intervals
  g.h.lo[h,] <- g.h[h,] - 1.96 * sqrt(hi.var)
  g.h.hi[h,] <- g.h[h,] + 1.96 * sqrt(hi.var)
  
  for (n in 1:N[h]){
  
    y.group.gene <- y.group[y.group$gene == unique(y.group$gene)[n],]
    Nobs <- nrow(y.group.gene)
    
    # Parameters for gene-level function
    params.n <- hyperpar[3:4,h]
    
    # Marginal covariance of g.star
    Kn.starstar <- Kh.starstar + C.squaredexp(t.star, t.star, params.n[1], params.n[2])
    
    Kn.star <- Kh.star
    # Select the indexes corresponding to the actual gene
    idx.gene <- ((n-1)*Nobs + 1):(n*Nobs)
    
    upd = C.squaredexp(t.star, y.group.gene$time, params.n[1], params.n[2])
    Kn.star[,idx.gene] <- Kn.star[,idx.gene] + upd
    
    # Posterior mean of h function
    g.hn[iter,] <- Kn.star %*% Sigmainv %*% y.group$log2exp

    # Posterior covariance matrix of hi.est
    hi.covmat <- Kn.starstar - Kn.star %*% Sigmainv %*% t(Kn.star)
    hi.var <- diag(hi.covmat)
    
    # 95% confidence intervals
    g.hn.lo[iter,] <- g.hn[iter,] - 1.96 * sqrt(hi.var)
    g.hn.hi[iter,] <- g.hn[iter,] + 1.96 * sqrt(hi.var)
    
    iter <- iter + 1
    
    for (r in 1:3){
      
      y.group.gene.rep <- y.group.gene[y.group.gene$replicate == unique(y.group.gene$replicate)[r],]
      Nobs <- nrow(y.group.gene.rep)
      
      # Parameters for gene-level function
      params.r <- hyperpar[5:6,h]
      
      # Marginal covariance of g.star
      Kr.starstar <- Kn.starstar + C.squaredexp(t.star, t.star, params.r[1], params.r[2])
      
      Kr.star <- Kn.star
      # Select the indexes corresponding to the actual gene
      idx.gene.rep <- idx.gene[((r-1)*Nobs + 1):(r*Nobs)]
      
      upd = C.squaredexp(t.star, y.group.gene.rep$time, params.r[1], params.r[2])
      Kr.star[,idx.gene.rep] <- Kr.star[,idx.gene.rep] + upd
      
      # Posterior mean of h function
      g.hnr[iter2,] <- Kn.star %*% Sigmainv %*% y.group$log2exp
      
      # Posterior covariance matrix of hi.est
      hi.covmat <- Kn.starstar - Kn.star %*% Sigmainv %*% t(Kn.star)
      hi.var <- diag(hi.covmat)
      
      # 95% confidence intervals
      g.hnr.lo[iter2,] <- g.hnr[iter2,] - 1.96 * sqrt(hi.var)
      g.hnr.hi[iter2,] <- g.hnr[iter2,] + 1.96 * sqrt(hi.var)
      
      iter2 <- iter2 + 1
    }
  }
}


par(mar = c(2,2,2,2), cex = 0.7, mfrow = c(5, 3))
iter <- 1
cols <- c('firebrick2','dodgerblue','goldenrod2')
for (h in 1:H){
  y.group <- y[y$group == unique(y$group)[h],]
  for (n in 1:N[h]){
    y.group.gene <- y.group[y.group$gene == unique(y.group$gene)[n],]
    for (r in 1:3){
      y.group.gene.rep <- y.group.gene[y.group.gene$replicate == unique(y.group.gene$replicate)[r],]
      plot(y.group.gene.rep$time, y.group.gene.rep$log2exp, type = 'p', col = 'black', bg = rep(unique(cols[as.numeric(y.group$gene)]), each = 3)[iter], main = bquote('Group'~.(h)*', Gene'~.(n)*', Rep.'~.(r)), xlab = 'Time', ylab = 'Gene Expression', cex = 1.2, pch = 21, ylim = c(min(y$log2exp), max(y$log2exp)))
      lines(t.star, g.hnr[iter,], lty = 2, lwd = 2)
      polygon(c(t.star, rev(t.star)), c(g.hnr.lo[iter,], rev(g.hnr.hi[iter,])), col = rgb(0.83,0.83,0.83,0.5), border = rgb(0.83, 0.83, 0.83, 0.5))
      iter <- iter + 1
    }
  }
}


cols <- rev(colorRampPalette(brewer.pal(8, "Accent"))(14))
par(mar=c(4,4,2,2), mfrow = c(1, 3), cex = 0.9, family = 'Palatino')
for (h in 1:H){
  y.group <- y[y$group == unique(y$group)[h],]
  plot(y.group$time, y.group$log2exp, type = 'p', col = 'black', bg = cols[as.numeric(y.group$gene)], xlab = 'Time', ylab = 'Gene Expression', main = bquote(Group~.(h)), cex = 1.2, pch = 21, ylim = c(min(y$log2exp), max(y$log2exp)))
  lines(t.star, g.h[h,], lty = 2, lwd = 2)
  polygon(c(t.star, rev(t.star)), c(g.h.lo[h,], rev(g.h.hi[h,])), col = rgb(0.83,0.83,0.83,0.5), border = rgb(0.83, 0.83, 0.83, 0.5))
}


h <- 1
par(mar=c(4,4,2,2), mfrow = c(2, 2), family = 'Palatino', cex = 0.9)
iter <- 1
cols <- c('firebrick2','dodgerblue','goldenrod2')
y.group <- y[y$group == unique(y$group)[h],]
for (n in 1:N[h]){
  y.group.gene <- y.group[y.group$gene == unique(y.group$gene)[n],]
  plot(y.group.gene$time, y.group.gene$log2exp, type = 'p', col = 'black', bg = cols[y.group.gene$replicate], xlab = 'Time', ylab = 'Gene Expression', cex = 1.2, pch = 21, ylim = c(min(y$log2exp), max(y$log2exp)))
  lines(t.star, g.hn[iter,], lty = 2, lwd = 2)
  polygon(c(t.star, rev(t.star)), c(g.hn.lo[iter,], rev(g.hn.hi[iter,])), col = rgb(0.83,0.83,0.83,0.5), border = rgb(0.83, 0.83, 0.83, 0.5))
  iter <- iter + 1
}

  