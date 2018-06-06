# Bayesian inference for logistic regression

rm(list=ls())
source('Functions-BayesLogReg.R')
source('Functions-Print.R')
library(plotrix)
library(ks)
library(mvtnorm)

M = 1e5
shift = .5 # .1, .5, 1. Ex = .5

###############################################################################
# Data
load('VEER1.Rdata')
n = nrow(X); Intercept = rep(1, n)
X = cbind(Intercept, X); d = ncol(X); 
data = list(X=X, Y=Y)

###############################################################################
# Results
load('../res/VEER1-prior.Rdata')
theta.grid = seq(-20, 20, by=.01)

###############################################################################
# All methods : Confidence intervals 
load(paste0('../res/VEER1-IS-prior-M', M,'-shift', shift, '.Rdata'))
load(paste0('../res/VEER1-IS-VB-M', M, '-shift', shift, '.Rdata'))
load(paste0('../res/VEER1-IS-MLE-M', M, '-shift', shift, '.Rdata'))
load(paste0('../res/VEER1-MH-full-M', M, '-shift', shift, '.Rdata'))
CI.prior; CI.VB; CI.MLE; CI.MH.full
MLE.var = summary(MLE)$cov.scaled

# IS from approximate posterior
pdf('../figs/ISproposal-density.pdf')
par(mfrow=c(2, 2), mex=.5)
sapply(1:4, function(j){
   beta.grid = seq(min(beta.sample[, j]), max(beta.sample[, j]), length.out=500)
   plot(beta.grid, dnorm(beta.grid, mean=VB$Parms.vb$mean[j], sd=sqrt(VB$Parms.vb$Variance[j, j])), col=2, lwd=2, type='l', 
        xlab='', ylab='', main=colnames(X)[j])
   # VB.weight = IS.VB$weight.sample / sum(IS.VB$weight.sample)
   # weighted.hist(IS.VB$beta.sample[, j], w=VB.weight, add=T, border=2)
   lines(beta.grid, dnorm(beta.grid, mean=beta.prior.mean[j], sd=sqrt(beta.prior.var[j, j])), col=4, lwd=2)
   lines(beta.grid, dnorm(beta.grid, mean=MLE$coefficients[j], sd=sqrt(MLE.var[j, j])), col=3, lwd=2)
   lines(density(beta.sample[, j]), lwd=2)
})
dev.off()

# IS from approximate posterior: 2D density
pdf('../figs/ISproposal-density2D.pdf')
prob = .01; cont = c(10, 50, 90)
par(mfrow=c(2, 2), mex=.5)
for (j in 2:3){
   for(k in (j+1):4){
      KDE = kde(beta.sample[, c(j, k)])
      plot(KDE, main='', xlab=colnames(X)[j], ylab=colnames(X)[j], lwd=2, xlim=quantile(beta.sample[, j], probs=c(prob, 1-prob)), 
           ylim=quantile(beta.sample[, k], probs=c(prob, 1-prob)), cont=cont)
      abline(v=mean(beta.sample[, j]), h=mean(beta.sample[, k]), lty=2, lwd=1)
      plotmixt(mus=VB$Parms.vb$mean[c(j, k)], Sigmas=VB$Parms.vb$Variance[c(j, k), c(j, k)], props=1, 
               col=2, lwd=2, add=T, cont=cont)
      abline(v=VB$Parms.vb$mean[j], h=VB$Parms.vb$mean[k], col=2, lty=2, lwd=1)
   }
}
dev.off()

###############################################################################
# MH : Sampling path

# Acceptance rate
shift.list = c(.1, .5, 1)
rate = sapply(shift.list, function(s){
   load(paste0('../res/VEER1-MH-full-M', M, '-shift', s, '.Rdata'))
   MH.full$accept / MH.full$M.MH
})
Res = rbind(round(rate, 3))
colnames(Res) = shift.list
rownames(Res) = 'acceptance rate'
F_PrintTab(Res)

# par(mfrow=c(4, 1), mex=.25)
# for(j in 1:4){
#    plot(MH.full$beta.path[, j], xlab='', ylab='', main='', type='l')
# }

pdf(paste0('../figs/MH-path-10shift', round(10*shift), '.pdf'))
par(mfrow=c(4, 1), mex=.25)
for(j in 1:4){
   plot(beta.sample[, j], xlab='', ylab='', main='', type='l')
}
dev.off()
 
# # Autocorrelation check
# par(mfrow=c(2, 2), mex=.3, pch=20)
# invisible(sapply(1:4, function(j){plot(MH.full$beta.path[1:(nrow(MH.full$beta.path)-1), j], 
#                                        MH.full$beta.path[2:nrow(MH.full$beta.path), j], main='', xlab='', ylab='')
# }))

# Autocorrelation check
pdf(paste0('../figs/MH-autocorrelation-10shift', round(10*shift), '.pdf'))
par(mfrow=c(2, 2), mex=.3, pch=20)
invisible(sapply(1:4, function(j){plot(beta.sample[1:(beta.nb-1), j], 
                                       beta.sample[2:beta.nb, j], main='', xlab='', ylab='')
}))
dev.off()

###############################################################################
# MH full  : Inference
load(paste0('../res/VEER1-MH-full-M', M, '-shift', shift, '.Rdata'))
beta.nb = nrow(beta.sample)
# Posterior histograms
pdf('../figs/posterior-hist.pdf')
par(mfrow=c(2, 2))
phi = list()
sapply(1:4, function(j){
   H = hist(beta.sample[, j], breaks=sqrt(beta.nb), xlab='', ylab='', main=colnames(X)[j])
   phi[[j]] <<- dnorm(theta.grid, mean=beta.prior.mean[j], sd=sqrt(beta.prior.var[j, j]))
   lines(theta.grid, beta.nb*mean(diff(H$breaks))*phi[[j]], lwd=2, col=4)
   abline(v=0, lty=2, lwd=2)
   })
dev.off()

# Posterior density
pdf('../figs/posterior-density.pdf')
par(mfrow=c(2, 2), mex=.5)
sapply(1:4, function(j){
   plot(density(beta.sample[, j]), xlab='', ylab='', main=colnames(X)[j], lwd=2)
   lines(theta.grid, phi[[j]], lwd=2, col=4)
   abline(v=0, lty=2, lwd=2)
   })
dev.off()

# Posterior estimates
post.estim = matrix(0, 4, 4)
colnames(post.estim) = c('post.mean', 'post.mode', 'lower.CI', 'upper.CI')
rownames(post.estim) = colnames(X)[1:4]
sapply(1:4, function(j){
   posterior = density(beta.sample[, j])
   post.estim[j, ] <<- c(mean(beta.sample[, j]), posterior$x[which.max(posterior$y)], CI.MH.full[, j])
})
F_PrintTab(post.estim)

# Posterior inference for a contrast
delta.sample = matrix(beta.sample[, 2] - beta.sample[, 3], beta.nb, 1)
pdf('../figs/posterior-density-delta.pdf')
par(mfrow=c(1, 1))
delta.post.density = density(delta.sample)
plot(delta.post.density, xlab='', ylab='', main='delta', lwd=2); 
abline(v=0, lty=2, lwd=2)
dev.off()
post.delta = matrix(0, 1, 4)
colnames(post.delta) = c('post.mean', 'post.mode', 'lower.CI', 'upper.CI')
rownames(post.delta) = 'delta'
post.delta[1, ] = c(mean(delta.sample), delta.post.density$x[which.max(delta.post.density$y)], 
                    as.vector(WeightedCI(IS=list(beta.sample=delta.sample, weight.sample=rep(1, beta.nb)/beta.nb))))
F_PrintTab(post.delta)

###############################################################################
# MH : Model comparison and averaging
load(paste0('../res/VEER1-MH-comp-M', M, '-shift', shift, '.Rdata'))
par(mfrow=c(1, 1), pch=20)
i0 = 1; x0 = X[i0, ]
pY = Esp.prob0 = Var.prob0 = rep(0, d)
for (k in 1:d){
   logpY = logLikelihood(data=list(X=as.matrix(X[, 1:k]), Y=Y), 
                         H.sample=list(beta=beta.sample.comp[[k]]))
   logpY.mean = mean(logpY)
   logpY = logpY - logpY.mean
   pY[k] = exp(logpY.mean) * mean(exp(logpY))
   if(k > 1){
      prob0 = plogis(beta.sample.comp[[k]]%*%x0[1:k])
   }else{
      prob0 = plogis(beta.sample.comp[[k]]*x0[1])
   }
   Esp.prob0[k] = mean(prob0); Var.prob0[k] = var(prob0)
}
pM.Y = pY / sum(pY)
plot(0:(d-1), pM.Y, type='b', xlab='Model')

as.vector(x0)

Sd.prob0 = sqrt(Var.prob0)
pM.Y = pY / sum(pY)
prob0.avg = cbind(round(log(pY), 2), round(pM.Y, 4), round(Esp.prob0, 3), round(Sd.prob0, 3))
colnames(prob0.avg) = c('logpY|M', 'pM|Y', 'Esp.prob0', 'Sd.prob0')
rownames(prob0.avg) = paste0('$M_{', 0:(d-1), '}$')
F_PrintTab(prob0.avg)

logpY.marg = log(sum(pY)/length(pY))
Esp.prob0.avg = pM.Y %*% Esp.prob0
Var.prob0.avg = pM.Y %*% Var.prob0 + (pM.Y %*% Esp.prob0^2) - Esp.prob0.avg^2
Sd.prob0.avg = sqrt(Var.prob0.avg)
paste0(round(logpY.marg, 2), ' & & ', round(Esp.prob0.avg, 3), ' & ', round(Sd.prob0.avg, 3))

###############################################################################
# MH : Data combination
pdf('../figs/posterior-combine.pdf')
par(mfrow=c(2, 2))
sapply(1:4, function(j){
   load(paste0('../res/VEER1-MH-full-M', M, '-shift', shift, '.Rdata'))
   plot(density(beta.sample[, j]), xlab='', ylab='', main=colnames(X)[j], lwd=2)
   load(paste0('../res/VEER1-MH-half-M', M, '-shift', shift, '.Rdata'))
   lines(density(beta.sample[, j]), col=2, lwd=2)
   lines(theta.grid, phi[[j]], lwd=2, col=4)
   abline(v=0, lty=2, lwd=2)
})
dev.off()

