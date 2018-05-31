# Bayesian inference for logistic regression

rm(list=ls())
source('Functions-BayesLogReg.R')
source('Functions-Print.R')

M = 1e5

###############################################################################
# Data
load('VEER1.Rdata')
n = nrow(X); Intercept = rep(1, n)
X = cbind(Intercept, X); d = ncol(X); 
data = list(X=X, Y=Y)

###############################################################################
# Results
load(paste0('../res/VEER1-IS-prior-M', M,'.Rdata'))
load(paste0('../res/VEER1-IS-VB-M', M,'.Rdata'))
load(paste0('../res/VEER1-IS-MLE-M', M,'.Rdata'))
load(paste0('../res/VEER1-MH-full-M', M,'.Rdata'))
load(paste0('../res/VEER1-MH-comp-M', M,'.Rdata'))

###############################################################################
# All methods : Confidence intervals 
CI.prior; CI.VB; CI.MLE; CI.MH.full

###############################################################################
# MH : Autocorrelation check
par(mfrow=c(d, d), mex=.3, pch=20)
invisible(sapply(1:d, function(j){sapply(1:d, function(k){
   plot(beta.sample[1:(beta.nb-1), j], beta.sample[2:beta.nb, k], main='', xlab='', ylab='')
})}))

###############################################################################
# MH : Inference
# Posterior histograms
pdf('../figs/posterior-hist.pdf')
par(mfrow=c(2, 2))
sapply(1:4, function(j){hist(beta.sample[, j], breaks=sqrt(M/2), xlab='', ylab='', main=colnames(X)[j])
   abline(v=0, col=2, lwd=2)})
dev.off()

# Posterior density
pdf('../figs/posterior-density.pdf')
par(mfrow=c(2, 2))
sapply(1:4, function(j){plot(density(beta.sample[, j]), xlab='', ylab='', main=colnames(X)[j])
   abline(v=0, col=2, lwd=2)})
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
plot(delta.post.density, xlab='', ylab='', main='delta'); abline(v=0, col=2, lwd=2)
dev.off()
post.delta = matrix(0, 1, 4)
colnames(post.delta) = c('post.mean', 'post.mode', 'lower.CI', 'upper.CI')
rownames(post.delta) = 'delta'
post.delta[1, ] = c(mean(delta.sample), delta.post.density$x[which.max(delta.post.density$y)], 
                    as.vector(WeightedCI(IS=list(beta.sample=delta.sample, weight.sample=rep(1, beta.nb)/beta.nb))))
F_PrintTab(post.delta)

###############################################################################
# MH : Model comparison
par(mfrow=c(1, 1))
pY = sapply(1:d, function(j){
   logpY = logLikelihood(data=list(X=as.matrix(X[, 1:j]), Y=Y), 
                         H.sample=list(beta=beta.sample.comp[[j]]))
   logpY.mean = mean(logpY)
   logpY = logpY - logpY.mean
   exp(logpY.mean) * mean(exp(logpY))
})
plot(0:(d-1), pY, type='b', xlab='Model')
