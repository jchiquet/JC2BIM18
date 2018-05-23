# Bayesian inference for logistic regression

rm(list=ls())
library(mvtnorm); library(Hmisc)
source('Functions-LogReg.R')

###############################################################################
# Data
# Simulations
n = 100; d = 5; 
beta = (1:d)-(d+1)/2
X = matrix(rnorm(n*d), n, d); prob = plogis(X%*%beta)
Y = rbinom(n, 1, prob)

# Van't Veer
library(cancerdata)
data(VEER1); Y = 1*(VEER1@phenoData@data$class=='DM'); X = as.matrix(t(VEER1@assayData$exprs)); dim(X)
# data(VIJVER1); Y = 1*(VIJVER1@phenoData@data$class=='DM'); X = as.matrix(t(VIJVER1@assayData$exprs)); dim(X)
Remove = as.vector(which(is.na(colSums(X))))
X = X[, -Remove]; dim(X)
p.val = sapply(1:ncol(X), function(j){wilcox.test(X[, j], Y, var.equal=F)$p.value})
X = X[, which(order(p.val) <= 15)]

# Dimensions
n = nrow(X); d = ncol(X); data = list(X=X, Y=Y)

# Prior & MC parms
M = 1e4
beta.prior.mean = rep(0, d); beta.prior.var = 1e2*diag(d)
beta.shift.var = 1e-1*diag(d)

###############################################################################
# Functions
SampleIS <- function(beta.prop.mean, beta.prop.var, beta.prior.mean, beta.prior.var, M, data=data){
   # beta.prop.mean = beta.VB.mean; beta.prop.var = beta.VB.var; 
   # beta.prop.mean = beta.prior.mean; beta.prop.var = beta.prior.var; 
   beta.sample = rmvnorm(M, mean=beta.prop.mean, sigma=beta.prop.var)
   loglik.sample = likelihood(data, H.sample=list(beta=beta.sample))
   logprior.sample = logPrior(H.sample=list(beta=beta.sample), HyperParms=list(mean=beta.prior.mean, Variance=beta.prior.var))
   logprop.sample = logApproxPost(H.sample=list(beta=beta.sample), HyperParms=list(mean=beta.prop.mean, Variance=beta.prop.var))
   logweight.sample = loglik.sample+logprior.sample-logprop.sample
   logweight.sample[which(logweight.sample < -100)] = -100; logweight.sample[which(logweight.sample > +100)] = +100
   weight.sample = exp(logweight.sample - mean(logweight.sample)); 
   weight.sample = weight.sample / sum(weight.sample)
   return(list(beta.sample=beta.sample, weight.sample=weight.sample, ESS=1/sum(weight.sample^2)))
}
WeightedCI <- function(IS, probs=c(.025, .975)){
   d = ncol(IS$beta.sample); prob.nb = length(probs); CI = matrix(0, prob.nb, d)
   sapply(1:d, function(j){
      beta.order = order(IS$beta.sample[, j])
      weight.sort = IS$weight.sample[beta.order]
      weight.cum = cumsum(weight.sort)
      beta.sort = IS$beta.sample[beta.order, j]
      sapply(1:prob.nb, function(k){CI[k, j] <<- beta.sort[which.min(abs(weight.cum-probs[k]))]})
   })
   return(CI)
}

###############################################################################
# Importance sampling 
par(mfrow=c(3, 1))
# IS from prior
IS.prior = SampleIS(beta.prior.mean, beta.prior.var, beta.prior.mean, beta.prior.var, M, data=data)
hist(log10(IS.prior$weight.sample), breaks=sqrt(M), main=round(IS.prior$ESS, 2))
CI.prior = WeightedCI(IS.prior); print(CI.prior)

# IS from VB
VB = VB(data=list(Y=Y, X=X), HyperParms=list(mean.prior=beta.prior.mean, Variance.prior=beta.prior.var))
beta.VB.mean = VB$Parms.vb$mean; beta.VB.var = VB$Parms.vb$Variance
IS.VB = SampleIS(beta.VB.mean, beta.VB.var, beta.prior.mean, beta.prior.var, M, data=data)
hist(log10(IS.VB$weight.sample), breaks=sqrt(M), main=round(IS.VB$ESS, 2))
CI.VB = WeightedCI(IS.VB); print(CI.VB)

# IS from MLE
MLE = glm(Y ~ -1 + X, family='binomial'); beta.MLE.mean = MLE$coefficients; beta.MLE.var = vcov(MLE)
IS.MLE = SampleIS(beta.MLE.mean, beta.MLE.var, beta.prior.mean, beta.prior.var, M, data=data)
hist(log10(IS.MLE$weight.sample), breaks=sqrt(M), main=round(IS.MLE$ESS, 2))
CI.MLE = WeightedCI(IS.MLE); print(CI.MLE)

###############################################################################
# Metropolis-Hastings

MH = round(1.2*M)
beta.tmp = matrix(0, 1, d) #rmvnorm(1, mean=beta.prior.mean, sigma=beta.prior.var)
logpY.tmp = likelihood(data, H.sample=list(beta=beta.tmp))
logprior.tmp = logPrior(H.sample=list(beta=beta.tmp), HyperParms=list(mean=beta.prior.mean, Variance=beta.prior.var))
beta.path = matrix(0, MH, d); iter = 0; sample = 0
par(mfrow=c(d, 1), mex=.3)
while(sample < MH){
   iter = iter+1
   if(iter%%round(sqrt(MH))==0){cat(sample, '/', iter, ' ', sep='')}
   beta.prop = beta.tmp + rmvnorm(1, sigma=beta.shift.var)
   logpY.prop = likelihood(data, H.sample=list(beta=beta.prop))
   logprior.prop = logPrior(H.sample=list(beta=beta.prop), HyperParms=list(mean=beta.prior.mean, Variance=beta.prior.var))
   logratio = logpY.prop + logprior.prop - logpY.tmp - logprior.tmp
   if(logratio < -100){ratio = 0}else if(logratio > +100){ratio = 1}else{ratio = min(c(1,exp(logratio)))}
   if(runif(1) < ratio){
      sample = sample+1; beta.path[sample, ] = beta.prop
      beta.tmp = beta.prop; logpY.tmp = logpY.prop; logprior.tmp = logprior.prop
   }
   # if((iter%%round(sqrt(MH))==0) & (sample>0)){
   #    for (j in 1:d){plot(beta.path[1:sample, j], main='', xlab='', ylab='', type='l')}
   # }
}
cat('acceptance =', sample/iter, '\n')
for (j in 1:d){plot(beta.path[, j], main='', xlab='', ylab='', type='l')}

# Remove burn-in period
beta.path = beta.path[(MH-M+1):MH, ]
beta.sample = beta.path[seq(1, M, by=10), ]
beta.nb = nrow(beta.sample)
par(mfrow=c(d, d), mex=.3, pch=20)
invisible(sapply(1:d, function(j){sapply(1:d, function(k){
   plot(beta.sample[1:(beta.nb-1), j], beta.sample[2:beta.nb, k], main='', xlab='', ylab='')
   })}))
CI.MH = WeightedCI(IS=list(beta.sample=beta.sample, weight.sample=rep(1, beta.nb)/beta.nb))

CI.prior; CI.VB; CI.MLE; CI.MH
