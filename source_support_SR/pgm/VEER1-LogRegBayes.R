# Bayesian inference for logistic regression

rm(list=ls())
library(mvtnorm); library(Hmisc)
source('Functions-BayesLogReg.R')
source('Functions-Print.R')

###############################################################################
# Data
load('VEER1.Rdata')
n = nrow(X); Intercept = rep(1, n)
X = cbind(Intercept, X); d = ncol(X); 
data = list(X=X, Y=Y)

# Prior & MC parms
M = 1e4
beta.prior.mean = rep(0, d); beta.prior.var = 1e2*diag(d)
shift = 1
beta.shift.var = shift*diag(d)
save(beta.prior.mean, beta.prior.var, file='../res/VEER1-prior.Rdata')

###############################################################################
# Importance sampling 
par(mfrow=c(3, 1))

# IS from prior
if(file.exists(paste0('../res/VEER1-IS-prior-M', M, '-shift', shift, '.Rdata'))==F){
   IS.prior = SampleIS(beta.prior.mean, beta.prior.var, beta.prior.mean, beta.prior.var, M, data=data)
   hist(log10(IS.prior$weight.sample), breaks=sqrt(M), main=round(IS.prior$ESS, 2))
   CI.prior = WeightedCI(IS.prior); print(CI.prior)
   save(IS.prior, CI.prior, file=paste0('../res/VEER1-IS-prior-M', M, '-shift', shift, '.Rdata'))
}

# IS from VB
if(file.exists(paste0('../res/VEER1-IS-VB-M', M, '-shift', shift, '.Rdata'))==F){
   VB = VB(data=list(Y=Y, X=X), HyperParms=list(mean.prior=beta.prior.mean, Variance.prior=beta.prior.var))
   beta.VB.mean = VB$Parms.vb$mean; beta.VB.var = VB$Parms.vb$Variance
   IS.VB = SampleIS(beta.VB.mean, beta.VB.var, beta.prior.mean, beta.prior.var, M, data=data)
   hist(log10(IS.VB$weight.sample), breaks=sqrt(M), main=round(IS.VB$ESS, 2))
   CI.VB = WeightedCI(IS.VB); print(CI.VB)
   save(VB, IS.VB, CI.VB, file=paste0('../res/VEER1-IS-VB-M', M, '-shift', shift, '.Rdata'))
}

# IS from MLE
if(file.exists(paste0('../res/VEER1-IS-MLE-M', M, '-shift', shift, '.Rdata'))==F){
   MLE = glm(Y ~ -1 + X, family='binomial'); beta.MLE.mean = MLE$coefficients; beta.MLE.var = vcov(MLE)
   IS.MLE = SampleIS(beta.MLE.mean, beta.MLE.var, beta.prior.mean, beta.prior.var, M, data=data)
   hist(log10(IS.MLE$weight.sample), breaks=sqrt(M), main=round(IS.MLE$ESS, 2))
   CI.MLE = WeightedCI(IS.MLE); print(CI.MLE)
   save(MLE, IS.MLE, CI.MLE, file=paste0('../res/VEER1-IS-MLE-M', M, '-shift', shift, '.Rdata'))
}

# ###############################################################################
# # Metropolis-Hastings : for code display
# p = d
# mu.prior = rep(0, p); Sigma.prior = 100*diag(p); Sigma.shift = .5*diag(p)
# theta.sample = matrix(0, M, p)
# logprior.cur = dmvnorm(theta.sample[1, ], mean=mu.prior, sigma=Sigma.prior, log=T)
# prob.cur = plogis(X%*%theta.sample[1, ])
# loglik.cur = sum(dbinom(Y, 1, prob.cur, log=T))
# for (m in 2:M){
#    theta.tmp = rmvnorm(1, mean=theta.sample[m-1, ], sigma=Sigma.shift)[1, ]
#    logprior.tmp = dmvnorm(theta.tmp, mean=mu.prior, sigma=Sigma.prior, log=T)
#    prob.tmp = plogis(X%*%theta.tmp)
#    loglik.tmp = sum(dbinom(Y, 1, prob.tmp, log=T))
#    alpha = exp(logprior.tmp + loglik.tmp - logprior.cur - loglik.cur)
#    if(runif(1) < alpha){
#       theta.sample[m, ] = theta.tmp; logprior.cur = logprior.tmp; loglik.cur = loglik.tmp
#    }else{theta.sample[m, ] = theta.sample[m-1, ]}
# }

###############################################################################
# Metropolis-Hastings : full model
if(file.exists(paste0('../res/VEER1-MH-full-M', M, '-shift', shift, '.Rdata'))==F){
   MH.full = MH(data, beta.prior.mean, beta.prior.var, beta.shift.var, M=M)
   cat('acceptance =', MH.full$sample/MH.full$iter, '\n')
   for (j in 1:d){plot(MH.full$beta.path[, j], main='', xlab='', ylab='', type='l')}
   # Remove burn-in period
   beta.path = MH.full$beta.path[(MH.full$M.MH-M+1):MH.full$M.MH, ]
   beta.sample = beta.path[seq(1, M, by=10), ]
   beta.nb = nrow(beta.sample)
   CI.MH.full = WeightedCI(IS=list(beta.sample=beta.sample, weight.sample=rep(1, beta.nb)/beta.nb))
   save(MH.full, beta.sample, beta.nb, CI.MH.full, file=paste0('../res/VEER1-MH-full-M', M, '-shift', shift, '.Rdata'))
}

###############################################################################
# Metropolis-Hastings : full model, half data
if(file.exists(paste0('../res/VEER1-MH-half-M', M, '-shift', shift, '.Rdata'))==F){
   sample1 = sort(sample(1:n, round(n/2)))
   sample2 = (1:n)[-sample1]
   MH.half = MH(data=list(X=X[sample1, ], Y=Y[sample1]), beta.prior.mean, beta.prior.var, beta.shift.var, M=M)
   cat('acceptance =', MH.half$sample/MH.half$iter, '\n')
   for (j in 1:d){plot(MH.half$beta.path[, j], main='', xlab='', ylab='', type='l')}
   # Remove burn-in period
   beta.path = MH.half$beta.path[(MH.half$M.MH-M+1):MH.half$M.MH, ]
   beta.sample = beta.path[seq(1, M, by=10), ]
   beta.nb = nrow(beta.sample)
   CI.MH.half = WeightedCI(IS=list(beta.sample=beta.sample, weight.sample=rep(1, beta.nb)/beta.nb))
   save(MH.half, beta.sample, beta.nb, CI.MH.half, sample1, file=paste0('../res/VEER1-MH-half-M', M, '-shift', shift, '.Rdata'))
}

###############################################################################
# Metropolis-Hastings : model comparison
if(file.exists(paste0('../res/VEER1-MH-comp-M', M, '-shift', shift, '.Rdata'))==F){
   MH.comp = beta.sample.comp = beta.nb.comp = list()
   for (j in 1:d){
      cat('Model ', j, ':\n')
      MH.comp[[j]] = MH(data=list(X=as.matrix(X[, 1:j]), Y=Y), 
                        beta.prior.mean=as.vector(beta.prior.mean[1:j]), 
                        beta.prior.var=as.matrix(beta.prior.var[1:j, 1:j]), 
                        beta.shift.var=as.matrix(beta.shift.var[1:j, 1:j]), M=M)
      cat('acceptance =', MH.comp[[j]]$sample/MH.comp[[j]]$iter, '\n')
      # Remove burn-in period
      beta.path = as.matrix(MH.comp[[j]]$beta.path[(MH.comp[[j]]$M.MH-M+1):MH.comp[[j]]$M.MH, ])
      beta.sample.comp[[j]] = beta.path[seq(1, M, by=10), ]
      beta.nb.comp[[j]] = nrow(beta.sample.comp[[j]])
   }
   save(MH.comp, beta.sample.comp, beta.nb.comp, file=paste0('../res/VEER1-MH-comp-M', M, '-shift', shift, '.Rdata'))
}

###############################################################################
# ABC
epsilon = .1; iter = sample = 0; beta.sample = matrix(0, M, d)
while(sample < M){
   iter = iter + 1
   beta.tmp = rmvnorm(1, mean=beta.prior.mean, sigma=beta.prior.var)[1, ]
   prob.tmp = plogis(X %*% beta.tmp)
   Y.tmp = rbinom(n, 1, prob.tmp)
   dist.tmp = sum(abs(Y.tmp - Y))/n
   if(dist.tmp < epsilon){
      sample = sample + 1
      beta.sample[sample, ] = beta.tmp
   }
   if (iter %% sqrt(M)==0){cat(sample, '/', iter, '')}
}
cat(sample, '/', iter, '\n')
