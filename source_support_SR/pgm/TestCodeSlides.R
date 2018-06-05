setwd('/home/robin/Bureau/JC2BIM18/source_support_SR/pgm')
library(mvtnorm)

load('VEER1.Rdata')
n = nrow(X); Intercept = rep(1, n)
X = cbind(Intercept, X); p = d = ncol(X); 
B = 1e5

mu.prior = rep(0, p); Sigma.prior = 100*diag(p); Sigma.shift = .5*diag(p)
theta.sample = matrix(0, B, p);
theta.cur = theta.sample[1, ]
logprior.cur = dmvnorm(theta.cur, mean=mu.prior, sigma=Sigma.prior, log=T)
prob.cur = plogis(X%*%theta.cur)
loglik.cur = sum(dbinom(Y, 1, prob.cur, log=T))
for (b in 2:B){
   theta.prop = rmvnorm(1, mean=theta.sample[b-1, ], sigma=Sigma.shift)[1, ]
   logprior.prop = dmvnorm(theta.prop, mean=mu.prior, sigma=Sigma.prior, log=T)
   prob.prop = plogis(X%*%theta.prop)
   loglik.prop = sum(dbinom(Y, 1, prob.prop, log=T))
   alpha = exp(logprior.prop + loglik.prop - logprior.cur - loglik.cur)
   if(runif(1) < alpha){
      theta.sample[b, ] = theta.cur = theta.prop
      logprior.cur = logprior.prop; loglik.cur = loglik.prop
   }else{
      theta.sample[b, ] = theta.sample[b-1, ]
   }
}

jump = which(theta.sample[1:(B-1), 1]!=theta.sample[2:B, 1])
accept = length(jump)
max(diff(jump))
accept/B
par(mfrow=c(4, 4), mex=.3)
for (j in 1:p){
   plot(theta.sample[, j], xlab='', ylab='', main='', type='l')
}
