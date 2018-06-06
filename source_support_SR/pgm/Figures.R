rm(list=ls())
library(plotrix)
library(ks)
source('Functions-Print.R')

###############################################################################
# Logistic regression
x = seq(-8, 8, by=.01)
a = 0; b = 1
pdf('../figs/logistic-curve.pdf')
par(mfrow=c(1, 1), mex=0.3)
plot(x, plogis(a + b*x), main='', xlab='', ylab='', lwd=4, col=2, type='l')
dev.off()

###############################################################################
# Beta binomial
u = seq(0, 1, by=.001)
p = 1/3; 
a.list = c(1, 10, 1, 5); b.list = c(1, 1, 10, 5); list.nb = length(a.list)
# Fixed n
n = 10; y = round(n*p)
for (i in 1:list.nb){
   a = a.list[i]; b = b.list[i]
   pdf(paste0('../figs/beta-binomial-a', a, '-b', b, '-n', n, '-p', round(100*p), '-zoomed.pdf'))
   par(mfrow=c(1, 1), mex=0.3)
   plot(u, dbeta(u, a+y, b+n-y), main='', xlab='', ylab='', lwd=4, 
        col=1, type='l', ylim=c(0, 5))
   lines(u, dbeta(u, a, b), lwd=4, col=4)
   abline(v=y/n, lty=2, lwd=4)
   dev.off()
}
# Increasing n
n.list = c(10, 100, 1000)
for (n in n.list){
   y = round(n*p)
   for (i in 1:list.nb){
      a = a.list[i]; b = b.list[i]
      pdf(paste0('../figs/beta-binomial-a', a, '-b', b, '-n', n, '-p', round(100*p), '.pdf'))
      par(mfrow=c(1, 1), mex=0.3)
      plot(u, dbeta(u, a+y, b+n-y), main='', xlab='', ylab='', lwd=4, 
           col=1, type='l', ylim=c(0, 30))
      lines(u, dbeta(u, a, b), lwd=4, col=4)
      abline(v=y/n, lty=2, lwd=4)
      dev.off()
   }
}

# ###############################################################################
# # Monte Carlo : Esp(exp(N(0, 10)))
# mu=0; sigma2 = 10; sigma = sqrt(sigma2)
# E = exp(mu+sigma2/2)
# # mean(exp(rnorm(M, mean=mu, sd=sigma)))
# M.list = 10^(3:6); M.nb = length(M.list); B = 1e2
# E.hat = matrix(0, B, M.nb)
# sapply(1:M.nb, function(i){
#    M = M.list[i]
#    sapply(1:B, function(b){
#       theta.sample = rnorm(M, mean=0, sd=sigma)
#       E.hat[b, i] <<- mean(exp(theta.sample))
#    })
# })
# 
# Res = cbind(round(c(apply(E.hat, 2, mean), E), 2), c(round(apply(E.hat, 2, sd), 2), '--'))
# colnames(Res) = c('mean', 'sd')
# rownames(Res) = c(M.list, 'truth')
# F_PrintTab(Res)
# 
# theta.grid = seq(-3*sigma, 3*sigma, by=.01)
# set.seed(3); theta.sample = rnorm(10, sd=sigma)
# phi = dnorm(theta.grid, sd=sigma)
# g = exp(theta.grid)
# pdf('../figs/EspLogNorm-MC.pdf')
# par(mfrow=c(1, 1), mex=.3, pch=20)
# plot(theta.grid, phi, lwd=4, col=4, type='l', xlab='', ylab='', ylim=c(0, 3))
# lines(theta.grid, g, col=2, lwd=4)
# points(theta.sample, rep(0, length(theta.sample)), cex=2)
# dev.off() 
# 
# pdf('../figs/EspLogNorm-MC-log.pdf')
# par(mfrow=c(1, 1), mex=.3, pch=20)
# plot(theta.grid, phi, lwd=4, col=4, type='l', xlab='', ylab='', log='y', ylim=c(1e-3, 1e3))
# lines(theta.grid, g, col=2, lwd=4)
# points(theta.sample, rep(1e-3, length(theta.sample)), cex=2)
# dev.off()

###############################################################################
# Importance sampling: Sampling from a beta
M = 1e3
y.grid = seq(0, 1, by=.01)
a.target = 6; b.target = 12; p.target = dbeta(y.grid, a.target, b.target)
a.prop = 1; b.prop = 1; p.prop = dbeta(y.grid, a.prop, b.prop)
y.prop = rbeta(M, a.prop, b.prop)
H = hist(y.prop, breaks=sqrt(M), xlim=c(0, 1))
bin.width = mean(diff(H$breaks))

pdf('../figs/ImportanceSampling-1.pdf')
par(mfrow=c(1, 1), mex=0.3)
plot(y.grid, bin.width*M*p.target, lwd=4, col=1, type='l', main='', xlab='', ylab='', xlim=c(0, 1), 
     ylim=c(0, 1.2*M*bin.width*max(p.target)))
dev.off()

pdf('../figs/ImportanceSampling-2.pdf')
par(mfrow=c(1, 1), mex=0.3)
plot(y.grid, bin.width*M*p.target, lwd=4, col=1, type='l', main='', xlab='', ylab='', xlim=c(0, 1), 
     ylim=c(0, 1.2*M*bin.width*max(p.target)))
lines(y.grid, bin.width*M*p.prop, lwd=4, col=4)
dev.off()

pdf('../figs/ImportanceSampling-3.pdf')
par(mfrow=c(1, 1), mex=0.3)
plot(y.grid, bin.width*M*p.target, lwd=4, col=1, type='l', main='', xlab='', ylab='', xlim=c(0, 1), 
     ylim=c(0, 1.2*M*bin.width*max(p.target)))
lines(y.grid, bin.width*M*p.prop, lwd=4, col=4)
hist(y.prop, breaks=sqrt(M), add=T, lwd=2, border=4)
dev.off()

W = dbeta(y.prop, a.target, b.target) / dbeta(y.prop, a.prop, b.prop)

pdf('../figs/ImportanceSampling-4.pdf')
par(mfrow=c(1, 1), mex=0.3)
plot(y.grid, bin.width*M*p.target, lwd=4, col=1, type='l', main='', xlab='', ylab='', xlim=c(0, 1), 
     ylim=c(0, 1.2*M*bin.width*max(p.target)))
hist(y.prop, breaks=sqrt(M), add=T, lwd=2, border=4)
weighted.hist(y.prop, w=W, breaks=H$breaks, add=T, border=1, lwd=2, ylab='', xaxis=F)
lines(y.grid, bin.width*M*p.prop, lwd=4, col=4)
lines(y.grid, bin.width*M*p.target, lwd=4, col=1)
dev.off()

###############################################################################
# Importance sampling: Influence of the proposal
M = 1e3; y.grid = seq(0, 1, by=.01); seed = 3
a.target = 6; b.target = 12; p.target = dbeta(y.grid, a.target, b.target)
a.list = c(1, 10, 1, 5); b.list = c(1, 1, 10, 5)
prop.nb = length(a.list)
par(mfrow=c(2, 2)); set.seed(seed)
for (i in 1:prop.nb){
   a.prop = a.list[i]; b.prop = b.list[i]
   y.prop = rbeta(M, a.prop, b.prop)
   p.prop = dbeta(y.grid, a.prop, b.prop)
   W = dbeta(y.prop, a.target, b.target) / dbeta(y.prop, a.prop, b.prop)
   ESS = mean(W)^2 / mean(W^2)   
   cat(a.target, b.target, a.prop, b.prop, ESS, '\n')
   H = hist(y.prop, breaks=sqrt(M), xlim=c(0, 1))
   bin.width = mean(diff(H$breaks))
   
   pdf(paste0('../figs/ImportanceSampling-Beta-a0', a.target, '-b0', b.target, '-a', a.prop, '-b', b.prop, '-seed', seed, '.pdf'))
   par(mfrow=c(1, 1), mex=0.3)
   hist(y.prop, breaks=sqrt(M), xlim=c(0, 1), 
        ylim=c(0, 1.5*bin.width*M*max(p.target)), ylab='', xlab='', 
        main=paste('ESS =', 100*round(ESS, 4), '%'), border=4, cex.main=2)
   weighted.hist(y.prop, w=W, add=T, breaks=H$breaks, xaxis=F)
   lines(y.grid, bin.width*M*p.prop, lwd=4, col=4)
   lines(y.grid, bin.width*M*p.target, lwd=4, col=1)
   dev.off()
}

###############################################################################
# Joint, marginal, conditional
mu = c(0, 0); Sigma = matrix(c(1, 1.25, 1.25, 2), 2, 2); cont = c(10, 50, 90)
x.grid = seq(-4, 4, by=.01); grid.lim = c(min(x.grid), max(x.grid))

# Joint
pdf('../figs/Joint.pdf')
par(mfrow=c(1 ,1), mex=.3)
plotmixt(mus=mu, Sigmas=Sigma, props=1, col=2, lwd=3, cont=cont, 
         xlim=grid.lim, ylim=grid.lim, main='', xlab='', ylab='', xaxt='n', yaxt='n')
abline(v=mu[1], h=mu[2], lty=2, col=2, lwd=2)
dev.off()

# Marginal
pdf('../figs/JointMarg.pdf')
par(mfrow=c(1 ,1), mex=.3)
plotmixt(mus=mu, Sigmas=Sigma, props=1, col=2, lwd=3, cont=cont, 
         xlim=grid.lim, ylim=grid.lim, main='', xlab='', ylab='', xaxt='n', yaxt='n')
abline(v=mu[1], h=mu[2], lty=2, col=2, lwd=2)
phi.marg = dnorm(x.grid, mean=mu[1], sd=sqrt(Sigma[1, 1])); 
coef.phi = 3
lines(x.grid, min(x.grid)+coef.phi*phi.marg, col=4, lwd=3); 
abline(v=mu[1], lty=2, col=4, lwd=2)
dev.off()

# Ref
pdf('../figs/JointMargRef.pdf')
par(mfrow=c(1 ,1), mex=.3)
plotmixt(mus=mu, Sigmas=Sigma, props=1, col=2, lwd=3, cont=cont, 
         xlim=grid.lim, ylim=grid.lim, main='', xlab='', ylab='', xaxt='n', yaxt='n')
abline(v=mu[1], h=mu[2], lty=2, col=2, lwd=2)
lines(x.grid, min(x.grid)+coef.phi*phi.marg, col=4, lwd=3); 
abline(v=mu[1], lty=2, col=4, lwd=2)
# Ref
y = 1.5; 
abline(h=y, lty=4, lwd=2)
dev.off()

# Conditional
pdf('../figs/JointMargCond.pdf')
par(mfrow=c(1 ,1), mex=.3)
plotmixt(mus=mu, Sigmas=Sigma, props=1, col=2, lwd=3, cont=cont, 
         xlim=grid.lim, ylim=grid.lim, main='', xlab='', ylab='', xaxt='n', yaxt='n')
abline(v=mu[1], h=mu[2], lty=2, col=2, lwd=2)
lines(x.grid, min(x.grid)+coef.phi*phi.marg, col=4, lwd=3); 
abline(v=mu[1], lty=2, col=4, lwd=2)
abline(h=y, lty=4, lwd=2)
# Conditional
mu.cond = mu[1] + (y-mu[2])*Sigma[1, 2]/Sigma[2, 2]
sigma2.cond = Sigma[1, 1] - Sigma[1, 2]*Sigma[2, 1]/Sigma[2, 2]
phi.cond = dnorm(x.grid, mean=mu.cond, sd=sqrt(sigma2.cond)); 
lines(x.grid, y+coef.phi*phi.cond, col=1, lwd=3); 
abline(v=mu.cond, lty=2, col=1, lwd=2)
dev.off()