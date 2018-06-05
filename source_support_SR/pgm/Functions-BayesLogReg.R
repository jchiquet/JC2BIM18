#############################################################################
library(mvtnorm)

##########################################################
logPrior <- function(H.sample, HyperParms){
   beta.sample=H.sample$beta
   mean.Prior = HyperParms$mean; 
   Variance.Prior = HyperParms$Variance;
   p = length(mean.Prior)
   if(p==1){beta.sample=matrix(beta.sample,ncol=1)}
   logPrior.sample = dmvnorm(beta.sample,mean.Prior, Variance.Prior,log=TRUE)
   return(logPrior.sample)
}

#############################################################################
# Log-likelihood for the logistic regression model
logLikelihood <- function(data, H.sample){
   Y <- data$Y
   X <- data$X
   p = dim(X)[2]
   beta = H.sample$beta
   if(p==1){beta=matrix(beta,ncol=1)}
   M = dim(beta)[1]
   n = length(Y)
   YY = matrix(Y,nrow=n,ncol=M,byrow=FALSE)
   logl = colSums(dbinom(YY, 1, plogis(X%*%t(beta)), log=T))
   return(logl)
}

########################################################################"
logApproxPost <- function(H.sample, HyperParms.ApproxPost){
   
   beta.sample=H.sample$beta
   
   mean.ApproxPost = HyperParms.ApproxPost$mean; 
   Variance.ApproxPost = HyperParms.ApproxPost$Variance;
   
   p = length(mean.ApproxPost)
   if(p==1){beta.sample=matrix(beta.sample,ncol=1)}
   logApproxPost.sample = dmvnorm(beta.sample,mean.ApproxPost, Variance.ApproxPost,log=TRUE)
   return(logApproxPost.sample)
}

###############################################################################
# Importance sampling
SampleIS <- function(beta.prop.mean, beta.prop.var, beta.prior.mean, beta.prior.var, M, data=data){
   # beta.prop.mean = beta.VB.mean; beta.prop.var = beta.VB.var; 
   # beta.prop.mean = beta.prior.mean; beta.prop.var = beta.prior.var; 
   beta.sample = rmvnorm(M, mean=beta.prop.mean, sigma=beta.prop.var)
   loglik.sample = logLikelihood(data, H.sample=list(beta=beta.sample))
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
# Variational Bayes algorithm for logistic regression (Jaakkola & Jordan, 2000)
VB <- function(data, HyperParms){

   X = data$X # Covariates
   Y = data$Y # Response
   n = length(Y)
   
   mean.prior = HyperParms$mean # prior mean vector for the coefficients
   Variance.prior = HyperParms$Variance # prior variance matrix for the coefficients
   XY = t(X) %*% (Y-.5)
   Precision.prior = solve(Variance.prior)
   
   # Intialization
   ResGLM = glm(Y~X-1,family = binomial(link="logit"))
   mean.vb = ResGLM$coefficients
   Variance.vb = summary(ResGLM)$cov.scaled
   Precision.vb = solve(Variance.vb)
      
   # 'EM' inference   
   tol = 1e-6; diff = 2*tol; iter = 0
   while (diff > tol){
      iter = iter+1
      Variance.vb = solve(Precision.vb)
      
      # 'E' step
      xsi = sqrt(diag(X %*% (mean.vb%*%t(mean.vb) + Variance.vb) %*% t(X)))
      lambda = tanh(xsi/2) / 4 / xsi
      
      # 'M' step
      Precision.tmp = Precision.prior + 2 * t(X) %*% diag(lambda) %*% X
      Variance.tmp = solve(Precision.tmp)
      mean.tmp = Variance.tmp %*% (Precision.prior %*% mean.prior + XY)
      
      # Test & update
      diff = max(max(abs(mean.vb-mean.tmp)), max(abs(Variance.vb-Variance.tmp)))
      mean.vb = mean.tmp; Variance.vb = Variance.tmp; Precision.vb = Precision.tmp
      
      #       print(t(mean.vb)); print(Variance.vb); 
      #       cat(diff, '')
   }
   cat('(', iter, ') ', sep='')
   logpY.vb = -sum(log(1+exp(-xsi))) - .5*sum(xsi) + sum(lambda*xsi^2) - 
      .5*t(mean.prior)%*%Precision.prior%*%mean.prior + 
      .5*t(mean.vb)%*%Precision.vb%*%mean.vb + 
      .5*(log(det(Variance.vb)/det(Variance.prior)))
   Parms.vb = list(mean = as.vector(mean.vb), Variance = Variance.vb, Precision = Precision.vb, logpY=logpY.vb)
   return(list(Parms.vb = Parms.vb, iter=iter,log.marg=logpY.vb))
}

###############################################################################
# Metropolis-Hastings algorithm for logistic regression
MH <- function(data, beta.prior.mean, beta.prior.var, beta.shift.var, M=1e3, coef.MH=1.2){
   M.MH = round(coef.MH*M); d = ncol(data$X)
   beta.cur = matrix(0, 1, d) #rmvnorm(1, mean=beta.prior.mean, sigma=beta.prior.var)
   logpY.cur = logLikelihood(data, H.sample=list(beta=beta.cur))
   logprior.cur = logPrior(H.sample=list(beta=beta.cur), HyperParms=list(mean=beta.prior.mean, Variance=beta.prior.var))
   beta.path = matrix(0, M.MH, d); beta.path[1, ] = beta.cur; accept = 0
   par(mfrow=c(d, 1), mex=.3)
   for(m in 1:M.MH){
      if(m%%round(sqrt(M.MH))==0){cat(accept, '/', m, ' ', sep='')}
      beta.prop = beta.cur + rmvnorm(1, sigma=beta.shift.var)
      logpY.prop = logLikelihood(data, H.sample=list(beta=beta.prop))
      logprior.prop = logPrior(H.sample=list(beta=beta.prop), HyperParms=list(mean=beta.prior.mean, Variance=beta.prior.var))
      logratio = logpY.prop + logprior.prop - logpY.cur - logprior.cur
      if(logratio < -100){ratio = 0}else if(logratio > +100){ratio = 1}else{ratio = min(c(1,exp(logratio)))}
      if(runif(1) < ratio){
         accept = accept+1; 
         beta.path[m, ] = beta.prop
         beta.cur = beta.prop; logpY.cur = logpY.prop; logprior.cur = logprior.prop
      }else{
         beta.path[m, ] = beta.cur
      }
   }
   return(list(beta.path=beta.path, accept=accept, M.MH=M.MH))
}
