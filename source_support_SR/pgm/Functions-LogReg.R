#############################################################################
library(mvtnorm)


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

#############################################################################
# Likelihood for the logistic regression model
likelihood <- function(data, H.sample){
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

#-----------------------------------------------------------------
# Importance sampling with a given init
#-----------------------------------------------------------------
rApproxPost <- function(M, HyperParms.ApproxPost){
  # Sampling theta
  mean.ApproxPost = HyperParms.ApproxPost$mean; 
  Variance.ApproxPost = HyperParms.ApproxPost$Variance;
  
  # Sampling gamma and pi 
  beta.sample = rmvnorm(M,mean.ApproxPost,Variance.ApproxPost)
  return(list(beta=beta.sample))
}

#-----------------------------------------------------------------
rPrior <- function(M, HyperParms){
  # Sampling beta
  mean.Prior = HyperParms$mean; 
  Variance.Prior = HyperParms$Variance;
  beta.sample = rmvnorm(M,mean.Prior,Variance.Prior)
  return(list(beta=beta.sample))
}


#-----------------------------------------------------------------
# Gibbs sampler
#-----------------------------------------------------------------
MCMC.Kernel <- function(data, H , rho, HyperParms.ApproxPost, HyperParms, op.save=FALSE, op.print=FALSE, Parms.MCMC=list(),op.SMC.classic){
  
  ### 
  # HyperParms : Param?tres de la loi a priori
  # HyperParms.ApproxPost : param?tres de la loi a posteriori approch?e
  # rho : param?tres de pond?ration post, post approch?e
  # B :  Number of it?rations
  
  
  Y = data$Y; 
  X = data$X; 
  
  
  
  mean.ApproxPost = HyperParms.ApproxPost$mean; 
  Variance.ApproxPost = HyperParms.ApproxPost$Variance;
  mean.Prior = HyperParms$mean; 
  Variance.Prior = HyperParms$Variance;
  
  tau = Parms.MCMC$tau
  Sigma =Parms.MCMC$Sigma
  B = Parms.MCMC$B
  
  p=dim(X)[2]
  seqbeta = matrix(0,B,p)
  
  beta=H$beta; 
  L1 = (1-rho)*dmvnorm(beta,mean.ApproxPost,Variance.ApproxPost,log=TRUE)
  L2 = rho*dmvnorm(beta,mean.Prior,Variance.Prior,log=TRUE)
  L3 = rho*sum(dbinom(Y, 1, plogis(X%*%matrix(beta,ncol=1)), log=T))
  LL = L1+L2+L3
  if(LL==-Inf){LL = -1e4}
  
  
  for (b in 1:B){
    if((op.print==TRUE)&(b%%5000==0)){print(c('MCMC Kernel : iteration',b))}
    
    beta_c = beta + rmvnorm(1,rep(0,p),tau*sample(c(1/10,1,10),1)*Sigma)
    #### proba acceptation
    L1c = (1-rho)*dmvnorm(beta_c,mean.ApproxPost,Variance.ApproxPost,log=TRUE)
    L2c = rho*dmvnorm(beta_c,mean.Prior,Variance.Prior,log=TRUE)
    L3c = rho*sum(dbinom(Y, 1, plogis(X%*%matrix(beta_c,ncol=1)), log=T))
    
    LLc = L1c+L2c+L3c;
    if(LLc==-Inf){LLc = -1e4}
    
    if(log(runif(1))<(LLc-LL)){
      beta=beta_c
      LL = LLc
    }
    
    if(op.save==TRUE){seqbeta[b,]=beta}
  }
  if(op.save==FALSE){return(list(beta=beta))}
  if(op.save==TRUE){return(list(seqbeta=seqbeta))}
}







