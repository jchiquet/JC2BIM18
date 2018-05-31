# Bayesian inference for logistic regression

rm(list=ls())
library(mvtnorm); library(Hmisc)
source('Functions-Print.R')
source('Functions-LogReg.R')

###############################################################################
# Data
load('VEER1.Rdata')
n = nrow(X); d = ncol(X); data = list(X=X, Y=Y)

# GLM
GLM = glm(Y ~ X, family=binomial)
summary(GLM)
F_PrintTab(summary(GLM)$coef)
