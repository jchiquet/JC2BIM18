# Formating Van't Veer dataset

rm(list=ls())
source('Functions-Print.R')

library(cancerdata)
data(VEER1); Y = 1*(VEER1@phenoData@data$class=='DM'); X = as.matrix(t(VEER1@assayData$exprs)); dim(X)
# data(VIJVER1); Y = 1*(VIJVER1@phenoData@data$class=='DM'); X = as.matrix(t(VIJVER1@assayData$exprs)); dim(X)
Remove = as.vector(which(is.na(colSums(X))))
X = X[, -Remove]; dim(X)
p.val = sapply(1:ncol(X), function(j){wilcox.test(X[, j], Y, var.equal=F)$p.value})
X = X[, which(order(p.val) <= 15)]
colnames(X) = gsub('_', '', colnames(X))
save(X, Y, file='VEER1.Rdata') 

tab = cbind(X[, 1:3], Y)[1:5, ]
colnames(tab) = c(colnames(X)[c(1:3)], 'Status')
rownames(tab) = 1:5
F_PrintTab(tab)
