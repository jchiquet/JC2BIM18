###############################################################################
# Logistic regression
x = seq(-3, 3, by=.01)
a = -2; b = 3
pdf('../figs/logistic-curve.pdf')
plot(x, plogis(a + b*x), main='', xlab='', ylab='', lwd=4, col=2, type='l')
dev.off()
