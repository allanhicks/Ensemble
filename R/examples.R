


if(F) {



# covariance of log variables
library(MASS)
mu <- c(100,20)
Sigma <- matrix(c(10,3,3,2),2,2)
Sigma
xy <- mvrnorm(n = 1000, mu, Sigma)
var(xy)
ExEy <- matrix(mu,nrow=length(mu),ncol=length(mu)) * t(matrix(mu,nrow=length(mu),ncol=length(mu)))
logSigma <- log(1+Sigma/(ExEy))   #Taylor expansion approximation: https://stats.stackexchange.com/questions/46874/covariance-of-transformed-random-variables/46930#46930
logxy <- mvrnorm(n=10000, log(mu), logSigma)
apply(exp(logxy),2,summary)
var(exp(logxy))

library(corpcor)
mu <- c(100,20,2)
Sigma <- matrix(c(10,5,2,5,3,2,2,2,3),3,3)
Sigma
xy <- mvrnorm(n = 1000, mu, Sigma)
var(xy)
ExEy <- matrix(mu,nrow=length(mu),ncol=length(mu)) * t(matrix(mu,nrow=length(mu),ncol=length(mu)))
logSigma <- log(1+Sigma/(ExEy))   #Taylor expansion approximation: https://stats.stackexchange.com/questions/46874/covariance-of-transformed-random-variables/46930#46930
logSigmaPD <- make.positive.definite(logSigma)
logxy <- mvrnorm(n=100000, log(mu), logSigmaPD)
apply(exp(logxy),2,summary)
var(exp(logxy))


source("C:\\GitHub\\Ensemble\\R\\createSample.R")
source("C:\\GitHub\\Ensemble\\R\\Utilities.R")
xxx <- createSample.fn(n=10000, means=mu, Sigma=Sigma, logSpace=FALSE, logBiasCorr=TRUE, lower=rep(-Inf, length = length(mean)), upper=rep( Inf, length = length(mean)), truncated=TRUE, fixPD=T, estAlpha=TRUE)
apply(xxx,2,summary)
round(var(xxx),1)

xxxlog <- createSample.fn(n=100000, means=mu, Sigma=Sigma, logSpace=TRUE, logBiasCorr=TRUE, lower=rep(-Inf, length = length(mean)), upper=rep( Inf, length = length(mean)), truncated=FALSE, fixPD=T, estAlpha=FALSE)
apply(xxxlog,2,summary)
round(var(xxxlog))


}
