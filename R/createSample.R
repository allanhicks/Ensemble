#' Create a multivariate normal sample
#' Checks for positive definite covariance matrix and makes a correction if not
#' @param n Sample size
#' @param means A vector of means
#' @param Sigma Variance-Covariance matrix
#' @param lower Lower bounds for the parameters (to sample from a truncate multivariate normal)
#' @param upper Upper bounds for the parameters (to sample from a truncate multivariate normal)
#' @param truncated Logical. Whether to use the truncated multivariate normal function or not. The truncated multivariate normal function may be slower.
#' @export
createSample.fn <- function(n, means, Sigma, lower=rep(-Inf, length = length(mean)), upper=rep( Inf, length = length(mean)), truncated=TRUE) {

	require(MASS)
	require(corpcor)
	#check if Sigma is positive-definite (uses corpcor package)
	#if not, remove all fixed parameters and check again

	fullSigma <- Sigma
	keep <- diag(Sigma)!=0
	Sigma <- Sigma[keep,keep]
	if(!is.positive.definite(Sigma)) {
			cat("Covariance matrix not positive definite, even after removing fixed parameters.\n")
			cat("Check your ADMB model to make sure that the Hessian was positive definite\n")
			cat("Or, use a smaller set of parameters to sample from (e.g., yrs=c(2015,2025)\n")
			stop("Exiting because can not sample\n")
	}
	fullMeans <- means
	means <- means[keep]
	fullLower <- lower
	lower <- lower[keep]
	lower[is.na(lower)] <- -Inf
	fullUpper <- upper
	upper <- upper[keep]
	upper[is.na(upper)] <- Inf

	if(truncated) {
		samp <- rtmvnorm.rejection(n, means, Sigma, lower=lower, upper=upper)
	} else {
		samp <- mvrnorm(n, means, Sigma)
	}

	#now, piece it back together because fixed values were removed.
	fullSamp <- matrix(NA,nrow=nrow(samp),ncol=ncol(fullSigma),dimnames=list(NULL,colnames(fullSigma)))
	fullSamp[,keep] <- samp
	fullSamp[,!keep] <- rep(fullMeans[!keep],each=nrow(samp))
	samp <- fullSamp

	return(samp)
}
