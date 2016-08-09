#' Create a multivariate normal sample
#' Checks for positive definite covariance matrix and makes a correction if not
#' @param n Sample size
#' @param means A vector of means
#' @param Sigma Variance-Covariance matrix
#' @export
createSample.fn <- function(n, means, Sigma) {

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
			stop("Exiting because can not sample\n")
	}
	fullMeans <- means
	means <- means[keep]

	samp <- mvrnorm(n, means, Sigma)

	#now, piece it back together because fixed values were removed.
	fullSamp <- matrix(NA,nrow=nrow(samp),ncol=ncol(fullSigma),dimnames=list(NULL,colnames(fullSigma)))
	fullSamp[,keep] <- samp
	fullSamp[,!keep] <- rep(fullMeans[!keep],each=nrow(samp))
	samp <- fullSamp

	return(samp)
}
