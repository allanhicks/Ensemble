#' Create a multivariate normal sample
#' Checks for positive definite covariance matrix and makes a correction if not
#' @param n Sample size
#' @param means A vector of means (always entered in normal space)
#' @param Sigma Variance-Covariance matrix (always entered in normal space)
#' @param logSpace logical: take the natural log of the variable before sampling, then exponentiating back
#' @param logBiasCorr logical: apply a bias correction if the mean in normal space (rather than the median) is provided for a lognormal sample
#' @param lower Lower bounds for the parameters (to sample from a truncate multivariate normal)
#' @param upper Upper bounds for the parameters (to sample from a truncate multivariate normal)
#' @param truncated Logical. Whether to use the truncated multivariate normal function or not. The truncated multivariate normal function may be slower.
#' @param fixPD fix the PD to the nearest positive definite, if it isn't
#' @param estAlpha estimate the acceptance rate if using the truncated normal
#' 
#' @authors Allan Hicks, Gwladys Lambert
#' @export
createSample.fn <- function(n, means, Sigma, logSpace=FALSE, logBiasCorr=TRUE, lower=rep(-Inf, length = length(mean)), upper=rep( Inf, length = length(mean)), truncated=TRUE, fixPD=T, estAlpha=TRUE) {

	require(MASS)
	require(corpcor)
	#require(Matrix)
	#check if Sigma is positive-definite (uses corpcor package)
	#if not, remove all fixed parameters and check again

	fullSigma <- Sigma
	keep <- diag(Sigma)!=0
	Sigma <- Sigma[keep,keep]
	if(!is.matrix(Sigma)) {
		Sigma <- as.matrix(Sigma)
	}

	fullMeans <- means
	means <- means[keep]
	fullLower <- lower
	lower <- lower[keep]
	lower[is.na(lower)] <- -Inf
	fullUpper <- upper
	upper <- upper[keep]
	upper[is.na(upper)] <- Inf

	if(logSpace) {
		#this is an approximation
		cat("WARNING: The Variance/Covariance matrix in log space is an approximation that may not be accurate for large variances\n")
		ExEy <- matrix(means,nrow=length(means),ncol=length(means)) * t(matrix(means,nrow=length(means),ncol=length(means)))

		if(logBiasCorr) {
			logVar <- log((sqrt(diag(Sigma))/means)^2 + 1)
			means <- log(means) - 0.5*logVar
		}else {
			means <- log(means)
		}

		Sigma <- log(1+Sigma/(ExEy))   #Taylor expansion approximation: https://stats.stackexchange.com/questions/46874/covariance-of-transformed-random-variables/46930#46930

	}
	if(!is.positive.definite(Sigma)) {
		if(fixPD) {
				cat("Fixed the covariance matrix to make it positive-definite\n")
				#Sigma <- as.matrix(nearPD(Sigma, corr = FALSE, keepDiag = FALSE)$mat)
				Sigma <- make.positive.definite(Sigma)
			}else{
				cat("Covariance matrix not positive definite, even after removing fixed parameters.\n")
				cat("Check your ADMB model to make sure that the Hessian was positive definite\n")
				cat("Or, use a smaller set of parameters to sample from (e.g., yrs=c(2015,2025)\n")
				stop("Exiting because can not sample\n")				
			}
	}
	if(truncated) {
		samp <- rtmvnorm.rejection(n, means, Sigma, lower=lower, upper=upper,estimatedAlpha=estAlpha)
	} else {
		samp <- mvrnorm(n, means, Sigma)
	}

	if(logSpace) {
		#transform back to normal space
		samp <- exp(samp)
	}
	#now, piece it back together because fixed values were removed.
	fullSamp <- matrix(NA,nrow=nrow(samp),ncol=ncol(fullSigma),dimnames=list(NULL,colnames(fullSigma)))
	fullSamp[,keep] <- samp
	fullSamp[,!keep] <- rep(fullMeans[!keep],each=nrow(samp))
	samp <- fullSamp

	return(samp)
}
