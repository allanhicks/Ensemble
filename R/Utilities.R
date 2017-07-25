#' To Read the admb .cor file
#' @param corFile   The location and name of the .cor file
#' @param preNames  The number of columns before the correlations are reported
#' @export
readCor.fn <- function(corFile, preNames=4) {
	xx <- readLines(corFile)
	xx <- xx[-1]  #remove the first line (determinant Hessian)
	nVars <- as.numeric(read.table(textConnection(xx[length(xx)]))[1])
	colNames <- (as.character(read.table(textConnection(xx[1]),header=F,stringsAsFactors=F)))
	xx <- xx[-1]  #remove column names line
	
	out <- as.data.frame(matrix(NA,ncol=length(colNames),nrow=nVars))
	names(out) <- colNames
	for(i in 1:nVars) {
		tmp <- read.table(textConnection(xx[i]))
		index <- as.numeric(tmp[1])
		nCells <- preNames+index
		out[i,1:nCells] <- as.numeric(tmp[1,])
		out[i,2] <- as.character(tmp[1,2])
	}	
	return(out)
}

#' convert a correlation matrix to a covariance matrix
#' @param corMat The correlation matrix
#' @param sds    A vector of the standard deviations associated with the diagonal of the correlation matrix
cor2cov.fn <- function(corMat,sds) {
	b <- sds %*% t(sds)
	out <- b * corMat
	return(out)
}

#' Creates a matrix X (n x k) with random realizations 
#' From a truncated multivariate normal distribution with k dimensions
#' Via rejection sampling from a multivariate normal distribution with the condition
#' lower <= Dx <= upper
#'
#' Adapted from the R package rtmvnorm version 1.4-9 
#' Written by Stefan Wilhelm (wilhelm@financial.com)
#' Posted on https://github.com/cran/tmvtnorm/blob/master/R/rtmvnorm.R
#'
#' @param n Number of realizations
#' @param mean Mean value vector (k x 1) of the normal distribution
#' @param sigma Covariance matrix (k x k) of the normal distribution
#' @param lower Lower bounds vector (k x 1) with lower <= mean <= upper
#' @param upper Upper bounds vector (k x 1) with lower <= x <= upper
rtmvnorm.rejection <- function(n, 
  						mean = rep(0, nrow(sigma)), 
  						sigma = diag(length(mean)), 
  						lower = rep(-Inf, length = length(mean)), 
  						upper = rep( Inf, length = length(mean))) {
	
	# Adapted from the R package rtmvnorm version 1.4-9 
	# Written by Stefan Wilhelm (wilhelm@financial.com)
	# Posted on https://github.com/cran/tmvtnorm/blob/master/R/rtmvnorm.R

	require(mvtnorm)
  	# No check of input parameters, checks are done in createSample.fn()
  
	# k = Dimension
	k <- length(mean)
  
	# Output matrix (n x k)
	Y <- matrix(NA, n, k)
  
	# Number of sample
  	numSamples <- n
  
  	# Number of accepted samples kept
  	numAcceptedSamplesTotal <- 0
  
  	# Acceptance rate (alpha) from the multivariate normal distribution
    alpha <- pmvnorm(lower=lower, upper=upper, mean=mean, sigma=sigma)
    if (alpha <= 0.01) warning(sprintf("Acceptance rate is very low (%s) and rejection sampling becomes inefficient. Consider using Gibbs sampling.", alpha))
    estimatedAlpha <- TRUE

    # Defer calculation of alpha. Assume for now that all samples will be accepted.
	# alpha <- 1
	# estimatedAlpha <- FALSE
  
    # Repeat from Multivariate NV and see how many samples remain after truncation
    while(numSamples > 0)
    {
		# Generate N / alpha samples from a multivariate normal distribution: If alpha is too low, rejection sampling becomes inefficient and N / alpha too large. Then create only N
		nproposals <- ifelse (numSamples/alpha > 1000000, numSamples, ceiling(max(numSamples/alpha,10)))
		X <- rmvnorm(nproposals, mean=mean, sigma=sigma)

		# Determine the proportion of samples after truncation
		# ind=apply(X, 1, function(x) all(x >= lower & x<=upper))
		ind <- logical(nproposals)
		for (i in 1:nproposals)
		{
	  	    ind[i] <- all(X[i,] >= lower & X[i,] <= upper)
		} 

		# Number of samples accepted in this run
		numAcceptedSamples <- sum(ind)
		    
		# If nothing has been accepted, then restart loop
		if (length(numAcceptedSamples) == 0 || numAcceptedSamples == 0) next   #halt the processing of the current iteration and advance the looping index

		if (!estimatedAlpha) {
			alpha <- numAcceptedSamples / nproposals
			if (alpha <= 0.01) warning(sprintf("Acceptance rate is very low (%s) and rejection sampling becomes inefficient. Consider using Gibbs sampling.", alpha))
		}

		#cat("nProposals=",nproposals,"numSamplesAccepted=",numAcceptedSamples," numSamplesNeeded = ",numSamples,"\n")
		numKeptSamples <- min(numAcceptedSamples, numSamples) #don't more than max desired
		Y[(numAcceptedSamplesTotal+1):(numAcceptedSamplesTotal+numKeptSamples),] <- X[which(ind)[1:numKeptSamples],]
		    
		# Number of total samples accepted
		numAcceptedSamplesTotal <- numAcceptedSamplesTotal + numAcceptedSamples

		# Number of remaining samples
		numSamples <- numSamples - numAcceptedSamples 
  	}
  	Y
}

#' A simple way to sample from multiple models with no correlation between parameters
#' @param totSamps Total number of samples across all models
#' @param means    The means of the parameter for each model
#' @param stds     The standard deviations of the parameter for each model
#' @param wts      The weight for each model.  A single value is equal weighting. Standardized to 1.
combineDistn.fn <- function(totSamps,means,stds,wts=1) {
	if(length(means)!=length(stds)) stop("means and stds not equal length\n")
	if(length(wts)==1) wts<-rep(wts,length(means))
	wts <- wts/sum(wts)

	out <- NULL
	for(i in 1:length(means)) {
		out <- c(out, rnorm(totSamps*wts[i], means, stds))
	}
	return(out)
}
