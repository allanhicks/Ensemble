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
