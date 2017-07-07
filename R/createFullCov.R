#' A function to create a covariance matrix from a full run (many parameters estimated) and a subset run (fewer parameters estimated)
#' It provides the correlations from the full run
#' The standard deviations (sd's) or variances are from the subset run for those estimated parameters, but the full run for the additional parameters
#' This is done so that one can use estimated sd's from a run they have more faith in than a full run, but
#'   correlations come from a full run and results in a positive definite covariance matrix (hopefully)
#' @param fullRun   SS output for the full run
#' @param subsetRun SS output for the subset run
#' @export
createFullCov.fn <- function(subsetRun, fullRun) {
	covFull <- getCov.fn(fullRun, vars=dimnames(fullRun$parameters)[[1]])
	covSub <- getCov.fn(subsetRun, vars=dimnames(subsetRun$parameters)[[1]])

	sdSub <- sqrt(diag(covSub))  #standard deviations for subset parameters
	sds <- sdSub
	#add in standard deviations for additional parameters
	sdFull <- sqrt(diag(covFull))
	if(length(sdSub) != length(sdFull)) stop("The models do not have the same parameter structure")
	addPars <- names(sdFull)[sdSub==0 & sdFull!=0]  #parameters with variance in full model only
	sds[addPars] <- sdFull[addPars]

	out <- getCov.fn(fullRun, vars=names(sdFull), sds=sds)
	return(out)
}