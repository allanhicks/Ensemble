#' Get the covariance matrix from a SS model for a specific set of parameters
#' Uses the SS version of the correlation matrix
#' @param ssModel A SS model loaded using SS_output from r4ss
#' @param vars    Names of the variables for the covariance matrix
#' @export
getCov.fn <- function(ssModel, vars) {
	#cat("Getting Covariance matrix. This may take a while. There are",length(vars),"rows.\nRow:")
	x <- ssModel$CoVar
	x <- x[x$label.i %in% vars,]
	x$label.j[x$label.j=="_"] <- x$label.i[x$label.j=="_"]
	x <- x[x$label.j %in% vars,]

	covar <- matrix(0,length(vars),length(vars),dimnames=list(vars,vars))

	for(i in vars) {
		indx <- x$label.i == i
		covar[i,x$label.j[indx]] <- x$corr[indx]
		covar[x$label.j[indx],i] <- x$corr[indx]
	}

	#saves sd's and convert to covariance matrix
	sds <- diag(covar)
	diag(covar) <- 1
	covar <- cor2cov.fn(covar,sds)
	rownames(covar) <- vars

	return(covar)
}
