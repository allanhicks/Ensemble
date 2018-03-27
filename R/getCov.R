#' Get the covariance matrix from a SS model for a specific set of parameters
#' Uses the SS version of the correlation matrix
#' @param ssModel A SS model loaded using SS_output from r4ss
#' @param vars    Names of the variables for the covariance matrix
#' @params sds    Standard deviations to put on the diagonal. Make sure to enter standard deviations and not variances. This is to allow the user to use correlations from a run, but impose different sd's (diagonal).
#' @export
getCov.fn <- function(ssModel, vars, sds=NULL) {
	#cat("Getting Covariance matrix. This may take a while. There are",length(vars),"rows.\nRow:")

	x <- ssModel$CoVar

	if(length(grep("SSB",vars)) > 0) {
		#change SPB to SSB to match newer r4ss versions and names of derived_quants
		x$label.i <- gsub("SPB","SSB",x$label.i)
		x$label.j <- gsub("SPB","SSB",x$label.j)		
	}

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
	if(is.null(sds)) {
		sds <- diag(covar)
	} else {
		if(length(sds) != length(diag(covar))) {stop("Supplied SDs not equal to number of parameters in covar matrix")}
	}
	diag(covar) <- 1
	covar <- cor2cov.fn(covar,sds)
	rownames(covar) <- vars

	return(covar)
}
