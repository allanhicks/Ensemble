#' Samples from each ensemble model 
#' @param models  A list of SS_ouput models
#' @param param   The name of the parameter to be combined (does not include year)
#' @param element The list element to get param from (e.g., derived_quants, parameters)
#' @param yrs     The yrs of the parameter to summarize.  If NULL, it will either include all years, or the param does not have a year.
#' @param totN    Total number of samples of combined models
#' @param wts     Weights applied to each model when combining. This is normalized, so can be any number. A single number means equal weights (default)
#' @param asList  Will keep the sample from each model as a list element.  By default it is FALSE, but will be changed to true if the variables are not the same for each model (which would occur if some are not estimated in a model)
#' @export
combineEnsemble.fn <- function(models, param="SPB", element="derived_quants", yrs=NULL, useCov=F, totN=1e6, wts=1, asList=FALSE) {

	if(length(wts)==1) wts<-rep(wts,length(models))
	wts <- wts/sum(wts)
	eachN <- totN * wts

	colNames <- switch(element,
		               derived_quants = c("Value","StdDev"),
		               parameters = c("Value","Parm_StDev"),
		               stop("'element' must be 'derived_quants' or 'parameters'.\n"))

	#get the values desired
	vals <- lapply(models,function(x,e,p,y,cn) {
								if(is.null(y)) {
									out <- getVals.fn(x[[e]], param=p, cols=cn)
								} else {
									out <- getVals.fn(x[[e]], param=paste(p,y,sep="_"),cols=cn)
								}
	        }, e=element, p=param, y=yrs, cn=colNames)

	lens <- unlist(lapply(vals,nrow))
	if(any(lens != lens[1])) {
		asList <- TRUE
	}
	#out <- combineDistn.fn(totN, vals) deprecated
	if(asList) {
		samps <- vector(length=length(models),mode="list")
	} else {
		samps <- NULL
	}
	for(i in 1:length(models)) {
		sig <- diag(nrow(vals[[i]]))
		diag(sig) <- vals[[i]][,colNames[2]]^2
		sig[is.na(sig)] <- 0
		rownames(sig) <- colnames(sig) <- rownames(vals[[i]])
		if(useCov) {
			#now fill in covariances if desired
			sig <- getCov.fn(models[[i]], rownames(vals[[i]]))
		}
		#test for positive definite and create multivariate normal sample
		tmp <- createSample.fn(eachN[i], vals[[i]][,colNames[1]],sig)
		rownames(tmp) <- paste("Model",i,1:eachN[i],sep="_")
		if(asList) {
			samps[[i]] <- tmp
			colnames(samps[[i]]) <- rownames(vals[[i]])
		} else {
			samps <- rbind(samps, tmp)			
		}
	}

	if(!asList) {
		colnames(samps) <- rownames(vals[[1]])		
	}

	return(samps)
}
