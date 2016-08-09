#' Gets the mean and standard deviation of parameters from SS output
#' @param x     The matrix of parameters or derived parameters from SS_output
#' @param param The name of the parameter
#' @param cols  The column names for the mean and standard deviation, respectively
getVals.fn <- function(x, param, cols=c("Value","StdDev")) {
	if(is.null(cols)) {cols <- 1:ncol(x)}
	if(is.null(param)) {param<-rownames(x)}

	if(length(param)>1) {
		out <- x[param,cols]
	}
	if(length(param)==1) {
		out <- x[grep(param,rownames(x),fixed=T),cols]
	}
	return(out)
}