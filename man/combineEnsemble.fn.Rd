% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/combineEnsemble.R
\name{combineEnsemble.fn}
\alias{combineEnsemble.fn}
\title{Samples from each ensemble model}
\usage{
combineEnsemble.fn(models, param = "SPB", element = "derived_quants",
  yrs = NULL, useCov = F, totN = 1e+06, wts = 1, asList = FALSE)
}
\arguments{
\item{models}{A list of SS_ouput models}

\item{param}{The name of the parameter to be combined (does not include year)}

\item{element}{The list element to get param from (e.g., derived_quants, parameters)}

\item{yrs}{The yrs of the parameter to summarize.  If NULL, it will either include all years, or the param does not have a year.}

\item{totN}{Total number of samples of combined models}

\item{wts}{Weights applied to each model when combining. This is normalized, so can be any number. A single number means equal weights (default)}

\item{asList}{Will keep the sample from each model as a list element.  By default it is FALSE, but will be changed to true if the variables are not the same for each model (which would occur if some are not estimated in a model)}
}
\description{
Samples from each ensemble model
}

