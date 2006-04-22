\name{coef.booltest}
\docType{methods}
\alias{coef.booltest}
\alias{coef-methods}
\alias{coef,ANY-method}
\alias{coef,booltest-method}
\title{Coefficients for Boolean Binary Response Models}
\description{
  This function implements a 'coef' S3 method for objects of class '"booltest"'.
  If the model was estimated by maximum likelihood, the maximum likelihood
  estimates are returned. If the model was estimated by Metropolis-Hastings or
  bootstrapped, the mean over the iterations for each estimated parameter is
  returned with a warning.
}
\section{Usage}{\describe{

\item{object = "booltest"}{coef(object)}
}}
\keyword{methods}
