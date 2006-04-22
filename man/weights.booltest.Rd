\name{weights.booltest}
\docType{methods}
\alias{weights.booltest}
\alias{weights-methods}
\alias{weights,ANY-method}
\alias{weights,booltest-method}
\title{Extract weights for Boolean Binary Response Models}
\description{
  This function implements a 'weights' S3 method for objects of class '"booltest"'.
  However, weights are not yet supported in boolean models, so a vector of ones
  is returned.
}
\section{Usage}{\describe{

\item{object = "booltest"}{weights(object)}
}}
\keyword{methods}
