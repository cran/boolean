\name{fitted.booltest}
\docType{methods}
\alias{fitted.booltest}
\alias{fitted.values.booltest}
\alias{fitted-methods}
\alias{fitted.values-methods}
\alias{fitted,ANY-method}
\alias{fitted.values,ANY-method}
\alias{fitted,booltest-method}
\alias{fitted.values,booltest-method}
\title{Fitted Values (Probabilities) for Boolean Binary Response Models}
\description{
  This function implements a 'fitted' (and an identical 'fitted.values')
  S3 method for objects of class '"booltest"'. Since these pertain to a
  binary response model, the 'fitted values' are probabilities.
  If \code{method="MH"} or \code{bootstrap > 0}, then the mean for each
  estimated parameter is used to create the fitted values.
}
\section{Usage}{\describe{

\item{object = "booltest"}{fitted(object)}
}}
\keyword{methods}
