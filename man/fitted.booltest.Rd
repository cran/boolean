\name{fitted.booltest}
\alias{fitted.booltest}
\alias{fitted.values.booltest}
\alias{fitted}
\alias{fitted.values}
\title{Fitted Values (Probabilities) for Boolean Binary Response Models}
\description{
  This function implements a 'fitted' (and an identical 'fitted.values')
  S3 method for objects of class 'booltest'. Since these pertain to a
  binary response model, the 'fitted values' are probabilities.
  If \code{method="MH"} or \code{bootstrap > 0} in the call to \code{\link{boolean}}, 
  then the mean for each estimated parameter is used to create the fitted values.
}
\usage{
## S3 method for class 'booltest'
\method{fitted}{booltest}(object, ...)
\method{fitted.values}{booltest}(object, ...)
}

\arguments{
  \item{object}{An object of class 'booltest'}
  \item{...}{Anything passed will be ignored}
}
\value{
  A numeric vector of probabilities.
}
\references{
  Braumoeller, Bear F. (2003) "Causal Complexity and the Study
  of Politics." \emph{Political Analysis} 11(3): 209-233.
}
\author{
Bear F. Braumoeller, Harvard University, \email{braumoeller.1@polisci.osu.edu},\cr
Ben Goodrich, Harvard University, \email{goodrich@fas.harvard.edu}, and\cr
Jacob Kline, Harvard University, \email{jkline@fas.harvard.edu}
}

\keyword{methods}

