\name{weights.booltest}
\alias{weights.booltest}
\alias{weights}
\title{Extract weights for Boolean Binary Response Models}
\description{
  This function implements a 'weights' S3 method for objects of class '"booltest"'.
  However, weights are not yet supported in boolean models, so a vector of ones
  is returned.
}
\usage{
## S3 method for class 'booltest'
\method{weights}{booltest}(object, ...)
}

\arguments{
  \item{object}{An object of class 'booltest'}
  \item{...}{Anything passed will be ignored}
}
\value{
  A numeric vector of ones.
}
\references{
  Braumoeller, Bear F. (2003) "Causal Complexity and the Study
  of Politics." \emph{Political Analysis} 11(3): 209-233.
}
\author{
Bear F. Braumoeller, Harvard University, \email{bfbraum@fas.harvard.edu},\cr
Ben Goodrich, Harvard University, \email{goodrich@fas.harvard.edu}, and\cr
Jacob Kline, Harvard University, \email{jkline@fas.harvard.edu}
}

\keyword{methods}

