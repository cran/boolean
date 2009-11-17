\name{coef.booltest}
\alias{coef.booltest}
\alias{coef}
\title{Coefficients for Boolean Binary Response Models}
\description{
  This function implements a 'coef' S3 method for objects of class 'booltest'.
  If the model was estimated by maximum likelihood, the maximum likelihood
  estimates are returned. If the model was estimated by Metropolis-Hastings or
  bootstrapped, the mean over the iterations for each estimated parameter is
  returned with a warning.
}
\usage{
## S3 method for class 'booltest'
\method{coef}{booltest}(object, ...)
}

\arguments{
  \item{object}{An object of class 'booltest'}
  \item{...}{Anything passed will be ignored.}
}
\value{
  A numeric vector with one element per estimated parameter.
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

