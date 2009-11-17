\name{summary.booltest}
\alias{summary.booltest}
\alias{summary}
\title{'summary' Method for Boolean Binary Response Models}
\description{
'summary' method for objects of class 'booltest'
}

\usage{
## S3 method for class 'booltest'
\method{summary}{booltest}(object, ...)
}

\arguments{
  \item{object}{An object of class 'booltest'}
  \item{...}{Further arguments}
}
\value{
  An object of class 'summary.booltest', a list with components
  \item{formula}{The formula from \code{\link{boolprep}}}
  \item{summary.mat}{A matrix of coefficients, standard errors, etc.}
  \item{log.likelihood}{The log-likelihood at the optimum, NULL if 
                        Metropolis-Hastings was used to estimate}
  \item{iterations}{The number of iterations, NULL if Metropolis-
                    Hastings was used to estimate.}
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

