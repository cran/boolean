\name{sim.booltest}
\alias{sim.booltest}
\title{'booltest' Method for Parameter Simulation Via 'zelig'}
\description{
These functions provide an interface to the \code{\link[Zelig]{sim}} function
in the Zelig library. They should not be called directly by the user.
}

\usage{
## S3 method for class 'booltest'
\method{sim}{booltest}(object, x, x1=NULL, num=c(1000, 100),
             prev = NULL, bootstrap = FALSE, bootfn=NULL, cond.data = NULL, ...)
}
\arguments{
  \item{object}{An object of class 'booltest'estimated by \code{\link[Zelig]{Zelig}}.}
  \item{x}{An object produced by \code{\link[Zelig]{setx}}.}
  \item{x1}{An optional, additional object produced by \code{\link[Zelig]{setx}}.}
  \item{num}{The number of simulations of the parameters.}
  \item{prev}{An optional matrix of parameter simulations supplied by the user.}
  \item{bootstrap}{Should be FALSE because bootstrapping is supported directly in the call
                   to \code{\link[Zelig]{Zelig}} for boolean models. See also 
                   \link{boolean}.}
  \item{bootfn}{Should be NULL because bootstrapping is supported directly in the call
                   to \code{\link[Zelig]{Zelig}} for boolean models. See also 
                   \link{boolean}.}
  \item{cond.data}{Should be NULL.}
  \item{...}{Should not be needed for a boolean model.}
}
\value{
  The param and qi methods are used internally by \code{sim}, which in turn produces 
  quantities of interest. See \code{\link[Zelig]{sim}} for details; the output is the
  same as for a probit or logit model estimated via \code{zelig}.
}
\author{
Bear F. Braumoeller, Harvard University, \email{bfbraum@fas.harvard.edu}, and\cr
Ben Goodrich, Harvard University, \email{goodrich@fas.harvard.edu}
}
\keyword{methods}


