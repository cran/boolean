\name{zelig.boolprep}
\alias{zelig.boolprep}
\title{'boolprep' Method for Estimation Via 'zelig'}
\description{
This function is an interface to the \code{\link[Zelig]{Zelig}} function
in the Zelig library. It should not be called directly by the user.
}

\usage{
## S3 method for class 'booltest'
\method{zelig}{boolprep}(formula, model, data, by, save.data = FALSE, cite, ...)
}

\arguments{
  \item{formula}{An object of class 'boolprep'}
  \item{data}{A dataframe containing all the relevant variables.}
  \item{model}{A character vector, "boolean" in this case.}
  \item{by}{A character vector with the name of a stratifying variable
            so that the model is estimated on one subset of observations
            at a time.}
  \item{save.data}{Should the data be saved in the resulting object, defaults
                   to FALSE}
  \item{cite}{I do not know what this does}
  \item{...}{Further arguments passed}
}
\value{
  An object of class 'booltest', see \link{boolean} for details.
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

