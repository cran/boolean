\name{setx.booltest}
\alias{setx.booltest}
\title{'booltest' Method for Post-Estimation Via 'zelig'}
\description{
This function is an interface to the \code{\link[Zelig]{setx}} function
in the Zelig library. It should not be called directly by the user.
}

\usage{
## S3 method for class 'booltest'
\method{setx}{booltest}(object, fn = list(numeric = mean, ordered = median, other = mode),
              data = NULL, cond = FALSE, counter = NULL, ...)  # good
}

\arguments{
  \item{object}{An object of class 'boolprep' estimated by \code{\link[Zelig]{Zelig}}.}
  \item{fn}{A list of functions; see \code{\link[Zelig]{setx}} for details.}
  \item{data}{An additional dataframe containing values for the explanatory variables
              or NULL to use the dataframe in the call to \code{\link[Zelig]{Zelig}}.}
  \item{cond}{A flag to indicate conditional prediction; 
              see \code{\link[Zelig]{setx}} for details.}
  \item{counter}{Should always be NULL.}
  \item{...}{Additional arguments to set explanatory variables to particular values.}
}
\value{
  A dataframe of explanatory variables set to the specified values; 
  see \code{\link[Zelig]{setx}} for details.
}
\author{
Bear F. Braumoeller, Harvard University, \email{bfbraum@fas.harvard.edu}, and\cr
Ben Goodrich, Harvard University, \email{goodrich@fas.harvard.edu}
}
\keyword{methods}

