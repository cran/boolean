\name{latent}
\alias{latent}
\title{Extract Probabilities from Each Path from Boolean Models}
\description{
  Boolean models involve more than one 'path', e.g. if y~a&(b|c), then
  a, b, and c are paths that are parameterized with explanatory variables.
  This function extracts estimates of a, b, and c, which are probabilities
  but not the probability that the dependent variable occurs.
}
\usage{
latent(object, probability = FALSE, invMills = FALSE)
}
\arguments{
  \item{object}{An object of class 'booltest'}
  \item{probability}{If TRUE, the linear predictor is transformed by the
     link function specified in the call to \code{\link{boolean}} to produce
     a vector of probabilities.}
  \item{invMills}{If TRUE, the inverse Mills ratio is calculated for paths
     where a probit link is specified in the call to \code{\link{boolean}}}.
  Note that if both probabilities and invMills are FALSE (the default), then
  the linear predictor for each path is returned.
}
\value{
  A numeric matrix with columns equal to the number of paths and rows equal to the
  number of observations. Whether the matrix contains the linear predictions,
  probabilities, or inverse Mills ratios depends on the call to \code{latent}.
}
\references{
  Braumoeller, Bear F. (2003) "Causal Complexity and the Study
  of Politics." \emph{Political Analysis} 11(3): 209-233.
  Nagler, Jonathon. (1994) "Scobit: An Alternative Estimator to Logit
  and Probit." \emph{American Journal of Political Science} 38(1): 230-255.
}
\author{
Bear F. Braumoeller, Harvard University, \email{braumoeller.1@polisci.osu.edu},\cr
Ben Goodrich, Harvard University, \email{goodrich@fas.harvard.edu}, and\cr
Jacob Kline, Harvard University, \email{jkline@fas.harvard.edu}
}
\keyword{models}