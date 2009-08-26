\name{effectplot}
\alias{effectplot}
\title{Visually Assess Nonlinearity in Effects of Explanatory Variables}
\description{
  This function plots the change in probability of the outcome in one
  explanatory variable as a function of a change in another explanatory variable.
}
\usage{
effectplot(z.out, variables, delta = 1, CI = 95, truehist = TRUE,
                       legend = TRUE, polygon = FALSE, ...)
}
\arguments{
  \item{z.out}{An object produced by \code{\link[Zelig]{Zelig}}.}
  \item{variables}{A character vector of length two. The first indicates which 
                   explanatory variable should constitute the x-axis of the plot.
                   The second indicates which explanatory variable should be 
                   perturbed by \code{delta}}
  \item{delta}{The amount by which to perturb the sample values of the second element 
               of \code{variable}.}
  \item{CI}{Confidence interval for (the change in) the probability of the outcome.}
  \item{truehist}{If TRUE, a histogram for \code{variable} is included in the plot.}
  \item{legend}{If TRUE, a legend is included in the plot.}
  \item{polygon}{If TRUE, then the confidence intervals as filled as ploygons.}
  \item{...}{Further arguments passed to \code{plot} or \code{\link[Zelig]{setx}}.}
}
\value{
  Nothing is returned.
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
\keyword{models}
