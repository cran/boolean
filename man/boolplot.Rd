\name{boolplot}
\alias{boolplot}
\title{Plot Probabilities from Boolean Models}
\description{
  This function plots the probability of the outcome as a function of
  an explanatory variable or the change in an explanatory variable.
}
\usage{
boolplot(z.out, variable, delta = 0, suppression.factor = FALSE, CI = 95,
         truehist = TRUE, legend = TRUE, plot.both = FALSE, polygon = FALSE, yscale = NULL, ...)
}
\arguments{
  \item{z.out}{An object produced by \code{\link[Zelig]{Zelig}}.}
  \item{variable}{A character vector indicating which explanatory variable
                  should constitute the x-axis of the plot.}
  \item{delta}{The amount by which to perturb the sample values of \code{variable}
               when plotting the change in probability as a function of the change
               in \code{variabel}. If \code{plot.both = TRUE}, then \code{delta} must
               be non-zero.}
  \item{suppression.factor}{If true, then the complement of the probability of the 
                            outcome constitutes the y-axis.}
  \item{CI}{Confidence interval for (the change in) the probability of the outcome.}
  \item{truehist}{If TRUE, a histogram for \code{variable} is included in the plot.}
  \item{legend}{If TRUE, a legend is included in the plot.}
  \item{plot.both}{If TRUE, the probability of the outcome is plotted as a function
                   of the level of \code{variable} and the change in the probability of
                   the outcome is plotted as a function of \code{variable} plus \code{delta}.}
  \item{polygon}{If TRUE, then the confidence intervals as filled as ploygons.}
  \item{yscale}{If not NULL, then it should be a numeric vector of length two indicating the
                lower and upper limits of the y-axis in the plot.}
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
