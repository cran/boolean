\name{booltest}
\alias{booltest}
\alias{booltest-class}
\title{ Boolean logit/probit objects }
\description{
  Boolean logit and probit procedures produce an object of class
  booltest.  Such an object contains slots @Calculus, @LogLik,
  @Variables, @Coefficients, @StandardErrors, @Iterations, @Hessian,
  @Gradient, @Zscore, @Probz, @Conf95lo, @Conf95hi, @pstructure, and
  @method (note that some slots may be left empty if the relevant
  information is not furnished by the maximizer).

  A booltest object may be summarized using summary().  More detailed
  output can be obtained by show().  The specialized functions
  \code{\link{boolfirst}} and \code{\link{boolprof}} produce plots of
  first differences and profile likelihoods, respectively.

  In addition, \code{plot()} produces a matrix of plots.  \code{plot()}
  takes a single argument, \code{panel}, that determines whether the
  matrix consists of first differences (\code{panel="boolfirst"}) or
  profile likelihoods (\code{panel="boolprof"}).  The default is first
  differences.
}
\references{  Braumoeller, Bear F. (2003) "Causal Complexity and the Study
  of Politics." \emph{Political Analysis} 11(3): 209-233.
}

\author{
Bear F. Braumoeller, Harvard University, \email{bfbraum@fas.harvard.edu} \cr
Jacob Kline, Harvard University, \email{jkline@fas.harvard.edu}
}
\seealso{ \code{\link{boolean}}, \code{\link{boolfirst}}, and \code{\link{boolprof}}.
}

\keyword{logic}
\keyword{optimize}
\keyword{distribution}
\keyword{models}
