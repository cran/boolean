\name{boolprof}
\alias{boolprof}

\title{Graph Profile Likelihoods for Boolean Binary Response Models}
\description{
  This command produces a graph of the profile likelihood for a
  given coefficient following the estimation of a Boolean binary
  response model.
}
\usage{
boolprof(object, gvar = NULL, range = NULL, M = 100)
}

\arguments{
  \item{object}{ The object created by the boolean command. }
  \item{gvar}{ A character string indicating the name of variable associated 
     with coefficient of interest. By default, \code{boolprof} plots the profile
     likelihood for each estimated parameter. Use of the default is
     encouraged. Otherwise, be sure to specify the name of the variable
     including the "id" that is tacked onto the end of the variable 
     name so as to identify which instance of the variable is of 
     interest, even if the variable appears only once.}
  \item{M}{ The number of points to plot. The default (100) is usually
     sufficient.}
  \item{range}{ Range of points over which graph is plotted (e.g.,
    \code{range=seq(0,1,0.01)}). If NULL (the default), \code{M} random draws
     from a multivariate normal distribution are plotted consistent 
     with the mode and Hessian for models estimated by maximum
     likelihood. For models estimated via bootstrapping or MCMC,
     \code{M} values are sampled randomly from the output.}
}

\value{
  Nothing is returned. A graph is plotted.
}

\references{ Braumoeller, Bear F. (2003) "Causal Complexity and the
     Study of Politics." \emph{Political Analysis} 11(3): 209-233. }

\author{
Bear F. Braumoeller, Harvard University, \email{braumoeller.1@polisci.osu.edu},\cr
Ben Goodrich, Harvard University, \email{goodrich@fas.harvard.edu}, and \cr
Jacob Kline, Harvard University, \email{jkline@fas.harvard.edu}
}

\seealso{ \code{\link{boolean}} and \code{\link{boolprep}}}

\examples{
set.seed(50)
x1<-rnorm(1000)
x2<-rnorm(1000)
x3<-rnorm(1000)
x4<-rnorm(1000)
x5<-rnorm(1000)
x6<-rnorm(1000)
y<-1-(1-pnorm(-2+0.33*x1+0.66*x2+1*x3)*1-(pnorm(1+1.5*x4-0.25*x5)*pnorm(1+0.2*x6)))
y <- y>runif(1000)
bp <- boolprep("(a|(b&c))", y, a = ~ x1 + x2 + x3, b = ~ x4 + x5, c = ~ x6)
answer <- boolean(bp, link = c("probit", "logit", "cloglog"))

## Plot profiles
boolprof(answer)

}
\keyword{logic}
\keyword{optimize}
\keyword{distribution}
\keyword{models}
