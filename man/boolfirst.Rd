\name{boolfirst}
\alias{boolfirst}
\title{ Graph first differences after Boolean logit/probit }
\description{
  This function graphs first differences for a given variable,
  with all other variables held at their means, following the
  execution of a \code{\link{boolean}} logit or probit command.
}
\usage{
boolfirst(object, gvar, instance = 0, range = 0)
}
\arguments{
  \item{object}{ The object created by the boolean command. }
  \item{gvar}{ Name of variable associated with coefficient of
     interest. }
  \item{instance}{ If variable occurs more than once in equation,
     permits user to specify which occurrence corresponds to the
     coefficient of interest. }
  \item{range}{ Range of points over which graph is plotted (e.g.,
  range=seq(0,1,0.01)) }
}
\references{ Braumoeller, Bear F. (2003) "Causal Complexity and the
     Study of Politics." \emph{Political Analysis} 11(3): 209-233. }

\author{
Bear F. Braumoeller, Harvard University, \email{bfbraum@fas.harvard.edu} \cr
Jacob Kline, Harvard University, \email{jkline@fas.harvard.edu}
}

\seealso{ \code{\link{boolean}}, \code{\link{boolprep}}, 
     and \code{\link{boolprof}}. }

\examples{
library("boolean")
set.seed(50)
x1<-rnorm(1000)
x2<-rnorm(1000)
x3<-rnorm(1000)
x4<-rnorm(1000)
x5<-rnorm(1000)
x6<-rnorm(1000)
e1<-rnorm(1000)/3
e2<-rnorm(1000)/3
e3<-rnorm(1000)/3
y<-1-(1-pnorm(-2+0.33*x1+0.66*x2+1*x3+e1)*1-(pnorm(1+1.5*x4-0.25*x5+e2)*pnorm(1+0.2*x6+e3)))
y <- y>runif(1000)
answer <- boolean(y ~( ((cons+x1+x2+x3)|((cons+x4+x5)&(cons+x6))) ), method="probit")
boolfirst(answer, "x4")
boolprof(answer, "x1", instance=2)
}
\keyword{logic}
\keyword{optimize}
\keyword{distribution}
\keyword{models}
