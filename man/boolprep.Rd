\name{boolprep}
\alias{boolprep}
\title{ Generate Functional Form for Boolean Logit or Probit }
\description{
  Takes as arguments the logical form of a Boolean model, the
  dependent variable, and the vectors of independent variables, and
  produces a functional form suitable for inclusion in a
  \code{\link{boolean}} logit or probit command.
}
\usage{
boolprep(form, depvar, a, b = "", c = "", d = "", e = "", f = "",
g = "", h = "", i = "", j = "", k = "", l = "", m = "", n = "",
o = "", p = "", q = "", r = "", s = "", t = "", u = "", v = "",
w = "", x = "", y = "", z = "", constant = TRUE)
}
\arguments{
  \item{form}{ Form of Boolean model, e.g., \code{a&b}, 
     \code{a&(b|c)}, etc. }
  \item{depvar}{ Dependent variable. }
  \item{a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z}{ Vectors
     of independent variables corresponding to \code{a}, \code{b}, etc.,
     in \code{form}. }
  \item{constant}{ Indicates whether constants should be placed in
     each vector of independent variables (note that the syntax
     requires that constant terms be explicitly included). }
}
\details{
  Because the structure of a Boolean logit or probit model can
  be arbitrarily complex, the layering of parentheses soon becomes
  tedious and unnecessarily difficult.  To make life easier, the
  \code{boolprep} command automates the task of creating the
  appropriate argument.
}
\value{
  A text string of the form \code{y ~ f(x1...xn)}.}
}
\references{  Braumoeller, Bear F. (2003) "Causal Complexity and the Study
  of Politics." \emph{Political Analysis} 11(3): 209-233.
}

\author{
Bear F. Braumoeller, Harvard University, \email{bfbraum@fas.harvard.edu} \cr
Jacob Kline, Harvard University, \email{jkline@fas.harvard.edu}
}

\seealso{ \code{\link{boolean}}, \code{\link{boolprof}}, 
     and \code{\link{boolfirst}}. }

\examples{
boolprep("a&b", "y", "x1+x2", "x3+x4+x5")

boolprep("a|(b&c)", "y", "x1+x2+x3", "x4+x5", "x6")

boolprep("(a|b|c)&(d|e)&f", "y", "x1+x2", "x3+x4+x5", "x6", "x7+x8", "x9+x10+x11", "x12+x13")
}
\keyword{logic}
\keyword{optimize}
\keyword{distribution}
\keyword{models}
