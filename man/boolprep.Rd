\name{boolprep}
\alias{boolprep}
\title{ Generate Functional Form for Boolean binary choice model }
\description{
  Takes as arguments the logical form of a Boolean model, the
  dependent variable, and sets of independent variables, and
  produces a functional form necessary for inclusion in a
  \code{\link{boolean}} command.
}
\usage{
boolprep(FORM, DEPVAR, ...)
}
\arguments{
  \item{FORM}{Form of Boolean model, e.g., \code{"a&b"}, 
     \code{"a&(b|c)"}, etc. Unlike previous versions of this package,
      you are not restricted to single letters. For example,
      \code{"(lucky|good)"} is now allowed as well. However, \code{FORM}
      must be a character string. }
  \item{DEPVAR}{Dependent variable, which must be a logical vector
      or a numeric vector with zeroes and ones only, or
      the name of such a dependent variable can be given as a 
      character string. }
  \item{...}{If there are n logical operators in FORM, users are
      required to pass n+1 arguments through the dots. These n+1
      arguments should be one-sided formulas and follow the same 
      syntax as in the glm() function (see \code{\link{formula}}), 
      such as \code{~ x1 + x2}. Note that if \code{FORM = "(lucky|good)"}, 
      then the required arguments would take a form such as 
      \code{lucky = ~x1, good = ~x2}. }
}
\details{
  Because the structure of a Boolean binary choice model can
  be arbitrarily complex, the layering of parentheses soon becomes
  tedious and unnecessarily difficult.  To make life easier, the
  \code{boolprep} command automates the task of creating the
  appropriate argument and checks for syntactical errors.
}
\value{
  An object of class boolprep with the structure of the model and the 
  appropriate formulas as list elements, none of which are expected to be 
  accessed by the user. However, \code{\link{boolean}} requires that a 
  boolprep object be passed to it.
}
\references{
  Braumoeller, Bear F. (2003) "Causal Complexity and the Study
  of Politics." \emph{Political Analysis} 11(3): 209-233.
}

\author{
Bear F. Braumoeller, Harvard University, \email{braumoeller.1@polisci.osu.edu}\cr
Ben Goodrich, Harvard University, \email{goodrich@fas.harvard.edu}, and\cr
Jacob Kline, Harvard University, \email{jkline@fas.harvard.edu}
}

\seealso{\code{\link{boolean}}}

\examples{
y  <- rbinom(1000, 1, 0.5)

bp <- boolprep("a & b", "y", a = ~x1+x2, b = ~x3+x4+x5)

bp <- boolprep("rich | (smart & diligent)", y, rich = ~ x1+x2+x3, 
               smart = ~ x4+x5, diligent = ~ x6)

bp <- boolprep("(a|b|c)&(d|e)&f", y, a=~x1+x2, b = ~x3+x4+x5, 
               c = ~x6, d = ~x7+x8, e = ~x9+x10+x11, f = ~x12+x13)
}
\keyword{logic}
\keyword{optimize}
\keyword{distribution}
\keyword{models}
