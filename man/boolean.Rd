\name{boolean}
\alias{boolean}
\alias{Boolean}
\title{ Partial-Observability Logit or Probit Models for Testing 
  Boolean Hypotheses }
\description{
  Boolean logit and probit are a family of partial-observability
  \emph{n}-variate models designed to permit researchers to model
  causal complexity, or multiple causal "paths" to a given outcome
  (e.g., a situation in which either a or b will produce y, one in
  which (a or b) and c produce y, and so forth).
}
\usage{
boolean(structure, method, maxoptions = "", optimizer="nlm",
        safety=1, bootstrap=FALSE, bootsize=100)
}
\arguments{
  \item{structure}{ 
  Structure of equation to be estimated, in standard y ~ f(x) form,
  using \code{&} to represent the Boolean operator "and" and \code{|}
  to represent the Boolean operator "or."  (Note that the syntax
  requires that constants be entered explicitly; see the entry for
  \code{\link{boolprep}} for details.)  Be sure to enter the correct
  functional form and balance parentheses; if in doubt, or just for
  convenience, use the \code{boolprep} command to prepare structure
  prior to estimation. }
  \item{method}{ Either "logit" or "probit". }
  \item{maxoptions}{ Maximization options (see \code{\link{nlm}} or
  \code{\link{optim}} for details). }
  \item{optimizer}{ Either "nlm" or "optim". }
  \item{safety}{ Number of search attempts.  The likelihood functions
  implied by Boolean procedures can become quite convoluted; in such
  cases, multiple searches from different starting points can be run.
  Works only when using \code{nlm}. }
  \item{bootstrap}{ If TRUE, bootstraps standard errors. }
  \item{bootsize}{ Number of iterations if bootstrap=TRUE. }
}
\details{
  Boolean permits estimation of Boolean logit and probit equations
  (see Braumoeller 2003 for derivation).  Boolean logit and probit
  model situations in which any number of antecedent conditions
  (a, b, c, ...) occur with probabilities that can be modeled using
  standard logit or probit curves, and the antecedent conditions
  combine in a manner described by Boolean logic to produce the
  dependent variable.  For example, theory might suggest that a binary
  dependent variable y occurs only when conditions a and b are present,
  that the probability of a's occurrence is influenced by variables
  x1...x4, and that the probability of b's occurrence is influenced by
  variables x5...x8.  Alternately, theory might suggest that y occurs
  if (a and b) or (c and d) are present (in other words, that a-d are
  INUS conditions; see Mackie 1965), and that the probability of the
  occurrence of each antecedent condition is influenced by some vector
  of independent variables.  Boolean logit and probit are designed for
  use in such situations.
}
\value{
  Returns an object of class booltest, with slots @Calculus, @LogLik,
  @Variables, @Coefficients, @StandardErrors, @Iterations, @Hessian,
  @Gradient, @Zscore, @Probz, @Conf95lo, @Conf95hi, @pstructure, and
  @method (note that some slots may be left empty if the relevant
  information is not furnished by the maximizer).
}
\references{  Braumoeller, Bear F. (2003) "Causal Complexity and the Study
  of Politics." \emph{Political Analysis} 11(3): 209-233.

  Mackie, J. L. (1965) "Causes and Conditions." \emph{American
  Philosophical Quarterly} 2(4): 245-264.}
\author{
Bear F. Braumoeller, Harvard University, \email{bfbraum@fas.harvard.edu} \cr
Jacob Kline, Harvard University, \email{jkline@fas.harvard.edu}
}
\note{ Examining profile likelihoods with \code{boolprof} is highly
  recommended.  Boolean logit and probit are partial observability models,
  which are generically starved for information; as a result, maximum
  likelihood estimation can encounter problems with plateaus in likelihood
  functions even with very large n. }

\seealso{ \code{\link{boolprep}} to prepare structure of equation,
  \code{\link{boolfirst}} to graph first differences after estimation, and
  \code{\link{boolprof}} to produce profile likelihoods after estimation. }

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

## Examine coefficients, standard errors, etc.
summary(answer)

## Examine "summary" output plus Hessian, gradient, etc.
show(answer)

## Plot first differences for model
plot(answer)

## Plot profiles
plot(answer, panel="boolprof")
}
\keyword{logic}
\keyword{optimize}
\keyword{distribution}
\keyword{models}
