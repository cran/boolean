\name{boolean}
\alias{boolean}
\alias{Boolean}
\title{ Partial-Observability Logit or Probit Models for Testing 
  Boolean Hypotheses }
\description{
  Boolean logit and probit are a family of partial-observability
  \emph{n}-variate models designed to permit researchers to model
  causal complexity, or multiple causal "paths" to a given
  outcome.
}
\usage{
boolean(structure, method, maxoptions = "", optimizer="nlm",
        safety=1, bootstrap=FALSE, bootsize=100, popsize=5000)
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
  \item{optimizer}{ Either "nlm", "optim", or "genoud". }
  \item{safety}{ Number of search attempts.  The likelihood functions
  implied by Boolean procedures can become quite convoluted; in such
  cases, multiple searches from different starting points can be run.
  Works only when using \code{nlm}. }
  \item{bootstrap}{ If TRUE, bootstraps standard errors. }
  \item{bootsize}{ Number of iterations if bootstrap=TRUE. }
  \item{popsize}{ Population size if optimizer=genoud. }
}
\details{
  Boolean permits estimation of Boolean logit and probit models
 (see Braumoeller 2003 for derivation).  Boolean logit and probit are a
 family of partial-observability n-variate models designed to permit
 researchers to model causal complexity, or multiple causal "paths" to a
 given outcome.  The various "paths" are modeled as latent dependent
 variables that are multiplied together in a manner determined by the
 logic of their (Boolean) interaction.  If, for example, we wanted to
 model a situation in which diet OR smoking causes heart failure, we
 would use one set of independent variables (caloric intake, fat intake,
 etc.) to predict the latent probability of diet-related coronary
 failure (y1*), use another set of variables (cigarettes smoked per day,
 exposure to second-hand smoke, etc.) to predict the latent probability
 of smoking-related coronary failure (y2*), and model the observed
 outcome (y, or coronary failure) as a function of the Boolean
 interaction of the two: Pr(y=1) = 1-([1-y1*] x [1-y2*]).  Independent
 variables that have an impact on both latent dependent variables can
 be included in both paths.  Any combination of ANDs and ORs can be
 posited, and the interaction of any number of latent dependent
 variables can be modeled, although the procedure becomes exponentially
 more data-intensive as the number of latent dependent variables
 increases.
}
\value{
  Returns an object of class booltest, with slots @Calculus, @LogLik,
  @Variables, @Coefficients, @StandardErrors, @Iterations, @Hessian,
  @Gradient, @Zscore, @Probz, @Conf95lo, @Conf95hi, @pstructure, and
  @method (note that some slots may be left empty if the relevant
  information is not furnished by the maximizer).
}
\references{
  Braumoeller, Bear F. (2003) "Causal Complexity and the Study
  of Politics." \emph{Political Analysis} 11(3): 209-233.
}
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
y<-1-(1-pnorm(-2+0.33*x1+0.66*x2+1*x3)*1-(pnorm(1+1.5*x4-0.25*x5)*pnorm(1+0.2*x6)))
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
