\name{boolean}
\alias{boolean}
\alias{Boolean}
\title{ Partial-Observability Binary Response Models for Testing 
  Boolean Hypotheses }
\description{
  Boolean binary response models are a family of partial-observability
  \emph{n}-variate models designed to permit researchers to model
  causal complexity, or multiple causal "paths" to a given
  outcome.
}
\usage{
boolean(formula, data, link = "logit", start.values = NULL,
        method = "nlm", bootstrap = 0, control = list(fnscale = -1), 
        weights = NULL, robust = FALSE, ...)
}
\arguments{
  \item{formula}{The object producted by \code{boolprep}; 
    see \code{\link{boolprep}} for details.}
  \item{data}{An optional dataframe that contains the variables specified
    in the call to \code{boolprep}.}
  \item{link}{A character string or character vector of length \emph{n} that 
    identifies the CDF(s) to be used for each causal path. If only one link function is given, 
    it is used for all \emph{n} causal paths.The string may be "cauchit",
    "cloglog", "log", "logit" (the default), "probit", "scobit", "scobitL", "scobitR", or 
    any combination thereof in a character vector of length \emph{n}.  Aside
    from "scobit", "scobitL", and "scobitR" these are the same as used in the \code{glm}
    function. For details, see \code{\link{binomial}}. For documentation of
    the scobit link function, see the Nagler reference below; "scobitR" is an alias for
    "scobit" and scobitL is a transformation of the scobit CDF that is left-skewed.}
  \item{start.values}{If not NULL (the default), a vector of starting values.
     Such a vector must have length equal to the number of estimated parameters
     and should be ordered so that the starting values for the coefficients come
     before the starting values for the ancillary parameter for any "scobit", "scobitL"
     or "scobitR" link functions. Otherwise, the parameters should be ordered left-to-right
     based on their order in the call to \code{\link{boolprep}}. The \code{start.values} 
     argument is useful for checking whether the results are robust to overdispersion of starting
     values. If NULL, the starting values will be created from a series of calls to 
     \code{\link{glm}} using the appropriate link function, assuming the ancillary parameter
     is one in the case of the various scobits.}
  \item{method}{ A character string or character vector of length two that specifies the 
     algorithm used to estimate the boolean model. If \code{bootstrap > 0} and \code{method}
     is of length two, then the first element of \code{method} is used to obtain the 
     maximum-likelihood estimates (presumably with a computationally intensive algorithm) 
     and the second element is used for the bootstraps (possibly using a faster algorithm).
     Possible choices are: "nlm" (the default) "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", 
     "SANN" (the preceding five represent the possibilities for \code{\link{optim}}), "constrOptim", 
     "genoud", or "MH" (short for Metropolis-Hastings without simulating annealing). 
     For "nlm", see \code{\link{nlm}}; for the optim-related methods, see \code{\link{optim}}; 
     for "constrOptim", see \code{\link{constrOptim}}; for genoud, see \code{\link[rgenoud]{genoud}}; 
     for "MH" (Metropolis-Hastings), see 
     \code{\link[MCMCpack]{MCMCmetrop1R}}. Additional arguments for each of these methods can 
     (and often should) be passed through the \code{...} as indicated in the respective documentations. 
     If \code{method = "constrOptim"}, the constrained optimization uses the BFGS algorithm. Note that 
     for previous versions of this package, the \code{method} argument did something very different.}
  \item{bootstrap}{The number of bootstraps to execute, which must be a non-negative integer.
     If zero (the default), standard errors are approximated using the Hessian.
     Note that bootstrapping is appropriate only if the \code{method} is something
     other than "MH".}
  \item{control}{A list, identical in usage with that for \code{optim}.}
  \item{weights}{A vector of weights with length equal to the number of valid observations, which
     is used to produce a weighted sum of the observation-specific log-likelihoods. If NULL (the
     default), a weight of one is used for each observation. Currently, only NULL is supported. }
  \item{robust}{If TRUE a sandwich-style variance-covariance matrix would be estimated. However,
     this option is not current supported.}
  \item{...}{Additional arguments to be passed to the procedure identified by the \code{method} 
     argument. Such arguments are often necessary to obtain reasonable results but are too numerous 
     to be listed here. However, it is not possible to pass a function that calculates the gradient.}
}
\details{
  Boolean permits estimation of Boolean binary response models
 (see Braumoeller 2003 for derivation), which are a
 family of partial-observability \emph{n}-variate models designed to permit
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
  Returns an object of class booltest, with the following list elements:
  \item{value}{The maximized value of the log-likelihood or NA, if \code{method="MH"}.}
  \item{counts}{The number of calls to the function and the gradient, or
    the number of MCMC iterations if \code{method="MH"}.}
  \item{convergence}{The convergence code for maximum-likelihood methods 
     or NA if \code{method = "MH"}, which does NOT indicate the
     presence, absence, or irrelevance of convergence of the Markov chain
     to the stationary distribution.}
  \item{message}{The convergence message for maximum-likelihood methods, (NA
     if \code{method = "nlm")} or NA if \code{method = "MH"}, which does NOT 
     indicate the presence, absence, or irrelevance of convergence of the Markov chain 
     to the stationary distribution.}
   \item{hessian}{The Hessian matrix, unless \code{method = "MH"}, in which case, 
     NULL. If \code{variance} is a "caught" error, examining the Hessian may 
     indicate where the problem is.}
   \item{coefficients}{Either a vector of maximum-likelihood estimates, or
     if \code{method = "MH"} or \code{bootstap > 0}, a matrix of MCMC
     or bootstrapped samples with rows equal to the number of samples and
     columns equal to the number of estimated parameters.}
  \item{variance}{A "try" of the variance-covariance matrix, unless 
     \code{method = "MH"}, in which case NA.}
   \item{boolean.call}{The call to the \code{boolean} function.}
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
\note{ Examining profile likelihoods with \code{boolprof} is highly
    recommended.  These are partial observability models are generically 
    starved for information; as a result, maximum likelihood estimation
    can encounter problems with plateaus in likelihood functions even
    with very large datasets. In principle, if \code{method = "BHHH"}, the
    Hessian can be approximated in circumstances where other maximization methods
    might fail to estimate the Hessian due to plateaus or discontinuities at
    the maximum. However, such estimates should be analyzed with caution.}

\seealso{ \code{\link{boolprep}} to prepare the structure of the model, and
  \code{\link{boolprof}} to produce profile likelihoods after estimation. }

\examples{
set.seed(50)
x1 <- rnorm(1000)
x2 <- rnorm(1000)
x3 <- rnorm(1000)
x4 <- rnorm(1000)
x5 <- rnorm(1000)
x6 <- rnorm(1000)
y<-1-(1-pnorm(-2+0.33*x1+0.66*x2+1*x3)*1-(pnorm(1+1.5*x4-0.25*x5)*pnorm(1+0.2*x6)))
y <- y>runif(1000)
bp <- boolprep("(a|(b&c))", "y", a = ~ x1 + x2 + x3, b = ~ x4 + x5, c = ~ x6)
answer <- boolean(bp, link = c("probit", "logit", "cloglog"), start.values = ## For speed
                  c(-1.750,  0.354,  0.698,  1.231,  1.473,  2.628, -0.452,  0.764,  0.173))

## Plot profiles
boolprof(answer)


## Examine coefficients
coef(answer)

## Examine coefficients, standard errors, etc.
summary(answer)

}
\keyword{logic}
\keyword{optimize}
\keyword{distribution}
\keyword{models}
