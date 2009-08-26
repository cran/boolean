#     This file is part of boolean, a program to estimate boolean models in R
#     Copyright 2003 -- 2009 Bear F. Braumoeller
#
#     Some portions of this code are derived from code in the Zelig package, which is
#     Copyright 2004 Kosuke Imai, Gary King and Olivia Lau and licensed under the GPL
#     version 2 or later.
#
#     boolean is free software: you can redistribute it and/or modify
#     it under the terms of the GNU Affero General Public License as published by
#     the Free Software Foundation, either version 2 of the License, or
#     (at your option) any later version.
#
#     boolean is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU Affero General Public License for more details.
#
#     You should have received a copy of the GNU Affero General Public License
#     along with boolean.  If not, see <http://www.gnu.org/licenses/>.


.onAttach <- function( ... ) {

param.booltest <<-
function(object, num, bootstrap = FALSE) {
  if(class(object)[2] != "ML") return(object$coefficients)
  if(!is.matrix(object$variance) || any(is.na(object$variance))) 
	stop("It is not possible to draw from a multivariate-normal sampling distribution because the variance-covariance matrix could not be estimated", call. = FALSE)
  vcv <- object$variance
  c.vcv <- chol(vcv)
  mat <- matrix(rnorm(ncol(vcv)*num), nrow = ncol(vcv), ncol = num)
  res <- t(t(c.vcv) %*% mat + object$coefficients)
  return(res)
}

qi.booltest <<- 
function(object, simpar = simpar, x = x, x1 = x1, y = NULL) {

prob <- function(par, X, map = MAP, predicted = FALSE, ...) {
  scobitMark <- ncol(X) + 1
  holder <- array(NA, dim = c(nrow(X), switch(is.matrix(par) + 1, nrow(map), c(nrow(map), nrow(par)))))
  if(!predicted){
    for(i in 1:nrow(map)) {
      SEQ <- map[i,1]:map[i,2]
      invlink <- boolean:::make.invlink(link[i])
      if(pmatch("scobit", link[i], nomatch = FALSE) > 0) {
        if(is.matrix(holder)) holder[,i] <- invlink(X[,SEQ] %*% par[SEQ], par[scobitMark])
        else holder[,i,] <- invlink(X[,SEQ] %*% t(as.matrix(par[,SEQ])), par[1,scobitMark])
        scobitMark <- scobitMark + 1
      }
      else {
        if(is.matrix(holder)) holder[,i] <- invlink(X[,SEQ] %*% par[SEQ])
        else holder[,i,] <- invlink(X[,SEQ] %*% t(par[,SEQ]))
      }
    } # END for loop
  }

  else {                                                                     # Obtain *predicted* probabilities
    num <- dim(holder)[1] * dim(holder)[3]                                   # The number of rbinom draws we need
    for(i in 1:nrow(MAP)) {                                                  # Otherwise the same as above
      SEQ <- map[i,1]:map[i,2]
      invlink <- boolean:::make.invlink(link[i])
      if(pmatch("scobit", link[i], nomatch = FALSE) > 0) {
        holder[,i,] <- rbinom(num, 1, invlink(X[,SEQ] %*% t(par[,SEQ]), par[1,scobitMark]))
        scobitMark <- scobitMark + 1
      }
      else holder[,i,] <- rbinom(num, 1, invlink(X[,SEQ] %*% t(par[,SEQ])))
    } # END for loop
  }

  probabilities <- bool.fun(holder)
  probabilities
}


  dta <- eval(object$boolean.call$data)
  admin <- boolean:::administrative(object$boolean.call$formula, dta)
  MAP <- admin$map

  if(is.null(object$boolean.call$link)) link <- formals(boolean)$link
  else {
    link <- as.character(object$boolean.call$link)
    link <- link[link != "c"]
  }
  if(length(link) == 1) link <- rep(as.character(link), nrow(admin[[3]]))

  bool.fun <- boolean:::make.bool.fun(object$boolean.call$formula)

  boolean.env <- new.env()
  environment(prob) <- boolean.env
  if(!is.null(y)) environment(y) <- boolean.env
  if(!is.null(x1)) environment(x1) <- boolean.env
  environment(x) <- boolean.env
  environment(MAP) <- boolean.env
  environment(link) <- boolean.env
  environment(bool.fun) <- boolean.env
  assign("%&%", boolean:::"%&%", envir = boolean.env)
  assign("%|%", boolean:::"%|%", envir = boolean.env)

  ev <- prob(par = simpar, X = x, map = MAP, predicted = FALSE)
  pr <- as.character(prob(par = simpar, X = x, map = MAP, predicted = TRUE))
  qi <- list(ev = ev, pr = pr)
  qi.name <- list(ev = "Expected Values: E(Y|X)",
                  pr = "Predicted Values: Y|X")
  if (!is.null(x1)){
    ev1 <- prob(par = simpar, X = x1, map = MAP, predicted = FALSE)
    qi$fd <- ev1-ev
    qi.name$fd <- "First Differences in Expected Values: E(Y|X1)-E(Y|X)"
    qi$rr <- ev1/ev
    qi.name$rr <- "Risk Ratios: P(Y=1|X1)/P(Y=1|X)"
  }
  if (!is.null(y)) {
    yvar <- matrix(rep(y, nrow(simpar)), nrow = nrow(simpar), byrow = TRUE)
    tmp.ev <- yvar - qi$ev
    tmp.pr <- yvar - as.integer(qi$pr)
    qi$ate.ev <- matrix(apply(tmp.ev, 1, mean), nrow = nrow(simpar))
    qi$ate.pr <- matrix(apply(tmp.pr, 1, mean), nrow = nrow(simpar))
    qi.name$ate.ev <- "Average Treatment Effect: Y - EV"
    qi.name$ate.pr <- "Average Treatment Effect: Y - PR"
  }
  list(qi=qi, qi.name=qi.name)
}

}