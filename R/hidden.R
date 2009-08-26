#     This file is part of boolean, a program to estimate boolean models in R
#     Copyright 2003 -- 2009 Bear F. Braumoeller
#
#     Some portions of this code are derived from code in the MCMCpack package, which is
#     Copyright 2003-2007 Andrew D. Martin and Kevin M. Quinn and licensed under the GPL
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


# Helper functions, consider moving these two  (to where? must be a Zelig thing)
"%&%" <- function(x,y) {x * y}                 ## The AND function
"%|%" <- function(x,y) {1 - (1 - x) * (1 - y)} ## The OR function

mlog <- function(x) {                                       # Takes the log, but won't return -Inf (might test on a Mac G5)
  ifelse(x < .Machine$double.xmax^-1.0498, -744.44, log(x)) # -744.44 looks like the smallest R can go without resorting to -Inf
}

mlog <- function(x) pmax(log(x), -744.44)                   # -744.44 looks like the smallest R can go without resorting to -Inf

make.invlink <- function (link) # check if this function has changed in base
{
    switch(link, logit = {
        linkinv <- function(eta) {
          thresh <- -log(.Machine$double.eps)
          eta <- pmin(pmax(eta, -thresh), thresh)
          exp(eta)/(1 + exp(eta))
        }
      }, probit = {
        linkinv <- function(eta) {
          thresh <- -qnorm(.Machine$double.eps)
          eta <- pmin(pmax(eta, -thresh), thresh)
          pnorm(eta)
        }
      }, cauchit = {
        linkinv <- function(eta) {
          thresh <- -qcauchy(.Machine$double.eps)
          eta <- pmin(pmax(eta, -thresh), thresh)
          pcauchy(eta)
        }
      }, cloglog = {
        linkinv <- function(eta) pmax(pmin(-expm1(-exp(eta)),
                                           1 - .Machine$double.eps), .Machine$double.eps)
      }, log = {
        linkinv <- function(eta) pmax(exp(eta), .Machine$double.eps)
      }, scobit = {
        linkinv <- function(eta, a) {
          thresh <- -log(.Machine$double.eps)
          eta <- pmin(pmax(eta, -thresh), thresh)
          (1 + exp(-eta))^(-exp(a))                         # Note parameterization of the exponent
        }
      }, scobitR = {                                        # Alias for scobit
        linkinv <- function(eta, a) {
          thresh <- -log(.Machine$double.eps)
          eta <- pmin(pmax(eta, -thresh), thresh)
          (1 + exp(-eta))^(-exp(a))                         # Note parameterization of the exponent
        }
      }, scobitL = {
        linkinv <- function(eta, a) {
          thresh <- -log(.Machine$double.eps)
          eta <- pmin(pmax(eta, -thresh), thresh)
          1 - (1 + exp(eta))^(-exp(a))                      # Note parameterization of the exponent
        }
      }, stop(sQuote(link), " link not recognized"))
    linkinv
}

parse.formula <- function(formula, data, intercept=TRUE) { # This MCMCpack function parses formulas

    # extract Y, X, and variable names for model formula and frame
  mt <- terms(formula, data=data)
  if(missing(data)) data <- sys.frame(sys.parent())
  mf <- match.call(expand.dots = FALSE)
  mf$intercept <- mf$justX <- NULL
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, sys.frame(sys.parent()))
  if (!intercept){
    attributes(mt)$intercept <- 0
  }

    # null model support
  X <- if (!is.empty.model(mt)) model.matrix(mt, mf, contrasts)
  X <- as.matrix(X)         # X matrix
  xvars <- dimnames(X)[[2]] # X variable names
  xobs  <- dimnames(X)[[1]] # X observation names
  Y <- as.character(formula[2])
  return(list(Y, X, xvars, xobs))
}

# This calls parse.formula for every mini-formula in the boolean statement and returns a list containing a big X matrix and
# a matrix indicating where the RHS data for each mini-GLM starts and stops in the big X matrix.
parse.list <- function(x = formula[-1], data) {
  x$depvar <- NULL
  X <- as.numeric(NULL)
  x.cols = matrix(NA, nrow = length(x), ncol = 2)
  for(i in 1:length(x)) {
    holder <- parse.formula(x[[i]], data)
    colnames(holder[[2]]) <- paste(holder[[3]], "_", names(x)[i], sep = "")
    X <- cbind(X, holder[[2]])

    if(i > 1) {
      x.cols[i,1] <- x.cols[i-1,2] + 1
      x.cols[i,2] <- x.cols[i,1] + ncol(holder[[2]]) - 1
    }
    else {
      x.cols[i,1] <- 1
      x.cols[1,2] <- ncol(holder[[2]])
    }
  }
  out <- list(x.lengths = x.cols, X = X)
  out
}

administrative <- function(formula, data, Weights = NULL) { # This is called by boolean() and processes formula
  temp.formula <- as.formula(formula[[1]])
  temp.formula[3] <- 1
  tt <- terms(temp.formula, data = data)
  if(!is.null(Weights) && !is.numeric(Weights)) Weights <- eval(parse(text = paste("data$", Weights, sep = "")))
#   mf <- model.frame(temp.formula, data = data, weights = Weights, na.action = na.pass)
  mf <- model.frame(temp.formula, data = data, na.action = na.pass)
  y <- model.response(mf, "numeric")
  Weights <- model.weights(mf)
  if(!is.null(Weights) && length(Weights) != length(y)) {
    cat("The weights vector must be of the same length as the other variables.\n")
    stop("Please respecify and call boolean() again.\n", call. = FALSE)
  }
  if(!is.null(Weights) && any(Weights <= 0)) {
    Weights[Weights <= 0] <- NA
    warning("Observations where weights <= 0 are dropped.")
  }
  if(is.null(Weights)) Weights <- rep(1, length(y))

  # Make big X matrix and deal with NAs in big X matrix
  old <- options("na.action")
  options(na.action = "na.pass")

  X <- parse.list(x = formula[-1], data)
  options(na.action = old[[1]])
  map <- X[[1]]
  X <- X[[2]]
  cc <- complete.cases(cbind(y, X, Weights))
  dropped <- sum(!cc)
  if(dropped > 0) { # Deal with rows that have NAs
    X <- X[cc,]
    y <- y[cc]
    Weights <- Weights[cc]
    cat(paste(dropped, "observations dropped due to missing data.\n", sep = " "))
  }

  if(unique(y)*1 != c(0,1) && unique(y)*1 != c(1,0)) {
    cat("Dependent variable must consist of some zeroes, some ones, and no other values.\n")
    stop("Please respecify and call boolean() again.\n", call. = FALSE)
  }


  out <- list(y = y, X = X, map = map, weights = Weights)
#   names(out)[1] <- formula$depvar
  return(out)
}

make.bool.fun <- function(formula) {
  structure.split <- strsplit(gsub(" ", "", formula$structure), split = NULL)[[1]]
  applied.formula <- as.character(NULL)
  count <- 1
  check <- c("&", "|")
  path <- as.character(NULL)
  for(i in (which(structure.split == "~") + 1):length(structure.split)) {
    if(structure.split[i] == "&") {
      if(structure.split[i-1] != ")") {
        applied.formula <- c(applied.formula, paste("x[", count, "]", sep = ""), "%&%")
        count <- count + 1
      }
      else applied.formula <- c(applied.formula, "%&%")
    }
    else if (structure.split[i] == "|") {
      if(structure.split[i-1] != ")") {
        applied.formula <- c(applied.formula, paste("x[", count, "]", sep = ""), "%|%")
        count <- count + 1
      }
      else applied.formula <- c(applied.formula, "%|%")
    }
    else if (structure.split[i] == "(") {
      applied.formula <- c(applied.formula, "(")
    }
    else if (structure.split[i] == ")") {
      if(structure.split[i - 1] != ")") {
        applied.formula <- c(applied.formula, paste("x[", count, "])", sep = ""))
        count <- count + 1
      }
      else applied.formula <- c(applied.formula, ")")
    }
    if(i == length(structure.split) && structure.split[i] != ")") {
      applied.formula <- c(applied.formula, paste("x[", count, "]", sep = ""))
    }
  }
  applied.formula <- gsub("[", "[,", applied.formula, fixed = TRUE)
  applied.formula <- c("switch((length(dim(x)) == 3) + 1,", applied.formula, ",",
                       gsub("]", ",]", applied.formula, fixed = TRUE))
  applied.formula <- paste(applied.formula, collapse = " ")
  eval(parse(text = paste("function(x) {", applied.formula, ")}", sep = "")))
}

llik <- function(par, NLM = FALSE, BHHH = FALSE, map = MAP, ...) {          # Calculates the log likelihood
  scobitMark <- ncol(X) + 1 # The alpha terms for scobit link functions (if any) start after the betas in par
  holder <- array(NA, dim = c(nrow(X), switch(is.matrix(par) + 1, nrow(map), c(nrow(map), nrow(par))))) # This holds probabilities
  for(i in 1:nrow(map)) {                                                   # Loop through the mini-GLMs
    SEQ <- map[i,1]:map[i,2]                                                # Mark relevant columns for the "ith" mini-GLM
    invlink <- make.invlink(link[i])                                        # Make the link function
    if(pmatch("scobit", link[i], nomatch = FALSE) > 0) {
      if(is.matrix(holder)) holder[,i] <- invlink(X[,SEQ] %*% par[SEQ], par[scobitMark])          # If necessary, use some type of scobit link,
      else holder[,i,] <- invlink(X[,SEQ] %*% t(par[,SEQ]), par[1,scobitMark])
      scobitMark <- scobitMark + 1                                          # and increment scobitMark
    }
    else {                                                                  # Invoke some one-parameter inverse link
      if(!is.matrix(par)) holder[,i] <- invlink(X[,SEQ] %*% par[SEQ])
      else holder[,i,] <- invlink(X[,SEQ] %*% t(par[,SEQ]))
    }
  } # END for loop

  probabilities <- bool.fun(holder)
  if(BHHH) return((y * mlog(probabilities) + (1 - y) * mlog(1 - probabilities))*weights)
  if(is.matrix(par)) out <- colSums((y * mlog(probabilities) + (1 - y) * mlog(1 - probabilities))*weights) ## CHECK
  else {
    out <- sum((y * mlog(probabilities) + (1 - y) * mlog(1 - probabilities))*weights)
    if(NLM) out <- out * -1
  }
  out
}


prob <- function(par, X, map = MAP, predicted = FALSE, ...) {
  scobitMark <- ncol(X) + 1
  holder <- array(NA, dim = c(nrow(X), switch(is.matrix(par) + 1, nrow(map), c(nrow(map), nrow(par)))))
  if(!predicted) {                                                         # If we just want probabilities, same as llik()
    for(i in 1:nrow(map)) {
      SEQ <- map[i,1]:map[i,2]
      invlink <- make.invlink(link[i])
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
      invlink <- make.invlink(link[i])
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
