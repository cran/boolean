# source("hidden.R")

# Helper functions, consider moving these two
"%&%" <- function(x,y) {x * y}                 ## The AND function
"%|%" <- function(x,y) {1 - (1 - x) * (1 - y)} ## The OR function

mlog <- function(x) {                                       # Takes the log, but won't return -Inf (might test on a Mac G5)
  ifelse(x < .Machine$double.xmax^-1.0498, -744.44, log(x)) # -744.44 looks like the smallest R can go without resorting to -Inf
}

make.invlink <- function (link) 
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
  mf <- model.frame(temp.formula, data = data, na.action = na.pass)
  y <- model.response(mf, "numeric")
#   Weights <- model.weights(mf)
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
  return(out)
}

make.bool.fun <- function(formula) {
  structure.split <- strsplit(gsub(" ", "", formula$structure, extended = FALSE), split = NULL)[[1]]
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
  applied.formula <- gsub("[", "[,", applied.formula, extended = FALSE, fixed = TRUE)
  applied.formula <- c("switch((length(dim(x)) == 3) + 1,", applied.formula, ",", 
                       gsub("]", ",]", applied.formula, extended = FALSE, fixed = TRUE))
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

"check.mcmc.parameters" <-
  function(burnin, mcmc, thin) { # This is a MCMCpack error-checking function
  
    if(mcmc %% thin != 0) {
      cat("Error: MCMC iterations not evenly divisible by thinning interval.\n")
      stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
    }
    if(mcmc < 0) {
      cat("Error: MCMC iterations negative.\n")
      stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE) 
    }
    if(burnin < 0) {
      cat("Error: Burnin iterations negative.\n")
      stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
    }
    if(thin < 0) {
      cat("Error: Thinning interval negative.\n")
      stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
    }    
    return(0)
  }


"vector.tune" <- function(mcmc.tune, K){ # This is a MCMCpack function that processes the tuning parameters
  if (max(is.na(mcmc.tune))){
    cat("Error: Vector tuning parameter cannot contain NAs.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)    
  }
  if (length(mcmc.tune) == 1){
    mcmc.tune <- rep(mcmc.tune, K)
  }
  if (length(mcmc.tune) != K){
    cat("Error: length(vector tuning parameter) != length(theta) or 1.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if (sum(mcmc.tune <= 0) != 0) {
    cat("Error: Vector tuning parameter cannot contain negative values.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if (length(mcmc.tune)==1){
    return(matrix(mcmc.tune, 1, 1))
  }
  else{
    return(diag(as.double(mcmc.tune)))
  }
}

check.offset <- function (args)
{
    if (sum(names(args) == "offset") == 1) {
        cat("Error: Offsets are currently not supported in MCMCpack.\n")
        stop("Please respecify and call ", calling.function(),
            " again.\n", call. = FALSE)
    }
    return(0)
}

"form.seeds" <-
   function(seed) {                      # This is a MCMCpack function that deals with the RNG seeds
      if(length(seed)==1) {
         if(is.na(seed)) seed <- 12345
         seed <- as.integer(seed)
         if(seed < 0) {
            cat("Error: Mersenne seed negative.\n")
            stop("Please respecify and call ", calling.function(), " again.",
              call.=FALSE)                       
         }
         seeds <- list(0, rep(seed,6), 0)
      }
      if(length(seed)==2) {
         if(!is.list(seed)) {
            cat("Error: List must be passed to use L'Ecuyer.\n")
            stop("Please respecify and call ", calling.function(), " again.",
              call.=FALSE)          
         }
         lec.seed <- seed[[1]]
         lec.substream <- as.integer(seed[[2]])
         if(is.na(lec.seed[1])) lec.seed <- rep(12345, 6)
         if(length(lec.seed) != 6) {
            cat("Error: L'Ecuyer seed not of length six.\n")
            stop("Please respecify and call ", calling.function(), " again.",
              call.=FALSE)          
         }
         if(!all(lec.seed >= 0))  {
             cat("Error: At least one L'Ecuyer seed negative.\n")
            stop("Please respecify and call ", calling.function(), " again.",
              call.=FALSE)          
         }
         if( max(lec.seed[1:3]) >= 4294967087){
           cat("Error: At least one of first three L'Ecuyer seeds\n")
           cat("  greater than or equal to 4294967087\n")
           stop("Please respecify and call ", calling.function(), " again.",
                call.=FALSE)          
         }
         if( all(lec.seed[1:3]) == 0 ){
           cat("Error: first three L'Ecuyer seeds == 0\n")
           stop("Please respecify and call ", calling.function(), " again.",
                call.=FALSE)          
         }
         if( max(lec.seed[4:6]) >= 4294944443){
           cat("Error: At least one of last three L'Ecuyer seeds\n")
           cat("  greater than or equal to 4294944443\n")
           stop("Please respecify and call ", calling.function(), " again.",
                call.=FALSE)          
         }         
         if( all(lec.seed[4:6]) == 0 ){
           cat("Error: last three L'Ecuyer seeds == 0\n")
           stop("Please respecify and call ", calling.function(), " again.",
                call.=FALSE)          
         }
         if(lec.substream < 1) {
            cat("Error: L'Ecuyer substream number not positive.\n")
            stop("Please respecify and call ", calling.function(), " again.",
                 call.=FALSE)               
         }
         seeds <- list(1, lec.seed, lec.substream) 
      }
      if(length(seed)>2) {
            cat("Error: Seed passed as length greater than two.\n")
            stop("Please respecify and call ", calling.function(), " again.",
              call.=FALSE)        
      }
      return(seeds)
   }
