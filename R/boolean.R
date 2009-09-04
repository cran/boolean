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

boolprep <- 
function(FORM, DEPVAR, ...) {
  if (!is.character(FORM)) FORM <- as.character(FORM)
  if(is.factor(DEPVAR) && nlevels(DEPVAR) == 2) depvar <- deparse(substitute(DEPVAR))
  if(all(is.numeric(DEPVAR)) | all(is.logical(DEPVAR))) depvar <- deparse(substitute(DEPVAR)) #as.character(match.call()$DEPVAR)
  else if(is.character(DEPVAR)) depvar <- DEPVAR
  else {
    cat("DEPVAR must be a binary variable or the name of a binary variable.\n")
    stop("Please respecify and call boolprep() again.\n", call. = FALSE)
  }
  form <- FORM
  form <- gsub(" ", "", form, extended = FALSE)
  form.split <- strsplit(form, split = NULL)[[1]]
  if(sum(form.split == "(") != sum(form.split == ")")) {
    cat(paste("Number of left parentheses does not match the number of right parentheses in ", FORM, "\n", sep = ""))
    stop("Please respecify and call boolprep() again\n", call. = FALSE)
  }
  check <- c("&", "|")
  for(i in 1:length(form.split)) {
    if(form.split[i] == "&" | form.split[i] == "|") {
      if(all(check != form.split[i])) {
        cat(paste("Insufficient number of nested parentheses in ", FORM, "\n", sep = ""))
        stop("Please check and call boolprep() again\n", call. = FALSE)
      }
      else check <- form.split[i]
    }
    else if(form.split[i] == ")" | form.split[i] == "(") check <- c("&", "|")
  }
  form.split <- form
  form.split <- gsub("(", "#", form.split, extended = FALSE)
  form.split <- gsub(")", "#", form.split, extended = FALSE)
  form.split <- gsub("&", "#", form.split, extended = FALSE)
  form.split <- gsub("|", "#", form.split, extended = FALSE)
  form.split <- strsplit(form.split, split = "#")[[1]]
  form.split <- form.split[form.split != ""]
  dots <- match.call(expand.dots = TRUE)[-(1:3)]
  if(!is.null(dots$constant) && is.logical(dots$constant)) {
    warning("See help(boolprep) for the new syntax")
    if(dots$constant == FALSE) stop("Please respecify and call boolprep() again\n", call. = FALSE)
    dots$constant <- NULL
  }
  if(length(dots) != length(form.split)) {
    cat(paste("Number of causal paths in", FORM, "does not match number of arguments in boolprep(...)\n", sep = " "))
    stop("Please respecify and call boolprep() again\n", call. = FALSE)
  }
  if(length(dots) == 1) {
    cat("You have entered a plain GLM instead of a boolean statement.\n")
    cat("If this is what you intended, use glm() or zelig()\n")
    stop("Otherwise, please respecify and call boolprep() again\n", call. = FALSE)
  }
  if(is.null(names(dots)) | any(names(dots) == "")) {
      cat("Some arguments in ... do not have names\n")
      cat("Make sure to use \nname_of_thing = ~ x1 + x2\n etc. when specifying arguments\n")
      stop("Please respecify and call boolprep() again\n", call. = FALSE)
  }

  structure <- paste(depvar, "~", FORM, sep = " ")
  out <- list(structure)
  for(i in 1:length(dots)) {
    if(names(dots)[i] != form.split[i]) {
      cat("The names of the formulas passed through the ... must match in spelling and order to the tokens in FORM\n")
      cat(paste(form.split[i], "and", names(dots)[i], "do not match\n", sep = " "))
      stop("Please respecify and call boolprep() again\n", call. = FALSE)
    }
    if(is.character(dots[[i]]) && all(strsplit(dots[[i]], split = NULL) != "~")) dots[[i]] <- paste("~", dots[[i]])
    out[[i + 1]] <- as.formula(dots[[i]])
  }
  names(out) <- c("structure", names(dots))
  out$depvar <- depvar
  class(out) <- "boolprep"
  out
}

boolean <- 
function(formula, data, link = "logit", start.values = NULL,
         method = "nlm", bootstrap = 0, control = list(fnscale = -1), weights = NULL, robust = FALSE, ...) {

  if(class(formula) != "boolprep") {
    cat("The formula argument must be an object created by boolprep()\n")
    stop("Please respecify and call boolean() again\n", call. = FALSE)
  }
  if(bootstrap < 0) {
    cat("bootstrap must be a non-negative integer\n")
    stop("Please respecify and call boolean() again\n", call. = FALSE)
  }
  if(bootstrap != as.integer(bootstrap)) {
    cat("bootstrap must be a non-negative integer\n")
    stop("Please respecify and call boolean() again\n", call. = FALSE)
  }
  if(bootstrap == 0 && length(method) == 2) {
    cat("Using a character vector of length *two* for method is sensible only if bootstrap > 0\n")
    cat("Either limit method to a single character string or specify bootsrap = 0\n")
    stop("Please respecify and call boolean() again\n", call. = FALSE)
  }
  if(!is.null(weights)) {
    warning("The weights argument is not currently supported")
    weights <- NULL
  }
  if(isTRUE(robust)) {
    warning("The robust argument is not currently supported")
    robust <- FALSE
  }

  if(missing(data)) data <- .GlobalEnv
  admin <- administrative(formula, data)
  y <- admin$y
  X <- admin$X
  MAP <- admin$map
  weights <- admin$weights

  if(length(link) == 1) link <- rep(link, nrow(MAP))
  if(length(link) != nrow(MAP)) {
    cat("The number of link functions must equal one or the number of tokens in the model\n")
    stop("Please respecify and call boolean() again\n", call. = FALSE)
  }

  npars <- ncol(X) + sum(pmatch(substr(link, 1, 6), "scobit", nomatch = FALSE, duplicates.ok = TRUE))
  if(is.null(start.values)) { # Based on plain GLMs
    start.values <- as.numeric(NULL)
    for(i in 2:length(formula)) {
      if(names(formula)[i] == "depvar") break
      temp.form  <- y ~ X[,MAP[i-1,1]:MAP[i-1,2]] -1
      temp.start <- coef(switch(link[i],
                                        "logit"   = glm(temp.form, binomial(link = "logit")),
                                        "probit"  = glm(temp.form, binomial(link = "probit")),
                                        "cloglog" = glm(temp.form, binomial(link = "cloglog")),
                                        "cauchit" = glm(temp.form, binomial(link = "cauchit")),
                                        "log"     = glm(temp.form, binomial(link = "log")),
                                                    glm(temp.form, binomial(link = "logit")))) ## covers scobit*
       start.values <- c(start.values, temp.start)
    }
    if(length(start.values) < npars) start.values <- c(start.values, rep(0, npars - length(start.values)))
  }
  else if(length(c(start.values)) != npars) {
    cat("Number of starting values does not equal the number of parameters\n")
    if(npars > ncol(X)) cat("Be sure to include starting values for the *log* of the ancillary parameters (e.g. zero)\n")
    stop("Please respecify and call boolean() again\n", call. = FALSE)
  }

  # Make the function to be evaluated
  bool.fun <- make.bool.fun(formula)

  # dots stuff
  dots <- list(...)

  # It is necessary to put everything into an environment, so that llik can "see" everything it needs
  boolean.env <- new.env()
  environment(llik) <- boolean.env
  environment(start.values) <- boolean.env
  environment(y) <- boolean.env
  environment(X) <- boolean.env
  environment(MAP) <- boolean.env
  environment(link) <- boolean.env
  environment(bool.fun) <- boolean.env
  environment(weights) <- boolean.env
  environment(dots) <- boolean.env

  method <- match.arg(method, c("Nelder-Mead", "BFGS", "CG", "SANN", "L-BFGS-B", "genoud", "nlm", "MH"), several.ok = TRUE)

  if(method[1] %in% c("Nelder-Mead", "BFGS", "CG", "SANN") | method[1] == "L-BFGS-B") {     # optim methods
    if(length(dots) > 0) {
      formals.temp <- formals(optim)
      for(i in 1:length(dots)) if(names(dots)[i] %in% names(formals.temp)) {
        formals.temp[which(names(formals.temp) == names(dots)[i])] <- dots[[i]]
      }
      formals(optim) <- formals.temp
    }
    CONTROL <- control
    CONTROL$fnscale <- -1
    METHOD <- method[1]
    out <- optim(start.values, llik, gr = NULL, hessian = TRUE, method = METHOD, control = CONTROL)
  }
  else if(method[1] == "constrOptim") {
    if(length(dots) > 0) {
      formals.temp <- formals(constrOptim)
      for(i in 1:length(dots)) if(names(dots)[i] %in% names(formals.temp)) {
        formals.temp[which(names(formals.temp) == names(dots)[i])] <- dots[[i]]
      }
      formals(constrOptim) <- formals.temp
    }
    CONTROL <- control
    CONTROL$fnscale <- -1
    out <- constrOptim(start.values, llik, grad = NULL, hessian = TRUE, method = "BFGS", control = CONTROL)
  }
  else if(method[1] == "genoud") {                      # Genetic algorithm optimization
    stopifnot(require(rgenoud))
    if(length(dots) > 0) {
      formals.temp <- formals(genoud)
      for(i in 1:length(dots)) if(names(dots)[i] %in% names(formals.temp)) {
        formals.temp[which(names(formals.temp) == names(dots)[i])] <- dots[[i]]
      }
      formals(genoud) <- formals.temp
    }
    out <- genoud(llik, npars, max = TRUE, starting.values = start.values, hessian = TRUE)
    out$counts <- out$generations
    out$generations <- NULL
    out$gradients <- NULL
    out$convergence <- NA
    out$message <- NA
  }
  else if(method[1] == "nlm") {                                  # If nlm() is used to *minimize* llik()
    if(length(dots) > 0) {
      formals.temp <- formals(nlm)
      for(i in 1:length(dots)) if(names(dots)[i] %in% names(formals.temp)) {
        formals.temp[which(names(formals.temp) == names(dots)[i])] <- dots[[i]]
      }
      formals(nlm) <- formals.temp
    }
    out <- nlm(llik, start.values, hessian = TRUE, NLM = TRUE)
    out$par <- out$estimate
    out$estimate <- NULL
    out$value <- out$minimum * -1
    out$minimum <- NULL
    out$hessian <- out$hessian * -1
    out$gradient <- NULL
    out$counts <- out$iterations
    out$iterations <- NULL
    out$convergence <- out$code
    out$code <- NULL
    out$message <- NA
  }
  else if(method[1] == "MH") {                                 # If Metropolis-Hastings is used to sample
    if(bootstrap > 0) {
      cat("Bootstapping and Metrolis-Hastings sampling are mutually exclusive\n")
      stop("Please respecify and call boolean() again.\n")
    }
    if(isTRUE(robust)) warning("Robust estimation of the variance-covariance is only possible with maximum-likelihood methods")
    stopifnot(require(MCMCpack))
    if(length(dots) > 0) {
      formals.temp <- formals(MCMCmetrop1R)
      for(i in 1:length(dots)) if(names(dots)[i] %in% names(formals.temp)) {
        formals.temp[which(names(formals.temp) == names(dots)[i])] <- dots[[i]]
      }
      formals(MCMCmetrop1R) <- formals.temp
    }
    CONTROL <- control
    CONTROL$fnscale <- -1
    if(length(dots) > 0) {
      formals.temp <- formals(optim)
      for(i in 1:length(dots)) if(names(dots)[i] %in% names(formals.temp)) {
        formals.temp[which(names(formals.temp) == names(dots)[i])] <- dots[[i]]
      }
      formals.temp$control <- CONTROL
      formals(optim) <- formals.temp
    }
    sample <- MCMCmetrop1R(llik, start.values, logfun = TRUE, y=admin$y, X=admin$X, MAP = admin$map)
    if(npars == ncol(X)) colnames(sample) <- colnames(X)
    else {
      varnames <- c(colnames(X), paste("log(alpha)_",
                     names(formula[-1])[pmatch(substr(link, 1, 6), "scobit", nomatch = FALSE, duplicates.ok = TRUE) == 1], sep = ""))
      names(out$par) <-  colnames(out$hessian) <- rownames(out$hessian) <- varnames
    }
    out <- list(coefficients = sample, variance = as.numeric(NA), convergence = as.numeric(NA), hessian = as.numeric(NA),
                message = as.numeric(NA), counts = attributes(sample)$burnin + attributes(sample)$mcmc)
  }
  else if(method[1] == "BHHH") {
    stop("'BHHH' method no longer supported")
#     if(length(dots) > 0) {
#       formals.temp <- formals(maxBHHH)
#       for(i in 1:length(dots)) if(names(dots)[i] %in% names(formals.temp)) {
#         formals.temp[which(names(formals.temp) == names(dots)[i])] <- dots[[i]]
#       }
#       formals(maxBHHH) <- formals.temp
#     }
#     out <- maxBHHH(llik, grad = NULL, hess = NULL, theta = start.values, BHHH = TRUE)
#     out$type <- NULL
#     out$gradient <- NULL
#     out$value <- out$maximum
#     out$maximum <- NULL
#     out$par <- out$estimate
#     out$estimate <- NULL
#     out$last.step <- NULL
#     out$activePar <- NULL
#     out$counts <- out$iterations
#     out$iterations <- NULL
#     out$convergence <- out$code
#     out$code <- NULL
  }
  if(method[1] != "MH") {
    attributes(out$par)$.Environment <- NULL
    if(npars == ncol(X)) names(out$par) <-  colnames(out$hessian) <- rownames(out$hessian) <- colnames(X)
    else {
      varnames <- c(colnames(X), paste("log(alpha)_",
                    names(formula[-1])[pmatch(substr(link, 1, 6), "scobit", nomatch = FALSE, duplicates.ok = TRUE) == 1], sep = ""))
      names(out$par) <-  colnames(out$hessian) <- rownames(out$hessian) <- varnames
    }
    out$coefficients <- out$par
    out$par <- NULL
    out$hessian <- out$hessian
    out$variance <- try(solve(-out$hessian), silent = TRUE)
#    if(isTRUE(robust)) ## Write this code
#    else
  }
  if(bootstrap > 0) {
    MLEs <- out$coefficients
    cat("The MLEs are:\n")
    print(cbind(MLEs))
    XX <- X
    yy <- y
    out$coefficients <- matrix(NA, nrow = bootstrap, ncol = length(MLEs))
    colnames(out$coefficients) <- names(MLEs)
    METHOD <- ifelse(length(method) == 1, method[1], method[2])
    if(length(dots) > 0) {
      formals.temp <- formals(optim)
      for(i in 1:length(dots)) if(names(dots)[i] %in% names(formals.temp)) {
        formals.temp[which(names(formals.temp) == names(dots)[i])] <- dots[[i]]
      }
      formals(optim) <- formals.temp
    }
    if(length(dots) > 0) {
      formals.temp <- formals(constrOptim)
      for(i in 1:length(dots)) if(names(dots)[i] %in% names(formals.temp)) {
        formals.temp[which(names(formals.temp) == names(dots)[i])] <- dots[[i]]
      }
      formals(constrOptim) <- formals.temp
    }
    if(length(dots) > 0) {
      formals.temp <- formals(genoud)
      for(i in 1:length(dots)) if(names(dots)[i] %in% names(formals.temp)) {
        formals.temp[which(names(formals.temp) == names(dots)[i])] <- dots[[i]]
      }
      formals(genoud) <- formals.temp
    }
    if(length(dots) > 0) {
      formals.temp <- formals(nlm)
      for(i in 1:length(dots)) if(names(dots)[i] %in% names(formals.temp)) {
        formals.temp[which(names(formals.temp) == names(dots)[i])] <- dots[[i]]
      }
      formals(nlm) <- formals.temp
    }
#     if(length(dots) > 0) {
#       formals.temp <- formals(maxBHHH)
#       for(i in 1:length(dots)) if(names(dots)[i] %in% names(formals.temp)) {
#         formals.temp[which(names(formals.temp) == names(dots)[i])] <- dots[[i]]
#       }
#       formals(maxBHHH) <- formals.temp
#     }
    CONTROL <- control
    CONTROL$fnscale <- -1

    for(i in 1:nrow(out$coefficients)) {
      sample.rows <- sample(1:nrow(XX), nrow(XX), replace = TRUE)
      X <- XX[c(sample.rows), ]
      y <- yy[c(sample.rows)]
      environment(X) <- boolean.env
      environment(y) <- boolean.env
      if(METHOD %in% c("Nelder-Mead", "BFGS", "CG", "SANN") | METHOD == "L-BFGS-B")
        out$coefficients[i,] <- optim(MLEs, llik, gr = NULL, hessian = FALSE, method = METHOD, control = CONTROL)$par
      else if (METHOD == "constrOptim")
        out$coefficients[i,] <- constrOptim(MLEs, llik, grad = NULL, hessian = FALSE, method = "BFGS", control = CONTROL)$par
      else if (METHOD == "genoud")
        out$coefficients[i,] <- genoud(llik, npars, max = TRUE, starting.values = MLEs, hessian = FALSE)$par
      else if (METHOD == "nlm")
        out$coefficients[i,] <- nlm(llik, MLEs, hessian = FALSE, NLM = TRUE)$estimate
#       else if (METHOD == "BHHH")
#         out$coefficients[i,] <- maxBHHH(llik, grad = NULL, hess = NULL, theta = MLEs)$estimate
      else {
        cat(paste(METHOD, "is an invalid option for the second element of method\n"))
        cat("Possible choices are: Nelder-Mead, BFGS, CG, L-BFGS-B, SANN, constrOptim, genoud, or nlm\n")
        stop("Please respecify and call boolean() again\n")
      }
      if(i %% 10 == 0) print(paste("Bootstrap iteration", i, "of", nrow(out$coefficients), sep = " "))
    }
  }
  out$boolean.call <- match.call(expand.dots = TRUE)
  out$boolean.call$formula <- eval(out$boolean.call$formula)
  out <- out[c("value", "coefficients", "variance", "hessian", "counts", "convergence", "message", "boolean.call")]
  if(method[1] == "MH") class(out) <- c("booltest", "MCMC")
  else if(bootstrap > 0) class(out) <- c("booltest", "BS")
  else class(out) <- c("booltest", "ML")
  return(out)
}

# Methods
coef.booltest <- 
function(object, ...) {
  if(class(object)[2] == "ML") return(object$coefficients)
  else if (class(object)[2] == "MCMC"){
    warning("These are the posterior means, which can only loosely be considered coefficients\n")
  }
  else if (class(object)[2] == "BS") {
    warning("These are the bootstrapped means, which can only loosely be considered coefficients\n")
  }
  return(colMeans(object$coefficients))
}

weights.booltest <- 
function(object, ...) {
    if(is.null(object$boolean.call$data)) object$boolean.call$data <- .GlobalEnv
  dta <- eval(object$boolean.call$data)
  admin <- administrative(object$boolean.call$formula, dta)
  weights <- admin$weights
  return(weights)
}

fitted.booltest <- fitted.values.booltest <- 
function(object, ...) {

  if(is.null(object$boolean.call$data)) object$boolean.call$data <- .GlobalEnv
  dta <- eval(object$boolean.call$data)
  admin <- administrative(object$boolean.call$formula, dta)
  y <- admin$y
  X <- admin$X
  MAP <- admin$map
  if(is.null(object$boolean.call$link)) link <- formals(boolean)$link
  else {
    link <- as.character(object$boolean.call$link)
    link <- link[link != "c"]
  }

  if(length(link) == 1) link <- rep(as.character(link), nrow(MAP))

  npars <- ncol(X) + sum(pmatch(substr(link, 1, 6), "scobit", nomatch = FALSE, duplicates.ok = TRUE))

  ## Make the function to be evaluated
  bool.fun <- make.bool.fun(eval(object$boolean.call$formula))

  # It is necessary to put everything into an environment, so that llik can "see" everything it needs
  boolean.env <- new.env()
  environment(prob) <- boolean.env
  environment(y) <- boolean.env
  environment(link) <- boolean.env
  environment(bool.fun) <- boolean.env

  ev <- prob(par = c(coef(object)), X, MAP)
  ev
}

summary.booltest <- 
function(object, ...) {
  out <- list(formula = object$boolean.call$formula)
  if(class(object)[2] == "MCMC")  {
    out$summary.mat <- summary(object$coefficients)
    out$log.likelihood <- NULL
    out$iterations <- NULL
  }
  else if(class(object)[2] == "BS") {
    mat <- as.matrix(cbind(Estimates = colMeans(object$coefficients), SEs = apply(object$coefficients, 2, FUN = sd)))
    mat <- cbind(mat, CI95Lo = apply(object$coefficients, 2, FUN = quantile, probs = .025),
                      CI95Hi = apply(object$coefficients, 2, FUN = quantile, probs = .975))
    out$summary.mat <- mat
    out$log.likelihood <- object$value
    out$iterations <- nrow(object$coefficients)
  }
  else { # ML models
    if(is.matrix(object$variance)) mat <- as.matrix(cbind(Estimates = object$coefficients, SEs = sqrt(diag(object$variance))))
    else mat <- as.matrix(cbind(Estimates = object$coefficients, SEs = as.numeric(NA)))
    mat <- cbind(mat, Zstat = mat[,1] / mat[,2])
    mat <- cbind(mat, ProbZ = pnorm(abs(mat[,3]), lower.tail = FALSE))
    mat <- cbind(mat, CI95Lo = mat[,1] - 1.96 * mat[,2], CI95Hi = mat[,1] + 1.96 * mat[,2])
    colnames(mat) <- c("Estimates", "SEs", "Zstat", "ProbZ", "CI95Lo", "CI95Hi")
    out$summary.mat <- mat
    out$log.likelihood <- object$value
    out$iterations <- object$counts[1]
  }
  class(out) <- "summary.booltest"
  return(out)
}

print.summary.booltest <- 
function(x, ...) {
  cat("\n")
  cat(x$formula$structure)
  cat("\nwhere\n")
  for(i in 2:length(x$formula)) cat(paste(names(x$formula)[i], x$formula[i], "\n", sep = " "))
  cat("\n")
  if(is.matrix(x$summary.mat)) {
    cat(paste("Log-likelihood = ", round(x$log.likelihood[[1]], 3), sep=""))
    cat("    ")
    cat(paste("Iterations = ", x$iterations, sep = ""))
    cat("\n")
  }
  print(x$summary.mat, 3)
}

latent <- 
function(object, probability = FALSE, invMills = FALSE) {
  if(isTRUE(probability) & isTRUE(invMills)) {
    cat("Either probability or invMills must be FALSE\n")
    stop("Please respecify and call latent again\n", call. = FALSE)
  }
  if(is.null(object$boolean.call$data)) object$boolean.call$data <- .GlobalEnv
  dta <- eval(object$boolean.call$data)
  admin <- administrative(object$boolean.call$formula, dta)
  y <- admin$y
  X <- admin$X
  map <- admin$map

  if(is.null(object$boolean.call$link)) link <- formals(boolean)$link
  else {
    link <- as.character(object$boolean.call$link)
    link <- link[link != "c"]
  }
  if(length(link) == 1) link <- rep(link, nrow(map))

  npars <- ncol(X) + sum(pmatch(substr(link, 1, 6), "scobit", nomatch = FALSE, duplicates.ok = TRUE))

  ## Make the function to be evaluated
  bool.fun <- make.bool.fun(object$boolean.call$formula)

  par <- object$coefficients
  scobitMark <- ncol(X) + 1
  holder <- array(NA, dim = c(nrow(X), switch(is.matrix(par) + 1, nrow(map), c(nrow(map), nrow(par)))))
  cnames <- names(object$boolean.call$formula)
  mark <- !(cnames %in% c("structure", "depvar"))
  cnames <- cnames[mark]
  colnames(holder) <- cnames
  for(i in 1:nrow(map)) {
    SEQ <- map[i,1]:map[i,2]
    invlink <- make.invlink(link[i])
    if(pmatch("scobit", link[i], nomatch = FALSE) > 0) {
      if(isTRUE(probability)) {
        if(is.matrix(holder)) holder[,i] <- invlink(X[,SEQ] %*% par[SEQ], par[scobitMark])
        else holder[,i,] <- invlink(X[,SEQ] %*% t(as.matrix(par[,SEQ])), par[1,scobitMark])
      }
      else if(isTRUE(invMills)) {
        if(is.matrix(holder)) holder[,i] <- NA
        else holder[,i,] <- NA
      }
      else {
        if(is.matrix(holder)) holder[,i] <- X[,SEQ] %*% par[SEQ]
        else holder[,i,] <- X[,SEQ] %*% t(as.matrix(par[,SEQ]))
      }
      scobitMark <- scobitMark + 1
    }
    else {
      if(isTRUE(probability)) {
        if(is.matrix(holder)) holder[,i] <- invlink(X[,SEQ] %*% par[SEQ])
        else holder[,i,] <- invlink(X[,SEQ] %*% t(as.matrix(par[,SEQ])))
      }
      else if(isTRUE(invMills)) {
        if(link[i] == "probit") {
          if(is.matrix(holder)) {
            eta <- X[,SEQ] %*% par[SEQ]
            holder[,i] <- dnorm(eta) / invlink(eta)
          }
          else {
            eta <- X[,SEQ,] %*% t(as.matrix(par[,SEQ]))
            holder[,i,] <- dnorm(eta) / invlink(eta)
          }
        }
        else {
          if(is.matrix(holder)) holder[,i] <- NA
          else holder[,i,] <- NA
          warning("Inverse Mills ratio is only estimatable for paths with a probit link")
        }
      }
      else {
        if(is.matrix(holder)) holder[,i] <- X[,SEQ] %*% par[SEQ]
        else holder[,i,] <- X[,SEQ] %*% t(par[,SEQ])
      }
    }
  }
  return(holder)
}

boolfirst <- 
function(...) {
  cat("The boolfirst() function has been deprecated and nothing is returned\n")
  cat("Please use the Zelig library to obtain first differences of probabilities\n")
}

boolprof <- 
function(object, gvar = NULL, range = NULL, M = 100) {
  if(is.null(object$boolean.call$data)) object$boolean.call$data <- .GlobalEnv
  dta <- eval(object$boolean.call$data)
  admin <- administrative(object$boolean.call$formula, dta)
  y <- admin$y
  X <- admin$X
  MAP <- admin$map

  if(is.null(object$boolean.call$link)) link <- formals(boolean)$link
  else {
    link <- as.character(object$boolean.call$link)
    link <- link[link != "c"]
  }
  if(length(link) == 1) link <- rep(as.character(link), nrow(MAP))

  npars <- ncol(X) + sum(pmatch(substr(link, 1, 6), "scobit", nomatch = FALSE, duplicates.ok = TRUE))

  ## Make the function to be evaluated
  bool.fun <- make.bool.fun(object$boolean.call$formula)

  # It is necessary to put everything into an environment, so that llik can "see" everything it needs
  boolean.env <- new.env()
  environment(llik) <- boolean.env
  environment(y) <- boolean.env
  environment(X) <- boolean.env
  environment(MAP) <- boolean.env
  environment(link) <- boolean.env
  environment(bool.fun) <- boolean.env
  weights <- eval(object$boolean.call$weights)
  if(length(weights) == 0) weights <- rep(1, nrow(X))
  environment(weights) <- boolean.env

  # Start unique part of boolprof() function
  if(is.null(range)) {
    if(class(object)[2] == "MCMC") {
      keep <- sample(1:nrow(object$coefficients), M, replace = FALSE, prob = NULL)
      simpar <- object$coefficients[c(keep),]
    }
    else if(class(object)[2] == "BS") {
      keep <- sample(1:nrow(object$coefficients), M, replace = FALSE, prob = NULL)
      simpar <- object$coefficients[c(keep),]
    }
    else {
      if(!is.matrix(object$variance)) {
        warning("Variance-covariance matrix indeterminate, consider speciying range, assuming plus or minus 1.0")
        simpar <- coef(object)
        simpar <- simpar + matrix(seq(from = -1, to = 1, length = M), nrow = nrow(simpar), ncol = M, byrow = TRUE)
      }
      else {
        vcv <- object$variance
        coeffs <- coef(object)
        simpar <- matrix(NA, nrow = M, ncol = length(coeffs))
        for(i in 1:ncol(simpar)) {
          SD <- sqrt(diag(vcv)[i])
          simpar[,i] <- coeffs[i] + seq(from = -2*SD, to = 2*SD, length = M)
        }
      }
    }
  }

  if(is.null(gvar)) {
    SIMPAR <- matrix(object$coefficients, nrow = M, ncol = npars, byrow = TRUE)
    devices <- ceiling(ncol(simpar)/9)
    i <- 1
    while(i <  devices) {
      get(getOption("device"))()
      i <- i + 1
    }
    par(mfrow = c(3,3))
    for(i in 1:ncol(simpar)) {
      if((i-1) %% 9 == 0) {
        devices <- devices - 1
        dev.set(devices)
        par(mfrow = c(3,3))
      }
      simpar.temp <- SIMPAR
      simpar.temp[,i] <- sort(simpar[,i])
      ll <- llik(par = simpar.temp)
      if(i <= ncol(X)) {
        plot(x = simpar.temp[,i], y = ll, type = "l", sub = colnames(simpar)[i],
             xlab = expression(hat(beta)), ylab = "Log-likelihood")
      }
      else {
        plot(x = exp(simpar.temp[,i]), y = ll, type = "l", sub = names(object$boolean.call$formula[-1])[i-ncol(X)-1], ## CHECK
             xlab = paste("ln(", expression(hat(alpha)), ")", sep = ""), ylab = "Log-likelihood")
      }
    }
    op <- par(mfrow = c(1,1))
  }
  else {
    simpar <- matrix(coef(object), nrow = length(range), ncol = npars, byrow = TRUE)
    colnames(simpar) <- rownames(coef(object))
    marker <- which(colnames(simpar) %in% gvar)
    if(!is.null(range)) {
      if(is.numeric(range)) simpar[,marker] <- range
      else {
        stopifnot(is.list(range))
        for(i in marker) simpar[,marker[i]] <- range[[i]]
      }
    }
    ll <- llik(par = simpar)
    par(mfrow = c(1,length(marker)))
    for(i in marker) {
      if(marker[i] <= ncol(X)) {
        plot(x = simpar[,marker], y = ll, type = "l",
               xlab = paste(expression(beta), "for", colnames(simpar)[marker], sep = " "), ylab = "Log-likelihood")
      }
      else {
        plot(x = exp(simpar[,marker]), y = ll, type = "l",
                          xlab = paste(expression(alpha), gvar, sep = " "),
             ylab = "Log-likelihood")
      }
    }
  }
}

boolplot <- 
function(z.out, variable, delta = 0, suppression.factor = FALSE, CI = 95,
         truehist = TRUE, legend = TRUE, plot.both = FALSE, polygon = FALSE, yscale = NULL, ...) {
  stopifnot(require(Zelig))
  if(class(z.out)[1] != "booltest") stop("boolplot() only works with boolean models")
  stopifnot(!is.null(z.out$zelig))
  if(delta == 0) {
        if(plot.both) stop("if plot.both = TRUE, then you must specify a nonzero delta to shift 'variable'")
  }
  else if(suppression.factor) stop("if suppression.factor = TRUE, delta must be zero")

  if(CI > 1) CI <- CI / 100

  dots <- list(...)
  x.out <- setx(z.out, fn = NULL)
  x.out.base <- setx(z.out)

  if(nrow(x.out.base) > 1) stopifnot(nrow(x.out.base) == nrow(x.out))
  x.out[,substr(colnames(x.out), start = 1, stop = nchar(variable)) != variable] <-
  x.out.base[,substr(colnames(x.out), start = 1, stop = nchar(variable)) != variable]

  if(length(dots) > 0) for(i in 1:length(dots))
    x.out[,substr(colnames(x.out), 1, nchar(names(dots)[i])) == names(dots)[i]] <- dots[[i]]

  xvar <- x.out[,as.logical(pmatch(substr(colnames(x.out), start = 1,
                 stop = nchar(variable)), table = variable, nomatch = FALSE))]
  x.out <- x.out[order(xvar),]
  xvar <- sort(xvar)
  xvar.unique <- !duplicated(xvar)
  xvar <- xvar[xvar.unique]
  x.out <- x.out[xvar.unique,]
  x.out.1 <- NULL
  if(any(delta != 0)) {
    x.out.1 <- x.out
    for(i in 1:ncol(x.out.1)) if(all(x.out.1[,i] == xvar)) x.out.1[,i] <- xvar + delta
  }
  s.out <- sim(z.out, x.out, x.out.1)
  if(isTRUE(suppression.factor)) s.out$qi$ev <- 1 - s.out$qi$ev
  if(all(delta == 0) | plot.both) {
    if(is.null(yscale)) yscale <- if(any(s.out$qi$fd < 0)) c(-1,1) else 0:1
    plot(x = xvar, y = rep(.5, length(xvar)), type = "n", ylim = yscale,
         xlab = variable, ylab = "")
    if(truehist) {
      stopifnot(require(MASS))
      par(new = TRUE)
      truehist(xvar, col = "white", border = "gray", xlab = "", ylab = "", ylim = yscale,  axes = FALSE)
    }
    if(!polygon) {
      lines(x = xvar, y = apply(s.out$qi$ev, 1, quantile, probs = (1 - CI) / 2), lty = "dashed")
      lines(x = xvar, y = apply(s.out$qi$ev, 1, quantile, probs = CI + (1 - CI) / 2), lty = "dashed")
    }
    else polygon(c(xvar, sort(xvar, decreasing = TRUE)), y = c(apply(s.out$qi$ev, 1, quantile, probs = CI + (1 - CI) / 2), rev(apply(s.out$qi$ev, 1, quantile, probs = (1 - CI) / 2))), col = "gray80", border = NA)
    lines(x = xvar, y = rowMeans(s.out$qi$ev), col = "red", pch = 20, type = if(truehist) "l" else "b")
    if(plot.both) {
      if(!polygon) {
        lines(x = xvar, y = apply(s.out$qi$fd, 1, quantile, probs = (1 - CI) / 2), lty = "dashed")
        lines(x = xvar, y = apply(s.out$qi$fd, 1, quantile, probs = CI + (1 - CI) / 2), lty = "dashed")
      }
      else polygon(c(xvar, sort(xvar, decreasing = TRUE)), y = c(apply(s.out$qi$fd, 1, quantile, probs = CI + (1 - CI) / 2), rev(apply(s.out$qi$fd, 1, quantile, probs = (1 - CI) / 2))), col = "gray80", border = NA)
      lines(x = xvar, y = rowMeans(s.out$qi$fd), col = "blue", pch = 20, lty = "dotted", type = if(truehist) "l" else "b")
      if(legend)
        legend(x = if(mean(s.out$qi$ev[1,]) < mean(s.out$qi$ev[nrow(s.out$qi$ev),]))
        "topleft" else "topright", legend = c(paste("Confidence interval (", CI, ")"),
        if(suppression.factor) "Expected value, Pr(y = 0)" else "Expected value, Pr(y = 1)",
        "Change in expected value"), lty = c(if(polygon) "solid" else "dashed", "solid", "dotted"), col = c(if(polygon) "gray80" else "black", "red", "blue"))
    }
    else if(legend)
      legend(x = if(mean(s.out$qi$ev[1,]) < mean(s.out$qi$ev[nrow(s.out$qi$ev),]))
         "topleft" else "topright", legend = c(paste("Confidence interval (", CI, ")"),
         if(suppression.factor) "Expected value, Pr(y = 0)" else "Expected value, Pr(y = 1)"), lty = c(if(polygon) "solid" else "dashed", "solid"), col = c(if(polygon) "gray80" else "black", "red"))
  }
  else {
    plot(x = xvar, y = rep(.5, length(xvar)), type = "n", ylim = 0:1,
         xlab = variable, ylab = "")
    if(truehist) {
      stopifnot(require(MASS))
      par(new = TRUE)
      truehist(xvar, col = "white", border = "gray", xlab = "", ylab = "", ylim = 0:1, axes = FALSE)
    }
    if(isTRUE(suppression.factor)) s.out$qi$fd <- 1 - s.out$qi$fd
    lines(x = xvar, y = apply(s.out$qi$fd, 1, quantile, probs = (1 - CI) / 2), lty = "dashed")
    lines(x = xvar, y = apply(s.out$qi$fd, 1, quantile, probs = CI + (1 - CI) / 2), lty = "dashed")
    lines(x = xvar, y = rowMeans(s.out$qi$fd), col = "blue", pch = 20, lty = "dotted", type = if(truehist) "l" else "b")
    if(legend) legend(x = if(mean(s.out$qi$ev[1,]) < mean(s.out$qi$ev[nrow(s.out$qi$ev),]))
         "topleft" else "topright", legend = c(paste("Confidence interval (", CI, ")"),
         if(suppression.factor) "Expected value, Pr(y = 0)" else "Expected value, Pr(y = 1)"), lty = c(if(polygon) "solid" else "dashed", "solid"), col = c(if(polygon) "gray80" else "black", "red"))
  }
}

effectplot <- 
function(z.out, variables, delta = 1, CI = 95, truehist = TRUE, 
                       legend = TRUE, polygon = FALSE, ...) {

  stopifnot(require(Zelig))
  if(class(z.out)[1] != "booltest") stop("boolplot() only works with boolean models")
  stopifnot(!is.null(z.out$zelig))
  if(CI > 1) CI <- CI / 100

  dots <- list(...)
  x.out <- setx(z.out, fn = NULL)
  x.out.base <- setx(z.out)

  if(nrow(x.out.base) > 1) stopifnot(nrow(x.out.base) == nrow(x.out))
  x.out[,substr(colnames(x.out), start = 1, stop = nchar(variables[1])) != variables[1]] <-
  x.out.base[,substr(colnames(x.out), start = 1, stop = nchar(variables[1])) != variables[1]]

  if(length(dots) > 0) for(i in 1:length(dots))
    x.out[,substr(colnames(x.out), 1, nchar(names(dots)[i])) == names(dots)[i]] <- dots[[i]]

  xvar <- x.out[,as.logical(pmatch(substr(colnames(x.out), start = 1,
                 stop = nchar(variables[1])), table = variables[1], nomatch = FALSE))]
  x.out <- x.out[order(xvar),]
  xvar <- sort(xvar)

  # Why was there no xvar.unique before?
  xvar.unique <- !duplicated(xvar)
  xvar <- xvar[xvar.unique]
  x.out <- x.out[xvar.unique,]

  x.out.1 <- x.out
  x.out.1[,substr(colnames(x.out), start = 1, stop = nchar(variables[1])) == variables[2]] <-
    x.out[,substr(colnames(x.out), start = 1, stop = nchar(variables[1])) == variables[2]] + delta

  s.out <- sim(z.out, x.out, x.out.1)
  yvar <- rowMeans(s.out$qi$fd)
  plot(x = xvar, y = yvar, type = "n", ylim = if(truehist) c(min(0, min(yvar)), 1) else NULL,
       xlab = variables[2], ylab = paste("Expected change in Pr(y) given a", delta, "unit change in", variables[1]))
  if(truehist) {
      stopifnot(require(MASS))
      par(new = TRUE)
      truehist(xvar, col = "white", border = "gray", xlab = "", ylab = "", ylim = c(min(0, min(yvar)), 1), axes = FALSE)
  }
  if(!polygon) {
      lines(x = xvar, y = apply(s.out$qi$fd, 1, quantile, probs = (1 - CI) / 2), lty = "dashed")
      lines(x = xvar, y = apply(s.out$qi$fd, 1, quantile, probs = CI + (1 - CI) / 2), lty = "dashed")
  }
  else polygon(c(xvar, sort(xvar, decreasing = TRUE)), y = c(apply(s.out$qi$fd, 1, quantile, probs = CI + (1 - CI) / 2), rev(apply(s.out$qi$fd, 1, quantile, probs = (1 - CI) / 2))), col = "gray80", border = NA)
  points(x = xvar, y = yvar, pch = 20, col = "red")
}


