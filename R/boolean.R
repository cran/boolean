## source("boolean.R")

boolprep <- function(FORM, DEPVAR, ...) {
  if (!is.character(FORM)) FORM <- as.character(FORM)
  if(all(is.numeric(DEPVAR)) | all(is.logical(DEPVAR))) depvar <- as.character(match.call()$DEPVAR)
  else if(is.character(DEPVAR)) depvar <- DEPVAR
  else {
    cat("DEPVAR must be a binary variable or the name of a binary variable.\n")
    stop("Please respecify and call boolprep() again.\n", call. = FALSE)
  }
  form <- FORM
  form <- gsub(" ", "", form, extended = FALSE)
  form.split <- strsplit(form, split = NULL)[[1]]
  if(sum(form.split == "(") != sum(form.split == ")")) {
    cat(paste("Number of left parentheses does not match the number of right parentheses in ", FORM, ".\n", sep = ""))
    stop("Please respecify and call boolprep() again.\n", call. = FALSE)
  }
  check <- c("&", "|")
  for(i in 1:length(form.split)) {
    if(form.split[i] == "&" | form.split[i] == "|") {
      if(all(check != form.split[i])) {
        cat(paste("Insufficient number of nested parentheses in ", FORM, ".\n", sep = ""))
        stop("Please check and call boolprep() again.\n", call. = FALSE)
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
    if(dots$constant == FALSE) stop("Please respecify and call boolprep() again.\n", call. = FALSE)
    dots$constant <- NULL
  }
  if(length(dots) != length(form.split)) {
    cat(paste("Number of causal paths in", FORM, "does not match number of arguments in boolprep(...).\n", sep = " "))
    stop("Please respecify and call boolprep() again.\n", call. = FALSE)
  }
  if(length(dots) == 1) {
    cat("You have entered a plain GLM instead of a boolean statement.\n")
    cat("If this is what you intended, use glm() or zelig().\n")
    stop("Otherwise, please respecify and call boolprep() again.\n", call. = FALSE)
  }
  if(is.null(names(dots)) | any(names(dots) == "")) {
      cat("Some arguments in ... do not have names.\n")
      cat("Make sure to use 'name_of_thing = ~ x1 + x2' etc. when specifying arguments.\n")
      stop("Please respecify and call boolprep() again.\n", call. = FALSE)
  }

  structure <- paste(depvar, "~", FORM, sep = " ")
  out <- list(structure)
  for(i in 1:length(dots)) {
    if(names(dots)[i] != form.split[i]) {
      cat("The names of the formulas passed through the ... must match in spelling and order to the tokens in FORM.\n")
      cat(paste(form.split[i], "and", names(dots)[i], "do not match.\n", sep = " "))
      stop("Please respecify and call boolprep() again.\n", call. = FALSE)
    }
    if(is.character(dots[[i]]) && all(strsplit(dots[[i]], split = NULL) != "~")) dots[[i]] <- paste("~", dots[[i]])
    out[[i + 1]] <- as.formula(dots[[i]])
  }
  names(out) <- c("structure", names(dots))
  class(out) <- "boolprep"
  out
}

boolean <- function(formula, data, link = "logit", start.values = NULL,
                    method = "nlm", bootstrap = 0, control = list(fnscale = -1), weights = NULL, robust = FALSE, ...) {
  if(class(formula) != "boolprep") {
    cat("The formula argument must be an object created by boolprep().\n")
    stop("Please respecify and call boolean() again.\n", call. = FALSE)
  }
  if(bootstrap < 0) {
    cat("bootstrap must be a non-negative integer.\n")
    stop("Please respecify and call boolean() again.\n", call. = FALSE)
  }
  if(bootstrap != as.integer(bootstrap)) {
    cat("bootstrap must be a non-negative integer.\n")
    stop("Please respecify and call boolean() again.\n", call. = FALSE)
  }
  if(bootstrap == 0 && length(method) == 2) {
    cat("Using a character vector of length *two* for method is sensible only if bootstrap > 0.\n")
    cat("Either limit method to a single character string or specify bootsrap = 0.\n")
    stop("Please respecify and call boolean() again.\n", call. = FALSE)
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
  y <- admin[[1]]
  X <- admin[[2]]
  MAP <- admin[[3]]
  weights <- admin[[4]]

  if(length(link) == 1) link <- rep(link, nrow(MAP))
  if(length(link) != nrow(MAP)) {
    cat("The number of link functions must equal one or the number of tokens in the model.\n")
    stop("Please respecify and call boolean() again.\n", call. = FALSE)
  }

  npars <- ncol(X) + sum(pmatch(substr(link, 1, 6), "scobit", nomatch = FALSE, duplicates.ok = TRUE))
  if(is.null(start.values)) { # Based on plain GLMs
    start.values <- as.numeric(NULL)
    for(i in 2:length(formula)) {
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
    cat("Number of starting values does not equal the number of parameters.\n")
    if(npars > ncol(X)) cat("Be sure to include starting values for the *log* of the ancillary parameters (e.g. zero).\n")
    stop("Please respecify and call boolean() again.\n", call. = FALSE)
  }

  ## Make the function to be evaluated
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

  if(method[1] %in% c("Nelder-Mead", "BFGS", "CG", "SANN") | method[1] == "L-BFGS-B") {                              # optim methods
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
    sample <- MCMCmetrop1R(llik, start.values, logfun = TRUE, y=admin[[1]], X=admin[[2]], MAP = admin[[3]])
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
    stopifnot(require(micEcon))
    if(length(dots) > 0) {
      formals.temp <- formals(maxBHHH)
      for(i in 1:length(dots)) if(names(dots)[i] %in% names(formals.temp)) {
        formals.temp[which(names(formals.temp) == names(dots)[i])] <- dots[[i]]
      }
      formals(maxBHHH) <- formals.temp
    }
    out <- maxBHHH(llik, grad = NULL, hess = NULL, theta = start.values, BHHH = TRUE)
    out$type <- NULL
    out$gradient <- NULL
    out$value <- out$maximum
    out$maximum <- NULL
    out$par <- out$estimate
    out$estimate <- NULL
    out$last.step <- NULL
    out$activePar <- NULL
    out$counts <- out$iterations
    out$iterations <- NULL
    out$convergence <- out$code
    out$code <- NULL
  }
  else {
    cat("Method must be a character string or character vector.\n")
    cat("Possible choices are: 'Nelder-Mead', 'BFGS', 'CG', 'L-BFGS-B', 'SANN', 'constrOptim', 'genoud', 'nlm', 'BHHH' or 'MH'.\n")
    stop("Please respecify and call boolean() again.\n")
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
    out$hessian <- out$hessian / mean(weights) ## CHECK
    if(isTRUE(robust)) {
      score <- matrix(NA, nrow = nrow(X), ncol = length(out$coefficients))
      for(i in 1:ncol(score)) score[,i] <- (llik(out$coefficients + 0.5 * 10^-6 * (1:ncol(score) == i), BHHH = TRUE)  - 
                                            llik(out$coefficients - 0.5 * 10^-6 * (1:ncol(score) == i), BHHH = TRUE)) / 10^-6
      B <- matrix(0, nrow = length(out$coefficients), ncol = length(out$coefficients))
      for(i in 1:nrow(score)) B <- B + outer(score[i,], score[i,]) * weights[i] ## CHECK
      B <- B / sum(weights) ## CHECK
      out$variance <- try(solve(out$hessian) %*% B %*% solve(out$hessian), silent = TRUE)
    }
    else out$variance <- try(solve(-out$hessian), silent = TRUE)
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
    if(length(dots) > 0) {
      formals.temp <- formals(maxBHHH)
      for(i in 1:length(dots)) if(names(dots)[i] %in% names(formals.temp)) {
        formals.temp[which(names(formals.temp) == names(dots)[i])] <- dots[[i]]
      }
      formals(maxBHHH) <- formals.temp
    }
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
      else if (METHOD == "BHHH")
        out$coefficients[i,] <- maxBHHH(llik, grad = NULL, hess = NULL, theta = MLEs)$estimate
      else {
        cat(paste(METHOD, "is an invalid option for the second element of method.\n"))
        cat("Possible choices are: 'Nelder-Mead', 'BFGS', 'CG', 'L-BFGS-B', 'SANN', 'constrOptim', 'genoud', 'nlm', or 'BHHH'.\n")
        stop("Please respecify and call boolean() again.\n")
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
coef.booltest <- function(object, ...) {
  if(class(object)[2] == "ML") return(cbind(object$coefficients))
  else if (class(object)[2] == "MCMC"){
    warning("These are the posterior means, which can only loosely be considered 'coefficients'.\n")
  }
  else if (class(object)[2] == "BS") {
    warning("These are the bootstrapped means, which can only loosely be considered 'coefficients'.\n")
  }
  return(cbind(colMeans(object$coefficients)))
}

weights.booltest <- function(object, ...) {
    if(is.null(object$boolean.call$data)) object$boolean.call$data <- .GlobalEnv
  dta <- eval(object$boolean.call$data)
  admin <- administrative(object$boolean.call$formula, dta)
  weights <- admin[[4]]
  return(weights)
}

fitted.booltest <- fitted.values.booltest <- function(object, ...) {

  if(is.null(object$boolean.call$data)) object$boolean.call$data <- .GlobalEnv
  dta <- eval(object$boolean.call$data)
  admin <- administrative(object$boolean.call$formula, dta)
  y <- admin[[1]]
  X <- admin[[2]]
  MAP <- admin[[3]]
  if(is.null(object$boolean.call$link)) link <- "logit"
  else link <- object$boolean.call$link
  if(length(link) == 1) link <- rep(as.character(link), nrow(admin[[3]]))
  else link <- as.character(link)[-1]

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

summary.booltest <- function(object, ...) {
  out <- list(formula = object$boolean.call$formula)
  if(class(object)[2] == "MCMC")  out$summary.mat <- summary(object$coefficients)
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

print.summary.booltest <- function(x, ...) {
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

latent <- function(object, probability = FALSE, invMills = FALSE) {
  if(isTRUE(probability) & isTRUE(invMills)) {
    cat("Either probability or invMills must be FALSE.\n")
    stop("Please respecify and call latent again.\n", call. = FALSE)
  }
  if(is.null(object$boolean.call$data)) object$boolean.call$data <- .GlobalEnv
  dta <- eval(object$boolean.call$data)
  admin <- administrative(object$boolean.call$formula, dta)
  y <- admin[[1]]
  X <- admin[[2]]
  map <- admin[[3]]

  if(is.null(object$boolean.call$link)) link <- formals(boolean)$link
  else link <- object$boolean.call$link
  if(length(link) == 1) link <- rep(link, nrow(admin[[3]]))

  npars <- ncol(X) + sum(pmatch(substr(link, 1, 6), "scobit", nomatch = FALSE, duplicates.ok = TRUE))

  ## Make the function to be evaluated
  bool.fun <- make.bool.fun(object$boolean.call$formula)

  par <- object$coefficients
  scobitMark <- ncol(X) + 1 
  holder <- array(NA, dim = c(nrow(X), switch(is.matrix(par) + 1, nrow(map), c(nrow(map), nrow(par)))))
  colnames(holder) <- names(object$boolean.call$formula)[-1]
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
          warning("Inverse Mills ratio is only estimatable for paths with a probit link.")
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

boolfirst <- function(...) {
  cat("The boolfirst() function has been deprecated and nothing is returned.\n")
  cat("Please use the Zelig library to obtain first differences of probabilities.\n")
}

boolprof <- function(object, gvar = NULL, range = NULL, M = 100) {
  if(is.null(object$boolean.call$data)) object$boolean.call$data <- .GlobalEnv
  dta <- eval(object$boolean.call$data)
  admin <- administrative(object$boolean.call$formula, dta)
  y <- admin[[1]]
  X <- admin[[2]]
  MAP <- admin[[3]]

  if(is.null(object$boolean.call$link)) link <- formals(boolean)$link
#  else link <- as.character(object$boolean.call$link)[-1]
  else {
    link <- as.character(object$boolean.call$link)
    link <- link[link != "c"]
  }
  if(length(link) == 1) link <- rep(as.character(link), nrow(admin[[3]]))

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
      c.vcv <- chol(object$variance)
      mat <- matrix(rnorm(ncol(c.vcv)*M), nrow = ncol(c.vcv), ncol = M)
      simpar <- t(c.vcv %*% mat + object$coefficients)
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
    if(is.null(range)) {
      cat("If gvar is specified, a range must be specified for that variable.\n")
      stop("Please respecify and call boolprof() again.\n", call. = FALSE)
    }
    simpar <- matrix(coef(object), nrow = length(range), ncol = npars, byrow = TRUE)
    colnames(simpar) <- rownames(coef(object))
    marker <- which(colnames(simpar) == gvar)
    simpar[,marker] <- range
    ll <- llik(par = simpar)
    if(marker <= ncol(X)) {
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
