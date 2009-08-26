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

zelig2boolean <- function(formula, model, data, M, ...) {
#   require(boolean)
  mf <- match.call(expand.dots = TRUE)
  mf$model <- mf$M <- NULL
  mf[[1]] <- as.name("boolean.wrapper")
  as.call(mf)
}

boolean.wrapper <- function(formula, data, ...) {
  kall <- match.call(expand.dots = TRUE)
  kall.copy <- kall
  kall.copy[[1]] <- as.name("boolean")

  structure <- gsub(" ", "", formula[[1]])
  structure.split <- strsplit(structure, split = NULL)[[1]]
  Y.NAME <- as.name(substr(structure, 1, which(structure.split == "~") -1))

  formula.copy <- formula
  formula <- paste(Y.NAME, "~", deparse(formula.copy[[2]][[2]]))
  for(i in 3:length(formula.copy)) if(names(formula.copy)[i] == "depvar") break
                                   else formula <- paste(formula, "+", deparse(formula.copy[[i]][[2]]))

  formula <- as.formula(formula)

define.data <- function(formula, data, x.name = "x",
                        y.name = "y", mf.name = "mf",
                        xlev.name = "xlev", lev.name = "lev") { # function was part of previous versions of Zelig
  mc <- match.call()
  mc$x.name <- mc$y.name <- mc$t.name <- mc$xlev.name <- mc$l.name <- NULL

  mc[[1]] <- as.name("model.frame")
  mf <- eval.parent(mc)
  assign(mf.name, mf, envir = parent.frame())

  Terms <- attr(mf, "terms")

  x <- model.matrix(Terms, mf)
  assign(x.name, x, envir = parent.frame())

  yvar <- attr(Terms, "response")
  xvars <- as.character(attr(Terms, "response"))[-1]
  xlev <- if (length(xvars) > 0) {
    xlev <- lapply(mf[xvars], levels)
    xlev[!sapply(xlev, is.null)]
  }
  assign(xlev.name, xlev, envir = parent.frame())

  y <- model.response(mf)
  assign(y.name, y, envir = parent.frame())

  lev <- levels(y)
  assign(lev.name, lev, envir = parent.frame())

  attr(Terms, "items") <- list(x = x.name, y = y.name,
                               mf = mf.name, xlev =
                               xlev.name, lev = lev.name)
  Terms
}

  Terms <- define.data(formula, data)

  formula <- formula.copy

  fit <- eval(kall.copy)

  lev <- eval(as.name(attr(Terms, "items")$lev))
  kall <- eval(as.name("kall"))
  xlev <- eval(as.name(attr(Terms, "items")$xlev))
  mf <- eval(as.name(attr(Terms, "items")$mf))
  x <- eval(as.name(attr(Terms, "items")$x))

  fit$xlevels <- xlev
  fit$lev <- lev
  fit$terms <- Terms


  attr(fit, "na.message") <- attr(mf, "na.message")

  if (!is.null(attr(mf, "na.action")))
    fit$na.action <- attr(mf, "na.action")

  fit
}

zelig.boolprep <-
function(formula, model, data, by, save.data = FALSE, cite, ...) {
    fn1 <- paste("zelig2", model, sep = "")
    fn2 <- paste("zelig3", model, sep = "")
    if (!exists(fn1))
        stop(model, " not supported. Type help.zelig(\"models\") to list supported models.")
    mf <- zelig.call <- match.call(expand.dots = TRUE)
    zelig.call[[1]] <- as.name("zelig")
    if (missing(by))
        by <- NULL
    N <- M <- 1
    object <- list()
    if ("mi" %in% class(data))
        M <- length(data)
    if (M > 1)
        dat <- data[[1]]
    else dat <- data
    if (!is.null(by)) {
        if (any(as.character(by) %in% c(formula[[2]], formula[[3]])))
            stop("the variable selected for subsetting cannot be called in the formula.")
        idx <- dat[, by]
        mf$by <- NULL
        lev <- sort(unique(idx))
        N <- length(lev)
    }
    mf <- do.call(fn1, list(formula, model, dat, N, ...))
    for (i in 1:N) {
        if (N > 1) {
            dat <- list()
            if (M > 1) {
                for (j in 1:M) dat[[j]] <- data[[j]][idx == lev[i],
                  ]
            }
            else dat <- data[idx == lev[i], ]
        }
        else dat <- data
        obj <- list()
        for (j in 1:M) {
            if (M > 1)
                d <- dat[[j]]
            else d <- dat
            if (is.data.frame(d)) {
#                 d <- d[complete.cases(model.frame(mf$formula,
#                   data = d, na.action = na.pass)), ]
                mf$data <- d
#                 res <- create.ZeligS4(eval(as.call(mf)))
		res <- eval(as.call(mf))
                if (exists(fn2))
                  res <- do.call(fn2, list(res = res, fcall = mf,
                    zcall = as.list(zelig.call)))
                res$call <- as.call(zelig.call)
                if (save.data)
                  res$zelig.data <- d
                res$zelig <- model
                if (M > 1)
                  obj[[j]] <- res
                else obj <- res
            }
        }
        if (M > 1)
            class(obj) <- "MI"
        if (N > 1)
            object[[i]] <- obj
        else object <- obj
    }
    if (N > 1) {
        class(object) <- "strata"
        names(object) <- lev
    }
    return(object)
}

setx.booltest <-
function (object, fn = list(numeric = mean, ordered = median, other = mode),
          data = NULL, cond = FALSE, counter = NULL, ...)  # good
{
    mc <- match.call()
    if (class(object)[1] == "MI")
        object <- object[[1]]
    mode <- function(x) {
        tb <- tapply(x, x, length)
        if (is.factor(x))
            value <- factor(unlist(labels(tb[seq(along = tb)[tb ==
                max(tb)]])), levels = levels(x))
        else if (is.logical(x))
            value <- as.logical(unlist(labels(tb[seq(along = tb)[tb ==
                max(tb)]])))
        else if (is.character(x))
            value <- as.character(unlist(labels(tb[seq(along = tb)[tb ==
                max(tb)]])))
        else stop(paste(vars[i], "is not a supported variable type."))
        if (length(value) > 1) {
            warning("There is more than one mode. The first level is selected.")
            value <- sort(value)[1]
        }
        return(value)
    }
    median.default <- median
    median <- function(x) {
        if (is.numeric(x))
            value <- median.default(x)
        else if (is.ordered(x))
            value <- factor(levels(x)[median.default(as.integer(x))],
                levels = levels(x))
        else stop("median cannot be calculated for this data type")
        return(value)
    }
    max.default <- max
    max <- function(x, na.rm = FALSE) {
        if (is.numeric(x))
            value <- max.default(x, na.rm = na.rm)
        else if (is.ordered(x))
            value <- factor(levels(x)[length(levels(x))], levels = levels(x))
        else stop("max cannot be calculated for this data type")
        return(value)
    }
    min.default <- min
    min <- function(x, na.rm = FALSE) {
        if (is.numeric(x))
            value <- min.default(x, na.rm = na.rm)
        else if (is.ordered(x))
            value <- factor(levels(x)[1], levels = levels(x))
        else stop("min cannot be calculated for this data type")
        return(value)
    }
    tt <- terms(object)
    tt.attr <- attributes(tt)
    env <- tt.attr$.Environment
    if (is.null(env))
        env <- parent.frame()
    if (is.null(data))
        if (nrow(as.data.frame(object$zelig.data)) > 0)
            dta <- object$data
        else dta <- eval(object$call$data, envir = env)
    else dta <- as.data.frame(data)
    mf <- model.frame(tt, data = dta, na.action = na.pass)
    if (any(class(tt) == "multiple"))
        vars <- unlist(c(attr(tt, "depVars"), attr(tt, "indVars")),
            use.names = FALSE)
    else vars <- all.vars(tt)
    if (!is.null(tt.attr$response) && tt.attr$response)
        resvars <- all.vars(tt.attr$variables[[1 + tt.attr$response]])
    else resvars <- NULL

    # Start modifications
    admin <- administrative(object$boolean.call$formula, dta)
    # dta <- as.data.frame(cbind(admin[[2]], admin[[1]]))
    # vars <- c(colnames(admin[[2]]), "vote")
    data <- as.data.frame(admin[[2]])
    # End Modifications


    #data <- dta[, names(dta) %in% vars, drop = FALSE]
    if (!is.null(counter)) {
        if (!any(counter == vars))
            stop("the variable specified for counter is not used in the model")
        treat <- data[, names(data) == counter]
        if (is.numeric(treat)) {
            data[treat == 1, names(data) == counter] <- 0
            data[treat == 0, names(data) == counter] <- 1
        }
        else if (is.factor(treat)) {
            lev <- levels(treat)
            if (length(lev) == 2) {
                treat <- as.numeric(treat) - 1
                data[treat == 1, names(data) == counter] <- lev[1]
                data[treat == 0, names(data) == counter] <- lev[2]
            }
            else stop("counter only takes a binary variable")
        }
        else if (is.logical(treat)) {
            treat <- as.numeric(treat)
            data[treat == 1, names(data) == counter] <- FALSE
            data[treat == 0, names(data) == counter] <- TRUE
        }
        else stop("not supported variable type for counter")
        if (!cond)
            stop("if counter is specified, cond must be TRUE")
    }
    if (cond) {
        if (is.null(data))
            stop("if cond = TRUE, you must specify the data frame.")
        if (is.null(mc$fn))
            fn <- NULL
        if (!is.null(fn)) {
            warning("when cond = TRUE, fn is coerced to NULL")
            fn <- NULL
        }
        maxl <- nrow(data)
    }
    else if (!is.null(fn)) {
        if (is.null(fn$numeric) || !is.function(fn$numeric)) {
            warning("fn$numeric coerced to mean().")
            fn$numeric <- mean
        }
        if (is.null(fn$ordered) || !is.function(fn$ordered) ||
            identical(mean, fn$ordered)) {
            warning("fn$ordered coreced to median().")
            fn$ordered <- median
        }
        else if (identical(min.default, fn$ordered))
            fn$ordered <- min
        else if (identical(max.default, fn$ordered))
            fn$ordered <- max
        else if (identical(median.default, fn$ordered))
            fn$ordered <- median
        if (is.null(fn$other) || !is.function(fn$other)) {
            warning("the only available fn for other is mode.")
            fn$other <- mode
        }
        for (i in 1:ncol(data)) {
            if (!(colnames(data)[i] %in% resvars)) {
                if (is.numeric(data[, i]))
                  value <- lapply(list(data[, i]), fn$numeric)[[1]]
                else if (is.ordered(data[, i]))
                  value <- lapply(list(data[, i]), fn$ordered)[[1]]
                else value <- lapply(list(data[, i]), fn$other)[[1]]
                data[, i] <- value
            }
        }
        maxl <- 1
    }
    else {
        maxl <- nrow(data)
    }
    opt <- vars[na.omit(pmatch(names(mc), vars))]
    if (length(opt) > 0)
        for (i in 1:length(opt)) {
            value <- eval(mc[[opt[i]]], envir = env)
            lv <- length(value)
            if (lv > 1)
                if (maxl == 1 || maxl == lv) {
                  maxl <- lv
                  data <- data[1:lv, , drop = FALSE]
                }
                else stop("vector inputs should have the same length.")
            if (is.factor(data[, opt[i]]))
                data[, opt[i]] <- list(factor(value, levels = levels(data[,
                  opt[i]])))
            else if (is.numeric(data[, opt[i]]))
                data[, opt[i]] <- list(as.numeric(value))
            else if (is.logical(data[, opt[i]]))
                data[, opt[i]] <- list(as.logical(value))
            else data[, opt[i]] <- list(value)
        }
    data <- data[1:maxl, , drop = FALSE]
    if (cond) {
        X <- model.frame(tt, data = dta)
        if (!is.null(counter)) {
            X <- list(treat = X[treat == 1, , drop = FALSE],
                control = X[treat == 0, , drop = FALSE])
            class(X$treat) <- class(X$control) <- c("data.frame",
                "cond")
            class(X) <- "setx.counter"
        }
        else class(X) <- c("data.frame", "cond")
    }
    else {
#         X <- as.data.frame(model.matrix(tt, data = data))
	X <- data
    }
    return(X)
}

sim.booltest <-
function(object, x, x1=NULL, num=c(1000, 100),
         prev = NULL, bootstrap = FALSE, bootfn=NULL, cond.data = NULL, ...) {
  if (class(object)[2] != "ML") num <- nrow(object$coefficients)
  else num <- num[1]

  simpar <- param.booltest(object, num, bootstrap)
  simqi <- qi.booltest(object, simpar = simpar, x = x, x1 = x1, y = NULL)
  c <- match.call()
  c$num <- num
  res <- list(x=x, x1=x1, call = c, zelig.call = object$call,
              par = simpar, qi=simqi$qi, qi.name=simqi$qi.name)
  class(res) <- "zelig"
  res
}

plot.zelig.boolean <-
function(x, xlab = "", user.par = FALSE, ...){
  Zelig:::plot.zelig.logit(x, xlab, user.par, ...)
}
