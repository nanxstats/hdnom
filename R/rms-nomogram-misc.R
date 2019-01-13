# datadist ---------------------------------------------------------------------

# The following code is derived from the code at
# https://github.com/harrelfe/rms/blob/master/R/datadist.s
# which is a part of the R package rms released under GPL,
# originally written by Frank E Harrell Jr.

rms_datadist <- function(
  ..., data, q.display, q.effect = c(.25, .75),
  adjto.cat = c("mode", "first"), n.unique = 10) {
  adjto.cat <- match.arg(adjto.cat)
  X <- list(...)

  argnames <- as.character(sys.call())[-1]

  if (inherits(x <- X[[1]], "datadist")) {
    Limits <- x$limits
    Values <- x$values
    X[[1]] <- NULL
    argnames <- argnames[-1]
  }
  else {
    Limits <- list()
    Values <- list()
  }

  if (is.data.frame(X[[1]])) {
    if (length(X) > 1) stop("when the first argument is a data frame, no other variables may be specified")
    X <- X[[1]]
  }

  else
    if (is.recursive(X[[1]]) &&
        length(Terms <- X[[1]]$terms) && length(D <- attr(Terms, "Design"))) {
      n <- D$name[D$assume != "interaction"]
      X <- list()
      if (missing(data)) {
        for (nm in n) X[[nm]] <- eval.parent(nm)
      } else
        if (length(names(data))) {
          j <- match(n, names(data), 0)
          if (any(j == 0)) {
            stop(paste(
              "variable(s)",
              paste(n[j == 0], collapse = " "),
              "in model not found on data=, \nwhich has variables",
              paste(names(data), collapse = " ")
            ))
          }
          for (nm in n) X[[nm]] <- data[[nm]]
        }
      else {
        for (nm in n) X[[nm]] <- get(nm, data)
      }
    }
  else {
    if (length(X) & !length(names(X))) names(X) <- argnames[1:length(X)]

    # NEED TO FIX: R has no database.object
    if (!missing(data)) {
      # This duplicative code is for efficiency for large data frames
      stop("program logic error")
      if (length(X)) {
        # if(is.numeric(data)) X <- c(X,database.object(data))
        # else
        X <- c(X, data)
      }
      else {
        # if(is.numeric(data)) X <- database.object(data)
        # else
        X <- data
      }
    }
  }
  nam <- names(X)
  p <- length(nam)
  if (p == 0) stop("you must specify individual variables or a data frame")

  maxl <- 0
  for (i in 1:p) {
    values <- NULL
    x <- X[[i]]
    if (is.character(x)) x <- as.factor(x)
    lx <- length(x)
    lev <- levels(x)
    ll <- length(lev)
    limits <- rep(NA, 5)
    if (is.matrix(x) | (i > 1 && lx != maxl)) {
      warning(paste(nam[i], "is a matrix or has incorrect length; ignored"))
    } else {
      if (ll && (ll < length(x))) values <- lev # if # levels=length(x) is ID variable
      # first look for ordered variable with numeric levels (scored() var)
      if (is.ordered(x) && all_numeric(lev)) {
        levx <- sort(as.numeric(lev))
        limits <- c(
          levx[1], levx[(ll + 1) / 2], levx[ll], levx[1], levx[ll],
          levx[1], levx[ll]
        )
        values <- levx
      }

      else if (ll) {
        adjto <- if (adjto.cat == "first") {
          lev[1]
        } else {
          tab <- table(x)
          (names(tab)[tab == max(tab)])[1]
        }
        limits <- factor(c(NA, adjto, NA, lev[1], lev[ll], lev[1], lev[ll]), levels = lev)
        # non-ordered categorical
      }
      else { # regular numeric variable
        clx <- setdiff(class(x), c("integer", "numeric"))
        # above prevents rounding of quantiles to integers
        y <- x[!is.na(x)]
        n <- length(y)
        if (n < 2) {
          stop(paste("fewer than 2 non-missing observations for", nam[i]))
        }
        values <- sort(unique(y))
        names(values) <- NULL
        nunique <- length(values)
        if (nunique < 2) {
          warning(paste(nam[i], "is constant"))
          limits <- rep(y[1], 7)
        }
        else {
          r <- range(values)
          limits[6:7] <- r
          if (nunique < 4) {
            q <- r
          } else {
            if (missing(q.display)) {
              q.display <- 10 / max(n, 200)
              q.display <- c(q.display, 1 - q.display)
            }
            q <- quantile(unclass(y), q.display)
          } # chron obj. not work here
          limits[4] <- q[1]
          limits[5] <- q[2]
          # check for very poorly distributed categorical numeric variable
          if (limits[4] == limits[5]) limits[4:5] <- r

          # use low category if binary var, middle if 3-level, median otherwise
          if (nunique < 3) {
            limits[2] <- values[1]
          } else
            if (nunique == 3) {
              limits[2] <- values[2]
            } else {
              limits[2] <- median(unclass(y))
            }

          if (nunique < 4) {
            q <- r
          } else {
            q <- quantile(unclass(y), q.effect)
          }
          limits[1] <- q[1]
          limits[3] <- q[2]
          if (limits[1] == limits[3]) limits[c(1, 3)] <- r
          if (nunique > n.unique) values <- NULL
          class(limits) <- clx
        }
      }
      Limits[[nam[i]]] <- limits
      if (length(values)) Values[[nam[i]]] <- values
      maxl <- max(maxl, lx)
    }
  }

  Limits <- structure(
    Limits,
    class = "data.frame",
    row.names = c(
      "Low:effect", "Adjust to",
      "High:effect", "Low:prediction",
      "High:prediction", "Low", "High"
    )
  )
  # data.frame(Limits) gives error with chron objects

  d <- list(limits = Limits, values = Values)
  class(d) <- "datadist"
  d
}

# check if all elements in a character vector are numeric
all_numeric <- function(x, what = c("test", "vector"), extras = c(".", "NA")) {
  what <- match.arg(what)
  x <- sub("[[:space:]]+$", "", x)
  x <- sub("^[[:space:]]+", "", x)
  xs <- x[x %nin% c("", extras)]
  if (!length(xs)) return(if (what == "test") FALSE else x)
  isnum <- suppressWarnings(!any(is.na(as.numeric(xs))))
  if (what == "test") {
    isnum
  } else if (isnum) {
    as.numeric(x)
  } else {
    x
  }
}

# misc -------------------------------------------------------------------------

# The following code is derived from the code at
# https://github.com/harrelfe/rms/blob/master/R/rmsMisc.s
# which is a part of the R package rms released under GPL,
# originally written by Frank E Harrell Jr.

# Returns a list such that variables with no = after them get the value NA
# for handling `...` arguments to nomogram().

rms_args <- function(.object, envir = parent.frame(2)) {
  if (length(.object) < 2) return(NULL)
  .names <- names(.object)[-1]
  # see if no variables given with = after their names
  if (!length(.names)) .names <- rep("", length(.object) - 1)
  .n <- length(.names)
  .vars <- sapply(.object, as.character)[-1]
  .res <- vector("list", .n)
  for (.i in 1:.n)
  {
    if (.names[.i] == "") {
      .names[.i] <- .vars[.i]
      .res[[.i]] <- NA
    }
    else {
      .res[[.i]] <- eval(.object[[.i + 1]], envir = envir)
    }
  }
  names(.res) <- .names
  .res
}

# Function to retrieve limits and values, from fit (if they are there)
# or from a datadist object. If need.all = FALSE and input is coming from
# datadist, insert columns with NAs for variables not defined.
#
# @param at essentially fit$Design
# @param x predictor matrix (non-zero varaible only) to be converted to datadist

get_lim <- function(at, x, allow.null = FALSE, need.all = TRUE) {
  nam <- at$name[at$assume != "interaction"]
  limits <- at$limits
  values <- at$values

  XDATADIST <- rms_datadist(as.data.frame(x))
  # XDATADIST <- .Options$datadist
  X <- lims <- vals <- NULL
  if (!is.null(XDATADIST)) {
    X <- XDATADIST
    lims <- X$limits
    if (is.null(lims)) stop("datadist conversion error in get_lim()")
    vals <- X$values
  }

  if ((length(X) + length(limits)) == 0) {
    if (allow.null) {
      lims <- list()
      for (nn in nam) lims[[nn]] <- rep(NA, 7)
      lims <- structure(
        lims,
        class = "data.frame",
        row.names = c(
          "Low:effect", "Adjust to", "High:effect", "Low:prediction",
          "High:prediction", "Low", "High"
        )
      )
      return(list(limits = lims, values = values))
    }
    stop("no datadist in effect now or during model fit")
  }

  na <- if (length(limits)) {
    sapply(limits, function(x) all(is.na(x)))
  } else {
    rep(TRUE, length(nam))
  }
  if (length(lims) && any(na)) {
    for (n in nam[na]) { # if() assumes NA stored in fit
      # for missing vars
      z <- limits[[n]]
      u <- if (match(n, names(lims), 0) > 0) lims[[n]] else NULL
      # This requires exact name match, not substring match
      if (is.null(u)) {
        if (need.all) {
          stop(paste(
            "variable", n,
            "does not have limits defined in fit or with datadist"
          ))
        } else {
          limits[[n]] <- rep(NA, 7)
        } # Added 28 Jul 94
      }
      else {
        limits[[n]] <- u
      }
    }
  }
  limits <- structure(
    limits,
    class = "data.frame",
    row.names = c(
      "Low:effect", "Adjust to", "High:effect", "Low:prediction",
      "High:prediction", "Low", "High"
    )
  )

  if (length(vals)) {
    values <- c(
      values,
      vals[match(names(vals), nam, 0) > 0 & match(names(vals), names(values), 0) == 0]
    )
  } # add in values from datadist corresponding to vars in model
  # not already defined for model

  list(limits = limits, values = values)
}

# Function to return limits for an individual variable, given an object
# created by get_lim

get_limi <- function(name, Limval, need.all = TRUE) {
  lim <- if (match(name, names(Limval$limits), 0) > 0) {
    Limval$limits[[name]]
  } else {
    NULL
  }
  if (is.null(Limval) || is.null(lim) || all(is.na(lim))) {
    if (need.all) {
      stop(paste(
        "no limits defined by datadist for variable",
        name
      ))
    }
    return(rep(NA, 7))
  }
  lim
}

# rms_levels
# Make each variable in an input data frame that is a
# factor variable in the model be a factor variable with
# the levels that were used in the model.  This is primarily
# so that row insertion will work right with <-[.data.frame
#
# at=Design attributes

rms_levels <- function(df, at) {
  ac <- at$assume.code
  for (nn in names(df))
  {
    j <- match(nn, at$name, 0)
    if (j > 0) {
      if ((ac[j] == 5 | ac[j] == 8) & length(lev <- at$parms[[nn]])) {
        df[[nn]] <- factor(df[[nn]], lev)
      }
    }
  }
  df
}

# Function to remove one or more terms from a model formula, using
# strictly character manipulation.  This handles problems such as
# [.terms removing offset() if you subset on anything
# For each character string in which, terms like string(...) are removed.

remove_formula_terms <- function(form, which = NULL, delete.response = FALSE) {
  if ("offset" %in% which) {
    form <- formula(terms(form)[TRUE])
    which <- setdiff(which, "offset")
  }
  # [.terms ignores offset variables.  Above logic handles nested () unlike
  # what is below
  form <- paste(deparse(form), collapse = "") # no string splitting
  if (delete.response) form <- gsub(".*~", "~", form)
  for (w in which) {
    pattern <- sprintf("\\+?[ ]*?%s\\(.*?\\)[ ]*?\\+{0,1}", w) # assume additive form
    form <- gsub(pattern, "", form)
  }
  as.formula(form)
}

# Function to return a list whose ith element contains indexes
# of all predictors related, indirectly or directly, to predictor i
# Predictor i and j are related indirectly if they are related to
# any predictors that interact
# Set type="direct" to only include factors interacting with i
# This function is used by nomogram.

related_predictors <- function(at, type = c("all", "direct")) {
  type <- match.arg(type)
  f <- sum(at$assume.code < 9)
  if (any(at$assume.code == 10)) stop("does not work with matrix factors")
  ia <- at$interactions
  x <- rep(NA, f)
  names(x) <- at$name[at$assume.code < 9]
  mode(x) <- "list"
  if (length(ia) == 0) {
    for (i in 1:f) x[[i]] <- integer(0)
    return(x)
  }
  for (i in 1:f)
  {
    r <- integer(0)
    for (j in 1:ncol(ia))
    {
      w <- ia[, j]
      if (any(w == i)) r <- c(r, w[w > 0 & w != i])
    }
    x[[i]] <- r
  }
  if (type == "direct") return(x)

  while (TRUE) {
    bigger <- FALSE
    for (j in 1:f)
    {
      xj <- x[[j]]
      y <- unlist(x[xj])
      y <- y[y != j]
      new <- unique(c(y, xj))
      bigger <- bigger | length(new) > length(xj)
      x[[j]] <- new
    }
    if (!bigger) break
  }
  x
}

# check values -----------------------------------------------------------------

# The following code is derived from the code at
# https://github.com/harrelfe/rms/blob/master/R/rms.trans.s
# which is a part of the R package rms released under GPL,
# originally written by Frank E Harrell Jr.

# If x is NA, returns list of possible values of factor i defined
#   in object f's attributes. For continuous factors, returns n values
# 	in default prediction range. Use n = 0 to return trio of effect
# 	limits. Use n<0 to return pretty(plotting range, nint = -n).
#   If type.range = "full" uses the full range instead of default plot rng.
# If x is not NA, checks that list to see that each value is allowable
# 	for the factor type, and returns x.
# The last argument is object returned from get_lim().
# The first argument is the Design list.

check_values <- function(f, i, x, n, limval, type.range = "plot") {
  as <- f$assume.code[i]
  name <- f$name[i]
  parms <- f$parms[[name]]
  isna <- length(x) == 1 && is.na(x)
  values <- limval$values[[name]]
  charval <- length(values) && is.character(values)
  if (isna & as != 7) {
    if (!length(limval) || match(name, dimnames(limval$limits)[[2]], 0) == 0 ||
        is.na(limval$limits["Adjust to", name])) {
      stop(paste("variable", name, "does not have limits defined by datadist"))
    }

    limits <- limval$limits[, name]
    lim <- if (type.range == "full") limits[6:7] else limits[4:5]
  }

  if (as < 5 | as == 6) {
    if (isna) {
      if (!length(values)) {
        if (n == 0) {
          x <- limits[1:3]
        } else {
          if (n > 0) {
            x <- seq(
              unclass(lim[1]), # handles chron
              unclass(lim[2]),
              length = n
            )
          } else {
            x <- pretty(unclass(lim[1:2]), n = -n)
          }
          class(x) <- class(lim)
        }
      } else {
        x <- values
      }
    } else {
      if (is.character(x) && !charval) {
        stop(paste(
          "character value not allowed for variable",
          name
        ))
      } # Allow any numeric value
      if (charval) {
        j <- match(x, values, 0)
        if (any(j == 0)) {
          stop(paste(
            "illegal values for categorical variable:",
            paste(x[j == 0], collapse = " "), "\nPossible levels:",
            paste(values, collapse = " ")
          ))
        }
      }
    }
  } else if (as == 5 | as == 8) {
    if (isna) {
      x <- parms
    } else {
      j <- match(x, parms, 0) # match converts x to char if needed
      if (any(j == 0)) {
        stop(paste(
          "illegal levels for categorical variable:",
          paste(x[j == 0], collapse = " "), "\nPossible levels:",
          paste(parms, collapse = " ")
        ))
      }
      x
    }
  }
  else if (as == 7) {
    if (isna) {
      x <- parms
    } else if (is.character(x)) {
      stop(paste(
        "character value not allowed for",
        "variable", name
      ))
    } else {
      j <- match(x, parms, 0)
      if (any(j == 0)) {
        stop(paste(
          "illegal levels for categorical variable:",
          paste(x[j == 0], collapse = " "), "\n", "Possible levels:",
          paste(parms, collapse = " ")
        ))
      }
    }
  }

  invisible(x)
}

# predict ----------------------------------------------------------------------

# The following code is derived from the code at
# https://github.com/harrelfe/rms/blob/master/R/predictrms.s
# which is a part of the R package rms released under GPL,
# originally written by Frank E Harrell Jr.

#' @importFrom stats model.extract model.frame model.matrix qnorm qt terms

rms_predict <- function(
  fit, newdata = NULL,
  type = c(
    "lp", "x", "data.frame", "terms", "cterms", "ccterms",
    "adjto", "adjto.data.frame", "model.frame"
  ),
  se.fit = FALSE, conf.int = FALSE,
  conf.type = c("mean", "individual", "simultaneous"),
  kint = NULL,
  na.action = na.keep, expand.na = TRUE,
  center.terms = type == "terms", ref.zero = FALSE, ...) {
  type <- match.arg(type)
  conf.type <- match.arg(conf.type)
  if (conf.type == "simultaneous") {
    # this requires library("multcomp") so not doing here...
    if (missing(newdata) || !length(newdata)) {
      stop('newdata must be given if conf.type = "simultaneous"')
    }
  }

  at <- fit$Design
  assume <- at$assume.code
  Limval <- get_lim(at, x = fit$x, allow.null = TRUE, need.all = FALSE)
  Values <- Limval$values
  non.ia <- assume != 9L
  non.strat <- assume != 8L
  f <- sum(non.ia)
  nstrata <- sum(assume == 8L)
  somex <- any(non.strat)
  rnam <- NULL
  cox <- inherits(fit, "cph")
  naa <- fit$na.action
  if (!expand.na) {
    naresid <- function(a, b) b
  } # don't really call naresid if drop NAs

  parms <- at$parms
  name <- at$name
  coeff <- fit$coefficients
  nrp <- num_intercepts(fit)
  nrpcoef <- num_intercepts(fit, "coef")
  if (!length(kint)) kint <- fit$interceptRef # orm or lrm otherwise NULL

  int.pres <- nrp > 0L

  assign <- fit$assign
  nama <- names(assign)[1L]
  asso <- 1 * (nama == "Intercept" | nama == "(Intercept)")

  Center <- if (cox) fit$center else 0.

  oldopts <- options(
    contrasts = c(
      factor = "contr.treatment",
      ordered = "contr.treatment"
    ), # was "contr.poly"
    Design.attr = at
  )

  # In SV4 options(two lists) causes problems
  on.exit({
    options(contrasts = oldopts$contrasts)
    options(Design.attr = NULL)
  })

  # Formula without resposne variable any offsets:
  formulano <- remove_formula_terms(
    fit$sformula,
    which = "offset",
    delete.response = TRUE
  )

  offset <- 0
  offpres <- FALSE
  # # offset is ignored for prediction (offset set to zero)
  # if(! missing(newdata) && length(newdata)) {
  #    offset <- model.offset(model.frame(remove_formula_terms(
  #      fit$sformula, delete.response = TRUE), newdata, na.action = na.action, ...
  #    ))
  # offpres <- length(offset) > 0
  # if(! offpres) offset <- 0
  # }

  Terms <- terms(formulano, specials = "strat")

  attr(Terms, "response") <- 0L
  attr(Terms, "intercept") <- 1L
  # Need intercept whenever design matrix is generated to get
  # current list of dummy variables for factor variables
  stra <- attr(Terms, "specials")$strat

  Terms.ns <- if (length(stra)) Terms[-stra] else Terms

  if (conf.int) {
    vconstant <- 0.
    if (conf.type == "individual") {
      vconstant <- fit$stats["Sigma"]^2
      if (is.na(vconstant)) {
        stop('conf.type="individual" requires that fit be from ols')
      }
    }
    zcrit <- if (length(idf <- fit$df.residual)) {
      qt((1. + conf.int) / 2., idf)
    } else {
      qnorm((1. + conf.int) / 2.)
    }
  }

  # Form design matrix for adjust-to values
  # Result of adj_to() is a model matrix with no intercept(s)
  adj_to <- function(type) {
    adjto <- list()
    ii <- 0L
    for (i in (1L:length(assume))[non.ia]) {
      ii <- ii + 1L
      xi <- get_limi(name[i], Limval, need.all = TRUE)[2L]
      if (assume[i] %in% c(5L, 8L)) {
        xi <- factor(xi, parms[[name[i]]])
      } else
        if (assume[i] == 7L) {
          stop("scored() is not implemented in hdnom")
          # xi <- scored(xi, name = name[i])
        } else
          if (assume[i] == 10L) {
            xi <- matrix(parms[[name[i]]], nrow = 1)
          } # matrx col medians
      adjto[[ii]] <- xi
    }
    names(adjto) <- name[non.ia]
    attr(adjto, "row.names") <- "1"
    class(adjto) <- "data.frame"
    if (type == "adjto.data.frame") return(adjto)
    adjto <- model.frame(Terms, adjto)
    adjto <- model.matrix(Terms.ns, adjto)[, -1, drop = FALSE]
    if (type == "adjto") {
      k <- (nrpcoef + 1L):length(coeff)
      nck <- names(coeff)[k]
      if (is.matrix(adjto)) {
        dimnames(adjto) <- list(dimnames(adjto)[[1L]], nck)
      } else {
        names(adjto) <- nck
      }
    }
    adjto
  }

  adjto <- NULL

  if (type %nin% c("adjto", "adjto.data.frame")) {
    X <- NULL
    if (missing(newdata) || !length(newdata)) {
      flp <- fit$linear.predictors
      if (type == "lp" && length(flp)) {
        LP <- naresid(naa, flp)
        if (int.pres) {
          lpkint <- attr(flp, "intercepts")
          if (!length(lpkint)) lpkint <- 1L
          if (length(kint) && kint != lpkint) {
            LP <- LP - coeff[lpkint] + coeff[kint]
          }
        }
        if (length(stra <- fit$strata)) {
          attr(LP, "strata") <- naresid(naa, stra)
        }
        if (!se.fit && !conf.int) {
          return(LP)
        } else
          if (length(fit$se.fit)) {
            if (nrp > 1L) {
              warning("se.fit is retrieved from the fit but it corresponded to kint")
            }
            retlist <- list(linear.predictors = LP)
            if (se.fit) retlist$se.fit <- naresid(naa, fit$se.fit)
            if (conf.int) {
              plminus <- zcrit * sqrt(retlist$se.fit^2 + vconstant)
              retlist$lower <- LP - plminus
              retlist$upper <- LP + plminus
            }
            return(retlist)
          }
      } # end type='lp'
      else
        if (type == "x") {
          return(
            structure(
              naresid(naa, fit$x),
              strata = if (length(stra <- fit$strata)) {
                naresid(naa, stra)
              } else {
                NULL
              }
            )
          )
        }
      X <- fit[["x"]]
      rnam <- dimnames(X)[[1]]
      if (!length(X)) {
        stop("newdata not given and fit did not store x")
      }
    } # end no newdata
    if (!length(X)) {
      if (!is.data.frame(newdata)) {
        if (is.list(newdata)) {
          loc <- if (!length(names(newdata))) 1L:f else name[assume != 9L]
          new <- matrix(double(1L), nrow = length(newdata[[1L]]), ncol = length(newdata))
          for (j in 1L:ncol(new)) new[, j] <- newdata[[loc[j]]]
          newdata <- new
        }
        if (!is.matrix(newdata)) newdata <- matrix(newdata, ncol = f)
        if (ncol(newdata) != f) {
          stop("# columns in newdata not= # factors in design")
        }
        X <- list()
        k <- 0L
        ii <- 0L
        for (i in (1L:length(assume))[non.ia]) {
          ii <- ii + 1L
          xi <- newdata[, ii]
          as <- assume[i]
          allna <- all(is.na(xi))
          if (as == 5L | as == 8L) {
            xi <- as.integer(xi)
            levels(xi) <- parms[[name[i]]]
            class(xi) <- "factor"
          }
          else if (as == 7L) {
            stop("scored() is not implemented in hdnom")
            # xi <- scored(xi, name = name[i])
          } else if (as == 10L) {
            if (i == 1) {
              ifact <- 1L
            } else {
              ifact <- 1L + sum(assume[1L:(i - 1L)] != 8L)
            }
            # Accounts for assign not being output for strata factors
            ncols <- length(assign[[ifact + asso]])
            if (allna) {
              xi <- matrix(double(1L), nrow = length(xi), ncol = ncols)
              for (j in 1L:ncol(xi)) xi[, j] <- parms[[name[i]]][j]
            }
            else {
              xi <- matrix(xi, nrow = length(xi), ncol = ncols)
            }
          }
          # Duplicate single value for all parts of matrix
          k <- k + 1L
          X[[k]] <- xi
        }
        names(X) <- name[non.ia]
        attr(X, "row.names") <- as.character(1L:nrow(newdata))
        class(X) <- "data.frame"
        newdata <- X
        # Note: data.frame() converts matrix variables to individual variables
        if (type == "data.frame") return(newdata)
      } # end !is.data.frame(newdata)
      else {
        # Need to convert any factors to have all levels in original fit
        # Otherwise, wrong dummy variables will be generated by model.matrix

        nm <- names(newdata)
        for (i in 1L:ncol(newdata)) {
          j <- match(nm[i], name)
          if (!is.na(j)) {
            asj <- assume[j]
            w <- newdata[, i]
            V <- NULL
            if (asj %in% c(5L, 7L, 8L) |
                (name[j] %in% names(Values) &&
                 length(V <- Values[[name[j]]]) && is.character(V))) {
              if (length(Pa <- parms[[name[j]]])) V <- Pa
              newdata[, i] <- factor(w, V)
              # Handles user specifying numeric values without quotes, that are levels
              ww <- is.na(newdata[, i]) & !is.na(unclass(w))
              if (any(ww)) {
                cat(
                  "Error in predictrms: Values in", names(newdata)[i],
                  "not in", V, ":\n"
                )
                print(as.character(w[ww]), quote = FALSE)
                stop()
              }
            }
          }
        }
      } # is.data.frame(newdata)
      X <- model.frame(Terms, newdata, na.action = na.action, ...)
      if (type == "model.frame") return(X)
      naa <- attr(X, "na.action")
      rnam <- row.names(X)

      strata <- list()
      nst <- 0
      ii <- 0
      for (i in 1L:ncol(X)) {
        ii <- ii + 1L
        xi <- X[[i]]
        asi <- attr(xi, "assume.code")
        as <- assume[ii]
        if (!length(asi) && as == 7L) {
          stop("scored() is not implemented in hdnom")
          # attr(X[, i], "contrasts") <- attr(scored(xi, name = name[ii]), "contrasts")
          if (length(xi) == 1L) warning("a bug in model.matrix can produce incorrect results\nwhen only one observation is being predicted for an ordered variable")
        }

        if (as == 8L) {
          nst <- nst + 1L
          ste <- paste(name[ii], parms[[name[ii]]], sep = "=")
          strata[[nst]] <- factor(ste[X[, i]], ste)
        }
      }
      X <- if (!somex) {
        NULL
      } else {
        model.matrix(Terms.ns, X)[, -1L, drop = FALSE]
      }

      if (nstrata > 0L) {
        names(strata) <- paste("S", 1L:nstrata, sep = "")
        strata <- interaction(as.data.frame(strata), drop = FALSE)
      }
    } # end !length(X)
    else {
      strata <- attr(X, "strata")
    }
  } # if(end adj.to adj.to.data.frame)
  if (somex) {
    cov <- matrix(1e-33, ncol = ncol(fit$x) + 1, nrow = ncol(fit$x) + 1)
    colnames(cov) <- c("Intercept", colnames(fit$x))
    rownames(cov) <- c("Intercept", colnames(fit$x))
    # replaces the original
    # cov <- vcov(fit, regcoef.only = TRUE, intercepts = kint)
    # since
    # vcov.lm <- function(object, ...) { so <- summary.lm(object); so$sigma^2 * so$cov.unscaled }
    # and sigma -> 0
    covnoint <- if (nrp == 0) {
      cov
    } else {
      cov <- matrix(1e-33, ncol = ncol(fit$x), nrow = ncol(fit$x))
      colnames(cov) <- colnames(fit$x)
      rownames(cov) <- colnames(fit$x)
      # vcov(fit, regcoef.only = TRUE, intercepts = "none")
    }
  }

  if (type %in% c("adjto.data.frame", "adjto")) return(adj_to(type))

  if (type == "x") {
    return(
      structure(
        naresid(naa, X),
        strata = if (nstrata > 0) naresid(naa, strata) else NULL,
        na.action = if (expand.na) NULL else naa
      )
    )
  }


  if (type == "lp") {
    if (somex) {
      xb <- matxv(X, coeff, kint = kint) - Center + offset
      names(xb) <- rnam
    }
    else {
      xb <- if (offpres) offset else numeric(0)
      if (nstrata > 0) attr(xb, "strata") <- naresid(naa, strata)
      return(structure(if (se.fit) {
        list(
          linear.predictors = xb,
          se.fit = rep(NA, length(xb))
        )
      } else {
        xb
      },
      na.action = if (expand.na) NULL else naa
      ))
    }
    xb <- naresid(naa, xb)
    if (nstrata > 0) attr(xb, "strata") <- naresid(naa, strata)
    ycenter <- if (ref.zero && somex) {
      if (!length(adjto)) adjto <- adj_to(type)
      matxv(adjto, coeff, kint = kint) - Center
    } else {
      0.
    }

    if (ref.zero || ((se.fit || conf.int) && somex)) {
      dx <- dim(X)
      n <- dx[1L]
      p <- dx[2L]
      if (cox && !ref.zero) X <- X - rep(fit$means, rep.int(n, p))
      if (ref.zero) {
        if (!length(adjto)) adjto <- adj_to(type)
        X <- X - rep(adjto, rep.int(n, p))
      }
      se <- drop(if (ref.zero || nrp == 0L) {
        sqrt(((X %*% covnoint) * X) %*% rep(1L, ncol(X)))
      } else {
        Xx <- cbind(Intercept = 1., X)
        sqrt(((Xx %*% cov) * Xx) %*% rep(1L, ncol(Xx)))
      })
      names(se) <- rnam

      sef <- naresid(naa, se)
      ww <- if (conf.int || se.fit) {
        if (se.fit) {
          list(linear.predictors = xb - ycenter, se.fit = sef)
        } else {
          list(linear.predictors = xb - ycenter)
        }
      }
      else {
        xb - ycenter
      }
      retlist <- structure(ww, na.action = if (expand.na) NULL else naa)
      if (conf.int) {
        if (conf.type == "simultaneous") {
          stop("conf.type = \"simultaneous\" is not implemented in hdnom")
          # num.intercepts.not.in.X <- length(coeff) - ncol(X)
          # u <- confint(
          #   multcomp::glht(
          #     fit,
          #     if (num.intercepts.not.in.X == 0L) X else Xx,
          #     df = if (length(idf)) idf else 0L
          #   ),
          #   level = conf.int
          # )$confint
          # retlist$lower <- u[, "lwr"]
          # retlist$upper <- u[, "upr"]
        } else {
          plminus <- zcrit * sqrt(sef^2 + vconstant)
          retlist$lower <- xb - plminus - ycenter
          retlist$upper <- xb + plminus - ycenter
        }
      }
      return(retlist)
    }
    else {
      return(structure(xb - ycenter, na.action = if (expand.na) NULL else naa))
    }
  } # end if type = "lp"

  if (type %in% c("terms", "cterms", "ccterms")) {
    if (!somex) {
      stop('type="terms" may not be given unless covariables present')
    }

    usevar <- if (type == "terms") non.strat else rep(TRUE, length(assume))
    fitted <- array(
      0, c(nrow(X), sum(usevar)),
      list(rnam, name[usevar])
    )
    if (se.fit) se <- fitted
    if (center.terms) {
      if (!length(adjto)) adjto <- adj_to(type)
      if (ncol(adjto) != ncol(X)) {
        if (dimnames(adjto)[[2L]][1L] %in% c("Intercept", "(Intercept)") &&
            dimnames(X)[[2L]][1L] %nin% c("Intercept", "(Intercept)")) {
          adjto <- adjto[, -1L, drop = FALSE]
        }
        if (ncol(adjto) != ncol(X)) stop("program logic error")
      }
      X <- sweep(X, 2L, adjto) # center columns
    }
    j <- 0L
    for (i in (1L:length(assume))[usevar]) {
      j <- j + 1L
      if (assume[i] != 8L) { # non-strat factor; otherwise leave fitted=0
        k <- assign[[j + asso]]
        num.intercepts.not.in.X <- length(coeff) - ncol(X)
        ko <- k - num.intercepts.not.in.X
        fitted[, j] <- matxv(X[, ko, drop = FALSE], coeff[k])
        if (se.fit) {
          se[, j] <-
            (((X[, ko, drop = FALSE] %*% cov[k, k, drop = FALSE]) *
                X[, ko, drop = FALSE]) %*% rep(1., length(ko)))^.5
        }
      }
    }
    if (type == "cterms") {
      # Combine all related interation terms with main effect terms
      w <- fitted[, non.ia, drop = FALSE] # non-interaction terms
      for (i in 1L:f) {
        ia <- interactions.containing(at, i)
        # subscripts of interaction terms related to predictor i
        if (length(ia)) w[, i] <- rowSums(fitted[, c(i, ia), drop = FALSE])
      }
      fitted <- w
    }

    if (type == "ccterms") {
      stop("type = \"ccterms\" is not implemented in hdnom")
      # z <- combineRelatedPredictors(at)
      # f <- length(z$names)
      # w <- matrix(NA, ncol = f, nrow = nrow(fitted))
      # colnames(w) <- sapply(z$names, paste, collapse = ", ")
      # for (i in 1L:f)
      #   w[, i] <- rowSums(fitted[, z$namesia[[i]], drop = FALSE])
      # fitted <- w
    }

    fitted <- structure(naresid(naa, fitted), strata = if (nstrata == 0) NULL else naresid(naa, strata))

    if (se.fit) {
      return(structure(list(fitted = fitted, se.fit = naresid(naa, se)), na.action = if (expand.na) NULL else naa))
    }
    else {
      return(structure(fitted, na.action = if (expand.na) NULL else naa))
    }
  }
}

"%nin%" <- function(x, table) match(x, table, nomatch = 0) == 0

na.keep <- function(mf) {
  w <- na.detail.response(mf)
  if (length(w)) {
    class(w) <- "keep"
  }

  attr(mf, "na.action") <- w
  mf
}

na.detail.response <- function(mf) {
  if (is.null(z <- .Options$na.detail.response) || !z) {
    return(NULL)
  }

  response <- model.extract(mf, response)
  if (is.null(response)) {
    return(NULL)
  }

  if (!is.matrix(response)) {
    response <- as.matrix(response)
  }

  GFUN <- options()$na.fun.response
  if (is.null(GFUN)) {
    GFUN <- function(x, ...) {
      if (is.matrix(x)) x <- x[, ncol(x)]
      x <- x[!is.na(x)]
      c(N = length(x), Mean = mean(x))
    }
  } else {
    GFUN <- eval.parent(as.name(GFUN))
  }

  w <- NULL
  nam <- names(mf)
  wnam <- NULL
  N <- nrow(mf)
  p <- ncol(mf)
  omit <- rep(FALSE, N)
  for (i in 2:p) {
    x <- mf[, i]
    if (is.matrix(x)) {
      x <- x[, 1]
    }

    isna <- is.na(x)
    omit <- omit | isna
    nmiss <- sum(isna)
    if (nmiss) {
      w <- cbind(w, GFUN(response[isna, ]))
      wnam <- c(wnam, paste(nam[i], "=NA", sep = ""))
    }

    n <- N - nmiss
    if (n) {
      w <- cbind(w, GFUN(response[!isna, ]))
      wnam <- c(wnam, paste(nam[i], "!=NA", sep = ""))
    }
  }

  # summarize responce for ANY x missing
  if (p > 2) {
    nmiss <- sum(omit)
    if (nmiss) {
      w <- cbind(w, GFUN(response[omit, ]))
      wnam <- c(wnam, "Any NA")
    }

    if (N - nmiss) {
      w <- cbind(w, GFUN(response[!omit, ]))
      wnam <- c(wnam, "No NA")
    }
  }

  dimnames(w)[[2]] <- wnam
  w
}


# multiply matrix by a vector

matxv <- function(a, b, kint = 1, bmat = FALSE) {
  bi <- attr(b, "intercepts")
  lbi <- length(bi)
  lkint <- length(kint)
  if (lkint > 1L) stop("kint must have length 0 or 1")

  if (bmat) {
    if (!is.matrix(a)) stop("a must be a matrix when b is a matrix")
    ca <- ncol(a)
    cb <- ncol(b)
    if (cb < ca) stop("number of columns in b must be >= number in a")
    if (cb == ca) return(a %*% t(b))
    excess <- cb - ca
    xx <- matrix(0, nrow = nrow(a), ncol = excess)
    if (lbi && lkint) {
      if (lbi != excess) {
        stop("b intercepts attribute has different length from number of excess elements in b")
      }
      bi <- round(bi)
      kint <- round(kint)
      if (!isTRUE(all.equal(sort(bi), sort(kint)))) {
        stop("b intercepts attribute do not match kint")
      }
      xx[] <- 1.
    }
    else if (lkint) {
      if (kint > excess) {
        stop("kint > number of excess elements in b")
      }
      xx[, kint] <- 1.
    }
    return(cbind(xx, a) %*% t(b))
  }

  if (!is.matrix(a)) {
    a <- if (length(b) == 1L) matrix(a, ncol = 1L) else matrix(a, nrow = 1L)
  }

  nc <- dim(a)[2]
  lb <- length(b)
  if (lb < nc) {
    stop(paste("columns in a (", nc, ") must be <= length of b (", length(b), ")", sep = ""))
  }

  if (nc == lb) return(drop(a %*% b))

  excess <- lb - nc
  if (lbi && lkint) {
    if (lbi != excess) {
      stop("b intercepts attribute has different length from number of excess elements in b")
    }
    bi <- round(bi)
    kint <- round(kint)
    if (!isTRUE(all.equal(sort(bi), sort(kint)))) {
      stop("b intercepts attribute do not match kint")
    }
    bkint <- b[1]
  }
  else if (lkint) {
    if (kint > excess) {
      stop("kint > number excess elements in b")
    }

    bkint <- b[kint]
  }
  else {
    bkint <- 0.
  }
  drop(bkint + (a %*% b[(lb - nc + 1L):lb]))
}
