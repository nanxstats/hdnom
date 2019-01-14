#' Print High-Dimensional Cox Model Objects
#'
#' Print information about high-dimensional Cox models.
#'
#' @param x Model object fitted by \code{hdcox.*()} functions.
#' @param ... Other parameters (not used).
#'
#' @method print hdcox.model
#'
#' @export
#'
#' @examples
#' library("survival")
#'
#' # Load imputed SMART data
#' data("smart")
#' x <- as.matrix(smart[, -c(1, 2)])
#' time <- smart$TEVENT
#' event <- smart$EVENT
#' y <- Surv(time, event)
#'
#' fit <- hdcox.lasso(x, y, nfolds = 5, rule = "lambda.1se", seed = 11)
#' print(fit)
print.hdcox.model <- function(x, ...) {
  if (!("hdcox.model" %in% class(x))) {
    stop('x must be of class "hdcox.model" fitted by hdcox.* functions')
  }

  model.type <- gsub("hdcox.model.", "", setdiff(class(x), "hdcox.model"))

  switch(

    model.type,

    lasso = {
      cat("High-Dimensional Cox Model Object\n")
      cat("Random seed:", x$"seed", "\n")
      cat("Model type: lasso\n")
      cat("Best lambda:", x$"lasso_best_lambda", "\n")
    },

    alasso = {
      cat("High-Dimensional Cox Model Object\n")
      cat("Random seed:", x$"seed", "\n")
      cat("Model type: adaptive lasso\n")
      cat("First step best lambda:", x$"ridge_best_lambda", "\n")
      cat("Second step best lambda:", x$"alasso_best_lambda", "\n")
    },

    enet = {
      cat("High-Dimensional Cox Model Object\n")
      cat("Random seed:", x$"seed", "\n")
      cat("Model type: elastic-net\n")
      cat("Best alpha:", x$"enet_best_alpha", "\n")
      cat("Best lambda:", x$"enet_best_lambda", "\n")
    },

    aenet = {
      cat("High-Dimensional Cox Model Object\n")
      cat("Random seed:", x$"seed", "\n")
      cat("Model type: adaptive elastic-net\n")
      cat("First step best alpha:", x$"enet_best_alpha", "\n")
      cat("First step best lambda:", x$"enet_best_lambda", "\n")
      cat("Second step best alpha:", x$"aenet_best_alpha", "\n")
      cat("Second step best lambda:", x$"aenet_best_lambda", "\n")
    },

    mcp = {
      cat("High-Dimensional Cox Model Object\n")
      cat("Random seed:", x$"seed", "\n")
      cat("Model type: MCP\n")
      cat("Best gamma:", x$"mcp_best_gamma", "\n")
      cat("Best lambda:", x$"mcp_best_lambda", "\n")
    },

    mnet = {
      cat("High-Dimensional Cox Model Object\n")
      cat("Random seed:", x$"seed", "\n")
      cat("Model type: Mnet\n")
      cat("Best gamma:", x$"mnet_best_gamma", "\n")
      cat("Best alpha:", x$"mnet_best_alpha", "\n")
      cat("Best lambda:", x$"mnet_best_lambda", "\n")
    },

    scad = {
      cat("High-Dimensional Cox Model Object\n")
      cat("Random seed:", x$"seed", "\n")
      cat("Model type: SCAD\n")
      cat("Best gamma:", x$"scad_best_gamma", "\n")
      cat("Best lambda:", x$"scad_best_lambda", "\n")
    },

    snet = {
      cat("High-Dimensional Cox Model Object\n")
      cat("Random seed:", x$"seed", "\n")
      cat("Model type: Snet\n")
      cat("Best gamma:", x$"snet_best_gamma", "\n")
      cat("Best alpha:", x$"snet_best_alpha", "\n")
      cat("Best lambda:", x$"snet_best_lambda", "\n")
    },

    flasso = {
      cat("High-Dimensional Cox Model Object\n")
      cat("Random seed:", x$"seed", "\n")
      cat("Model type: fused lasso\n")
      cat("Best lambda1:", x$"flasso_best_lambda1", "\n")
      cat("Best lambda2:", x$"flasso_best_lambda2", "\n")
    }
  )
}

#' Make Predictions from High-Dimensional Cox Models
#'
#' Predict overall survival probability at certain time points
#' from established Cox models.
#'
#' @param object Model object fitted by \code{hdcox.*()} functions.
#' @param x Data matrix used to fit the model.
#' @param y Response matrix made with \code{\link[survival]{Surv}}.
#' @param newx Matrix (with named columns) of new values for \code{x}
#' at which predictions are to be made.
#' @param pred.at Time point at which prediction should take place.
#' @param ... Other parameters (not used).
#'
#' @return A \code{nrow(newx) x length(pred.at)} matrix containing
#' overall survival probablity.
#'
#' @method predict hdcox.model
#'
#' @importFrom penalized survival
#' @importMethodsFrom penalized predict
#'
#' @export
#'
#' @examples
#' library("survival")
#'
#' # Load imputed SMART data
#' data("smart")
#' x <- as.matrix(smart[, -c(1, 2)])
#' time <- smart$TEVENT
#' event <- smart$EVENT
#' y <- Surv(time, event)
#'
#' fit <- hdcox.lasso(x, y, nfolds = 5, rule = "lambda.1se", seed = 11)
#' predict(fit, x, y, newx = x[101:105, ], pred.at = 1:10 * 365)
predict.hdcox.model <- function(object, x, y, newx, pred.at, ...) {
  if (!("hdcox.model" %in% class(object))) {
    stop('object must be of class "hdcox.model" fitted by hdcox.* functions')
  }

  model.type <- gsub("hdcox.model.", "", setdiff(class(object), "hdcox.model"))

  if (!("matrix" %in% class(newx))) stop("newx must be a matrix")

  time <- y[, 1L]
  event <- y[, 2L]

  obj.type <- switch(
    model.type,
    lasso = "glmnet", alasso = "glmnet",
    enet = "glmnet", aenet = "glmnet",
    mcp = "ncvreg", mnet = "ncvreg",
    scad = "ncvreg", snet = "ncvreg",
    flasso = "penalized"
  )

  switch(

    obj.type,

    glmnet = {
      lp <- predict(object[[paste0(model.type, "_model")]], x, type = "link")
      basesurv <- glmnet.basesurv(time, event, lp, pred.at)
      lpnew <- predict(object[[paste0(model.type, "_model")]], newx, type = "link")
      p <- exp(exp(lpnew) %*% -t(basesurv$"cumulative_base_hazard"))
    },

    ncvreg = {
      lp <- predict(object[[paste0(model.type, "_model")]], x, type = "link")
      basesurv <- ncvreg.basesurv(time, event, lp, pred.at)
      lpnew <- predict(object[[paste0(model.type, "_model")]], newx, type = "link")
      p <- exp(exp(lpnew) %*% -t(basesurv$"cumulative_base_hazard"))
      # # alternative method using ncvreg built-in prediction directly
      # # almost identical results, but sometimes produces NAs in practice
      # # e.g. pred.at = 1:10 * 365
      # survfun = predict(object[[paste0(model.type, '_model')]], newx, type = 'survival')
      # p = matrix(NA, nrow = nrow(newx), ncol = length(pred.at))
      # for (i in 1L:nrow(newx)) p[i, ] = sapply(pred.at, survfun[[i]])
    },

    penalized = {
      pred <- predict(object[[paste0(model.type, "_model")]], newx)
      p <- matrix(NA, nrow = nrow(newx), ncol = length(pred.at))
      for (i in 1L:length(pred.at)) p[, i] <- survival(pred, time = pred.at[i])
    }
  )

  colnames(p) <- as.character(pred.at)

  p
}

#' Extract Information of Selected Variables from High-Dimensional Cox Models
#'
#' Extract the names and type of selected variables from established
#' high-dimensional Cox models.
#'
#' @param object Model object fitted by \code{hdcox.*()} functions.
#' @param x Data matrix used to fit the model.
#'
#' @export hdnom.varinfo
#'
#' @return A list containing the index, name, type and range of the
#' selected variables.
#'
#' @examples
#' library("survival")
#'
#' # Load imputed SMART data
#' data("smart")
#' x <- as.matrix(smart[, -c(1, 2)])
#' time <- smart$TEVENT
#' event <- smart$EVENT
#' y <- Surv(time, event)
#'
#' # Fit Cox model with lasso penalty
#' fit <- hdcox.lasso(x, y, nfolds = 5, rule = "lambda.1se", seed = 11)
#' hdnom.varinfo(fit, x)
hdnom.varinfo <- function(object, x) {
  if (!("hdcox.model" %in% class(object))) {
    stop('object must be of class "hdcox.model" fitted by hdcox.* functions')
  }

  model.type <- gsub("hdcox.model.", "", setdiff(class(object), "hdcox.model"))

  obj.type <- switch(
    model.type,
    lasso = "glmnet", alasso = "glmnet",
    enet = "glmnet", aenet = "glmnet",
    mcp = "ncvreg", mnet = "ncvreg",
    scad = "ncvreg", snet = "ncvreg",
    flasso = "penalized"
  )

  switch(

    obj.type,

    glmnet = {
      nonzero_idx <- which(as.logical(abs(object[[paste0(model.type, "_model")]][["beta"]]) > .Machine$double.eps))
      nonzero_var <- rownames(object[[paste0(model.type, "_model")]][["beta"]])[nonzero_idx]
    },

    ncvreg = {
      nonzero_idx <- which(object[[paste0(model.type, "_model")]][["beta"]][-1L, ] > .Machine$double.eps)
      nonzero_var <- names(object[[paste0(model.type, "_model")]][["beta"]][-1L, ])[nonzero_idx]
    },

    penalized = {
      nonzero_idx <- which(object[[paste0(model.type, "_model")]]@"penalized" > .Machine$double.eps)
      nonzero_var <- colnames(x)[nonzero_idx]
    }
  )

  varinfo <- list("index" = NULL, "name" = NULL, "type" = NULL, "domain" = NULL)

  varinfo[["index"]] <- nonzero_idx
  varinfo[["name"]] <- nonzero_var
  nvar <- length(varinfo[["name"]])

  is.wholenumber <- function(
    x, tol = .Machine$double.eps^0.5) abs(x - round(x)) < tol

  # type: logical, categorical or continuous
  for (i in 1L:nvar) {
    var <- x[, varinfo[["name"]][i]]
    if (all(is.wholenumber(var)) & nlevels(as.factor(var)) == 2L) {
      varinfo[["type"]][i] <- "logical"
      varinfo[["domain"]][[i]] <- unique(var)
    } else if (all(is.wholenumber(var)) & nlevels(as.factor(var)) > 2L) {
      varinfo[["type"]][i] <- "categorical"
      varinfo[["domain"]][[i]] <- c(min(var), max(var))
    } else if (any(!is.wholenumber(var))) {
      varinfo[["type"]][i] <- "continuous"
      varinfo[["domain"]][[i]] <- c(min(var), max(var))
    } else {
      stop(paste0("unrecognized variable type: ", varinfo[["name"]][i]))
    }
  }

  class(varinfo) <- "hdnom.varinfo"

  varinfo
}
