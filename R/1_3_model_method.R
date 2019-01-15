#' Print high-dimensional Cox model objects
#'
#' Print high-dimensional Cox model objects
#'
#' @param x Model object.
#' @param ... Other parameters (not used).
#'
#' @method print hdnom.model
#'
#' @export
#'
#' @examples
#' data("smart")
#' x <- as.matrix(smart[, -c(1, 2)])
#' time <- smart$TEVENT
#' event <- smart$EVENT
#' y <- survival::Surv(time, event)
#'
#' fit <- fit_lasso(x, y, nfolds = 5, rule = "lambda.1se", seed = 11)
#' print(fit)
print.hdnom.model <- function(x, ...) {
  model <- x$model
  model_type <- x$type

  switch(

    model_type,

    lasso = {
      cat("High-Dimensional Cox Model Object\n")
      cat("Random seed:", x$"seed", "\n")
      cat("Model type: lasso\n")
      cat("Best lambda:", x$"lambda", "\n")
    },

    alasso = {
      cat("High-Dimensional Cox Model Object\n")
      cat("Random seed:", x$"seed", "\n")
      cat("Model type: adaptive lasso\n")
      cat("First step best lambda:", x$"lambda_init", "\n")
      cat("Second step best lambda:", x$"lambda", "\n")
    },

    enet = {
      cat("High-Dimensional Cox Model Object\n")
      cat("Random seed:", x$"seed", "\n")
      cat("Model type: elastic-net\n")
      cat("Best alpha:", x$"alpha", "\n")
      cat("Best lambda:", x$"lambda", "\n")
    },

    aenet = {
      cat("High-Dimensional Cox Model Object\n")
      cat("Random seed:", x$"seed", "\n")
      cat("Model type: adaptive elastic-net\n")
      cat("First step best alpha:", x$"alpha_init", "\n")
      cat("First step best lambda:", x$"lambda_init", "\n")
      cat("Second step best alpha:", x$"alpha", "\n")
      cat("Second step best lambda:", x$"lambda", "\n")
    },

    mcp = {
      cat("High-Dimensional Cox Model Object\n")
      cat("Random seed:", x$"seed", "\n")
      cat("Model type: MCP\n")
      cat("Best gamma:", x$"gamma", "\n")
      cat("Best lambda:", x$"lambda", "\n")
    },

    mnet = {
      cat("High-Dimensional Cox Model Object\n")
      cat("Random seed:", x$"seed", "\n")
      cat("Model type: Mnet\n")
      cat("Best gamma:", x$"gamma", "\n")
      cat("Best alpha:", x$"alpha", "\n")
      cat("Best lambda:", x$"lambda", "\n")
    },

    scad = {
      cat("High-Dimensional Cox Model Object\n")
      cat("Random seed:", x$"seed", "\n")
      cat("Model type: SCAD\n")
      cat("Best gamma:", x$"gamma", "\n")
      cat("Best lambda:", x$"lambda", "\n")
    },

    snet = {
      cat("High-Dimensional Cox Model Object\n")
      cat("Random seed:", x$"seed", "\n")
      cat("Model type: Snet\n")
      cat("Best gamma:", x$"gamma", "\n")
      cat("Best alpha:", x$"alpha", "\n")
      cat("Best lambda:", x$"lambda", "\n")
    },

    flasso = {
      cat("High-Dimensional Cox Model Object\n")
      cat("Random seed:", x$"seed", "\n")
      cat("Model type: fused lasso\n")
      cat("Best lambda1:", x$"lambda1", "\n")
      cat("Best lambda2:", x$"lambda2", "\n")
    }
  )

  invisible(x)
}

#' Make predictions from high-dimensional Cox models
#'
#' Predict overall survival probability at certain time points
#' from fitted Cox models.
#'
#' @param object Model object.
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
#' @method predict hdnom.model
#'
#' @importFrom penalized survival
#' @importMethodsFrom penalized predict
#'
#' @export
#'
#' @examples
#' data("smart")
#' x <- as.matrix(smart[, -c(1, 2)])
#' time <- smart$TEVENT
#' event <- smart$EVENT
#' y <- survival::Surv(time, event)
#'
#' fit <- fit_lasso(x, y, nfolds = 5, rule = "lambda.1se", seed = 11)
#' predict(fit, x, y, newx = x[101:105, ], pred.at = 1:10 * 365)
predict.hdnom.model <- function(object, x, y, newx, pred.at, ...) {
  model <- object$model
  model_type <- object$type

  if (!("matrix" %in% class(newx))) stop("newx must be a matrix")

  time <- y[, 1L]
  event <- y[, 2L]

  obj_type <- switch(
    model_type,
    lasso = "glmnet", alasso = "glmnet", enet = "glmnet", aenet = "glmnet",
    mcp = "ncvreg", mnet = "ncvreg", scad = "ncvreg", snet = "ncvreg",
    flasso = "penalized"
  )

  switch(

    obj_type,

    glmnet = {
      lp <- predict(model, x, type = "link")
      basesurv <- glmnet_basesurv(time, event, lp, pred.at)
      lpnew <- predict(model, newx, type = "link")
      p <- exp(exp(lpnew) %*% -t(basesurv$"cumulative_base_hazard"))
    },

    ncvreg = {
      lp <- predict(model, x, type = "link")
      basesurv <- ncvreg_basesurv(time, event, lp, pred.at)
      lpnew <- predict(model, newx, type = "link")
      p <- exp(exp(lpnew) %*% -t(basesurv$"cumulative_base_hazard"))
      # # alternative method using ncvreg built-in prediction directly
      # # almost identical results, but sometimes produces NAs in practice
      # # e.g. pred.at = 1:10 * 365
      # survfun = predict(model, newx, type = "survival")
      # p = matrix(NA, nrow = nrow(newx), ncol = length(pred.at))
      # for (i in 1L:nrow(newx)) p[i, ] = sapply(pred.at, survfun[[i]])
    },

    penalized = {
      pred <- predict(model, newx)
      p <- matrix(NA, nrow = nrow(newx), ncol = length(pred.at))
      for (i in 1L:length(pred.at)) p[, i] <- survival(pred, time = pred.at[i])
    }
  )

  colnames(p) <- as.character(pred.at)
  p
}

#' Extract information of selected variables from high-dimensional Cox models
#'
#' Extract the names and type of selected variables from fitted
#' high-dimensional Cox models.
#'
#' @param object Model object.
#' @param x Data matrix used to fit the model.
#'
#' @export infer_variable_type
#'
#' @return A list containing the index, name, type and range of the
#' selected variables.
#'
#' @examples
#' data("smart")
#' x <- as.matrix(smart[, -c(1, 2)])
#' time <- smart$TEVENT
#' event <- smart$EVENT
#' y <- survival::Surv(time, event)
#'
#' fit <- fit_lasso(x, y, nfolds = 5, rule = "lambda.1se", seed = 11)
#' infer_variable_type(fit, x)
infer_variable_type <- function(object, x) {
  model <- object$model
  model_type <- object$type

  obj_type <- switch(
    model_type,
    lasso = "glmnet", alasso = "glmnet", enet = "glmnet", aenet = "glmnet",
    mcp = "ncvreg", mnet = "ncvreg", scad = "ncvreg", snet = "ncvreg",
    flasso = "penalized"
  )

  switch(

    obj_type,

    glmnet = {
      nonzero_idx <- which(as.logical(abs(model$beta) > .Machine$double.eps))
      nonzero_var <- rownames(model$beta)[nonzero_idx]
    },

    ncvreg = {
      nonzero_idx <- which(model$beta[-1L, ] > .Machine$double.eps)
      nonzero_var <- names(model$beta[-1L, ])[nonzero_idx]
    },

    penalized = {
      nonzero_idx <- which(model@"penalized" > .Machine$double.eps)
      nonzero_var <- colnames(x)[nonzero_idx]
    }
  )

  res <- list("index" = NULL, "name" = NULL, "type" = NULL, "domain" = NULL)

  res$index <- nonzero_idx
  res$name <- nonzero_var
  nvar <- length(res$name)

  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) abs(x - round(x)) < tol

  # type: logical, categorical or continuous
  for (i in 1L:nvar) {
    var <- x[, res[["name"]][i]]
    if (all(is.wholenumber(var)) & nlevels(as.factor(var)) == 2L) {
      res[["type"]][i] <- "logical"
      res[["domain"]][[i]] <- unique(var)
    } else if (all(is.wholenumber(var)) & nlevels(as.factor(var)) > 2L) {
      res[["type"]][i] <- "categorical"
      res[["domain"]][[i]] <- c(min(var), max(var))
    } else if (any(!is.wholenumber(var))) {
      res[["type"]][i] <- "continuous"
      res[["domain"]][[i]] <- c(min(var), max(var))
    } else {
      stop(paste0("unrecognized variable type: ", res[["name"]][i]))
    }
  }

  class(res) <- "hdnom.variable.type"
  res
}
