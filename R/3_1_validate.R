#' Validate high-dimensional Cox models with time-dependent AUC
#'
#' Validate high-dimensional Cox models with time-dependent AUC
#'
#' @param x Matrix of training data used for fitting the model;
#' on which to run the validation.
#' @param time Survival time.
#' Must be of the same length with the number of rows as \code{x}.
#' @param event Status indicator, normally 0 = alive, 1 = dead.
#' Must be of the same length with the number of rows as \code{x}.
#' @param model.type Model type to validate. Could be one of \code{"lasso"},
#' \code{"alasso"}, \code{"flasso"}, \code{"enet"}, \code{"aenet"},
#' \code{"mcp"}, \code{"mnet"}, \code{"scad"}, or \code{"snet"}.
#' @param alpha Value of the elastic-net mixing parameter alpha for
#' \code{enet}, \code{aenet}, \code{mnet}, and \code{snet} models.
#' For \code{lasso}, \code{alasso}, \code{mcp}, and \code{scad} models,
#' please set \code{alpha = 1}.
#' \code{alpha=1}: lasso (l1) penalty; \code{alpha=0}: ridge (l2) penalty.
#' Note that for \code{mnet} and \code{snet} models,
#' \code{alpha} can be set to very close to 0 but not 0 exactly.
#' @param lambda Value of the penalty parameter lambda to use in the
#' model fits on the resampled data. From the fitted Cox model.
#' @param pen.factor Penalty factors to apply to each coefficient.
#' From the fitted \emph{adaptive lasso} or \emph{adaptive elastic-net} model.
#' @param gamma Value of the model parameter gamma for
#' MCP/SCAD/Mnet/Snet models.
#' @param lambda1 Value of the penalty parameter lambda1 for fused lasso model.
#' @param lambda2 Value of the penalty parameter lambda2 for fused lasso model.
#' @param method Validation method.
#' Could be \code{"bootstrap"}, \code{"cv"}, or \code{"repeated.cv"}.
#' @param boot.times Number of repetitions for bootstrap.
#' @param nfolds Number of folds for cross-validation and
#' repeated cross-validation.
#' @param rep.times Number of repeated times for repeated cross-validation.
#' @param tauc.type Type of time-dependent AUC.
#' Including \code{"CD"} proposed by Chambless and Diao (2006).,
#' \code{"SZ"} proposed by Song and Zhou (2008).,
#' \code{"UNO"} proposed by Uno et al. (2007).
#' @param tauc.time Numeric vector. Time points at which to evaluate
#' the time-dependent AUC.
#' @param seed A random seed for resampling.
#' @param trace Logical. Output the validation progress or not.
#' Default is \code{TRUE}.
#'
#' @export validate
#'
#' @references
#' Chambless, L. E. and G. Diao (2006).
#' Estimation of time-dependent area under the ROC curve for long-term
#' risk prediction.
#' \emph{Statistics in Medicine} 25, 3474--3486.
#'
#' Song, X. and X.-H. Zhou (2008).
#' A semiparametric approach for the covariate specific ROC curve with
#' survival outcome.
#' \emph{Statistica Sinica} 18, 947--965.
#'
#' Uno, H., T. Cai, L. Tian, and L. J. Wei (2007).
#' Evaluating prediction rules for t-year survivors with censored
#' regression models.
#' \emph{Journal of the American Statistical Association} 102, 527--537.
#'
#' @examples
#' library("survival")
#'
#' # Load imputed SMART data
#' data(smart)
#' x <- as.matrix(smart[, -c(1, 2)])[1:500, ]
#' time <- smart$TEVENT[1:500]
#' event <- smart$EVENT[1:500]
#' y <- Surv(time, event)
#'
#' # Fit penalized Cox model with lasso penalty
#' fit <- fit_lasso(x, y, nfolds = 5, rule = "lambda.1se", seed = 11)
#'
#' # Model validation by bootstrap with time-dependent AUC
#' # Normally boot.times should be set to 200 or more,
#' # we set it to 3 here only to save example running time.
#' val.boot <- validate(
#'   x, time, event,
#'   model.type = "lasso",
#'   alpha = 1, lambda = fit$lasso_best_lambda,
#'   method = "bootstrap", boot.times = 3,
#'   tauc.type = "UNO", tauc.time = seq(0.25, 2, 0.25) * 365,
#'   seed = 1010
#' )
#'
#' # Model validation by 5-fold cross-validation with time-dependent AUC
#' val.cv <- validate(
#'   x, time, event,
#'   model.type = "lasso",
#'   alpha = 1, lambda = fit$lasso_best_lambda,
#'   method = "cv", nfolds = 5,
#'   tauc.type = "UNO", tauc.time = seq(0.25, 2, 0.25) * 365,
#'   seed = 1010
#' )
#'
#' # Model validation by repeated cross-validation with time-dependent AUC
#' val.repcv <- validate(
#'   x, time, event,
#'   model.type = "lasso",
#'   alpha = 1, lambda = fit$lasso_best_lambda,
#'   method = "repeated.cv", nfolds = 5, rep.times = 3,
#'   tauc.type = "UNO", tauc.time = seq(0.25, 2, 0.25) * 365,
#'   seed = 1010
#' )
#'
#' # bootstrap-based discrimination curves has a very narrow band
#' print(val.boot)
#' summary(val.boot)
#' plot(val.boot)
#'
#' # k-fold cv provides a more strict evaluation than bootstrap
#' print(val.cv)
#' summary(val.cv)
#' plot(val.cv)
#'
#' # repeated cv provides similar results as k-fold cv
#' # but more robust than k-fold cv
#' print(val.repcv)
#' summary(val.repcv)
#' plot(val.repcv)
#' # # Test fused lasso, SCAD, and Mnet models ###
#' # library("hdnom")
#' # library("survival")
#' #
#' # # Load imputed SMART data
#' # data(smart)
#' # x = as.matrix(smart[, -c(1, 2)])[1:500,]
#' # time = smart$TEVENT[1:500]
#' # event = smart$EVENT[1:500]
#' # y = Surv(time, event)
#' #
#' # set.seed(1010)
#' # val.boot = validate(
#' #   x, time, event, model.type = "flasso",
#' #   lambda1 = 5, lambda2 = 2,
#' #   method = "bootstrap", boot.times = 10,
#' #   tauc.type = "UNO", tauc.time = seq(0.25, 2, 0.25) * 365,
#' #   seed = 1010)
#' #
#' # val.cv = validate(
#' #   x, time, event, model.type = "scad",
#' #   gamma = 3.7, alpha = 1, lambda = 0.05,
#' #   method = "cv", nfolds = 5,
#' #   tauc.type = "UNO", tauc.time = seq(0.25, 2, 0.25) * 365,
#' #   seed = 1010)
#' #
#' # val.repcv = validate(
#' #   x, time, event, model.type = "mnet",
#' #   gamma = 3, alpha = 0.3, lambda = 0.05,
#' #   method = "repeated.cv", nfolds = 5, rep.times = 3,
#' #   tauc.type = "UNO", tauc.time = seq(0.25, 2, 0.25) * 365,
#' #   seed = 1010)
#' #
#' # print(val.boot)
#' # summary(val.boot)
#' # plot(val.boot)
#' #
#' # print(val.cv)
#' # summary(val.cv)
#' # plot(val.cv)
#' #
#' # print(val.repcv)
#' # summary(val.repcv)
#' # plot(val.repcv)
validate <- function(
  x, time, event,
  model.type =
    c(
      "lasso", "alasso", "flasso", "enet", "aenet",
      "mcp", "mnet", "scad", "snet"
    ),
  alpha, lambda, pen.factor = NULL, gamma,
  lambda1, lambda2,
  method = c("bootstrap", "cv", "repeated.cv"),
  boot.times = NULL, nfolds = NULL, rep.times = NULL,
  tauc.type = c("CD", "SZ", "UNO"), tauc.time,
  seed = 1001, trace = TRUE) {
  model.type <- match.arg(model.type)
  method <- match.arg(method)
  tauc.type <- match.arg(tauc.type)

  set.seed(seed)

  if (method == "bootstrap") {
    if (!is.null(nfolds) || !is.null(rep.times)) {
      stop('nfolds and rep.times must be NULL when method = "bootstrap"')
    }

    if (is.null(boot.times)) stop("please specify boot.times")

    # generate bootstrap sample index
    samp_mat <- matrix(NA, nrow(x), boot.times)
    for (i in 1L:boot.times) samp_mat[, i] <- sample(1L:nrow(x), replace = TRUE)

    # bootstrap validation main loop
    tauc <- vector("list", boot.times)
    for (i in 1L:ncol(samp_mat)) {
      if (trace) cat("Start bootstrap sample", i, "\n")

      samp_idx <- samp_mat[, i]
      x_tr <- x[samp_idx, , drop = FALSE]
      time_tr <- time[samp_idx]
      event_tr <- event[samp_idx]
      y_tr <- Surv(time_tr, event_tr)
      x_te <- x # use original dataset as test set
      y_te <- Surv(time, event) # use original dataset as test set

      if (model.type %in% c("lasso", "alasso", "enet", "aenet")) {
        tauc[[i]] <-
          glmnet_validate_tauc(
            x_tr = x_tr, x_te = x_te, y_tr = y_tr, y_te = y_te,
            alpha = alpha, lambda = lambda, pen.factor = pen.factor,
            tauc.type = tauc.type, tauc.time = tauc.time
          )
      }

      if (model.type %in% c("mcp", "mnet", "scad", "snet")) {
        tauc[[i]] <-
          ncvreg_validate_tauc(
            x_tr = x_tr, x_te = x_te, y_tr = y_tr, y_te = y_te,
            model.type = model.type,
            gamma = gamma, alpha = alpha, lambda = lambda,
            tauc.type = tauc.type, tauc.time = tauc.time
          )
      }

      if (model.type %in% c("flasso")) {
        tauc[[i]] <-
          penalized_validate_tauc(
            x_tr = x_tr, x_te = x_te, y_tr = y_tr, y_te = y_te,
            lambda1 = lambda1, lambda2 = lambda2,
            tauc.type = tauc.type, tauc.time = tauc.time
          )
      }
    }
  } else if (method == "cv") {
    if (!is.null(boot.times) || !is.null(rep.times)) {
      stop('boot.times and rep.times must be NULL when method = "cv"')
    }

    if (is.null(nfolds)) stop("please specify nfolds")

    # generate cross-validation sample index
    row_x <- nrow(x)
    samp_idx <- sample(rep_len(1L:nfolds, row_x))

    # cross-validation main loop
    tauc <- vector("list", nfolds)
    for (i in 1L:nfolds) {
      if (trace) cat("Start fold", i, "\n")

      x_tr <- x[samp_idx != i, , drop = FALSE]
      time_tr <- time[samp_idx != i]
      event_tr <- event[samp_idx != i]
      y_tr <- Surv(time_tr, event_tr)
      x_te <- x[samp_idx == i, , drop = FALSE]
      time_te <- time[samp_idx == i]
      event_te <- event[samp_idx == i]
      y_te <- Surv(time_te, event_te)

      if (model.type %in% c("lasso", "alasso", "enet", "aenet")) {
        tauc[[i]] <-
          glmnet_validate_tauc(
            x_tr = x_tr, x_te = x_te, y_tr = y_tr, y_te = y_te,
            alpha = alpha, lambda = lambda, pen.factor = pen.factor,
            tauc.type = tauc.type, tauc.time = tauc.time
          )
      }

      if (model.type %in% c("mcp", "mnet", "scad", "snet")) {
        tauc[[i]] <-
          ncvreg_validate_tauc(
            x_tr = x_tr, x_te = x_te, y_tr = y_tr, y_te = y_te,
            model.type = model.type,
            gamma = gamma, alpha = alpha, lambda = lambda,
            tauc.type = tauc.type, tauc.time = tauc.time
          )
      }

      if (model.type %in% c("flasso")) {
        tauc[[i]] <-
          penalized_validate_tauc(
            x_tr = x_tr, x_te = x_te, y_tr = y_tr, y_te = y_te,
            lambda1 = lambda1, lambda2 = lambda2,
            tauc.type = tauc.type, tauc.time = tauc.time
          )
      }
    }
  } else if (method == "repeated.cv") {
    if (!is.null(boot.times)) {
      stop('boot.times must be NULL when method = "repeated.cv"')
    }

    if (is.null(nfolds) || is.null(rep.times)) {
      stop("please specify nfolds and rep.times")
    }

    # generate repeated cross-validation sample index list
    row_x <- nrow(x)
    samp_idx <- vector("list", rep.times)
    for (k in 1L:rep.times) samp_idx[[k]] <- sample(rep_len(1L:nfolds, row_x))

    # repeated cross-validation main loop
    tauc <- vector("list", rep.times)
    for (k in 1L:rep.times) tauc[[k]] <- vector("list", nfolds)

    for (j in 1L:rep.times) {
      for (i in 1L:nfolds) {
        if (trace) cat("Start repeat round", j, "fold", i, "\n")

        x_tr <- x[samp_idx[[j]] != i, , drop = FALSE]
        time_tr <- time[samp_idx[[j]] != i]
        event_tr <- event[samp_idx[[j]] != i]
        y_tr <- Surv(time_tr, event_tr)
        x_te <- x[samp_idx[[j]] == i, , drop = FALSE]
        time_te <- time[samp_idx[[j]] == i]
        event_te <- event[samp_idx[[j]] == i]
        y_te <- Surv(time_te, event_te)

        if (model.type %in% c("lasso", "alasso", "enet", "aenet")) {
          tauc[[j]][[i]] <-
            glmnet_validate_tauc(
              x_tr = x_tr, x_te = x_te, y_tr = y_tr, y_te = y_te,
              alpha = alpha, lambda = lambda, pen.factor = pen.factor,
              tauc.type = tauc.type, tauc.time = tauc.time
            )
        }

        if (model.type %in% c("mcp", "mnet", "scad", "snet")) {
          tauc[[j]][[i]] <-
            ncvreg_validate_tauc(
              x_tr = x_tr, x_te = x_te, y_tr = y_tr, y_te = y_te,
              model.type = model.type,
              gamma = gamma, alpha = alpha, lambda = lambda,
              tauc.type = tauc.type, tauc.time = tauc.time
            )
        }

        if (model.type %in% c("flasso")) {
          tauc[[j]][[i]] <-
            penalized_validate_tauc(
              x_tr = x_tr, x_te = x_te, y_tr = y_tr, y_te = y_te,
              lambda1 = lambda1, lambda2 = lambda2,
              tauc.type = tauc.type, tauc.time = tauc.time
            )
        }
      }
    }
  } else {
    stop('method must be one of "bootstrap", cv", or "repeated.cv"')
  }

  switch(

    method,

    bootstrap = {
      if (model.type %in% c("lasso", "alasso", "enet", "aenet")) {
        class(tauc) <- c(
          "hdnom.validate",
          "glmnet.validate.bootstrap"
        )
        attr(tauc, "model.type") <- model.type
        attr(tauc, "alpha") <- alpha
        attr(tauc, "lambda") <- lambda
        attr(tauc, "pen.factor") <- pen.factor
        attr(tauc, "boot.times") <- boot.times
        attr(tauc, "tauc.type") <- tauc.type
        attr(tauc, "tauc.time") <- tauc.time
        attr(tauc, "seed") <- seed
      }

      if (model.type %in% c("mcp", "mnet", "scad", "snet")) {
        class(tauc) <- c(
          "hdnom.validate",
          "ncvreg.validate.bootstrap"
        )
        attr(tauc, "model.type") <- model.type
        attr(tauc, "gamma") <- gamma
        attr(tauc, "alpha") <- alpha
        attr(tauc, "lambda") <- lambda
        attr(tauc, "boot.times") <- boot.times
        attr(tauc, "tauc.type") <- tauc.type
        attr(tauc, "tauc.time") <- tauc.time
        attr(tauc, "seed") <- seed
      }

      if (model.type %in% c("flasso")) {
        class(tauc) <- c(
          "hdnom.validate",
          "penalized.validate.bootstrap"
        )
        attr(tauc, "model.type") <- model.type
        attr(tauc, "lambda1") <- lambda1
        attr(tauc, "lambda2") <- lambda2
        attr(tauc, "boot.times") <- boot.times
        attr(tauc, "tauc.type") <- tauc.type
        attr(tauc, "tauc.time") <- tauc.time
        attr(tauc, "seed") <- seed
      }
    },

    cv = {
      if (model.type %in% c("lasso", "alasso", "enet", "aenet")) {
        class(tauc) <- c(
          "hdnom.validate",
          "glmnet.validate.cv"
        )
        attr(tauc, "model.type") <- model.type
        attr(tauc, "alpha") <- alpha
        attr(tauc, "lambda") <- lambda
        attr(tauc, "pen.factor") <- pen.factor
        attr(tauc, "nfolds") <- nfolds
        attr(tauc, "tauc.type") <- tauc.type
        attr(tauc, "tauc.time") <- tauc.time
        attr(tauc, "seed") <- seed
      }

      if (model.type %in% c("mcp", "mnet", "scad", "snet")) {
        class(tauc) <- c(
          "hdnom.validate",
          "ncvreg.validate.cv"
        )
        attr(tauc, "model.type") <- model.type
        attr(tauc, "gamma") <- gamma
        attr(tauc, "alpha") <- alpha
        attr(tauc, "lambda") <- lambda
        attr(tauc, "nfolds") <- nfolds
        attr(tauc, "tauc.type") <- tauc.type
        attr(tauc, "tauc.time") <- tauc.time
        attr(tauc, "seed") <- seed
      }

      if (model.type %in% c("flasso")) {
        class(tauc) <- c(
          "hdnom.validate",
          "penalized.validate.cv"
        )
        attr(tauc, "model.type") <- model.type
        attr(tauc, "lambda1") <- lambda1
        attr(tauc, "lambda2") <- lambda2
        attr(tauc, "nfolds") <- nfolds
        attr(tauc, "tauc.type") <- tauc.type
        attr(tauc, "tauc.time") <- tauc.time
        attr(tauc, "seed") <- seed
      }
    },

    repeated.cv = {
      if (model.type %in% c("lasso", "alasso", "enet", "aenet")) {
        class(tauc) <- c(
          "hdnom.validate",
          "glmnet.validate.repeated.cv"
        )
        attr(tauc, "model.type") <- model.type
        attr(tauc, "alpha") <- alpha
        attr(tauc, "lambda") <- lambda
        attr(tauc, "pen.factor") <- pen.factor
        attr(tauc, "nfolds") <- nfolds
        attr(tauc, "rep.times") <- rep.times
        attr(tauc, "tauc.type") <- tauc.type
        attr(tauc, "tauc.time") <- tauc.time
        attr(tauc, "seed") <- seed
      }

      if (model.type %in% c("mcp", "mnet", "scad", "snet")) {
        class(tauc) <- c(
          "hdnom.validate",
          "ncvreg.validate.repeated.cv"
        )
        attr(tauc, "model.type") <- model.type
        attr(tauc, "gamma") <- gamma
        attr(tauc, "alpha") <- alpha
        attr(tauc, "lambda") <- lambda
        attr(tauc, "nfolds") <- nfolds
        attr(tauc, "rep.times") <- rep.times
        attr(tauc, "tauc.type") <- tauc.type
        attr(tauc, "tauc.time") <- tauc.time
        attr(tauc, "seed") <- seed
      }

      if (model.type %in% c("flasso")) {
        class(tauc) <- c(
          "hdnom.validate",
          "penalized.validate.repeated.cv"
        )
        attr(tauc, "model.type") <- model.type
        attr(tauc, "lambda1") <- lambda1
        attr(tauc, "lambda2") <- lambda2
        attr(tauc, "nfolds") <- nfolds
        attr(tauc, "rep.times") <- rep.times
        attr(tauc, "tauc.type") <- tauc.type
        attr(tauc, "tauc.time") <- tauc.time
        attr(tauc, "seed") <- seed
      }
    }
  )

  tauc
}
