#' Calibrate high-dimensional Cox models
#'
#' Calibrate high-dimensional Cox models
#'
#' @param x Matrix of training data used for fitting the model;
#' on which to run the calibration.
#' @param time Survival time.
#' Must be of the same length with the number of rows as \code{x}.
#' @param event Status indicator, normally 0 = alive, 1 = dead.
#' Must be of the same length with the number of rows as \code{x}.
#' @param model.type Model type to calibrate. Could be one of \code{"lasso"},
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
#' model fits on the resampled data. From the Cox model you have built.
#' @param pen.factor Penalty factors to apply to each coefficient.
#' From the built \emph{adaptive lasso} or \emph{adaptive elastic-net} model.
#' @param gamma Value of the model parameter gamma for
#' MCP/SCAD/Mnet/Snet models.
#' @param lambda1 Value of the penalty parameter lambda1 for fused lasso model.
#' @param lambda2 Value of the penalty parameter lambda2 for fused lasso model.
#' @param method Calibration method.
#' Options including \code{"fitting"}, \code{"bootstrap"}, \code{"cv"},
#' and \code{"repeated.cv"}.
#' @param boot.times Number of repetitions for bootstrap.
#' @param nfolds Number of folds for cross-validation and
#' repeated cross-validation.
#' @param rep.times Number of repeated times for repeated cross-validation.
#' @param pred.at Time point at which calibration should take place.
#' @param ngroup Number of groups to be formed for calibration.
#' @param seed A random seed for resampling.
#' @param trace Logical. Output the calibration progress or not.
#' Default is \code{TRUE}.
#'
#' @importFrom survival Surv
#' @importFrom stats quantile
#' @importFrom stats median
#'
#' @export calibrate
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
#' fit <- fit_lasso(x, y, nfolds = 5, rule = "lambda.1se", seed = 11)
#'
#' # Model calibration by fitting the original data directly
#' cal.fitting <- calibrate(
#'   x, time, event,
#'   model.type = "lasso",
#'   alpha = 1, lambda = fit$lasso_best_lambda,
#'   method = "fitting",
#'   pred.at = 365 * 9, ngroup = 5,
#'   seed = 1010
#' )
#'
#' # Model calibration by 5-fold cross-validation
#' cal.cv <- calibrate(
#'   x, time, event,
#'   model.type = "lasso",
#'   alpha = 1, lambda = fit$lasso_best_lambda,
#'   method = "cv", nfolds = 5,
#'   pred.at = 365 * 9, ngroup = 5,
#'   seed = 1010
#' )
#'
#' print(cal.fitting)
#' summary(cal.fitting)
#' plot(cal.fitting)
#'
#' print(cal.cv)
#' summary(cal.cv)
#' plot(cal.cv)
#'
#' # ### Testing fused lasso, SCAD, and Mnet models ###
#' # library("hdnom")
#' # library("survival")
#' #
#' # # Load imputed SMART data
#' # data(smart)
#' # x = as.matrix(smart[, -c(1, 2)])[1:500, ]
#' # time = smart$TEVENT[1:500]
#' # event = smart$EVENT[1:500]
#' # y = Surv(time, event)
#' #
#' # set.seed(1010)
#' # cal.fitting = calibrate(
#' #   x, time, event, model.type = "flasso",
#' #   lambda1 = 5, lambda2 = 2,
#' #   method = "fitting",
#' #   pred.at = 365 * 9, ngroup = 5,
#' #   seed = 1010)
#' #
#' # cal.boot = calibrate(
#' #   x, time, event, model.type = "scad",
#' #   gamma = 3.7, alpha = 1, lambda = 0.03,
#' #   method = "bootstrap", boot.times = 10,
#' #   pred.at = 365 * 9, ngroup = 5,
#' #   seed = 1010)
#' #
#' # cal.cv = calibrate(
#' #   x, time, event, model.type = "mnet",
#' #   gamma = 3, alpha = 0.3, lambda = 0.03,
#' #   method = "cv", nfolds = 5,
#' #   pred.at = 365 * 9, ngroup = 5,
#' #   seed = 1010)
#' #
#' # cal.repcv = calibrate(
#' #   x, time, event, model.type = "flasso",
#' #   lambda1 = 5, lambda2 = 2,
#' #   method = "repeated.cv", nfolds = 5, rep.times = 3,
#' #   pred.at = 365 * 9, ngroup = 5,
#' #   seed = 1010)
#' #
#' # print(cal.fitting)
#' # summary(cal.fitting)
#' # plot(cal.fitting)
#' #
#' # print(cal.boot)
#' # summary(cal.boot)
#' # plot(cal.boot)
#' #
#' # print(cal.cv)
#' # summary(cal.cv)
#' # plot(cal.cv)
#' #
#' # print(cal.repcv)
#' # summary(cal.repcv)
#' # plot(cal.repcv)
calibrate <- function(
  x, time, event,
  model.type =
    c(
      "lasso", "alasso", "flasso", "enet", "aenet",
      "mcp", "mnet", "scad", "snet"
    ),
  alpha, lambda, pen.factor = NULL, gamma,
  lambda1, lambda2,
  method = c("fitting", "bootstrap", "cv", "repeated.cv"),
  boot.times = NULL, nfolds = NULL, rep.times = NULL,
  pred.at, ngroup = 5,
  seed = 1001, trace = TRUE) {
  model.type <- match.arg(model.type)
  method <- match.arg(method)
  if (length(pred.at) != 1L) stop("pred.at should only contain 1 time point")

  set.seed(seed)

  if (method == "fitting") {
    if (!is.null(boot.times) || !is.null(nfolds) || !is.null(rep.times)) {
      stop('boot.times, nfolds, and rep.times must be NULL when method = "fitting"')
    }

    if (trace) cat("Start fitting ...\n")

    if (model.type %in% c("lasso", "alasso", "enet", "aenet")) {
      pred_list <- glmnet_calibrate_surv_prob_pred(
        x_tr = x, x_te = x, y_tr = Surv(time, event),
        alpha = alpha, lambda = lambda, pen.factor = pen.factor,
        pred.at = pred.at
      )
    }

    if (model.type %in% c("mcp", "mnet", "scad", "snet")) {
      pred_list <- ncvreg_calibrate_surv_prob_pred(
        x_tr = x, x_te = x, y_tr = Surv(time, event),
        model.type = model.type,
        gamma = gamma, alpha = alpha, lambda = lambda,
        pred.at = pred.at
      )
    }

    if (model.type %in% c("flasso")) {
      pred_list <- penalized_calibrate_surv_prob_pred(
        x_tr = x, x_te = x, y_tr = Surv(time, event),
        lambda1 = lambda1, lambda2 = lambda2,
        pred.at = pred.at
      )
    }

    pred_prob <- rep(NA, nrow(x))
    for (i in 1L:length(pred_prob)) pred_prob[i] <- pred_list$p[i, pred_list$idx]
    grp <- cut(pred_prob, quantile(pred_prob, seq(0, 1, 1 / ngroup)), labels = 1L:ngroup)

    pred_prob_median <- tapply(pred_prob, grp, median)

    true_prob <- calibrate_surv_prob_true(
      pred_prob, grp, time, event, pred.at, ngroup
    )

    prob <- cbind(pred_prob_median, true_prob)
    colnames(prob)[1L] <- "Predicted"
  } else if (method == "bootstrap") {
    if (!is.null(nfolds) || !is.null(rep.times)) {
      stop('nfolds and rep.times must be NULL when method = "bootstrap"')
    }

    if (is.null(boot.times)) stop("please specify boot.times")

    # generate bootstrap sample index
    samp_mat <- matrix(NA, nrow(x), boot.times)
    for (i in 1L:boot.times) samp_mat[, i] <- sample(1L:nrow(x), replace = TRUE)

    # bootstrap validation main loop
    pred_list_list <- vector("list", boot.times)
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
        pred_list_list[[i]] <- glmnet_calibrate_surv_prob_pred(
          x_tr = x_tr, x_te = x_te, y_tr = y_tr,
          alpha = alpha, lambda = lambda, pen.factor = pen.factor,
          pred.at = pred.at
        )
      }

      if (model.type %in% c("mcp", "mnet", "scad", "snet")) {
        pred_list_list[[i]] <- ncvreg_calibrate_surv_prob_pred(
          x_tr = x_tr, x_te = x_te, y_tr = y_tr,
          model.type = model.type,
          alpha = alpha, lambda = lambda, gamma = gamma,
          pred.at = pred.at
        )
      }

      if (model.type %in% c("flasso")) {
        pred_list_list[[i]] <- penalized_calibrate_surv_prob_pred(
          x_tr = x_tr, x_te = x_te, y_tr = y_tr,
          lambda1 = lambda1, lambda2 = lambda2,
          pred.at = pred.at
        )
      }
    }

    # get prediction from each bootstrap round
    pred_list_list_tmp <- vector("list", boot.times)
    for (i in 1L:boot.times) {
      pred_list_list_tmp[[i]] <- pred_list_list[[i]]$p[, pred_list_list[[i]]$idx]
    }

    # get mean of predicted probability across bootstrap rounds
    pred_prob <- Reduce("+", pred_list_list_tmp) / length(pred_list_list_tmp)
    grp <- cut(pred_prob, quantile(pred_prob, seq(0, 1, 1 / ngroup)), labels = 1L:ngroup)

    pred_prob_median <- tapply(pred_prob, grp, median)

    true_prob <- calibrate_surv_prob_true(
      pred_prob, grp, time, event, pred.at, ngroup
    )

    prob <- cbind(pred_prob_median, true_prob)
    colnames(prob)[1L] <- "Predicted"
  } else if (method == "cv") {
    if (!is.null(boot.times) || !is.null(rep.times)) {
      stop('boot.times and rep.times must be NULL when method = "cv"')
    }

    if (is.null(nfolds)) stop("please specify nfolds")

    # generate cross-validation sample index
    row_x <- nrow(x)
    samp_idx <- sample(rep_len(1L:nfolds, row_x))

    # cross-validation main loop
    pred_list_list <- vector("list", nfolds)
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
      idx_te <- which(samp_idx == i)

      if (model.type %in% c("lasso", "alasso", "enet", "aenet")) {
        pred_list_list[[i]] <- glmnet_calibrate_surv_prob_pred(
          x_tr = x_tr, x_te = x_te, y_tr = y_tr,
          alpha = alpha, lambda = lambda, pen.factor = pen.factor,
          pred.at = pred.at
        )
      }

      if (model.type %in% c("mcp", "mnet", "scad", "snet")) {
        pred_list_list[[i]] <- ncvreg_calibrate_surv_prob_pred(
          x_tr = x_tr, x_te = x_te, y_tr = y_tr,
          model.type = model.type,
          alpha = alpha, lambda = lambda, gamma = gamma,
          pred.at = pred.at
        )
      }

      if (model.type %in% c("flasso")) {
        pred_list_list[[i]] <- penalized_calibrate_surv_prob_pred(
          x_tr = x_tr, x_te = x_te, y_tr = y_tr,
          lambda1 = lambda1, lambda2 = lambda2,
          pred.at = pred.at
        )
      }

      rownames(pred_list_list[[i]]$p) <- idx_te
    }

    # get prediction from each fold
    pred_list_list_tmp <- vector("list", nfolds)
    for (i in 1L:nfolds) {
      pred_list_list_tmp[[i]] <- pred_list_list[[i]]$p[, pred_list_list[[i]]$idx]
    }

    # combine and sort all predictions
    pred_prob_unsorted <- unlist(pred_list_list_tmp)
    pred_prob <- pred_prob_unsorted[order(as.integer(names(pred_prob_unsorted)))]
    grp <- cut(pred_prob, quantile(pred_prob, seq(0, 1, 1 / ngroup)), labels = 1L:ngroup)

    pred_prob_median <- tapply(pred_prob, grp, median)

    true_prob <- calibrate_surv_prob_true(
      pred_prob, grp, time, event, pred.at, ngroup
    )

    prob <- cbind(pred_prob_median, true_prob)
    colnames(prob)[1L] <- "Predicted"
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
    pred_list_list <- vector("list", rep.times)
    for (k in 1L:rep.times) pred_list_list[[k]] <- vector("list", nfolds)

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
        idx_te <- which(samp_idx[[j]] == i)

        if (model.type %in% c("lasso", "alasso", "enet", "aenet")) {
          pred_list_list[[j]][[i]] <- glmnet_calibrate_surv_prob_pred(
            x_tr = x_tr, x_te = x_te, y_tr = y_tr,
            alpha = alpha, lambda = lambda, pen.factor = pen.factor,
            pred.at = pred.at
          )
        }

        if (model.type %in% c("mcp", "mnet", "scad", "snet")) {
          pred_list_list[[j]][[i]] <- ncvreg_calibrate_surv_prob_pred(
            x_tr = x_tr, x_te = x_te, y_tr = y_tr,
            model.type = model.type,
            alpha = alpha, lambda = lambda, gamma = gamma,
            pred.at = pred.at
          )
        }

        if (model.type %in% c("flasso")) {
          pred_list_list[[j]][[i]] <- penalized_calibrate_surv_prob_pred(
            x_tr = x_tr, x_te = x_te, y_tr = y_tr,
            lambda1 = lambda1, lambda2 = lambda2,
            pred.at = pred.at
          )
        }

        rownames(pred_list_list[[j]][[i]]$p) <- idx_te
      }
    }

    # get prediction from each fold for each repeat
    pred_list_list_tmp <- vector("list", rep.times)
    for (i in 1L:rep.times) pred_list_list_tmp[[i]] <- vector("list", nfolds)

    for (j in 1L:rep.times) {
      for (i in 1L:nfolds) {
        pred_list_list_tmp[[j]][[i]] <-
          pred_list_list[[j]][[i]]$p[, pred_list_list[[j]][[i]]$idx]
      }
    }

    # combine and sort all predictions for each repeat time
    pred_prob_unsorted <- vector("list", rep.times)
    for (i in 1L:rep.times) pred_prob_unsorted[[i]] <- unlist(pred_list_list_tmp[[i]])
    pred_prob_tmp <- vector("list", rep.times)
    for (i in 1L:rep.times) {
      pred_prob_tmp[[i]] <-
        pred_prob_unsorted[[i]][order(as.integer(names(pred_prob_unsorted[[i]])))]
    }

    # get mean of predicted probability across repeat times
    pred_prob <- Reduce("+", pred_prob_tmp) / length(pred_prob_tmp)
    grp <- cut(pred_prob, quantile(pred_prob, seq(0, 1, 1 / ngroup)), labels = 1L:ngroup)

    pred_prob_median <- tapply(pred_prob, grp, median)

    true_prob <- calibrate_surv_prob_true(
      pred_prob, grp, time, event, pred.at, ngroup
    )

    prob <- cbind(pred_prob_median, true_prob)
    colnames(prob)[1L] <- "Predicted"
  } else {
    stop('method must be one of "fitting", "bootstrap", "cv", or "repeated.cv"')
  }

  switch(

    method,

    fitting = {
      if (model.type %in% c("lasso", "alasso", "enet", "aenet")) {
        class(prob) <- c(
          "hdnom.calibrate",
          "glmnet.calibrate.fitting"
        )
        attr(prob, "alpha") <- alpha
        attr(prob, "lambda") <- lambda
        attr(prob, "pen.factor") <- pen.factor
      }

      if (model.type %in% c("mcp", "mnet", "scad", "snet")) {
        class(prob) <- c(
          "hdnom.calibrate",
          "ncvreg.calibrate.fitting"
        )
        attr(prob, "alpha") <- alpha
        attr(prob, "lambda") <- lambda
        attr(prob, "gamma") <- gamma
      }

      if (model.type %in% c("flasso")) {
        class(prob) <- c(
          "hdnom.calibrate",
          "penalized.calibrate.fitting"
        )
        attr(prob, "lambda1") <- lambda1
        attr(prob, "lambda2") <- lambda2
      }
    },

    bootstrap = {
      if (model.type %in% c("lasso", "alasso", "enet", "aenet")) {
        class(prob) <- c(
          "hdnom.calibrate",
          "glmnet.calibrate.bootstrap"
        )
        attr(prob, "alpha") <- alpha
        attr(prob, "lambda") <- lambda
        attr(prob, "pen.factor") <- pen.factor
        attr(prob, "boot.times") <- boot.times
      }

      if (model.type %in% c("mcp", "mnet", "scad", "snet")) {
        class(prob) <- c(
          "hdnom.calibrate",
          "ncvreg.calibrate.bootstrap"
        )
        attr(prob, "alpha") <- alpha
        attr(prob, "lambda") <- lambda
        attr(prob, "gamma") <- gamma
        attr(prob, "boot.times") <- boot.times
      }

      if (model.type %in% c("flasso")) {
        class(prob) <- c(
          "hdnom.calibrate",
          "penalized.calibrate.bootstrap"
        )
        attr(prob, "lambda1") <- lambda1
        attr(prob, "lambda2") <- lambda2
        attr(prob, "boot.times") <- boot.times
      }
    },

    cv = {
      if (model.type %in% c("lasso", "alasso", "enet", "aenet")) {
        class(prob) <- c(
          "hdnom.calibrate",
          "glmnet.calibrate.cv"
        )
        attr(prob, "alpha") <- alpha
        attr(prob, "lambda") <- lambda
        attr(prob, "pen.factor") <- pen.factor
        attr(prob, "nfolds") <- nfolds
      }

      if (model.type %in% c("mcp", "mnet", "scad", "snet")) {
        class(prob) <- c(
          "hdnom.calibrate",
          "ncvreg.calibrate.cv"
        )
        attr(prob, "alpha") <- alpha
        attr(prob, "lambda") <- lambda
        attr(prob, "gamma") <- gamma
        attr(prob, "nfolds") <- nfolds
      }

      if (model.type %in% c("flasso")) {
        class(prob) <- c(
          "hdnom.calibrate",
          "penalized.calibrate.cv"
        )
        attr(prob, "lambda1") <- lambda1
        attr(prob, "lambda2") <- lambda2
        attr(prob, "nfolds") <- nfolds
      }
    },

    repeated.cv = {
      if (model.type %in% c("lasso", "alasso", "enet", "aenet")) {
        class(prob) <- c(
          "hdnom.calibrate",
          "glmnet.calibrate.repeated.cv"
        )
        attr(prob, "alpha") <- alpha
        attr(prob, "lambda") <- lambda
        attr(prob, "pen.factor") <- pen.factor
        attr(prob, "nfolds") <- nfolds
        attr(prob, "rep.times") <- rep.times
      }

      if (model.type %in% c("mcp", "mnet", "scad", "snet")) {
        class(prob) <- c(
          "hdnom.calibrate",
          "ncvreg.calibrate.repeated.cv"
        )
        attr(prob, "alpha") <- alpha
        attr(prob, "lambda") <- lambda
        attr(prob, "gamma") <- gamma
        attr(prob, "nfolds") <- nfolds
        attr(prob, "rep.times") <- rep.times
      }

      if (model.type %in% c("flasso")) {
        class(prob) <- c(
          "hdnom.calibrate",
          "penalized.calibrate.repeated.cv"
        )
        attr(prob, "lambda1") <- lambda1
        attr(prob, "lambda2") <- lambda2
        attr(prob, "nfolds") <- nfolds
        attr(prob, "rep.times") <- rep.times
      }
    }
  )

  attr(prob, "model.type") <- model.type
  attr(prob, "pred.at") <- pred.at
  attr(prob, "ngroup") <- ngroup
  attr(prob, "risk.group") <- grp
  attr(prob, "surv.time") <- time
  attr(prob, "surv.event") <- event
  attr(prob, "seed") <- seed

  prob
}
