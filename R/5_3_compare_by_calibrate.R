#' Compare high-dimensional Cox models by model calibration
#'
#' Compare high-dimensional Cox models by model calibration
#'
#' @param x Matrix of training data used for fitting the model;
#' on which to run the calibration.
#' @param time Survival time.
#' Must be of the same length with the number of rows as \code{x}.
#' @param event Status indicator, normally 0 = alive, 1 = dead.
#' Must be of the same length with the number of rows as \code{x}.
#' @param model.type Model types to compare. Could be at least two of
#' \code{"lasso"}, \code{"alasso"}, \code{"flasso"}, \code{"enet"},
#' \code{"aenet"}, \code{"mcp"}, \code{"mnet"}, \code{"scad"},
#' or \code{"snet"}.
#' @param method Calibration method.
#' Could be \code{"bootstrap"}, \code{"cv"}, or \code{"repeated.cv"}.
#' @param boot.times Number of repetitions for bootstrap.
#' @param nfolds Number of folds for cross-validation and
#' repeated cross-validation.
#' @param rep.times Number of repeated times for repeated cross-validation.
#' @param pred.at Time point at which calibration should take place.
#' @param ngroup Number of groups to be formed for calibration.
#' @param seed A random seed for cross-validation fold division.
#' @param trace Logical. Output the calibration progress or not.
#' Default is \code{TRUE}.
#'
#' @export compare_by_calibrate
#'
#' @examples
#' data(smart)
#' x <- as.matrix(smart[, -c(1, 2)])
#' time <- smart$TEVENT
#' event <- smart$EVENT
#'
#' # Compare lasso and adaptive lasso by 5-fold cross-validation
#' cmp.cal.cv <- compare_by_calibrate(
#'   x, time, event,
#'   model.type = c("lasso", "alasso"),
#'   method = "fitting",
#'   pred.at = 365 * 9, ngroup = 5, seed = 1001
#' )
#'
#' print(cmp.cal.cv)
#' summary(cmp.cal.cv)
#' plot(cmp.cal.cv)
compare_by_calibrate <- function(
    x, time, event,
    model.type = c(
      "lasso", "alasso", "flasso", "enet", "aenet",
      "mcp", "mnet", "scad", "snet"
    ),
    method = c("fitting", "bootstrap", "cv", "repeated.cv"),
    boot.times = NULL, nfolds = NULL, rep.times = NULL,
    pred.at, ngroup = 5,
    seed = 1001, trace = TRUE) {
  method <- match.arg(method)
  if (length(pred.at) != 1L) stop("pred.at should only contain 1 time point")

  if (!all(model.type %in% c(
    "lasso", "alasso", "flasso", "enet", "aenet",
    "mcp", "mnet", "scad", "snet"
  ))) {
    stop("Unknown model type(s) specified")
  }

  nmodel <- length(model.type)
  problist <- vector("list", nmodel)

  # Sanity check for arguments
  if (method == "fitting") {
    if (!is.null(boot.times) || !is.null(nfolds) || !is.null(rep.times)) {
      stop('boot.times, nfolds, and rep.times must be NULL when method = "fitting"')
    }
  }

  if (method == "bootstrap") {
    if (!is.null(nfolds) || !is.null(rep.times)) {
      stop('nfolds and rep.times must be NULL when method = "bootstrap"')
    }
    if (is.null(boot.times)) stop("please specify boot.times")
  }

  if (method == "cv") {
    if (!is.null(boot.times) || !is.null(rep.times)) {
      stop('boot.times and rep.times must be NULL when method = "cv"')
    }
    if (is.null(nfolds)) stop("please specify nfolds")
  }

  if (method == "repeated.cv") {
    if (!is.null(boot.times)) {
      stop('boot.times must be NULL when method = "repeated.cv"')
    }
    if (is.null(nfolds) || is.null(rep.times)) {
      stop("please specify nfolds and rep.times")
    }
  }

  for (i in 1L:nmodel) {
    if (trace) cat("Starting model", i, ":", model.type[i], "\n")

    switch(model.type[i],
      lasso = {
        cvfit <- fit_lasso(
          x, Surv(time, event),
          nfolds = 5L,
          rule = "lambda.1se", seed = seed
        )

        problist[[i]] <- hdnom::calibrate(
          x, time, event,
          model.type = "lasso",
          alpha = 1, lambda = cvfit$"lambda",
          method = method, boot.times = boot.times, nfolds = nfolds, rep.times = rep.times,
          pred.at = pred.at, ngroup = ngroup,
          seed = seed, trace = trace
        )
      },
      alasso = {
        cvfit <- fit_alasso(
          x, Surv(time, event),
          nfolds = 5L,
          rule = "lambda.1se", seed = rep(seed, 2)
        )

        problist[[i]] <- hdnom::calibrate(
          x, time, event,
          model.type = "alasso",
          alpha = 1, lambda = cvfit$"lambda", pen.factor = cvfit$"pen_factor",
          method = method, boot.times = boot.times, nfolds = nfolds, rep.times = rep.times,
          pred.at = pred.at, ngroup = ngroup,
          seed = seed, trace = trace
        )
      },
      flasso = {
        cvfit <- fit_flasso(x, Surv(time, event), nfolds = 5L, seed = seed)

        problist[[i]] <- hdnom::calibrate(
          x, time, event,
          model.type = "flasso",
          lambda1 = cvfit$"lambda1",
          lambda2 = cvfit$"lambda2",
          method = method, boot.times = boot.times, nfolds = nfolds, rep.times = rep.times,
          pred.at = pred.at, ngroup = ngroup,
          seed = seed, trace = trace
        )
      },
      enet = {
        cvfit <- fit_enet(
          x, Surv(time, event),
          nfolds = 5L,
          alphas = c(0.1, 0.25, 0.5, 0.75, 0.9), # to reduce computation time
          rule = "lambda.1se", seed = seed
        )

        problist[[i]] <- hdnom::calibrate(
          x, time, event,
          model.type = "enet",
          alpha = cvfit$"alpha", lambda = cvfit$"lambda",
          method = method, boot.times = boot.times, nfolds = nfolds, rep.times = rep.times,
          pred.at = pred.at, ngroup = ngroup,
          seed = seed, trace = trace
        )
      },
      aenet = {
        cvfit <- fit_aenet(
          x, Surv(time, event),
          nfolds = 5L,
          alphas = c(0.1, 0.25, 0.5, 0.75, 0.9), # to reduce computation time
          rule = "lambda.1se", seed = rep(seed, 2)
        )

        problist[[i]] <- hdnom::calibrate(
          x, time, event,
          model.type = "aenet",
          alpha = cvfit$"alpha", lambda = cvfit$"lambda", pen.factor = cvfit$"pen_factor",
          method = method, boot.times = boot.times, nfolds = nfolds, rep.times = rep.times,
          pred.at = pred.at, ngroup = ngroup,
          seed = seed, trace = trace
        )
      },
      mcp = {
        cvfit <- fit_mcp(x, Surv(time, event), nfolds = 5L, seed = seed)

        problist[[i]] <- hdnom::calibrate(
          x, time, event,
          model.type = "mcp",
          alpha = 1, gamma = cvfit$"gamma", lambda = cvfit$"lambda",
          method = method, boot.times = boot.times, nfolds = nfolds, rep.times = rep.times,
          pred.at = pred.at, ngroup = ngroup,
          seed = seed, trace = trace
        )
      },
      mnet = {
        cvfit <- fit_mnet(
          x, Surv(time, event),
          nfolds = 5L,
          alphas = c(0.1, 0.25, 0.5, 0.75, 0.9), # to reduce computation time
          seed = seed
        )

        problist[[i]] <- hdnom::calibrate(
          x, time, event,
          model.type = "mnet",
          alpha = cvfit$"alpha", gamma = cvfit$"gamma", lambda = cvfit$"lambda",
          method = method, boot.times = boot.times, nfolds = nfolds, rep.times = rep.times,
          pred.at = pred.at, ngroup = ngroup,
          seed = seed, trace = trace
        )
      },
      scad = {
        cvfit <- fit_scad(x, Surv(time, event), nfolds = 5L, seed = seed)

        problist[[i]] <- hdnom::calibrate(
          x, time, event,
          model.type = "scad",
          alpha = 1, gamma = cvfit$"gamma", lambda = cvfit$"lambda",
          method = method, boot.times = boot.times, nfolds = nfolds, rep.times = rep.times,
          pred.at = pred.at, ngroup = ngroup,
          seed = seed, trace = trace
        )
      },
      snet = {
        cvfit <- fit_snet(
          x, Surv(time, event),
          nfolds = 5L,
          alphas = c(0.1, 0.25, 0.5, 0.75, 0.9), # to reduce computation time
          seed = seed
        )

        problist[[i]] <- hdnom::calibrate(
          x, time, event,
          model.type = "snet",
          alpha = cvfit$"alpha", gamma = cvfit$"gamma", lambda = cvfit$"lambda",
          method = method, boot.times = boot.times, nfolds = nfolds, rep.times = rep.times,
          pred.at = pred.at, ngroup = ngroup,
          seed = seed, trace = trace
        )
      }
    )
  }

  names(problist) <- model.type
  class(problist) <- c("hdnom.compare.calibrate")

  problist
}
