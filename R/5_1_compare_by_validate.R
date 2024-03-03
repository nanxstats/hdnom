#' Compare high-dimensional Cox models by model validation
#'
#' Compare high-dimensional Cox models by model validation
#'
#' @param x Matrix of training data used for fitting the model;
#' on which to run the validation.
#' @param time Survival time.
#' Must be of the same length with the number of rows as \code{x}.
#' @param event Status indicator, normally 0 = alive, 1 = dead.
#' Must be of the same length with the number of rows as \code{x}.
#' @param model.type Model types to compare. Could be at least two
#' of \code{"lasso"}, \code{"alasso"}, \code{"flasso"}, \code{"enet"},
#' \code{"aenet"}, \code{"mcp"}, \code{"mnet"}, \code{"scad"},
#' or \code{"snet"}.
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
#' @param seed A random seed for cross-validation fold division.
#' @param trace Logical. Output the validation progress or not.
#' Default is \code{TRUE}.
#'
#' @export compare_by_validate
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
#' data(smart)
#' x <- as.matrix(smart[, -c(1, 2)])[1:1000, ]
#' time <- smart$TEVENT[1:1000]
#' event <- smart$EVENT[1:1000]
#'
#' # Compare lasso and adaptive lasso by 5-fold cross-validation
#' cmp.val.cv <- compare_by_validate(
#'   x, time, event,
#'   model.type = c("lasso", "alasso"),
#'   method = "cv", nfolds = 5, tauc.type = "UNO",
#'   tauc.time = seq(0.25, 2, 0.25) * 365, seed = 1001
#' )
#'
#' print(cmp.val.cv)
#' summary(cmp.val.cv)
#' plot(cmp.val.cv)
#' plot(cmp.val.cv, interval = TRUE)
compare_by_validate <- function(
    x, time, event,
    model.type = c(
      "lasso", "alasso", "flasso", "enet", "aenet",
      "mcp", "mnet", "scad", "snet"
    ),
    method = c("bootstrap", "cv", "repeated.cv"),
    boot.times = NULL, nfolds = NULL, rep.times = NULL,
    tauc.type = c("CD", "SZ", "UNO"), tauc.time,
    seed = 1001, trace = TRUE) {
  method <- match.arg(method)
  tauc.type <- match.arg(tauc.type)

  if (!all(model.type %in% c(
    "lasso", "alasso", "flasso", "enet", "aenet",
    "mcp", "mnet", "scad", "snet"
  ))) {
    stop("Unknown model type(s) specified")
  }

  nmodel <- length(model.type)
  tauclist <- vector("list", nmodel)

  # Sanity check for arguments
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

        tauclist[[i]] <- hdnom::validate(
          x, time, event,
          model.type = "lasso",
          alpha = 1, lambda = cvfit$"lambda",
          method = method, boot.times = boot.times, nfolds = nfolds, rep.times = rep.times,
          tauc.type = tauc.type, tauc.time = tauc.time,
          seed = seed, trace = trace
        )
      },
      alasso = {
        cvfit <- fit_alasso(
          x, Surv(time, event),
          nfolds = 5L,
          rule = "lambda.1se", seed = rep(seed, 2)
        )

        tauclist[[i]] <- hdnom::validate(
          x, time, event,
          model.type = "alasso",
          alpha = 1, lambda = cvfit$"lambda", pen.factor = cvfit$"pen_factor",
          method = method, boot.times = boot.times, nfolds = nfolds, rep.times = rep.times,
          tauc.type = tauc.type, tauc.time = tauc.time,
          seed = seed, trace = trace
        )
      },
      flasso = {
        cvfit <- fit_flasso(x, Surv(time, event), nfolds = 5L, seed = seed)

        tauclist[[i]] <- hdnom::validate(
          x, time, event,
          model.type = "flasso",
          lambda1 = cvfit$"lambda1",
          lambda2 = cvfit$"lambda2",
          method = method, boot.times = boot.times, nfolds = nfolds, rep.times = rep.times,
          tauc.type = tauc.type, tauc.time = tauc.time,
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

        tauclist[[i]] <- hdnom::validate(
          x, time, event,
          model.type = "enet",
          alpha = cvfit$"alpha", lambda = cvfit$"lambda",
          method = method, boot.times = boot.times, nfolds = nfolds, rep.times = rep.times,
          tauc.type = tauc.type, tauc.time = tauc.time,
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

        tauclist[[i]] <- hdnom::validate(
          x, time, event,
          model.type = "aenet",
          alpha = cvfit$"alpha", lambda = cvfit$"lambda", pen.factor = cvfit$"pen_factor",
          method = method, boot.times = boot.times, nfolds = nfolds, rep.times = rep.times,
          tauc.type = tauc.type, tauc.time = tauc.time,
          seed = seed, trace = trace
        )
      },
      mcp = {
        cvfit <- fit_mcp(x, Surv(time, event), nfolds = 5L, seed = seed)

        tauclist[[i]] <- hdnom::validate(
          x, time, event,
          model.type = "mcp",
          alpha = 1, gamma = cvfit$"gamma", lambda = cvfit$"lambda",
          method = method, boot.times = boot.times, nfolds = nfolds, rep.times = rep.times,
          tauc.type = tauc.type, tauc.time = tauc.time,
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

        tauclist[[i]] <- hdnom::validate(
          x, time, event,
          model.type = "mnet",
          alpha = cvfit$"alpha", gamma = cvfit$"gamma", lambda = cvfit$"lambda",
          method = method, boot.times = boot.times, nfolds = nfolds, rep.times = rep.times,
          tauc.type = tauc.type, tauc.time = tauc.time,
          seed = seed, trace = trace
        )
      },
      scad = {
        cvfit <- fit_scad(x, Surv(time, event), nfolds = 5L, seed = seed)

        tauclist[[i]] <- hdnom::validate(
          x, time, event,
          model.type = "scad",
          alpha = 1, gamma = cvfit$"gamma", lambda = cvfit$"lambda",
          method = method, boot.times = boot.times, nfolds = nfolds, rep.times = rep.times,
          tauc.type = tauc.type, tauc.time = tauc.time,
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

        tauclist[[i]] <- hdnom::validate(
          x, time, event,
          model.type = "snet",
          alpha = cvfit$"alpha", gamma = cvfit$"gamma", lambda = cvfit$"lambda",
          method = method, boot.times = boot.times, nfolds = nfolds, rep.times = rep.times,
          tauc.type = tauc.type, tauc.time = tauc.time,
          seed = seed, trace = trace
        )
      }
    )
  }

  names(tauclist) <- model.type
  class(tauclist) <- c("hdnom.compare.validate")

  tauclist
}
