#' Externally calibrate high-dimensional Cox models
#'
#' Externally calibrate high-dimensional Cox models
#'
#' @param object Model object fitted by \code{hdnom::fit_*()}.
#' @param x Matrix of training data used for fitting the model.
#' @param time Survival time of the training data.
#' Must be of the same length with the number of rows as \code{x}.
#' @param event Status indicator of the training data,
#' normally 0 = alive, 1 = dead.
#' Must be of the same length with the number of rows as \code{x}.
#' @param x_new Matrix of predictors for the external validation data.
#' @param time_new Survival time of the external validation data.
#' Must be of the same length with the number of rows as \code{x_new}.
#' @param event_new Status indicator of the external validation data,
#' normally 0 = alive, 1 = dead.
#' Must be of the same length with the number of rows as \code{x_new}.
#' @param pred.at Time point at which external calibration should take place.
#' @param ngroup Number of groups to be formed for external calibration.
#'
#' @importFrom survival Surv
#' @importFrom stats quantile
#' @importFrom stats median
#'
#' @export calibrate_external
#'
#' @examples
#' library("survival")
#'
#' # Load imputed SMART data
#' data(smart)
#' # Use the first 1000 samples as training data
#' # (the data used for internal validation)
#' x <- as.matrix(smart[, -c(1, 2)])[1:1000, ]
#' time <- smart$TEVENT[1:1000]
#' event <- smart$EVENT[1:1000]
#'
#' # Take the next 1000 samples as external calibration data
#' # In practice, usually use data collected in other studies
#' x_new <- as.matrix(smart[, -c(1, 2)])[1001:2000, ]
#' time_new <- smart$TEVENT[1001:2000]
#' event_new <- smart$EVENT[1001:2000]
#'
#' # Fit Cox model with lasso penalty
#' fit <- fit_lasso(
#'   x, Surv(time, event),
#'   nfolds = 5, rule = "lambda.min", seed = 1001
#' )
#'
#' # External calibration
#' cal.ext <- calibrate_external(
#'   fit, x, time, event,
#'   x_new, time_new, event_new,
#'   pred.at = 365 * 5, ngroup = 5
#' )
#'
#' print(cal.ext)
#' summary(cal.ext)
#' plot(cal.ext, xlim = c(0.6, 1), ylim = c(0.6, 1))
#' # # Test fused lasso, MCP, and Snet models
#' # data(smart)
#' # # Use first 500 samples as training data
#' # # (the data used for internal validation)
#' # x <- as.matrix(smart[, -c(1, 2)])[1:500, ]
#' # time <- smart$TEVENT[1:500]
#' # event <- smart$EVENT[1:500]
#' #
#' # # Take 1000 samples as external validation data.
#' # # In practice, usually use data collected in other studies.
#' # x_new <- as.matrix(smart[, -c(1, 2)])[1001:2000, ]
#' # time_new <- smart$TEVENT[1001:2000]
#' # event_new <- smart$EVENT[1001:2000]
#' #
#' # flassofit <- fit_flasso(x, survival::Surv(time, event), nfolds = 5, seed = 11)
#' # scadfit <- fit_mcp(x, survival::Surv(time, event), nfolds = 5, seed = 11)
#' # mnetfit <- fit_snet(x, survival::Surv(time, event), nfolds = 5, seed = 11)
#' #
#' # cal.ext1 <- calibrate_external(
#' #   flassofit, x, time, event,
#' #   x_new, time_new, event_new,
#' #   pred.at = 365 * 5, ngroup = 5)
#' #
#' # cal.ext2 <- calibrate_external(
#' #   scadfit, x, time, event,
#' #   x_new, time_new, event_new,
#' #   pred.at = 365 * 5, ngroup = 5)
#' #
#' # cal.ext3 <- calibrate_external(
#' #   mnetfit, x, time, event,
#' #   x_new, time_new, event_new,
#' #   pred.at = 365 * 5, ngroup = 5)
#' #
#' # print(cal.ext1)
#' # summary(cal.ext1)
#' # plot(cal.ext1)
#' #
#' # print(cal.ext2)
#' # summary(cal.ext2)
#' # plot(cal.ext2)
#' #
#' # print(cal.ext3)
#' # summary(cal.ext3)
#' # plot(cal.ext3)
calibrate_external <- function(
    object, x, time, event,
    x_new, time_new, event_new,
    pred.at, ngroup = 5) {
  if (!("hdnom.model" %in% class(object))) {
    stop('object must be of class "hdnom.model"')
  }

  model <- object$model
  model_type <- object$type

  if (!("matrix" %in% class(x_new))) stop("x_new must be a matrix")

  if (length(pred.at) != 1L) stop("pred.at should only contain 1 time point")

  if (model_type %in% c("lasso", "alasso", "enet", "aenet")) {
    pred_list <- glmnet_calibrate_external_surv_prob_pred(
      model, x, x_new, Surv(time, event), pred.at
    )
  }

  if (model_type %in% c("mcp", "mnet", "scad", "snet")) {
    pred_list <- ncvreg_calibrate_external_surv_prob_pred(
      model, x, x_new, Surv(time, event), pred.at
    )
  }

  if (model_type %in% c("flasso")) {
    pred_list <- penalized_calibrate_external_surv_prob_pred(
      model, x, x_new, Surv(time, event), pred.at
    )
  }

  pred_prob <- rep(NA, nrow(x_new))
  for (i in seq_along(pred_prob)) pred_prob[i] <- pred_list$p[i, pred_list$idx]
  grp <- cut(pred_prob, quantile(pred_prob, seq(0, 1, 1 / ngroup)), labels = 1L:ngroup)

  pred_prob_median <- tapply(pred_prob, grp, median)

  true_prob <- calibrate_external_surv_prob_true(
    pred_prob, grp, time_new, event_new, pred.at, ngroup
  )

  prob <- cbind(pred_prob_median, true_prob)
  colnames(prob)[1L] <- "Predicted"

  if (model_type %in% c("lasso", "alasso", "enet", "aenet")) {
    class(prob) <- c("hdnom.calibrate.external", "glmnet.calibrate.external")
  }

  if (model_type %in% c("mcp", "mnet", "scad", "snet")) {
    class(prob) <- c("hdnom.calibrate.external", "ncvreg.calibrate.external")
  }

  if (model_type %in% c("flasso")) {
    class(prob) <- c("hdnom.calibrate.external", "penalized.calibrate.external")
  }

  attr(prob, "model.type") <- model_type
  attr(prob, "pred.at") <- pred.at
  attr(prob, "ngroup") <- ngroup
  attr(prob, "risk.group") <- grp
  attr(prob, "surv.time") <- time_new
  attr(prob, "surv.event") <- event_new

  prob
}
