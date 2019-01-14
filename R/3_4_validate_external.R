#' Externally validate high-dimensional Cox models with time-dependent AUC
#'
#' Externally validate high-dimensional Cox models with time-dependent AUC
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
#' @param tauc.type Type of time-dependent AUC.
#' Including \code{"CD"} proposed by Chambless and Diao (2006).,
#' \code{"SZ"} proposed by Song and Zhou (2008).,
#' \code{"UNO"} proposed by Uno et al. (2007).
#' @param tauc.time Numeric vector. Time points at which to evaluate
#' the time-dependent AUC.
#'
#' @export validate_external
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
#' # Use the first 1000 samples as training data
#' # (the data used for internal validation)
#' x <- as.matrix(smart[, -c(1, 2)])[1:1000, ]
#' time <- smart$TEVENT[1:1000]
#' event <- smart$EVENT[1:1000]
#'
#' # Take the next 1000 samples as external validation data
#' # In practice, usually use data collected in other studies
#' x_new <- as.matrix(smart[, -c(1, 2)])[1001:2000, ]
#' time_new <- smart$TEVENT[1001:2000]
#' event_new <- smart$EVENT[1001:2000]
#'
#' # Fit Cox model with lasso penalty
#' fit <- fit_lasso(
#'   x, Surv(time, event),
#'   nfolds = 5, rule = "lambda.1se", seed = 11
#' )
#'
#' # External validation with time-dependent AUC
#' val.ext <- validate_external(
#'   fit, x, time, event,
#'   x_new, time_new, event_new,
#'   tauc.type = "UNO",
#'   tauc.time = seq(0.25, 2, 0.25) * 365
#' )
#'
#' print(val.ext)
#' summary(val.ext)
#' plot(val.ext)
#'
#' # ### Testing fused lasso, MCP and Snet models ###
#' # library("survival")
#' #
#' # # Load imputed SMART data
#' # data(smart)
#' # # Use first 600 samples as training data
#' # # (the data used for internal validation)
#' # x = as.matrix(smart[, -c(1, 2)])[1:600, ]
#' # time = smart$TEVENT[1:600]
#' # event = smart$EVENT[1:600]
#' #
#' # # Take 500 samples as external validation data.
#' # # In practice, usually use data collected in other studies.
#' # x_new = as.matrix(smart[, -c(1, 2)])[1001:1500, ]
#' # time_new = smart$TEVENT[1001:1500]
#' # event_new = smart$EVENT[1001:1500]
#' #
#' # flassofit = fit_flasso(x, Surv(time, event), nfolds = 5, seed = 11)
#' # scadfit = fit_mcp(x, Surv(time, event), nfolds = 5, seed = 11)
#' # mnetfit = fit_snet(x, Surv(time, event), nfolds = 5, seed = 11)
#' #
#' # val.ext1 = validate_external(
#' #   flassofit, x, time, event,
#' #   x_new, time_new, event_new,
#' #   tauc.type = "UNO",
#' #   tauc.time = seq(0.25, 2, 0.25) * 365)
#' #
#' # val.ext2 = validate_external(
#' #   scadfit, x, time, event,
#' #   x_new, time_new, event_new,
#' #   tauc.type = "CD",
#' #   tauc.time = seq(0.25, 2, 0.25) * 365)
#' #
#' # val.ext3 = validate_external(
#' #   mnetfit, x, time, event,
#' #   x_new, time_new, event_new,
#' #   tauc.type = "SZ",
#' #   tauc.time = seq(0.25, 2, 0.25) * 365)
#' #
#' # print(val.ext1)
#' # summary(val.ext1)
#' # plot(val.ext1)
#' #
#' # print(val.ext2)
#' # summary(val.ext2)
#' # plot(val.ext2)
#' #
#' # print(val.ext3)
#' # summary(val.ext3)
#' # plot(val.ext3)
validate_external <- function(
  object, x, time, event,
  x_new, time_new, event_new,
  tauc.type = c("CD", "SZ", "UNO"), tauc.time) {
  if (!("hdnom.model" %in% class(object))) {
    stop('object must be of class "hdnom.model"')
  }

  model.type <- gsub("hdnom.model.", "", setdiff(class(object), "hdnom.model"))
  model_object <- object[[paste0(model.type, "_model")]]

  if (!("matrix" %in% class(x_new))) stop("x_new must be a matrix")

  tauc.type <- match.arg(tauc.type)

  x_tr <- x
  time_tr <- time
  event_tr <- event
  y_tr <- Surv(time_tr, event_tr)
  x_te <- x_new
  time_te <- time_new
  event_te <- event_new
  y_te <- Surv(time_te, event_te)

  if (model.type %in% c("lasso", "alasso", "enet", "aenet")) {
    tauc <- glmnet_validate_external_tauc(
      object = model_object, x_tr = x_tr, x_te = x_te, y_tr = y_tr, y_te = y_te,
      tauc.type = tauc.type, tauc.time = tauc.time
    )
  }

  if (model.type %in% c("mcp", "mnet", "scad", "snet")) {
    tauc <- ncvreg_validate_external_tauc(
      object = model_object, x_tr = x_tr, x_te = x_te, y_tr = y_tr, y_te = y_te,
      tauc.type = tauc.type, tauc.time = tauc.time
    )
  }

  if (model.type %in% c("flasso")) {
    tauc <- penalized_validate_external_tauc(
      object = model_object, x_tr = x_tr, x_te = x_te, y_tr = y_tr, y_te = y_te,
      tauc.type = tauc.type, tauc.time = tauc.time
    )
  }

  if (model.type %in% c("lasso", "alasso", "enet", "aenet")) {
    class(tauc) <- c(
      "hdnom.validate.external",
      "glmnet.validate.external"
    )
    attr(tauc, "model.type") <- model.type
    attr(tauc, "tauc.type") <- tauc.type
    attr(tauc, "tauc.time") <- tauc.time
  }

  if (model.type %in% c("mcp", "mnet", "scad", "snet")) {
    class(tauc) <- c(
      "hdnom.validate.external",
      "ncvreg.validate.external"
    )
    attr(tauc, "model.type") <- model.type
    attr(tauc, "tauc.type") <- tauc.type
    attr(tauc, "tauc.time") <- tauc.time
  }

  if (model.type %in% c("flasso")) {
    class(tauc) <- c(
      "hdnom.validate.external",
      "penalized.validate.external"
    )
    attr(tauc, "model.type") <- model.type
    attr(tauc, "tauc.type") <- tauc.type
    attr(tauc, "tauc.time") <- tauc.time
  }

  tauc
}
