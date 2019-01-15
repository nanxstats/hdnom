#' Construct nomogram ojects for high-dimensional Cox models
#'
#' Construct nomograms ojects for high-dimensional Cox models
#'
#' @param object Model object fitted by `hdnom::fit_*()` functions.
#' @param x Matrix of training data used for fitting the model.
#' @param time Survival time. Must be of the same length with
#' the number of rows as \code{x}.
#' @param event Status indicator, normally 0 = alive, 1 = dead.
#' Must be of the same length with the number of rows as \code{x}.
#' @param pred.at Time point at which to plot nomogram prediction axis.
#' @param fun.at Function values to label on axis.
#' @param funlabel Label for \code{fun} axis.
#'
#' @note The nomogram visualizes the model under the automatically
#' selected "optimal" hyperparameters (e.g. lambda, alpha, gamma).
#'
#' @export as_nomogram
#'
#' @importFrom stats coef as.formula
#'
#' @examples
#' data(smart)
#' x <- as.matrix(smart[, -c(1, 2)])
#' time <- smart$TEVENT
#' event <- smart$EVENT
#' y <- survival::Surv(time, event)
#'
#' fit <- fit_lasso(x, y, nfolds = 5, rule = "lambda.1se", seed = 11)
#'
#' nom <- as_nomogram(
#'   fit, x, time, event, pred.at = 365 * 2,
#'   funlabel = "2-Year Overall Survival Probability"
#' )
#'
#' print(nom)
#' plot(nom)
as_nomogram <- function(object, x, time, event, pred.at = NULL, fun.at = NULL, funlabel = NULL) {
  model <- object$model
  model_type <- object$type

  if (nrow(x) != length(time) || nrow(x) != length(event)) {
    stop("Number of x rows and length of time/event did not match")
  }

  if (is.null(pred.at)) stop("Missing argument pred.at")

  # Convert hdnom model to nomogram object
  nomogram_object <- convert_model(model, x)

  # Compute survival curves
  if (model_type %in% c("lasso", "alasso", "enet", "aenet")) {
    if (!all(c("coxnet", "glmnet") %in% class(model))) {
      stop('model class must be "glmnet" and "coxnet"')
    }

    if (length(model$"lambda") != 1L) {
      stop("There should be one and only one lambda in the model")
    }

    idx_ones <- which(event == 1L)
    survtime_ones <- time[idx_ones]
    names(survtime_ones) <- idx_ones
    survtime_ones <- sort(survtime_ones)
    survtime_at <- survtime_ones[which(survtime_ones > pred.at)[1L] - 1L]
    survtime_at_idx <- names(survtime_at)

    survcurve <- glmnet_survcurve(model, time, event, x, survtime_ones)
  }

  if (model_type %in% c("mcp", "mnet", "scad", "snet")) {
    if (!all(c("ncvsurv", "ncvreg") %in% class(model))) {
      stop('model class must be "ncvreg" and "ncvsurv"')
    }

    idx_ones <- which(event == 1L)
    survtime_ones <- time[idx_ones]
    names(survtime_ones) <- idx_ones
    survtime_ones <- sort(survtime_ones)
    survtime_at <- survtime_ones[which(survtime_ones > pred.at)[1L] - 1L]
    survtime_at_idx <- names(survtime_at)

    survcurve <- ncvreg_survcurve(model, time, event, x, survtime_ones)
  }

  if (model_type %in% c("flasso")) {
    if (!("penfit" %in% class(model))) {
      stop('model class must be "penfit"')
    }

    idx_ones <- which(event == 1L)
    survtime_ones <- time[idx_ones]
    names(survtime_ones) <- idx_ones
    survtime_ones <- sort(survtime_ones)
    survtime_at <- survtime_ones[which(survtime_ones > pred.at)[1L] - 1L]
    survtime_at_idx <- names(survtime_at)

    survcurve <- penalized_survcurve(model, time, event, x, survtime_ones)
  }

  # Compute baseline harzard
  baseline <- exp(
    log(survcurve$p[1L, which(colnames(survcurve$p) == survtime_at_idx)]) /
      exp(survcurve$lp[1L])
  )
  bhfun <- function(z) baseline^exp(z)

  # Set prediction time points
  if (is.null(fun.at)) fun.at <- c(0.05, 0.2, 0.4, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99)
  if (is.null(funlabel)) funlabel <- paste("Overall Survival Probability at Time", pred.at)

  nom <- list(
    "nomogram" = nomogram_object,
    "survcurve" = survcurve,
    "bhfun" = bhfun,
    "pred.at" = pred.at,
    "fun.at" = fun.at,
    "funlabel" = funlabel
  )

  class(nom) <- "hdnom.nomogram"
  nom
}
