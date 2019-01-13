#' Nomograms for High-Dimensional Cox Models
#'
#' Nomograms for High-Dimensional Cox Models
#'
#' @param object Fitted model object.
#' @param model.type Fitted model type. Could be one of \code{"lasso"},
#' \code{"alasso"}, \code{"flasso"}, \code{"enet"}, \code{"aenet"},
#' \code{"mcp"}, \code{"mnet"}, \code{"scad"}, or \code{"snet"}.
#' @param x Matrix of training data used for fitting the model.
#' @param time Survival time.
#' Must be of the same length with the number of rows as \code{x}.
#' @param event Status indicator, normally 0 = alive, 1 = dead.
#' Must be of the same length with the number of rows as \code{x}.
#' @param pred.at Time point at which to plot nomogram prediction axis.
#' @param fun.at Function values to label on axis.
#' @param funlabel Label for \code{fun} axis.
#'
#' @note We will try to use the value of the selected penalty
#' parameter (e.g. lambda, alpha) in the model object fitted by
#' \code{\link[glmnet]{glmnet}}, \code{\link[ncvreg]{ncvsurv}}, or
#' \code{\link[penalized]{penalized}}. The selected variables under
#' the penalty parameter will be used to build the nomogram
#' and make predictions. This means, if you use routines other
#' than \code{hdcox.*} functions to fit the model, please
#' make sure that there is only one set of necessary parameters
#' (e.g. only a single lambda should be in glmnet objects intead
#' of a vector of multiple lambdas) in the fitted model object.
#'
#' @export hdnom.nomogram
#'
#' @importFrom stats coef as.formula
#'
#' @examples
#' library("hdnom")
#' library("survival")
#'
#' # Load imputed SMART data
#' data(smart)
#' x <- as.matrix(smart[, -c(1, 2)])
#' time <- smart$TEVENT
#' event <- smart$EVENT
#' y <- Surv(time, event)
#'
#' # Fit penalized Cox model with lasso penalty
#' fit <- hdcox.lasso(x, y, nfolds = 5, rule = "lambda.1se", seed = 11)
#'
#' # Generate hdnom.nomogram objects and plot nomogram
#' nom <- hdnom.nomogram(
#'   fit$lasso_model,
#'   model.type = "lasso",
#'   x, time, event, pred.at = 365 * 2,
#'   funlabel = "2-Year Overall Survival Probability"
#' )
#'
#' print(nom)
#' plot(nom)
hdnom.nomogram <- function(
  object,
  model.type =
    c(
      "lasso", "alasso", "flasso", "enet", "aenet",
      "mcp", "mnet", "scad", "snet"
    ),
  x, time, event,
  pred.at = NULL, fun.at = NULL, funlabel = NULL) {

  # input parameter sanity check
  model.type <- match.arg(model.type)

  if (nrow(x) != length(time) || nrow(x) != length(event)) {
    stop("Number of x rows and length of time/event did not match")
  }

  if (is.null(pred.at)) stop("Missing argument pred.at")

  # convert hdcox models to nomogram object
  nomogram.obj <- convert.model(object, x)

  # compute survival curves
  if (model.type %in% c("lasso", "alasso", "enet", "aenet")) {
    if (!all(c("coxnet", "glmnet") %in% class(object))) {
      stop('object class must be "glmnet" and "coxnet"')
    }

    if (length(object$"lambda") != 1L) {
      stop("There should be one and only one lambda in the model object")
    }

    idx_ones <- which(event == 1L)
    survtime_ones <- time[idx_ones]
    names(survtime_ones) <- idx_ones
    survtime_ones <- sort(survtime_ones)
    survtime_at <- survtime_ones[which(survtime_ones > pred.at)[1L] - 1L]
    survtime_at_idx <- names(survtime_at)

    survcurve <- glmnet.survcurve(
      object = object, time = time, event = event,
      x = x, survtime = survtime_ones
    )
  }

  if (model.type %in% c("mcp", "mnet", "scad", "snet")) {
    if (!all(c("ncvsurv", "ncvreg") %in% class(object))) {
      stop('object class must be "ncvreg" and "ncvsurv"')
    }

    idx_ones <- which(event == 1L)
    survtime_ones <- time[idx_ones]
    names(survtime_ones) <- idx_ones
    survtime_ones <- sort(survtime_ones)
    survtime_at <- survtime_ones[which(survtime_ones > pred.at)[1L] - 1L]
    survtime_at_idx <- names(survtime_at)

    survcurve <- ncvreg.survcurve(
      object = object, time = time, event = event,
      x = x, survtime = survtime_ones
    )
  }

  if (model.type %in% c("flasso")) {
    if (!("penfit" %in% class(object))) {
      stop('object class must be "penfit"')
    }

    idx_ones <- which(event == 1L)
    survtime_ones <- time[idx_ones]
    names(survtime_ones) <- idx_ones
    survtime_ones <- sort(survtime_ones)
    survtime_at <- survtime_ones[which(survtime_ones > pred.at)[1L] - 1L]
    survtime_at_idx <- names(survtime_at)

    survcurve <- penalized.survcurve(
      object = object, time = time, event = event,
      x = x, survtime = survtime_ones
    )
  }

  # compute baseline harzard
  baseline <- exp(
    log(survcurve$p[1L, which(colnames(survcurve$p) == survtime_at_idx)]) /
      exp(survcurve$lp[1L])
  )
  bhfun <- function(z) baseline^exp(z)

  # set prediction time points
  if (is.null(fun.at)) {
    fun.at <- c(0.05, 0.2, 0.4, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99)
  }

  if (is.null(funlabel)) {
    funlabel <- paste("Overall Survival Probability at Time", pred.at)
  }

  nom <- list(
    "nomogram.obj" = nomogram.obj,
    "survcurve" = survcurve,
    "bhfun" = bhfun,
    "pred.at" = pred.at,
    "fun.at" = fun.at,
    "funlabel" = funlabel
  )

  class(nom) <- "hdnom.nomogram"

  nom
}

#' Print Nomograms Objects Generated by hdnom.nomogram
#'
#' Print Nomograms Objects Generated by hdnom.nomogram
#'
#' @param x An object returned by \code{\link{hdnom.nomogram}}.
#' @param ... Other parameters.
#'
#' @method print hdnom.nomogram
#'
#' @export
#'
#' @examples
#' NULL
print.hdnom.nomogram <- function(x, ...) {
  if (class(x) != "hdnom.nomogram") {
    stop('object class must be "hdnom.nomogram"')
  }

  nom.raw <- nomogram.raw(
    fit = x$"nomogram.obj", fun = x$"bhfun",
    fun.at = x$"fun.at", funlabel = x$"funlabel",
    lp = TRUE, vnames = "labels", ...
  )

  print(nom.raw)
}

#' Plot Nomogram Objects Generated by hdnom.nomogram
#'
#' Plot Nomogram Objects Generated by hdnom.nomogram
#'
#' @param x An object returned by \code{\link{hdnom.nomogram}}.
#' @param ... Other parameters.
#'
#' @method plot hdnom.nomogram
#'
#' @export
#'
#' @importFrom graphics plot
#'
#' @examples
#' NULL
plot.hdnom.nomogram <- function(x, ...) {
  if (class(x) != "hdnom.nomogram") {
    stop('object class must be "hdnom.nomogram"')
  }

  nom.raw <- nomogram.raw(
    fit = x$"nomogram.obj", fun = x$"bhfun",
    fun.at = x$"fun.at", funlabel = x$"funlabel",
    lp = TRUE, vnames = "labels", ...
  )

  plot(nom.raw)
}
