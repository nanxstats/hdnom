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
#' @param ddist Data frame version of x, made by \code{\link[rms]{datadist}}.
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
#' @importFrom rms ols nomogram
#' @importFrom stats coef as.formula
#'
#' @examples
#' library("rms")
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
#' # Prepare data needed for plotting nomogram
#' x.df <- as.data.frame(x)
#' dd <- datadist(x.df)
#' options(datadist = "dd")
#'
#' # Generate hdnom.nomogram objects and plot nomogram
#' nom <- hdnom.nomogram(
#'   fit$lasso_model,
#'   model.type = "lasso",
#'   x, time, event, x.df, pred.at = 365 * 2,
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
  x, time, event, ddist,
  pred.at = NULL, fun.at = NULL, funlabel = NULL) {
  model.type <- match.arg(model.type)

  if (nrow(x) != length(time) || nrow(x) != length(event)) {
    stop("Number of x rows and length of time/event did not match")
  }

  if (is.null(pred.at)) stop("Missing argument pred.at")

  if (model.type %in% c("lasso", "alasso", "enet", "aenet")) {
    if (!all(c("coxnet", "glmnet") %in% class(object))) {
      stop('object class must be "glmnet" and "coxnet"')
    }

    if (length(object$"lambda") != 1L) {
      stop("There should be one and only one lambda in the model object")
    }

    glmnet_pred_lp <- as.vector(
      predict(object, newx = x, s = object$"lambda", type = "link")
    )

    all_vars <- rownames(object$"beta")
    selected_vars <-
      all_vars[which(as.vector(
        abs(coef(object, s = object$"lambda")) > .Machine$double.eps
      ))]
    ols_formula <- paste(
      "glmnet_pred_lp ~",
      paste(paste("`", selected_vars, "`", sep = ""), collapse = " + ")
    )
    ols_fit <- ols(
      as.formula(ols_formula),
      data = ddist, sigma = 1, x = TRUE, y = TRUE
    )

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

    ncvreg_pred_lp <- predict(object, X = x, type = "link")

    all_vars <- rownames(object$"beta")
    selected_vars <-
      all_vars[which(as.vector(abs(coef(object)) > .Machine$double.eps))]
    ols_formula <- paste(
      "ncvreg_pred_lp ~",
      paste(paste("`", selected_vars, "`", sep = ""), collapse = " + ")
    )
    ols_fit <- ols(
      as.formula(ols_formula),
      data = ddist, sigma = 1, x = TRUE, y = TRUE
    )

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

    penalized_pred_lp <- as.numeric(object@lin.pred)

    all_vars <- colnames(x)
    selected_vars <-
      all_vars[which(abs(object@penalized) > .Machine$double.eps)]
    ols_formula <- paste(
      "penalized_pred_lp ~",
      paste(paste("`", selected_vars, "`", sep = ""), collapse = " + ")
    )
    ols_fit <- ols(
      as.formula(ols_formula),
      data = ddist, sigma = 1, x = TRUE, y = TRUE
    )

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

  baseline <- exp(
    log(survcurve$p[1L, which(colnames(survcurve$p) == survtime_at_idx)]) /
      exp(survcurve$lp[1L])
  )
  bhfun <- function(z) baseline^exp(z)

  if (is.null(fun.at)) {
    fun.at <- c(0.05, 0.2, 0.4, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99)
  }

  if (is.null(funlabel)) {
    funlabel <- paste("Overall Survival Probability at Time", pred.at)
  }

  nom <- list(
    "ols_fit" = ols_fit,
    "survcurve" = survcurve,
    "bhfun" = bhfun,
    "pred.at" = pred.at,
    "fun.at" = fun.at,
    "funlabel" = funlabel
  )

  class(nom) <- "hdnom.nomogram"

  nom
}

#' Survival Curve Prediction for glmnet Objects
#'
#' Derived from c060::predictProb.coxnet
#'
#' @param object \code{glmnet} model object
#' @param time Survival time
#' @param event Status indicator
#' @param x Predictor matrix
#' @param survtime Survival time to evaluate
#'
#' @return list containing predicted survival probabilities and
#' linear predictors for all samples
#'
#' @export glmnet.survcurve
#'
#' @examples
#' NULL
glmnet.survcurve <- function(object, time, event, x, survtime) {
  lp <- as.numeric(predict(
    object,
    newx = data.matrix(x),
    s = object$"lambda", type = "link"
  ))
  basesurv <- glmnet.basesurv(time, event, lp, sort(survtime))
  p <- exp(exp(lp) %*% (-t(basesurv$cumulative_base_hazard)))
  colnames(p) <- names(sort(survtime))

  if (nrow(p) != nrow(x) || ncol(p) != length(survtime)) {
    stop("Prediction error when estimating baseline hazard")
  }

  list("p" = p, "lp" = lp)
}

#' Breslow Baseline Hazard Estimator for glmnet Objects
#'
#' Derived from \code{peperr:::basesurv} and \code{gbm::basehaz.gbm}.
#'
#' @param time Survival time
#' @param event Status indicator
#' @param lp Linear predictors
#' @param times.eval Survival time to evaluate
#' @param centered Should we center the survival curve?
#' See \code{\link[survival]{basehaz}} for details.
#'
#' @importFrom stats approx
#'
#' @export glmnet.basesurv
#'
#' @return list containing cumulative base hazard
#'
#' @examples
#' NULL
glmnet.basesurv <- function(
  time, event, lp,
  times.eval = NULL, centered = FALSE) {
  if (is.null(times.eval)) times.eval <- sort(unique(time))

  t.unique <- sort(unique(time[event == 1L]))
  alpha <- length(t.unique)

  for (i in 1L:length(t.unique)) {
    alpha[i] <- sum(
      time[event == 1L] == t.unique[i]
    ) / sum(exp(lp[time >= t.unique[i]]))
  }

  obj <- approx(t.unique, cumsum(alpha), yleft = 0, xout = times.eval, rule = 2)

  if (centered) obj$y <- obj$y * exp(mean(lp))
  obj$z <- exp(-obj$y)

  names(obj) <- c("times", "cumulative_base_hazard", "base_surv")

  obj
}

#' Survival Curve Prediction for ncvreg Objects
#'
#' Derived from c060::predictProb.coxnet
#'
#' @param object \code{ncvreg} model object
#' @param time Survival time
#' @param event Status indicator
#' @param x Predictor matrix
#' @param survtime Survival time to evaluate
#'
#' @export ncvreg.survcurve
#'
#' @return list containing predicted survival probabilities and
#' linear predictors for all samples
#'
#' @examples
#' NULL
ncvreg.survcurve <- function(object, time, event, x, survtime) {
  lp <- as.numeric(predict(object, X = data.matrix(x), type = "link"))
  basesurv <- ncvreg.basesurv(time, event, lp, sort(survtime))
  p <- exp(exp(lp) %*% (-t(basesurv$cumulative_base_hazard)))
  colnames(p) <- names(sort(survtime))

  if (nrow(p) != nrow(x) || ncol(p) != length(survtime)) {
    stop("Prediction error when estimating baseline hazard")
  }

  list("p" = p, "lp" = lp)
}

#' Breslow Baseline Hazard Estimator for ncvreg Objects
#'
#' Derived from \code{peperr:::basesurv} and \code{gbm::basehaz.gbm}.
#'
#' @param time Survival time
#' @param event Status indicator
#' @param lp Linear predictors
#' @param times.eval Survival time to evaluate
#' @param centered Should we center the survival curve?
#' See \code{\link[survival]{basehaz}} for details.
#'
#' @importFrom stats approx
#'
#' @export ncvreg.basesurv
#'
#' @return list containing cumulative base hazard
#'
#' @examples
#' NULL
ncvreg.basesurv <- function(
  time, event, lp,
  times.eval = NULL, centered = FALSE) {
  if (is.null(times.eval)) times.eval <- sort(unique(time))

  t.unique <- sort(unique(time[event == 1L]))
  alpha <- length(t.unique)

  for (i in 1L:length(t.unique)) {
    alpha[i] <- sum(
      time[event == 1L] == t.unique[i]
    ) / sum(exp(lp[time >= t.unique[i]]))
  }

  obj <- approx(t.unique, cumsum(alpha), yleft = 0, xout = times.eval, rule = 2)

  if (centered) obj$y <- obj$y * exp(mean(lp))
  obj$z <- exp(-obj$y)

  names(obj) <- c("times", "cumulative_base_hazard", "base_surv")

  obj
}

#' Survival Curve Prediction for "penalized" Objects
#'
#' Derived from c060::predictProb.coxnet
#'
#' @param object \code{penalized} model object
#' @param time Survival time
#' @param event Status indicator
#' @param x Predictor matrix
#' @param survtime Survival time to evaluate
#'
#' @export penalized.survcurve
#'
#' @return list containing predicted survival probabilities and
#' linear predictors for all samples
#'
#' @examples
#' NULL
penalized.survcurve <- function(object, time, event, x, survtime) {
  lp <- as.numeric(object@lin.pred)
  basesurv <- penalized.basesurv(time, event, lp, sort(survtime))
  p <- exp(exp(lp) %*% (-t(basesurv$cumulative_base_hazard)))
  colnames(p) <- names(sort(survtime))

  if (nrow(p) != nrow(x) || ncol(p) != length(survtime)) {
    stop("Prediction error when estimating baseline hazard")
  }

  list("p" = p, "lp" = lp)
}

#' Breslow Baseline Hazard Estimator for "penalized" Objects
#'
#' Derived from \code{peperr:::basesurv} and \code{gbm::basehaz.gbm}.
#'
#' @param time Survival time
#' @param event Status indicator
#' @param lp Linear predictors
#' @param times.eval Survival time to evaluate
#' @param centered Should we center the survival curve?
#' See \code{\link[survival]{basehaz}} for details.
#'
#' @importFrom stats approx
#'
#' @export penalized.basesurv
#'
#' @return list containing cumulative base hazard
#'
#' @examples
#' NULL
penalized.basesurv <- function(
  time, event, lp,
  times.eval = NULL, centered = FALSE) {
  if (is.null(times.eval)) times.eval <- sort(unique(time))

  t.unique <- sort(unique(time[event == 1L]))
  alpha <- length(t.unique)

  for (i in 1L:length(t.unique)) {
    alpha[i] <- sum(
      time[event == 1L] == t.unique[i]
    ) / sum(exp(lp[time >= t.unique[i]]))
  }

  obj <- approx(t.unique, cumsum(alpha), yleft = 0, xout = times.eval, rule = 2)

  if (centered) obj$y <- obj$y * exp(mean(lp))
  obj$z <- exp(-obj$y)

  names(obj) <- c("times", "cumulative_base_hazard", "base_surv")

  obj
}

#' Print Nomograms Objects Generated by hdnom.nomogram
#'
#' Print Nomograms Objects Generated by hdnom.nomogram
#'
#' @param x An object returned by \code{\link{hdnom.nomogram}}.
#' @param ... Other parameters for \code{\link[rms]{nomogram}}.
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

  nom <- nomogram(
    fit = x$"ols_fit", fun = x$"bhfun",
    fun.at = x$"fun.at", funlabel = x$"funlabel",
    lp = TRUE, vnames = "labels", ...
  )

  print(nom)
}

#' Plot Nomogram Objects Generated by hdnom.nomogram
#'
#' Plot Nomogram Objects Generated by hdnom.nomogram
#'
#' @param x An object returned by \code{\link{hdnom.nomogram}}.
#' @param ... Other parameters for \code{\link[rms]{nomogram}}.
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

  nom <- nomogram(
    fit = x$"ols_fit", fun = x$"bhfun",
    fun.at = x$"fun.at", funlabel = x$"funlabel",
    lp = TRUE, vnames = "labels", ...
  )

  plot(nom)
}
