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
