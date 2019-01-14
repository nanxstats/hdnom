#' Survival curve prediction for glmnet objects
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
#' @export glmnet_survcurve
#'
#' @examples
#' NULL
glmnet_survcurve <- function(object, time, event, x, survtime) {
  lp <- as.numeric(predict(
    object,
    newx = data.matrix(x),
    s = object$"lambda", type = "link"
  ))
  basesurv <- glmnet_basesurv(time, event, lp, sort(survtime))
  p <- exp(exp(lp) %*% (-t(basesurv$cumulative_base_hazard)))
  colnames(p) <- names(sort(survtime))

  if (nrow(p) != nrow(x) || ncol(p) != length(survtime)) {
    stop("Prediction error when estimating baseline hazard")
  }

  list("p" = p, "lp" = lp)
}

#' Breslow baseline hazard estimator for glmnet objects
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
#' @export glmnet_basesurv
#'
#' @return list containing cumulative base hazard
#'
#' @examples
#' NULL
glmnet_basesurv <- function(
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

#' Survival curve prediction for ncvreg objects
#'
#' Derived from c060::predictProb.coxnet
#'
#' @param object \code{ncvreg} model object
#' @param time Survival time
#' @param event Status indicator
#' @param x Predictor matrix
#' @param survtime Survival time to evaluate
#'
#' @export ncvreg_survcurve
#'
#' @return list containing predicted survival probabilities and
#' linear predictors for all samples
#'
#' @examples
#' NULL
ncvreg_survcurve <- function(object, time, event, x, survtime) {
  lp <- as.numeric(predict(object, X = data.matrix(x), type = "link"))
  basesurv <- ncvreg_basesurv(time, event, lp, sort(survtime))
  p <- exp(exp(lp) %*% (-t(basesurv$cumulative_base_hazard)))
  colnames(p) <- names(sort(survtime))

  if (nrow(p) != nrow(x) || ncol(p) != length(survtime)) {
    stop("Prediction error when estimating baseline hazard")
  }

  list("p" = p, "lp" = lp)
}

#' Breslow baseline hazard estimator for ncvreg objects
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
#' @export ncvreg_basesurv
#'
#' @return list containing cumulative base hazard
#'
#' @examples
#' NULL
ncvreg_basesurv <- function(
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

#' Survival curve prediction for penfit objects
#'
#' Derived from c060::predictProb.coxnet
#'
#' @param object \code{penalized} model object
#' @param time Survival time
#' @param event Status indicator
#' @param x Predictor matrix
#' @param survtime Survival time to evaluate
#'
#' @export penalized_survcurve
#'
#' @return list containing predicted survival probabilities and
#' linear predictors for all samples
#'
#' @examples
#' NULL
penalized_survcurve <- function(object, time, event, x, survtime) {
  lp <- as.numeric(object@lin.pred)
  basesurv <- penalized_basesurv(time, event, lp, sort(survtime))
  p <- exp(exp(lp) %*% (-t(basesurv$cumulative_base_hazard)))
  colnames(p) <- names(sort(survtime))

  if (nrow(p) != nrow(x) || ncol(p) != length(survtime)) {
    stop("Prediction error when estimating baseline hazard")
  }

  list("p" = p, "lp" = lp)
}

#' Breslow baseline hazard estimator for penfit objects
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
#' @export penalized_basesurv
#'
#' @return list containing cumulative base hazard
#'
#' @examples
#' NULL
penalized_basesurv <- function(
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
