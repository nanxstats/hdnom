#' Compute glmnet predicted survival probabilities for external calibration
#'
#' @importFrom glmnet glmnet
#' @importFrom stats predict
#'
#' @return list containing predicted survival probability
#'
#' @keywords internal
glmnet_calibrate_external_surv_prob_pred <- function(object, x_tr, x_te, y_tr, pred.at) {
  lp <- as.numeric(predict(object, newx = data.matrix(x_tr), type = "link"))
  lpnew <- as.numeric(predict(object, newx = data.matrix(x_te), type = "link"))

  time_tr <- y_tr[, 1L]
  event_tr <- y_tr[, 2L]
  idx_ones <- which(event_tr == 1L)
  if (length(idx_ones) == 0L) {
    stop("No 1 events in the training fold, please try other random seeds")
  }
  survtime_ones <- time_tr[idx_ones]
  names(survtime_ones) <- idx_ones
  survtime_ones <- sort(survtime_ones)

  basesurv <- glmnet_basesurv(time_tr, event_tr, lp, survtime_ones)
  p <- exp(exp(lpnew) %*% (-t(basesurv$cumulative_base_hazard)))

  if (nrow(p) != nrow(x_te) || ncol(p) != length(survtime_ones)) {
    stop("Prediction error when estimating baseline hazard")
  }

  idx <- length(which(survtime_ones <= pred.at))

  list("p" = p, "idx" = idx)
}

#' Compute ncvreg predicted survival probabilities for external calibration
#'
#' @importFrom stats predict
#'
#' @return list containing predicted survival probability
#'
#' @keywords internal
ncvreg_calibrate_external_surv_prob_pred <- function(object, x_tr, x_te, y_tr, pred.at) {
  lp <- as.numeric(predict(object, X = data.matrix(x_tr), type = "link"))
  lpnew <- as.numeric(predict(object, X = data.matrix(x_te), type = "link"))

  time_tr <- y_tr[, 1L]
  event_tr <- y_tr[, 2L]
  idx_ones <- which(event_tr == 1L)
  if (length(idx_ones) == 0L) {
    stop("No 1 events in the training fold, please try other random seeds")
  }
  survtime_ones <- time_tr[idx_ones]
  names(survtime_ones) <- idx_ones
  survtime_ones <- sort(survtime_ones)

  basesurv <- ncvreg_basesurv(time_tr, event_tr, lp, survtime_ones)
  p <- exp(exp(lpnew) %*% (-t(basesurv$cumulative_base_hazard)))

  if (nrow(p) != nrow(x_te) || ncol(p) != length(survtime_ones)) {
    stop("Prediction error when estimating baseline hazard")
  }

  idx <- length(which(survtime_ones <= pred.at))

  list("p" = p, "idx" = idx)
}

#' Compute penfit predicted survival probabilities for external calibration
#'
#' @importFrom penalized penalized
#' @importFrom stats predict
#'
#' @return list containing predicted survival probability
#'
#' @keywords internal
penalized_calibrate_external_surv_prob_pred <- function(object, x_tr, x_te, y_tr, pred.at) {
  lp <- as.vector(data.matrix(x_tr) %*% as.matrix(object@"penalized"))
  lpnew <- as.vector(data.matrix(x_te) %*% as.matrix(object@"penalized"))

  time_tr <- y_tr[, 1L]
  event_tr <- y_tr[, 2L]
  idx_ones <- which(event_tr == 1L)
  if (length(idx_ones) == 0L) {
    stop("No 1 events in the training fold, please try other random seeds")
  }
  survtime_ones <- time_tr[idx_ones]
  names(survtime_ones) <- idx_ones
  survtime_ones <- sort(survtime_ones)

  basesurv <- penalized_basesurv(time_tr, event_tr, lp, survtime_ones)
  p <- exp(exp(lpnew) %*% (-t(basesurv$cumulative_base_hazard)))

  if (nrow(p) != nrow(x_te) || ncol(p) != length(survtime_ones)) {
    stop("Prediction error when estimating baseline hazard")
  }

  idx <- length(which(survtime_ones <= pred.at))

  list("p" = p, "idx" = idx)
}

#' Compute Kaplan-Meier estimated survival probabilities for external calibration
#'
#' @importFrom survival survfit
#' @importFrom survival Surv
#'
#' @return list
#'
#' @keywords internal
calibrate_external_surv_prob_true <- function(
    pred_prob, grp,
    time_new, event_new,
    pred.at, ngroup) {
  true_prob <- matrix(NA, ncol = 3L, nrow = ngroup)
  colnames(true_prob) <- c("Observed", "Lower 95%", "Upper 95%")

  for (i in 1L:ngroup) {
    time_grp <- time_new[which(grp == i)]
    event_grp <- event_new[which(grp == i)]
    km <- survfit(Surv(time_grp, event_grp) ~ 1, type = "kaplan-meier")
    idx <- which(km$time > pred.at)[1L] - 1L
    km_pred_at <- km$surv[idx]
    ll_pred_at <- km$lower[idx]
    ul_pred_at <- km$upper[idx]
    true_prob[i, ] <- c(km_pred_at, ll_pred_at, ul_pred_at)
  }

  true_prob
}
