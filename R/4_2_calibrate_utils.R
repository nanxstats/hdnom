#' Compute glmnet predicted survival probabilities for calibration
#'
#' @importFrom glmnet glmnet
#' @importFrom stats predict
#'
#' @return list containing predicted survival probability
#'
#' @keywords internal
glmnet_calibrate_surv_prob_pred <- function(
  x_tr, x_te, y_tr,
  alpha, lambda, pen.factor,
  pred.at) {
  if (is.null(pen.factor)) {
    object <- glmnet(
      x = x_tr, y = y_tr, family = "cox",
      alpha = alpha, lambda = lambda
    )
  } else {
    object <- glmnet(
      x = x_tr, y = y_tr, family = "cox",
      alpha = alpha, lambda = lambda,
      penalty.factor = pen.factor
    )
  }

  lp <- as.numeric(
    predict(object, newx = data.matrix(x_tr), s = lambda, type = "link")
  )
  lpnew <- as.numeric(
    predict(object, newx = data.matrix(x_te), s = lambda, type = "link")
  )

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

#' Compute ncvreg predicted survival probabilities for calibration
#'
#' @importFrom ncvreg ncvsurv
#' @importFrom stats predict
#'
#' @return list containing predicted survival probability
#'
#' @keywords internal
ncvreg_calibrate_surv_prob_pred <- function(
  x_tr, x_te, y_tr,
  model.type,
  alpha, lambda, gamma,
  pred.at) {
  if (model.type == "mcp") {
    object <- ncvreg::ncvsurv(
      X = x_tr, y = y_tr,
      penalty = "MCP", gamma = gamma,
      alpha = 1, lambda = lambda
    )
  }

  if (model.type == "mnet") {
    object <- ncvreg::ncvsurv(
      X = x_tr, y = y_tr,
      penalty = "MCP", gamma = gamma,
      alpha = alpha, lambda = lambda
    )
  }

  if (model.type == "scad") {
    object <- ncvreg::ncvsurv(
      X = x_tr, y = y_tr,
      penalty = "SCAD", gamma = gamma,
      alpha = 1, lambda = lambda
    )
  }

  if (model.type == "snet") {
    object <- ncvreg::ncvsurv(
      X = x_tr, y = y_tr,
      penalty = "SCAD", gamma = gamma,
      alpha = alpha, lambda = lambda
    )
  }

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

#' Compute penfit predicted survival probabilities for calibration
#'
#' @importFrom penalized penalized
#' @importFrom stats predict
#'
#' @return list containing predicted survival probability
#'
#' @keywords internal
penalized_calibrate_surv_prob_pred <- function(
  x_tr, x_te, y_tr,
  lambda1, lambda2,
  pred.at) {
  object <- penalized(
    response = y_tr, penalized = x_tr,
    lambda1 = lambda1, lambda2 = lambda2,
    maxiter = 25, epsilon = 1e-3, # for faster convergence, consistent with `fit_flasso()`
    fusedl = TRUE, standardize = TRUE, model = "cox"
  )

  lp <- as.vector(data.matrix(x_tr) %*% as.matrix(object@penalized))
  lpnew <- as.vector(data.matrix(x_te) %*% as.matrix(object@penalized))

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

#' Compute Kaplan-Meier estimated survival probabilities for calibration
#'
#' @importFrom survival survfit
#' @importFrom survival Surv
#'
#' @return list
#'
#' @keywords internal
calibrate_surv_prob_true <- function(
  pred_prob, grp,
  time, event,
  pred.at, ngroup) {
  true_prob <- matrix(NA, ncol = 3L, nrow = ngroup)
  colnames(true_prob) <- c("Observed", "Lower 95%", "Upper 95%")

  for (i in 1L:ngroup) {
    time_grp <- time[which(grp == i)]
    event_grp <- event[which(grp == i)]
    km <- survfit(Surv(time_grp, event_grp) ~ 1, type = "kaplan-meier")
    idx <- which(km$time > pred.at)[1L] - 1L
    km_pred_at <- km$surv[idx]
    ll_pred_at <- km$lower[idx]
    ul_pred_at <- km$upper[idx]
    true_prob[i, ] <- c(km_pred_at, ll_pred_at, ul_pred_at)
  }

  true_prob
}
