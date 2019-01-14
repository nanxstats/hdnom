#' Log-rank test for internal calibration and external calibration results
#'
#' Log-rank test for internal calibration and external calibration results
#'
#' @param object An object returned by \code{\link{calibrate}} or
#' \code{\link{calibrate_external}}.
#'
#' @importFrom survival survdiff
#' @importFrom survival Surv
#' @importFrom stats formula
#'
#' @export logrank_test
#'
#' @examples
#' library("survival")
#'
#' # Load imputed SMART data
#' data("smart")
#'
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
#'   nfolds = 5, rule = "lambda.1se", seed = 11
#' )
#'
#' # Internal calibration
#' cal.int <- calibrate(
#'   x, time, event,
#'   model.type = "lasso",
#'   alpha = 1, lambda = fit$lasso_best_lambda,
#'   method = "cv", nfolds = 5,
#'   pred.at = 365 * 9, ngroup = 3
#' )
#'
#' logrank_test(cal.int)
#'
#' # External calibration
#' cal.ext <- calibrate_external(
#'   fit, x, time, event,
#'   x_new, time_new, event_new,
#'   pred.at = 365 * 5, ngroup = 3
#' )
#'
#' logrank_test(cal.ext)
logrank_test <- function(object) {
  if (!(any(c("hdnom.calibrate", "hdnom.calibrate.external") %in% class(object)))) {
    stop('object class must be "hdnom.calibrate" or "hdnom.calibrate.external"')
  }

  time <- attr(object, "surv.time")
  event <- attr(object, "surv.event")
  grp <- attr(object, "risk.group")

  sdiff <- survdiff(formula("Surv(time, event) ~ grp"))
  pval <- pchisq(sdiff$"chisq", length(sdiff$"n") - 1L, lower.tail = FALSE)
  sdiff$"pval" <- pval

  sdiff
}
