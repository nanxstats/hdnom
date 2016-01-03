#' Log-Rank Test for Internal Calibration and External Calibration Results
#'
#' Log-Rank Test for Internal Calibration and External Calibration Results
#'
#' @param object An object returned by \code{\link{hdnom.calibrate}} or
#' \code{\link{hdnom.external.calibrate}}.
#'
#' @importFrom survival survdiff
#' @importFrom survival Surv
#' @importFrom stats formula
#'
#' @export hdnom.logrank
#'
#' @examples
#' library("survival")
#'
#' # Load imputed SMART data
#' data("smart")
#'
#' # Use the first 1000 samples as training data
#' # (the data used for internal validation)
#' x = as.matrix(smart[, -c(1, 2)])[1:1000, ]
#' time = smart$TEVENT[1:1000]
#' event = smart$EVENT[1:1000]
#'
#' # Take the next 1000 samples as external calibration data
#' # In practice, usually use data collected in other studies
#' x_new = as.matrix(smart[, -c(1, 2)])[1001:2000, ]
#' time_new = smart$TEVENT[1001:2000]
#' event_new = smart$EVENT[1001:2000]
#'
#' # Fit Cox model by lasso penalization
#' lassofit = hdcox.lasso(x, Surv(time, event), nfolds = 5, rule = "lambda.1se", seed = 11)
#'
#' ### Internal calibration
#' cal.int = hdnom.calibrate(x, time, event, model.type = "lasso",
#'                           alpha = 1, lambda = lassofit$'lasso_best_lambda',
#'                           method = "cv", nfolds = 5,
#'                           pred.at = 365 * 9, ngroup = 3)
#'
#' hdnom.logrank(cal.int)
#'
#' ### External calibration
#' cal.ext =
#'   hdnom.external.calibrate(lassofit, x, time, event,
#'                            x_new, time_new, event_new,
#'                            pred.at = 365 * 5, ngroup = 3)
#'
#' hdnom.logrank(cal.ext)
hdnom.logrank = function(object) {

  if (!(any(c('hdnom.calibrate', 'hdnom.external.calibrate') %in% class(object))))
    stop('object class must be "hdnom.calibrate" or "hdnom.external.calibrate"')

  time = attr(object, 'surv.time')
  event = attr(object, 'surv.event')
  grp = attr(object, 'risk.group')

  sdiff = survdiff(formula('Surv(time, event) ~ grp'))
  pval = pchisq(sdiff$'chisq', length(sdiff$'n') - 1L, lower.tail = FALSE)
  sdiff$'pval' = pval

  return(sdiff)

}
