#' Externally Calibrate High-Dimensional Cox Models
#'
#' Externally Calibrate High-Dimensional Cox Models
#'
#' @param object Model object fitted by \code{hdcox.*()} functions.
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
#' @param pred.at Time point at which external calibration should take place.
#' @param ngroup Number of groups to be formed for external calibration.
#'
#' @importFrom survival Surv
#' @importFrom stats quantile
#' @importFrom stats median
#'
#' @export hdnom.external.calibrate
#'
#' @examples
#' library("survival")
#'
#' # Load imputed SMART data
#' data(smart)
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
#' # External calibration
#' cal.ext =
#'   hdnom.external.calibrate(lassofit, x, time, event,
#'                            x_new, time_new, event_new,
#'                            pred.at = 365 * 5, ngroup = 5)
#'
#' print(cal.ext)
#' summary(cal.ext)
#' plot(cal.ext, xlim = c(0.6, 1), ylim = c(0.6, 1))
#'
# ### Testing fused lasso, MCP, and Snet models ###
# library("survival")
#
# # Load imputed SMART data
# data(smart)
# # Use first 500 samples as training data
# # (the data used for internal validation)
# x = as.matrix(smart[, -c(1, 2)])[1:500, ]
# time = smart$TEVENT[1:500]
# event = smart$EVENT[1:500]
#
# # Take 1000 samples as external validation data.
# # In practice, usually use data collected in other studies.
# x_new = as.matrix(smart[, -c(1, 2)])[1001:2000, ]
# time_new = smart$TEVENT[1001:2000]
# event_new = smart$EVENT[1001:2000]
#
# flassofit = hdcox.flasso(x, Surv(time, event), nfolds = 5, seed = 11)
# scadfit = hdcox.mcp(x, Surv(time, event), nfolds = 5, seed = 11)
# mnetfit = hdcox.snet(x, Surv(time, event), nfolds = 5, seed = 11)
#
# cal.ext1 =
#   hdnom.external.calibrate(flassofit, x, time, event,
#                            x_new, time_new, event_new,
#                            pred.at = 365 * 5, ngroup = 5)
#
# cal.ext2 =
#   hdnom.external.calibrate(scadfit, x, time, event,
#                            x_new, time_new, event_new,
#                            pred.at = 365 * 5, ngroup = 5)
#
# cal.ext3 =
#   hdnom.external.calibrate(mnetfit, x, time, event,
#                            x_new, time_new, event_new,
#                            pred.at = 365 * 5, ngroup = 5)
#
# print(cal.ext1)
# summary(cal.ext1)
# plot(cal.ext1)
#
# print(cal.ext2)
# summary(cal.ext2)
# plot(cal.ext2)
#
# print(cal.ext3)
# summary(cal.ext3)
# plot(cal.ext3)
hdnom.external.calibrate = function(object, x, time, event,
                                    x_new, time_new, event_new,
                                    pred.at, ngroup = 5) {

  if (!('hdcox.model' %in% class(object)))
    stop('object must be of class "hdcox.model" fitted by hdcox.* functions')

  model.type = gsub('hdcox.model.', '', setdiff(class(object), 'hdcox.model'))
  model_object = object[[paste0(model.type, '_model')]]

  if (!('matrix' %in% class(x_new))) stop('x_new must be a matrix')

  if (length(pred.at) != 1L) stop('pred.at should only contain 1 time point')

  if (model.type %in% c('lasso', 'alasso', 'enet', 'aenet')) {
    pred_list = glmnet.external.calibrate.internal.pred(
      object = model_object, x_tr = x, x_te = x_new, y_tr = Surv(time, event),
      pred.at = pred.at
    )
  }

  if (model.type %in% c('mcp', 'mnet', 'scad', 'snet')) {
    pred_list = ncvreg.external.calibrate.internal.pred(
      object = model_object, x_tr = x, x_te = x_new, y_tr = Surv(time, event),
      pred.at = pred.at
    )
  }

  if (model.type %in% c('flasso')) {
    pred_list = penalized.external.calibrate.internal.pred(
      object = model_object, x_tr = x, x_te = x_new, y_tr = Surv(time, event),
      pred.at = pred.at
    )
  }

  pred_prob = rep(NA, nrow(x_new))
  for (i in 1L:length(pred_prob)) pred_prob[i] = pred_list$p[i, pred_list$idx]
  grp = cut(pred_prob, quantile(pred_prob, seq(0, 1, 1/ngroup)), labels = 1L:ngroup)

  pred_prob_median = tapply(pred_prob, grp, median)

  true_prob = hdnom.external.calibrate.internal.true(
    pred_prob, grp, time_new, event_new, pred.at, ngroup
  )

  prob = cbind(pred_prob_median, true_prob)
  colnames(prob)[1L] = 'Predicted'

  if (model.type %in% c('lasso', 'alasso', 'enet', 'aenet'))
    class(prob) = c('hdnom.external.calibrate', 'glmnet.external.calibrate')

  if (model.type %in% c('mcp', 'mnet', 'scad', 'snet'))
    class(prob) = c('hdnom.external.calibrate', 'ncvreg.external.calibrate')

  if (model.type %in% c('flasso'))
    class(prob) = c('hdnom.external.calibrate', 'penalized.external.calibrate')

  attr(prob, 'model.type') = model.type
  attr(prob, 'pred.at')    = pred.at
  attr(prob, 'ngroup')     = ngroup
  attr(prob, 'risk.group') = grp
  attr(prob, 'surv.time') = time_new
  attr(prob, 'surv.event') = event_new

  prob

}

#' Compute glmnet Predicted Survival Probabilities for External Calibration
#'
#' @importFrom glmnet glmnet
#' @importFrom stats predict
#'
#' @return list containing predicted survival probability
#'
#' @keywords internal
glmnet.external.calibrate.internal.pred = function(object, x_tr, x_te, y_tr,
                                                   alpha, lambda, pen.factor,
                                                   pred.at) {

  lp = as.numeric(predict(object, newx = data.matrix(x_tr), type = 'link'))
  lpnew = as.numeric(predict(object, newx = data.matrix(x_te), type = 'link'))

  time_tr = y_tr[, 1L]
  event_tr = y_tr[, 2L]
  idx_ones = which(event_tr == 1L)
  if (length(idx_ones) == 0L)
    stop('No 1 events in the training fold, please try other random seeds')
  survtime_ones = time_tr[idx_ones]
  names(survtime_ones) = idx_ones
  survtime_ones = sort(survtime_ones)

  basesurv = glmnet.basesurv(time_tr, event_tr, lp, survtime_ones)
  p = exp(exp(lpnew) %*% (-t(basesurv$cumulative_base_hazard)))

  if (nrow(p) != nrow(x_te) || ncol(p) != length(survtime_ones))
    stop('Prediction error when estimating baseline hazard')

  idx = length(which(survtime_ones <= pred.at))

  list('p' = p, 'idx' = idx)

}

#' Compute ncvreg Predicted Survival Probabilities for External Calibration
#'
#' @importFrom ncvreg ncvsurv
#' @importFrom stats predict
#'
#' @return list containing predicted survival probability
#'
#' @keywords internal
ncvreg.external.calibrate.internal.pred = function(object, x_tr, x_te, y_tr,
                                                   pred.at) {

  lp = as.numeric(predict(object, X = data.matrix(x_tr), type = 'link'))
  lpnew = as.numeric(predict(object, X = data.matrix(x_te), type = 'link'))

  time_tr = y_tr[, 1L]
  event_tr = y_tr[, 2L]
  idx_ones = which(event_tr == 1L)
  if (length(idx_ones) == 0L)
    stop('No 1 events in the training fold, please try other random seeds')
  survtime_ones = time_tr[idx_ones]
  names(survtime_ones) = idx_ones
  survtime_ones = sort(survtime_ones)

  basesurv = ncvreg.basesurv(time_tr, event_tr, lp, survtime_ones)
  p = exp(exp(lpnew) %*% (-t(basesurv$cumulative_base_hazard)))

  if (nrow(p) != nrow(x_te) || ncol(p) != length(survtime_ones))
    stop('Prediction error when estimating baseline hazard')

  idx = length(which(survtime_ones <= pred.at))

  list('p' = p, 'idx' = idx)

}

#' Compute "penalized" Predicted Survival Probabilities for External Calibration
#'
#' @importFrom penalized penalized
#' @importFrom stats predict
#'
#' @return list containing predicted survival probability
#'
#' @keywords internal
penalized.external.calibrate.internal.pred = function(object, x_tr, x_te, y_tr,
                                                      pred.at) {

  lp = as.vector(data.matrix(x_tr) %*% as.matrix(object@'penalized'))
  lpnew = as.vector(data.matrix(x_te) %*% as.matrix(object@'penalized'))

  time_tr = y_tr[, 1L]
  event_tr = y_tr[, 2L]
  idx_ones = which(event_tr == 1L)
  if (length(idx_ones) == 0L)
    stop('No 1 events in the training fold, please try other random seeds')
  survtime_ones = time_tr[idx_ones]
  names(survtime_ones) = idx_ones
  survtime_ones = sort(survtime_ones)

  basesurv = penalized.basesurv(time_tr, event_tr, lp, survtime_ones)
  p = exp(exp(lpnew) %*% (-t(basesurv$cumulative_base_hazard)))

  if (nrow(p) != nrow(x_te) || ncol(p) != length(survtime_ones))
    stop('Prediction error when estimating baseline hazard')

  idx = length(which(survtime_ones <= pred.at))

  list('p' = p, 'idx' = idx)

}

#' Compute Kaplan-Meier Estimated Survival Probabilities for External Calibration
#'
#' @importFrom survival survfit
#' @importFrom survival Surv
#'
#' @return list
#'
#' @keywords internal
hdnom.external.calibrate.internal.true = function(pred_prob, grp,
                                                  time_new, event_new,
                                                  pred.at, ngroup) {

  true_prob = matrix(NA, ncol = 3L, nrow = ngroup)
  colnames(true_prob) = c("Observed", "Lower 95%", "Upper 95%")

  for (i in 1L:ngroup) {
    time_grp = time_new[which(grp == i)]
    event_grp = event_new[which(grp == i)]
    km = survfit(Surv(time_grp, event_grp) ~ 1, type = 'kaplan-meier')
    idx = which(km$time > pred.at)[1L] - 1L
    km_pred_at = km$surv[idx]
    ll_pred_at = km$lower[idx]
    ul_pred_at = km$upper[idx]
    true_prob[i, ] = c(km_pred_at, ll_pred_at, ul_pred_at)
  }

  return(true_prob)

}

#' Print External Calibration Results
#'
#' Print External Calibration Results
#'
#' @param x an object returned by \code{\link{hdnom.external.calibrate}}.
#' @param ... other parameters (not used).
#'
#' @method print hdnom.external.calibrate
#'
#' @export
#'
#' @examples
#' NULL
print.hdnom.external.calibrate = function(x, ...) {

  if (!('hdnom.external.calibrate' %in% class(x)))
    stop('object class must be "hdnom.external.calibrate"')

  method = setdiff(class(x), 'hdnom.external.calibrate')

  cat('High-Dimensional Cox Model External Calibration Object\n')
  cat('Model type:', attr(x, 'model.type'), '\n')
  cat('Calibration time point:', attr(x, 'pred.at'), '\n')
  cat('Number of groups formed for calibration:', attr(x, 'ngroup'), '\n')

}

#' Summary of External Calibration Results
#'
#' Summary of External Calibration Results
#'
#' @param object an object returned by \code{\link{hdnom.external.calibrate}}.
#' @param ... other parameters (not used).
#'
#' @method summary hdnom.external.calibrate
#'
#' @export
#'
#' @examples
#' NULL
summary.hdnom.external.calibrate = function(object, ...) {

  if (!('hdnom.external.calibrate' %in% class(object)))
    stop('object class must be "hdnom.external.calibrate"')

  method = setdiff(class(object), 'hdnom.external.calibrate')

  attr(object, 'model.type') = NULL
  attr(object, 'pred.at')    = NULL
  attr(object, 'ngroup')     = NULL
  attr(object, 'risk.group') = NULL
  attr(object, 'surv.time')  = NULL
  attr(object, 'surv.event') = NULL

  cat('  External Calibration Summary Table\n')
  class(object) = 'matrix'
  print(object)

}

#' Plot External Calibration Results
#'
#' Plot External Calibration Results
#'
#' @param x an object returned by \code{\link{hdnom.external.calibrate}}.
#' @param xlim x axis limits of the plot.
#' @param ylim y axis limits of the plot.
#' @param col.pal Color palette to use. Possible values are
#' \code{"JCO"}, \code{"Lancet"}, \code{"NPG"}, and \code{"AAAS"}.
#' Default is \code{"JCO"}.
#' @param ... other parameters for \code{plot}.
#'
#' @method plot hdnom.external.calibrate
#'
#' @export
#'
#' @importFrom ggplot2 ggplot aes_string geom_errorbar
#' geom_line geom_point geom_abline xlab ylab theme_bw
#'
#' @examples
#' NULL
plot.hdnom.external.calibrate =
  function(x, xlim = c(0, 1), ylim = c(0, 1),
           col.pal = c('JCO', 'Lancet', 'NPG', 'AAAS'), ...) {

    if (!('hdnom.external.calibrate' %in% class(x)))
      stop('object class must be "hdnom.external.calibrate"')

    df = data.frame('pre' = x[, 'Predicted'], 'obs' = x[, 'Observed'],
                    'll' = x[, 'Lower 95%'], 'ul' = x[, 'Upper 95%'])

    col.pal = match.arg(col.pal)
    col_pal = switch (
      col.pal,
      JCO   = palette.jco()[1], Lancet = palette.lancet()[1],
      NPG   = palette.npg()[1], AAAS   = palette.aaas()[1])

    ggplot(df, aes_string(x = 'pre', y = 'obs',
                          xmin = xlim[1L], xmax = xlim[2L],
                          ymin = ylim[1L], ymax = ylim[2L])) +
      geom_abline(slope = 1, intercept = 0, colour = 'grey') +
      geom_errorbar(aes_string(ymin = 'll', ymax = 'ul'), colour = col_pal) +
      geom_line(colour = col_pal) +
      geom_point(size = 3, colour = col_pal) +
      xlab('Predicted Survival Probability') +
      ylab('Observed Survival Probability') +
      theme_bw()

  }
