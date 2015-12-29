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
#' @param lambda Value of the penalty parameter lambda in
#' \code{\link[glmnet]{glmnet}} or \code{\link[ncvreg]{ncvsurv}}.
#' Required except when \code{model.type == "flasso"}.
#' We will use the selected variables at the provided \code{lambda} to
#' build the nomogram, and make predictions.
#' See the example for choosing a proper lambda value extracted
#' from cross-validation results.
#' @param pred.at Time point at which to plot nomogram prediction axis.
#' @param fun.at Function values to label on axis.
#' @param funlabel Label for \code{fun} axis.
#'
#' @export hdnom.nomogram
#'
#' @importFrom rms ols nomogram
#' @importFrom stats coef as.formula
#'
#' @examples
#' library("glmnet")
#' library("survival")
#' library("rms")
#'
#' # Load imputed SMART data
#' data(smart)
#' x = as.matrix(smart[, -c(1, 2)])
#' time = smart$TEVENT
#' event = smart$EVENT
#' x.df = as.data.frame(x)
#' dd = datadist(x.df)
#' options(datadist = "dd")
#'
#' # Fit penalized Cox model (lasso penalty) with glmnet
#' set.seed(1010)
#' cvfit = cv.glmnet(x, Surv(time, event), family = "cox", nfolds = 10)
#' fit = glmnet(x, Surv(time, event), family = "cox")
#'
#' # Generate hdnom.nomogram objects and plot nomogram
#' nom = hdnom.nomogram(fit, model.type = "lasso", x, time, event, x.df,
#'                      lambda = cvfit$lambda.1se, pred.at = 365 * 2,
#'                      funlabel = "2-Year Overall Survival Probability")
#'
#' print(nom)
#' plot(nom)
#'
# ### Testing fused lasso models ###
# library("penalized")
# library("survival")
# library("rms")
#
# # Load imputed SMART data
# data(smart)
# x = as.matrix(smart[, -c(1, 2)])
# x = x[1:500, ]
# time = smart$TEVENT[1:500]
# event = smart$EVENT[1:500]
# x.df = as.data.frame(x)
# dd = datadist(x.df)
# options(datadist = "dd")
#
# # Fit penalized Cox model (fused lasso penalty) with penalized package
# set.seed(1010)
# cvfit = optL1(response = Surv(time, event), penalized = x, fusedl = TRUE,
#               standardize = TRUE, model = 'cox', fold = 3, trace = TRUE)
#
# # Generate hdnom.nomogram objects and plot nomogram
# nom = hdnom.nomogram(cvfit$fullfit, model.type = "flasso",
#                      x, time, event, x.df, pred.at = 365 * 2,
#                      funlabel = "2-Year Overall Survival Probability")
#
# print(nom)
# plot(nom)
#
# ### Testing SCAD models ###
# library("ncvreg")
#
# cvfit = cv.ncvsurv(x, Surv(time, event),
#                    model = 'cox', penalty = 'SCAD',
#                    alpha = 1, nfolds = 5,
#                    seed = 1010, trace = TRUE)
# fit = ncvsurv(x, Surv(time, event), model = 'cox', penalty = 'SCAD', alpha = 1)
#
# nom = hdnom.nomogram(fit, model.type = "scad", x, time, event, x.df,
#                      lambda = cvfit$lambda.min, pred.at = 365 * 2,
#                      funlabel = "2-Year Overall Survival Probability")
#
# print(nom)
# plot(nom)
#
# ### Testing Mnet models ###
# cvfit = cv.ncvsurv(x, Surv(time, event),
#                    model = 'cox', penalty = 'MCP',
#                    alpha = 0.5, nfolds = 5,
#                    seed = 1010, trace = TRUE)
# fit = ncvsurv(x, Surv(time, event), model = 'cox', penalty = 'MCP', alpha = 0.5)
#
# # Generate hdnom.nomogram objects and plot nomogram
# nom = hdnom.nomogram(fit, model.type = "mnet", x, time, event, x.df,
#                      lambda = cvfit$lambda.min, pred.at = 365 * 2,
#                      funlabel = "2-Year Overall Survival Probability")
#
# print(nom)
# plot(nom)
hdnom.nomogram = function(object,
                          model.type = c('lasso', 'alasso', 'flasso',
                                         'enet', 'aenet',
                                         'mcp', 'mnet',
                                         'scad', 'snet'),
                          x, time, event, ddist,
                          lambda = NULL, pred.at = NULL,
                          fun.at = NULL, funlabel = NULL) {

  model.type = match.arg(model.type)

  if (nrow(x) != length(time) || nrow(x) != length(event))
    stop('Number of x rows and length of time/event did not match')

  if (is.null(pred.at)) stop('Missing argument pred.at')

  if (model.type %in% c('lasso', 'alasso', 'enet', 'aenet')) {

    if (!all(c('coxnet', 'glmnet') %in% class(object)))
      stop('object class must be "glmnet" and "coxnet"')

    if (is.null(lambda)) stop('Missing argument lambda')

    glmnet_pred_lp = as.vector(predict(object, newx = x, s = lambda, type = 'link'))

    all_vars = rownames(object$beta)
    selected_vars =
      all_vars[which(as.vector(abs(coef(object, s = lambda)) > .Machine$double.eps))]
    ols_formula = paste('glmnet_pred_lp ~',
                        paste(paste('`', selected_vars, '`', sep = ''),
                              collapse = ' + '))
    ols_fit = ols(as.formula(ols_formula), data = ddist,
                  sigma = 1, x = TRUE, y = TRUE)

    idx_ones = which(event == 1L)
    survtime_ones = time[idx_ones]
    names(survtime_ones) = idx_ones
    survtime_ones = sort(survtime_ones)
    survtime_at = survtime_ones[which(survtime_ones > pred.at)[1L] - 1L]
    survtime_at_idx = names(survtime_at)

    survcurve = glmnet.survcurve(object = object, time = time, event = event,
                                 x = x, survtime = survtime_ones, lambda = lambda)

  }

  if (model.type %in% c('mcp', 'mnet', 'scad', 'snet')) {

    if (!all(c('ncvsurv', 'ncvreg') %in% class(object)))
      stop('object class must be "ncvreg" and "ncvsurv"')

    if (is.null(lambda)) stop('Missing argument lambda')

    ncvreg_pred_lp = predict(object, X = x, type = 'link')

    all_vars = rownames(object$beta)
    selected_vars =
      all_vars[which(as.vector(abs(coef(object)) > .Machine$double.eps))]
    ols_formula = paste('ncvreg_pred_lp ~',
                        paste(paste('`', selected_vars, '`', sep = ''),
                              collapse = ' + '))
    ols_fit = ols(as.formula(ols_formula), data = ddist,
                  sigma = 1, x = TRUE, y = TRUE)

    idx_ones = which(event == 1L)
    survtime_ones = time[idx_ones]
    names(survtime_ones) = idx_ones
    survtime_ones = sort(survtime_ones)
    survtime_at = survtime_ones[which(survtime_ones > pred.at)[1L] - 1L]
    survtime_at_idx = names(survtime_at)

    survcurve = ncvreg.survcurve(object = object, time = time, event = event,
                                 x = x, survtime = survtime_ones, lambda = lambda)

  }

  if (model.type %in% c('flasso')) {

    if (!('penfit' %in% class(object)))
      stop('object class must be "penfit"')

    lambda = object@lambda1
    penalized_pred_lp = as.numeric(object@lin.pred)

    all_vars = colnames(x)
    selected_vars =
      all_vars[which(abs(object@penalized) > .Machine$double.eps)]
    ols_formula = paste('penalized_pred_lp ~',
                        paste(paste('`', selected_vars, '`', sep = ''),
                              collapse = ' + '))
    ols_fit = ols(as.formula(ols_formula), data = ddist,
                  sigma = 1, x = TRUE, y = TRUE)

    idx_ones = which(event == 1L)
    survtime_ones = time[idx_ones]
    names(survtime_ones) = idx_ones
    survtime_ones = sort(survtime_ones)
    survtime_at = survtime_ones[which(survtime_ones > pred.at)[1L] - 1L]
    survtime_at_idx = names(survtime_at)

    survcurve = penalized.survcurve(object = object, time = time, event = event,
                                    x = x, survtime = survtime_ones)

  }

  baseline = exp(
    log(survcurve$p[1L, which(colnames(survcurve$p) == survtime_at_idx)]) /
      exp(survcurve$lp[1L]))
  bhfun = function(z) baseline^exp(z)

  if (is.null(fun.at))
    fun.at = c(0.05, 0.2, 0.4, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99)

  if (is.null(funlabel))
    funlabel = paste('Overall Survival Probability at Time', pred.at)

  nom = list('ols_fit' = ols_fit, 'lambda' = lambda, 'survcurve' = survcurve,
             'bhfun' = bhfun, 'pred.at' = pred.at,
             'fun.at' = fun.at, 'funlabel' = funlabel)

  class(nom) = 'hdnom.nomogram'

  nom

}

#' Survival Curve Prediction for glmnet Objects
#'
#' Derived from c060::predictProb.coxnet
#'
#' @return list containing predicted survival probabilities and
#' linear predictors for all samples
#'
#' @keywords internal
glmnet.survcurve = function(object, time, event, x, survtime, lambda) {

  lp = as.numeric(predict(object, newx = data.matrix(x),
                          s = lambda, type = 'link'))
  basesurv = glmnet.basesurv(time, event, lp, sort(survtime))
  p = exp(exp(lp) %*% (-t(basesurv$cumulative_base_hazard)))
  colnames(p) = names(sort(survtime))

  if (nrow(p) != nrow(x) || ncol(p) != length(survtime))
    stop('Prediction error when estimating baseline hazard')

  list('p' = p, 'lp' = lp)

}

#' Breslow Baseline Hazard Estimator for glmnet Objects
#'
#' Derived from \code{peperr:::basesurv} and \code{gbm::basehaz.gbm}.
#'
#' @importFrom stats approx
#'
#' @return list containing cumulative base hazard
#'
#' @keywords internal
glmnet.basesurv = function(time, event, lp,
                           times.eval = NULL, centered = FALSE) {

  if (is.null(times.eval)) times.eval = sort(unique(time))

  t.unique = sort(unique(time[event == 1L]))
  alpha = length(t.unique)

  for (i in 1L:length(t.unique)) {
    alpha[i] = sum(time[event == 1L] ==
                     t.unique[i])/sum(exp(lp[time >= t.unique[i]]))
  }

  obj = approx(t.unique, cumsum(alpha), yleft = 0, xout = times.eval, rule = 2)

  if (centered) obj$y = obj$y * exp(mean(lp))
  obj$z = exp(-obj$y)

  names(obj) = c('times', 'cumulative_base_hazard', 'base_surv')

  obj

}

#' Survival Curve Prediction for ncvreg Objects
#'
#' Derived from c060::predictProb.coxnet
#'
#' @return list containing predicted survival probabilities and
#' linear predictors for all samples
#'
#' @keywords internal
ncvreg.survcurve = function(object, time, event, x, survtime, lambda) {

  lp = as.numeric(predict(object, X = data.matrix(x), type = 'link'))
  basesurv = ncvreg.basesurv(time, event, lp, sort(survtime))
  p = exp(exp(lp) %*% (-t(basesurv$cumulative_base_hazard)))
  colnames(p) = names(sort(survtime))

  if (nrow(p) != nrow(x) || ncol(p) != length(survtime))
    stop('Prediction error when estimating baseline hazard')

  list('p' = p, 'lp' = lp)

}

#' Breslow Baseline Hazard Estimator for ncvreg Objects
#'
#' Derived from \code{peperr:::basesurv} and \code{gbm::basehaz.gbm}.
#'
#' @importFrom stats approx
#'
#' @return list containing cumulative base hazard
#'
#' @keywords internal
ncvreg.basesurv = function(time, event, lp,
                           times.eval = NULL, centered = FALSE) {

  if (is.null(times.eval)) times.eval = sort(unique(time))

  t.unique = sort(unique(time[event == 1L]))
  alpha = length(t.unique)

  for (i in 1L:length(t.unique)) {
    alpha[i] = sum(time[event == 1L] ==
                     t.unique[i])/sum(exp(lp[time >= t.unique[i]]))
  }

  obj = approx(t.unique, cumsum(alpha), yleft = 0, xout = times.eval, rule = 2)

  if (centered) obj$y = obj$y * exp(mean(lp))
  obj$z = exp(-obj$y)

  names(obj) = c('times', 'cumulative_base_hazard', 'base_surv')

  obj

}

#' Survival Curve Prediction for "penalized" Objects
#'
#' Derived from c060::predictProb.coxnet
#'
#' @return list containing predicted survival probabilities and
#' linear predictors for all samples
#'
#' @keywords internal
penalized.survcurve = function(object, time, event, x, survtime) {

  lp = as.numeric(object@lin.pred)
  basesurv = penalized.basesurv(time, event, lp, sort(survtime))
  p = exp(exp(lp) %*% (-t(basesurv$cumulative_base_hazard)))
  colnames(p) = names(sort(survtime))

  if (nrow(p) != nrow(x) || ncol(p) != length(survtime))
    stop('Prediction error when estimating baseline hazard')

  list('p' = p, 'lp' = lp)

}

#' Breslow Baseline Hazard Estimator for "penalized" Objects
#'
#' Derived from \code{peperr:::basesurv} and \code{gbm::basehaz.gbm}.
#'
#' @importFrom stats approx
#'
#' @return list containing cumulative base hazard
#'
#' @keywords internal
penalized.basesurv = function(time, event, lp,
                              times.eval = NULL, centered = FALSE) {

  if (is.null(times.eval)) times.eval = sort(unique(time))

  t.unique = sort(unique(time[event == 1L]))
  alpha = length(t.unique)

  for (i in 1L:length(t.unique)) {
    alpha[i] = sum(time[event == 1L] ==
                     t.unique[i])/sum(exp(lp[time >= t.unique[i]]))
  }

  obj = approx(t.unique, cumsum(alpha), yleft = 0, xout = times.eval, rule = 2)

  if (centered) obj$y = obj$y * exp(mean(lp))
  obj$z = exp(-obj$y)

  names(obj) = c('times', 'cumulative_base_hazard', 'base_surv')

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
print.hdnom.nomogram = function(x, ...) {

  if (class(x) != 'hdnom.nomogram')
    stop('object class must be "hdnom.nomogram"')

  nom = nomogram(fit = x$ols_fit, fun = x$bhfun,
                 fun.at = x$fun.at, funlabel = x$funlabel,
                 lp = TRUE, vnames = 'labels', ...)

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
plot.hdnom.nomogram = function(x, ...) {

  if (class(x) != 'hdnom.nomogram')
    stop('object class must be "hdnom.nomogram"')

  nom = nomogram(fit = x$ols_fit, fun = x$bhfun,
                 fun.at = x$fun.at, funlabel = x$funlabel,
                 lp = TRUE, vnames = 'labels', ...)

  plot(nom)

}
