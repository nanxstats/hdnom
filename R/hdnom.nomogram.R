#' Nomograms for High-Dimensional Cox models
#'
#' Nomograms for High-Dimensional Cox models
#'
#' @param object Fitted \code{glmnet} model object.
#' @param x Matrix of training data used for the \code{glmnet} object.
#' @param time Survival time.
#' Must be of the same length with the number of rows as \code{x}.
#' @param event Status indicator, normally 0 = alive, 1 = dead.
#' Must be of the same length with the number of rows as \code{x}.
#' @param ddist Data frame version of x, made by \code{\link[rms]{datadist}}.
#' @param lambda Value of the penalty parameter lambda in
#' \code{\link[glmnet]{glmnet}}.
#' We will use the selected variables at the provided \code{s} to
#' build the nomogram, and make predictions.
#' See the example for choosing a proper lambda value extracted
#' from cross-validation results.
#' @param pred.at Time point at which to plot nomogram prediction axis.
#' @param fun.at Function values to label on axis.
#' @param funlabel Label for \code{fun} axis.
#' @param ... Other arguments for \code{\link[rms]{nomogram}}.
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
#' nom = hdnom.nomogram(fit, x, time, event, x.df,
#'                      lambda = cvfit$lambda.1se, pred.at = 365 * 2,
#'                      funlabel = "2-Year Overall Survival Probability")
#'
#' print(nom)
#' plot(nom)
hdnom.nomogram = function(object, x, time, event, ddist,
                          lambda = NULL, pred.at = NULL,
                          fun.at = NULL, funlabel = NULL) {

  if (!all(c('coxnet', 'glmnet') %in% class(object)))
    stop('object class must be "glmnet" and "coxnet"')

  if (nrow(x) != length(time) || nrow(x) != length(event))
    stop('Number of x rows and length of time/event did not match')

  if (is.null(lambda)) stop('Missing argument lambda')

  if (is.null(lambda)) stop('Missing argument pred.at')

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

#' Survival curve prediction for glmnet objects
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

#' Breslow baseline hazard estimator
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

#' Print nomograms objects generated by hdnom.nomogram
#'
#' Print nomograms objects generated by hdnom.nomogram
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

#' Plot nomogram objects generated by hdnom.nomogram
#'
#' Plot nomogram objects generated by hdnom.nomogram
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
