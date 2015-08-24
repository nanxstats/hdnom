#' Nomograms for Cox models fitted with glmnet
#'
#' Nomograms for Cox models fitted with glmnet
#'
#' @param object Fitted \code{"glmnet"} model object.
#' @param x Matrix of training data used for the \code{glmnet} object.
#' @param time Survival time.
#' Must be of the same length with the number of rows as \code{x}.
#' @param event Status indicator, normally 0 = alive, 1 = dead.
#' Must be of the same length with the number of rows as \code{x}.
#' @param ddist Data frame version of x, made by \code{\link[rms]{datadist}}.
#' @param s Value of the penalty parameter lambda in
#' \code{\link[glmnet]{glmnet}}.
#' We will use the selected variables at the provided \code{s} to
#' build the nomogram, and make predictions.
#' See the example for choosing a proper lambda value extracted
#' from cross-validation results.
#' @param pred.at Time point at which to draw nomogram prediction axis.
#' @param fun.at Function values to label on axis.
#' @param funlabel Label for \code{fun} axis.
#' @param ... other arguments for \code{\link[rms]{nomogram}}.
#'
#' @export glmnet.nomogram
#'
#' @importFrom rms ols nomogram
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
#' # Generate glmnet.nomogram objects and draw nomogram
#' nom = glmnet.nomogram(fit, x, time, event, x.df,
#'                       s = cvfit$lambda.1se, pred.at = 365 * 2,
#'                       funlabel = "2-Year Overall Survival Probability")
#'
#' print(nom)
#' plot(nom)
glmnet.nomogram = function(object, x, time, event, ddist,
                           s = NULL, pred.at = NULL,
                           fun.at = NULL, funlabel = NULL) {

  if (!all(c('coxnet', 'glmnet') %in% class(object)))
    stop('object class must be "glmnet" and "coxnet"')

  if (nrow(x) != length(time) || nrow(x) != length(event))
    stop('Number of x rows and length of time/event did not match')

  if (is.null(s)) stop('Missing argument s')

  if (is.null(s)) stop('Missing argument pred.at')

  glmnet_pred_lp = as.vector(predict(object, newx = x, s = s, type = 'link'))

  all_vars = rownames(object$beta)
  selected_vars =
    all_vars[which(as.vector(abs(coef(object, s = s)) > .Machine$double.eps))]
  ols_formula = paste('glmnet_pred_lp ~',
                      paste(selected_vars, collapse = ' + '))
  ols_fit = ols(as.formula(ols_formula), data = ddist,
                sigma = 1, x = TRUE, y = TRUE)

  idx_ones = which(event == 1L)
  survtime_ones = time[idx_ones]
  names(survtime_ones) = idx_ones
  survtime_ones = sort(survtime_ones)
  survtime_at = survtime_ones[which(survtime_ones > pred.at)[1L] - 1L]
  survtime_at_idx = names(survtime_at)

  survcurve = survcurve.glmnet(object = object, time = time, event = event,
                               x = x, survtime = survtime_ones, lambda = s)
  baseline = exp(
    log(survcurve$p[1L, which(colnames(survcurve$p) == survtime_at_idx)]) /
      exp(survcurve$lp[1L]))
  bhfun = function(z) baseline^exp(z)

  if (is.null(fun.at))
    fun.at = c(0.05, 0.2, 0.4, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99)

  if (is.null(funlabel))
    funlabel = paste('Overall Survival Probability at Time', pred.at)

  nom = list('ols_fit' = ols_fit, 's' = s, 'bhfun' = bhfun,
             'pred.at' = pred.at, 'fun.at' = fun.at, 'funlabel' = funlabel)

  class(nom) = 'glmnet.nomogram'

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
survcurve.glmnet = function(object, time, event, x, survtime, lambda) {

  lp = as.numeric(predict(object, newx = data.matrix(x),
                          s = lambda, type = 'link'))
  basesurv = basesurv.glmnet(time, event, lp, sort(survtime))
  p = exp(exp(lp) %*% (-t(basesurv$cumulative_base_hazard)))
  colnames(p) = names(sort(survtime))

  if (nrow(p) != nrow(x) || ncol(p) != length(survtime))
    stop('Prediction error when estimating baseline hazard')

  list('p' = p, 'lp' = lp)

}

#' Breslow estimator for baseline hazard
#'
#' Derived from \code{peperr:::basesurv} and \code{gbm::basehaz.gbm}.
#'
#' @return list containing cumulative base hazard
#'
#' @keywords internal
basesurv.glmnet = function(time, event, lp,
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

#' Plot nomogram objects generated with glmnet.nomogram
#'
#' Plot nomogram objects generated with glmnet.nomogram
#'
#' @param x a \code{"glmnet.nomogram"} object.
#' @param ... other parameters for \code{\link[rms]{nomogram}}.
#'
#' @method plot glmnet.nomogram
#'
#' @export
#'
#' @examples
#' NULL
plot.glmnet.nomogram = function(x, ...) {

  if (class(x) != 'glmnet.nomogram')
    stop('object class must be "glmnet.nomogram"')

  nom = nomogram(fit = x$ols_fit, fun = x$bhfun,
                 fun.at = x$fun.at, funlabel = x$funlabel,
                 lp = TRUE, vnames = 'labels', ...)

  plot(nom)

}

#' Print nomograms objects generated with glmnet.nomogram
#'
#' Print nomograms objects generated with glmnet.nomogram
#'
#' @param x a \code{"glmnet.nomogram"} object.
#' @param ... other parameters for \code{\link[rms]{nomogram}}.
#'
#' @method print glmnet.nomogram
#'
#' @export
#'
#' @examples
#' NULL
print.glmnet.nomogram = function(x, ...) {

  if (class(x) != 'glmnet.nomogram')
    stop('object class must be "glmnet.nomogram"')

  nom = nomogram(fit = x$ols_fit, fun = x$bhfun,
                 fun.at = x$fun.at, funlabel = x$funlabel,
                 lp = TRUE, vnames = 'labels', ...)

  print(nom)

}
