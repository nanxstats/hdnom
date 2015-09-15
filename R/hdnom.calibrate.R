#' Calibrate High-Dimensional Cox models
#'
#' Calibrate High-Dimensional Cox models
#'
#' @param x Matrix of training data used for the \code{glmnet} object;
#' on which to run the calibration.
#' @param time Survival time.
#' Must be of the same length with the number of rows as \code{x}.
#' @param event Status indicator, normally 0 = alive, 1 = dead.
#' Must be of the same length with the number of rows as \code{x}.
#' @param alpha Value of the elastic-net mixing parameter alpha in
#' glmnet. \code{alpha = 1}: lasso; \code{alpha = 0}: ridge.
#' From the Cox model you have built.
#' @param lambda Value of the penalty parameter lambda to use in the
#' glmnet fits on the resampled data. From the Cox model you have built.
#' @param pen.factor Penalty factors to apply to each coefficient.
#' From the built \emph{adaptive lasso} or \emph{adaptive elastic-net} model.
#' @param method Calibration method.
#' Options including \code{"fitting"}, \code{"bootstrap"}, \code{"cv"},
#' and \code{"repeated.cv"}.
#' @param boot.times Number of repetitions for bootstrap.
#' @param nfolds Number of folds for cross-validation and
#' repeated cross-validation.
#' @param rep.times Number of repeated times for repeated cross-validation.
#' @param pred.at Time point at which calibration should take place.
#' @param ngroup Number of groups to be formed for calibration.
#' @param trace Logical. Print trace or not. Default is \code{TRUE}.
#'
#' @importFrom survival Surv
#' @importFrom stats quantile
#' @importFrom stats median
#'
#' @export hdnom.calibrate
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
#'
#' # Fit penalized Cox model (lasso penalty) with glmnet
#' set.seed(1010)
#' cvfit = cv.glmnet(x, Surv(time, event), family = "cox", nfolds = 5)
#'
#' # Model calibration by fitting the original data directly
#' cal.fitting = hdnom.calibrate(x, time, event,
#'                               alpha = 1, lambda = cvfit$lambda.1se,
#'                               method = "fitting",
#'                               pred.at = 365 * 9, ngroup = 5)
#'
#' # Model calibration by bootstrap
#' # Normally boot.times should be set to 200 or more,
#' # we set it to 3 here only to save example running time.
#' cal.boot = hdnom.calibrate(x, time, event,
#'                            alpha = 1, lambda = cvfit$lambda.1se,
#'                            method = "bootstrap", boot.times = 3,
#'                            pred.at = 365 * 9, ngroup = 5)
#'
#' # Model calibration by 10-fold cross-validation
#' cal.cv = hdnom.calibrate(x, time, event,
#'                          alpha = 1, lambda = cvfit$lambda.1se,
#'                          method = "cv", nfolds = 5,
#'                          pred.at = 365 * 9, ngroup = 5)
#'
#' # Model calibration by repeated cross-validation
#' cal.repcv = hdnom.calibrate(x, time, event,
#'                             alpha = 1, lambda = cvfit$lambda.1se,
#'                             method = "repeated.cv", nfolds = 5, rep.times = 3,
#'                             pred.at = 365 * 9, ngroup = 5)
#'
#' print(cal.fitting)
#' summary(cal.fitting)
#' plot(cal.fitting)
#'
#' print(cal.boot)
#' summary(cal.boot)
#' plot(cal.boot)
#'
#' print(cal.cv)
#' summary(cal.cv)
#' plot(cal.cv)
#'
#' print(cal.repcv)
#' summary(cal.repcv)
#' plot(cal.repcv)
hdnom.calibrate = function(x, time, event,
                           alpha, lambda, pen.factor = NULL,
                           method = c('fitting', 'bootstrap', 'cv', 'repeated.cv'),
                           boot.times = NULL, nfolds = NULL, rep.times = NULL,
                           pred.at, ngroup = 5,
                           trace = TRUE) {

  method = match.arg(method)
  if (length(pred.at) != 1L) stop('pred.at should only contain 1 time point')

  if (method == 'fitting') {

    if (!is.null(boot.times) || !is.null(nfolds) || !is.null(rep.times))
      stop('boot.times, nfolds, and rep.times must be NULL when method = "fitting"')

    if (trace) cat('Start fitting ...\n')

    pred_list = glmnet.calibrate.internal.pred(
      x_tr = x, x_te = x, y_tr = Surv(time, event), y_te = Surv(time, event),
      alpha = alpha, lambda = lambda, pen.factor = pen.factor,
      pred.at = pred.at
    )

    pred_prob = rep(NA, nrow(x))
    for (i in 1L:length(pred_prob)) pred_prob[i] = pred_list$p[i, pred_list$idx]
    grp = cut(pred_prob, quantile(pred_prob, seq(0, 1, 1/ngroup)), labels = 1L:ngroup)

    pred_prob_median = tapply(pred_prob, grp, median)

    true_prob = glmnet.calibrate.internal.true(
      pred_prob, grp, time, event, pred.at, ngroup
    )

    prob = cbind(pred_prob_median, true_prob)
    colnames(prob)[1L] = 'Predicted'

  } else if (method == 'bootstrap') {

    if (!is.null(nfolds) || !is.null(rep.times))
      stop('nfolds and rep.times must be NULL when method = "bootstrap"')

    if (is.null(boot.times)) stop('please specify boot.times')

    # generate bootstrap sample index
    samp_mat = matrix(NA, nrow(x), boot.times)
    for(i in 1L:boot.times) samp_mat[, i] = sample(1L:nrow(x), replace = TRUE)

    # bootstrap validation main loop
    pred_list_list = vector('list', boot.times)
    for (i in 1L:ncol(samp_mat)) {

      if (trace) cat('Start bootstrap sample', i, '\n')

      samp_idx = samp_mat[, i]
      x_tr = x[samp_idx, ]
      time_tr = time[samp_idx]
      event_tr = event[samp_idx]
      y_tr = Surv(time_tr, event_tr)
      x_te = x  # use original dataset as test set
      y_te = Surv(time, event)  # use original dataset as test set

      pred_list_list[[i]] = glmnet.calibrate.internal.pred(
        x_tr = x_tr, x_te = x_te, y_tr = y_tr, y_te = y_te,
        alpha = alpha, lambda = lambda, pen.factor = pen.factor,
        pred.at = pred.at
      )

    }

    # get prediction from each bootstrap round
    pred_list_list_tmp = vector('list', boot.times)
    for (i in 1L:boot.times) {
      pred_list_list_tmp[[i]] = pred_list_list[[i]]$p[, pred_list_list[[i]]$idx]
    }

    # get mean of predicted probability across bootstrap rounds
    pred_prob = Reduce('+', pred_list_list_tmp)/length(pred_list_list_tmp)
    grp = cut(pred_prob, quantile(pred_prob, seq(0, 1, 1/ngroup)), labels = 1L:ngroup)

    pred_prob_median = tapply(pred_prob, grp, median)

    true_prob = glmnet.calibrate.internal.true(
      pred_prob, grp, time, event, pred.at, ngroup
    )

    prob = cbind(pred_prob_median, true_prob)
    colnames(prob)[1L] = 'Predicted'

  } else if (method == 'cv') {

    if (!is.null(boot.times) || !is.null(rep.times))
      stop('boot.times and rep.times must be NULL when method = "cv"')

    if (is.null(nfolds)) stop('please specify nfolds')

    # generate cross-validation sample index
    row_x = nrow(x)
    samp_idx = sample(rep_len(1L:nfolds, row_x))

    # cross-validation main loop
    pred_list_list = vector('list', nfolds)
    for (i in 1L:nfolds) {

      if (trace) cat('Start fold', i, '\n')

      x_tr = x[samp_idx != i, ]
      time_tr = time[samp_idx != i]
      event_tr = event[samp_idx != i]
      y_tr = Surv(time_tr, event_tr)
      x_te  = x[samp_idx == i, ]
      time_te = time[samp_idx == i]
      event_te = event[samp_idx == i]
      y_te = Surv(time_te, event_te)
      idx_te = which(samp_idx == i)

      pred_list_list[[i]] = glmnet.calibrate.internal.pred(
        x_tr = x_tr, x_te = x_te, y_tr = y_tr, y_te = y_te,
        alpha = alpha, lambda = lambda, pen.factor = pen.factor,
        pred.at = pred.at
      )

      rownames(pred_list_list[[i]]$p) = idx_te

    }

    # get prediction from each fold
    pred_list_list_tmp = vector('list', nfolds)
    for (i in 1L:nfolds) {
      pred_list_list_tmp[[i]] = pred_list_list[[i]]$p[, pred_list_list[[i]]$idx]
    }

    # combine and sort all predictions
    pred_prob_unsorted = unlist(pred_list_list_tmp)
    pred_prob = pred_prob_unsorted[order(as.integer(names(pred_prob_unsorted)))]
    grp = cut(pred_prob, quantile(pred_prob, seq(0, 1, 1/ngroup)), labels = 1L:ngroup)

    pred_prob_median = tapply(pred_prob, grp, median)

    true_prob = glmnet.calibrate.internal.true(
      pred_prob, grp, time, event, pred.at, ngroup
    )

    prob = cbind(pred_prob_median, true_prob)
    colnames(prob)[1L] = 'Predicted'

  } else if (method == 'repeated.cv') {

    if (!is.null(boot.times))
      stop('boot.times must be NULL when method = "repeated.cv"')

    if (is.null(nfolds) || is.null(rep.times))
      stop('please specify nfolds and rep.times')

    # generate repeated cross-validation sample index list
    row_x = nrow(x)
    samp_idx = vector('list', rep.times)
    for (k in 1L:rep.times) samp_idx[[k]] = sample(rep_len(1L:nfolds, row_x))

    # repeated cross-validation main loop
    pred_list_list = vector('list', rep.times)
    for (k in 1L:rep.times) pred_list_list[[k]] = vector('list', nfolds)

    for (j in 1L:rep.times) {
      for (i in 1L:nfolds) {

        if (trace) cat('Start repeat round', j, 'fold', i, '\n')

        x_tr = x[samp_idx[[j]] != i, ]
        time_tr = time[samp_idx[[j]] != i]
        event_tr = event[samp_idx[[j]] != i]
        y_tr = Surv(time_tr, event_tr)
        x_te  = x[samp_idx[[j]] == i, ]
        time_te = time[samp_idx[[j]] == i]
        event_te = event[samp_idx[[j]] == i]
        y_te = Surv(time_te, event_te)
        idx_te = which(samp_idx[[j]] == i)

        pred_list_list[[j]][[i]] = glmnet.calibrate.internal.pred(
          x_tr = x_tr, x_te = x_te, y_tr = y_tr, y_te = y_te,
          alpha = alpha, lambda = lambda, pen.factor,
          pred.at = pred.at
        )

        rownames(pred_list_list[[j]][[i]]$p) = idx_te

      }
    }

    # get prediction from each fold for each repeat
    pred_list_list_tmp = vector('list', rep.times)
    for (i in 1L:rep.times) pred_list_list_tmp[[i]] = vector('list', nfolds)

    for (j in 1L:rep.times) {
      for (i in 1L:nfolds) {
        pred_list_list_tmp[[j]][[i]] =
          pred_list_list[[j]][[i]]$p[, pred_list_list[[j]][[i]]$idx]
      }
    }

    # combine and sort all predictions for each repeat time
    pred_prob_unsorted = vector('list', rep.times)
    for (i in 1L:rep.times) pred_prob_unsorted[[i]] = unlist(pred_list_list_tmp[[i]])
    pred_prob_tmp = vector('list', rep.times)
    for (i in 1L:rep.times) {
      pred_prob_tmp[[i]] =
        pred_prob_unsorted[[i]][order(as.integer(names(pred_prob_unsorted[[i]])))]
    }

    # get mean of predicted probability across repeat times
    pred_prob = Reduce('+', pred_prob_tmp)/length(pred_prob_tmp)
    grp = cut(pred_prob, quantile(pred_prob, seq(0, 1, 1/ngroup)), labels = 1L:ngroup)

    pred_prob_median = tapply(pred_prob, grp, median)

    true_prob = glmnet.calibrate.internal.true(
      pred_prob, grp, time, event, pred.at, ngroup
    )

    prob = cbind(pred_prob_median, true_prob)
    colnames(prob)[1L] = 'Predicted'

  } else {
    stop('method must be one of "fitting", "bootstrap", "cv", or "repeated.cv"')
  }

  switch(method,
         fitting = {
           class(prob) = c('hdnom.calibrate',
                           'glmnet.calibrate.fitting')
           attr(prob, 'alpha') = alpha
           attr(prob, 'lambda') = lambda
           attr(prob, 'pen.factor') = pen.factor
           attr(prob, 'pred.at') = pred.at
           attr(prob, 'ngroup') = ngroup
         },
         bootstrap = {
           class(prob) = c('hdnom.calibrate',
                           'glmnet.calibrate.bootstrap')
           attr(prob, 'alpha')      = alpha
           attr(prob, 'lambda')     = lambda
           attr(prob, 'pen.factor') = pen.factor
           attr(prob, 'pred.at')    = pred.at
           attr(prob, 'ngroup')     = ngroup
           attr(prob, 'boot.times') = boot.times
         },
         cv = {
           class(prob) = c('hdnom.calibrate',
                           'glmnet.calibrate.cv')
           attr(prob, 'alpha')   = alpha
           attr(prob, 'lambda')  = lambda
           attr(prob, 'pen.factor') = pen.factor
           attr(prob, 'pred.at') = pred.at
           attr(prob, 'ngroup')  = ngroup
           attr(prob, 'nfolds')  = nfolds
         },
         repeated.cv = {
           class(prob) = c('hdnom.calibrate',
                           'glmnet.calibrate.repeated.cv')
           attr(prob, 'alpha')      = alpha
           attr(prob, 'lambda')     = lambda
           attr(prob, 'pen.factor') = pen.factor
           attr(prob, 'pred.at')    = pred.at
           attr(prob, 'ngroup')     = ngroup
           attr(prob, 'nfolds')     = nfolds
           attr(prob, 'rep.times')  = rep.times
         }
  )

  prob

}

#' Compute predicted survival probabilities for calibration
#'
#' @importFrom glmnet glmnet
#' @importFrom stats predict
#'
#' @return list containing predicted survival probability
#'
#' @keywords internal
glmnet.calibrate.internal.pred = function(x_tr, x_te, y_tr, y_te,
                                          alpha, lambda, pen.factor,
                                          pred.at) {

  if (is.null(pen.factor)) {
    object = glmnet(x = x_tr, y = y_tr, family = 'cox',
                    alpha = alpha, lambda = lambda)
  } else {
    object = glmnet(x = x_tr, y = y_tr, family = 'cox',
                    alpha = alpha, lambda = lambda,
                    penalty.factor = pen.factor)
  }

  lp = as.numeric(predict(object, newx = data.matrix(x_te),
                          s = lambda, type = 'link'))

  time_te = y_te[, 1L]
  event_te = y_te[, 2L]
  idx_ones = which(event_te == 1L)
  if (length(idx_ones) == 0L)
    stop('No 1 events in the testing fold, please set another random seed.')
  survtime_ones = time_te[idx_ones]
  names(survtime_ones) = idx_ones
  survtime_ones = sort(survtime_ones)

  basesurv = glmnet.basesurv(time_te, event_te, lp, survtime_ones)
  p = exp(exp(lp) %*% (-t(basesurv$cumulative_base_hazard)))

  if (nrow(p) != nrow(x_te) || ncol(p) != length(survtime_ones))
    stop('Prediction error when estimating baseline hazard')

  idx = length(which(survtime_ones <= pred.at))

  list('p' = p, 'idx' = idx)

}

#' Compute Kaplan-Meier estimated survival probabilities for calibration
#'
#' @importFrom survival survfit
#' @importFrom survival Surv
#'
#' @return list
#'
#' @keywords internal
glmnet.calibrate.internal.true = function(pred_prob, grp,
                                          time, event,
                                          pred.at, ngroup) {

  true_prob = matrix(NA, ncol = 3L, nrow = ngroup)
  colnames(true_prob) = c("Observed", "Lower 95%", "Upper 95%")

  for (i in 1L:ngroup) {
    time_grp = time[which(grp == i)]
    event_grp = event[which(grp == i)]
    km = survfit(Surv(time_grp, event_grp) ~ 1, type = 'kaplan-meier')
    idx = which(km$time > pred.at)[1L] - 1L
    km_pred_at = km$surv[idx]
    ll_pred_at = km$lower[idx]
    ul_pred_at = km$upper[idx]
    true_prob[i, ] = c(km_pred_at, ll_pred_at, ul_pred_at)
  }

  return(true_prob)

}

#' Print calibration result generated by hdnom.calibrate
#'
#' Print calibration result generated by hdnom.calibrate
#'
#' @param x an object returned by \code{\link{hdnom.calibrate}}.
#' @param ... other parameters (not used).
#'
#' @method print hdnom.calibrate
#'
#' @export
#'
#' @examples
#' NULL
print.hdnom.calibrate = function (x, ...) {

  if (!('hdnom.calibrate' %in% class(x)))
    stop('object class must be "hdnom.calibrate"')

  method = setdiff(class(x), 'hdnom.calibrate')

  if (method == 'glmnet.calibrate.fitting') {

    cat('High-Dimensional Cox Model Calibration Object\n')
    cat('Calibration method: fitting\n')
    cat('glmnet model alpha:', attr(x, 'alpha'), '\n')
    cat('glmnet model lambda:', attr(x, 'lambda'), '\n')
    if (is.null(attr(x, 'pen.factor'))) {
      cat('glmnet model penalty factor: not specified\n')
    } else {
      cat('glmnet model penalty factor: specified\n')
    }
    cat('Calibration time point:', attr(x, 'pred.at'), '\n')
    cat('Number of groups formed for calibration:', attr(x, 'ngroup'), '\n')

  } else if  (method == 'glmnet.calibrate.bootstrap') {

    cat('High-Dimensional Cox Model Calibration Object\n')
    cat('Calibration method: bootstrap\n')
    cat('Bootstrap samples:', attr(x, 'boot.times'), '\n')
    cat('glmnet model alpha:', attr(x, 'alpha'), '\n')
    cat('glmnet model lambda:', attr(x, 'lambda'), '\n')
    if (is.null(attr(x, 'pen.factor'))) {
      cat('glmnet model penalty factor: not specified\n')
    } else {
      cat('glmnet model penalty factor: specified\n')
    }
    cat('Calibration time point:', attr(x, 'pred.at'), '\n')
    cat('Number of groups formed for calibration:', attr(x, 'ngroup'), '\n')

  } else if  (method == 'glmnet.calibrate.cv') {

    cat('High-Dimensional Cox Model Calibration Object\n')
    cat('Calibration method: k-fold cross-validation\n')
    cat('Cross-validation folds:', attr(x, 'nfolds'), '\n')
    cat('glmnet model alpha:', attr(x, 'alpha'), '\n')
    cat('glmnet model lambda:', attr(x, 'lambda'), '\n')
    if (is.null(attr(x, 'pen.factor'))) {
      cat('glmnet model penalty factor: not specified\n')
    } else {
      cat('glmnet model penalty factor: specified\n')
    }
    cat('Calibration time point:', attr(x, 'pred.at'), '\n')
    cat('Number of groups formed for calibration:', attr(x, 'ngroup'), '\n')

  } else if  (method == 'glmnet.calibrate.repeated.cv') {

    cat('High-Dimensional Cox Model Calibration Object\n')
    cat('Calibration method: repeated cross-validation\n')
    cat('Cross-validation folds:', attr(x, 'nfolds'), '\n')
    cat('Cross-validation repeated times:', attr(x, 'rep.times'), '\n')
    cat('glmnet model alpha:', attr(x, 'alpha'), '\n')
    cat('glmnet model lambda:', attr(x, 'lambda'), '\n')
    if (is.null(attr(x, 'pen.factor'))) {
      cat('glmnet model penalty factor: not specified\n')
    } else {
      cat('glmnet model penalty factor: specified\n')
    }
    cat('Calibration time point:', attr(x, 'pred.at'), '\n')
    cat('Number of groups formed for calibration:', attr(x, 'ngroup'), '\n')

  } else {

    stop('hdnom.calibrate object is not valid')

  }

}

#' Summary calibration result generated by hdnom.calibrate
#'
#' Summary calibration result generated by hdnom.calibrate
#'
#' @param object an object returned by \code{\link{hdnom.calibrate}}.
#' @param ... other parameters (not used).
#'
#' @method summary hdnom.calibrate
#'
#' @export
#'
#' @examples
#' NULL
summary.hdnom.calibrate = function(object, ...) {

  if (!('hdnom.calibrate' %in% class(object)))
    stop('object class must be "hdnom.calibrate"')

  method = setdiff(class(object), 'hdnom.calibrate')

  if (method == 'glmnet.calibrate.fitting') {

    attr(object, 'alpha')      = NULL
    attr(object, 'lambda')     = NULL
    attr(object, 'pen.factor') = NULL
    attr(object, 'pred.at')    = NULL
    attr(object, 'ngroup')     = NULL

  } else if  (method == 'glmnet.calibrate.bootstrap') {

    attr(object, 'boot.times') = NULL
    attr(object, 'alpha')      = NULL
    attr(object, 'lambda')     = NULL
    attr(object, 'pen.factor') = NULL
    attr(object, 'pred.at')    = NULL
    attr(object, 'ngroup')     = NULL

  } else if  (method == 'glmnet.calibrate.cv') {

    attr(object, 'nfolds')     = NULL
    attr(object, 'alpha')      = NULL
    attr(object, 'lambda')     = NULL
    attr(object, 'pen.factor') = NULL
    attr(object, 'pred.at')    = NULL
    attr(object, 'ngroup')     = NULL

  } else if  (method == 'glmnet.calibrate.repeated.cv') {

    attr(object, 'nfolds')     = NULL
    attr(object, 'rep.times')  = NULL
    attr(object, 'alpha')      = NULL
    attr(object, 'lambda')     = NULL
    attr(object, 'pen.factor') = NULL
    attr(object, 'pred.at')    = NULL
    attr(object, 'ngroup')     = NULL

  } else {

    stop('hdnom.calibrate object is not valid')

  }

  cat('  Calibration Summary Table\n')
  class(object) = 'matrix'
  print(object)

}

#' Plot calibration curves
#'
#' Plot calibration curves
#'
#' @param x an object returned by \code{\link{hdnom.calibrate}}.
#' @param xlim x axis limits of the plot.
#' @param ylim y axis limits of the plot.
#' @param ... other parameters for \code{plot}.
#'
#' @method plot hdnom.calibrate
#'
#' @export
#'
#' @importFrom graphics arrows abline
#'
#' @examples
#' NULL
plot.hdnom.calibrate = function(x, xlim = c(0, 1), ylim = c(0, 1), ...) {

  if (!('hdnom.calibrate' %in% class(x)))
    stop('object class must be "hdnom.calibrate"')

  pre = x[, 'Predicted']
  obs = x[, 'Observed']
  ll  = x[, 'Lower 95%']
  ul  = x[, 'Upper 95%']

  plot(pre, obs, xlim = xlim, ylim = ylim,
       xlab = 'Predicted Survival Probability',
       ylab = 'Observed Survival Probability',
       pch = 16, type = 'b', ...)
  arrows(x0 = pre, y0 = ll, y1 = ul,
         angle = 90, code = 3, length = 0.05, lwd = 1)
  abline(a = 0, b = 1, lty = 2)

}
