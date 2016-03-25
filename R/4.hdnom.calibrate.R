#' Calibrate High-Dimensional Cox Models
#'
#' Calibrate High-Dimensional Cox Models
#'
#' @param x Matrix of training data used for fitting the model;
#' on which to run the calibration.
#' @param time Survival time.
#' Must be of the same length with the number of rows as \code{x}.
#' @param event Status indicator, normally 0 = alive, 1 = dead.
#' Must be of the same length with the number of rows as \code{x}.
#' @param model.type Model type to calibrate. Could be one of \code{"lasso"},
#' \code{"alasso"}, \code{"flasso"}, \code{"enet"}, \code{"aenet"},
#' \code{"mcp"}, \code{"mnet"}, \code{"scad"}, or \code{"snet"}.
#' @param alpha Value of the elastic-net mixing parameter alpha for
#' \code{enet}, \code{aenet}, \code{mnet}, and \code{snet} models.
#' For \code{lasso}, \code{alasso}, \code{mcp}, and \code{scad} models,
#' please set \code{alpha = 1}.
#' \code{alpha=1}: lasso (l1) penalty; \code{alpha=0}: ridge (l2) penalty.
#' Note that for \code{mnet} and \code{snet} models,
#' \code{alpha} can be set to very close to 0 but not 0 exactly.
#' @param lambda Value of the penalty parameter lambda to use in the
#' model fits on the resampled data. From the Cox model you have built.
#' @param pen.factor Penalty factors to apply to each coefficient.
#' From the built \emph{adaptive lasso} or \emph{adaptive elastic-net} model.
#' @param gamma Value of the model parameter gamma for MCP/SCAD/Mnet/Snet models.
#' @param method Calibration method.
#' Options including \code{"fitting"}, \code{"bootstrap"}, \code{"cv"},
#' and \code{"repeated.cv"}.
#' @param boot.times Number of repetitions for bootstrap.
#' @param nfolds Number of folds for cross-validation and
#' repeated cross-validation.
#' @param rep.times Number of repeated times for repeated cross-validation.
#' @param pred.at Time point at which calibration should take place.
#' @param ngroup Number of groups to be formed for calibration.
#' @param seed A random seed for resampling.
#' @param trace Logical. Output the calibration progress or not.
#' Default is \code{TRUE}.
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
#'
#' # Load imputed SMART data
#' data(smart)
#' x = as.matrix(smart[, -c(1, 2)])[1:1000, ]
#' time = smart$TEVENT[1:1000]
#' event = smart$EVENT[1:1000]
#'
#' # Fit penalized Cox model (lasso penalty) with glmnet
#' set.seed(1010)
#' cvfit = cv.glmnet(x, Surv(time, event), family = "cox", nfolds = 5)
#'
#' # Model calibration by fitting the original data directly
#' cal.fitting = hdnom.calibrate(x, time, event, model.type = "lasso",
#'                               alpha = 1, lambda = cvfit$lambda.1se,
#'                               method = "fitting",
#'                               pred.at = 365 * 9, ngroup = 5,
#'                               seed = 1010)
#'
#' # Model calibration by bootstrap
#' # Normally boot.times should be set to 200 or more,
#' # we set it to 3 here only to save example running time.
#' cal.boot = hdnom.calibrate(x, time, event, model.type = "lasso",
#'                            alpha = 1, lambda = cvfit$lambda.1se,
#'                            method = "bootstrap", boot.times = 3,
#'                            pred.at = 365 * 9, ngroup = 5,
#'                            seed = 1010)
#'
#' # Model calibration by 5-fold cross-validation
#' cal.cv = hdnom.calibrate(x, time, event, model.type = "lasso",
#'                          alpha = 1, lambda = cvfit$lambda.1se,
#'                          method = "cv", nfolds = 5,
#'                          pred.at = 365 * 9, ngroup = 5,
#'                          seed = 1010)
#'
#' # Model calibration by repeated cross-validation
#' cal.repcv = hdnom.calibrate(x, time, event, model.type = "lasso",
#'                             alpha = 1, lambda = cvfit$lambda.1se,
#'                             method = "repeated.cv", nfolds = 3, rep.times = 3,
#'                             pred.at = 365 * 9, ngroup = 5,
#'                             seed = 1010)
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
# ### Testing fused lasso, SCAD, and Mnet models ###
# library("survival")
# library("rms")
#
# # Load imputed SMART data
# data(smart)
# x = as.matrix(smart[, -c(1, 2)])[1:500, ]
# time = smart$TEVENT[1:500]
# event = smart$EVENT[1:500]
# y = Surv(time, event)
#
# set.seed(1010)
# cal.fitting = hdnom.calibrate(x, time, event, model.type = "flasso",
#                               lambda = 60,
#                               method = "fitting",
#                               pred.at = 365 * 9, ngroup = 5,
#                               seed = 1010)
#
# cal.boot = hdnom.calibrate(x, time, event, model.type = "scad",
#                            gamma = 3.7, alpha = 1, lambda = 0.03,
#                            method = "bootstrap", boot.times = 10,
#                            pred.at = 365 * 9, ngroup = 5,
#                            seed = 1010)
#
# cal.cv = hdnom.calibrate(x, time, event, model.type = "mnet",
#                          gamma = 3, alpha = 0.3, lambda = 0.03,
#                          method = "cv", nfolds = 5,
#                          pred.at = 365 * 9, ngroup = 5,
#                          seed = 1010)
#
# cal.repcv = hdnom.calibrate(x, time, event, model.type = "flasso",
#                             lambda = 60,
#                             method = "repeated.cv", nfolds = 5, rep.times = 3,
#                             pred.at = 365 * 9, ngroup = 5,
#                             seed = 1010)
#
# print(cal.fitting)
# summary(cal.fitting)
# plot(cal.fitting)
#
# print(cal.boot)
# summary(cal.boot)
# plot(cal.boot)
#
# print(cal.cv)
# summary(cal.cv)
# plot(cal.cv)
#
# print(cal.repcv)
# summary(cal.repcv)
# plot(cal.repcv)
hdnom.calibrate = function(x, time, event,
                           model.type = c('lasso', 'alasso', 'flasso',
                                          'enet', 'aenet',
                                          'mcp', 'mnet',
                                          'scad', 'snet'),
                           alpha, lambda, pen.factor = NULL, gamma,
                           method = c('fitting', 'bootstrap', 'cv', 'repeated.cv'),
                           boot.times = NULL, nfolds = NULL, rep.times = NULL,
                           pred.at, ngroup = 5,
                           seed = 1001, trace = TRUE) {

  model.type = match.arg(model.type)
  method = match.arg(method)
  if (length(pred.at) != 1L) stop('pred.at should only contain 1 time point')

  set.seed(seed)

  if (method == 'fitting') {

    if (!is.null(boot.times) || !is.null(nfolds) || !is.null(rep.times))
      stop('boot.times, nfolds, and rep.times must be NULL when method = "fitting"')

    if (trace) cat('Start fitting ...\n')

    if (model.type %in% c('lasso', 'alasso', 'enet', 'aenet')) {
      pred_list = glmnet.calibrate.internal.pred(
        x_tr = x, x_te = x, y_tr = Surv(time, event),
        alpha = alpha, lambda = lambda, pen.factor = pen.factor,
        pred.at = pred.at
      )
    }

    if (model.type %in% c('mcp', 'mnet', 'scad', 'snet')) {
      pred_list = ncvreg.calibrate.internal.pred(
        x_tr = x, x_te = x, y_tr = Surv(time, event),
        model.type = model.type,
        gamma = gamma, alpha = alpha, lambda = lambda,
        pred.at = pred.at
      )
    }

    if (model.type %in% c('flasso')) {
      pred_list = penalized.calibrate.internal.pred(
        x_tr = x, x_te = x, y_tr = Surv(time, event),
        lambda = lambda,
        pred.at = pred.at
      )
    }

    pred_prob = rep(NA, nrow(x))
    for (i in 1L:length(pred_prob)) pred_prob[i] = pred_list$p[i, pred_list$idx]
    grp = cut(pred_prob, quantile(pred_prob, seq(0, 1, 1/ngroup)), labels = 1L:ngroup)

    pred_prob_median = tapply(pred_prob, grp, median)

    true_prob = hdnom.calibrate.internal.true(
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
      x_tr = x[samp_idx, , drop = FALSE]
      time_tr = time[samp_idx]
      event_tr = event[samp_idx]
      y_tr = Surv(time_tr, event_tr)
      x_te = x  # use original dataset as test set
      y_te = Surv(time, event)  # use original dataset as test set

      if (model.type %in% c('lasso', 'alasso', 'enet', 'aenet')) {
        pred_list_list[[i]] = glmnet.calibrate.internal.pred(
          x_tr = x_tr, x_te = x_te, y_tr = y_tr,
          alpha = alpha, lambda = lambda, pen.factor = pen.factor,
          pred.at = pred.at
        )
      }

      if (model.type %in% c('mcp', 'mnet', 'scad', 'snet')) {
        pred_list_list[[i]] = ncvreg.calibrate.internal.pred(
          x_tr = x_tr, x_te = x_te, y_tr = y_tr,
          model.type = model.type,
          alpha = alpha, lambda = lambda, gamma = gamma,
          pred.at = pred.at
        )
      }

      if (model.type %in% c('flasso')) {
        pred_list_list[[i]] = penalized.calibrate.internal.pred(
          x_tr = x_tr, x_te = x_te, y_tr = y_tr,
          lambda = lambda,
          pred.at = pred.at
        )
      }

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

    true_prob = hdnom.calibrate.internal.true(
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

      x_tr = x[samp_idx != i, , drop = FALSE]
      time_tr = time[samp_idx != i]
      event_tr = event[samp_idx != i]
      y_tr = Surv(time_tr, event_tr)
      x_te  = x[samp_idx == i, , drop = FALSE]
      time_te = time[samp_idx == i]
      event_te = event[samp_idx == i]
      y_te = Surv(time_te, event_te)
      idx_te = which(samp_idx == i)

      if (model.type %in% c('lasso', 'alasso', 'enet', 'aenet')) {
        pred_list_list[[i]] = glmnet.calibrate.internal.pred(
          x_tr = x_tr, x_te = x_te, y_tr = y_tr,
          alpha = alpha, lambda = lambda, pen.factor = pen.factor,
          pred.at = pred.at
        )
      }

      if (model.type %in% c('mcp', 'mnet', 'scad', 'snet')) {
        pred_list_list[[i]] = ncvreg.calibrate.internal.pred(
          x_tr = x_tr, x_te = x_te, y_tr = y_tr,
          model.type = model.type,
          alpha = alpha, lambda = lambda, gamma = gamma,
          pred.at = pred.at
        )
      }

      if (model.type %in% c('flasso')) {
        pred_list_list[[i]] = penalized.calibrate.internal.pred(
          x_tr = x_tr, x_te = x_te, y_tr = y_tr,
          lambda = lambda,
          pred.at = pred.at
        )
      }

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

    true_prob = hdnom.calibrate.internal.true(
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

        x_tr = x[samp_idx[[j]] != i, , drop = FALSE]
        time_tr = time[samp_idx[[j]] != i]
        event_tr = event[samp_idx[[j]] != i]
        y_tr = Surv(time_tr, event_tr)
        x_te  = x[samp_idx[[j]] == i, , drop = FALSE]
        time_te = time[samp_idx[[j]] == i]
        event_te = event[samp_idx[[j]] == i]
        y_te = Surv(time_te, event_te)
        idx_te = which(samp_idx[[j]] == i)

        if (model.type %in% c('lasso', 'alasso', 'enet', 'aenet')) {
          pred_list_list[[j]][[i]] = glmnet.calibrate.internal.pred(
            x_tr = x_tr, x_te = x_te, y_tr = y_tr,
            alpha = alpha, lambda = lambda, pen.factor = pen.factor,
            pred.at = pred.at
          )
        }

        if (model.type %in% c('mcp', 'mnet', 'scad', 'snet')) {
          pred_list_list[[j]][[i]] = ncvreg.calibrate.internal.pred(
            x_tr = x_tr, x_te = x_te, y_tr = y_tr,
            model.type = model.type,
            alpha = alpha, lambda = lambda, gamma = gamma,
            pred.at = pred.at
          )
        }

        if (model.type %in% c('flasso')) {
          pred_list_list[[j]][[i]] = penalized.calibrate.internal.pred(
            x_tr = x_tr, x_te = x_te, y_tr = y_tr,
            lambda = lambda,
            pred.at = pred.at
          )
        }

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

    true_prob = hdnom.calibrate.internal.true(
      pred_prob, grp, time, event, pred.at, ngroup
    )

    prob = cbind(pred_prob_median, true_prob)
    colnames(prob)[1L] = 'Predicted'

  } else {
    stop('method must be one of "fitting", "bootstrap", "cv", or "repeated.cv"')
  }

  switch(method,

         fitting = {

           if (model.type %in% c('lasso', 'alasso', 'enet', 'aenet')) {
             class(prob) = c('hdnom.calibrate',
                             'glmnet.calibrate.fitting')
             attr(prob, 'alpha') = alpha
             attr(prob, 'lambda') = lambda
             attr(prob, 'pen.factor') = pen.factor
           }

           if (model.type %in% c('mcp', 'mnet', 'scad', 'snet')) {
             class(prob) = c('hdnom.calibrate',
                             'ncvreg.calibrate.fitting')
             attr(prob, 'alpha') = alpha
             attr(prob, 'lambda') = lambda
             attr(prob, 'gamma') = gamma
           }

           if (model.type %in% c('flasso')) {
             class(prob) = c('hdnom.calibrate',
                             'penalized.calibrate.fitting')
             attr(prob, 'lambda') = lambda
           }

         },

         bootstrap = {

           if (model.type %in% c('lasso', 'alasso', 'enet', 'aenet')) {
             class(prob) = c('hdnom.calibrate',
                             'glmnet.calibrate.bootstrap')
             attr(prob, 'alpha')      = alpha
             attr(prob, 'lambda')     = lambda
             attr(prob, 'pen.factor') = pen.factor
             attr(prob, 'boot.times') = boot.times
           }

           if (model.type %in% c('mcp', 'mnet', 'scad', 'snet')) {
             class(prob) = c('hdnom.calibrate',
                             'ncvreg.calibrate.bootstrap')
             attr(prob, 'alpha')      = alpha
             attr(prob, 'lambda')     = lambda
             attr(prob, 'gamma')      = gamma
             attr(prob, 'boot.times') = boot.times
           }

           if (model.type %in% c('flasso')) {
             class(prob) = c('hdnom.calibrate',
                             'penalized.calibrate.bootstrap')
             attr(prob, 'lambda')     = lambda
             attr(prob, 'boot.times') = boot.times
           }

         },

         cv = {

           if (model.type %in% c('lasso', 'alasso', 'enet', 'aenet')) {
             class(prob) = c('hdnom.calibrate',
                             'glmnet.calibrate.cv')
             attr(prob, 'alpha')      = alpha
             attr(prob, 'lambda')     = lambda
             attr(prob, 'pen.factor') = pen.factor
             attr(prob, 'nfolds')     = nfolds
           }

           if (model.type %in% c('mcp', 'mnet', 'scad', 'snet')) {
             class(prob) = c('hdnom.calibrate',
                             'ncvreg.calibrate.cv')
             attr(prob, 'alpha')      = alpha
             attr(prob, 'lambda')     = lambda
             attr(prob, 'gamma')      = gamma
             attr(prob, 'nfolds')     = nfolds
           }

           if (model.type %in% c('flasso')) {
             class(prob) = c('hdnom.calibrate',
                             'penalized.calibrate.cv')
             attr(prob, 'lambda')     = lambda
             attr(prob, 'nfolds')     = nfolds
           }

         },

         repeated.cv = {

           if (model.type %in% c('lasso', 'alasso', 'enet', 'aenet')) {
             class(prob) = c('hdnom.calibrate',
                             'glmnet.calibrate.repeated.cv')
             attr(prob, 'alpha')      = alpha
             attr(prob, 'lambda')     = lambda
             attr(prob, 'pen.factor') = pen.factor
             attr(prob, 'nfolds')     = nfolds
             attr(prob, 'rep.times')  = rep.times
           }

           if (model.type %in% c('mcp', 'mnet', 'scad', 'snet')) {
             class(prob) = c('hdnom.calibrate',
                             'ncvreg.calibrate.repeated.cv')
             attr(prob, 'alpha')      = alpha
             attr(prob, 'lambda')     = lambda
             attr(prob, 'gamma')      = gamma
             attr(prob, 'nfolds')     = nfolds
             attr(prob, 'rep.times')  = rep.times
           }

           if (model.type %in% c('flasso')) {
             class(prob) = c('hdnom.calibrate',
                             'penalized.calibrate.repeated.cv')
             attr(prob, 'lambda')     = lambda
             attr(prob, 'nfolds')     = nfolds
             attr(prob, 'rep.times')  = rep.times
           }

         }
  )

  attr(prob, 'model.type') = model.type
  attr(prob, 'pred.at')    = pred.at
  attr(prob, 'ngroup')     = ngroup
  attr(prob, 'risk.group') = grp
  attr(prob, 'surv.time')  = time
  attr(prob, 'surv.event') = event
  attr(prob, 'seed')       = seed

  prob

}

#' Compute glmnet Predicted Survival Probabilities for Calibration
#'
#' @importFrom glmnet glmnet
#' @importFrom stats predict
#'
#' @return list containing predicted survival probability
#'
#' @keywords internal
glmnet.calibrate.internal.pred = function(x_tr, x_te, y_tr,
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

  lp = as.numeric(predict(object, newx = data.matrix(x_tr), s = lambda, type = 'link'))
  lpnew = as.numeric(predict(object, newx = data.matrix(x_te), s = lambda, type = 'link'))

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

#' Compute ncvreg Predicted Survival Probabilities for Calibration
#'
#' @importFrom ncvreg ncvsurv
#' @importFrom stats predict
#'
#' @return list containing predicted survival probability
#'
#' @keywords internal
ncvreg.calibrate.internal.pred = function(x_tr, x_te, y_tr,
                                          model.type,
                                          alpha, lambda, gamma,
                                          pred.at) {

  if (model.type == 'mcp') {
    object = ncvsurv(X = x_tr, y = y_tr, model = 'cox',
                     penalty = 'MCP', gamma = gamma,
                     alpha = 1, lambda = lambda)
  }

  if (model.type == 'mnet') {
    object = ncvsurv(X = x_tr, y = y_tr, model = 'cox',
                     penalty = 'MCP', gamma = gamma,
                     alpha = alpha, lambda = lambda)
  }

  if (model.type == 'scad') {
    object = ncvsurv(X = x_tr, y = y_tr, model = 'cox',
                     penalty = 'SCAD', gamma = gamma,
                     alpha = 1, lambda = lambda)
  }

  if (model.type == 'snet') {
    object = ncvsurv(X = x_tr, y = y_tr, model = 'cox',
                     penalty = 'SCAD', gamma = gamma,
                     alpha = alpha, lambda = lambda)
  }

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

#' Compute "penalized" Predicted Survival Probabilities for Calibration
#'
#' @importFrom penalized penalized
#' @importFrom stats predict
#'
#' @return list containing predicted survival probability
#'
#' @keywords internal
penalized.calibrate.internal.pred = function(x_tr, x_te, y_tr,
                                             lambda,
                                             pred.at) {

  object = penalized(response = y_tr, penalized = x_tr,
                     lambda1 = lambda, lambda2 = 0,
                     fusedl = TRUE, standardize = TRUE, model = 'cox')

  lp = as.vector(data.matrix(x_tr) %*% as.matrix(object@penalized))
  lpnew = as.vector(data.matrix(x_te) %*% as.matrix(object@penalized))

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

#' Compute Kaplan-Meier Estimated Survival Probabilities for Calibration
#'
#' @importFrom survival survfit
#' @importFrom survival Surv
#'
#' @return list
#'
#' @keywords internal
hdnom.calibrate.internal.true = function(pred_prob, grp,
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

#' Print Calibration Results
#'
#' Print Calibration Results
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
print.hdnom.calibrate = function(x, ...) {

  if (!('hdnom.calibrate' %in% class(x)))
    stop('object class must be "hdnom.calibrate"')

  method = setdiff(class(x), 'hdnom.calibrate')

  switch(method,

         glmnet.calibrate.fitting = {
           cat('High-Dimensional Cox Model Calibration Object\n')
           cat('Random seed:', attr(x, 'seed'), '\n')
           cat('Calibration method: fitting\n')
           cat('Model type:', attr(x, 'model.type'), '\n')
           cat('glmnet model alpha:', attr(x, 'alpha'), '\n')
           cat('glmnet model lambda:', attr(x, 'lambda'), '\n')
           if (is.null(attr(x, 'pen.factor'))) {
             cat('glmnet model penalty factor: not specified\n')
           } else {
             cat('glmnet model penalty factor: specified\n')
           }
           cat('Calibration time point:', attr(x, 'pred.at'), '\n')
           cat('Number of groups formed for calibration:', attr(x, 'ngroup'), '\n')
         },

         glmnet.calibrate.bootstrap = {
           cat('High-Dimensional Cox Model Calibration Object\n')
           cat('Random seed:', attr(x, 'seed'), '\n')
           cat('Calibration method: bootstrap\n')
           cat('Bootstrap samples:', attr(x, 'boot.times'), '\n')
           cat('Model type:', attr(x, 'model.type'), '\n')
           cat('glmnet model alpha:', attr(x, 'alpha'), '\n')
           cat('glmnet model lambda:', attr(x, 'lambda'), '\n')
           if (is.null(attr(x, 'pen.factor'))) {
             cat('glmnet model penalty factor: not specified\n')
           } else {
             cat('glmnet model penalty factor: specified\n')
           }
           cat('Calibration time point:', attr(x, 'pred.at'), '\n')
           cat('Number of groups formed for calibration:', attr(x, 'ngroup'), '\n')
         },

         glmnet.calibrate.cv = {
           cat('High-Dimensional Cox Model Calibration Object\n')
           cat('Random seed:', attr(x, 'seed'), '\n')
           cat('Calibration method: k-fold cross-validation\n')
           cat('Cross-validation folds:', attr(x, 'nfolds'), '\n')
           cat('Model type:', attr(x, 'model.type'), '\n')
           cat('glmnet model alpha:', attr(x, 'alpha'), '\n')
           cat('glmnet model lambda:', attr(x, 'lambda'), '\n')
           if (is.null(attr(x, 'pen.factor'))) {
             cat('glmnet model penalty factor: not specified\n')
           } else {
             cat('glmnet model penalty factor: specified\n')
           }
           cat('Calibration time point:', attr(x, 'pred.at'), '\n')
           cat('Number of groups formed for calibration:', attr(x, 'ngroup'), '\n')
         },

         glmnet.calibrate.repeated.cv = {
           cat('High-Dimensional Cox Model Calibration Object\n')
           cat('Random seed:', attr(x, 'seed'), '\n')
           cat('Calibration method: repeated cross-validation\n')
           cat('Cross-validation folds:', attr(x, 'nfolds'), '\n')
           cat('Cross-validation repeated times:', attr(x, 'rep.times'), '\n')
           cat('Model type:', attr(x, 'model.type'), '\n')
           cat('glmnet model alpha:', attr(x, 'alpha'), '\n')
           cat('glmnet model lambda:', attr(x, 'lambda'), '\n')
           if (is.null(attr(x, 'pen.factor'))) {
             cat('glmnet model penalty factor: not specified\n')
           } else {
             cat('glmnet model penalty factor: specified\n')
           }
           cat('Calibration time point:', attr(x, 'pred.at'), '\n')
           cat('Number of groups formed for calibration:', attr(x, 'ngroup'), '\n')
         },

         ncvreg.calibrate.fitting = {
           cat('High-Dimensional Cox Model Calibration Object\n')
           cat('Random seed:', attr(x, 'seed'), '\n')
           cat('Calibration method: fitting\n')
           cat('Model type:', attr(x, 'model.type'), '\n')
           cat('ncvreg model gamma:', attr(x, 'gamma'), '\n')
           cat('ncvreg model alpha:', attr(x, 'alpha'), '\n')
           cat('ncvreg model lambda:', attr(x, 'lambda'), '\n')
           cat('Calibration time point:', attr(x, 'pred.at'), '\n')
           cat('Number of groups formed for calibration:', attr(x, 'ngroup'), '\n')
         },

         ncvreg.calibrate.bootstrap = {
           cat('High-Dimensional Cox Model Calibration Object\n')
           cat('Random seed:', attr(x, 'seed'), '\n')
           cat('Calibration method: bootstrap\n')
           cat('Bootstrap samples:', attr(x, 'boot.times'), '\n')
           cat('Model type:', attr(x, 'model.type'), '\n')
           cat('ncvreg model gamma:', attr(x, 'gamma'), '\n')
           cat('ncvreg model alpha:', attr(x, 'alpha'), '\n')
           cat('ncvreg model lambda:', attr(x, 'lambda'), '\n')
           cat('Calibration time point:', attr(x, 'pred.at'), '\n')
           cat('Number of groups formed for calibration:', attr(x, 'ngroup'), '\n')
         },

         ncvreg.calibrate.cv = {
           cat('High-Dimensional Cox Model Calibration Object\n')
           cat('Random seed:', attr(x, 'seed'), '\n')
           cat('Calibration method: k-fold cross-validation\n')
           cat('Cross-validation folds:', attr(x, 'nfolds'), '\n')
           cat('Model type:', attr(x, 'model.type'), '\n')
           cat('ncvreg model gamma:', attr(x, 'gamma'), '\n')
           cat('ncvreg model alpha:', attr(x, 'alpha'), '\n')
           cat('ncvreg model lambda:', attr(x, 'lambda'), '\n')
           cat('Calibration time point:', attr(x, 'pred.at'), '\n')
           cat('Number of groups formed for calibration:', attr(x, 'ngroup'), '\n')
         },

         ncvreg.calibrate.repeated.cv = {
           cat('High-Dimensional Cox Model Calibration Object\n')
           cat('Random seed:', attr(x, 'seed'), '\n')
           cat('Calibration method: repeated cross-validation\n')
           cat('Cross-validation folds:', attr(x, 'nfolds'), '\n')
           cat('Cross-validation repeated times:', attr(x, 'rep.times'), '\n')
           cat('Model type:', attr(x, 'model.type'), '\n')
           cat('ncvreg model alpha:', attr(x, 'gamma'), '\n')
           cat('ncvreg model alpha:', attr(x, 'alpha'), '\n')
           cat('ncvreg model lambda:', attr(x, 'lambda'), '\n')
           cat('Calibration time point:', attr(x, 'pred.at'), '\n')
           cat('Number of groups formed for calibration:', attr(x, 'ngroup'), '\n')
         },

         penalized.calibrate.fitting = {
           cat('High-Dimensional Cox Model Calibration Object\n')
           cat('Random seed:', attr(x, 'seed'), '\n')
           cat('Calibration method: fitting\n')
           cat('Model type:', attr(x, 'model.type'), '\n')
           cat('Fused lasso model lambda:', attr(x, 'lambda'), '\n')
           cat('Calibration time point:', attr(x, 'pred.at'), '\n')
           cat('Number of groups formed for calibration:', attr(x, 'ngroup'), '\n')
         },

         penalized.calibrate.bootstrap = {
           cat('High-Dimensional Cox Model Calibration Object\n')
           cat('Random seed:', attr(x, 'seed'), '\n')
           cat('Calibration method: bootstrap\n')
           cat('Bootstrap samples:', attr(x, 'boot.times'), '\n')
           cat('Model type:', attr(x, 'model.type'), '\n')
           cat('Fused lasso model lambda:', attr(x, 'lambda'), '\n')
           cat('Calibration time point:', attr(x, 'pred.at'), '\n')
           cat('Number of groups formed for calibration:', attr(x, 'ngroup'), '\n')
         },

         penalized.calibrate.cv = {
           cat('High-Dimensional Cox Model Calibration Object\n')
           cat('Random seed:', attr(x, 'seed'), '\n')
           cat('Calibration method: k-fold cross-validation\n')
           cat('Cross-validation folds:', attr(x, 'nfolds'), '\n')
           cat('Model type:', attr(x, 'model.type'), '\n')
           cat('Fused lasso model lambda:', attr(x, 'lambda'), '\n')
           cat('Calibration time point:', attr(x, 'pred.at'), '\n')
           cat('Number of groups formed for calibration:', attr(x, 'ngroup'), '\n')
         },

         penalized.calibrate.repeated.cv = {
           cat('High-Dimensional Cox Model Calibration Object\n')
           cat('Random seed:', attr(x, 'seed'), '\n')
           cat('Calibration method: repeated cross-validation\n')
           cat('Cross-validation folds:', attr(x, 'nfolds'), '\n')
           cat('Cross-validation repeated times:', attr(x, 'rep.times'), '\n')
           cat('Model type:', attr(x, 'model.type'), '\n')
           cat('Fused lasso model lambda:', attr(x, 'lambda'), '\n')
           cat('Calibration time point:', attr(x, 'pred.at'), '\n')
           cat('Number of groups formed for calibration:', attr(x, 'ngroup'), '\n')
         }

  )

}

#' Summary of Calibration Results
#'
#' Summary of Calibration Results
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

  switch(method,

         glmnet.calibrate.fitting = {
           attr(object, 'alpha')      = NULL
           attr(object, 'lambda')     = NULL
           attr(object, 'pen.factor') = NULL
         },

         glmnet.calibrate.bootstrap = {
           attr(object, 'boot.times') = NULL
           attr(object, 'alpha')      = NULL
           attr(object, 'lambda')     = NULL
           attr(object, 'pen.factor') = NULL

         },

         glmnet.calibrate.cv = {
           attr(object, 'nfolds')     = NULL
           attr(object, 'alpha')      = NULL
           attr(object, 'lambda')     = NULL
           attr(object, 'pen.factor') = NULL
         },

         glmnet.calibrate.repeated.cv = {
           attr(object, 'nfolds')     = NULL
           attr(object, 'rep.times')  = NULL
           attr(object, 'alpha')      = NULL
           attr(object, 'lambda')     = NULL
           attr(object, 'pen.factor') = NULL
         },

         ncvreg.calibrate.fitting = {
           attr(object, 'alpha')      = NULL
           attr(object, 'lambda')     = NULL
           attr(object, 'gamma')      = NULL
         },

         ncvreg.calibrate.bootstrap = {
           attr(object, 'alpha')      = NULL
           attr(object, 'lambda')     = NULL
           attr(object, 'gamma')      = NULL
           attr(object, 'boot.times') = NULL
         },

         ncvreg.calibrate.cv = {
           attr(object, 'alpha')      = NULL
           attr(object, 'lambda')     = NULL
           attr(object, 'gamma')      = NULL
           attr(object, 'nfolds')     = NULL
         },

         ncvreg.calibrate.repeated.cv = {
           attr(object, 'alpha')      = NULL
           attr(object, 'lambda')     = NULL
           attr(object, 'gamma')      = NULL
           attr(object, 'nfolds')     = NULL
           attr(object, 'rep.times')  = NULL
         },

         penalized.calibrate.fitting = {
           attr(object, 'lambda')     = NULL
         },

         penalized.calibrate.bootstrap = {
           attr(object, 'lambda')     = NULL
           attr(object, 'boot.times') = NULL
         },

         penalized.calibrate.cv = {
           attr(object, 'lambda')     = NULL
           attr(object, 'nfolds')     = NULL
         },

         penalized.calibrate.repeated.cv = {
           attr(object, 'lambda')     = NULL
           attr(object, 'nfolds')     = NULL
           attr(object, 'rep.times')  = NULL
         }

  )

  attr(object, 'model.type') = NULL
  attr(object, 'pred.at')    = NULL
  attr(object, 'ngroup')     = NULL
  attr(object, 'risk.group') = NULL
  attr(object, 'surv.time')  = NULL
  attr(object, 'surv.event') = NULL
  attr(object, 'seed')       = NULL

  cat('  Calibration Summary Table\n')
  class(object) = 'matrix'
  print(object)

}

#' Plot Calibration Results
#'
#' Plot Calibration Results
#'
#' @param x an object returned by \code{\link{hdnom.calibrate}}.
#' @param xlim x axis limits of the plot.
#' @param ylim y axis limits of the plot.
#' @param col.pal Color palette to use. Possible values are
#' \code{"JCO"}, \code{"Lancet"}, \code{"NPG"}, and \code{"AAAS"}.
#' Default is \code{"JCO"}.
#' @param ... other parameters for \code{plot}.
#'
#' @method plot hdnom.calibrate
#'
#' @export
#'
#' @importFrom ggplot2 ggplot aes_string geom_errorbar
#' geom_line geom_point geom_abline xlab ylab theme_bw
#'
#' @examples
#' NULL
plot.hdnom.calibrate =
  function(x, xlim = c(0, 1), ylim = c(0, 1),
           col.pal = c('JCO', 'Lancet', 'NPG', 'AAAS'), ...) {

    if (!('hdnom.calibrate' %in% class(x)))
      stop('object class must be "hdnom.calibrate"')

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
