#' Validate Cox models fitted with glmnet using time-dependent AUC
#'
#' Validate Cox models fitted with glmnet using time-dependent AUC
#'
#' @param x Matrix of training data used for the \code{glmnet} object;
#' on which to run the validation.
#' @param time Survival time.
#' Must be of the same length with the number of rows as \code{x}.
#' @param event Status indicator, normally 0 = alive, 1 = dead.
#' Must be of the same length with the number of rows as \code{x}.
#' @param alpha Value of the elastic-net mixing parameter alpha in
#' glmnet. \code{alpha=1}: lasso; \code{alpha=0}: ridge.
#' From the Cox model you have built.
#' @param lambda Value of the penalty parameter lambda to use in the
#' glmnet fits on the resampled data. From the Cox model you have built.
#' @param method Validation method.
#' Could be \code{"bootstrap"}, \code{"cv"}, or \code{"repeated.cv"}.
#' @param boot.times Number of repetitions for bootstrap.
#' @param nfolds Number of folds for cross-validation and
#' repeated cross-validation.
#' @param rep.times Number of repeated times for repeated cross-validation.
#' @param tauc.type Type of time-dependent AUC.
#' Including \code{"CD"} proposed by Chambless and Diao (2006).,
#' \code{"SZ"} proposed by Song and Zhou (2008).,
#' \code{"UNO"} proposed by Uno et al. (2007).
#' @param tauc.time Numeric vector. Time points at which to evaluate
#' the time-dependent AUC.
#' @param trace Logical. Print trace or not. Default is \code{TRUE}.
#'
#' @export glmnet.validate
#'
#' @references
#' Chambless, L. E. and G. Diao (2006).
#' Estimation of time-dependent area under the ROC curve for long-term
#' risk prediction.
#' \emph{Statistics in Medicine} 25, 3474--3486.
#'
#' Song, X. and X.-H. Zhou (2008).
#' A semiparametric approach for the covariate specific ROC curve with
#' survival outcome.
#' \emph{Statistica Sinica} 18, 947--965.
#'
#' Uno, H., T. Cai, L. Tian, and L. J. Wei (2007).
#' Evaluating prediction rules for t-year survivors with censored
#' regression models.
#' \emph{Journal of the American Statistical Association} 102, 527--537.
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
#' cvfit = cv.glmnet(x, Surv(time, event), family = "cox", nfolds = 10)
#' fit = glmnet(x, Surv(time, event), lambda = cvfit$lambda.1se, family = "cox")
#'
#' # Model validation by bootstrap with time-dependent AUC
#' val.boot = glmnet.validate(x, time, event,
#'                            alpha = 1, lambda = cvfit$lambda.1se,
#'                            method = "bootstrap", boot.times = 20,
#'                            tauc.type = "UNO", tauc.time = seq(0.25, 2, 0.25) * 365)
#'
#' # Model validation by 10-fold cross-validation with time-dependent AUC
#' val.cv = glmnet.validate(x, time, event,
#'                          alpha = 1, lambda = cvfit$lambda.1se,
#'                          method = "cv", nfolds = 10,
#'                          tauc.type = "UNO", tauc.time = seq(0.25, 2, 0.25) * 365)
#'
#' # Model validation by repeated cross-validation with time-dependent AUC
#' val.repcv = glmnet.validate(x, time, event,
#'                             alpha = 1, lambda = cvfit$lambda.1se,
#'                             method = "repeated.cv", nfolds = 10, rep.times = 5,
#'                             tauc.type = "UNO", tauc.time = seq(0.25, 2, 0.25) * 365)
#'
#' # bootstrap-based discrimination curves has a very narrow band
#' print(val.boot)
#' summary(val.boot)
#' plot(val.boot, ylim = c(0.4, 0.8))
#'
#' # k-fold cv provides a more strict evaluation than bootstrap
#' print(val.cv)
#' summary(val.cv)
#' plot(val.cv, ylim = c(0.4, 0.8))
#'
#' # repeated cv provides similar results as k-fold cv
#' # but more stable than k-fold cv
#' print(val.repcv)
#' summary(val.repcv)
#' plot(val.repcv, ylim = c(0.4, 0.8))
glmnet.validate = function (x, time, event,
                            alpha, lambda,
                            method = c('bootstrap', 'cv', 'repeated.cv'),
                            boot.times = NULL, nfolds = NULL, rep.times = NULL,
                            tauc.type = c("CD", "SZ", "UNO"), tauc.time,
                            trace = TRUE) {

  method = match.arg(method)
  tauc.type = match.arg(tauc.type)

  if (method == 'bootstrap') {

    if (!is.null(nfolds) || !is.null(rep.times))
      stop('nfolds and rep.times must be NULL when method = "bootstrap"')

    if (is.null(boot.times)) stop('please specify boot.times')

    # generate bootstrap sample index
    samp_mat = matrix(NA, nrow(x), boot.times)
    for(i in 1L:boot.times) samp_mat[, i] = sample(1L:nrow(x), replace = TRUE)

    # bootstrap validation main loop
    tauc = vector('list', boot.times)
    for (i in 1L:ncol(samp_mat)) {

      if (trace) cat('Start bootstrap sample', i, '\n')

      samp_idx = samp_mat[, i]
      x_tr = x[samp_idx, ]
      time_tr = time[samp_idx]
      event_tr = event[samp_idx]
      y_tr = Surv(time_tr, event_tr)
      x_te = x  # use original dataset as test set
      y_te = Surv(time, event)  # use original dataset as test set

      tauc[[i]] =
        glmnet.validate.internal(
          x_tr = x_tr, x_te = x_te, y_tr = y_tr, y_te = y_te,
          alpha = alpha, lambda = lambda,
          tauc.type = tauc.type, tauc.time = tauc.time
        )

    }

  } else if (method == 'cv') {

    if (!is.null(boot.times) || !is.null(rep.times))
      stop('boot.times and rep.times must be NULL when method = "cv"')

    if (is.null(nfolds)) stop('please specify nfolds')

    # generate cross-validation sample index
    row_x = nrow(x)
    samp_idx = sample(rep_len(1L:nfolds, row_x))

    # cross-validation main loop
    tauc = vector('list', nfolds)
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

      tauc[[i]] =
        glmnet.validate.internal(
          x_tr = x_tr, x_te = x_te, y_tr = y_tr, y_te = y_te,
          alpha = alpha, lambda = lambda,
          tauc.type = tauc.type, tauc.time = tauc.time
        )

    }

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
    tauc = vector('list', rep.times)
    for (k in 1L:rep.times) tauc[[k]] = vector('list', nfolds)

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

        tauc[[j]][[i]] =
          glmnet.validate.internal(
            x_tr = x_tr, x_te = x_te, y_tr = y_tr, y_te = y_te,
            alpha = alpha, lambda = lambda,
            tauc.type = tauc.type, tauc.time = tauc.time
          )

      }
    }

  } else {
    stop('method must be one of "bootstrap", "cv", or "repeated.cv"')
  }

  switch(method,
         bootstrap = {
           class(tauc) = c('glmnet.validate',
                           'glmnet.validate.bootstrap')
           attr(tauc, 'boot.times') = boot.times
           attr(tauc, 'tauc.type')  = tauc.type
           attr(tauc, 'tauc.time')  = tauc.time
         },
         cv = {
           class(tauc) = c('glmnet.validate',
                           'glmnet.validate.cv')
           attr(tauc, 'nfolds')    = nfolds
           attr(tauc, 'tauc.type') = tauc.type
           attr(tauc, 'tauc.time') = tauc.time
         },
         repeated.cv = {
           class(tauc) = c('glmnet.validate',
                           'glmnet.validate.repeated.cv')
           attr(tauc, 'nfolds')    = nfolds
           attr(tauc, 'rep.times') = rep.times
           attr(tauc, 'tauc.type') = tauc.type
           attr(tauc, 'tauc.time') = tauc.time
         }
  )

  tauc

}

#' Compute validation measures for glmnet objects by bootstrap
#'
#' @importFrom survAUC AUC.cd AUC.sh AUC.uno
#' @importFrom glmnet glmnet
#' @importFrom survival Surv
#'
#' @return time-dependent AUC (tAUC) value
#'
#' @keywords internal
glmnet.validate.internal = function(x_tr, x_te, y_tr, y_te,
                                    alpha, lambda,
                                    tauc.type, tauc.time) {

  samp_fit = glmnet(x = x_tr, y = y_tr, family = 'cox',
                    alpha = alpha, lambda = lambda)

  lp_tr = as.vector(predict(samp_fit, newx = x_tr, type = 'link'))
  lp_te = as.vector(predict(samp_fit, newx = x_te, type = 'link'))

  tauc_list = switch(tauc.type,
                     CD = {
                       AUC.cd(Surv.rsp = y_tr, Surv.rsp.new = y_te,
                              lp = lp_tr, lpnew = lp_te,
                              times = tauc.time)
                     },
                     SZ = {
                       AUC.sh(Surv.rsp = y_tr, Surv.rsp.new = y_te,
                              lp = lp_tr, lpnew = lp_te,
                              times = tauc.time)
                     },
                     UNO = {
                       AUC.uno(Surv.rsp = y_tr, Surv.rsp.new = y_te,
                               lpnew = lp_te,
                               times = tauc.time)
                     }
  )

  tauc_list

}

#' Print validation result generated with glmnet.validate
#'
#' Print validation result generated with glmnet.validate
#'
#' @param x a \code{"glmnet.validate"} object.
#' @param ... other parameters (not used).
#'
#' @method print glmnet.validate
#'
#' @export
#'
#' @examples
#' NULL
print.glmnet.validate = function (x, ...) {

  if (!('glmnet.validate' %in% class(x)))
    stop('object class must be "glmnet.validate"')

  method = setdiff(class(x), 'glmnet.validate')

  if (method == 'glmnet.validate.bootstrap') {

    cat('High-dimensional Cox model validation object\n')
    cat('Validation method: bootstrap\n')
    cat('Bootstrap samples:', attr(x, 'boot.times'), '\n')
    cat('Time-dependent AUC type:', attr(x, 'tauc.type'), '\n')
    cat('Evaluation time points for tAUC:', attr(x, 'tauc.time'))

  } else if (method == 'glmnet.validate.cv') {

    cat('High-dimensional Cox model validation object\n')
    cat('Validation method: k-fold cross-validation\n')
    cat('Cross-validation folds:', attr(x, 'nfolds'), '\n')
    cat('Time-dependent AUC type:', attr(x, 'tauc.type'), '\n')
    cat('Evaluation time points for tAUC:', attr(x, 'tauc.time'))

  } else if (method == 'glmnet.validate.repeated.cv') {

    cat('High-dimensional Cox model validation object\n')
    cat('Validation method: repeated cross-validation\n')
    cat('Cross-validation folds:', attr(x, 'nfolds'), '\n')
    cat('Cross-validation repeated times:', attr(x, 'rep.times'), '\n')
    cat('Time-dependent AUC type:', attr(x, 'tauc.type'), '\n')
    cat('Evaluation time points for tAUC:', attr(x, 'tauc.time'))

  } else {

    stop('glmnet.validate object is not valid')

  }

}

#' Summary validation result generated with glmnet.validate
#'
#' Summary validation result generated with glmnet.validate
#'
#' @param object a \code{"glmnet.validate"} object.
#' @param silent Print summary table header or not.
#' Default is \code{FALSE}.
#' @param ... other parameters (not used).
#'
#' @method summary glmnet.validate
#'
#' @export
#'
#' @examples
#' NULL
summary.glmnet.validate = function (object, silent = FALSE, ...) {

  if (!('glmnet.validate' %in% class(object)))
    stop('object class must be "glmnet.validate"')

  method = setdiff(class(object), 'glmnet.validate')

  if (method == 'glmnet.validate.bootstrap') {

    boot.times = attr(object, 'boot.times')
    tauc.time = attr(object, 'tauc.time')
    aucmat = matrix(NA, ncol = length(tauc.time), nrow = boot.times)
    for (i in 1L:boot.times) aucmat[i, ] = object[[i]]$auc
    summary_mat = rbind(apply(aucmat, 2, mean), apply(aucmat, 2, quantile))
    rownames(summary_mat) = c('Mean', 'Min', '0.25 Qt.',
                              'Median', '0.75 Qt.', 'Max')
    colnames(summary_mat) = tauc.time

  } else if (method == 'glmnet.validate.cv') {

    nfolds = attr(object, 'nfolds')
    tauc.time = attr(object, 'tauc.time')
    aucmat = matrix(NA, ncol = length(tauc.time), nrow = nfolds)
    for (i in 1L:nfolds) aucmat[i, ] = object[[i]]$auc
    summary_mat = rbind(apply(aucmat, 2, mean), apply(aucmat, 2, quantile))
    rownames(summary_mat) = c('Mean', 'Min', '0.25 Qt.',
                              'Median', '0.75 Qt.', 'Max')
    colnames(summary_mat) = tauc.time

  } else if (method == 'glmnet.validate.repeated.cv') {

    nfolds = attr(object, 'nfolds')
    rep.times = attr(object, 'rep.times')
    tauc.time = attr(object, 'tauc.time')
    auclist = vector('list', rep.times)
    for (i in 1L:rep.times) {
      auclist[[i]] = matrix(NA, ncol = length(tauc.time), nrow = nfolds)
    }
    for (i in 1L:rep.times) {
      for (j in 1L:nfolds) {
        auclist[[i]][j, ] = object[[i]][[j]]$auc
      }
    }

    summary_list = vector('list', rep.times)
    for (i in 1L:rep.times) {
      summary_list[[i]] = rbind(apply(auclist[[i]], 2, mean),
                                apply(auclist[[i]], 2, quantile))
    }

    summary_mat = Reduce('+', summary_list)/length(summary_list)
    rownames(summary_mat) = c('Mean of Mean', 'Mean of Min',
                              'Mean of 0.25 Qt.', 'Mean of Median',
                              'Mean of 0.75 Qt.', 'Mean of Max')
    colnames(summary_mat) = tauc.time

    if (!silent) {
      cat('Note: for repeated CV, we evaluated quantile statistic tables for\n')
      cat('each CV repeat, then calculated element-wise mean across all tables.\n')
    }

  } else {
    stop('glmnet.validate object is not valid')
  }

  if (!silent)
    cat('Time-dependent AUC summary over all bootstrap runs at evaluation time points\n')

  return(summary_mat)

}

#' Plot optimism-corrected time-dependent discrimination curves
#'
#' Plot optimism-corrected time-dependent discrimination curves
#'
#' @param x a \code{"glmnet.validate"} object.
#' @param ylim y axis limits of the plot.
#' @param xlab title for the x axis.
#' @param ylab title for the y axis.
#' @param ... other parameters for \code{plot}.
#'
#' @method plot glmnet.validate
#'
#' @export
#'
#' @examples
#' NULL
plot.glmnet.validate = function (x, ylim = c(0.5, 1), xlab = 'Time',
                                 ylab = 'Time-dependent Area under ROC',
                                 ...) {

  if (!('glmnet.validate' %in% class(x)))
    stop('object class must be "glmnet.validate"')

  mat = summary(x, silent = TRUE)
  tauc_time = attr(x, 'tauc.time')
  tauc_mean = mat[1L, ]
  tauc_median = mat[4L, ]
  tauc_q25 = mat[3L, ]
  tauc_q75 = mat[5L, ]

  # two panels, one for plot, one for legend
  layout(rbind(1, 2), heights = c(7, 1))

  plot(tauc_time, tauc_median, type = 'l',
       xlab = xlab, ylab = ylab, ylim = ylim)
  polygon(c(tauc_time, rev(tauc_time)), c(tauc_q25, rev(tauc_q75)),
          col = 'grey85', border = FALSE)
  lines(tauc_time, tauc_median, lty = 5, lwd = 2)
  lines(tauc_time, tauc_mean, lty = 1, lwd = 2)
  lines(tauc_time, tauc_q75, lty = 3, lwd = 1)
  lines(tauc_time, tauc_q25, lty = 3, lwd = 1)

  for (i in seq(0.1, 1, 0.1)) abline(h = i, lty = 3)

  mar = par('mar')
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend('center', 'groups',
         c('Median', 'Mean', '25th/75th Quantiles'),
         lwd = c(2, 2, 1), lty = c(5, 1, 3),
         bty = 'n', horiz = TRUE, text.width = c(0.1, 0.1, 0.1))
  par(mar = mar)

}
