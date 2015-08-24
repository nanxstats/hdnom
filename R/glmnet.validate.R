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
#' NULL
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
           class(tauc) = c('glmnet.validate', 'glmnet.validate.bootstrap')
         },
         cv = {
           class(tauc) = c('glmnet.validate', 'glmnet.validate.cv')
         },
         repeated.cv = {
           class(tauc) = c('glmnet.validate', 'glmnet.validate.repeated.cv')
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

#' Plot validation result generated with glmnet.validate
#'
#' Plot validation result generated with glmnet.validate
#'
#' @param x a \code{"glmnet.validate"} object.
#' @param ... other parameters for \code{plot}.
#'
#' @method plot glmnet.validate
#'
#' @export
#'
#' @examples
#' NULL
plot.glmnet.validate = function (x, ...) {

  if (class(x) != 'glmnet.validate')
    stop('object class must be "glmnet.validate"')

  val = x

  plot(val, ...)

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

  if (class(x) != 'glmnet.validate')
    stop('object class must be "glmnet.validate"')

  val = x

  print(val)

}
