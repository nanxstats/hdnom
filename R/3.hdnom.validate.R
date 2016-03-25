#' Validate High-Dimensional Cox Models with Time-Dependent AUC
#'
#' Validate High-Dimensional Cox Models with Time-Dependent AUC
#'
#' @param x Matrix of training data used for fitting the model;
#' on which to run the validation.
#' @param time Survival time.
#' Must be of the same length with the number of rows as \code{x}.
#' @param event Status indicator, normally 0 = alive, 1 = dead.
#' Must be of the same length with the number of rows as \code{x}.
#' @param model.type Model type to validate. Could be one of \code{"lasso"},
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
#' model fits on the resampled data. From the built Cox model.
#' @param pen.factor Penalty factors to apply to each coefficient.
#' From the built \emph{adaptive lasso} or \emph{adaptive elastic-net} model.
#' @param gamma Value of the model parameter gamma for
#' MCP/SCAD/Mnet/Snet models.
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
#' @param seed A random seed for resampling.
#' @param trace Logical. Output the validation progress or not.
#' Default is \code{TRUE}.
#'
#' @export hdnom.validate
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
#'
#' # Load imputed SMART data
#' data(smart)
#' x = as.matrix(smart[, -c(1, 2)])[1:500, ]
#' time = smart$TEVENT[1:500]
#' event = smart$EVENT[1:500]
#'
#' # Fit penalized Cox model (lasso penalty) with glmnet
#' set.seed(1010)
#' cvfit = cv.glmnet(x, Surv(time, event), family = "cox", nfolds = 5)
#'
#' # Model validation by bootstrap with time-dependent AUC
#' # Normally boot.times should be set to 200 or more,
#' # we set it to 3 here only to save example running time.
#' val.boot = hdnom.validate(x, time, event, model.type = "lasso",
#'                           alpha = 1, lambda = cvfit$lambda.1se,
#'                           method = "bootstrap", boot.times = 3,
#'                           tauc.type = "UNO", tauc.time = seq(0.25, 2, 0.25) * 365,
#'                           seed = 1010)
#'
#' # Model validation by 5-fold cross-validation with time-dependent AUC
#' val.cv = hdnom.validate(x, time, event, model.type = "lasso",
#'                         alpha = 1, lambda = cvfit$lambda.1se,
#'                         method = "cv", nfolds = 5,
#'                         tauc.type = "UNO", tauc.time = seq(0.25, 2, 0.25) * 365,
#'                         seed = 1010)
#'
#' # Model validation by repeated cross-validation with time-dependent AUC
#' val.repcv = hdnom.validate(x, time, event, model.type = "lasso",
#'                            alpha = 1, lambda = cvfit$lambda.1se,
#'                            method = "repeated.cv", nfolds = 5, rep.times = 3,
#'                            tauc.type = "UNO", tauc.time = seq(0.25, 2, 0.25) * 365,
#'                            seed = 1010)
#'
#' # bootstrap-based discrimination curves has a very narrow band
#' print(val.boot)
#' summary(val.boot)
#' plot(val.boot)
#'
#' # k-fold cv provides a more strict evaluation than bootstrap
#' print(val.cv)
#' summary(val.cv)
#' plot(val.cv)
#'
#' # repeated cv provides similar results as k-fold cv
#' # but more robust than k-fold cv
#' print(val.repcv)
#' summary(val.repcv)
#' plot(val.repcv)
#'
# ### Testing fused lasso, SCAD, and Mnet models ###
# library("survival")
# library("rms")
#
# # Load imputed SMART data
# data(smart)
# x = as.matrix(smart[, -c(1, 2)])[1:500,]
# time = smart$TEVENT[1:500]
# event = smart$EVENT[1:500]
# y = Surv(time, event)
#
# set.seed(1010)
# val.boot = hdnom.validate(x, time, event, model.type = "flasso",
#                           lambda = 60,
#                           method = "bootstrap", boot.times = 10,
#                           tauc.type = "UNO", tauc.time = seq(0.25, 2, 0.25) * 365,
#                           seed = 1010)
#
# val.cv = hdnom.validate(x, time, event, model.type = "scad",
#                         gamma = 3.7, alpha = 1, lambda = 0.05,
#                         method = "cv", nfolds = 5,
#                         tauc.type = "UNO", tauc.time = seq(0.25, 2, 0.25) * 365,
#                         seed = 1010)
#
# val.repcv = hdnom.validate(x, time, event, model.type = "mnet",
#                            gamma = 3, alpha = 0.3, lambda = 0.05,
#                            method = "repeated.cv", nfolds = 5, rep.times = 3,
#                            tauc.type = "UNO", tauc.time = seq(0.25, 2, 0.25) * 365,
#                            seed = 1010)
#
# print(val.boot)
# summary(val.boot)
# plot(val.boot)
#
# print(val.cv)
# summary(val.cv)
# plot(val.cv)
#
# print(val.repcv)
# summary(val.repcv)
# plot(val.repcv)
hdnom.validate = function(x, time, event,
                          model.type = c('lasso', 'alasso', 'flasso',
                                         'enet', 'aenet',
                                         'mcp', 'mnet',
                                         'scad', 'snet'),
                          alpha, lambda, pen.factor = NULL, gamma,
                          method = c('bootstrap', 'cv', 'repeated.cv'),
                          boot.times = NULL, nfolds = NULL, rep.times = NULL,
                          tauc.type = c("CD", "SZ", "UNO"), tauc.time,
                          seed = 1001, trace = TRUE) {

  model.type = match.arg(model.type)
  method = match.arg(method)
  tauc.type = match.arg(tauc.type)

  set.seed(seed)

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
      x_tr = x[samp_idx, , drop = FALSE]
      time_tr = time[samp_idx]
      event_tr = event[samp_idx]
      y_tr = Surv(time_tr, event_tr)
      x_te = x  # use original dataset as test set
      y_te = Surv(time, event)  # use original dataset as test set

      if (model.type %in% c('lasso', 'alasso', 'enet', 'aenet')) {
        tauc[[i]] =
          glmnet.validate.internal(
            x_tr = x_tr, x_te = x_te, y_tr = y_tr, y_te = y_te,
            alpha = alpha, lambda = lambda, pen.factor = pen.factor,
            tauc.type = tauc.type, tauc.time = tauc.time
          )
      }

      if (model.type %in% c('mcp', 'mnet', 'scad', 'snet')) {
        tauc[[i]] =
          ncvreg.validate.internal(
            x_tr = x_tr, x_te = x_te, y_tr = y_tr, y_te = y_te,
            model.type = model.type,
            gamma = gamma, alpha = alpha, lambda = lambda,
            tauc.type = tauc.type, tauc.time = tauc.time
          )
      }

      if (model.type %in% c('flasso')) {
        tauc[[i]] =
          penalized.validate.internal(
            x_tr = x_tr, x_te = x_te, y_tr = y_tr, y_te = y_te,
            lambda = lambda,
            tauc.type = tauc.type, tauc.time = tauc.time
          )
      }

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

      x_tr = x[samp_idx != i, , drop = FALSE]
      time_tr = time[samp_idx != i]
      event_tr = event[samp_idx != i]
      y_tr = Surv(time_tr, event_tr)
      x_te  = x[samp_idx == i, , drop = FALSE]
      time_te = time[samp_idx == i]
      event_te = event[samp_idx == i]
      y_te = Surv(time_te, event_te)

      if (model.type %in% c('lasso', 'alasso', 'enet', 'aenet')) {
        tauc[[i]] =
          glmnet.validate.internal(
            x_tr = x_tr, x_te = x_te, y_tr = y_tr, y_te = y_te,
            alpha = alpha, lambda = lambda, pen.factor = pen.factor,
            tauc.type = tauc.type, tauc.time = tauc.time
          )
      }

      if (model.type %in% c('mcp', 'mnet', 'scad', 'snet')) {
        tauc[[i]] =
          ncvreg.validate.internal(
            x_tr = x_tr, x_te = x_te, y_tr = y_tr, y_te = y_te,
            model.type = model.type,
            gamma = gamma, alpha = alpha, lambda = lambda,
            tauc.type = tauc.type, tauc.time = tauc.time
          )
      }

      if (model.type %in% c('flasso')) {
        tauc[[i]] =
          penalized.validate.internal(
            x_tr = x_tr, x_te = x_te, y_tr = y_tr, y_te = y_te,
            lambda = lambda,
            tauc.type = tauc.type, tauc.time = tauc.time
          )
      }

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

        x_tr = x[samp_idx[[j]] != i, , drop = FALSE]
        time_tr = time[samp_idx[[j]] != i]
        event_tr = event[samp_idx[[j]] != i]
        y_tr = Surv(time_tr, event_tr)
        x_te  = x[samp_idx[[j]] == i, , drop = FALSE]
        time_te = time[samp_idx[[j]] == i]
        event_te = event[samp_idx[[j]] == i]
        y_te = Surv(time_te, event_te)

        if (model.type %in% c('lasso', 'alasso', 'enet', 'aenet')) {
          tauc[[j]][[i]] =
            glmnet.validate.internal(
              x_tr = x_tr, x_te = x_te, y_tr = y_tr, y_te = y_te,
              alpha = alpha, lambda = lambda, pen.factor = pen.factor,
              tauc.type = tauc.type, tauc.time = tauc.time
            )
        }

        if (model.type %in% c('mcp', 'mnet', 'scad', 'snet')) {
          tauc[[j]][[i]] =
            ncvreg.validate.internal(
              x_tr = x_tr, x_te = x_te, y_tr = y_tr, y_te = y_te,
              model.type = model.type,
              gamma = gamma, alpha = alpha, lambda = lambda,
              tauc.type = tauc.type, tauc.time = tauc.time
            )
        }

        if (model.type %in% c('flasso')) {
          tauc[[j]][[i]] =
            penalized.validate.internal(
              x_tr = x_tr, x_te = x_te, y_tr = y_tr, y_te = y_te,
              lambda = lambda,
              tauc.type = tauc.type, tauc.time = tauc.time
            )
        }

      }
    }

  } else {
    stop('method must be one of "bootstrap", cv", or "repeated.cv"')
  }

  switch(method,
         bootstrap = {

           if (model.type %in% c('lasso', 'alasso', 'enet', 'aenet')) {
             class(tauc) = c('hdnom.validate',
                             'glmnet.validate.bootstrap')
             attr(tauc, 'model.type') = model.type
             attr(tauc, 'alpha')      = alpha
             attr(tauc, 'lambda')     = lambda
             attr(tauc, 'pen.factor') = pen.factor
             attr(tauc, 'boot.times') = boot.times
             attr(tauc, 'tauc.type')  = tauc.type
             attr(tauc, 'tauc.time')  = tauc.time
             attr(tauc, 'seed')       = seed
           }

           if (model.type %in% c('mcp', 'mnet', 'scad', 'snet')) {
             class(tauc) = c('hdnom.validate',
                             'ncvreg.validate.bootstrap')
             attr(tauc, 'model.type') = model.type
             attr(tauc, 'gamma')      = gamma
             attr(tauc, 'alpha')      = alpha
             attr(tauc, 'lambda')     = lambda
             attr(tauc, 'boot.times') = boot.times
             attr(tauc, 'tauc.type')  = tauc.type
             attr(tauc, 'tauc.time')  = tauc.time
             attr(tauc, 'seed')       = seed
           }

           if (model.type %in% c('flasso')) {
             class(tauc) = c('hdnom.validate',
                             'penalized.validate.bootstrap')
             attr(tauc, 'model.type') = model.type
             attr(tauc, 'lambda')     = lambda
             attr(tauc, 'boot.times') = boot.times
             attr(tauc, 'tauc.type')  = tauc.type
             attr(tauc, 'tauc.time')  = tauc.time
             attr(tauc, 'seed')       = seed
           }

         },

         cv = {

           if (model.type %in% c('lasso', 'alasso', 'enet', 'aenet')) {
             class(tauc) = c('hdnom.validate',
                             'glmnet.validate.cv')
             attr(tauc, 'model.type') = model.type
             attr(tauc, 'alpha')      = alpha
             attr(tauc, 'lambda')     = lambda
             attr(tauc, 'pen.factor') = pen.factor
             attr(tauc, 'nfolds')     = nfolds
             attr(tauc, 'tauc.type')  = tauc.type
             attr(tauc, 'tauc.time')  = tauc.time
             attr(tauc, 'seed')       = seed
           }

           if (model.type %in% c('mcp', 'mnet', 'scad', 'snet')) {
             class(tauc) = c('hdnom.validate',
                             'ncvreg.validate.cv')
             attr(tauc, 'model.type') = model.type
             attr(tauc, 'gamma')      = gamma
             attr(tauc, 'alpha')      = alpha
             attr(tauc, 'lambda')     = lambda
             attr(tauc, 'nfolds')     = nfolds
             attr(tauc, 'tauc.type')  = tauc.type
             attr(tauc, 'tauc.time')  = tauc.time
             attr(tauc, 'seed')       = seed
           }

           if (model.type %in% c('flasso')) {
             class(tauc) = c('hdnom.validate',
                             'penalized.validate.cv')
             attr(tauc, 'model.type') = model.type
             attr(tauc, 'lambda')     = lambda
             attr(tauc, 'nfolds')     = nfolds
             attr(tauc, 'tauc.type')  = tauc.type
             attr(tauc, 'tauc.time')  = tauc.time
             attr(tauc, 'seed')       = seed
           }

         },

         repeated.cv = {

           if (model.type %in% c('lasso', 'alasso', 'enet', 'aenet')) {
             class(tauc) = c('hdnom.validate',
                             'glmnet.validate.repeated.cv')
             attr(tauc, 'model.type') = model.type
             attr(tauc, 'alpha')      = alpha
             attr(tauc, 'lambda')     = lambda
             attr(tauc, 'pen.factor') = pen.factor
             attr(tauc, 'nfolds')     = nfolds
             attr(tauc, 'rep.times')  = rep.times
             attr(tauc, 'tauc.type')  = tauc.type
             attr(tauc, 'tauc.time')  = tauc.time
             attr(tauc, 'seed')       = seed
           }

           if (model.type %in% c('mcp', 'mnet', 'scad', 'snet')) {
             class(tauc) = c('hdnom.validate',
                             'ncvreg.validate.repeated.cv')
             attr(tauc, 'model.type') = model.type
             attr(tauc, 'gamma')      = gamma
             attr(tauc, 'alpha')      = alpha
             attr(tauc, 'lambda')     = lambda
             attr(tauc, 'nfolds')     = nfolds
             attr(tauc, 'rep.times')  = rep.times
             attr(tauc, 'tauc.type')  = tauc.type
             attr(tauc, 'tauc.time')  = tauc.time
             attr(tauc, 'seed')       = seed
           }

           if (model.type %in% c('flasso')) {
             class(tauc) = c('hdnom.validate',
                             'penalized.validate.repeated.cv')
             attr(tauc, 'model.type') = model.type
             attr(tauc, 'lambda')     = lambda
             attr(tauc, 'nfolds')     = nfolds
             attr(tauc, 'rep.times')  = rep.times
             attr(tauc, 'tauc.type')  = tauc.type
             attr(tauc, 'tauc.time')  = tauc.time
             attr(tauc, 'seed')       = seed
           }

         }

  )

  tauc

}

#' Compute Validation Measures for glmnet Objects
#'
#' @importFrom survAUC AUC.cd AUC.sh AUC.uno
#' @importFrom glmnet glmnet
#' @importFrom survival Surv
#'
#' @return time-dependent AUC (tAUC) value
#'
#' @keywords internal
glmnet.validate.internal = function(x_tr, x_te, y_tr, y_te,
                                    alpha, lambda, pen.factor,
                                    tauc.type, tauc.time) {

  if (is.null(pen.factor)) {
    samp_fit = glmnet(x = x_tr, y = y_tr, family = 'cox',
                      alpha = alpha, lambda = lambda)
  } else {
    samp_fit = glmnet(x = x_tr, y = y_tr, family = 'cox',
                      alpha = alpha, lambda = lambda,
                      penalty.factor = pen.factor)
  }

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

#' Compute Validation Measures for ncvreg Model Objects
#'
#' @importFrom survAUC AUC.cd AUC.sh AUC.uno
#' @importFrom ncvreg ncvsurv
#' @importFrom survival Surv
#'
#' @return time-dependent AUC (tAUC) value
#'
#' @keywords internal
ncvreg.validate.internal = function(x_tr, x_te, y_tr, y_te, model.type,
                                    gamma, alpha, lambda,
                                    tauc.type, tauc.time) {

  if (model.type == 'mcp') {
    samp_fit = ncvsurv(X = x_tr, y = y_tr, model = 'cox',
                       penalty = 'MCP', gamma = gamma,
                       alpha = 1, lambda = lambda)
  }

  if (model.type == 'mnet') {
    samp_fit = ncvsurv(X = x_tr, y = y_tr, model = 'cox',
                       penalty = 'MCP', gamma = gamma,
                       alpha = alpha, lambda = lambda)
  }

  if (model.type == 'scad') {
    samp_fit = ncvsurv(X = x_tr, y = y_tr, model = 'cox',
                       penalty = 'SCAD', gamma = gamma,
                       alpha = 1, lambda = lambda)
  }

  if (model.type == 'snet') {
    samp_fit = ncvsurv(X = x_tr, y = y_tr, model = 'cox',
                       penalty = 'SCAD', gamma = gamma,
                       alpha = alpha, lambda = lambda)
  }

  lp_tr = as.vector(predict(samp_fit, X = x_tr, type = 'link'))
  lp_te = as.vector(predict(samp_fit, X = x_te, type = 'link'))

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

#' Compute Validation Measures for "penalized" Model Objects
#'
#' @importFrom survAUC AUC.cd AUC.sh AUC.uno
#' @importFrom penalized penalized
#' @importFrom survival Surv
#'
#' @return time-dependent AUC (tAUC) value
#'
#' @keywords internal
penalized.validate.internal = function(x_tr, x_te, y_tr, y_te,
                                       lambda,
                                       tauc.type, tauc.time) {

  samp_fit = penalized(response = y_tr, penalized = x_tr,
                       lambda1 = lambda, lambda2 = 0,
                       fusedl = TRUE, standardize = TRUE, model = 'cox')

  lp_tr = as.vector(samp_fit@lin.pred)
  lp_te = as.vector(x_te %*% as.matrix(samp_fit@penalized))

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

#' Print Validation Results
#'
#' Print Validation Results
#'
#' @param x An object returned by \code{\link{hdnom.validate}}.
#' @param ... Other parameters (not used).
#'
#' @method print hdnom.validate
#'
#' @export
#'
#' @examples
#' NULL
print.hdnom.validate = function(x, ...) {

  if (!('hdnom.validate' %in% class(x)))
    stop('object class must be "hdnom.validate"')

  method = setdiff(class(x), 'hdnom.validate')

  switch(method,

         glmnet.validate.bootstrap = {
           cat('High-Dimensional Cox Model Validation Object\n')
           cat('Random seed:', attr(x, 'seed'), '\n')
           cat('Validation method: bootstrap\n')
           cat('Bootstrap samples:', attr(x, 'boot.times'), '\n')
           cat('Model type:', attr(x, 'model.type'), '\n')
           cat('glmnet model alpha:', attr(x, 'alpha'), '\n')
           cat('glmnet model lambda:', attr(x, 'lambda'), '\n')
           if (is.null(attr(x, 'pen.factor'))) {
             cat('glmnet model penalty factor: not specified\n')
           } else {
             cat('glmnet model penalty factor: specified\n')
           }
           cat('Time-dependent AUC type:', attr(x, 'tauc.type'), '\n')
           cat('Evaluation time points for tAUC:', attr(x, 'tauc.time'))
         },

         glmnet.validate.cv = {
           cat('High-Dimensional Cox Model Validation Object\n')
           cat('Random seed:', attr(x, 'seed'), '\n')
           cat('Validation method: k-fold cross-validation\n')
           cat('Cross-validation folds:', attr(x, 'nfolds'), '\n')
           cat('Model type:', attr(x, 'model.type'), '\n')
           cat('glmnet model alpha:', attr(x, 'alpha'), '\n')
           cat('glmnet model lambda:', attr(x, 'lambda'), '\n')
           if (is.null(attr(x, 'pen.factor'))) {
             cat('glmnet model penalty factor: not specified\n')
           } else {
             cat('glmnet model penalty factor: specified\n')
           }
           cat('Time-dependent AUC type:', attr(x, 'tauc.type'), '\n')
           cat('Evaluation time points for tAUC:', attr(x, 'tauc.time'))
         },

         glmnet.validate.repeated.cv = {
           cat('High-Dimensional Cox Model Validation Object\n')
           cat('Random seed:', attr(x, 'seed'), '\n')
           cat('Validation method: repeated cross-validation\n')
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
           cat('Time-dependent AUC type:', attr(x, 'tauc.type'), '\n')
           cat('Evaluation time points for tAUC:', attr(x, 'tauc.time'))
         },

         ncvreg.validate.bootstrap = {
           cat('High-Dimensional Cox Model Validation Object\n')
           cat('Random seed:', attr(x, 'seed'), '\n')
           cat('Validation method: bootstrap\n')
           cat('Bootstrap samples:', attr(x, 'boot.times'), '\n')
           cat('Model type:', attr(x, 'model.type'), '\n')
           cat('ncvreg model gamma:', attr(x, 'gamma'), '\n')
           cat('ncvreg model alpha:', attr(x, 'alpha'), '\n')
           cat('ncvreg model lambda:', attr(x, 'lambda'), '\n')
           cat('Time-dependent AUC type:', attr(x, 'tauc.type'), '\n')
           cat('Evaluation time points for tAUC:', attr(x, 'tauc.time'))
         },

         ncvreg.validate.cv = {
           cat('High-Dimensional Cox Model Validation Object\n')
           cat('Random seed:', attr(x, 'seed'), '\n')
           cat('Validation method: k-fold cross-validation\n')
           cat('Cross-validation folds:', attr(x, 'nfolds'), '\n')
           cat('Model type:', attr(x, 'model.type'), '\n')
           cat('ncvreg model gamma:', attr(x, 'gamma'), '\n')
           cat('ncvreg model alpha:', attr(x, 'alpha'), '\n')
           cat('ncvreg model lambda:', attr(x, 'lambda'), '\n')
           cat('Time-dependent AUC type:', attr(x, 'tauc.type'), '\n')
           cat('Evaluation time points for tAUC:', attr(x, 'tauc.time'))
         },

         ncvreg.validate.repeated.cv = {
           cat('High-Dimensional Cox Model Validation Object\n')
           cat('Random seed:', attr(x, 'seed'), '\n')
           cat('Validation method: repeated cross-validation\n')
           cat('Cross-validation folds:', attr(x, 'nfolds'), '\n')
           cat('Cross-validation repeated times:', attr(x, 'rep.times'), '\n')
           cat('Model type:', attr(x, 'model.type'), '\n')
           cat('ncvreg model gamma:', attr(x, 'gamma'), '\n')
           cat('ncvreg model alpha:', attr(x, 'alpha'), '\n')
           cat('ncvreg model lambda:', attr(x, 'lambda'), '\n')
           cat('Time-dependent AUC type:', attr(x, 'tauc.type'), '\n')
           cat('Evaluation time points for tAUC:', attr(x, 'tauc.time'))
         },

         penalized.validate.bootstrap = {
           cat('High-Dimensional Cox Model Validation Object\n')
           cat('Random seed:', attr(x, 'seed'), '\n')
           cat('Validation method: bootstrap\n')
           cat('Bootstrap samples:', attr(x, 'boot.times'), '\n')
           cat('Model type:', attr(x, 'model.type'), '\n')
           cat('Fused lasso model lambda:', attr(x, 'lambda'), '\n')
           cat('Time-dependent AUC type:', attr(x, 'tauc.type'), '\n')
           cat('Evaluation time points for tAUC:', attr(x, 'tauc.time'))
         },

         penalized.validate.cv = {
           cat('High-Dimensional Cox Model Validation Object\n')
           cat('Random seed:', attr(x, 'seed'), '\n')
           cat('Validation method: k-fold cross-validation\n')
           cat('Cross-validation folds:', attr(x, 'nfolds'), '\n')
           cat('Model type:', attr(x, 'model.type'), '\n')
           cat('Fused lasso model lambda:', attr(x, 'lambda'), '\n')
           cat('Time-dependent AUC type:', attr(x, 'tauc.type'), '\n')
           cat('Evaluation time points for tAUC:', attr(x, 'tauc.time'))
         },

         penalized.validate.repeated.cv = {
           cat('High-Dimensional Cox Model Validation Object\n')
           cat('Random seed:', attr(x, 'seed'), '\n')
           cat('Validation method: repeated cross-validation\n')
           cat('Cross-validation folds:', attr(x, 'nfolds'), '\n')
           cat('Cross-validation repeated times:', attr(x, 'rep.times'), '\n')
           cat('Model type:', attr(x, 'model.type'), '\n')
           cat('Fused lasso model lambda:', attr(x, 'lambda'), '\n')
           cat('Time-dependent AUC type:', attr(x, 'tauc.type'), '\n')
           cat('Evaluation time points for tAUC:', attr(x, 'tauc.time'))
         }

  )

}

#' Summary of Validation Results
#'
#' Summary of Validation Results
#'
#' @param object An object \code{\link{hdnom.validate}}.
#' @param silent Print summary table header or not,
#' default is \code{FALSE}.
#' @param ... Other parameters (not used).
#'
#' @method summary hdnom.validate
#'
#' @export
#'
#' @examples
#' NULL
summary.hdnom.validate = function(object, silent = FALSE, ...) {

  if (!('hdnom.validate' %in% class(object)))
    stop('object class must be "hdnom.validate"')

  method = setdiff(class(object), 'hdnom.validate')

  if (grepl('validate.bootstrap', method)) {

    boot.times = attr(object, 'boot.times')
    tauc.time = attr(object, 'tauc.time')
    aucmat = matrix(NA, ncol = length(tauc.time), nrow = boot.times)
    for (i in 1L:boot.times) aucmat[i, ] = object[[i]]$auc
    summary_mat = rbind(apply(aucmat, 2, mean), apply(aucmat, 2, quantile))
    rownames(summary_mat) = c('Mean', 'Min', '0.25 Qt.',
                              'Median', '0.75 Qt.', 'Max')
    colnames(summary_mat) = tauc.time

  } else if (grepl('validate.cv', method)) {

    nfolds = attr(object, 'nfolds')
    tauc.time = attr(object, 'tauc.time')
    aucmat = matrix(NA, ncol = length(tauc.time), nrow = nfolds)
    for (i in 1L:nfolds) aucmat[i, ] = object[[i]]$auc
    summary_mat = rbind(apply(aucmat, 2, mean), apply(aucmat, 2, quantile))
    rownames(summary_mat) = c('Mean', 'Min', '0.25 Qt.',
                              'Median', '0.75 Qt.', 'Max')
    colnames(summary_mat) = tauc.time

  } else if (grepl('validate.repeated.cv', method)) {

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
    stop('hdnom.validate object is not valid')
  }

  if (!silent)
    cat('Time-Dependent AUC Summary at Evaluation Time Points\n')

  return(summary_mat)

}

#' Plot Optimism-Corrected Time-Dependent Discrimination Curves for Validation
#'
#' Plot Optimism-Corrected Time-Dependent Discrimination Curves for Validation
#'
#' @param x An object returned by \code{\link{hdnom.validate}}.
#' @param col.pal Color palette to use. Possible values are
#' \code{"JCO"}, \code{"Lancet"}, \code{"NPG"}, and \code{"AAAS"}.
#' Default is \code{"JCO"}.
#' @param ... Other parameters (not used).
#'
#' @method plot hdnom.validate
#'
#' @export
#'
#' @importFrom ggplot2 ggplot aes_string geom_point geom_line geom_point
#' geom_ribbon scale_x_continuous scale_fill_manual scale_colour_manual
#' theme_bw theme ylab
#'
#' @examples
#' NULL
plot.hdnom.validate =
  function(x, col.pal = c('JCO', 'Lancet', 'NPG', 'AAAS'),
           ...) {

    if (!('hdnom.validate' %in% class(x)))
      stop('object class must be "hdnom.validate"')

    df = as.data.frame(t(summary(x, silent = TRUE)))
    tauc_time = attr(x, 'tauc.time')

    # special processing for repeated cv
    if (any(grepl(pattern = 'validate.repeated.cv', class(x))))
      names(df) = sapply(strsplit(names(df), 'Mean of '), '[', 2L)

    df[, 'Time'] = tauc_time
    names(df)[which(names(df) == '0.25 Qt.')] = 'Qt25'
    names(df)[which(names(df) == '0.75 Qt.')] = 'Qt75'

    col.pal = match.arg(col.pal)
    col_pal = switch (
      col.pal,
      JCO   = palette.jco()[1], Lancet = palette.lancet()[1],
      NPG   = palette.npg()[1], AAAS   = palette.aaas()[1])

    ggplot(data = df, aes_string(x = 'Time', y = 'Mean')) +
      geom_point(colour = col_pal) +
      geom_line(colour = col_pal) +
      geom_point(data = df, aes_string(x = 'Time', y = 'Median'),
                 colour = col_pal) +
      geom_line(data = df, aes_string(x = 'Time', y = 'Median'),
                colour = col_pal, linetype = 'dashed') +
      geom_ribbon(data = df, aes_string(ymin = 'Qt25', ymax = 'Qt75'),
                  linetype = 0, alpha = 0.2) +
      geom_ribbon(data = df, aes_string(ymin = 'Min', ymax = 'Max'),
                  linetype = 0, alpha = 0.1) +
      scale_x_continuous(breaks = df$'Time') +
      theme_bw() +
      theme(legend.position = 'none') +
      ylab('Area under ROC')

  }
