#' Externally Validate High-Dimensional Cox Models with Time-Dependent AUC
#'
#' Externally Validate High-Dimensional Cox Models with Time-Dependent AUC
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
#' @param tauc.type Type of time-dependent AUC.
#' Including \code{"CD"} proposed by Chambless and Diao (2006).,
#' \code{"SZ"} proposed by Song and Zhou (2008).,
#' \code{"UNO"} proposed by Uno et al. (2007).
#' @param tauc.time Numeric vector. Time points at which to evaluate
#' the time-dependent AUC.
#'
#' @export hdnom.external.validate
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
#' # Take the next 1000 samples as external validation data
#' # In practice, usually use data collected in other studies
#' x_new = as.matrix(smart[, -c(1, 2)])[1001:2000, ]
#' time_new = smart$TEVENT[1001:2000]
#' event_new = smart$EVENT[1001:2000]
#'
#' # Fit Cox model by lasso penalization
#' lassofit = hdcox.lasso(x, Surv(time, event), nfolds = 5, rule = "lambda.1se", seed = 11)
#'
#' # External validation with time-dependent AUC
#' val.ext =
#'   hdnom.external.validate(lassofit, x, time, event,
#'                           x_new, time_new, event_new,
#'                           tauc.type = "UNO",
#'                           tauc.time = seq(0.25, 2, 0.25) * 365)
#'
#' print(val.ext)
#' summary(val.ext)
#' plot(val.ext)
#'
# ### Testing fused lasso, MCP and Snet models ###
# library("survival")
#
# # Load imputed SMART data
# data(smart)
# # Use first 600 samples as training data
# # (the data used for internal validation)
# x = as.matrix(smart[, -c(1, 2)])[1:600, ]
# time = smart$TEVENT[1:600]
# event = smart$EVENT[1:600]
#
# # Take 500 samples as external validation data.
# # In practice, usually use data collected in other studies.
# x_new = as.matrix(smart[, -c(1, 2)])[1001:1500, ]
# time_new = smart$TEVENT[1001:1500]
# event_new = smart$EVENT[1001:1500]
#
# flassofit = hdcox.flasso(x, Surv(time, event), nfolds = 5, seed = 11)
# scadfit = hdcox.mcp(x, Surv(time, event), nfolds = 5, seed = 11)
# mnetfit = hdcox.snet(x, Surv(time, event), nfolds = 5, seed = 11)
#
# val.ext1 = hdnom.external.validate(flassofit, x, time, event,
#                                    x_new, time_new, event_new,
#                                    tauc.type = "UNO",
#                                    tauc.time = seq(0.25, 2, 0.25) * 365)
#
# val.ext2 = hdnom.external.validate(scadfit, x, time, event,
#                                    x_new, time_new, event_new,
#                                    tauc.type = "CD",
#                                    tauc.time = seq(0.25, 2, 0.25) * 365)
#
# val.ext3 = hdnom.external.validate(mnetfit, x, time, event,
#                                    x_new, time_new, event_new,
#                                    tauc.type = "SZ",
#                                    tauc.time = seq(0.25, 2, 0.25) * 365)
#
# print(val.ext1)
# summary(val.ext1)
# plot(val.ext1)
#
# print(val.ext2)
# summary(val.ext2)
# plot(val.ext2)
#
# print(val.ext3)
# summary(val.ext3)
# plot(val.ext3)
hdnom.external.validate = function(object, x, time, event,
                                   x_new, time_new, event_new,
                                   tauc.type = c("CD", "SZ", "UNO"), tauc.time) {

  if (!('hdcox.model' %in% class(object)))
    stop('object must be of class "hdcox.model" fitted by hdcox.* functions')

  model.type = gsub('hdcox.model.', '', setdiff(class(object), 'hdcox.model'))
  model_object = object[[paste0(model.type, '_model')]]

  if (!('matrix' %in% class(x_new))) stop('x_new must be a matrix')

  tauc.type = match.arg(tauc.type)

  x_tr = x
  time_tr = time
  event_tr = event
  y_tr = Surv(time_tr, event_tr)
  x_te  = x_new
  time_te = time_new
  event_te = event_new
  y_te = Surv(time_te, event_te)

  if (model.type %in% c('lasso', 'alasso', 'enet', 'aenet')) {
    tauc =
      glmnet.external.validate.internal(
        object = model_object, x_tr = x_tr, x_te = x_te, y_tr = y_tr, y_te = y_te,
        tauc.type = tauc.type, tauc.time = tauc.time
      )
  }

  if (model.type %in% c('mcp', 'mnet', 'scad', 'snet')) {
    tauc =
      ncvreg.external.validate.internal(
        object = model_object, x_tr = x_tr, x_te = x_te, y_tr = y_tr, y_te = y_te,
        tauc.type = tauc.type, tauc.time = tauc.time
      )
  }

  if (model.type %in% c('flasso')) {
    tauc =
      penalized.external.validate.internal(
        object = model_object, x_tr = x_tr, x_te = x_te, y_tr = y_tr, y_te = y_te,
        tauc.type = tauc.type, tauc.time = tauc.time
      )
  }

  if (model.type %in% c('lasso', 'alasso', 'enet', 'aenet')) {
    class(tauc) = c('hdnom.external.validate',
                    'glmnet.external.validate')
    attr(tauc, 'model.type') = model.type
    attr(tauc, 'tauc.type')  = tauc.type
    attr(tauc, 'tauc.time')  = tauc.time
  }

  if (model.type %in% c('mcp', 'mnet', 'scad', 'snet')) {
    class(tauc) = c('hdnom.external.validate',
                    'ncvreg.external.validate')
    attr(tauc, 'model.type') = model.type
    attr(tauc, 'tauc.type')  = tauc.type
    attr(tauc, 'tauc.time')  = tauc.time
  }

  if (model.type %in% c('flasso')) {
    class(tauc) = c('hdnom.external.validate',
                    'penalized.external.validate')
    attr(tauc, 'model.type') = model.type
    attr(tauc, 'tauc.type')  = tauc.type
    attr(tauc, 'tauc.time')  = tauc.time
  }

  tauc

}

#' Compute External Validation Measures for glmnet Objects
#'
#' @importFrom survAUC AUC.cd AUC.sh AUC.uno
#' @importFrom glmnet glmnet
#' @importFrom survival Surv
#'
#' @return time-dependent AUC (tAUC) value
#'
#' @keywords internal
glmnet.external.validate.internal = function(object, x_tr, x_te, y_tr, y_te,
                                             tauc.type, tauc.time) {

  lp_tr = as.vector(predict(object, newx = x_tr, type = 'link'))
  lp_te = as.vector(predict(object, newx = x_te, type = 'link'))

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

#' Compute External Validation Measures for ncvreg Model Objects
#'
#' @importFrom survAUC AUC.cd AUC.sh AUC.uno
#' @importFrom ncvreg ncvsurv
#' @importFrom survival Surv
#'
#' @return time-dependent AUC (tAUC) value
#'
#' @keywords internal
ncvreg.external.validate.internal = function(object, x_tr, x_te, y_tr, y_te,
                                             tauc.type, tauc.time) {

  lp_tr = as.vector(predict(object, X = x_tr, type = 'link'))
  lp_te = as.vector(predict(object, X = x_te, type = 'link'))

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

#' Compute External Validation Measures for "penalized" Model Objects
#'
#' @importFrom survAUC AUC.cd AUC.sh AUC.uno
#' @importFrom penalized penalized
#' @importFrom survival Surv
#'
#' @return time-dependent AUC (tAUC) value
#'
#' @keywords internal
penalized.external.validate.internal = function(object, x_tr, x_te, y_tr, y_te,
                                                tauc.type, tauc.time) {

  lp_tr = as.vector(object@'lin.pred')
  lp_te = as.vector(x_te %*% as.matrix(object@'penalized'))

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

#' Print External Validation Results
#'
#' Print External Validation Results
#'
#' @param x An object returned by \code{\link{hdnom.external.validate}}.
#' @param ... Other parameters (not used).
#'
#' @method print hdnom.external.validate
#'
#' @export
#'
#' @examples
#' NULL
print.hdnom.external.validate = function(x, ...) {

  if (!('hdnom.external.validate' %in% class(x)))
    stop('object class must be "hdnom.external.validate"')

  method = setdiff(class(x), 'hdnom.external.validate')

  switch(method,

         glmnet.external.validate = {
           cat('High-Dimensional Cox Model External Validation Object\n')
           cat('Model type:', attr(x, 'model.type'), '\n')
           cat('Time-dependent AUC type:', attr(x, 'tauc.type'), '\n')
           cat('Evaluation time points for tAUC:', attr(x, 'tauc.time'))
         },

         ncvreg.external.validate = {
           cat('High-Dimensional Cox Model External Validation Object\n')
           cat('Model type:', attr(x, 'model.type'), '\n')
           cat('Time-dependent AUC type:', attr(x, 'tauc.type'), '\n')
           cat('Evaluation time points for tAUC:', attr(x, 'tauc.time'))
         },

         penalized.external.validate = {
           cat('High-Dimensional Cox Model External Validation Object\n')
           cat('Model type:', attr(x, 'model.type'), '\n')
           cat('Time-dependent AUC type:', attr(x, 'tauc.type'), '\n')
           cat('Evaluation time points for tAUC:', attr(x, 'tauc.time'))
         }

  )

}

#' Summary of External Validation Results
#'
#' Summary of External Validation Results
#'
#' @param object An object \code{\link{hdnom.external.validate}}.
#' @param silent Print summary table header or not,
#' default is \code{FALSE}.
#' @param ... Other parameters (not used).
#'
#' @method summary hdnom.external.validate
#'
#' @export
#'
#' @examples
#' NULL
summary.hdnom.external.validate = function(object, silent = FALSE, ...) {

  if (!('hdnom.external.validate' %in% class(object)))
    stop('object class must be "hdnom.external.validate"')

  tauc.time = attr(object, 'tauc.time')
  aucmat = matrix(NA, ncol = length(tauc.time), nrow = 1L)
  aucmat[1L, ] = object$'auc'
  rownames(aucmat) = 'AUC'
  colnames(aucmat) = tauc.time

  if (!silent)
    cat('Time-Dependent AUC Summary at Evaluation Time Points\n')

  return(aucmat)

}

#' Plot Time-Dependent Discrimination Curves for External Validation
#'
#' Plot Time-Dependent Discrimination Curves for External Validation
#'
#' @param x An object returned by \code{\link{hdnom.external.validate}}.
#' @param col.pal Color palette to use. Possible values are
#' \code{"JCO"}, \code{"Lancet"}, \code{"NPG"}, and \code{"AAAS"}.
#' Default is \code{"JCO"}.
#' @param ... Other parameters (not used).
#'
#' @method plot hdnom.external.validate
#'
#' @export
#'
#' @importFrom ggplot2 ggplot aes_string geom_point geom_line geom_point
#' geom_ribbon scale_x_continuous scale_fill_manual scale_colour_manual
#' theme_bw theme ylab
#'
#' @examples
#' NULL
plot.hdnom.external.validate =
  function(x, col.pal = c('JCO', 'Lancet', 'NPG', 'AAAS'), ...) {

  if (!('hdnom.external.validate' %in% class(x)))
    stop('object class must be "hdnom.external.validate"')

  df = as.data.frame(t(summary(x, silent = TRUE)))
  tauc_time = attr(x, 'tauc.time')

  df[, 'Time'] = tauc_time

  col.pal = match.arg(col.pal)
  col_pal = switch (
    col.pal,
    JCO   = palette.jco()[1], Lancet = palette.lancet()[1],
    NPG   = palette.npg()[1], AAAS   = palette.aaas()[1])

  ggplot(data = df, aes_string(x = 'Time', y = 'AUC')) +
    geom_point(colour = col_pal) +
    geom_line(colour = col_pal) +
    scale_x_continuous(breaks = df$'Time') +
    theme_bw() +
    theme(legend.position = 'none') +
    ylab('Area under ROC')

}
