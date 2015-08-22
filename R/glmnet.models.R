#' Automatic alpha tuning function by k-fold cross-validation
#'
#' @return best model object and best alpha
#'
#' @importFrom glmnet cv.glmnet
#'
#' @keywords internal
tuneAlpha = function(..., alphas) {

  `%xdopar%` = foreach::`%dopar%`

  # Run models
  modelList <- foreach::foreach(alpha = alphas) %xdopar% {
    cv.glmnet(..., alpha = alpha)
  }

  # Choose best model by best alpha # TODO: choose best model by tAUC instead of cvm
  errors = unlist(lapply(modelList, function(x) min(sqrt(x$cvm))))

  return(list('best.model' = modelList[[which.min(errors)]],
              'best.alpha' = alphas[which.min(errors)]))

}

#' Adaptive elastic-net model selection for high-dimensional Cox models
#'
#' Automatic adaptive elastic-net model selection for high-dimensional
#' Cox models, evaluated by time-dependent AUC (tAUC).
#'
#' @param x data matrix
#' @param y response matrix made by \code{\link[survival]{Surv}}
#'
#' @importFrom glmnet glmnet
#'
#' @export hdcox.aenet
#'
#' @examples
#' set.seed(1001)
hdcox.aenet = function(x, y, nfolds = 5L, alphas = seq(0.05, 0.95, 0.05)) {

  # Tuning alpha for the both two stages of adaptive enet estimation

  enet_y = tuneAlpha(x = x, y = y, family = 'cox',
                     nfolds = nfolds, alphas = alphas)

  # fit the model on all the data use the parameters got by CV
  best_alpha_enet  = enet_y$best.alpha
  best_lambda_enet = enet_y$best.model$lambda.min
  enet_all = glmnet(x, y, family = 'cox',
                    lambda = best_lambda_enet,
                    alpha  = best_alpha_enet)

  bhat = as.matrix(enet_all$beta)
  if(all(bhat == 0)) bhat = rep(.Machine$double.eps * 2, length(bhat))

  # adaptive penalty
  adpen = (1/pmax(abs(bhat), .Machine$double.eps))

  aenet_y = tuneAlpha(x, y, family = 'cox', nfolds = nfolds, alphas = alphas,
                      exclude = which(bhat == 0), penalty.factor = adpen)

  # fit the model on all the data use the parameters got by CV
  best_alpha_aenet  = aenet_y$best.alpha
  best_lambda_aenet = aenet_y$best.model$lambda.min
  aenet_all = glmnet(x, y, family = 'cox',
                     exclude = which(bhat == 0),
                     lambda = best_lambda_aenet,
                     alpha  = best_alpha_aenet,
                     penalty.factor = adpen)

  coxaenet_model = list('best_alpha_enet' = best_alpha_enet,
                        'best_lambda_enet' = best_lambda_enet,
                        'enet_all' = enet_all,
                        'best_alpha_aenet' = best_alpha_aenet,
                        'best_lambda_aenet' = best_lambda_aenet,
                        'aenet_all' = aenet_all)

  return(coxaenet_model)

}

# # Example
# # how to calculate time-dependent AUC
# library("glmnet")
# library("survival")
# library("doParallel")
# registerDoParallel(detectCores())
#
# data("smart")
# x.tr = as.matrix(smart[1:2900, 3:29])
# x.te = as.matrix(smart[2901:3873, 3:29])
# y.tr = Surv(smart[1:2900, 1], smart[1:2900, 2])
# y.te = Surv(smart[2901:3873, 1], smart[2901:3873, 2])
#
# set.seed(1001)
# allfit = hdcox.aenet(x.tr, y.tr)
# fit = allfit$aenet_all
#
# library('survAUC')
# lp.tr = as.vector(predict(fit, newx = x.tr, type = 'link'))
# lp.te = as.vector(predict(fit, newx = x.te, type = 'link'))
# times = seq(0.25, 2, 0.25) * 365
#
# cd = AUC.cd(y.tr, y.te, lp.tr, lp.te, times)$auc
# sh = AUC.sh(y.tr, y.te, lp.tr, lp.te, times)$auc
# uno = AUC.uno(y.tr, y.te, lp.te, times)$auc
#
# tauc = data.frame(
#   tAUC.type = rep(c('AUC.cd', 'AUC.sh', 'AUC.uno'), each = length(times)),
#   tAUC.value = c(cd, sh, uno),
#   time.month = rep(times/365*12, 3))
#
# # line plot
# library('ggplot2')
# cbPalette = c("#E69F00", "#009E73", "#0072B2")
# ggplot(data = tauc, aes(x = time.month, y = tAUC.value,
#                         group = tAUC.type, colour = tAUC.type)) +
#   geom_line(size = 1) +
#   geom_point(size = 3) +
#   scale_colour_manual(values = cbPalette) +
#   theme_bw()

#' Adaptive lasso model selection for high-dimensional Cox models
#'
#' Automatic adaptive lasso model selection for high-dimensional
#' Cox models, evaluated by time-dependent AUC (tAUC).
#'
#' @param x data matrix
#' @param y response matrix made by \code{\link[survival]{Surv}}
#'
#' @export hdcox.alasso
#'
#' @examples
#' set.seed(1001)
hdcox.alasso = function(x, y, nfolds = 5L, alphas = seq(0.05, 0.95, 0.05)) {
  NULL
}

#' Elastic-net model selection for high-dimensional Cox models
#'
#' Automatic elastic-net model selection for high-dimensional
#' Cox models, evaluated by time-dependent AUC (tAUC).
#'
#' @param x data matrix
#' @param y response matrix made by \code{\link[survival]{Surv}}
#'
#' @export hdcox.enet
#'
#' @examples
#' set.seed(1001)
hdcox.enet = function(x, y, nfolds = 5L, alphas = seq(0.05, 0.95, 0.05)) {
  NULL
}

#' Lasso model selection for high-dimensional Cox models
#'
#' Automatic lasso model selection for high-dimensional
#' Cox models, evaluated by time-dependent AUC (tAUC).
#'
#' @param x data matrix
#' @param y response matrix made by \code{\link[survival]{Surv}}
#'
#' @export hdcox.lasso
#'
#' @examples
#' set.seed(1001)
hdcox.lasso = function(x, y, nfolds = 5L, alphas = seq(0.05, 0.95, 0.05)) {
  NULL
}
