# These functions are used to select models based on penalized partial
# likelihood, not tAUC, since it would be too complicated and unnecessary
# if we use tAUC when tuning alpha and lambda simultaneously.
#
# We should use these functions to build a model first, then validate
# the model using glmnet.validation and glmnet.calibrate functions.

#' Automatic alpha tuning function by k-fold cross-validation
#'
#' @return best model object and best alpha
#'
#' @importFrom glmnet cv.glmnet
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#'
#' @keywords internal
tuneAlpha = function(..., alphas, seed) {

  # Run models
  modelList <- foreach(alphas = alphas) %dopar% {
    set.seed(seed)
    cv.glmnet(..., alpha = alphas)
  }

  # Choose model for best alpha first (then lambda)
  # Criterion: penalized partial likelihood
  errors = unlist(lapply(modelList, function(x) min(sqrt(x$cvm))))

  return(list('best.model' = modelList[[which.min(errors)]],
              'best.alpha' = alphas[which.min(errors)]))

}

#' Adaptive elastic-net model selection for high-dimensional Cox models
#'
#' Automatic adaptive elastic-net model selection for high-dimensional
#' Cox models, evaluated by penalized partial-likelihood.
#'
#' @param x data matrix
#' @param y response matrix made by \code{\link[survival]{Surv}}
#' @param nfolds fold numbers of cross-validation
#' @param alphas alphas to tune in \code{\link[glmnet]{cv.glmnet}}
#' @param rule Model selection criterion, \code{"lambda.min"} or
#' \code{"lambda.1se"}. See \code{\link[glmnet]{cv.glmnet}}
#' for details.
#' @param seeds two random seeds for cross-validation fold division
#' in the two steps
#'
#' @importFrom glmnet glmnet
#'
#' @export hdcox.aenet
#'
#' @examples
#' library("glmnet")
#' library("survival")
#' library("rms")
#' library("doParallel")
#' registerDoParallel(detectCores())
#'
#' # Load imputed SMART data
#' data("smart")
#' x = as.matrix(smart[, -c(1, 2)])
#' time = smart$TEVENT
#' event = smart$EVENT
#' y = Surv(time, event)
#'
#' # Fit Cox model by adaptive elastic-net penalization
#' aenetfit = hdcox.aenet(x, y, nfolds = 10, rule = 'lambda.1se', seeds = c(5, 7))
#'
#' # Prepare data for glmnet.nomogram
#' x.df = as.data.frame(x)
#' dd = datadist(x.df)
#' options(datadist = "dd")
#'
#' # Generate glmnet.nomogram objects and draw nomogram
#' nom = glmnet.nomogram(aenetfit$aenet_all, x, time, event, x.df,
#'                       s = aenetfit$best_lambda_aenet, pred.at = 365 * 2,
#'                       funlabel = "2-Year Overall Survival Probability")
#'
#' plot(nom)
hdcox.aenet = function(x, y, nfolds = 5L, alphas = seq(0.05, 0.95, 0.05),
                       rule = c('lambda.min', 'lambda.1se'),
                       seeds = c(1001, 1002)) {

  rule = match.arg(rule)

  # Tuning alpha for the both two stages of adaptive enet estimation
  enet_y = tuneAlpha(x, y, family = 'cox',
                     nfolds = nfolds, alphas = alphas,
                     seed = seeds[1L])

  # fit the model on all the data use the parameters got by CV
  best_alpha_enet  = enet_y$best.alpha

  if (rule == 'lambda.min') {
    best_lambda_enet = enet_y$best.model$lambda.min
  } else if (rule == 'lambda.1se') {
    best_lambda_enet = enet_y$best.model$lambda.1se
  }

  enet_all = glmnet(x, y, family = 'cox',
                    lambda = best_lambda_enet,
                    alpha  = best_alpha_enet)

  bhat = as.matrix(enet_all$beta)
  if(all(bhat == 0)) bhat = rep(.Machine$double.eps * 2, length(bhat))

  # adaptive penalty
  adpen = (1/pmax(abs(bhat), .Machine$double.eps))

  aenet_y = tuneAlpha(x, y, family = 'cox', nfolds = nfolds, alphas = alphas,
                      exclude = which(bhat == 0), penalty.factor = adpen,
                      seed = seeds[2L])

  # fit the model on all the data use the parameters got by CV
  best_alpha_aenet  = aenet_y$best.alpha
  best_lambda_aenet = aenet_y$best.model$lambda.min
  aenet_all = glmnet(x, y, family = 'cox',
                     exclude = which(bhat == 0),
                     lambda  = best_lambda_aenet,
                     alpha   = best_alpha_aenet,
                     penalty.factor = adpen)

  coxaenet_model = list('best_alpha_enet' = best_alpha_enet,
                        'best_lambda_enet' = best_lambda_enet,
                        'enet_all' = enet_all,
                        'best_alpha_aenet' = best_alpha_aenet,
                        'best_lambda_aenet' = best_lambda_aenet,
                        'aenet_all' = aenet_all)

  class(coxaenet_model) = 'hdcox.model.aenet'

  return(coxaenet_model)

}

#' Adaptive lasso model selection for high-dimensional Cox models
#'
#' Automatic adaptive lasso model selection for high-dimensional
#' Cox models, evaluated by penalized partial-likelihood.
#'
#' @param x data matrix
#' @param y response matrix made by \code{\link[survival]{Surv}}
#' @param nfolds fold numbers of cross-validation
#' @param rule Model selection criterion, \code{"lambda.min"} or
#' \code{"lambda.1se"}. See \code{\link[glmnet]{cv.glmnet}}
#' for details.
#' @param seeds two random seeds for cross-validation fold division
#' in the two steps
#'
#' @export hdcox.alasso
#'
#' @examples
#' library("glmnet")
#' library("survival")
#' library("rms")
#' library("doParallel")
#' registerDoParallel(detectCores())
#'
#' # Load imputed SMART data
#' data("smart")
#' x = as.matrix(smart[, -c(1, 2)])
#' time = smart$TEVENT
#' event = smart$EVENT
#' y = Surv(time, event)
#'
#' # Fit Cox model by adaptive lasso penalization
#' alassofit = hdcox.alasso(x, y, nfolds = 10, rule = 'lambda.1se', seeds = c(7, 11))
#'
#' # Prepare data for glmnet.nomogram
#' x.df = as.data.frame(x)
#' dd = datadist(x.df)
#' options(datadist = "dd")
#'
#' # Generate glmnet.nomogram objects and draw nomogram
#' nom = glmnet.nomogram(alassofit$alasso_all, x, time, event, x.df,
#'                       s = alassofit$best_lambda_alasso, pred.at = 365 * 2,
#'                       funlabel = "2-Year Overall Survival Probability")
#'
#' plot(nom)
hdcox.alasso = function(x, y, nfolds = 5L,
                        rule = c('lambda.min', 'lambda.1se'),
                        seeds = c(1001, 1002)) {

  # Tuning lambda for the both two stages of adaptive lasso estimation
  set.seed(seeds[1L])
  lasso_y = cv.glmnet(x, y, family = 'cox', nfolds = nfolds, alpha = 0)

  # fit the model on all the data use the parameters got by CV
  if (rule == 'lambda.min') {
    best_lambda_lasso = lasso_y$lambda.min
  } else if (rule == 'lambda.1se') {
    best_lambda_lasso = lasso_y$lambda.1se
  }

  lasso_all = glmnet(x, y, family = 'cox',
                     lambda = best_lambda_lasso, alpha = 0)

  bhat = as.matrix(lasso_all$beta)
  if(all(bhat == 0)) bhat = rep(.Machine$double.eps * 2, length(bhat))

  # adaptive penalty
  adpen = (1/pmax(abs(bhat), .Machine$double.eps))

  set.seed(seeds[2L])
  alasso_y = cv.glmnet(x, y, family = 'cox', nfolds = nfolds, alpha = 1,
                       penalty.factor = adpen)

  # fit the model on all the data use the parameters got by CV
  if (rule == 'lambda.min') {
    best_lambda_alasso = alasso_y$lambda.min
  } else if (rule == 'lambda.1se') {
    best_lambda_alasso = alasso_y$lambda.1se
  }

  alasso_all = glmnet(x, y, family = 'cox', lambda = best_lambda_alasso,
                      alpha = 1, penalty.factor = adpen)

  coxalasso_model = list('best_lambda_lasso' = best_lambda_lasso,
                         'lasso_all' = lasso_all,
                         'best_lambda_alasso' = best_lambda_alasso,
                         'alasso_all' = alasso_all)

  class(coxalasso_model) = 'hdcox.model.alasso'

  return(coxalasso_model)

}

#' Elastic-net model selection for high-dimensional Cox models
#'
#' Automatic elastic-net model selection for high-dimensional
#' Cox models, evaluated by penalized partial-likelihood.
#'
#' @param x data matrix
#' @param y response matrix made by \code{\link[survival]{Surv}}
#' @param nfolds fold numbers of cross-validation
#' @param alphas alphas to tune in \code{\link[glmnet]{cv.glmnet}}
#' @param rule Model selection criterion, \code{"lambda.min"} or
#' \code{"lambda.1se"}. See \code{\link[glmnet]{cv.glmnet}}
#' for details.
#' @param seed a random seed for cross-validation fold division
#'
#' @export hdcox.enet
#'
#' @examples
#' library("glmnet")
#' library("survival")
#' library("rms")
#' library("doParallel")
#' registerDoParallel(detectCores())
#'
#' # Load imputed SMART data
#' data("smart")
#' x = as.matrix(smart[, -c(1, 2)])
#' time = smart$TEVENT
#' event = smart$EVENT
#' y = Surv(time, event)
#'
#' # Fit Cox model by adaptive elastic-net penalization
#' enetfit = hdcox.enet(x, y, nfolds = 10, rule = 'lambda.1se', seed = 11)
#'
#' # Prepare data for glmnet.nomogram
#' x.df = as.data.frame(x)
#' dd = datadist(x.df)
#' options(datadist = "dd")
#'
#' # Generate glmnet.nomogram objects and draw nomogram
#' nom = glmnet.nomogram(enetfit$enet_all, x, time, event, x.df,
#'                       s = enetfit$best_lambda_enet, pred.at = 365 * 2,
#'                       funlabel = "2-Year Overall Survival Probability")
#'
#' plot(nom)
hdcox.enet = function(x, y, nfolds = 5L, alphas = seq(0.05, 0.95, 0.05),
                      rule = c('lambda.min', 'lambda.1se'), seed = 1001) {

  enet_y = tuneAlpha(x, y, family = 'cox',
                     nfolds = nfolds, alphas = alphas,
                     seed = seed)

  # fit the model on all the data use the parameters got by CV
  best_alpha_enet  = enet_y$best.alpha

  if (rule == 'lambda.min') {
    best_lambda_enet = enet_y$best.model$lambda.min
  } else if (rule == 'lambda.1se') {
    best_lambda_enet = enet_y$best.model$lambda.1se
  }

  enet_all = glmnet(x, y, family = 'cox',
                    lambda = best_lambda_enet,
                    alpha  = best_alpha_enet)

  coxenet_model = list('best_alpha_enet' = best_alpha_enet,
                       'best_lambda_enet' = best_lambda_enet,
                       'enet_all' = enet_all)

  class(coxenet_model) = 'hdcox.model.enet'

  return(coxenet_model)

}

#' Lasso model selection for high-dimensional Cox models
#'
#' Automatic lasso model selection for high-dimensional
#' Cox models, evaluated by penalized partial-likelihood.
#'
#' @param x data matrix
#' @param y response matrix made by \code{\link[survival]{Surv}}
#' @param nfolds fold numbers of cross-validation
#' @param rule Model selection criterion, \code{"lambda.min"} or
#' \code{"lambda.1se"}. See \code{\link[glmnet]{cv.glmnet}}
#' for details.
#' @param seed a random seed for cross-validation fold division
#'
#' @export hdcox.lasso
#'
#' @examples
#' library("glmnet")
#' library("survival")
#' library("rms")
#' library("doParallel")
#' registerDoParallel(detectCores())
#'
#' # Load imputed SMART data
#' data("smart")
#' x = as.matrix(smart[, -c(1, 2)])
#' time = smart$TEVENT
#' event = smart$EVENT
#' y = Surv(time, event)
#'
#' # Fit Cox model by adaptive lasso penalization
#' lassofit = hdcox.lasso(x, y, nfolds = 10, rule = 'lambda.1se', seed = 11)
#'
#' # Prepare data for glmnet.nomogram
#' x.df = as.data.frame(x)
#' dd = datadist(x.df)
#' options(datadist = "dd")
#'
#' # Generate glmnet.nomogram objects and draw nomogram
#' nom = glmnet.nomogram(lassofit$lasso_all, x, time, event, x.df,
#'                       s = lassofit$best_lambda_lasso, pred.at = 365 * 2,
#'                       funlabel = "2-Year Overall Survival Probability")
#'
#' plot(nom)
hdcox.lasso = function(x, y, nfolds = 5L,
                       rule = c('lambda.min', 'lambda.1se'),
                       seed = 1001) {

  set.seed(seed)
  lasso_y = cv.glmnet(x, y, family = 'cox', nfolds = nfolds, alpha = 1)

  # fit the model on all the data use the parameters got by CV
  if (rule == 'lambda.min') {
    best_lambda_lasso = lasso_y$lambda.min
  } else if (rule == 'lambda.1se') {
    best_lambda_lasso = lasso_y$lambda.1se
  }

  lasso_all = glmnet(x, y, family = 'cox',
                     lambda = best_lambda_lasso, alpha = 1)

  coxlasso_model = list('best_lambda_lasso' = best_lambda_lasso,
                        'lasso_all' = lasso_all)

  class(coxlasso_model) = 'hdcox.model.lasso'

  return(coxlasso_model)

}
