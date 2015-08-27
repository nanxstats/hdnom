# These functions are used to select models based on penalized partial
# likelihood, not tAUC, since it would be too complicated and unnecessary
# if we use tAUC when tuning alpha and lambda simultaneously.
#
# We should use these functions to build a model first, then validate
# the model using hdnom.validation and hdnom.calibrate functions.

#' Automatic alpha tuning function by k-fold cross-validation
#'
#' @return best model object and best alpha
#'
#' @importFrom glmnet cv.glmnet
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#'
#' @keywords internal
glmnet.tune.alpha = function(..., alphas, seed) {

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
#' @param x Data matrix.
#' @param y Response matrix made with \code{\link[survival]{Surv}}.
#' @param nfolds Fold numbers of cross-validation.
#' @param alphas Alphas to tune in \code{\link[glmnet]{cv.glmnet}}.
#' @param rule Model selection criterion, \code{"lambda.min"} or
#' \code{"lambda.1se"}. See \code{\link[glmnet]{cv.glmnet}}
#' for details.
#' @param seeds Two random seeds for cross-validation fold division
#' in two estimation steps.
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
#' \donttest{aenetfit = hdcox.aenet(x, y, nfolds = 10, rule = "lambda.1se", seeds = c(5, 7))
#'
#' # Prepare data for hdnom.nomogram
#' x.df = as.data.frame(x)
#' dd = datadist(x.df)
#' options(datadist = "dd")
#'
#' # Generate hdnom.nomogram objects and plot nomogram
#' nom = hdnom.nomogram(aenetfit$aenet_model, x, time, event, x.df,
#'                      lambda = aenetfit$aenet_best_lambda, pred.at = 365 * 2,
#'                      funlabel = "2-Year Overall Survival Probability")
#'
#' plot(nom)}
# Skipped testing the above code since parallel computation is required
hdcox.aenet = function(x, y, nfolds = 5L, alphas = seq(0.05, 0.95, 0.05),
                       rule = c('lambda.min', 'lambda.1se'),
                       seeds = c(1001, 1002)) {

  rule = match.arg(rule)

  # Tuning alpha for the both two stages of adaptive enet estimation
  enet_y = glmnet.tune.alpha(x, y, family = 'cox',
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

  aenet_y = glmnet.tune.alpha(x, y, family = 'cox', nfolds = nfolds,
                              exclude = which(bhat == 0),
                              penalty.factor = adpen,
                              alphas = alphas,
                              seed = seeds[2L])

  # fit the model on all the data use the parameters got by CV
  best_alpha_aenet  = aenet_y$best.alpha

  if (rule == 'lambda.min') {
    best_lambda_aenet = aenet_y$best.model$lambda.min
  } else if (rule == 'lambda.1se') {
    best_lambda_aenet = aenet_y$best.model$lambda.1se
  }

  aenet_all = glmnet(x, y, family = 'cox',
                     exclude = which(bhat == 0),
                     lambda  = best_lambda_aenet,
                     penalty.factor = adpen,
                     alpha   = best_alpha_aenet)

  adpen_vec = as.vector(adpen)
  adpen_name = rownames(adpen)
  names(adpen_vec) = adpen_name

  coxaenet_model = list('enet_best_alpha' = best_alpha_enet,
                        'enet_best_lambda' = best_lambda_enet,
                        'enet_model' = enet_all,
                        'aenet_best_alpha' = best_alpha_aenet,
                        'aenet_best_lambda' = best_lambda_aenet,
                        'aenet_model' = aenet_all,
                        'pen_factor' = adpen_vec)

  class(coxaenet_model) = 'hdcox.model.aenet'

  return(coxaenet_model)

}

#' Adaptive lasso model selection for high-dimensional Cox models
#'
#' Automatic adaptive lasso model selection for high-dimensional
#' Cox models, evaluated by penalized partial-likelihood.
#'
#' @param x Data matrix.
#' @param y Response matrix made by \code{\link[survival]{Surv}}.
#' @param nfolds Fold numbers of cross-validation.
#' @param rule Model selection criterion, \code{"lambda.min"} or
#' \code{"lambda.1se"}. See \code{\link[glmnet]{cv.glmnet}}
#' for details.
#' @param seeds Two random seeds for cross-validation fold division
#' in two estimation steps.
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
#' \donttest{alassofit = hdcox.alasso(x, y, nfolds = 10, rule = "lambda.1se", seeds = c(7, 11))
#'
#' # Prepare data for hdnom.nomogram
#' x.df = as.data.frame(x)
#' dd = datadist(x.df)
#' options(datadist = "dd")
#'
#' # Generate hdnom.nomogram objects and plot nomogram
#' nom = hdnom.nomogram(alassofit$alasso_model, x, time, event, x.df,
#'                      lambda = alassofit$alasso_best_lambda, pred.at = 365 * 2,
#'                      funlabel = "2-Year Overall Survival Probability")
#'
#' plot(nom)}
# Skipped testing the above code since parallel computation is required
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

  adpen_vec = as.vector(adpen)
  adpen_name = rownames(adpen)
  names(adpen_vec) = adpen_name

  coxalasso_model = list('ridge_best_lambda' = best_lambda_lasso,
                         'ridge_model' = lasso_all,
                         'alasso_best_lambda' = best_lambda_alasso,
                         'alasso_model' = alasso_all,
                         'pen_factor' = adpen_vec)

  class(coxalasso_model) = 'hdcox.model.alasso'

  return(coxalasso_model)

}

#' Elastic-net model selection for high-dimensional Cox models
#'
#' Automatic elastic-net model selection for high-dimensional
#' Cox models, evaluated by penalized partial-likelihood.
#'
#' @param x Data matrix.
#' @param y Response matrix made by \code{\link[survival]{Surv}}.
#' @param nfolds Fold numbers of cross-validation.
#' @param alphas Alphas to tune in \code{\link[glmnet]{cv.glmnet}}.
#' @param rule Model selection criterion, \code{"lambda.min"} or
#' \code{"lambda.1se"}. See \code{\link[glmnet]{cv.glmnet}}
#' for details.
#' @param seed A random seed for cross-validation fold division.
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
#' \donttest{enetfit = hdcox.enet(x, y, nfolds = 10, rule = "lambda.1se", seed = 11)
#'
#' # Prepare data for hdnom.nomogram
#' x.df = as.data.frame(x)
#' dd = datadist(x.df)
#' options(datadist = "dd")
#'
#' # Generate hdnom.nomogram objects and plot nomogram
#' nom = hdnom.nomogram(enetfit$enet_model, x, time, event, x.df,
#'                      lambda = enetfit$enet_best_lambda, pred.at = 365 * 2,
#'                      funlabel = "2-Year Overall Survival Probability")
#'
#' plot(nom)}
# Skipped testing the above code since parallel computation is required
hdcox.enet = function(x, y, nfolds = 5L, alphas = seq(0.05, 0.95, 0.05),
                      rule = c('lambda.min', 'lambda.1se'), seed = 1001) {

  enet_y = glmnet.tune.alpha(x, y, family = 'cox',
                             nfolds = nfolds, alphas = alphas,
                             seed = seed)

  # fit the model on all the data use the parameters got by CV
  best_alpha_enet = enet_y$best.alpha

  if (rule == 'lambda.min') {
    best_lambda_enet = enet_y$best.model$lambda.min
  } else if (rule == 'lambda.1se') {
    best_lambda_enet = enet_y$best.model$lambda.1se
  }

  enet_all = glmnet(x, y, family = 'cox',
                    lambda = best_lambda_enet,
                    alpha  = best_alpha_enet)

  coxenet_model = list('enet_best_alpha' = best_alpha_enet,
                       'enet_best_lambda' = best_lambda_enet,
                       'enet_model' = enet_all)

  class(coxenet_model) = 'hdcox.model.enet'

  return(coxenet_model)

}

#' Lasso model selection for high-dimensional Cox models
#'
#' Automatic lasso model selection for high-dimensional
#' Cox models, evaluated by penalized partial-likelihood.
#'
#' @param x Data matrix.
#' @param y Response matrix made by \code{\link[survival]{Surv}}.
#' @param nfolds Fold numbers of cross-validation.
#' @param rule Model selection criterion, \code{"lambda.min"} or
#' \code{"lambda.1se"}. See \code{\link[glmnet]{cv.glmnet}}
#' for details.
#' @param seed A random seed for cross-validation fold division.
#'
#' @export hdcox.lasso
#'
#' @examples
#' library("glmnet")
#' library("survival")
#' library("rms")
#'
#' # Load imputed SMART data
#' data("smart")
#' x = as.matrix(smart[, -c(1, 2)])
#' time = smart$TEVENT
#' event = smart$EVENT
#' y = Surv(time, event)
#'
#' # Fit Cox model by adaptive lasso penalization
#' lassofit = hdcox.lasso(x, y, nfolds = 10, rule = "lambda.1se", seed = 11)
#'
#' # Prepare data for hdnom.nomogram
#' x.df = as.data.frame(x)
#' dd = datadist(x.df)
#' options(datadist = "dd")
#'
#' # Generate hdnom.nomogram objects and plot nomogram
#' nom = hdnom.nomogram(lassofit$lasso_model, x, time, event, x.df,
#'                      lambda = lassofit$lasso_best_lambda, pred.at = 365 * 2,
#'                      funlabel = "2-Year Overall Survival Probability")
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

  coxlasso_model = list('lasso_best_lambda' = best_lambda_lasso,
                        'lasso_model' = lasso_all)

  class(coxlasso_model) = 'hdcox.model.lasso'

  return(coxlasso_model)

}
