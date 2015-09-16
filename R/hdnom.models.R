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
glmnet.tune.alpha = function(..., alphas, seed, parallel) {

  if (!parallel) {
    modelList = vector('list', length(alphas))
    for (i in 1L:length(alphas)) {
      set.seed(seed)
      modelList[[i]] = cv.glmnet(..., alpha = alphas[i])
    }
  } else {
    modelList <- foreach(alphas = alphas) %dopar% {
      set.seed(seed)
      cv.glmnet(..., alpha = alphas)
    }
  }

  # Choose model for best lambda first (then alpha)
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
#' @param parallel Logical. Enable parallel parameter tuning or not,
#' default is {FALSE}. To enable parallel tuning, load the
#' \code{doParallel} package and run \code{registerDoParallel()}
#' with the number of CPU cores before calling this function.
#'
#' @importFrom glmnet glmnet
#'
#' @export hdcox.aenet
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
#' # To enable parallel parameter tuning, first run:
#' # library("doParallel")
#' # registerDoParallel(detectCores())
#' # then set hdcox.aenet(..., parallel = TRUE).
#'
#' # Fit Cox model by adaptive elastic-net penalization
#' aenetfit = hdcox.aenet(x, y, nfolds = 3, alphas = c(0.3, 0.7),
#'                        rule = "lambda.1se", seeds = c(5, 7))
#'
#' # Prepare data for hdnom.nomogram
#' x.df = as.data.frame(x)
#' dd = datadist(x.df)
#' options(datadist = "dd")
#'
#' # Generate hdnom.nomogram objects and plot nomogram
#' nom = hdnom.nomogram(aenetfit$aenet_model, model.type = 'aenet',
#'                      x, time, event, x.df,
#'                      lambda = aenetfit$aenet_best_lambda, pred.at = 365 * 2,
#'                      funlabel = "2-Year Overall Survival Probability")
#'
#' plot(nom)
hdcox.aenet = function(x, y, nfolds = 5L, alphas = seq(0.05, 0.95, 0.05),
                       rule = c('lambda.min', 'lambda.1se'),
                       seeds = c(1001, 1002),
                       parallel = FALSE) {

  rule = match.arg(rule)

  # Tuning alpha for the both two stages of adaptive enet estimation
  enet_y = glmnet.tune.alpha(x, y, family = 'cox',
                             nfolds = nfolds, alphas = alphas,
                             seed = seeds[1L], parallel = parallel)

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
                              seed = seeds[2L],
                              parallel = parallel)

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
#'
#' # Load imputed SMART data
#' data("smart")
#' x = as.matrix(smart[, -c(1, 2)])
#' time = smart$TEVENT
#' event = smart$EVENT
#' y = Surv(time, event)
#'
#' # Fit Cox model by adaptive lasso penalization
#' alassofit = hdcox.alasso(x, y, nfolds = 3, rule = "lambda.1se", seeds = c(7, 11))
#'
#' # Prepare data for hdnom.nomogram
#' x.df = as.data.frame(x)
#' dd = datadist(x.df)
#' options(datadist = "dd")
#'
#' # Generate hdnom.nomogram objects and plot nomogram
#' nom = hdnom.nomogram(alassofit$alasso_model, model.type = 'alasso',
#'                      x, time, event, x.df,
#'                      lambda = alassofit$alasso_best_lambda, pred.at = 365 * 2,
#'                      funlabel = "2-Year Overall Survival Probability")
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
#' @param parallel Logical. Enable parallel parameter tuning or not,
#' default is {FALSE}. To enable parallel tuning, load the
#' \code{doParallel} package and run \code{registerDoParallel()}
#' with the number of CPU cores before calling this function.
#'
#' @export hdcox.enet
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
#' # To enable parallel parameter tuning, first run:
#' # library("doParallel")
#' # registerDoParallel(detectCores())
#' # then set hdcox.enet(..., parallel = TRUE).
#'
#' # Fit Cox model by adaptive elastic-net penalization
#' enetfit = hdcox.enet(x, y, nfolds = 3, alphas = c(0.3, 0.7),
#'                      rule = "lambda.1se", seed = 11)
#'
#' # Prepare data for hdnom.nomogram
#' x.df = as.data.frame(x)
#' dd = datadist(x.df)
#' options(datadist = "dd")
#'
#' # Generate hdnom.nomogram objects and plot nomogram
#' nom = hdnom.nomogram(enetfit$enet_model, model.type = 'enet',
#'                      x, time, event, x.df,
#'                      lambda = enetfit$enet_best_lambda, pred.at = 365 * 2,
#'                      funlabel = "2-Year Overall Survival Probability")
#'
#' plot(nom)
hdcox.enet = function(x, y, nfolds = 5L, alphas = seq(0.05, 0.95, 0.05),
                      rule = c('lambda.min', 'lambda.1se'),
                      seed = 1001, parallel = FALSE) {

  enet_y = glmnet.tune.alpha(x, y, family = 'cox',
                             nfolds = nfolds, alphas = alphas,
                             seed = seed, parallel = parallel)

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
#' lassofit = hdcox.lasso(x, y, nfolds = 5, rule = "lambda.1se", seed = 11)
#'
#' # Prepare data for hdnom.nomogram
#' x.df = as.data.frame(x)
#' dd = datadist(x.df)
#' options(datadist = "dd")
#'
#' # Generate hdnom.nomogram objects and plot nomogram
#' nom = hdnom.nomogram(lassofit$lasso_model, model.type = 'lasso',
#'                      x, time, event, x.df,
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

#' Fused lasso model selection for high-dimensional Cox models
#'
#' Automatic fused lasso model selection for high-dimensional
#' Cox models, evaluated by penalized partial-likelihood.
#'
#' @param x Data matrix.
#' @param y Response matrix made by \code{\link[survival]{Surv}}.
#' @param nfolds Fold numbers of cross-validation.
#' @param seed A random seed for cross-validation fold division.
#' @param trace Output the cross-validation parameter tuning
#' progress or not. Default is \code{FALSE}.
#'
#' @note The cross-validation procedure used in this function is the
#' so-called approximated cross-validation provided by the \code{penalized}
#' package. Be careful dealing with the results since they might be more
#' optimistic than fully fitting the model. This cross-validation method
#' is more suitable for datasets with larger number of observations,
#' and a higher number of cross-validation folds.
#'
#' @importFrom penalized optL1
#'
#' @export hdcox.flasso
#'
#' @examples
#' library("penalized")
#' library("survival")
#' library("rms")
#'
#' # Load imputed SMART data; only use the first 120 samples
#' data("smart")
#' x = as.matrix(smart[, -c(1, 2)])[1:120, ]
#' time = smart$TEVENT[1:120]
#' event = smart$EVENT[1:120]
#' y = Surv(time, event)
#'
#' # Fit Cox model by fused lasso penalization
#' flassofit = hdcox.flasso(x, y, nfolds = 3, seed = 11, trace = TRUE)
#'
#' # Prepare data for hdnom.nomogram
#' x.df = as.data.frame(x)
#' dd = datadist(x.df)
#' options(datadist = "dd")
#'
#' # Generate hdnom.nomogram objects and plot nomogram
#' nom = hdnom.nomogram(flassofit$flasso_model, model.type = 'flasso',
#'                      x, time, event, x.df,
#'                      lambda = flassofit$flasso_best_lambda, pred.at = 365 * 2,
#'                      funlabel = "2-Year Overall Survival Probability")
#'
#' plot(nom)
hdcox.flasso = function(x, y, nfolds = 5L,
                        seed = 1001, trace = FALSE) {

  set.seed(seed)
  flasso_all = optL1(response = y, penalized = x, fusedl = TRUE,
                     standardize = TRUE, model = 'cox', fold = nfolds,
                     trace = trace)

  coxflasso_model = list('flasso_best_lambda' = flasso_all$lambda,
                         'flasso_model' = flasso_all$fullfit)

  class(coxflasso_model) = 'hdcox.model.flasso'

  return(coxflasso_model)

}

#' Automatic MCP/SCAD gamma tuning function by k-fold cross-validation
#'
#' @return best model object and best gamma
#'
#' @importFrom ncvreg cv.ncvsurv
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#'
#' @keywords internal
ncvreg.tune.gamma = function(..., gammas, seed, parallel) {

  if (!parallel) {
    modelList = vector('list', length(gammas))
    for (i in 1L:length(gammas)) {
      set.seed(seed)
      modelList[[i]] = cv.ncvsurv(..., gamma = gammas[i])
    }
  } else {
    modelList <- foreach(gammas = gammas) %dopar% {
      set.seed(seed)
      cv.ncvsurv(..., gamma = gammas)
    }
  }

  # Choose model for best lambda first (then gamma)
  # Criterion: penalized partial likelihood
  errors = unlist(lapply(modelList, function(x) min(sqrt(x$cve))))

  return(list('best.model' = modelList[[which.min(errors)]],
              'best.gamma' = gammas[which.min(errors)]))

}

#' MCP model selection for high-dimensional Cox models
#'
#' Automatic MCP model selection for high-dimensional
#' Cox models, evaluated by penalized partial-likelihood.
#'
#' @param x Data matrix.
#' @param y Response matrix made by \code{\link[survival]{Surv}}.
#' @param nfolds Fold numbers of cross-validation.
#' @param gammas Gammas to tune in \code{\link[ncvreg]{cv.ncvsurv}}.
#' @param seed A random seed for cross-validation fold division.
#' @param trace Output the cross-validation parameter tuning
#' progress or not. Default is \code{FALSE}.
#' @param parallel Logical. Enable parallel parameter tuning or not,
#' default is {FALSE}. To enable parallel tuning, load the
#' \code{doParallel} package and run \code{registerDoParallel()}
#' with the number of CPU cores before calling this function.
#'
#' @importFrom ncvreg ncvsurv
#'
#' @export hdcox.mcp
#'
#' @examples
#' library("ncvreg")
#' library("survival")
#' library("rms")
#'
#' # Load imputed SMART data; only use the first 120 samples
#' data("smart")
#' x = as.matrix(smart[, -c(1, 2)])[1:120, ]
#' time = smart$TEVENT[1:120]
#' event = smart$EVENT[1:120]
#' y = Surv(time, event)
#'
#' # Fit Cox model by MCP penalization
#' mcpfit = hdcox.mcp(x, y, nfolds = 3, gammas = c(2.1, 3))
#'
#' # Prepare data for hdnom.nomogram
#' x.df = as.data.frame(x)
#' dd = datadist(x.df)
#' options(datadist = "dd")
#'
#' # Generate hdnom.nomogram objects and plot nomogram
#' nom = hdnom.nomogram(mcpfit$mcp_model, model.type = 'mcp', x, time, event, x.df,
#'                      lambda = mcpfit$mcp_best_lambda, pred.at = 365 * 2,
#'                      funlabel = "2-Year Overall Survival Probability")
#' plot(nom)
hdcox.mcp = function(x, y, nfolds = 5L, gammas = c(1.01, 1.7, 3, 100),
                     seed = 1001, trace = FALSE, parallel = FALSE) {

  mcp_y = ncvreg.tune.gamma(x, y, model = 'cox', penalty = 'MCP', alpha = 1,
                            nfolds = nfolds, gammas = gammas, seed = seed,
                            trace = trace, parallel = parallel)

  mcp_best_gamma  = mcp_y$best.gamma
  mcp_best_lambda = mcp_y$best.model$lambda.min

  # fit the model on all the data use the parameters got by CV
  mcp_all = ncvsurv(x, y, model = 'cox', penalty = 'MCP', alpha = 1,
                    gamma = mcp_best_gamma, lambda = mcp_best_lambda)

  coxmcp_model = list('mcp_best_gamma' = mcp_best_gamma,
                      'mcp_best_lambda' = mcp_best_lambda,
                      'mcp_model' = mcp_all)

  class(coxmcp_model) = 'hdcox.model.mcp'

  return(coxmcp_model)

}

#' Automatic Mnet/Snet gamma and alpha tuning function by k-fold cross-validation
#'
#' @return best model object, best gamma, and best alpha
#'
#' @importFrom ncvreg cv.ncvsurv
#' @importFrom foreach %dopar%
#' @importFrom foreach %:%
#' @importFrom foreach foreach
#'
#' @keywords internal
ncvreg.tune.gamma.alpha = function(..., gammas, alphas, seed, parallel) {

  if (!parallel) {

    modelList = vector('list', length(gammas))
    for (k in 1L:length(modelList)) {
      modelList[[k]] = vector('list', length(alphas))
    }

    for (i in 1L:length(gammas)) {
      for (j in 1L:length(alphas)) {
        set.seed(seed)
        modelList[[i]][[j]] =
          cv.ncvsurv(..., gamma = gammas[i], alpha = alphas[j])
      }
    }

    simplemodelList = unlist(modelList, recursive = FALSE)

  } else {

    modelList <- foreach(gammas = gammas) %:%
      foreach(alphas = alphas) %dopar% {
        set.seed(seed)
        cv.ncvsurv(..., gamma = gammas, alpha = alphas)
      }

    simplemodelList = unlist(modelList, recursive = FALSE)

  }

  # Choose model for best lambda first (then gamma/alpha)
  # Criterion: penalized partial likelihood
  errors = unlist(lapply(simplemodelList,
                         function(x) min(sqrt(x$cve))))

  return(list('best.model' = simplemodelList[[which.min(errors)]],
              'best.gamma' = simplemodelList[[which.min(errors)]]$fit$gamma,
              'best.alpha' = simplemodelList[[which.min(errors)]]$fit$alpha))

}

#' Mnet model selection for high-dimensional Cox models
#'
#' Automatic Mnet model selection for high-dimensional
#' Cox models, evaluated by penalized partial-likelihood.
#'
#' @param x Data matrix.
#' @param y Response matrix made by \code{\link[survival]{Surv}}.
#' @param nfolds Fold numbers of cross-validation.
#' @param gammas Gammas to tune in \code{\link[ncvreg]{cv.ncvsurv}}.
#' @param alphas Alphas to tune in \code{\link[ncvreg]{cv.ncvsurv}}.
#' @param seed A random seed for cross-validation fold division.
#' @param trace Output the cross-validation parameter tuning
#' progress or not. Default is \code{FALSE}.
#' @param parallel Logical. Enable parallel parameter tuning or not,
#' default is {FALSE}. To enable parallel tuning, load the
#' \code{doParallel} package and run \code{registerDoParallel()}
#' with the number of CPU cores before calling this function.
#'
#' @importFrom ncvreg ncvsurv
#'
#' @export hdcox.mnet
#'
#' @examples
#' library("ncvreg")
#' library("survival")
#' library("rms")
#'
#' # Load imputed SMART data; only use the first 120 samples
#' data("smart")
#' x = as.matrix(smart[, -c(1, 2)])[1:120, ]
#' time = smart$TEVENT[1:120]
#' event = smart$EVENT[1:120]
#' y = Surv(time, event)
#'
#' # Fit Cox model by Mnet penalization
#' mnetfit = hdcox.mnet(x, y, nfolds = 3, gammas = 3, alphas = c(0.3, 0.8))
#'
#' # Prepare data for hdnom.nomogram
#' x.df = as.data.frame(x)
#' dd = datadist(x.df)
#' options(datadist = "dd")
#'
#' # Generate hdnom.nomogram objects and plot nomogram
#' nom = hdnom.nomogram(mnetfit$mnet_model, model.type = 'mnet',
#'                      x, time, event, x.df,
#'                      lambda = mnetfit$mnet_best_lambda,
#'                      pred.at = 365 * 2,
#'                      funlabel = "2-Year Overall Survival Probability")
#'
#' plot(nom)
hdcox.mnet = function(x, y, nfolds = 5L, gammas = c(1.01, 1.7, 3, 100),
                      alphas = seq(0.05, 0.95, 0.05),
                      seed = 1001, trace = FALSE, parallel = FALSE) {

  mnet_y = ncvreg.tune.gamma.alpha(x, y, model = 'cox', penalty = 'MCP',
                                   nfolds = nfolds,
                                   gammas = gammas, alphas = alphas,
                                   seed = seed, trace = trace,
                                   parallel = parallel)

  mnet_best_gamma  = mnet_y$best.gamma
  mnet_best_alpha  = mnet_y$best.alpha
  mnet_best_lambda = mnet_y$best.model$lambda.min

  # fit the model on all the data use the parameters got by CV
  mnet_all = ncvsurv(x, y, model = 'cox', penalty = 'MCP',
                     gamma = mnet_best_gamma,
                     alpha = mnet_best_alpha,
                     lambda = mnet_best_lambda)

  coxmnet_model = list('mnet_best_gamma'  = mnet_best_gamma,
                       'mnet_best_alpha'  = mnet_best_alpha,
                       'mnet_best_lambda' = mnet_best_lambda,
                       'mnet_model'       = mnet_all)

  class(coxmnet_model) = 'hdcox.model.mnet'

  return(coxmnet_model)

}

#' SCAD model selection for high-dimensional Cox models
#'
#' Automatic SCAD model selection for high-dimensional
#' Cox models, evaluated by penalized partial-likelihood.
#'
#' @param x Data matrix.
#' @param y Response matrix made by \code{\link[survival]{Surv}}.
#' @param nfolds Fold numbers of cross-validation.
#' @param gammas Gammas to tune in \code{\link[ncvreg]{cv.ncvsurv}}.
#' @param seed A random seed for cross-validation fold division.
#' @param trace Output the cross-validation parameter tuning
#' progress or not. Default is \code{FALSE}.
#' @param parallel Logical. Enable parallel parameter tuning or not,
#' default is {FALSE}. To enable parallel tuning, load the
#' \code{doParallel} package and run \code{registerDoParallel()}
#' with the number of CPU cores before calling this function.
#'
#' @importFrom ncvreg ncvsurv
#'
#' @export hdcox.scad
#'
#' @examples
#' library("ncvreg")
#' library("survival")
#' library("rms")
#'
#' # Load imputed SMART data; only use the first 120 samples
#' data("smart")
#' x = as.matrix(smart[, -c(1, 2)])[1:120, ]
#' time = smart$TEVENT[1:120]
#' event = smart$EVENT[1:120]
#' y = Surv(time, event)
#'
#' # Fit Cox model by SCAD penalization
#' scadfit = hdcox.scad(x, y, nfolds = 3, gammas = c(3.7, 5))
#'
#' # Prepare data for hdnom.nomogram
#' x.df = as.data.frame(x)
#' dd = datadist(x.df)
#' options(datadist = "dd")
#'
#' # Generate hdnom.nomogram objects and plot nomogram
#' nom = hdnom.nomogram(scadfit$scad_model, model.type = 'scad', x, time, event, x.df,
#'                      lambda = scadfit$scad_best_lambda, pred.at = 365 * 2,
#'                      funlabel = "2-Year Overall Survival Probability")
#'
#' plot(nom)
hdcox.scad = function(x, y, nfolds = 5L, gammas = c(2.01, 2.3, 3.7, 200),
                      seed = 1001, trace = FALSE, parallel = FALSE) {

  scad_y = ncvreg.tune.gamma(x, y, model = 'cox', penalty = 'SCAD', alpha = 1,
                             nfolds = nfolds, gammas = gammas, seed = seed,
                             trace = trace, parallel = parallel)

  scad_best_gamma  = scad_y$best.gamma
  scad_best_lambda = scad_y$best.model$lambda.min

  # fit the model on all the data use the parameters got by CV
  scad_all = ncvsurv(x, y, model = 'cox', penalty = 'SCAD', alpha = 1,
                     gamma = scad_best_gamma, lambda = scad_best_lambda)

  coxscad_model = list('scad_best_gamma' = scad_best_gamma,
                       'scad_best_lambda' = scad_best_lambda,
                       'scad_model' = scad_all)

  class(coxscad_model) = 'hdcox.model.scad'

  return(coxscad_model)

}

#' Snet model selection for high-dimensional Cox models
#'
#' Automatic Snet model selection for high-dimensional
#' Cox models, evaluated by penalized partial-likelihood.
#'
#' @param x Data matrix.
#' @param y Response matrix made by \code{\link[survival]{Surv}}.
#' @param nfolds Fold numbers of cross-validation.
#' @param gammas Gammas to tune in \code{\link[ncvreg]{cv.ncvsurv}}.
#' @param alphas Alphas to tune in \code{\link[ncvreg]{cv.ncvsurv}}.
#' @param seed A random seed for cross-validation fold division.
#' @param trace Output the cross-validation parameter tuning
#' progress or not. Default is \code{FALSE}.
#' @param parallel Logical. Enable parallel parameter tuning or not,
#' default is {FALSE}. To enable parallel tuning, load the
#' \code{doParallel} package and run \code{registerDoParallel()}
#' with the number of CPU cores before calling this function.
#'
#' @importFrom ncvreg ncvsurv
#'
#' @export hdcox.snet
#'
#' @examples
#' library("ncvreg")
#' library("survival")
#' library("rms")
#'
#' # Load imputed SMART data; only use the first 120 samples
#' data("smart")
#' x = as.matrix(smart[, -c(1, 2)])[1:120, ]
#' time = smart$TEVENT[1:120]
#' event = smart$EVENT[1:120]
#' y = Surv(time, event)
#'
#' # Fit Cox model by Snet penalization
#' snetfit = hdcox.snet(x, y, nfolds = 3, gammas = 3.7, alphas = c(0.3, 0.8))
#'
#' # Prepare data for hdnom.nomogram
#' x.df = as.data.frame(x)
#' dd = datadist(x.df)
#' options(datadist = "dd")
#'
#' # Generate hdnom.nomogram objects and plot nomogram
#' nom = hdnom.nomogram(snetfit$snet_model, model.type = 'snet',
#'                      x, time, event, x.df,
#'                      lambda = snetfit$snet_best_lambda,
#'                      pred.at = 365 * 2,
#'                      funlabel = "2-Year Overall Survival Probability")
#'
#' plot(nom)
hdcox.snet = function(x, y, nfolds = 5L, gammas = c(2.01, 2.3, 3.7, 200),
                      alphas = seq(0.05, 0.95, 0.05),
                      seed = 1001, trace = FALSE, parallel = FALSE) {

  snet_y = ncvreg.tune.gamma.alpha(x, y, model = 'cox', penalty = 'SCAD',
                                   nfolds = nfolds,
                                   gammas = gammas, alphas = alphas,
                                   seed = seed, trace = trace,
                                   parallel = parallel)

  snet_best_gamma  = snet_y$best.gamma
  snet_best_alpha  = snet_y$best.alpha
  snet_best_lambda = snet_y$best.model$lambda.min

  # fit the model on all the data use the parameters got by CV
  snet_all = ncvsurv(x, y, model = 'cox', penalty = 'SCAD',
                     gamma = snet_best_gamma,
                     alpha = snet_best_alpha,
                     lambda = snet_best_lambda)

  coxsnet_model = list('snet_best_gamma'  = snet_best_gamma,
                       'snet_best_alpha'  = snet_best_alpha,
                       'snet_best_lambda' = snet_best_lambda,
                       'snet_model'       = snet_all)

  class(coxsnet_model) = 'hdcox.model.snet'

  return(coxsnet_model)

}
