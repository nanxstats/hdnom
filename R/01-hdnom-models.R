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
    model_list = vector('list', length(alphas))
    for (i in 1L:length(alphas)) {
      set.seed(seed)
      model_list[[i]] = cv.glmnet(..., alpha = alphas[i])
    }
  } else {
    model_list <- foreach(alphas = alphas) %dopar% {
      set.seed(seed)
      cv.glmnet(..., alpha = alphas)
    }
  }

  # select model for best lambda first (then alpha)
  # criterion: penalized partial likelihood
  errors = unlist(lapply(model_list, function(x) min(sqrt(x$cvm))))

  list('best.model' = model_list[[which.min(errors)]],
       'best.alpha' = alphas[which.min(errors)])

}

#' Adaptive Elastic-Net Model Selection for High-Dimensional Cox Models
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
#' @param seed Two random seeds for cross-validation fold division
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
#' # Fit Cox model with adaptive elastic-net penalty
#' fit = hdcox.aenet(
#'   x, y, nfolds = 3, alphas = c(0.3, 0.7),
#'   rule = "lambda.1se", seed = c(5, 7))
#'
#' # Prepare data for hdnom.nomogram
#' x.df = as.data.frame(x)
#' dd = datadist(x.df)
#' options(datadist = "dd")
#'
#' # Generate hdnom.nomogram objects and plot nomogram
#' nom = hdnom.nomogram(
#'   fit$aenet_model, model.type = "aenet",
#'   x, time, event, x.df, pred.at = 365 * 2,
#'   funlabel = "2-Year Overall Survival Probability")
#'
#' plot(nom)
hdcox.aenet = function(
  x, y, nfolds = 5L, alphas = seq(0.05, 0.95, 0.05),
  rule = c('lambda.min', 'lambda.1se'),
  seed = c(1001, 1002),
  parallel = FALSE) {

  rule = match.arg(rule)

  # Tuning alpha for the both two stages of adaptive enet estimation
  enet_cv = glmnet.tune.alpha(
    x, y, family = 'cox',
    nfolds = nfolds, alphas = alphas,
    seed = seed[1L], parallel = parallel)

  # fit the model on all the data use the parameters got by CV
  best_alpha_enet = enet_cv$best.alpha

  if (rule == 'lambda.min') {
    best_lambda_enet = enet_cv$best.model$lambda.min
  } else if (rule == 'lambda.1se') {
    best_lambda_enet = enet_cv$best.model$lambda.1se
  }

  enet_full = glmnet(
    x, y, family = 'cox',
    lambda = best_lambda_enet,
    alpha  = best_alpha_enet)

  bhat = as.matrix(enet_full$beta)
  if(all(bhat == 0)) bhat = rep(.Machine$double.eps * 2, length(bhat))

  # adaptive penalty
  adpen = (1/pmax(abs(bhat), .Machine$double.eps))

  aenet_cv = glmnet.tune.alpha(
    x, y, family = 'cox', nfolds = nfolds,
    exclude = which(bhat == 0),
    penalty.factor = adpen,
    alphas = alphas,
    seed = seed[2L],
    parallel = parallel)

  # fit the model on all the data use the parameters got by CV
  best_alpha_aenet  = aenet_cv$best.alpha

  if (rule == 'lambda.min') {
    best_lambda_aenet = aenet_cv$best.model$lambda.min
  } else if (rule == 'lambda.1se') {
    best_lambda_aenet = aenet_cv$best.model$lambda.1se
  }

  aenet_full = glmnet(
    x, y, family = 'cox',
    exclude = which(bhat == 0),
    lambda  = best_lambda_aenet,
    penalty.factor = adpen,
    alpha   = best_alpha_aenet)

  if (aenet_full$df < 0.5)
    stop('Null model produced by the full fit (all coefficients are zero).
         Please try to tune rule, alphas, seed, nfolds, or increase sample size.')

  adpen_vec = as.vector(adpen)
  adpen_name = rownames(adpen)
  names(adpen_vec) = adpen_name

  coxaenet_model = list(
    'seed' = seed,
    'enet_best_alpha' = best_alpha_enet,
    'enet_best_lambda' = best_lambda_enet,
    'enet_model' = enet_full,
    'aenet_best_alpha' = best_alpha_aenet,
    'aenet_best_lambda' = best_lambda_aenet,
    'aenet_model' = aenet_full,
    'pen_factor' = adpen_vec)

  class(coxaenet_model) = c('hdcox.model', 'hdcox.model.aenet')

  coxaenet_model

}

#' Adaptive Lasso Model Selection for High-Dimensional Cox Models
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
#' @param seed Two random seeds for cross-validation fold division
#' in two estimation steps.
#'
#' @export hdcox.alasso
#'
#' @examples
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
#' # Fit Cox model with adaptive lasso penalty
#' fit = hdcox.alasso(x, y, nfolds = 3, rule = "lambda.1se", seed = c(7, 11))
#'
#' # Prepare data for hdnom.nomogram
#' x.df = as.data.frame(x)
#' dd = datadist(x.df)
#' options(datadist = "dd")
#'
#' # Generate hdnom.nomogram objects and plot nomogram
#' nom = hdnom.nomogram(
#'   fit$alasso_model, model.type = "alasso",
#'   x, time, event, x.df, pred.at = 365 * 2,
#'   funlabel = "2-Year Overall Survival Probability")
#'
#' plot(nom)
hdcox.alasso = function(
  x, y, nfolds = 5L,
  rule = c('lambda.min', 'lambda.1se'),
  seed = c(1001, 1002)) {

  # Tuning lambda for the both two stages of adaptive lasso estimation
  set.seed(seed[1L])
  lasso_cv = cv.glmnet(x, y, family = 'cox', nfolds = nfolds, alpha = 0)

  # fit the model on all the data use the parameters got by CV
  if (rule == 'lambda.min') {
    best_lambda_lasso = lasso_cv$lambda.min
  } else if (rule == 'lambda.1se') {
    best_lambda_lasso = lasso_cv$lambda.1se
  }

  lasso_full = glmnet(
    x, y, family = 'cox',
    lambda = best_lambda_lasso, alpha = 0)

  bhat = as.matrix(lasso_full$beta)
  if(all(bhat == 0)) bhat = rep(.Machine$double.eps * 2, length(bhat))

  # adaptive penalty
  adpen = (1/pmax(abs(bhat), .Machine$double.eps))

  set.seed(seed[2L])
  alasso_cv = cv.glmnet(
    x, y, family = 'cox', nfolds = nfolds, alpha = 1,
    penalty.factor = adpen)

  # fit the model on all the data use the parameters got by CV
  if (rule == 'lambda.min') {
    best_lambda_alasso = alasso_cv$lambda.min
  } else if (rule == 'lambda.1se') {
    best_lambda_alasso = alasso_cv$lambda.1se
  }

  alasso_full = glmnet(
    x, y, family = 'cox', lambda = best_lambda_alasso,
    alpha = 1, penalty.factor = adpen)

  if (alasso_full$df < 0.5)
    stop('Null model produced by the full fit (all coefficients are zero).
         Please try to tune rule, seed, nfolds, or increase sample size.')

  adpen_vec = as.vector(adpen)
  adpen_name = rownames(adpen)
  names(adpen_vec) = adpen_name

  coxalasso_model = list(
    'seed' = seed,
    'ridge_best_lambda' = best_lambda_lasso,
    'ridge_model' = lasso_full,
    'alasso_best_lambda' = best_lambda_alasso,
    'alasso_model' = alasso_full,
    'pen_factor' = adpen_vec)

  class(coxalasso_model) = c('hdcox.model', 'hdcox.model.alasso')

  coxalasso_model

}

#' Elastic-Net Model Selection for High-Dimensional Cox Models
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
#' # Fit Cox model with elastic-net penalty
#' fit = hdcox.enet(x, y, nfolds = 3, alphas = c(0.3, 0.7),
#'                  rule = "lambda.1se", seed = 11)
#'
#' # Prepare data for hdnom.nomogram
#' x.df = as.data.frame(x)
#' dd = datadist(x.df)
#' options(datadist = "dd")
#'
#' # Generate hdnom.nomogram objects and plot nomogram
#' nom = hdnom.nomogram(
#'   fit$enet_model, model.type = "enet",
#'   x, time, event, x.df, pred.at = 365 * 2,
#'   funlabel = "2-Year Overall Survival Probability")
#'
#' plot(nom)
hdcox.enet = function(
  x, y, nfolds = 5L, alphas = seq(0.05, 0.95, 0.05),
  rule = c('lambda.min', 'lambda.1se'),
  seed = 1001, parallel = FALSE) {

  enet_cv = glmnet.tune.alpha(
    x, y, family = 'cox',
    nfolds = nfolds, alphas = alphas,
    seed = seed, parallel = parallel)

  # fit the model on all the data use the parameters got by CV
  best_alpha_enet = enet_cv$best.alpha

  if (rule == 'lambda.min') {
    best_lambda_enet = enet_cv$best.model$lambda.min
  } else if (rule == 'lambda.1se') {
    best_lambda_enet = enet_cv$best.model$lambda.1se
  }

  enet_full = glmnet(
    x, y, family = 'cox',
    lambda = best_lambda_enet,
    alpha  = best_alpha_enet)

  if (enet_full$df < 0.5)
    stop('Null model produced by the full fit (all coefficients are zero).
         Please try to tune rule, alphas, seed, nfolds, or increase sample size.')

  coxenet_model = list(
    'seed' = seed,
    'enet_best_alpha' = best_alpha_enet,
    'enet_best_lambda' = best_lambda_enet,
    'enet_model' = enet_full)

  class(coxenet_model) = c('hdcox.model', 'hdcox.model.enet')

  coxenet_model

}

#' Lasso Model Selection for High-Dimensional Cox Models
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
#' # Fit Cox model with lasso penalty
#' fit = hdcox.lasso(x, y, nfolds = 5, rule = "lambda.1se", seed = 11)
#'
#' # Prepare data for hdnom.nomogram
#' x.df = as.data.frame(x)
#' dd = datadist(x.df)
#' options(datadist = "dd")
#'
#' # Generate hdnom.nomogram objects and plot nomogram
#' nom = hdnom.nomogram(
#'   fit$lasso_model, model.type = "lasso",
#'   x, time, event, x.df, pred.at = 365 * 2,
#'   funlabel = "2-Year Overall Survival Probability")
#'
#' plot(nom)
hdcox.lasso = function(
  x, y, nfolds = 5L,
  rule = c('lambda.min', 'lambda.1se'),
  seed = 1001) {

  set.seed(seed)
  lasso_cv = cv.glmnet(x, y, family = 'cox', nfolds = nfolds, alpha = 1)

  # fit the model on all the data use the parameters got by CV
  if (rule == 'lambda.min') {
    best_lambda_lasso = lasso_cv$lambda.min
  } else if (rule == 'lambda.1se') {
    best_lambda_lasso = lasso_cv$lambda.1se
  }

  lasso_full = glmnet(
    x, y, family = 'cox',
    lambda = best_lambda_lasso, alpha = 1)

  if (lasso_full$df < 0.5)
    stop('Null model produced by the full fit (all coefficients are zero).
         Please try to tune rule, seed, nfolds, or increase sample size.')

  coxlasso_model = list(
    'seed' = seed,
    'lasso_best_lambda' = best_lambda_lasso,
    'lasso_model' = lasso_full)

  class(coxlasso_model) = c('hdcox.model', 'hdcox.model.lasso')

  coxlasso_model

}

#' Automatic lambda tuning function for fused lasso by k-fold cross-validation
#'
#' @return best model object, best lambda1, and best lambda2
#'
#' @importFrom penalized cvl
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#'
#' @keywords internal
penalized.tune.lambda = function(..., lambda1, lambda2, seed, trace, parallel) {

  nlambda1 = length(lambda1)
  nlambda2 = length(lambda2)

  if (!parallel) {

    model_list = vector('list', nlambda1)
    for (k in 1L:nlambda1) model_list[[k]] = vector('list', nlambda2)

    for (i in 1L:nlambda1) {
      for (j in 1L:nlambda2) {
        set.seed(seed)
        if (trace) cat('Starting: lambda1 =', lambda1[i], 'lambda2 =', lambda2[j], '\n')
        model_list[[i]][[j]] =
          penalized::cvl(..., lambda1 = lambda1[i], lambda2 = lambda2[j])
      }
    }

    # store lambda combinations
    for (i in 1L:nlambda1) {
      for (j in 1L:nlambda2) {
        model_list[[i]][[j]][['lambda']] =
          c('lambda1' = lambda1[i], 'lambda2' = lambda2[j])
      }
    }

    simple_model_list = unlist(model_list, recursive = FALSE)

  } else {

    model_list <- foreach(lambda1 = lambda1) %:%
      foreach(lambda2 = lambda2) %dopar% {
        set.seed(seed)
        penalized::cvl(..., lambda1 = lambda1, lambda2 = lambda2)
      }

    # store lambda combinations
    for (i in 1L:nlambda1) {
      for (j in 1L:nlambda2) {
        model_list[[i]][[j]][['lambda']] =
          c('lambda1' = lambda1[i], 'lambda2' = lambda2[j])
      }
    }

    simple_model_list = unlist(model_list, recursive = FALSE)

  }

  # choose model for best lambda combination
  # criterion: cross-validated likelihood
  max_cvl = which.max(unlist(sapply(simple_model_list, '[', 'cvl')))

  list('best.model' = simple_model_list[[max_cvl]],
       'best.lambda1' = simple_model_list[[max_cvl]][['lambda']]['lambda1'],
       'best.lambda2' = simple_model_list[[max_cvl]][['lambda']]['lambda2'])

}

#' Fused Lasso Model Selection for High-Dimensional Cox Models
#'
#' Automatic fused lasso model selection for high-dimensional
#' Cox models, evaluated by cross-validated likelihood.
#'
#' @param x Data matrix.
#' @param y Response matrix made by \code{\link[survival]{Surv}}.
#' @param nfolds Fold numbers of cross-validation.
#' @param lambda1 Vector of lambda1 candidates.
#' Default is \code{0.001, 0.05, 0.5, 1, 5}.
#' @param lambda2 Vector of lambda2 candidates.
#' Default is \code{0.001, 0.01, 0.5}.
#' @param maxiter The maximum number of iterations allowed.
#' Default is \code{25}.
#' @param epsilon The convergence criterion.
#' Default is \code{1e-3}.
#' @param seed A random seed for cross-validation fold division.
#' @param trace Output the cross-validation parameter tuning
#' progress or not. Default is \code{FALSE}.
#' @param parallel Logical. Enable parallel parameter tuning or not,
#' default is {FALSE}. To enable parallel tuning, load the
#' \code{doParallel} package and run \code{registerDoParallel()}
#' with the number of CPU cores before calling this function.
#' @param ... other parameters to \code{\link[penalized]{cvl}}
#' and \code{\link[penalized]{penalized}}.
#'
#' @note The cross-validation procedure used in this function is the
#' \emph{approximated cross-validation} provided by the \code{penalized}
#' package. Be careful dealing with the results since they might be more
#' optimistic than a traditional CV procedure. This cross-validation
#' method is more suitable for datasets with larger number of observations,
#' and a higher number of cross-validation folds.
#'
#' @importFrom penalized penalized
#'
#' @export hdcox.flasso
#'
#' @examples
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
#' # Fit Cox model with fused lasso penalty
#' fit = hdcox.flasso(x, y,
#'   lambda1 = c(1, 10), lambda2 = c(0.01),
#'   nfolds = 3, seed = 11)
#'
#' # Prepare data for hdnom.nomogram
#' x.df = as.data.frame(x)
#' dd = datadist(x.df)
#' options(datadist = "dd")
#'
#' # Generate hdnom.nomogram objects and plot nomogram
#' nom = hdnom.nomogram(
#'   fit$flasso_model, model.type = "flasso",
#'   x, time, event, x.df, pred.at = 365 * 2,
#'   funlabel = "2-Year Overall Survival Probability")
#'
#' plot(nom)
hdcox.flasso = function(
  x, y, nfolds = 5L,
  lambda1 = c(0.001, 0.05, 0.5, 1, 5),
  lambda2 = c(0.001, 0.01, 0.5),
  maxiter = 25, epsilon = 1e-3,
  seed = 1001, trace = FALSE, parallel = FALSE, ...) {

  if (trace) cat('Starting cross-validation...\n')
  flasso_cv = penalized.tune.lambda(
    response = y, penalized = x, fold = nfolds,
    lambda1 = lambda1, lambda2 = lambda2,
    maxiter = maxiter, epsilon = epsilon,
    seed = seed, trace = trace, parallel = parallel,
    fusedl = TRUE, standardize = TRUE, model = 'cox', ...)

  flasso_best_lambda1  = flasso_cv$'best.lambda1'
  flasso_best_lambda2  = flasso_cv$'best.lambda2'

  # fit the model on all the data use the parameters got by CV
  if (trace) cat('Fitting fused lasso model with full data...\n')
  flasso_full = penalized(
    response = y, penalized = x,
    lambda1 = flasso_best_lambda1,
    lambda2 = flasso_best_lambda2,
    maxiter = maxiter, epsilon = epsilon,
    trace = trace,
    fusedl = TRUE, standardize = FALSE, model = 'cox', ...)

  if (all(abs(flasso_full@penalized) < .Machine$double.eps))
    stop('Null model produced by the full fit (all coefficients are zero).
         Please try changing the seed, nfolds, or increase sample size.')

  coxflasso_model = list(
    'seed' = seed,
    'flasso_best_lambda1' = flasso_best_lambda1,
    'flasso_best_lambda2' = flasso_best_lambda2,
    'flasso_model' = flasso_full)

  class(coxflasso_model) = c('hdcox.model', 'hdcox.model.flasso')

  coxflasso_model

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
ncvreg.tune.gamma = function(..., gammas, eps, max.iter, seed, parallel) {

  if (!parallel) {
    model_list = vector('list', length(gammas))
    for (i in 1L:length(gammas)) {
      set.seed(seed)
      model_list[[i]] = cv.ncvsurv(..., gamma = gammas[i],
                                   eps = eps, max.iter = max.iter)
    }
  } else {
    model_list <- foreach(gammas = gammas) %dopar% {
      set.seed(seed)
      cv.ncvsurv(..., gamma = gammas,
                 eps = eps, max.iter = max.iter)
    }
  }

  # select model for best lambda first (then gamma)
  # criterion: penalized partial likelihood
  errors = unlist(lapply(model_list, function(x) min(sqrt(x$cve))))

  list('best.model' = model_list[[which.min(errors)]],
       'best.gamma' = gammas[which.min(errors)])

}

#' MCP Model Selection for High-Dimensional Cox Models
#'
#' Automatic MCP model selection for high-dimensional
#' Cox models, evaluated by penalized partial-likelihood.
#'
#' @param x Data matrix.
#' @param y Response matrix made by \code{\link[survival]{Surv}}.
#' @param nfolds Fold numbers of cross-validation.
#' @param gammas Gammas to tune in \code{\link[ncvreg]{cv.ncvsurv}}.
#' @param eps Convergence threshhold.
#' @param max.iter Maximum number of iterations.
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
#' library("survival")
#' library("rms")
#'
#' # Load imputed SMART data; only use the first 150 samples
#' data("smart")
#' x = as.matrix(smart[, -c(1, 2)])[1:150, ]
#' time = smart$TEVENT[1:150]
#' event = smart$EVENT[1:150]
#' y = Surv(time, event)
#'
#' # Fit Cox model with MCP penalty
#' fit = hdcox.mcp(x, y, nfolds = 3, gammas = c(2.1, 3), seed = 1001)
#'
#' # Prepare data for hdnom.nomogram
#' x.df = as.data.frame(x)
#' dd = datadist(x.df)
#' options(datadist = "dd")
#'
#' # Generate hdnom.nomogram objects and plot nomogram
#' nom = hdnom.nomogram(
#'   fit$mcp_model, model.type = "mcp",
#'   x, time, event, x.df, pred.at = 365 * 2,
#'   funlabel = "2-Year Overall Survival Probability")
#'
#' plot(nom)
hdcox.mcp = function(
  x, y, nfolds = 5L, gammas = c(1.01, 1.7, 3, 100),
  eps = 1e-4, max.iter = 10000L,
  seed = 1001, trace = FALSE, parallel = FALSE) {

  mcp_cv = ncvreg.tune.gamma(
    x, y, penalty = 'MCP', alpha = 1,
    nfolds = nfolds, gammas = gammas,
    eps = eps, max.iter = max.iter,
    seed = seed, trace = trace, parallel = parallel)

  mcp_best_gamma  = mcp_cv$best.gamma
  mcp_best_lambda = mcp_cv$best.model$lambda.min

  # fit the model on all the data use the parameters got by CV
  mcp_full =
    ncvreg::ncvsurv(
      x, y, penalty = 'MCP', alpha = 1,
      gamma = mcp_best_gamma, lambda = mcp_best_lambda,
      eps = eps, max.iter = max.iter)

  # deal with null models, thanks for the suggestion from Patrick Breheny
  if (all(abs(mcp_full$beta[-1L, ]) < .Machine$double.eps))
    stop('Null model produced by the full fit (all coefficients are zero).
         Please try to tune gammas, seed, nfolds, or increase sample size.')

  coxmcp_model = list(
    'seed' = seed,
    'mcp_best_gamma' = mcp_best_gamma,
    'mcp_best_lambda' = mcp_best_lambda,
    'mcp_model' = mcp_full)

  class(coxmcp_model) = c('hdcox.model', 'hdcox.model.mcp')

  coxmcp_model

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
ncvreg.tune.gamma.alpha = function(
  ..., gammas, alphas, eps, max.iter, seed, parallel) {

  if (!parallel) {

    model_list = vector('list', length(gammas))
    for (k in 1L:length(model_list)) {
      model_list[[k]] = vector('list', length(alphas))
    }

    for (i in 1L:length(gammas)) {
      for (j in 1L:length(alphas)) {
        set.seed(seed)
        model_list[[i]][[j]] =
          cv.ncvsurv(..., gamma = gammas[i], alpha = alphas[j],
                     eps = eps, max.iter = max.iter)
      }
    }

    simple_model_list = unlist(model_list, recursive = FALSE)

  } else {

    model_list <- foreach(gammas = gammas) %:%
      foreach(alphas = alphas) %dopar% {
        set.seed(seed)
        cv.ncvsurv(..., gamma = gammas, alpha = alphas,
                   eps = eps, max.iter = max.iter)
      }

    simple_model_list = unlist(model_list, recursive = FALSE)

  }

  # select model for best lambda first (then gamma/alpha)
  # criterion: penalized partial likelihood
  errors = unlist(lapply(simple_model_list,
                         function(x) min(sqrt(x$cve))))

  list('best.model' = simple_model_list[[which.min(errors)]],
       'best.gamma' = simple_model_list[[which.min(errors)]]$fit$gamma,
       'best.alpha' = simple_model_list[[which.min(errors)]]$fit$alpha)

}

#' Mnet Model Selection for High-Dimensional Cox Models
#'
#' Automatic Mnet model selection for high-dimensional
#' Cox models, evaluated by penalized partial-likelihood.
#'
#' @param x Data matrix.
#' @param y Response matrix made by \code{\link[survival]{Surv}}.
#' @param nfolds Fold numbers of cross-validation.
#' @param gammas Gammas to tune in \code{\link[ncvreg]{cv.ncvsurv}}.
#' @param alphas Alphas to tune in \code{\link[ncvreg]{cv.ncvsurv}}.
#' @param eps Convergence threshhold.
#' @param max.iter Maximum number of iterations.
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
#' # Fit Cox model with Mnet penalty
#' fit = hdcox.mnet(
#'   x, y, nfolds = 3,
#'   gammas = 3, alphas = c(0.3, 0.8),
#'   max.iter = 15000, seed = 1010)
#'
#' # Prepare data for hdnom.nomogram
#' x.df = as.data.frame(x)
#' dd = datadist(x.df)
#' options(datadist = "dd")
#'
#' # Generate hdnom.nomogram objects and plot nomogram
#' nom = hdnom.nomogram(
#'   fit$mnet_model, model.type = "mnet",
#'   x, time, event, x.df, pred.at = 365 * 2,
#'   funlabel = "2-Year Overall Survival Probability")
#'
#' plot(nom)
hdcox.mnet = function(
  x, y, nfolds = 5L,
  gammas = c(1.01, 1.7, 3, 100),
  alphas = seq(0.05, 0.95, 0.05),
  eps = 1e-4, max.iter = 10000L,
  seed = 1001, trace = FALSE, parallel = FALSE) {

  mnet_cv = ncvreg.tune.gamma.alpha(
    x, y, penalty = 'MCP',
    nfolds = nfolds,
    gammas = gammas, alphas = alphas,
    eps = eps, max.iter = max.iter,
    seed = seed, trace = trace,
    parallel = parallel)

  mnet_best_gamma  = mnet_cv$best.gamma
  mnet_best_alpha  = mnet_cv$best.alpha
  mnet_best_lambda = mnet_cv$best.model$lambda.min

  # fit the model on all the data use the parameters got by CV
  mnet_full =
    ncvreg::ncvsurv(
      x, y, penalty = 'MCP',
      gamma = mnet_best_gamma,
      alpha = mnet_best_alpha,
      lambda = mnet_best_lambda,
      eps = eps, max.iter = max.iter)

  if (all(abs(mnet_full$beta[-1L, ]) < .Machine$double.eps))
    stop('Null model produced by the full fit (all coefficients are zero).
         Please try to tune gammas, alphas, seed, nfolds, or increase sample size.')

  coxmnet_model = list(
    'seed' = seed,
    'mnet_best_gamma'  = mnet_best_gamma,
    'mnet_best_alpha'  = mnet_best_alpha,
    'mnet_best_lambda' = mnet_best_lambda,
    'mnet_model'       = mnet_full)

  class(coxmnet_model) = c('hdcox.model', 'hdcox.model.mnet')

  coxmnet_model

}

#' SCAD Model Selection for High-Dimensional Cox Models
#'
#' Automatic SCAD model selection for high-dimensional
#' Cox models, evaluated by penalized partial-likelihood.
#'
#' @param x Data matrix.
#' @param y Response matrix made by \code{\link[survival]{Surv}}.
#' @param nfolds Fold numbers of cross-validation.
#' @param gammas Gammas to tune in \code{\link[ncvreg]{cv.ncvsurv}}.
#' @param eps Convergence threshhold.
#' @param max.iter Maximum number of iterations.
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
#' # Fit Cox model with SCAD penalty
#' fit = hdcox.scad(
#'   x, y, nfolds = 3, gammas = c(3.7, 5),
#'   max.iter = 15000, seed = 1010)
#'
#' # Prepare data for hdnom.nomogram
#' x.df = as.data.frame(x)
#' dd = datadist(x.df)
#' options(datadist = "dd")
#'
#' # Generate hdnom.nomogram objects and plot nomogram
#' nom = hdnom.nomogram(
#'   fit$scad_model, model.type = "scad",
#'   x, time, event, x.df, pred.at = 365 * 2,
#'   funlabel = "2-Year Overall Survival Probability")
#'
#' plot(nom)
hdcox.scad = function(
  x, y, nfolds = 5L,
  gammas = c(2.01, 2.3, 3.7, 200),
  eps = 1e-4, max.iter = 10000L,
  seed = 1001, trace = FALSE, parallel = FALSE) {

  scad_cv = ncvreg.tune.gamma(
    x, y, penalty = 'SCAD', alpha = 1,
    nfolds = nfolds, gammas = gammas,
    eps = eps, max.iter = max.iter,
    seed = seed, trace = trace, parallel = parallel)

  scad_best_gamma  = scad_cv$best.gamma
  scad_best_lambda = scad_cv$best.model$lambda.min

  # fit the model on all the data use the parameters got by CV
  scad_full = ncvreg::ncvsurv(
    x, y, penalty = 'SCAD', alpha = 1,
    gamma = scad_best_gamma, lambda = scad_best_lambda,
    eps = eps, max.iter = max.iter)

  if (all(abs(scad_full$beta[-1L, ]) < .Machine$double.eps))
    stop('Null model produced by the full fit (all coefficients are zero).
         Please try to tune gammas, seed, nfolds, or increase sample size.')

  coxscad_model = list(
    'seed' = seed,
    'scad_best_gamma' = scad_best_gamma,
    'scad_best_lambda' = scad_best_lambda,
    'scad_model' = scad_full)

  class(coxscad_model) = c('hdcox.model', 'hdcox.model.scad')

  coxscad_model

}

#' Snet Model Selection for High-Dimensional Cox Models
#'
#' Automatic Snet model selection for high-dimensional
#' Cox models, evaluated by penalized partial-likelihood.
#'
#' @param x Data matrix.
#' @param y Response matrix made by \code{\link[survival]{Surv}}.
#' @param nfolds Fold numbers of cross-validation.
#' @param gammas Gammas to tune in \code{\link[ncvreg]{cv.ncvsurv}}.
#' @param alphas Alphas to tune in \code{\link[ncvreg]{cv.ncvsurv}}.
#' @param eps Convergence threshhold.
#' @param max.iter Maximum number of iterations.
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
#' # Fit Cox model with Snet penalty
#' fit = hdcox.snet(
#'   x, y, nfolds = 3,
#'   gammas = 3.7, alphas = c(0.3, 0.8),
#'   max.iter = 15000, seed = 1010)
#'
#' # Prepare data for hdnom.nomogram
#' x.df = as.data.frame(x)
#' dd = datadist(x.df)
#' options(datadist = "dd")
#'
#' # Generate hdnom.nomogram objects and plot nomogram
#' nom = hdnom.nomogram(
#'   fit$snet_model, model.type = "snet",
#'   x, time, event, x.df, pred.at = 365 * 2,
#'   funlabel = "2-Year Overall Survival Probability")
#'
#' plot(nom)
hdcox.snet = function(
  x, y, nfolds = 5L,
  gammas = c(2.01, 2.3, 3.7, 200),
  alphas = seq(0.05, 0.95, 0.05),
  eps = 1e-4, max.iter = 10000L,
  seed = 1001, trace = FALSE, parallel = FALSE) {

  snet_cv = ncvreg.tune.gamma.alpha(
    x, y, penalty = 'SCAD',
    nfolds = nfolds,
    gammas = gammas, alphas = alphas,
    eps = eps, max.iter = max.iter,
    seed = seed, trace = trace,
    parallel = parallel)

  snet_best_gamma  = snet_cv$best.gamma
  snet_best_alpha  = snet_cv$best.alpha
  snet_best_lambda = snet_cv$best.model$lambda.min

  # fit the model on all the data use the parameters got by CV
  snet_full = ncvreg::ncvsurv(
    x, y, penalty = 'SCAD',
    gamma = snet_best_gamma,
    alpha = snet_best_alpha,
    lambda = snet_best_lambda,
    eps = eps, max.iter = max.iter)

  if (all(abs(snet_full$beta[-1L, ]) < .Machine$double.eps))
    stop('Null model produced by the full fit (all coefficients are zero).
         Please try to tune gammas, alphas, seed, nfolds, or increase sample size.')

  coxsnet_model = list(
    'seed' = seed,
    'snet_best_gamma'  = snet_best_gamma,
    'snet_best_alpha'  = snet_best_alpha,
    'snet_best_lambda' = snet_best_lambda,
    'snet_model'       = snet_full)

  class(coxsnet_model) = c('hdcox.model', 'hdcox.model.snet')

  coxsnet_model

}



#' Make Predictions from High-Dimensional Cox Models
#'
#' Predict overall survival probability at certain time points
#' from established Cox models.
#'
#' @param object Model object fitted by \code{hdcox.*()} functions.
#' @param x Data matrix used to fit the model.
#' @param y Response matrix made with \code{\link[survival]{Surv}}.
#' @param newx Matrix (with named columns) of new values for \code{x}
#' at which predictions are to be made.
#' @param pred.at Time point at which prediction should take place.
#' @param ... Other parameters (not used).
#'
#' @return A \code{nrow(newx) x length(pred.at)} matrix containing
#' overall survival probablity.
#'
#' @method predict hdcox.model
#'
#' @importFrom penalized survival
#' @importMethodsFrom penalized predict
#'
#' @export
#'
#' @examples
#' library("survival")
#'
#' # Load imputed SMART data
#' data("smart")
#' x = as.matrix(smart[, -c(1, 2)])
#' time = smart$TEVENT
#' event = smart$EVENT
#' y = Surv(time, event)
#'
#' fit = hdcox.lasso(x, y, nfolds = 5, rule = "lambda.1se", seed = 11)
#' predict(fit, x, y, newx = x[101:105, ], pred.at = 1:10 * 365)
predict.hdcox.model = function(object, x, y, newx, pred.at, ...) {

  if (!('hdcox.model' %in% class(object)))
    stop('object must be of class "hdcox.model" fitted by hdcox.* functions')

  model.type = gsub('hdcox.model.', '', setdiff(class(object), 'hdcox.model'))

  if (!('matrix' %in% class(newx))) stop('newx must be a matrix')

  time  = y[, 1L]
  event = y[, 2L]

  obj.type = switch(
    model.type,
    lasso  = 'glmnet', alasso = 'glmnet',
    enet   = 'glmnet', aenet  = 'glmnet',
    mcp    = 'ncvreg', mnet   = 'ncvreg',
    scad   = 'ncvreg', snet   = 'ncvreg',
    flasso = 'penalized'
  )

  switch(

    obj.type,

    glmnet = {
      lp = predict(object[[paste0(model.type, '_model')]], x, type = 'link')
      basesurv = glmnet.basesurv(time, event, lp, pred.at)
      lpnew = predict(object[[paste0(model.type, '_model')]], newx, type = 'link')
      p = exp(exp(lpnew) %*% -t(basesurv$'cumulative_base_hazard'))
    },

    ncvreg = {
      lp = predict(object[[paste0(model.type, '_model')]], x, type = 'link')
      basesurv = ncvreg.basesurv(time, event, lp, pred.at)
      lpnew = predict(object[[paste0(model.type, '_model')]], newx, type = 'link')
      p = exp(exp(lpnew) %*% -t(basesurv$'cumulative_base_hazard'))
      # # alternative method using ncvreg built-in prediction directly
      # # almost identical results, but sometimes produces NAs in practice
      # # e.g. pred.at = 1:10 * 365
      # survfun = predict(object[[paste0(model.type, '_model')]], newx, type = 'survival')
      # p = matrix(NA, nrow = nrow(newx), ncol = length(pred.at))
      # for (i in 1L:nrow(newx)) p[i, ] = sapply(pred.at, survfun[[i]])
    },

    penalized = {
      pred = predict(object[[paste0(model.type, '_model')]], newx)
      p = matrix(NA, nrow = nrow(newx), ncol = length(pred.at))
      for (i in 1L:length(pred.at)) p[, i] = survival(pred, time = pred.at[i])
    }

  )

  colnames(p) = as.character(pred.at)

  p

}

#' Extract Information of Selected Variables from High-Dimensional Cox Models
#'
#' Extract the names and type of selected variables from established
#' high-dimensional Cox models.
#'
#' @param object Model object fitted by \code{hdcox.*()} functions.
#' @param x Data matrix used to fit the model.
#'
#' @export hdnom.varinfo
#'
#' @return A list containing the index, name, type and range of the
#' selected variables.
#'
#' @examples
#' library("survival")
#'
#' # Load imputed SMART data
#' data("smart")
#' x = as.matrix(smart[, -c(1, 2)])
#' time = smart$TEVENT
#' event = smart$EVENT
#' y = Surv(time, event)
#'
#' # Fit Cox model with lasso penalty
#' fit = hdcox.lasso(x, y, nfolds = 5, rule = "lambda.1se", seed = 11)
#' hdnom.varinfo(fit, x)
hdnom.varinfo = function(object, x) {

  if (!('hdcox.model' %in% class(object)))
    stop('object must be of class "hdcox.model" fitted by hdcox.* functions')

  model.type = gsub('hdcox.model.', '', setdiff(class(object), 'hdcox.model'))

  obj.type = switch(
    model.type,
    lasso  = 'glmnet', alasso = 'glmnet',
    enet   = 'glmnet', aenet  = 'glmnet',
    mcp    = 'ncvreg', mnet   = 'ncvreg',
    scad   = 'ncvreg', snet   = 'ncvreg',
    flasso = 'penalized'
  )

  switch(

    obj.type,

    glmnet = {
      nonzero_idx = which(as.logical(abs(object[[paste0(model.type, '_model')]][['beta']]) > .Machine$double.eps))
      nonzero_var = rownames(object[[paste0(model.type, '_model')]][['beta']])[nonzero_idx]
    },

    ncvreg = {
      nonzero_idx = which(object[[paste0(model.type, '_model')]][['beta']][-1L, ] > .Machine$double.eps)
      nonzero_var = names(object[[paste0(model.type, '_model')]][['beta']][-1L, ])[nonzero_idx]
    },

    penalized = {
      nonzero_idx = which(object[[paste0(model.type, '_model')]]@'penalized' > .Machine$double.eps)
      nonzero_var = colnames(x)[nonzero_idx]
    }

  )

  varinfo = list('index' = NULL, 'name' = NULL, 'type' = NULL, 'domain' = NULL)

  varinfo[['index']] = nonzero_idx
  varinfo[['name']]  = nonzero_var
  nvar = length(varinfo[['name']])

  is.wholenumber = function(
    x, tol = .Machine$double.eps^0.5) abs(x - round(x)) < tol

  # type: logical, categorical or continuous
  for (i in 1L:nvar) {
    var = x[, varinfo[['name']][i]]
    if (all(is.wholenumber(var)) & nlevels(as.factor(var)) == 2L) {
      varinfo[['type']][i] = 'logical'
      varinfo[['domain']][[i]] = unique(var)
    } else if (all(is.wholenumber(var)) & nlevels(as.factor(var)) > 2L) {
      varinfo[['type']][i] = 'categorical'
      varinfo[['domain']][[i]] = c(min(var), max(var))
    } else if (any(!is.wholenumber(var))) {
      varinfo[['type']][i] = 'continuous'
      varinfo[['domain']][[i]] = c(min(var), max(var))
    } else {
      stop(paste0('unrecognized variable type: ', varinfo[['name']][i]))
    }
  }

  class(varinfo) = 'hdnom.varinfo'

  varinfo

}

#' Print High-Dimensional Cox Model Objects
#'
#' Print information about high-dimensional Cox models.
#'
#' @param x Model object fitted by \code{hdcox.*()} functions.
#' @param ... Other parameters (not used).
#'
#' @method print hdcox.model
#'
#' @export
#'
#' @examples
#' library("survival")
#'
#' # Load imputed SMART data
#' data("smart")
#' x = as.matrix(smart[, -c(1, 2)])
#' time = smart$TEVENT
#' event = smart$EVENT
#' y = Surv(time, event)
#'
#' fit = hdcox.lasso(x, y, nfolds = 5, rule = "lambda.1se", seed = 11)
#' print(fit)
print.hdcox.model = function(x, ...) {

  if (!('hdcox.model' %in% class(x)))
    stop('x must be of class "hdcox.model" fitted by hdcox.* functions')

  model.type = gsub('hdcox.model.', '', setdiff(class(x), 'hdcox.model'))

  switch(

    model.type,

    lasso = {
      cat('High-Dimensional Cox Model Object\n')
      cat('Random seed:', x$'seed', '\n')
      cat('Model type: lasso\n')
      cat('Best lambda:', x$'lasso_best_lambda', '\n')
    },

    alasso = {
      cat('High-Dimensional Cox Model Object\n')
      cat('Random seed:', x$'seed', '\n')
      cat('Model type: adaptive lasso\n')
      cat('First step best lambda:', x$'ridge_best_lambda', '\n')
      cat('Second step best lambda:', x$'alasso_best_lambda', '\n')
    },

    enet = {
      cat('High-Dimensional Cox Model Object\n')
      cat('Random seed:', x$'seed', '\n')
      cat('Model type: elastic-net\n')
      cat('Best alpha:', x$'enet_best_alpha', '\n')
      cat('Best lambda:', x$'enet_best_lambda', '\n')
    },

    aenet = {
      cat('High-Dimensional Cox Model Object\n')
      cat('Random seed:', x$'seed', '\n')
      cat('Model type: adaptive elastic-net\n')
      cat('First step best alpha:', x$'enet_best_alpha', '\n')
      cat('First step best lambda:', x$'enet_best_lambda', '\n')
      cat('Second step best alpha:', x$'aenet_best_alpha', '\n')
      cat('Second step best lambda:', x$'aenet_best_lambda', '\n')
    },

    mcp = {
      cat('High-Dimensional Cox Model Object\n')
      cat('Random seed:', x$'seed', '\n')
      cat('Model type: MCP\n')
      cat('Best gamma:', x$'mcp_best_gamma', '\n')
      cat('Best lambda:', x$'mcp_best_lambda', '\n')
    },

    mnet = {
      cat('High-Dimensional Cox Model Object\n')
      cat('Random seed:', x$'seed', '\n')
      cat('Model type: Mnet\n')
      cat('Best gamma:', x$'mnet_best_gamma', '\n')
      cat('Best alpha:', x$'mnet_best_alpha', '\n')
      cat('Best lambda:', x$'mnet_best_lambda', '\n')
    },

    scad = {
      cat('High-Dimensional Cox Model Object\n')
      cat('Random seed:', x$'seed', '\n')
      cat('Model type: SCAD\n')
      cat('Best gamma:', x$'scad_best_gamma', '\n')
      cat('Best lambda:', x$'scad_best_lambda', '\n')
    },

    snet = {
      cat('High-Dimensional Cox Model Object\n')
      cat('Random seed:', x$'seed', '\n')
      cat('Model type: Snet\n')
      cat('Best gamma:', x$'snet_best_gamma', '\n')
      cat('Best alpha:', x$'snet_best_alpha', '\n')
      cat('Best lambda:', x$'snet_best_lambda', '\n')
    },

    flasso = {
      cat('High-Dimensional Cox Model Object\n')
      cat('Random seed:', x$'seed', '\n')
      cat('Model type: fused lasso\n')
      cat('Best lambda1:', x$'flasso_best_lambda1', '\n')
      cat('Best lambda2:', x$'flasso_best_lambda2', '\n')
    }

  )

}
