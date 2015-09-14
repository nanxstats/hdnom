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
#' # Load imputed SMART data
#' data("smart")
#' x = as.matrix(smart[, -c(1, 2)])
#' time = smart$TEVENT
#' event = smart$EVENT
#' y = Surv(time, event)
#'
#' # Fit Cox model by MCP penalization
#' mcpfit = hdcox.mcp(x, y, nfolds = 3, gammas = c(1.7, 3), trace = TRUE)
#'
#' # Prepare data for hdnom.nomogram
#' x.df = as.data.frame(x)
#' dd = datadist(x.df)
#' options(datadist = "dd")
#'
#' # Generate hdnom.nomogram objects and plot nomogram
#' nom = hdnom.nomogram(mcpfit$mcp_model, x, time, event, x.df,
#'                      gamma = mcpfit$mcp_best_gamma,
#'                      lambda = mcpfit$mcp_best_lambda,
#'                      pred.at = 365 * 2,
#'                      funlabel = "2-Year Overall Survival Probability")
#'
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
#' # Load imputed SMART data
#' data("smart")
#' x = as.matrix(smart[, -c(1, 2)])
#' time = smart$TEVENT
#' event = smart$EVENT
#' y = Surv(time, event)
#'
#' # Fit Cox model by Mnet penalization
#' mnetfit = hdcox.mnet(x, y, nfolds = 3, gammas = c(1.7, 3),
#'                      alphas = c(0.3, 0.8), trace = TRUE)
#'
#' # Prepare data for hdnom.nomogram
#' x.df = as.data.frame(x)
#' dd = datadist(x.df)
#' options(datadist = "dd")
#'
#' # Generate hdnom.nomogram objects and plot nomogram
#' nom = hdnom.nomogram(mnetfit$mnet_model, x, time, event, x.df,
#'                      gamma = mnetfit$mnet_best_gamma,
#'                      alpha = mnetfit$mnet_best_alpha,
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
