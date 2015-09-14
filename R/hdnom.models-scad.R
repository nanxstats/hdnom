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
#' # Load imputed SMART data
#' data("smart")
#' x = as.matrix(smart[, -c(1, 2)])
#' time = smart$TEVENT
#' event = smart$EVENT
#' y = Surv(time, event)
#'
#' # Fit Cox model by SCAD penalization
#' scadfit = hdcox.scad(x, y, nfolds = 3, gammas = c(2.3, 3.7), trace = TRUE)
#'
#' # Prepare data for hdnom.nomogram
#' x.df = as.data.frame(x)
#' dd = datadist(x.df)
#' options(datadist = "dd")
#'
#' # Generate hdnom.nomogram objects and plot nomogram
#' nom = hdnom.nomogram(scadfit$scad_model, x, time, event, x.df,
#'                      gamma = scadfit$scad_best_gamma,
#'                      lambda = scadfit$scad_best_lambda,
#'                      pred.at = 365 * 2,
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
#' # Load imputed SMART data
#' data("smart")
#' x = as.matrix(smart[, -c(1, 2)])
#' time = smart$TEVENT
#' event = smart$EVENT
#' y = Surv(time, event)
#'
#' # Fit Cox model by Snet penalization
#' snetfit = hdcox.snet(x, y, nfolds = 3, gammas = c(2.01, 3.7),
#'                      alphas = c(0.3, 0.8), trace = TRUE)
#'
#' # Prepare data for hdnom.nomogram
#' x.df = as.data.frame(x)
#' dd = datadist(x.df)
#' options(datadist = "dd")
#'
#' # Generate hdnom.nomogram objects and plot nomogram
#' nom = hdnom.nomogram(snetfit$snet_model, x, time, event, x.df,
#'                      gamma = snetfit$snet_best_gamma,
#'                      alpha = snetfit$snet_best_alpha,
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
