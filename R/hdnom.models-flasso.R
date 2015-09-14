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
#' # Load imputed SMART data
#' data("smart")
#' x = as.matrix(smart[, -c(1, 2)])
#' time = smart$TEVENT
#' event = smart$EVENT
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
#' nom = hdnom.nomogram(flassofit$flasso_model, x, time, event, x.df,
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
