#' Nomograms for Cox models fitted with glmnet
#'
#' Nomograms for Cox models fitted with glmnet
#'
#' @param Fitted \code{"glmnet"} model object.
#' @param x Matrix of training data used for the \code{glmnet} object.
#' @param ddist Data frame version of x, made by \code{datadist()}
#' @param s Value of the penalty parameter lambda in \code{glmnet}.
#' We will use the selected variables at the provided \code{s} to
#' build the nomogram, and make predictions.
#' See the example for choosing a proper lambda value extracted
#' from cross-validation results.
#' @param ... other arguments for \code{nomogram}.
#'
#' @export glmnet.nomogram
#'
#' @importFrom rms ols nomogram
#'
#' @examples
#' TBA
glmnet.nomogram = function(object, x, ddist, s, ...) {

  if (!all(c('coxnet', 'glmnet') %in% class(object)))
    stop('object class must be "glmnet" and "coxnet"')

  glmnet_pred_lp = as.vector(predict(object, newx = x, s = s, type = 'link'))

  all_vars = rownames(object$beta)
  selected_vars = all_vars[which(abs(coef(fit, s = s)) > .Machine$double.eps)]
  ols_formula = paste('glmnet_pred_lp ~',
                      paste(selected_vars, collapse = ' + '))
  ols_fit = ols(as.formula(ols_formula), data = ddist,
                sigma = 1, x = TRUE, y = TRUE)

  nom = nomogram(ols_fit, ...)

}

library('glmnet')
library('survival')
library('rms')

# Original SMART data not truncated
smart = read.table('~/Desktop/hdnomo-support/smartdata/SMARTs.tsv', header = TRUE, sep = '\t')
smart = na.omit(smart)
x = as.matrix(smart[, -c(1, 2)])
time = smart[, 1]
status = smart[, 2]

x.df = as.data.frame(x)
dd = datadist(x.df)
options(datadist = "dd")

set.seed(1010)
cvfit = cv.glmnet(x, Surv(time, status), family = "cox", nfolds = 5)
fit = glmnet(x, Surv(time, status), family = "cox")
nom = glmnet.nomogram(fit, x, ddist = x.df, s = cvfit$lambda.min)
plot(nom)
