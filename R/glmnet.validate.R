#' Validate Cox models fitted with glmnet using time-dependent AUC
#'
#' Validate Cox models fitted with glmnet using time-dependent AUC
#'
#' @param nomogram object returned by \code{\link{glmnet.nomogram}}.
#' @param x Matrix of training data used for the \code{glmnet} object;
#' on which to run the validation.
#' @param time Survival time.
#' Must be of the same length with the number of rows as \code{x}.
#' @param event Status indicator, normally 0 = alive, 1 = dead.
#' Must be of the same length with the number of rows as \code{x}.
#' @param s Value of the penalty parameter lambda in glmnet. By default,
#' we use the same penalty parameter as the nomogram object.
#' @param method Validation method. "bootstrap", "cv", or "repeated.cv".
#' @param B Number of repetitions for bootstrap.
#' @param nfolds Number of folds for cross-validation and
#' repeated cross-validation.
#' @param rep.times Number of repeated times for repeated cross-validation.
#'
#' @export glmnet.validate
#'
#' @examples
#' set.seed(1002)
glmnet.validate = function (nomogram, x, time, event, s = nomogram$s,
                            method = c('bootstrap', 'cv', 'repeated.cv'),
                            B = NULL, nfolds = NULL, rep.times = NULL) {

  method = match.arg(method)

  if (method == 'bootstrap') {

    if (!is.null(nfolds) || !is.null(rep.times))
      stop('nfolds and rep.times must be NULL when method = "bootstrap"')

    if (is.null(B)) B = 200L

    samp_mat = matrix(NA, nrow(x), B)
    for(i in 1L:B) samp_mat[, i] = sample(1L:nrow(x), replace = TRUE)



  } else if (method == 'cv') {

    if (!is.null(B) || !is.null(rep.times))
      stop('B and rep.times must be NULL when method = "cv"')

    if (is.null(nfolds)) nfolds = 5L

  } else if (method == 'repeated.cv') {

    if (!is.null(B)) stop('B must be NULL when method = "repeated.cv"')

    if (is.null(nfolds)) nfolds = 5L
    if (is.null(rep.times)) rep.times = 200L

  } else {
    stop('method must be one of "bootstrap", "cv", or "repeated.cv"')
  }

  val = list()

  class(val) = 'glmnet.validate'

  val

}

#' Compute validation measures for glmnet objects
#'
#' @return time-dependent AUC (tAUC) value
#'
#' @importFrom survival Surv
#' @importFrom glmnet glmnet
#'
#' @keywords internal
glmnet.validate.internal = function(samp_idx, nomogram, x, time, event, s) {

  samp_data = x[samp_idx, ]
  samp_fit = glmnet(samp_data, Surv(time, event), family = 'cox')
  lp.boot = samp_fit$x %*% as.matrix(samp_fit$coefficients) - mean(samp_fit$x %*% as.matrix(samp_fit$coefficients))
  lp.test = (nomogram$svy.cox)$x %*% as.matrix(samp_fit$coefficients) - mean((nomogram$svy.cox)$x %*% as.matrix(samp_fit$coefficients))
  cindex.train = 1 - rcorr.cens(lp.boot, Surv(samp_data$survival, samp_data$surv_cens))[[1]]
  cindex.test = 1 - rcorr.cens(lp.test, Surv(x$survival, x$surv_cens))[[1]]
  return(cindex.train - cindex.test)

}

#' Print validation result objects generated with glmnet.validate
#'
#' Print validation result objects generated with glmnet.validate
#'
#' @param object \code{"glmnet.validate"} object.
#'
#' @method print glmnet.validate
#'
#' @export
#'
#' @examples
#' NULL
print.glmnet.validate = function (object) {

  if (class(object) != 'glmnet.validate')
    stop('object class must be "glmnet.validate"')

  val = object

  print(val)

}

#' Plot validation result objects generated with glmnet.validate
#'
#' Plot validation result objects generated with glmnet.validate
#'
#' @param object \code{"glmnet.validate"} object.
#'
#' @method plot glmnet.validate
#'
#' @export
#'
#' @examples
#' NULL
plot.glmnet.validate = function (object) {

  if (class(object) != 'glmnet.validate')
    stop('object class must be "glmnet.validate"')

  val = object

  plot(val)

}
