#' Calibrate Cox models fitted with glmnet
#'
#' Calibrate Cox models fitted with glmnet
#'
#' @param x Matrix of training data used for the \code{glmnet} object;
#' on which to run the calibration.
#' @param time Survival time.
#' Must be of the same length with the number of rows as \code{x}.
#' @param event Status indicator, normally 0 = alive, 1 = dead.
#' Must be of the same length with the number of rows as \code{x}.
#' @param alpha Value of the elastic-net mixing parameter alpha in
#' glmnet. \code{alpha=1}: lasso; \code{alpha=0}: ridge.
#' From the Cox model you have built.
#' @param lambda Value of the penalty parameter lambda to use in the
#' glmnet fits on the resampled data. From the Cox model you have built.
#' @param method Calibration method.
#' Options including \code{"fitting"}, \code{"bootstrap"}, \code{"cv"},
#' and \code{"repeated.cv"}.
#' @param boot.times Number of repetitions for bootstrap.
#' @param nfolds Number of folds for cross-validation and
#' repeated cross-validation.
#' @param rep.times Number of repeated times for repeated cross-validation.
#' @param pred.at Time point at which calibration should take place.
#' @param ngroup number of groups to be formed for calibration
#'
#' @export glmnet.calibrate
#'
#' @examples
#' NULL
glmnet.calibrate = function (x, time, event,
                             alpha, lambda,
                             method = c('fitting', 'bootstrap', 'cv', 'repeated.cv'),
                             boot.times = NULL, nfolds = NULL, rep.times = NULL,
                             pred.at, ngroup = 5) {

  # It is better to use loess to estimate a smooth nonparametric calibration curve.
  cal = list()

  class(cal) = 'glmnet.calibrate'

  cal

}

#' Plot calibration curves
#'
#' Plot calibration curves
#'
#' @param x a \code{"glmnet.calibrate"} object.
#' @param ... other parameters for \code{plot}.
#'
#' @method plot glmnet.calibrate
#'
#' @export
#'
#' @examples
#' NULL
plot.glmnet.calibrate = function (x, ...) {

  if (class(x) != 'glmnet.calibrate')
    stop('object class must be "glmnet.calibrate"')

  cal = x

  plot(cal, ...)

}

#' Print calibration result generated with glmnet.calibrate
#'
#' Print calibration result generated with glmnet.calibrate
#'
#' @param x a \code{"glmnet.calibrate"} object.
#' @param ... other parameters (not used).
#'
#' @method print glmnet.calibrate
#'
#' @export
#'
#' @examples
#' NULL
print.glmnet.calibrate = function (x, ...) {

  if (class(x) != 'glmnet.calibrate')
    stop('object class must be "glmnet.calibrate"')

  cal = x

  print(cal, ...)

}
