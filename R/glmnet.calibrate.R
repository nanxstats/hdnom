#' Calibrate Cox models fitted with glmnet
#'
#' Calibrate Cox models fitted with glmnet
#'
#' @export glmnet.calibrate
#'
#' @examples
#' set.seed(1003)
glmnet.calibrate = function () {

  cal = list()

  class(cal) = 'glmnet.calibrate'

  cal

}

#' Plot calibration result objects generated with glmnet.calibrate
#'
#' Plot calibration result objects generated with glmnet.calibrate
#'
#' @param object \code{"glmnet.calibrate"} object.
#'
#' @method plot glmnet.calibrate
#'
#' @export
#'
#' @examples
#' NULL
plot.glmnet.calibrate = function (object) {

  if (class(object) != 'glmnet.calibrate')
    stop('object class must be "glmnet.calibrate"')

  cal = object

  plot(cal)

}
