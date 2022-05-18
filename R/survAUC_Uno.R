#' AUC estimator proposed by Uno et al.
#'
#' Uno's estimator of cumulative/dynamic AUC for right-censored time-to-event data
#'
#' @param Surv.rsp A \code{Surv(.,.)} object containing to the outcome of the training data.
#' @param Surv.rsp.new A \code{Surv(.,.)} object containing the outcome of the test data.
#' @param lpnew The vector of predictors obtained from the test data.
#' @param times A vector of time points at which to evaluate AUC.
#' @param savesensspec A logical specifying whether sensitivities and specificities
#' should be saved.
#'
#' @details The \code{sens.uno} and \code{spec.uno} functions implement the estimators of
#' time-dependent true and false positive rates proposed in Section 5.1 of Uno et al. (2007).
#'
#' The \code{AUC.uno} function implements the estimator of cumulative/dynamic AUC
#' that is based on the TPR and FPR estimators proposed by
#' Uno et al. (2007). It is given by the area(s) under the time-dependent
#' ROC curve(s) estimated by \code{sens.sh} and \code{spec.sh}. The \code{iauc}
#' summary measure is given by the integral of AUC on
#' [0, max(\code{times})] (weighted by the estimated probability density of
#' the time-to-event outcome).
#'
#' Uno's estimators are based on inverse-probability-of-censoring
#' weights and do not assume a specific working model for deriving the predictor
#' \code{lpnew}. It is assumed, however, that there is a one-to-one
#' relationship between the predictor and the expected survival times conditional
#' on the predictor. Note that the estimators implemented in \code{sens.uno},
#' \code{spec.uno} and \code{AUC.uno} are restricted to situations
#' where the random censoring assumption holds.
#'
#' @returns \code{AUC.uno} returns an object of class \code{survAUC}. Specifically,
#' \code{AUC.uno} returns a list with the following components:
#' \item{auc}{The cumulative/dynamic AUC estimates (evaluated at \code{times}).}
#' \item{times}{The vector of time points at which AUC is evaluated.}
#' \item{iauc}{The summary measure of AUC.}
#' \code{sens.uno} and \code{spec.uno} return matrices of dimensions \code{times} x
#' \code{(lpnew + 1)}. The elements of these matrices are the sensitivity and
#' specificity estimates for each threshold of \code{lpnew} and for each time point
#' specified in \code{times}.
#'
#' @references Uno, H., T. Cai, L. Tian, and L. J. Wei (2007).
#' Evaluating prediction rules for t-year survivors with censored regression models.
#' \emph{Journal of the American Statistical Association} \bold{102}, 527--537.
#'
#' @examples
#' library(survival)
#'
#' TR <- ovarian[1:16, ]
#' TE <- ovarian[17:26, ]
#' train.fit <- coxph(Surv(futime, fustat) ~ age, x = TRUE, y = TRUE, method = "breslow", data = TR)
#'
#' lpnew <- predict(train.fit, newdata = TE)
#' Surv.rsp <- Surv(TR$futime, TR$fustat)
#' Surv.rsp.new <- Surv(TE$futime, TE$fustat)
#' times <- seq(10, 1000, 10)
#'
#' AUC_Uno <- AUC.uno(Surv.rsp, Surv.rsp.new, lpnew, times)
#' names(AUC_Uno)
#' AUC_Uno$iauc
#'
#' @noRd
AUC.uno <- function(Surv.rsp, Surv.rsp.new, lpnew, times, savesensspec = FALSE) {
  thresh <- my.sort(unique(lpnew))
  n_th <- length(thresh)
  n_t <- length(times)

  # Sensitivity, Specificity, and AUC
  auc.uno <- .C(
    "auc_uno",
    as.numeric(vector("numeric", length = n_t)),
    as.numeric(0),
    as.numeric(vector("numeric", length = n_t * (n_th + 1)) + 1),
    as.numeric(vector("numeric", length = n_t * (n_th + 1))),
    as.numeric(Surv.rsp[, 1]),
    as.numeric(1 - Surv.rsp[, 2]),
    as.numeric(thresh),
    as.numeric(times),
    as.numeric(lpnew),
    as.numeric(Surv.rsp.new[, 1]),
    as.numeric(Surv.rsp.new[, 2]),
    as.integer(n_th),
    as.integer(n_t),
    as.integer(dim(Surv.rsp.new)[1]),
    as.integer(dim(Surv.rsp)[1]),
    PACKAGE = "hdnom"
  )
  if (!savesensspec) {
    erg <- list(auc = auc.uno[[1]], times = auc.uno[[8]], iauc = auc.uno[[2]])
  } else {
    erg <- list(
      auc = auc.uno[[1]], times = auc.uno[[8]], iauc = auc.uno[[2]],
      sens = matrix(auc.uno[[3]], n_t, n_th + 1),
      spec = matrix(auc.uno[[4]], n_t, n_th + 1),
      thresh = auc.uno[[7]]
    )
  }
  class(erg) <- "survAUC"
  erg
}
