#' Compute validation measures for glmnet objects
#'
#' @importFrom glmnet glmnet
#' @importFrom survival Surv
#'
#' @return time-dependent AUC (tAUC) value
#'
#' @keywords internal
glmnet_validate_tauc <- function(
  x_tr, x_te, y_tr, y_te,
  alpha, lambda, pen.factor,
  tauc.type, tauc.time) {
  if (is.null(pen.factor)) {
    samp_fit <- glmnet(
      x = x_tr, y = y_tr, family = "cox",
      alpha = alpha, lambda = lambda
    )
  } else {
    samp_fit <- glmnet(
      x = x_tr, y = y_tr, family = "cox",
      alpha = alpha, lambda = lambda,
      penalty.factor = pen.factor
    )
  }

  lp_tr <- as.vector(predict(samp_fit, newx = x_tr, type = "link"))
  lp_te <- as.vector(predict(samp_fit, newx = x_te, type = "link"))

  tauc_list <- switch(

    tauc.type,

    CD = {
      AUC.cd(
        Surv.rsp = y_tr, Surv.rsp.new = y_te,
        lp = lp_tr, lpnew = lp_te,
        times = tauc.time
      )
    },

    SZ = {
      AUC.sh(
        Surv.rsp = y_tr, Surv.rsp.new = y_te,
        lp = lp_tr, lpnew = lp_te,
        times = tauc.time
      )
    },

    UNO = {
      AUC.uno(
        Surv.rsp = y_tr, Surv.rsp.new = y_te,
        lpnew = lp_te,
        times = tauc.time
      )
    }
  )

  tauc_list
}

#' Compute validation measures for ncvreg model objects
#'
#' @importFrom ncvreg ncvsurv
#' @importFrom survival Surv
#'
#' @return time-dependent AUC (tAUC) value
#'
#' @keywords internal
ncvreg_validate_tauc <- function(
  x_tr, x_te, y_tr, y_te, model.type,
  gamma, alpha, lambda,
  tauc.type, tauc.time) {
  if (model.type == "mcp") {
    samp_fit <- ncvreg::ncvsurv(
      X = x_tr, y = y_tr,
      penalty = "MCP", gamma = gamma,
      alpha = 1, lambda = lambda
    )
  }

  if (model.type == "mnet") {
    samp_fit <- ncvreg::ncvsurv(
      X = x_tr, y = y_tr,
      penalty = "MCP", gamma = gamma,
      alpha = alpha, lambda = lambda
    )
  }

  if (model.type == "scad") {
    samp_fit <- ncvreg::ncvsurv(
      X = x_tr, y = y_tr,
      penalty = "SCAD", gamma = gamma,
      alpha = 1, lambda = lambda
    )
  }

  if (model.type == "snet") {
    samp_fit <- ncvreg::ncvsurv(
      X = x_tr, y = y_tr,
      penalty = "SCAD", gamma = gamma,
      alpha = alpha, lambda = lambda
    )
  }

  lp_tr <- as.vector(predict(samp_fit, X = x_tr, type = "link"))
  lp_te <- as.vector(predict(samp_fit, X = x_te, type = "link"))

  tauc_list <- switch(

    tauc.type,

    CD = {
      AUC.cd(
        Surv.rsp = y_tr, Surv.rsp.new = y_te,
        lp = lp_tr, lpnew = lp_te,
        times = tauc.time
      )
    },

    SZ = {
      AUC.sh(
        Surv.rsp = y_tr, Surv.rsp.new = y_te,
        lp = lp_tr, lpnew = lp_te,
        times = tauc.time
      )
    },

    UNO = {
      AUC.uno(
        Surv.rsp = y_tr, Surv.rsp.new = y_te,
        lpnew = lp_te,
        times = tauc.time
      )
    }
  )

  tauc_list
}

#' Compute validation measures for penfit model objects
#'
#' @importFrom penalized penalized
#' @importFrom survival Surv
#'
#' @return time-dependent AUC (tAUC) value
#'
#' @keywords internal
penalized_validate_tauc <- function(
  x_tr, x_te, y_tr, y_te,
  lambda1, lambda2,
  tauc.type, tauc.time) {
  samp_fit <- penalized(
    response = y_tr, penalized = x_tr,
    lambda1 = lambda1, lambda2 = lambda2,
    maxiter = 25, epsilon = 1e-3, # for faster convergence, consistent with `fit_flasso()`
    fusedl = TRUE, standardize = TRUE, model = "cox"
  )

  lp_tr <- as.vector(samp_fit@lin.pred)
  lp_te <- as.vector(x_te %*% as.matrix(samp_fit@penalized))

  tauc_list <- switch(

    tauc.type,

    CD = {
      AUC.cd(
        Surv.rsp = y_tr, Surv.rsp.new = y_te,
        lp = lp_tr, lpnew = lp_te,
        times = tauc.time
      )
    },

    SZ = {
      AUC.sh(
        Surv.rsp = y_tr, Surv.rsp.new = y_te,
        lp = lp_tr, lpnew = lp_te,
        times = tauc.time
      )
    },

    UNO = {
      AUC.uno(
        Surv.rsp = y_tr, Surv.rsp.new = y_te,
        lpnew = lp_te,
        times = tauc.time
      )
    }
  )

  tauc_list
}
