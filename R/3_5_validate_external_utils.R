#' Compute external validation measures for glmnet objects
#'
#' @importFrom glmnet glmnet
#' @importFrom survival Surv
#'
#' @return time-dependent AUC (tAUC) value
#'
#' @keywords internal
glmnet_validate_external_tauc <- function(
    object, x_tr, x_te, y_tr, y_te,
    tauc.type, tauc.time) {
  lp_tr <- as.vector(predict(object, newx = x_tr, type = "link"))
  lp_te <- as.vector(predict(object, newx = x_te, type = "link"))

  tauc_list <- switch(tauc.type,
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

#' Compute external validation measures for ncvreg model objects
#'
#' @importFrom survival Surv
#'
#' @return time-dependent AUC (tAUC) value
#'
#' @keywords internal
ncvreg_validate_external_tauc <- function(
    object, x_tr, x_te, y_tr, y_te,
    tauc.type, tauc.time) {
  lp_tr <- as.vector(predict(object, X = x_tr, type = "link"))
  lp_te <- as.vector(predict(object, X = x_te, type = "link"))

  tauc_list <- switch(tauc.type,
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

#' Compute external validation measures for penfit model objects
#'
#' @importFrom penalized penalized
#' @importFrom survival Surv
#'
#' @return time-dependent AUC (tAUC) value
#'
#' @keywords internal
penalized_validate_external_tauc <- function(
    object, x_tr, x_te, y_tr, y_te,
    tauc.type, tauc.time) {
  lp_tr <- as.vector(object@"lin.pred")
  lp_te <- as.vector(x_te %*% as.matrix(object@"penalized"))

  tauc_list <- switch(tauc.type,
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
