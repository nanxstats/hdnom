# Functions for model selection with CV + penalized partial likelihood
# or AIC/BIC/EBIC. We use these functions to build a model first,
# then validate the model using validate and calibrate functions.

#' Automatic alpha tuning function by k-fold cross-validation
#'
#' @return best model object and best alpha
#'
#' @importFrom glmnet cv.glmnet
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#'
#' @keywords internal
glmnet_tune_alpha <- function(..., alphas, seed, parallel) {
  if (!parallel) {
    model_list <- vector("list", length(alphas))
    for (i in 1L:length(alphas)) {
      set.seed(seed)
      model_list[[i]] <- cv.glmnet(..., alpha = alphas[i])
    }
  } else {
    model_list <- foreach(alphas = alphas) %dopar% {
      set.seed(seed)
      cv.glmnet(..., alpha = alphas)
    }
  }

  # select model for best lambda first (then alpha)
  # criterion: penalized partial likelihood
  errors <- unlist(lapply(model_list, function(x) min(sqrt(x$cvm))))

  list(
    "best.model" = model_list[[which.min(errors)]],
    "best.alpha" = alphas[which.min(errors)]
  )
}

#' Automatic MCP/SCAD gamma tuning function by k-fold cross-validation
#'
#' @return best model object and best gamma
#'
#' @importFrom ncvreg cv.ncvsurv
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#'
#' @keywords internal
ncvreg_tune_gamma <- function(..., gammas, eps, max.iter, seed, parallel) {
  if (!parallel) {
    model_list <- vector("list", length(gammas))
    for (i in 1L:length(gammas)) {
      set.seed(seed)
      model_list[[i]] <- cv.ncvsurv(
        ...,
        gamma = gammas[i],
        eps = eps, max.iter = max.iter
      )
    }
  } else {
    model_list <- foreach(gammas = gammas) %dopar% {
      set.seed(seed)
      cv.ncvsurv(
        ...,
        gamma = gammas,
        eps = eps, max.iter = max.iter
      )
    }
  }

  # select model for best lambda first (then gamma)
  # criterion: penalized partial likelihood
  errors <- unlist(lapply(model_list, function(x) min(sqrt(x$cve))))

  list(
    "best.model" = model_list[[which.min(errors)]],
    "best.gamma" = gammas[which.min(errors)]
  )
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
ncvreg_tune_gamma_alpha <- function(
  ..., gammas, alphas, eps, max.iter, seed, parallel) {
  if (!parallel) {
    model_list <- vector("list", length(gammas))
    for (k in 1L:length(model_list)) {
      model_list[[k]] <- vector("list", length(alphas))
    }

    for (i in 1L:length(gammas)) {
      for (j in 1L:length(alphas)) {
        set.seed(seed)
        model_list[[i]][[j]] <-
          cv.ncvsurv(
            ...,
            gamma = gammas[i], alpha = alphas[j],
            eps = eps, max.iter = max.iter
          )
      }
    }

    simple_model_list <- unlist(model_list, recursive = FALSE)
  } else {
    model_list <- foreach(gammas = gammas) %:%
      foreach(alphas = alphas) %dopar% {
        set.seed(seed)
        cv.ncvsurv(
          ...,
          gamma = gammas, alpha = alphas,
          eps = eps, max.iter = max.iter
        )
      }

    simple_model_list <- unlist(model_list, recursive = FALSE)
  }

  # select model for best lambda first (then gamma/alpha)
  # criterion: penalized partial likelihood
  errors <- unlist(lapply(
    simple_model_list,
    function(x) min(sqrt(x$cve))
  ))

  list(
    "best.model" = simple_model_list[[which.min(errors)]],
    "best.gamma" = simple_model_list[[which.min(errors)]]$fit$gamma,
    "best.alpha" = simple_model_list[[which.min(errors)]]$fit$alpha
  )
}

#' Automatic lambda tuning function for fused lasso by k-fold cross-validation
#'
#' @return best model object, best lambda1, and best lambda2
#'
#' @importFrom penalized cvl
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#'
#' @keywords internal
penalized_tune_lambda <- function(..., lambda1, lambda2, seed, trace, parallel) {
  nlambda1 <- length(lambda1)
  nlambda2 <- length(lambda2)

  if (!parallel) {
    model_list <- vector("list", nlambda1)
    for (k in 1L:nlambda1) model_list[[k]] <- vector("list", nlambda2)

    for (i in 1L:nlambda1) {
      for (j in 1L:nlambda2) {
        set.seed(seed)
        if (trace) cat("Starting: lambda1 =", lambda1[i], "lambda2 =", lambda2[j], "\n")
        model_list[[i]][[j]] <-
          penalized::cvl(..., lambda1 = lambda1[i], lambda2 = lambda2[j])
      }
    }

    # store lambda combinations
    for (i in 1L:nlambda1) {
      for (j in 1L:nlambda2) {
        model_list[[i]][[j]][["lambda"]] <-
          c("lambda1" = lambda1[i], "lambda2" = lambda2[j])
      }
    }

    simple_model_list <- unlist(model_list, recursive = FALSE)
  } else {
    model_list <- foreach(lambda1 = lambda1) %:%
      foreach(lambda2 = lambda2) %dopar% {
        set.seed(seed)
        penalized::cvl(..., lambda1 = lambda1, lambda2 = lambda2)
      }

    # store lambda combinations
    for (i in 1L:nlambda1) {
      for (j in 1L:nlambda2) {
        model_list[[i]][[j]][["lambda"]] <-
          c("lambda1" = lambda1[i], "lambda2" = lambda2[j])
      }
    }

    simple_model_list <- unlist(model_list, recursive = FALSE)
  }

  # choose model for best lambda combination
  # criterion: cross-validated likelihood
  max_cvl <- which.max(unlist(sapply(simple_model_list, "[", "cvl")))

  list(
    "best.model" = simple_model_list[[max_cvl]],
    "best.lambda1" = simple_model_list[[max_cvl]][["lambda"]]["lambda1"],
    "best.lambda2" = simple_model_list[[max_cvl]][["lambda"]]["lambda2"]
  )
}
