#' Model selection for high-dimensional Cox models with lasso penalty
#'
#' Automatic model selection for high-dimensional Cox models
#' with lasso penalty, evaluated by penalized partial-likelihood.
#'
#' @param x Data matrix.
#' @param y Response matrix made by \code{\link[survival]{Surv}}.
#' @param nfolds Fold numbers of cross-validation.
#' @param rule Model selection criterion, \code{"lambda.min"} or
#' \code{"lambda.1se"}. See \code{\link[glmnet]{cv.glmnet}}
#' for details.
#' @param seed A random seed for cross-validation fold division.
#'
#' @export fit_lasso
#'
#' @examples
#' library("survival")
#'
#' # Load imputed SMART data
#' data("smart")
#' x <- as.matrix(smart[, -c(1, 2)])
#' time <- smart$TEVENT
#' event <- smart$EVENT
#' y <- Surv(time, event)
#'
#' # Fit Cox model with lasso penalty
#' fit <- fit_lasso(x, y, nfolds = 5, rule = "lambda.1se", seed = 11)
#'
#' nom <- as_nomogram(
#'   fit$lasso_model,
#'   model.type = "lasso",
#'   x, time, event, pred.at = 365 * 2,
#'   funlabel = "2-Year Overall Survival Probability"
#' )
#'
#' plot(nom)
fit_lasso <- function(
  x, y, nfolds = 5L,
  rule = c("lambda.min", "lambda.1se"),
  seed = 1001) {
  set.seed(seed)
  lasso_cv <- cv.glmnet(x, y, family = "cox", nfolds = nfolds, alpha = 1)

  # fit the model on all the data use the parameters got by CV
  if (rule == "lambda.min") {
    best_lambda_lasso <- lasso_cv$lambda.min
  } else if (rule == "lambda.1se") {
    best_lambda_lasso <- lasso_cv$lambda.1se
  }

  lasso_full <- glmnet(
    x, y,
    family = "cox",
    lambda = best_lambda_lasso, alpha = 1
  )

  if (lasso_full$df < 0.5) {
    stop("Null model produced by the full fit (all coefficients are zero). Please try to tune rule, seed, nfolds, or increase sample size.")
  }

  coxlasso_model <- list(
    "seed" = seed,
    "lasso_best_lambda" = best_lambda_lasso,
    "lasso_model" = lasso_full
  )

  class(coxlasso_model) <- c("hdnom.model", "hdnom.model.lasso")

  coxlasso_model
}

#' Model selection for high-dimensional Cox models with adaptive lasso penalty
#'
#' Automatic model selection for high-dimensional Cox models
#' with adaptive lasso penalty, evaluated by penalized partial-likelihood.
#'
#' @param x Data matrix.
#' @param y Response matrix made by \code{\link[survival]{Surv}}.
#' @param nfolds Fold numbers of cross-validation.
#' @param rule Model selection criterion, \code{"lambda.min"} or
#' \code{"lambda.1se"}. See \code{\link[glmnet]{cv.glmnet}}
#' for details.
#' @param seed Two random seeds for cross-validation fold division
#' in two estimation steps.
#'
#' @export fit_alasso
#'
#' @examples
#' library("survival")
#'
#' # Load imputed SMART data
#' data("smart")
#' x <- as.matrix(smart[, -c(1, 2)])
#' time <- smart$TEVENT
#' event <- smart$EVENT
#' y <- Surv(time, event)
#'
#' # Fit Cox model with adaptive lasso penalty
#' fit <- fit_alasso(x, y, nfolds = 3, rule = "lambda.1se", seed = c(7, 11))
#'
#' nom <- as_nomogram(
#'   fit$alasso_model,
#'   model.type = "alasso",
#'   x, time, event, pred.at = 365 * 2,
#'   funlabel = "2-Year Overall Survival Probability"
#' )
#'
#' plot(nom)
fit_alasso <- function(
  x, y, nfolds = 5L,
  rule = c("lambda.min", "lambda.1se"),
  seed = c(1001, 1002)) {

  # Tuning lambda for the both two stages of adaptive lasso estimation
  set.seed(seed[1L])
  lasso_cv <- cv.glmnet(x, y, family = "cox", nfolds = nfolds, alpha = 0)

  # fit the model on all the data use the parameters got by CV
  if (rule == "lambda.min") {
    best_lambda_lasso <- lasso_cv$lambda.min
  } else if (rule == "lambda.1se") {
    best_lambda_lasso <- lasso_cv$lambda.1se
  }

  lasso_full <- glmnet(
    x, y,
    family = "cox",
    lambda = best_lambda_lasso, alpha = 0
  )

  bhat <- as.matrix(lasso_full$beta)
  if (all(bhat == 0)) bhat <- rep(.Machine$double.eps * 2, length(bhat))

  # adaptive penalty
  adpen <- (1 / pmax(abs(bhat), .Machine$double.eps))

  set.seed(seed[2L])
  alasso_cv <- cv.glmnet(
    x, y,
    family = "cox", nfolds = nfolds, alpha = 1,
    penalty.factor = adpen
  )

  # fit the model on all the data use the parameters got by CV
  if (rule == "lambda.min") {
    best_lambda_alasso <- alasso_cv$lambda.min
  } else if (rule == "lambda.1se") {
    best_lambda_alasso <- alasso_cv$lambda.1se
  }

  alasso_full <- glmnet(
    x, y,
    family = "cox", lambda = best_lambda_alasso,
    alpha = 1, penalty.factor = adpen
  )

  if (alasso_full$df < 0.5) {
    stop("Null model produced by the full fit (all coefficients are zero). Please try to tune rule, seed, nfolds, or increase sample size.")
  }

  adpen_vec <- as.vector(adpen)
  adpen_name <- rownames(adpen)
  names(adpen_vec) <- adpen_name

  coxalasso_model <- list(
    "seed" = seed,
    "ridge_best_lambda" = best_lambda_lasso,
    "ridge_model" = lasso_full,
    "alasso_best_lambda" = best_lambda_alasso,
    "alasso_model" = alasso_full,
    "pen_factor" = adpen_vec
  )

  class(coxalasso_model) <- c("hdnom.model", "hdnom.model.alasso")

  coxalasso_model
}

#' Model selection for high-dimensional Cox models with elastic-net penalty
#'
#' Automatic model selection for high-dimensional Cox models
#' with elastic-net penalty, evaluated by penalized partial-likelihood.
#'
#' @param x Data matrix.
#' @param y Response matrix made by \code{\link[survival]{Surv}}.
#' @param nfolds Fold numbers of cross-validation.
#' @param alphas Alphas to tune in \code{\link[glmnet]{cv.glmnet}}.
#' @param rule Model selection criterion, \code{"lambda.min"} or
#' \code{"lambda.1se"}. See \code{\link[glmnet]{cv.glmnet}}
#' for details.
#' @param seed A random seed for cross-validation fold division.
#' @param parallel Logical. Enable parallel parameter tuning or not,
#' default is {FALSE}. To enable parallel tuning, load the
#' \code{doParallel} package and run \code{registerDoParallel()}
#' with the number of CPU cores before calling this function.
#'
#' @export fit_enet
#'
#' @examples
#' library("survival")
#'
#' # Load imputed SMART data
#' data("smart")
#' x <- as.matrix(smart[, -c(1, 2)])
#' time <- smart$TEVENT
#' event <- smart$EVENT
#' y <- Surv(time, event)
#'
#' # To enable parallel parameter tuning, first run:
#' # library("doParallel")
#' # registerDoParallel(detectCores())
#' # then set fit_enet(..., parallel = TRUE).
#'
#' # Fit Cox model with elastic-net penalty
#' fit <- fit_enet(x, y,
#'   nfolds = 3, alphas = c(0.3, 0.7),
#'   rule = "lambda.1se", seed = 11
#' )
#'
#' nom <- as_nomogram(
#'   fit$enet_model,
#'   model.type = "enet",
#'   x, time, event, pred.at = 365 * 2,
#'   funlabel = "2-Year Overall Survival Probability"
#' )
#'
#' plot(nom)
fit_enet <- function(
  x, y, nfolds = 5L, alphas = seq(0.05, 0.95, 0.05),
  rule = c("lambda.min", "lambda.1se"),
  seed = 1001, parallel = FALSE) {
  enet_cv <- glmnet_tune_alpha(
    x, y,
    family = "cox",
    nfolds = nfolds, alphas = alphas,
    seed = seed, parallel = parallel
  )

  # fit the model on all the data use the parameters got by CV
  best_alpha_enet <- enet_cv$best.alpha

  if (rule == "lambda.min") {
    best_lambda_enet <- enet_cv$best.model$lambda.min
  } else if (rule == "lambda.1se") {
    best_lambda_enet <- enet_cv$best.model$lambda.1se
  }

  enet_full <- glmnet(
    x, y,
    family = "cox",
    lambda = best_lambda_enet,
    alpha = best_alpha_enet
  )

  if (enet_full$df < 0.5) {
    stop("Null model produced by the full fit (all coefficients are zero). Please try to tune rule, alphas, seed, nfolds, or increase sample size.")
  }

  coxenet_model <- list(
    "seed" = seed,
    "enet_best_alpha" = best_alpha_enet,
    "enet_best_lambda" = best_lambda_enet,
    "enet_model" = enet_full
  )

  class(coxenet_model) <- c("hdnom.model", "hdnom.model.enet")

  coxenet_model
}

#' Model selection for high-dimensional Cox models with adaptive elastic-net penalty
#'
#' Automatic model selection for high-dimensional Cox models
#' with adaptive elastic-net penalty, evaluated by penalized partial-likelihood.
#'
#' @param x Data matrix.
#' @param y Response matrix made with \code{\link[survival]{Surv}}.
#' @param nfolds Fold numbers of cross-validation.
#' @param alphas Alphas to tune in \code{\link[glmnet]{cv.glmnet}}.
#' @param rule Model selection criterion, \code{"lambda.min"} or
#' \code{"lambda.1se"}. See \code{\link[glmnet]{cv.glmnet}}
#' for details.
#' @param seed Two random seeds for cross-validation fold division
#' in two estimation steps.
#' @param parallel Logical. Enable parallel parameter tuning or not,
#' default is {FALSE}. To enable parallel tuning, load the
#' \code{doParallel} package and run \code{registerDoParallel()}
#' with the number of CPU cores before calling this function.
#'
#' @importFrom glmnet glmnet
#'
#' @export fit_aenet
#'
#' @examples
#' library("survival")
#'
#' # Load imputed SMART data
#' data("smart")
#' x <- as.matrix(smart[, -c(1, 2)])
#' time <- smart$TEVENT
#' event <- smart$EVENT
#' y <- Surv(time, event)
#'
#' # To enable parallel parameter tuning, first run:
#' # library("doParallel")
#' # registerDoParallel(detectCores())
#' # then set fit_aenet(..., parallel = TRUE).
#'
#' # Fit Cox model with adaptive elastic-net penalty
#' fit <- fit_aenet(
#'   x, y,
#'   nfolds = 3, alphas = c(0.3, 0.7),
#'   rule = "lambda.1se", seed = c(5, 7)
#' )
#'
#' nom <- as_nomogram(
#'   fit$aenet_model,
#'   model.type = "aenet",
#'   x, time, event, pred.at = 365 * 2,
#'   funlabel = "2-Year Overall Survival Probability"
#' )
#'
#' plot(nom)
fit_aenet <- function(
  x, y, nfolds = 5L, alphas = seq(0.05, 0.95, 0.05),
  rule = c("lambda.min", "lambda.1se"),
  seed = c(1001, 1002),
  parallel = FALSE) {
  rule <- match.arg(rule)

  # Tuning alpha for the both two stages of adaptive enet estimation
  enet_cv <- glmnet_tune_alpha(
    x, y,
    family = "cox",
    nfolds = nfolds, alphas = alphas,
    seed = seed[1L], parallel = parallel
  )

  # fit the model on all the data use the parameters got by CV
  best_alpha_enet <- enet_cv$best.alpha

  if (rule == "lambda.min") {
    best_lambda_enet <- enet_cv$best.model$lambda.min
  } else if (rule == "lambda.1se") {
    best_lambda_enet <- enet_cv$best.model$lambda.1se
  }

  enet_full <- glmnet(
    x, y,
    family = "cox",
    lambda = best_lambda_enet,
    alpha = best_alpha_enet
  )

  bhat <- as.matrix(enet_full$beta)
  if (all(bhat == 0)) bhat <- rep(.Machine$double.eps * 2, length(bhat))

  # adaptive penalty
  adpen <- (1 / pmax(abs(bhat), .Machine$double.eps))

  aenet_cv <- glmnet_tune_alpha(
    x, y,
    family = "cox", nfolds = nfolds,
    exclude = which(bhat == 0),
    penalty.factor = adpen,
    alphas = alphas,
    seed = seed[2L],
    parallel = parallel
  )

  # fit the model on all the data use the parameters got by CV
  best_alpha_aenet <- aenet_cv$best.alpha

  if (rule == "lambda.min") {
    best_lambda_aenet <- aenet_cv$best.model$lambda.min
  } else if (rule == "lambda.1se") {
    best_lambda_aenet <- aenet_cv$best.model$lambda.1se
  }

  aenet_full <- glmnet(
    x, y,
    family = "cox",
    exclude = which(bhat == 0),
    lambda = best_lambda_aenet,
    penalty.factor = adpen,
    alpha = best_alpha_aenet
  )

  if (aenet_full$df < 0.5) {
    stop("Null model produced by the full fit (all coefficients are zero). Please try to tune rule, alphas, seed, nfolds, or increase sample size.")
  }

  adpen_vec <- as.vector(adpen)
  adpen_name <- rownames(adpen)
  names(adpen_vec) <- adpen_name

  coxaenet_model <- list(
    "seed" = seed,
    "enet_best_alpha" = best_alpha_enet,
    "enet_best_lambda" = best_lambda_enet,
    "enet_model" = enet_full,
    "aenet_best_alpha" = best_alpha_aenet,
    "aenet_best_lambda" = best_lambda_aenet,
    "aenet_model" = aenet_full,
    "pen_factor" = adpen_vec
  )

  class(coxaenet_model) <- c("hdnom.model", "hdnom.model.aenet")

  coxaenet_model
}

#' Model selection for high-dimensional Cox models with SCAD penalty
#'
#' Automatic model selection for high-dimensional Cox models
#' with SCAD penalty, evaluated by penalized partial-likelihood.
#'
#' @param x Data matrix.
#' @param y Response matrix made by \code{\link[survival]{Surv}}.
#' @param nfolds Fold numbers of cross-validation.
#' @param gammas Gammas to tune in \code{\link[ncvreg]{cv.ncvsurv}}.
#' @param eps Convergence threshhold.
#' @param max.iter Maximum number of iterations.
#' @param seed A random seed for cross-validation fold division.
#' @param trace Output the cross-validation parameter tuning
#' progress or not. Default is \code{FALSE}.
#' @param parallel Logical. Enable parallel parameter tuning or not,
#' default is {FALSE}. To enable parallel tuning, load the
#' \code{doParallel} package and run \code{registerDoParallel()}
#' with the number of CPU cores before calling this function.
#'
#' @importFrom ncvreg ncvsurv
#'
#' @export fit_scad
#'
#' @examples
#' library("survival")
#'
#' # Load imputed SMART data; only use the first 120 samples
#' data("smart")
#' x <- as.matrix(smart[, -c(1, 2)])[1:120, ]
#' time <- smart$TEVENT[1:120]
#' event <- smart$EVENT[1:120]
#' y <- Surv(time, event)
#'
#' # Fit Cox model with SCAD penalty
#' fit <- fit_scad(
#'   x, y,
#'   nfolds = 3, gammas = c(3.7, 5),
#'   max.iter = 15000, seed = 1010
#' )
#'
#' nom <- as_nomogram(
#'   fit$scad_model,
#'   model.type = "scad",
#'   x, time, event, pred.at = 365 * 2,
#'   funlabel = "2-Year Overall Survival Probability"
#' )
#'
#' plot(nom)
fit_scad <- function(
  x, y, nfolds = 5L,
  gammas = c(2.01, 2.3, 3.7, 200),
  eps = 1e-4, max.iter = 10000L,
  seed = 1001, trace = FALSE, parallel = FALSE) {
  scad_cv <- ncvreg_tune_gamma(
    x, y,
    penalty = "SCAD", alpha = 1,
    nfolds = nfolds, gammas = gammas,
    eps = eps, max.iter = max.iter,
    seed = seed, trace = trace, parallel = parallel
  )

  scad_best_gamma <- scad_cv$best.gamma
  scad_best_lambda <- scad_cv$best.model$lambda.min

  # fit the model on all the data use the parameters got by CV
  scad_full <- ncvreg::ncvsurv(
    x, y,
    penalty = "SCAD", alpha = 1,
    gamma = scad_best_gamma, lambda = scad_best_lambda,
    eps = eps, max.iter = max.iter
  )

  if (all(abs(scad_full$beta[-1L, ]) < .Machine$double.eps)) {
    stop("Null model produced by the full fit (all coefficients are zero). Please try to tune gammas, seed, nfolds, or increase sample size.")
  }

  coxscad_model <- list(
    "seed" = seed,
    "scad_best_gamma" = scad_best_gamma,
    "scad_best_lambda" = scad_best_lambda,
    "scad_model" = scad_full
  )

  class(coxscad_model) <- c("hdnom.model", "hdnom.model.scad")

  coxscad_model
}

#' Model selection for high-dimensional Cox models with Snet penalty
#'
#' Automatic model selection for high-dimensional Cox models
#' with Snet penalty, evaluated by penalized partial-likelihood.
#'
#' @param x Data matrix.
#' @param y Response matrix made by \code{\link[survival]{Surv}}.
#' @param nfolds Fold numbers of cross-validation.
#' @param gammas Gammas to tune in \code{\link[ncvreg]{cv.ncvsurv}}.
#' @param alphas Alphas to tune in \code{\link[ncvreg]{cv.ncvsurv}}.
#' @param eps Convergence threshhold.
#' @param max.iter Maximum number of iterations.
#' @param seed A random seed for cross-validation fold division.
#' @param trace Output the cross-validation parameter tuning
#' progress or not. Default is \code{FALSE}.
#' @param parallel Logical. Enable parallel parameter tuning or not,
#' default is {FALSE}. To enable parallel tuning, load the
#' \code{doParallel} package and run \code{registerDoParallel()}
#' with the number of CPU cores before calling this function.
#'
#' @importFrom ncvreg ncvsurv
#'
#' @export fit_snet
#'
#' @examples
#' library("survival")
#'
#' # Load imputed SMART data; only use the first 120 samples
#' data("smart")
#' x <- as.matrix(smart[, -c(1, 2)])[1:120, ]
#' time <- smart$TEVENT[1:120]
#' event <- smart$EVENT[1:120]
#' y <- Surv(time, event)
#'
#' # Fit Cox model with Snet penalty
#' fit <- fit_snet(
#'   x, y,
#'   nfolds = 3,
#'   gammas = 3.7, alphas = c(0.3, 0.8),
#'   max.iter = 15000, seed = 1010
#' )
#'
#' nom <- as_nomogram(
#'   fit$snet_model,
#'   model.type = "snet",
#'   x, time, event, pred.at = 365 * 2,
#'   funlabel = "2-Year Overall Survival Probability"
#' )
#'
#' plot(nom)
fit_snet <- function(
  x, y, nfolds = 5L,
  gammas = c(2.01, 2.3, 3.7, 200),
  alphas = seq(0.05, 0.95, 0.05),
  eps = 1e-4, max.iter = 10000L,
  seed = 1001, trace = FALSE, parallel = FALSE) {
  snet_cv <- ncvreg_tune_gamma_alpha(
    x, y,
    penalty = "SCAD",
    nfolds = nfolds,
    gammas = gammas, alphas = alphas,
    eps = eps, max.iter = max.iter,
    seed = seed, trace = trace,
    parallel = parallel
  )

  snet_best_gamma <- snet_cv$best.gamma
  snet_best_alpha <- snet_cv$best.alpha
  snet_best_lambda <- snet_cv$best.model$lambda.min

  # fit the model on all the data use the parameters got by CV
  snet_full <- ncvreg::ncvsurv(
    x, y,
    penalty = "SCAD",
    gamma = snet_best_gamma,
    alpha = snet_best_alpha,
    lambda = snet_best_lambda,
    eps = eps, max.iter = max.iter
  )

  if (all(abs(snet_full$beta[-1L, ]) < .Machine$double.eps)) {
    stop("Null model produced by the full fit (all coefficients are zero). Please try to tune gammas, alphas, seed, nfolds, or increase sample size.")
  }

  coxsnet_model <- list(
    "seed" = seed,
    "snet_best_gamma" = snet_best_gamma,
    "snet_best_alpha" = snet_best_alpha,
    "snet_best_lambda" = snet_best_lambda,
    "snet_model" = snet_full
  )

  class(coxsnet_model) <- c("hdnom.model", "hdnom.model.snet")

  coxsnet_model
}

#' Model selection for high-dimensional Cox models with MCP penalty
#'
#' Automatic model selection for high-dimensional Cox models
#' with MCP penalty, evaluated by penalized partial-likelihood.
#'
#' @param x Data matrix.
#' @param y Response matrix made by \code{\link[survival]{Surv}}.
#' @param nfolds Fold numbers of cross-validation.
#' @param gammas Gammas to tune in \code{\link[ncvreg]{cv.ncvsurv}}.
#' @param eps Convergence threshhold.
#' @param max.iter Maximum number of iterations.
#' @param seed A random seed for cross-validation fold division.
#' @param trace Output the cross-validation parameter tuning
#' progress or not. Default is \code{FALSE}.
#' @param parallel Logical. Enable parallel parameter tuning or not,
#' default is {FALSE}. To enable parallel tuning, load the
#' \code{doParallel} package and run \code{registerDoParallel()}
#' with the number of CPU cores before calling this function.
#'
#' @importFrom ncvreg ncvsurv
#'
#' @export fit_mcp
#'
#' @examples
#' library("survival")
#'
#' # Load imputed SMART data; only use the first 150 samples
#' data("smart")
#' x <- as.matrix(smart[, -c(1, 2)])[1:150, ]
#' time <- smart$TEVENT[1:150]
#' event <- smart$EVENT[1:150]
#' y <- Surv(time, event)
#'
#' # Fit Cox model with MCP penalty
#' fit <- fit_mcp(x, y, nfolds = 3, gammas = c(2.1, 3), seed = 1001)
#'
#' nom <- as_nomogram(
#'   fit$mcp_model,
#'   model.type = "mcp",
#'   x, time, event, pred.at = 365 * 2,
#'   funlabel = "2-Year Overall Survival Probability"
#' )
#'
#' plot(nom)
fit_mcp <- function(
  x, y, nfolds = 5L, gammas = c(1.01, 1.7, 3, 100),
  eps = 1e-4, max.iter = 10000L,
  seed = 1001, trace = FALSE, parallel = FALSE) {
  mcp_cv <- ncvreg_tune_gamma(
    x, y,
    penalty = "MCP", alpha = 1,
    nfolds = nfolds, gammas = gammas,
    eps = eps, max.iter = max.iter,
    seed = seed, trace = trace, parallel = parallel
  )

  mcp_best_gamma <- mcp_cv$best.gamma
  mcp_best_lambda <- mcp_cv$best.model$lambda.min

  # fit the model on all the data use the parameters got by CV
  mcp_full <-
    ncvreg::ncvsurv(
      x, y,
      penalty = "MCP", alpha = 1,
      gamma = mcp_best_gamma, lambda = mcp_best_lambda,
      eps = eps, max.iter = max.iter
    )

  # deal with null models, thanks for the suggestion from Patrick Breheny
  if (all(abs(mcp_full$beta[-1L, ]) < .Machine$double.eps)) {
    stop("Null model produced by the full fit (all coefficients are zero). Please try to tune gammas, seed, nfolds, or increase sample size.")
  }

  coxmcp_model <- list(
    "seed" = seed,
    "mcp_best_gamma" = mcp_best_gamma,
    "mcp_best_lambda" = mcp_best_lambda,
    "mcp_model" = mcp_full
  )

  class(coxmcp_model) <- c("hdnom.model", "hdnom.model.mcp")

  coxmcp_model
}

#' Model selection for high-dimensional Cox models with Mnet penalty
#'
#' Automatic model selection for high-dimensional Cox models
#' with Mnet penalty, evaluated by penalized partial-likelihood.
#'
#' @param x Data matrix.
#' @param y Response matrix made by \code{\link[survival]{Surv}}.
#' @param nfolds Fold numbers of cross-validation.
#' @param gammas Gammas to tune in \code{\link[ncvreg]{cv.ncvsurv}}.
#' @param alphas Alphas to tune in \code{\link[ncvreg]{cv.ncvsurv}}.
#' @param eps Convergence threshhold.
#' @param max.iter Maximum number of iterations.
#' @param seed A random seed for cross-validation fold division.
#' @param trace Output the cross-validation parameter tuning
#' progress or not. Default is \code{FALSE}.
#' @param parallel Logical. Enable parallel parameter tuning or not,
#' default is {FALSE}. To enable parallel tuning, load the
#' \code{doParallel} package and run \code{registerDoParallel()}
#' with the number of CPU cores before calling this function.
#'
#' @importFrom ncvreg ncvsurv
#'
#' @export fit_mnet
#'
#' @examples
#' library("survival")
#'
#' # Load imputed SMART data; only use the first 120 samples
#' data("smart")
#' x <- as.matrix(smart[, -c(1, 2)])[1:120, ]
#' time <- smart$TEVENT[1:120]
#' event <- smart$EVENT[1:120]
#' y <- Surv(time, event)
#'
#' # Fit Cox model with Mnet penalty
#' fit <- fit_mnet(
#'   x, y,
#'   nfolds = 3,
#'   gammas = 3, alphas = c(0.3, 0.8),
#'   max.iter = 15000, seed = 1010
#' )
#'
#' nom <- as_nomogram(
#'   fit$mnet_model,
#'   model.type = "mnet",
#'   x, time, event, pred.at = 365 * 2,
#'   funlabel = "2-Year Overall Survival Probability"
#' )
#'
#' plot(nom)
fit_mnet <- function(
  x, y, nfolds = 5L,
  gammas = c(1.01, 1.7, 3, 100),
  alphas = seq(0.05, 0.95, 0.05),
  eps = 1e-4, max.iter = 10000L,
  seed = 1001, trace = FALSE, parallel = FALSE) {
  mnet_cv <- ncvreg_tune_gamma_alpha(
    x, y,
    penalty = "MCP",
    nfolds = nfolds,
    gammas = gammas, alphas = alphas,
    eps = eps, max.iter = max.iter,
    seed = seed, trace = trace,
    parallel = parallel
  )

  mnet_best_gamma <- mnet_cv$best.gamma
  mnet_best_alpha <- mnet_cv$best.alpha
  mnet_best_lambda <- mnet_cv$best.model$lambda.min

  # fit the model on all the data use the parameters got by CV
  mnet_full <-
    ncvreg::ncvsurv(
      x, y,
      penalty = "MCP",
      gamma = mnet_best_gamma,
      alpha = mnet_best_alpha,
      lambda = mnet_best_lambda,
      eps = eps, max.iter = max.iter
    )

  if (all(abs(mnet_full$beta[-1L, ]) < .Machine$double.eps)) {
    stop("Null model produced by the full fit (all coefficients are zero). Please try to tune gammas, alphas, seed, nfolds, or increase sample size.")
  }

  coxmnet_model <- list(
    "seed" = seed,
    "mnet_best_gamma" = mnet_best_gamma,
    "mnet_best_alpha" = mnet_best_alpha,
    "mnet_best_lambda" = mnet_best_lambda,
    "mnet_model" = mnet_full
  )

  class(coxmnet_model) <- c("hdnom.model", "hdnom.model.mnet")

  coxmnet_model
}

#' Model selection for high-dimensional Cox models with fused lasso penalty
#'
#' Automatic model selection for high-dimensional Cox models
#' with fused lasso penalty, evaluated by cross-validated likelihood.
#'
#' @param x Data matrix.
#' @param y Response matrix made by \code{\link[survival]{Surv}}.
#' @param nfolds Fold numbers of cross-validation.
#' @param lambda1 Vector of lambda1 candidates.
#' Default is \code{0.001, 0.05, 0.5, 1, 5}.
#' @param lambda2 Vector of lambda2 candidates.
#' Default is \code{0.001, 0.01, 0.5}.
#' @param maxiter The maximum number of iterations allowed.
#' Default is \code{25}.
#' @param epsilon The convergence criterion.
#' Default is \code{1e-3}.
#' @param seed A random seed for cross-validation fold division.
#' @param trace Output the cross-validation parameter tuning
#' progress or not. Default is \code{FALSE}.
#' @param parallel Logical. Enable parallel parameter tuning or not,
#' default is {FALSE}. To enable parallel tuning, load the
#' \code{doParallel} package and run \code{registerDoParallel()}
#' with the number of CPU cores before calling this function.
#' @param ... other parameters to \code{\link[penalized]{cvl}}
#' and \code{\link[penalized]{penalized}}.
#'
#' @note The cross-validation procedure used in this function is the
#' \emph{approximated cross-validation} provided by the \code{penalized}
#' package. Be careful dealing with the results since they might be more
#' optimistic than a traditional CV procedure. This cross-validation
#' method is more suitable for datasets with larger number of observations,
#' and a higher number of cross-validation folds.
#'
#' @importFrom penalized penalized
#'
#' @export fit_flasso
#'
#' @examples
#' library("survival")
#'
#' # Load imputed SMART data; only use the first 120 samples
#' data("smart")
#' x <- as.matrix(smart[, -c(1, 2)])[1:120, ]
#' time <- smart$TEVENT[1:120]
#' event <- smart$EVENT[1:120]
#' y <- Surv(time, event)
#'
#' # Fit Cox model with fused lasso penalty
#' fit <- fit_flasso(x, y,
#'   lambda1 = c(1, 10), lambda2 = c(0.01),
#'   nfolds = 3, seed = 11
#' )
#'
#' nom <- as_nomogram(
#'   fit$flasso_model,
#'   model.type = "flasso",
#'   x, time, event, pred.at = 365 * 2,
#'   funlabel = "2-Year Overall Survival Probability"
#' )
#'
#' plot(nom)
fit_flasso <- function(
  x, y, nfolds = 5L,
  lambda1 = c(0.001, 0.05, 0.5, 1, 5),
  lambda2 = c(0.001, 0.01, 0.5),
  maxiter = 25, epsilon = 1e-3,
  seed = 1001, trace = FALSE, parallel = FALSE, ...) {
  if (trace) cat("Starting cross-validation...\n")
  flasso_cv <- penalized_tune_lambda(
    response = y, penalized = x, fold = nfolds,
    lambda1 = lambda1, lambda2 = lambda2,
    maxiter = maxiter, epsilon = epsilon,
    seed = seed, trace = trace, parallel = parallel,
    fusedl = TRUE, standardize = TRUE, model = "cox", ...
  )

  flasso_best_lambda1 <- flasso_cv$"best.lambda1"
  flasso_best_lambda2 <- flasso_cv$"best.lambda2"

  # fit the model on all the data use the parameters got by CV
  if (trace) cat("Fitting fused lasso model with full data...\n")
  flasso_full <- penalized(
    response = y, penalized = x,
    lambda1 = flasso_best_lambda1,
    lambda2 = flasso_best_lambda2,
    maxiter = maxiter, epsilon = epsilon,
    trace = trace,
    fusedl = TRUE, standardize = FALSE, model = "cox", ...
  )

  if (all(abs(flasso_full@penalized) < .Machine$double.eps)) {
    stop("Null model produced by the full fit (all coefficients are zero). Please try changing the seed, nfolds, or increase sample size.")
  }

  coxflasso_model <- list(
    "seed" = seed,
    "flasso_best_lambda1" = flasso_best_lambda1,
    "flasso_best_lambda2" = flasso_best_lambda2,
    "flasso_model" = flasso_full
  )

  class(coxflasso_model) <- c("hdnom.model", "hdnom.model.flasso")

  coxflasso_model
}
