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
#' data("smart")
#' x <- as.matrix(smart[, -c(1, 2)])
#' time <- smart$TEVENT
#' event <- smart$EVENT
#' y <- survival::Surv(time, event)
#'
#' fit <- fit_lasso(x, y, nfolds = 5, rule = "lambda.1se", seed = 11)
#'
#' nom <- as_nomogram(
#'   fit, x, time, event, pred.at = 365 * 2,
#'   funlabel = "2-Year Overall Survival Probability"
#' )
#'
#' plot(nom)
fit_lasso <- function(
  x, y, nfolds = 5L,
  rule = c("lambda.min", "lambda.1se"),
  seed = 1001) {
  call <- match.call()
  rule <- match.arg(rule)

  set.seed(seed)
  lasso_cv <- cv.glmnet(x, y, family = "cox", nfolds = nfolds, alpha = 1)

  if (rule == "lambda.min") {
    lambda_opt <- lasso_cv$lambda.min
  } else if (rule == "lambda.1se") {
    lambda_opt <- lasso_cv$lambda.1se
  }

  lasso_full <- glmnet(
    x, y,
    family = "cox",
    lambda = lambda_opt, alpha = 1
  )

  if (lasso_full$df < 0.5) {
    stop("Null model produced by the full fit (all coefficients are zero). Please try to tune rule, seed, nfolds, or increase sample size.")
  }

  model <- list(
    "model" = lasso_full,
    "lambda" = lambda_opt,
    "type" = "lasso",
    "seed" = seed,
    "call" = call
  )

  class(model) <- c("hdnom.model", "hdnom.model.lasso")
  model
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
#' data("smart")
#' x <- as.matrix(smart[, -c(1, 2)])
#' time <- smart$TEVENT
#' event <- smart$EVENT
#' y <- survival::Surv(time, event)
#'
#' fit <- fit_alasso(x, y, nfolds = 3, rule = "lambda.1se", seed = c(7, 11))
#'
#' nom <- as_nomogram(
#'   fit, x, time, event, pred.at = 365 * 2,
#'   funlabel = "2-Year Overall Survival Probability"
#' )
#'
#' plot(nom)
fit_alasso <- function(
  x, y, nfolds = 5L,
  rule = c("lambda.min", "lambda.1se"),
  seed = c(1001, 1002)) {
  call <- match.call()
  rule <- match.arg(rule)

  # Tune lambda for the both two stages of adaptive lasso estimation
  set.seed(seed[1L])
  lasso_cv <- cv.glmnet(x, y, family = "cox", nfolds = nfolds, alpha = 0)

  if (rule == "lambda.min") {
    lambda_opt_init <- lasso_cv$lambda.min
  } else if (rule == "lambda.1se") {
    lambda_opt_init <- lasso_cv$lambda.1se
  }

  lasso_full <- glmnet(
    x, y,
    family = "cox",
    lambda = lambda_opt_init, alpha = 0
  )

  bhat <- as.matrix(lasso_full$beta)
  if (all(bhat == 0)) bhat <- rep(.Machine$double.eps * 2, length(bhat))

  # Compute adaptive penalty
  adpen <- (1 / pmax(abs(bhat), .Machine$double.eps))

  set.seed(seed[2L])
  alasso_cv <- cv.glmnet(
    x, y,
    family = "cox", nfolds = nfolds, alpha = 1,
    penalty.factor = adpen
  )

  if (rule == "lambda.min") {
    lambda_opt <- alasso_cv$lambda.min
  } else if (rule == "lambda.1se") {
    lambda_opt <- alasso_cv$lambda.1se
  }

  alasso_full <- glmnet(
    x, y,
    family = "cox", lambda = lambda_opt,
    alpha = 1, penalty.factor = adpen
  )

  if (alasso_full$df < 0.5) {
    stop("Null model produced by the full fit (all coefficients are zero). Please try to tune rule, seed, nfolds, or increase sample size.")
  }

  adpen_vec <- as.vector(adpen)
  adpen_name <- rownames(adpen)
  names(adpen_vec) <- adpen_name

  model <- list(
    "model" = alasso_full,
    "lambda" = lambda_opt,
    "model_init" = lasso_full,
    "lambda_init" = lambda_opt_init,
    "pen_factor" = adpen_vec,
    "type" = "alasso",
    "seed" = seed,
    "call" = call
  )

  class(model) <- c("hdnom.model", "hdnom.model.alasso")
  model
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
#' data("smart")
#' x <- as.matrix(smart[, -c(1, 2)])
#' time <- smart$TEVENT
#' event <- smart$EVENT
#' y <- survival::Surv(time, event)
#'
#' # To enable parallel parameter tuning, first run:
#' # library("doParallel")
#' # registerDoParallel(detectCores())
#' # then set fit_enet(..., parallel = TRUE).
#'
#' fit <- fit_enet(x, y,
#'   nfolds = 3, alphas = c(0.3, 0.7),
#'   rule = "lambda.1se", seed = 11
#' )
#'
#' nom <- as_nomogram(
#'   fit, x, time, event, pred.at = 365 * 2,
#'   funlabel = "2-Year Overall Survival Probability"
#' )
#'
#' plot(nom)
fit_enet <- function(
  x, y, nfolds = 5L, alphas = seq(0.05, 0.95, 0.05),
  rule = c("lambda.min", "lambda.1se"),
  seed = 1001, parallel = FALSE) {
  call <- match.call()
  rule <- match.arg(rule)

  enet_cv <- glmnet_tune_alpha(
    x, y,
    family = "cox",
    nfolds = nfolds, alphas = alphas,
    seed = seed, parallel = parallel
  )

  alpha_opt <- enet_cv$best.alpha

  if (rule == "lambda.min") {
    lambda_opt <- enet_cv$best.model$lambda.min
  } else if (rule == "lambda.1se") {
    lambda_opt <- enet_cv$best.model$lambda.1se
  }

  enet_full <- glmnet(
    x, y,
    family = "cox",
    lambda = lambda_opt,
    alpha = alpha_opt
  )

  if (enet_full$df < 0.5) {
    stop("Null model produced by the full fit (all coefficients are zero). Please try to tune rule, alphas, seed, nfolds, or increase sample size.")
  }

  model <- list(
    "model" = enet_full,
    "alpha" = alpha_opt,
    "lambda" = lambda_opt,
    "type" = "enet",
    "seed" = seed,
    "call" = call
  )

  class(model) <- c("hdnom.model", "hdnom.model.enet")
  model
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
#' data("smart")
#' x <- as.matrix(smart[, -c(1, 2)])
#' time <- smart$TEVENT
#' event <- smart$EVENT
#' y <- survival::Surv(time, event)
#'
#' # To enable parallel parameter tuning, first run:
#' # library("doParallel")
#' # registerDoParallel(detectCores())
#' # then set fit_aenet(..., parallel = TRUE).
#'
#' fit <- fit_aenet(
#'   x, y,
#'   nfolds = 3, alphas = c(0.3, 0.7),
#'   rule = "lambda.1se", seed = c(5, 7)
#' )
#'
#' nom <- as_nomogram(
#'   fit, x, time, event, pred.at = 365 * 2,
#'   funlabel = "2-Year Overall Survival Probability"
#' )
#'
#' plot(nom)
fit_aenet <- function(
  x, y, nfolds = 5L, alphas = seq(0.05, 0.95, 0.05),
  rule = c("lambda.min", "lambda.1se"),
  seed = c(1001, 1002),
  parallel = FALSE) {
  call <- match.call()
  rule <- match.arg(rule)

  # Tune alpha for the both two stages of adaptive enet estimation
  enet_cv <- glmnet_tune_alpha(
    x, y,
    family = "cox",
    nfolds = nfolds, alphas = alphas,
    seed = seed[1L], parallel = parallel
  )

  alpha_opt_init <- enet_cv$best.alpha

  if (rule == "lambda.min") {
    lambda_opt_init <- enet_cv$best.model$lambda.min
  } else if (rule == "lambda.1se") {
    lambda_opt_init <- enet_cv$best.model$lambda.1se
  }

  enet_full <- glmnet(
    x, y,
    family = "cox",
    lambda = lambda_opt_init,
    alpha = alpha_opt_init
  )

  bhat <- as.matrix(enet_full$beta)
  if (all(bhat == 0)) bhat <- rep(.Machine$double.eps * 2, length(bhat))

  # Compute adaptive penalty
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

  alpha_opt <- aenet_cv$best.alpha

  if (rule == "lambda.min") {
    lambda_opt <- aenet_cv$best.model$lambda.min
  } else if (rule == "lambda.1se") {
    lambda_opt <- aenet_cv$best.model$lambda.1se
  }

  aenet_full <- glmnet(
    x, y,
    family = "cox",
    exclude = which(bhat == 0),
    lambda = lambda_opt,
    alpha = alpha_opt,
    penalty.factor = adpen
  )

  if (aenet_full$df < 0.5) {
    stop("Null model produced by the full fit (all coefficients are zero). Please try to tune rule, alphas, seed, nfolds, or increase sample size.")
  }

  adpen_vec <- as.vector(adpen)
  adpen_name <- rownames(adpen)
  names(adpen_vec) <- adpen_name

  model <- list(
    "model" = aenet_full,
    "alpha" = alpha_opt,
    "lambda" = lambda_opt,
    "model_init" = enet_full,
    "alpha_init" = alpha_opt_init,
    "lambda_init" = lambda_opt_init,
    "pen_factor" = adpen_vec,
    "type" = "aenet",
    "seed" = seed,
    "call" = call
  )

  class(model) <- c("hdnom.model", "hdnom.model.aenet")
  model
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
#' data("smart")
#' x <- as.matrix(smart[, -c(1, 2)])[1:120, ]
#' time <- smart$TEVENT[1:120]
#' event <- smart$EVENT[1:120]
#' y <- survival::Surv(time, event)
#'
#' fit <- fit_scad(
#'   x, y,
#'   nfolds = 3, gammas = c(3.7, 5),
#'   max.iter = 15000, seed = 1010
#' )
#'
#' nom <- as_nomogram(
#'   fit, x, time, event, pred.at = 365 * 2,
#'   funlabel = "2-Year Overall Survival Probability"
#' )
#'
#' plot(nom)
fit_scad <- function(
  x, y, nfolds = 5L,
  gammas = c(2.01, 2.3, 3.7, 200),
  eps = 1e-4, max.iter = 10000L,
  seed = 1001, trace = FALSE, parallel = FALSE) {
  call <- match.call()

  scad_cv <- ncvreg_tune_gamma(
    x, y,
    penalty = "SCAD", alpha = 1,
    nfolds = nfolds, gammas = gammas,
    eps = eps, max.iter = max.iter,
    seed = seed, trace = trace, parallel = parallel
  )

  gamma_opt <- scad_cv$best.gamma
  lambda_opt <- scad_cv$best.model$lambda.min

  scad_full <- ncvreg::ncvsurv(
    x, y,
    penalty = "SCAD", alpha = 1,
    gamma = gamma_opt, lambda = lambda_opt,
    eps = eps, max.iter = max.iter
  )

  if (all(abs(scad_full$beta[-1L, ]) < .Machine$double.eps)) {
    stop("Null model produced by the full fit (all coefficients are zero). Please try to tune gammas, seed, nfolds, or increase sample size.")
  }

  model <- list(
    "model" = scad_full,
    "gamma" = gamma_opt,
    "lambda" = lambda_opt,
    "type" = "scad",
    "seed" = seed,
    "call" = call
  )

  class(model) <- c("hdnom.model", "hdnom.model.scad")
  model
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
#' data("smart")
#' x <- as.matrix(smart[, -c(1, 2)])[1:120, ]
#' time <- smart$TEVENT[1:120]
#' event <- smart$EVENT[1:120]
#' y <- survival::Surv(time, event)
#'
#' fit <- fit_snet(
#'   x, y,
#'   nfolds = 3,
#'   gammas = 3.7, alphas = c(0.3, 0.8),
#'   max.iter = 15000, seed = 1010
#' )
#'
#' nom <- as_nomogram(
#'   fit, x, time, event, pred.at = 365 * 2,
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
  call <- match.call()

  snet_cv <- ncvreg_tune_gamma_alpha(
    x, y,
    penalty = "SCAD",
    nfolds = nfolds,
    gammas = gammas, alphas = alphas,
    eps = eps, max.iter = max.iter,
    seed = seed, trace = trace,
    parallel = parallel
  )

  gamma_opt <- snet_cv$best.gamma
  alpha_opt <- snet_cv$best.alpha
  lambda_opt <- snet_cv$best.model$lambda.min

  snet_full <- ncvreg::ncvsurv(
    x, y,
    penalty = "SCAD",
    gamma = gamma_opt,
    alpha = alpha_opt,
    lambda = lambda_opt,
    eps = eps, max.iter = max.iter
  )

  if (all(abs(snet_full$beta[-1L, ]) < .Machine$double.eps)) {
    stop("Null model produced by the full fit (all coefficients are zero). Please try to tune gammas, alphas, seed, nfolds, or increase sample size.")
  }

  model <- list(
    "model" = snet_full,
    "gamma" = gamma_opt,
    "alpha" = alpha_opt,
    "lambda" = lambda_opt,
    "type" = "snet",
    "seed" = seed,
    "call" = call
  )

  class(model) <- c("hdnom.model", "hdnom.model.snet")
  model
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
#' data("smart")
#' x <- as.matrix(smart[, -c(1, 2)])[1:150, ]
#' time <- smart$TEVENT[1:150]
#' event <- smart$EVENT[1:150]
#' y <- survival::Surv(time, event)
#'
#' fit <- fit_mcp(x, y, nfolds = 3, gammas = c(2.1, 3), seed = 1001)
#'
#' nom <- as_nomogram(
#'   fit, x, time, event, pred.at = 365 * 2,
#'   funlabel = "2-Year Overall Survival Probability"
#' )
#'
#' plot(nom)
fit_mcp <- function(
  x, y, nfolds = 5L, gammas = c(1.01, 1.7, 3, 100),
  eps = 1e-4, max.iter = 10000L,
  seed = 1001, trace = FALSE, parallel = FALSE) {
  call <- match.call()

  mcp_cv <- ncvreg_tune_gamma(
    x, y,
    penalty = "MCP", alpha = 1,
    nfolds = nfolds, gammas = gammas,
    eps = eps, max.iter = max.iter,
    seed = seed, trace = trace, parallel = parallel
  )

  gamma_opt <- mcp_cv$best.gamma
  lambda_opt <- mcp_cv$best.model$lambda.min

  mcp_full <-
    ncvreg::ncvsurv(
      x, y,
      penalty = "MCP", alpha = 1,
      gamma = gamma_opt,
      lambda = lambda_opt,
      eps = eps, max.iter = max.iter
    )

  # Deal with null models, thanks for the suggestion from Prof. Patrick Breheny
  if (all(abs(mcp_full$beta[-1L, ]) < .Machine$double.eps)) {
    stop("Null model produced by the full fit (all coefficients are zero). Please try to tune gammas, seed, nfolds, or increase sample size.")
  }

  model <- list(
    "model" = mcp_full,
    "gamma" = gamma_opt,
    "lambda" = lambda_opt,
    "type" = "mcp",
    "seed" = seed,
    "call" = call
  )

  class(model) <- c("hdnom.model", "hdnom.model.mcp")
  model
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
#' data("smart")
#' x <- as.matrix(smart[, -c(1, 2)])[1:120, ]
#' time <- smart$TEVENT[1:120]
#' event <- smart$EVENT[1:120]
#' y <- survival::Surv(time, event)
#'
#' fit <- fit_mnet(
#'   x, y,
#'   nfolds = 3,
#'   gammas = 3, alphas = c(0.3, 0.8),
#'   max.iter = 15000, seed = 1010
#' )
#'
#' nom <- as_nomogram(
#'   fit, x, time, event, pred.at = 365 * 2,
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
  call <- match.call()

  mnet_cv <- ncvreg_tune_gamma_alpha(
    x, y,
    penalty = "MCP",
    nfolds = nfolds,
    gammas = gammas, alphas = alphas,
    eps = eps, max.iter = max.iter,
    seed = seed, trace = trace,
    parallel = parallel
  )

  gamma_opt <- mnet_cv$best.gamma
  alpha_opt <- mnet_cv$best.alpha
  lambda_opt <- mnet_cv$best.model$lambda.min

  mnet_full <- ncvreg::ncvsurv(
    x, y,
    penalty = "MCP",
    gamma = gamma_opt,
    alpha = alpha_opt,
    lambda = lambda_opt,
    eps = eps, max.iter = max.iter
  )

  if (all(abs(mnet_full$beta[-1L, ]) < .Machine$double.eps)) {
    stop("Null model produced by the full fit (all coefficients are zero). Please try to tune gammas, alphas, seed, nfolds, or increase sample size.")
  }

  model <- list(
    "model" = mnet_full,
    "gamma" = gamma_opt,
    "alpha" = alpha_opt,
    "lambda" = lambda_opt,
    "type" = "mnet",
    "seed" = seed,
    "call" = call
  )

  class(model) <- c("hdnom.model", "hdnom.model.mnet")
  model
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
#' data("smart")
#' x <- as.matrix(smart[, -c(1, 2)])[1:120, ]
#' time <- smart$TEVENT[1:120]
#' event <- smart$EVENT[1:120]
#' y <- survival::Surv(time, event)
#'
#' fit <- fit_flasso(x, y,
#'   lambda1 = c(1, 10), lambda2 = c(0.01),
#'   nfolds = 3, seed = 11
#' )
#'
#' nom <- as_nomogram(
#'   fit, x, time, event, pred.at = 365 * 2,
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
  call <- match.call()

  if (trace) cat("Starting cross-validation...\n")
  flasso_cv <- penalized_tune_lambda(
    response = y, penalized = x, fold = nfolds,
    lambda1 = lambda1, lambda2 = lambda2,
    maxiter = maxiter, epsilon = epsilon,
    seed = seed, trace = trace, parallel = parallel,
    fusedl = TRUE, standardize = TRUE, model = "cox", ...
  )

  lambda1_opt <- flasso_cv$"best.lambda1"
  lambda2_opt <- flasso_cv$"best.lambda2"

  if (trace) cat("Fitting fused lasso model with full data...\n")
  flasso_full <- penalized(
    response = y, penalized = x,
    lambda1 = lambda1_opt,
    lambda2 = lambda2_opt,
    maxiter = maxiter, epsilon = epsilon,
    trace = trace,
    fusedl = TRUE, standardize = FALSE, model = "cox", ...
  )

  if (all(abs(flasso_full@penalized) < .Machine$double.eps)) {
    stop("Null model produced by the full fit (all coefficients are zero). Please try changing the seed, nfolds, or increase sample size.")
  }

  model <- list(
    "model" = flasso_full,
    "lambda1" = lambda1_opt,
    "lambda2" = lambda2_opt,
    "type" = "flasso",
    "seed" = seed,
    "call" = call
  )

  class(model) <- c("hdnom.model", "hdnom.model.flasso")
  model
}
