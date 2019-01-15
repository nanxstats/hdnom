# convert model object to nomogram.raw() input model object
#
# @param model glmnet, ncvreg or penfit model object, e.g. fit$model
# @param x predictor matrix (original, with all variables)
#
# @return converted model object, can be used as input
# for as_nomogram_raw() to construct nomogram object
convert_model <- function(model, x) {
  is_glmnet <- function(model) "glmnet" %in% class(model)
  is_ncvreg <- function(model) "ncvreg" %in% class(model)
  is_penalized <- function(model) "penfit" %in% class(model)

  res <- list(
    "coefficients" = NULL, "non.slopes" = NULL,
    "linear.predictors" = NULL, "fitted.values" = NULL, "Design" = NULL,
    "x" = NULL, "stats" = NULL, "df.residual" = NULL,
    "assign" = NULL, "sformula" = NULL
  )

  # data required by nomogram.raw()

  # $coefficients
  if (is_glmnet(model) | is_ncvreg(model)) res$coefficients <- c("Intercept" = 0.0, model$beta[as.vector(model$beta) != 0, 1L])
  if (is_penalized(model)) res$coefficients <- c("Intercept" = 0.0, model@penalized[model@penalized != 0])

  # $non.slopes
  res$non.slopes <- 1

  # $linear.predictors
  if (is_glmnet(model)) {
    res$linear.predictors <- as.vector(predict(model, newx = x, s = model$"lambda", type = "link"))
  }
  if (is_ncvreg(model)) {
    res$linear.predictors <- predict(model, X = x, type = "link")
  }
  if (is_penalized(model)) {
    res$linear.predictors <- as.numeric(model@lin.pred)
  }

  # $fitted.values
  res$fitted.values <- res$linear.predictors

  # $Design
  res$Design <- list(
    "name" = NULL, "label" = NULL, "units" = NULL, "colnames" = NULL,
    "mmcolnames" = NULL, "assume" = NULL, "assume.code" = NULL,
    "parms" = NULL, "limits" = NULL, "values" = NULL,
    "nonlinear" = NULL, "interactions" = NULL
  )

  # focus on non-zero variables
  if (is_glmnet(model) | is_ncvreg(model)) idx_nzv <- which(as.vector(model$beta) != 0)
  if (is_penalized(model)) idx_nzv <- which(model@penalized != 0)
  x_nzv <- as.data.frame(x[, idx_nzv, drop = FALSE])
  n_nzv <- length(idx_nzv)

  # $Design$limits
  res$Design$limits <- rms_datadist(x_nzv)$limits

  # $Design$... others
  res$Design$name <- names(x_nzv)
  res$Design$label <- names(x_nzv)
  res$Design$units <- rep("", n_nzv)
  names(res$Design$units) <- names(x_nzv)
  res$Design$colnames <- names(x_nzv)
  res$Design$mmcolnames <- names(x_nzv)
  res$Design$assume <- rep("asis", n_nzv)
  res$Design$assume.code <- rep(1L, n_nzv)
  res$Design$parms <- list()
  res$Design$values <- structure(list(), .Names = character(0))
  res$Design$nonlinear <- vector("list", n_nzv)
  names(res$Design$nonlinear) <- names(x_nzv)
  for (i in 1L:n_nzv) res$Design$nonlinear[[i]] <- FALSE

  # data required by rms_predict()

  # $x
  res$x <- x[, idx_nzv, drop = FALSE]

  # $stats["Sigma"]
  res$stats <- c("Sigma" = 1)

  # $df.residual
  res$df.residual <- nrow(x) - n_nzv - 1L

  # $assign
  names_nzv <- colnames(x)[idx_nzv]
  res$assign <- as.list(2L:(n_nzv + 1L))
  names(res$assign) <- names_nzv

  # $sformula
  res$sformula <- as.formula(paste("reponse_name_to_remove ~ ", paste(names_nzv, collapse = " + ")))

  res
}
