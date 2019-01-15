#' Print calibration results
#'
#' Print calibration results
#'
#' @param x An object returned by \code{\link{calibrate}}.
#' @param ... Other parameters (not used).
#'
#' @method print hdnom.calibrate
#'
#' @export
#'
#' @examples
#' NULL
print.hdnom.calibrate <- function(x, ...) {
  method <- setdiff(class(x), "hdnom.calibrate")

  switch(

    method,

    glmnet.calibrate.fitting = {
      cat("High-Dimensional Cox Model Calibration Object\n")
      cat("Random seed:", attr(x, "seed"), "\n")
      cat("Calibration method: fitting\n")
      cat("Model type:", attr(x, "model.type"), "\n")
      cat("glmnet model alpha:", attr(x, "alpha"), "\n")
      cat("glmnet model lambda:", attr(x, "lambda"), "\n")
      if (is.null(attr(x, "pen.factor"))) {
        cat("glmnet model penalty factor: not specified\n")
      } else {
        cat("glmnet model penalty factor: specified\n")
      }
      cat("Calibration time point:", attr(x, "pred.at"), "\n")
      cat("Number of groups formed for calibration:", attr(x, "ngroup"), "\n")
    },

    glmnet.calibrate.bootstrap = {
      cat("High-Dimensional Cox Model Calibration Object\n")
      cat("Random seed:", attr(x, "seed"), "\n")
      cat("Calibration method: bootstrap\n")
      cat("Bootstrap samples:", attr(x, "boot.times"), "\n")
      cat("Model type:", attr(x, "model.type"), "\n")
      cat("glmnet model alpha:", attr(x, "alpha"), "\n")
      cat("glmnet model lambda:", attr(x, "lambda"), "\n")
      if (is.null(attr(x, "pen.factor"))) {
        cat("glmnet model penalty factor: not specified\n")
      } else {
        cat("glmnet model penalty factor: specified\n")
      }
      cat("Calibration time point:", attr(x, "pred.at"), "\n")
      cat("Number of groups formed for calibration:", attr(x, "ngroup"), "\n")
    },

    glmnet.calibrate.cv = {
      cat("High-Dimensional Cox Model Calibration Object\n")
      cat("Random seed:", attr(x, "seed"), "\n")
      cat("Calibration method: k-fold cross-validation\n")
      cat("Cross-validation folds:", attr(x, "nfolds"), "\n")
      cat("Model type:", attr(x, "model.type"), "\n")
      cat("glmnet model alpha:", attr(x, "alpha"), "\n")
      cat("glmnet model lambda:", attr(x, "lambda"), "\n")
      if (is.null(attr(x, "pen.factor"))) {
        cat("glmnet model penalty factor: not specified\n")
      } else {
        cat("glmnet model penalty factor: specified\n")
      }
      cat("Calibration time point:", attr(x, "pred.at"), "\n")
      cat("Number of groups formed for calibration:", attr(x, "ngroup"), "\n")
    },

    glmnet.calibrate.repeated.cv = {
      cat("High-Dimensional Cox Model Calibration Object\n")
      cat("Random seed:", attr(x, "seed"), "\n")
      cat("Calibration method: repeated cross-validation\n")
      cat("Cross-validation folds:", attr(x, "nfolds"), "\n")
      cat("Cross-validation repeated times:", attr(x, "rep.times"), "\n")
      cat("Model type:", attr(x, "model.type"), "\n")
      cat("glmnet model alpha:", attr(x, "alpha"), "\n")
      cat("glmnet model lambda:", attr(x, "lambda"), "\n")
      if (is.null(attr(x, "pen.factor"))) {
        cat("glmnet model penalty factor: not specified\n")
      } else {
        cat("glmnet model penalty factor: specified\n")
      }
      cat("Calibration time point:", attr(x, "pred.at"), "\n")
      cat("Number of groups formed for calibration:", attr(x, "ngroup"), "\n")
    },

    ncvreg.calibrate.fitting = {
      cat("High-Dimensional Cox Model Calibration Object\n")
      cat("Random seed:", attr(x, "seed"), "\n")
      cat("Calibration method: fitting\n")
      cat("Model type:", attr(x, "model.type"), "\n")
      cat("ncvreg model gamma:", attr(x, "gamma"), "\n")
      cat("ncvreg model alpha:", attr(x, "alpha"), "\n")
      cat("ncvreg model lambda:", attr(x, "lambda"), "\n")
      cat("Calibration time point:", attr(x, "pred.at"), "\n")
      cat("Number of groups formed for calibration:", attr(x, "ngroup"), "\n")
    },

    ncvreg.calibrate.bootstrap = {
      cat("High-Dimensional Cox Model Calibration Object\n")
      cat("Random seed:", attr(x, "seed"), "\n")
      cat("Calibration method: bootstrap\n")
      cat("Bootstrap samples:", attr(x, "boot.times"), "\n")
      cat("Model type:", attr(x, "model.type"), "\n")
      cat("ncvreg model gamma:", attr(x, "gamma"), "\n")
      cat("ncvreg model alpha:", attr(x, "alpha"), "\n")
      cat("ncvreg model lambda:", attr(x, "lambda"), "\n")
      cat("Calibration time point:", attr(x, "pred.at"), "\n")
      cat("Number of groups formed for calibration:", attr(x, "ngroup"), "\n")
    },

    ncvreg.calibrate.cv = {
      cat("High-Dimensional Cox Model Calibration Object\n")
      cat("Random seed:", attr(x, "seed"), "\n")
      cat("Calibration method: k-fold cross-validation\n")
      cat("Cross-validation folds:", attr(x, "nfolds"), "\n")
      cat("Model type:", attr(x, "model.type"), "\n")
      cat("ncvreg model gamma:", attr(x, "gamma"), "\n")
      cat("ncvreg model alpha:", attr(x, "alpha"), "\n")
      cat("ncvreg model lambda:", attr(x, "lambda"), "\n")
      cat("Calibration time point:", attr(x, "pred.at"), "\n")
      cat("Number of groups formed for calibration:", attr(x, "ngroup"), "\n")
    },

    ncvreg.calibrate.repeated.cv = {
      cat("High-Dimensional Cox Model Calibration Object\n")
      cat("Random seed:", attr(x, "seed"), "\n")
      cat("Calibration method: repeated cross-validation\n")
      cat("Cross-validation folds:", attr(x, "nfolds"), "\n")
      cat("Cross-validation repeated times:", attr(x, "rep.times"), "\n")
      cat("Model type:", attr(x, "model.type"), "\n")
      cat("ncvreg model alpha:", attr(x, "gamma"), "\n")
      cat("ncvreg model alpha:", attr(x, "alpha"), "\n")
      cat("ncvreg model lambda:", attr(x, "lambda"), "\n")
      cat("Calibration time point:", attr(x, "pred.at"), "\n")
      cat("Number of groups formed for calibration:", attr(x, "ngroup"), "\n")
    },

    penalized.calibrate.fitting = {
      cat("High-Dimensional Cox Model Calibration Object\n")
      cat("Random seed:", attr(x, "seed"), "\n")
      cat("Calibration method: fitting\n")
      cat("Model type:", attr(x, "model.type"), "\n")
      cat("Fused lasso model lambda1:", attr(x, "lambda1"), "\n")
      cat("Fused lasso model lambda2:", attr(x, "lambda2"), "\n")
      cat("Calibration time point:", attr(x, "pred.at"), "\n")
      cat("Number of groups formed for calibration:", attr(x, "ngroup"), "\n")
    },

    penalized.calibrate.bootstrap = {
      cat("High-Dimensional Cox Model Calibration Object\n")
      cat("Random seed:", attr(x, "seed"), "\n")
      cat("Calibration method: bootstrap\n")
      cat("Bootstrap samples:", attr(x, "boot.times"), "\n")
      cat("Model type:", attr(x, "model.type"), "\n")
      cat("Fused lasso model lambda1:", attr(x, "lambda1"), "\n")
      cat("Fused lasso model lambda2:", attr(x, "lambda2"), "\n")
      cat("Calibration time point:", attr(x, "pred.at"), "\n")
      cat("Number of groups formed for calibration:", attr(x, "ngroup"), "\n")
    },

    penalized.calibrate.cv = {
      cat("High-Dimensional Cox Model Calibration Object\n")
      cat("Random seed:", attr(x, "seed"), "\n")
      cat("Calibration method: k-fold cross-validation\n")
      cat("Cross-validation folds:", attr(x, "nfolds"), "\n")
      cat("Model type:", attr(x, "model.type"), "\n")
      cat("Fused lasso model lambda1:", attr(x, "lambda1"), "\n")
      cat("Fused lasso model lambda2:", attr(x, "lambda2"), "\n")
      cat("Calibration time point:", attr(x, "pred.at"), "\n")
      cat("Number of groups formed for calibration:", attr(x, "ngroup"), "\n")
    },

    penalized.calibrate.repeated.cv = {
      cat("High-Dimensional Cox Model Calibration Object\n")
      cat("Random seed:", attr(x, "seed"), "\n")
      cat("Calibration method: repeated cross-validation\n")
      cat("Cross-validation folds:", attr(x, "nfolds"), "\n")
      cat("Cross-validation repeated times:", attr(x, "rep.times"), "\n")
      cat("Model type:", attr(x, "model.type"), "\n")
      cat("Fused lasso model lambda1:", attr(x, "lambda1"), "\n")
      cat("Fused lasso model lambda2:", attr(x, "lambda2"), "\n")
      cat("Calibration time point:", attr(x, "pred.at"), "\n")
      cat("Number of groups formed for calibration:", attr(x, "ngroup"), "\n")
    }
  )

  invisible(x)
}

#' Summary of calibration results
#'
#' Summary of calibration results
#'
#' @param object An object returned by \code{\link{calibrate}}.
#' @param ... Other parameters (not used).
#'
#' @method summary hdnom.calibrate
#'
#' @export
#'
#' @examples
#' NULL
summary.hdnom.calibrate <- function(object, ...) {
  method <- setdiff(class(object), "hdnom.calibrate")

  x <- object

  switch(

    method,

    glmnet.calibrate.fitting = {
      attr(x, "alpha") <- NULL
      attr(x, "lambda") <- NULL
      attr(x, "pen.factor") <- NULL
    },

    glmnet.calibrate.bootstrap = {
      attr(x, "boot.times") <- NULL
      attr(x, "alpha") <- NULL
      attr(x, "lambda") <- NULL
      attr(x, "pen.factor") <- NULL
    },

    glmnet.calibrate.cv = {
      attr(x, "nfolds") <- NULL
      attr(x, "alpha") <- NULL
      attr(x, "lambda") <- NULL
      attr(x, "pen.factor") <- NULL
    },

    glmnet.calibrate.repeated.cv = {
      attr(x, "nfolds") <- NULL
      attr(x, "rep.times") <- NULL
      attr(x, "alpha") <- NULL
      attr(x, "lambda") <- NULL
      attr(x, "pen.factor") <- NULL
    },

    ncvreg.calibrate.fitting = {
      attr(x, "alpha") <- NULL
      attr(x, "lambda") <- NULL
      attr(x, "gamma") <- NULL
    },

    ncvreg.calibrate.bootstrap = {
      attr(x, "alpha") <- NULL
      attr(x, "lambda") <- NULL
      attr(x, "gamma") <- NULL
      attr(x, "boot.times") <- NULL
    },

    ncvreg.calibrate.cv = {
      attr(x, "alpha") <- NULL
      attr(x, "lambda") <- NULL
      attr(x, "gamma") <- NULL
      attr(x, "nfolds") <- NULL
    },

    ncvreg.calibrate.repeated.cv = {
      attr(x, "alpha") <- NULL
      attr(x, "lambda") <- NULL
      attr(x, "gamma") <- NULL
      attr(x, "nfolds") <- NULL
      attr(x, "rep.times") <- NULL
    },

    penalized.calibrate.fitting = {
      attr(x, "lambda1") <- NULL
      attr(x, "lambda2") <- NULL
    },

    penalized.calibrate.bootstrap = {
      attr(x, "lambda1") <- NULL
      attr(x, "lambda2") <- NULL
      attr(x, "boot.times") <- NULL
    },

    penalized.calibrate.cv = {
      attr(x, "lambda1") <- NULL
      attr(x, "lambda2") <- NULL
      attr(x, "nfolds") <- NULL
    },

    penalized.calibrate.repeated.cv = {
      attr(x, "lambda1") <- NULL
      attr(x, "lambda2") <- NULL
      attr(x, "nfolds") <- NULL
      attr(x, "rep.times") <- NULL
    }
  )

  attr(x, "model.type") <- NULL
  attr(x, "pred.at") <- NULL
  attr(x, "ngroup") <- NULL
  attr(x, "risk.group") <- NULL
  attr(x, "surv.time") <- NULL
  attr(x, "surv.event") <- NULL
  attr(x, "seed") <- NULL

  cat("  Calibration Summary Table\n")
  class(x) <- "matrix"
  print(x)

  invisible(x)
}

#' Plot calibration results
#'
#' Plot calibration results
#'
#' @param x An object returned by \code{\link{calibrate}}.
#' @param xlim x axis limits of the plot.
#' @param ylim y axis limits of the plot.
#' @param col.pal Color palette to use. Possible values are
#' \code{"JCO"}, \code{"Lancet"}, \code{"NPG"}, and \code{"AAAS"}.
#' Default is \code{"JCO"}.
#' @param ... Other parameters for \code{plot}.
#'
#' @method plot hdnom.calibrate
#'
#' @export
#'
#' @importFrom ggplot2 ggplot aes_string geom_errorbar
#' geom_line geom_point geom_abline xlab ylab
#'
#' @examples
#' NULL
plot.hdnom.calibrate <- function(
  x, xlim = c(0, 1), ylim = c(0, 1),
  col.pal = c("JCO", "Lancet", "NPG", "AAAS"), ...) {
  df <- data.frame(
    "pre" = x[, "Predicted"], "obs" = x[, "Observed"],
    "ll" = x[, "Lower 95%"], "ul" = x[, "Upper 95%"]
  )

  col.pal <- match.arg(col.pal)
  col_pal <- switch(
    col.pal,
    JCO = palette_jco()[1], Lancet = palette_lancet()[1],
    NPG = palette_npg()[1], AAAS = palette_aaas()[1]
  )

  ggplot(
    df,
    aes_string(
      x = "pre", y = "obs",
      xmin = xlim[1L], xmax = xlim[2L],
      ymin = ylim[1L], ymax = ylim[2L]
    )
  ) +
    geom_abline(slope = 1, intercept = 0, colour = "grey") +
    geom_errorbar(aes_string(ymin = "ll", ymax = "ul"), colour = col_pal) +
    geom_line(colour = col_pal) +
    geom_point(size = 3, colour = col_pal) +
    xlab("Predicted Survival Probability") +
    ylab("Observed Survival Probability") +
    theme_hdnom()
}
