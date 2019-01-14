#' Print validation results
#'
#' Print validation results
#'
#' @param x An object returned by \code{\link{validate}}.
#' @param ... Other parameters (not used).
#'
#' @method print hdnom.validate
#'
#' @export
#'
#' @examples
#' NULL
print.hdnom.validate <- function(x, ...) {
  method <- setdiff(class(x), "hdnom.validate")

  switch(

    method,

    glmnet.validate.bootstrap = {
      cat("High-Dimensional Cox Model Validation Object\n")
      cat("Random seed:", attr(x, "seed"), "\n")
      cat("Validation method: bootstrap\n")
      cat("Bootstrap samples:", attr(x, "boot.times"), "\n")
      cat("Model type:", attr(x, "model.type"), "\n")
      cat("glmnet model alpha:", attr(x, "alpha"), "\n")
      cat("glmnet model lambda:", attr(x, "lambda"), "\n")
      if (is.null(attr(x, "pen.factor"))) {
        cat("glmnet model penalty factor: not specified\n")
      } else {
        cat("glmnet model penalty factor: specified\n")
      }
      cat("Time-dependent AUC type:", attr(x, "tauc.type"), "\n")
      cat("Evaluation time points for tAUC:", attr(x, "tauc.time"))
    },

    glmnet.validate.cv = {
      cat("High-Dimensional Cox Model Validation Object\n")
      cat("Random seed:", attr(x, "seed"), "\n")
      cat("Validation method: k-fold cross-validation\n")
      cat("Cross-validation folds:", attr(x, "nfolds"), "\n")
      cat("Model type:", attr(x, "model.type"), "\n")
      cat("glmnet model alpha:", attr(x, "alpha"), "\n")
      cat("glmnet model lambda:", attr(x, "lambda"), "\n")
      if (is.null(attr(x, "pen.factor"))) {
        cat("glmnet model penalty factor: not specified\n")
      } else {
        cat("glmnet model penalty factor: specified\n")
      }
      cat("Time-dependent AUC type:", attr(x, "tauc.type"), "\n")
      cat("Evaluation time points for tAUC:", attr(x, "tauc.time"))
    },

    glmnet.validate.repeated.cv = {
      cat("High-Dimensional Cox Model Validation Object\n")
      cat("Random seed:", attr(x, "seed"), "\n")
      cat("Validation method: repeated cross-validation\n")
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
      cat("Time-dependent AUC type:", attr(x, "tauc.type"), "\n")
      cat("Evaluation time points for tAUC:", attr(x, "tauc.time"))
    },

    ncvreg.validate.bootstrap = {
      cat("High-Dimensional Cox Model Validation Object\n")
      cat("Random seed:", attr(x, "seed"), "\n")
      cat("Validation method: bootstrap\n")
      cat("Bootstrap samples:", attr(x, "boot.times"), "\n")
      cat("Model type:", attr(x, "model.type"), "\n")
      cat("ncvreg model gamma:", attr(x, "gamma"), "\n")
      cat("ncvreg model alpha:", attr(x, "alpha"), "\n")
      cat("ncvreg model lambda:", attr(x, "lambda"), "\n")
      cat("Time-dependent AUC type:", attr(x, "tauc.type"), "\n")
      cat("Evaluation time points for tAUC:", attr(x, "tauc.time"))
    },

    ncvreg.validate.cv = {
      cat("High-Dimensional Cox Model Validation Object\n")
      cat("Random seed:", attr(x, "seed"), "\n")
      cat("Validation method: k-fold cross-validation\n")
      cat("Cross-validation folds:", attr(x, "nfolds"), "\n")
      cat("Model type:", attr(x, "model.type"), "\n")
      cat("ncvreg model gamma:", attr(x, "gamma"), "\n")
      cat("ncvreg model alpha:", attr(x, "alpha"), "\n")
      cat("ncvreg model lambda:", attr(x, "lambda"), "\n")
      cat("Time-dependent AUC type:", attr(x, "tauc.type"), "\n")
      cat("Evaluation time points for tAUC:", attr(x, "tauc.time"))
    },

    ncvreg.validate.repeated.cv = {
      cat("High-Dimensional Cox Model Validation Object\n")
      cat("Random seed:", attr(x, "seed"), "\n")
      cat("Validation method: repeated cross-validation\n")
      cat("Cross-validation folds:", attr(x, "nfolds"), "\n")
      cat("Cross-validation repeated times:", attr(x, "rep.times"), "\n")
      cat("Model type:", attr(x, "model.type"), "\n")
      cat("ncvreg model gamma:", attr(x, "gamma"), "\n")
      cat("ncvreg model alpha:", attr(x, "alpha"), "\n")
      cat("ncvreg model lambda:", attr(x, "lambda"), "\n")
      cat("Time-dependent AUC type:", attr(x, "tauc.type"), "\n")
      cat("Evaluation time points for tAUC:", attr(x, "tauc.time"))
    },

    penalized.validate.bootstrap = {
      cat("High-Dimensional Cox Model Validation Object\n")
      cat("Random seed:", attr(x, "seed"), "\n")
      cat("Validation method: bootstrap\n")
      cat("Bootstrap samples:", attr(x, "boot.times"), "\n")
      cat("Model type:", attr(x, "model.type"), "\n")
      cat("Fused lasso model lambda1:", attr(x, "lambda1"), "\n")
      cat("Fused lasso model lambda2:", attr(x, "lambda2"), "\n")
      cat("Time-dependent AUC type:", attr(x, "tauc.type"), "\n")
      cat("Evaluation time points for tAUC:", attr(x, "tauc.time"))
    },

    penalized.validate.cv = {
      cat("High-Dimensional Cox Model Validation Object\n")
      cat("Random seed:", attr(x, "seed"), "\n")
      cat("Validation method: k-fold cross-validation\n")
      cat("Cross-validation folds:", attr(x, "nfolds"), "\n")
      cat("Model type:", attr(x, "model.type"), "\n")
      cat("Fused lasso model lambda1:", attr(x, "lambda1"), "\n")
      cat("Fused lasso model lambda2:", attr(x, "lambda2"), "\n")
      cat("Time-dependent AUC type:", attr(x, "tauc.type"), "\n")
      cat("Evaluation time points for tAUC:", attr(x, "tauc.time"))
    },

    penalized.validate.repeated.cv = {
      cat("High-Dimensional Cox Model Validation Object\n")
      cat("Random seed:", attr(x, "seed"), "\n")
      cat("Validation method: repeated cross-validation\n")
      cat("Cross-validation folds:", attr(x, "nfolds"), "\n")
      cat("Cross-validation repeated times:", attr(x, "rep.times"), "\n")
      cat("Model type:", attr(x, "model.type"), "\n")
      cat("Fused lasso model lambda1:", attr(x, "lambda1"), "\n")
      cat("Fused lasso model lambda2:", attr(x, "lambda2"), "\n")
      cat("Time-dependent AUC type:", attr(x, "tauc.type"), "\n")
      cat("Evaluation time points for tAUC:", attr(x, "tauc.time"))
    }
  )

  invisible(x)
}

#' Summary of validation results
#'
#' Summary of validation results
#'
#' @param object A \code{\link{validate}} object.
#' @param silent Print summary table header or not,
#' default is \code{FALSE}.
#' @param ... Other parameters (not used).
#'
#' @method summary hdnom.validate
#'
#' @export
#'
#' @examples
#' NULL
summary.hdnom.validate <- function(object, silent = FALSE, ...) {
  method <- setdiff(class(object), "hdnom.validate")

  if (grepl("validate.bootstrap", method)) {
    boot.times <- attr(object, "boot.times")
    tauc.time <- attr(object, "tauc.time")
    aucmat <- matrix(NA, ncol = length(tauc.time), nrow = boot.times)
    for (i in 1L:boot.times) aucmat[i, ] <- object[[i]]$auc
    summary_mat <- rbind(apply(aucmat, 2, mean), apply(aucmat, 2, quantile))
    rownames(summary_mat) <- c(
      "Mean", "Min", "0.25 Qt.",
      "Median", "0.75 Qt.", "Max"
    )
    colnames(summary_mat) <- tauc.time
  } else if (grepl("validate.cv", method)) {
    nfolds <- attr(object, "nfolds")
    tauc.time <- attr(object, "tauc.time")
    aucmat <- matrix(NA, ncol = length(tauc.time), nrow = nfolds)
    for (i in 1L:nfolds) aucmat[i, ] <- object[[i]]$auc
    summary_mat <- rbind(apply(aucmat, 2, mean), apply(aucmat, 2, quantile))
    rownames(summary_mat) <- c(
      "Mean", "Min", "0.25 Qt.",
      "Median", "0.75 Qt.", "Max"
    )
    colnames(summary_mat) <- tauc.time
  } else if (grepl("validate.repeated.cv", method)) {
    nfolds <- attr(object, "nfolds")
    rep.times <- attr(object, "rep.times")
    tauc.time <- attr(object, "tauc.time")
    auclist <- vector("list", rep.times)
    for (i in 1L:rep.times) {
      auclist[[i]] <- matrix(NA, ncol = length(tauc.time), nrow = nfolds)
    }
    for (i in 1L:rep.times) {
      for (j in 1L:nfolds) {
        auclist[[i]][j, ] <- object[[i]][[j]]$auc
      }
    }

    summary_list <- vector("list", rep.times)
    for (i in 1L:rep.times) {
      summary_list[[i]] <- rbind(
        apply(auclist[[i]], 2, mean),
        apply(auclist[[i]], 2, quantile)
      )
    }

    summary_mat <- Reduce("+", summary_list) / length(summary_list)
    rownames(summary_mat) <- c(
      "Mean of Mean", "Mean of Min",
      "Mean of 0.25 Qt.", "Mean of Median",
      "Mean of 0.75 Qt.", "Mean of Max"
    )
    colnames(summary_mat) <- tauc.time

    if (!silent) {
      cat("Note: for repeated CV, we evaluated quantile statistic tables for\n")
      cat("each CV repeat, then calculated element-wise mean across all tables.\n")
    }
  } else {
    stop("hdnom.validate object is not valid")
  }

  if (!silent) cat("Time-Dependent AUC Summary at Evaluation Time Points\n")
  print(summary_mat)

  invisible(summary_mat)
}

#' Plot optimism-corrected time-dependent discrimination curves for validation
#'
#' Plot optimism-corrected time-dependent discrimination curves for validation
#'
#' @param x An object returned by \code{\link{validate}}.
#' @param col.pal Color palette to use. Possible values are
#' \code{"JCO"}, \code{"Lancet"}, \code{"NPG"}, and \code{"AAAS"}.
#' Default is \code{"JCO"}.
#' @param ylim Range of y coordinates. For example, \code{c(0.5, 1)}.
#' @param ... Other parameters (not used).
#'
#' @method plot hdnom.validate
#'
#' @export
#'
#' @importFrom ggplot2 ggplot aes_string geom_point geom_line geom_point
#' geom_ribbon scale_x_continuous scale_fill_manual scale_colour_manual
#' theme_bw theme ylab coord_cartesian
#'
#' @examples
#' NULL
plot.hdnom.validate <- function(
  x, col.pal = c("JCO", "Lancet", "NPG", "AAAS"), ylim = NULL, ...) {
  df <- as.data.frame(t(summary(x, silent = TRUE)))
  tauc_time <- attr(x, "tauc.time")

  # special processing for repeated cv
  if (any(grepl(pattern = "validate.repeated.cv", class(x)))) {
    names(df) <- sapply(strsplit(names(df), "Mean of "), "[", 2L)
  }

  df[, "Time"] <- tauc_time
  names(df)[which(names(df) == "0.25 Qt.")] <- "Qt25"
  names(df)[which(names(df) == "0.75 Qt.")] <- "Qt75"

  col.pal <- match.arg(col.pal)
  col_pal <- switch(
    col.pal,
    JCO = palette_jco()[1], Lancet = palette_lancet()[1],
    NPG = palette_npg()[1], AAAS = palette_aaas()[1]
  )

  ggplot(data = df, aes_string(x = "Time", y = "Mean")) +
    geom_point(colour = col_pal) +
    geom_line(colour = col_pal) +
    geom_point(
      data = df, aes_string(x = "Time", y = "Median"),
      colour = col_pal
    ) +
    geom_line(
      data = df, aes_string(x = "Time", y = "Median"),
      colour = col_pal, linetype = "dashed"
    ) +
    geom_ribbon(
      data = df, aes_string(ymin = "Qt25", ymax = "Qt75"),
      linetype = 0, alpha = 0.2
    ) +
    geom_ribbon(
      data = df, aes_string(ymin = "Min", ymax = "Max"),
      linetype = 0, alpha = 0.1
    ) +
    scale_x_continuous(breaks = df$"Time") +
    coord_cartesian(ylim = ylim) +
    theme_bw() +
    theme(legend.position = "none") +
    ylab("Area under ROC")
}
