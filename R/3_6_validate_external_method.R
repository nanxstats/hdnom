#' Print external validation results
#'
#' Print external validation results
#'
#' @param x An object returned by \code{\link{validate_external}}.
#' @param ... Other parameters (not used).
#'
#' @method print hdnom.validate.external
#'
#' @export
#'
#' @examples
#' NULL
print.hdnom.validate.external <- function(x, ...) {
  method <- setdiff(class(x), "hdnom.validate.external")

  switch(

    method,

    glmnet.validate.external = {
      cat("High-Dimensional Cox Model External Validation Object\n")
      cat("Model type:", attr(x, "model.type"), "\n")
      cat("Time-dependent AUC type:", attr(x, "tauc.type"), "\n")
      cat("Evaluation time points for tAUC:", attr(x, "tauc.time"))
    },

    ncvreg.validate.external = {
      cat("High-Dimensional Cox Model External Validation Object\n")
      cat("Model type:", attr(x, "model.type"), "\n")
      cat("Time-dependent AUC type:", attr(x, "tauc.type"), "\n")
      cat("Evaluation time points for tAUC:", attr(x, "tauc.time"))
    },

    penalized.validate.external = {
      cat("High-Dimensional Cox Model External Validation Object\n")
      cat("Model type:", attr(x, "model.type"), "\n")
      cat("Time-dependent AUC type:", attr(x, "tauc.type"), "\n")
      cat("Evaluation time points for tAUC:", attr(x, "tauc.time"))
    }
  )

  invisible(x)
}

#' Summary of external validation results
#'
#' Summary of external validation results
#'
#' @param object An object returned by \code{\link{validate_external}}.
#' @param silent Print summary table header or not,
#' default is \code{FALSE}.
#' @param ... Other parameters (not used).
#'
#' @method summary hdnom.validate.external
#'
#' @export
#'
#' @examples
#' NULL
summary.hdnom.validate.external <- function(object, silent = FALSE, ...) {
  tauc.time <- attr(object, "tauc.time")
  aucmat <- matrix(NA, ncol = length(tauc.time), nrow = 1L)
  aucmat[1L, ] <- object$"auc"
  rownames(aucmat) <- "AUC"
  colnames(aucmat) <- tauc.time

  if (!silent) cat("Time-Dependent AUC Summary at Evaluation Time Points\n")
  print(aucmat)

  invisible(aucmat)
}

#' Plot time-dependent discrimination curves for external validation
#'
#' Plot time-dependent discrimination curves for external validation
#'
#' @param x An object returned by \code{\link{validate_external}}.
#' @param col.pal Color palette to use. Possible values are
#' \code{"JCO"}, \code{"Lancet"}, \code{"NPG"}, and \code{"AAAS"}.
#' Default is \code{"JCO"}.
#' @param ylim Range of y coordinates. For example, \code{c(0.5, 1)}.
#' @param ... Other parameters (not used).
#'
#' @method plot hdnom.validate.external
#'
#' @export
#'
#' @importFrom ggplot2 ggplot aes_string geom_point geom_line geom_point
#' geom_ribbon scale_x_continuous scale_fill_manual scale_colour_manual
#' theme_bw theme ylab coord_cartesian
#'
#' @examples
#' NULL
plot.hdnom.validate.external <- function(
  x, col.pal = c("JCO", "Lancet", "NPG", "AAAS"), ylim = NULL, ...) {
  df <- as.data.frame(t(summary(x, silent = TRUE)))
  tauc_time <- attr(x, "tauc.time")

  df[, "Time"] <- tauc_time

  col.pal <- match.arg(col.pal)
  col_pal <- switch(
    col.pal,
    JCO = palette_jco()[1], Lancet = palette_lancet()[1],
    NPG = palette_npg()[1], AAAS = palette_aaas()[1]
  )

  ggplot(data = df, aes_string(x = "Time", y = "AUC")) +
    geom_point(colour = col_pal) +
    geom_line(colour = col_pal) +
    scale_x_continuous(breaks = df$"Time") +
    coord_cartesian(ylim = ylim) +
    theme_bw() +
    theme(legend.position = "none") +
    ylab("Area under ROC")
}
