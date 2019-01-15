#' Print external calibration results
#'
#' Print external calibration results
#'
#' @param x An object returned by \code{\link{calibrate_external}}.
#' @param ... Other parameters (not used).
#'
#' @method print hdnom.calibrate.external
#'
#' @export
#'
#' @examples
#' NULL
print.hdnom.calibrate.external <- function(x, ...) {
  method <- setdiff(class(x), "hdnom.calibrate.external")

  cat("High-Dimensional Cox Model External Calibration Object\n")
  cat("Model type:", attr(x, "model.type"), "\n")
  cat("Calibration time point:", attr(x, "pred.at"), "\n")
  cat("Number of groups formed for calibration:", attr(x, "ngroup"), "\n")

  invisible(x)
}

#' Summary of external calibration results
#'
#' Summary of external calibration results
#'
#' @param object An object returned by \code{\link{calibrate_external}}.
#' @param ... Other parameters (not used).
#'
#' @method summary hdnom.calibrate.external
#'
#' @export
#'
#' @examples
#' NULL
summary.hdnom.calibrate.external <- function(object, ...) {
  method <- setdiff(class(object), "hdnom.calibrate.external")

  x <- object

  attr(x, "model.type") <- NULL
  attr(x, "pred.at") <- NULL
  attr(x, "ngroup") <- NULL
  attr(x, "risk.group") <- NULL
  attr(x, "surv.time") <- NULL
  attr(x, "surv.event") <- NULL

  cat("  External Calibration Summary Table\n")
  class(x) <- "matrix"
  print(x)

  invisible(x)
}

#' Plot external calibration results
#'
#' Plot external calibration results
#'
#' @param x An object returned by \code{\link{calibrate_external}}.
#' @param xlim x axis limits of the plot.
#' @param ylim y axis limits of the plot.
#' @param col.pal Color palette to use. Possible values are
#' \code{"JCO"}, \code{"Lancet"}, \code{"NPG"}, and \code{"AAAS"}.
#' Default is \code{"JCO"}.
#' @param ... Other parameters for \code{plot}.
#'
#' @method plot hdnom.calibrate.external
#'
#' @export
#'
#' @importFrom ggplot2 ggplot aes_string geom_errorbar
#' geom_line geom_point geom_abline xlab ylab
#'
#' @examples
#' NULL
plot.hdnom.calibrate.external <- function(
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
