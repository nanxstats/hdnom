#' Print model comparison by calibration results
#'
#' Print model comparison by calibration results
#'
#' @param x An object returned by \code{\link{compare_by_calibrate}}.
#' @param ... Other parameters (not used).
#'
#' @method print hdnom.compare.calibrate
#'
#' @export
#'
#' @examples
#' NULL
print.hdnom.compare.calibrate <- function(x, ...) {
  for (i in seq_along(x)) {
    print(x[[i]])
    cat("\n")
  }
  invisible(x)
}

#' Summary of model comparison by calibration results
#'
#' Summary of model comparison by calibration results
#'
#' @param object An object returned by \code{\link{compare_by_calibrate}}.
#' @param ... Other parameters (not used).
#'
#' @method summary hdnom.compare.calibrate
#'
#' @export
#'
#' @examples
#' NULL
summary.hdnom.compare.calibrate <- function(object, ...) {
  for (i in seq_along(object)) {
    cat("  Model type:", names(object)[i], "\n")
    summary(object[[i]])
    cat("\n")
  }
  invisible(object)
}

#' Plot model comparison by calibration results
#'
#' Plot model comparison by calibration results
#'
#' @param x An object returned by \code{\link{compare_by_calibrate}}.
#' @param xlim x axis limits of the plot.
#' @param ylim y axis limits of the plot.
#' @param col.pal Color palette to use. Possible values are
#' \code{"JCO"}, \code{"Lancet"}, \code{"NPG"}, and \code{"AAAS"}.
#' Default is \code{"JCO"}.
#' @param ... Other parameters (not used).
#'
#' @method plot hdnom.compare.calibrate
#'
#' @export
#'
#' @importFrom ggplot2 ggplot aes_string geom_errorbar
#' geom_line geom_point geom_abline scale_colour_manual
#' xlab ylab
#'
#' @examples
#' NULL
plot.hdnom.compare.calibrate <- function(
    x, xlim = c(0, 1), ylim = c(0, 1),
    col.pal = c("JCO", "Lancet", "NPG", "AAAS"), ...) {
  n <- length(x)
  dflist <- vector("list", n)

  for (i in 1L:n) {
    dflist[[i]] <- data.frame(
      "pre" = x[[i]][, "Predicted"],
      "obs" = x[[i]][, "Observed"],
      "ll" = x[[i]][, "Lower 95%"],
      "ul" = x[[i]][, "Upper 95%"]
    )
    dflist[[i]][, "Model"] <- names(x)[i]
  }

  df <- Reduce("rbind", dflist)

  col.pal <- match.arg(col.pal)
  col_pal <- switch(col.pal,
    JCO = palette_jco(),
    Lancet = palette_lancet(),
    NPG = palette_npg(),
    AAAS = palette_aaas()
  )

  ggplot(
    df,
    aes_string(
      x = "pre", y = "obs",
      xmin = xlim[1L], xmax = xlim[2L],
      ymin = ylim[1L], ymax = ylim[2L],
      colour = "Model"
    )
  ) +
    geom_errorbar(aes_string(ymin = "ll", ymax = "ul"), width = 0.02) +
    geom_line() +
    geom_point(size = 2) +
    geom_abline(slope = 1, intercept = 0, colour = "grey") +
    xlab("Predicted Survival Probability") +
    ylab("Observed Survival Probability") +
    scale_colour_manual(values = col_pal) +
    theme_hdnom()
}
