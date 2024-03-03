#' Print model comparison by validation results
#'
#' Print model comparison by validation results
#'
#' @param x An object returned by \code{\link{compare_by_validate}}.
#' @param ... Other parameters (not used).
#'
#' @method print hdnom.compare.validate
#'
#' @export
#'
#' @examples
#' NULL
print.hdnom.compare.validate <- function(x, ...) {
  for (i in seq_along(x)) {
    print(x[[i]])
    cat("\n\n")
  }
  invisible(x)
}

#' Summary of model comparison by validation results
#'
#' Summary of model comparison by validation results
#'
#' @param object An object \code{\link{compare_by_validate}}.
#' @param silent Print summary table header or not,
#' default is \code{FALSE}.
#' @param ... Other parameters (not used).
#'
#' @method summary hdnom.compare.validate
#'
#' @export
#'
#' @examples
#' NULL
summary.hdnom.compare.validate <- function(object, silent = FALSE, ...) {
  for (i in seq_along(object)) {
    cat("Model type:", names(object)[i], "\n")
    print(summary(object[[i]], silent = TRUE))
    cat("\n")
  }
  invisible(object)
}

#' Plot model comparison by validation results
#'
#' Plot model comparison by validation results
#'
#' @param x An object returned by \code{\link{compare_by_validate}}.
#' @param interval Show maximum, minimum, 0.25, and 0.75 quantiles of
#' time-dependent AUC as ribbons? Default is \code{FALSE}.
#' @param col.pal Color palette to use. Possible values are
#' \code{"JCO"}, \code{"Lancet"}, \code{"NPG"}, and \code{"AAAS"}.
#' Default is \code{"JCO"}.
#' @param ylim Range of y coordinates. For example, \code{c(0.5, 1)}.
#' @param ... Other parameters (not used).
#'
#' @method plot hdnom.compare.validate
#'
#' @importFrom ggplot2 ggplot aes_string geom_point geom_line
#' scale_x_continuous scale_colour_manual ylab coord_cartesian
#'
#' @export
#'
#' @examples
#' NULL
plot.hdnom.compare.validate <- function(
    x, interval = FALSE,
    col.pal = c("JCO", "Lancet", "NPG", "AAAS"),
    ylim = NULL, ...) {
  n <- length(x)
  dflist <- vector("list", n)

  for (i in 1L:n) {
    dflist[[i]] <- as.data.frame(t(summary(x[[i]], silent = TRUE)))
    tauc_time <- attr(x[[i]], "tauc.time")

    # special processing for repeated cv
    if (any(grepl(pattern = "validate.repeated.cv", class(x[[i]])))) {
      names(dflist[[i]]) <- sapply(strsplit(names(dflist[[i]]), "Mean of "), "[", 2L)
    }

    dflist[[i]][, "Time"] <- tauc_time
    dflist[[i]][, "Model"] <- names(x)[i]
    names(dflist[[i]])[which(names(dflist[[i]]) == "0.25 Qt.")] <- "Qt25"
    names(dflist[[i]])[which(names(dflist[[i]]) == "0.75 Qt.")] <- "Qt75"
  }

  df <- Reduce("rbind", dflist)

  col.pal <- match.arg(col.pal)
  col_pal <- switch(col.pal,
    JCO = palette_jco(),
    Lancet = palette_lancet(),
    NPG = palette_npg(),
    AAAS = palette_aaas()
  )

  if (!interval) {
    ggplot(
      data = df,
      aes_string(x = "Time", y = "Mean", colour = "Model", fill = "Model")
    ) +
      geom_point() +
      geom_line() +
      geom_point(
        data = df,
        aes_string(x = "Time", y = "Median", colour = "Model", fill = "Model")
      ) +
      geom_line(
        data = df,
        aes_string(x = "Time", y = "Median", colour = "Model"),
        linetype = "dashed"
      ) +
      scale_x_continuous(breaks = df$"Time") +
      scale_colour_manual(values = col_pal) +
      coord_cartesian(ylim = ylim) +
      theme_hdnom() +
      ylab("Area under ROC")
  } else {
    ggplot(
      data = df,
      aes_string(x = "Time", y = "Mean", colour = "Model", fill = "Model")
    ) +
      geom_point() +
      geom_line() +
      geom_point(
        data = df,
        aes_string(x = "Time", y = "Median", colour = "Model", fill = "Model")
      ) +
      geom_line(
        data = df,
        aes_string(x = "Time", y = "Median", colour = "Model"),
        linetype = "dashed"
      ) +
      geom_ribbon(
        data = df,
        aes_string(ymin = "Qt25", ymax = "Qt75", colour = "Model", fill = "Model"),
        linetype = 0, alpha = 0.1
      ) +
      geom_ribbon(
        data = df,
        aes_string(ymin = "Min", ymax = "Max", colour = "Model", fill = "Model"),
        linetype = 0, alpha = 0.05
      ) +
      scale_x_continuous(breaks = df$"Time") +
      scale_colour_manual(values = col_pal) +
      scale_fill_manual(values = col_pal) +
      coord_cartesian(ylim = ylim) +
      theme_hdnom() +
      ylab("Area under ROC")
  }
}
