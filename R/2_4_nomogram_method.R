#' Print nomograms objects
#'
#' Print nomograms objects
#'
#' @param x An object returned by \code{\link{as_nomogram}}.
#' @param ... Other parameters.
#'
#' @method print hdnom.nomogram
#'
#' @export
#'
#' @examples
#' NULL
print.hdnom.nomogram <- function(x, ...) {
  print(as_nomogram_raw(
    fit = x$"nomogram", fun = x$"bhfun",
    fun.at = x$"fun.at", funlabel = x$"funlabel",
    lp = TRUE, vnames = "labels", ...
  ))
  invisible(x)
}

#' Plot nomogram objects
#'
#' Plot nomogram objects
#'
#' @param x An object returned by \code{\link{as_nomogram}}.
#' @param ... Other parameters.
#'
#' @method plot hdnom.nomogram
#'
#' @export
#'
#' @importFrom graphics plot
#'
#' @examples
#' NULL
plot.hdnom.nomogram <- function(x, ...) {
  plot(as_nomogram_raw(
    fit = x$"nomogram", fun = x$"bhfun",
    fun.at = x$"fun.at", funlabel = x$"funlabel",
    lp = TRUE, vnames = "labels", ...
  ))
  invisible(x)
}
