#' @noRd
my.sort <- function(x, index.return = FALSE, decreasing = FALSE) {
  if (decreasing) {
    x <- -x
  }
  y <- sort(x, method = "quick", index.return = index.return)
  if (decreasing) {
    y$x <- -y$x
  }
  y
}
