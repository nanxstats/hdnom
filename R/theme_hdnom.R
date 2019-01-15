#' Plot theme (ggplot2) for hdnom
#'
#' @param base_size base font size
#'
#' @export theme_hdnom
#'
#' @importFrom ggplot2 margin

theme_hdnom <- function(base_size = 14) {
  theme_bw(base_size = base_size) +
    theme(
      text = element_text(size = base_size),
      plot.title = element_text(face = "bold", size = base_size * 1.3, hjust = 0.5),
      axis.title = element_text(face = "plain"),
      axis.title.x = element_text(size = base_size * 1.25, vjust = -2),
      axis.title.y = element_text(size = base_size * 1.25, angle = 90, vjust = 4),
      axis.text = element_text(size = base_size, color = "#000000"),
      axis.ticks = element_line(),
      axis.ticks.length = unit(2, "mm"),
      axis.line = element_line(colour = "#000000"),
      legend.key = element_rect(colour = NA),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.key.size = unit(2, "mm"),
      legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
      legend.title = element_text(face = "plain", size = base_size),
      plot.background = element_rect(colour = NA),
      plot.margin = unit(c(10, 5, 5, 5), "mm"),
      panel.background = element_rect(colour = NA),
      panel.border = element_rect(colour = NA),
      panel.grid.major = element_line(colour = "#F0F0F0"),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
      strip.text = element_text(face = "bold", size = base_size)
    )
}
