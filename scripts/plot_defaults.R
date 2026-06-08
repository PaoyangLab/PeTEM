PETEM_THEME_BASE_SIZE <- 16
PETEM_PLOT_TITLE_SIZE <- 14
PETEM_AXIS_TITLE_SIZE <- 14
PETEM_AXIS_TEXT_SIZE <- 12
PETEM_LEGEND_TITLE_SIZE <- 12
PETEM_LEGEND_TEXT_SIZE <- 12
PETEM_ANNOTATION_TEXT_SIZE <- 4

PETEM_BASE_CEX_AXIS <- 1.2
PETEM_BASE_CEX_LABEL <- 1.3
PETEM_BASE_CEX_NOTE <- 1.0
PETEM_BASE_LINE_WIDTH <- 4

petem_theme_bw <- function(base_size = PETEM_THEME_BASE_SIZE) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(size = PETEM_PLOT_TITLE_SIZE, face = "bold", hjust = 0.5),
      axis.title = element_text(size = PETEM_AXIS_TITLE_SIZE, face = "bold"),
      axis.text = element_text(size = PETEM_AXIS_TEXT_SIZE, face = "bold"),
      legend.title = element_text(size = PETEM_LEGEND_TITLE_SIZE, face = "bold"),
      legend.text = element_text(size = PETEM_LEGEND_TEXT_SIZE)
    )
}

petem_theme_classic <- function(base_size = PETEM_THEME_BASE_SIZE) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(size = PETEM_PLOT_TITLE_SIZE, face = "bold", hjust = 0.5),
      axis.title = element_text(size = PETEM_AXIS_TITLE_SIZE, face = "bold"),
      axis.text = element_text(size = PETEM_AXIS_TEXT_SIZE, face = "bold"),
      legend.title = element_text(size = PETEM_LEGEND_TITLE_SIZE, face = "bold"),
      legend.text = element_text(size = PETEM_LEGEND_TEXT_SIZE)
    )
}

petem_base_par <- function(extra = list()) {
  defaults <- list(
    cex.axis = PETEM_BASE_CEX_AXIS,
    cex.lab = PETEM_BASE_CEX_LABEL,
    cex.main = PETEM_BASE_CEX_LABEL,
    font.axis = 2,
    font.lab = 2
  )
  do.call(par, c(defaults, extra))
}
