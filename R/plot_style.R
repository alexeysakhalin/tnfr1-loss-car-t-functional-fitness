
library(ggplot2)

theme_pub <- function() {
  theme_classic(base_size = 16) +
    theme(
      axis.title = element_text(size = 19, face = "bold"),
      axis.text  = element_text(size = 16, face = "bold"),

      strip.text = element_text(size = 16, face = "bold"),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 1),

      plot.title = element_text(
        size = 18,
        face = "bold",
        hjust = 0.5   #
      ),

      legend.title = element_text(size = 14, face = "bold"),
      legend.text  = element_text(size = 13),

      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 1),

      plot.margin = margin(15, 15, 25, 15)
    )
}

save_ggplot_pair <- function(plot_object, png_path, width, height, dpi = 600) {
  dir.create(dirname(png_path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    filename = png_path,
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = "white"
  )
  tiff_path <- paste0(tools::file_path_sans_ext(png_path), ".tiff")
  ggplot2::ggsave(
    filename = tiff_path,
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    compression = "lzw",
    bg = "white"
  )
  invisible(c(png = png_path, tiff = tiff_path))
}
