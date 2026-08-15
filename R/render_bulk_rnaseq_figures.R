#!/usr/bin/env Rscript

# Render the validated bulk RNA-seq manuscript panels in one R process and
# record the package namespaces that were present in that same process.

figure_scripts <- c(
  file.path("R", "03_figure_1_B_C.R"),
  file.path("R", "04_figure_2_B_suppl_S2D.R")
)

missing_scripts <- figure_scripts[!file.exists(figure_scripts)]
if (length(missing_scripts)) {
  stop(
    "Bulk RNA-seq figure script(s) not found: ",
    paste(missing_scripts, collapse = ", "),
    ". Run this entry point from the repository root."
  )
}

for (script in figure_scripts) {
  message("Rendering ", script)
  source(script, echo = FALSE, chdir = FALSE, encoding = "UTF-8")
}

expected_render_namespaces <- c("ggplot2", "ggrepel", "scales", "VennDiagram")
missing_namespaces <- setdiff(expected_render_namespaces, loadedNamespaces())
if (length(missing_namespaces)) {
  stop(
    "Expected rendering namespace(s) were not loaded: ",
    paste(missing_namespaces, collapse = ", "),
    "."
  )
}

release_dir <- file.path("results", "release")
dir.create(release_dir, recursive = TRUE, showWarnings = FALSE)
session_path <- file.path(release_dir, "sessionInfo_bulk_figures.txt")
session_lines <- c(
  "Bulk RNA-seq manuscript figure rendering session",
  paste("Rendered scripts:", paste(figure_scripts, collapse = ", ")),
  "",
  capture.output(sessionInfo())
)
writeLines(session_lines, session_path, useBytes = TRUE)
message("Recorded rendering-session provenance at ", session_path)
