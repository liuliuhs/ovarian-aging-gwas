#!/usr/bin/env Rscript

# Figure 4E: pQTL-GWAS colocalization posterior probabilities.
# Only genes belonging to the 737-gene candidate set are retained.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

if (.Platform$OS.type == "windows") {
  invisible(suppressWarnings(try(Sys.setlocale("LC_CTYPE", "Chinese"), silent = TRUE)))
}

args <- commandArgs(trailingOnly = TRUE)

input_file <- if (length(args) >= 1L) {
  args[[1]]
} else {
  "C:/Users/92975/Documents/F4_runtime/F4D_pQTL_reference_figure_audit.tsv"
}

output_prefix <- if (length(args) >= 2L) {
  args[[2]]
} else {
  "C:/Users/92975/Documents/F4_output/Figure4E_pQTL_colocalization"
}

if (!file.exists(input_file)) stop("Input file not found: ", input_file)

# fread may fail on Chinese paths under some Windows locales.
input_read_path <- input_file
if (.Platform$OS.type == "windows" && grepl("[^ -~]", input_file)) {
  input_read_path <- file.path(tempdir(), "Figure4E_pQTL_input.tsv")
  copied <- file.copy(input_file, input_read_path, overwrite = TRUE)
  if (!isTRUE(copied)) stop("Could not create a temporary ASCII-path input copy.")
}

coloc_data <- fread(input_read_path, sep = "\t", header = TRUE, data.table = TRUE)
if (!identical(input_read_path, input_file)) unlink(input_read_path)

required_columns <- c("gene_symbol", "figure_PP_display", "in_S22_737")
missing_columns <- setdiff(required_columns, names(coloc_data))
if (length(missing_columns) > 0L) {
  stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
}

# This filter removes SF3B4 because it is not present in the 737-gene S22 set.
plot_data <- coloc_data[
  in_S22_737 %in% c(TRUE, "TRUE", 1L) &
    !is.na(gene_symbol) & !is.na(figure_PP_display)
]
plot_data[, pp4_percent := as.numeric(figure_PP_display)]
plot_data <- plot_data[is.finite(pp4_percent)]
plot_data <- unique(plot_data, by = "gene_symbol")
setorder(plot_data, -pp4_percent, gene_symbol)

if (nrow(plot_data) == 0L) stop("No eligible colocalization results were found.")

plot_data[, support := ifelse(pp4_percent >= 80, "PP4 >= 80%", "PP4 < 80%")]
plot_data[, gene_symbol := factor(gene_symbol, levels = rev(as.character(gene_symbol)))]

format_pp4 <- function(x) {
  ifelse(
    x >= 1,
    sprintf("%.1f", x),
    ifelse(
      x >= 0.01,
      sub("0+$", "", sub("\\.$", "", sprintf("%.3f", x))),
      format(x, scientific = TRUE, digits = 3, trim = TRUE)
    )
  )
}

plot_data[, value_label := format_pp4(pp4_percent)]
plot_data[, label_x := ifelse(pp4_percent < 3, 4.0, pp4_percent + 1.5)]

ovary_colour <- "#B83F82"
background_colour <- "#BDBDBD"

p <- ggplot(plot_data, aes(y = gene_symbol, x = pp4_percent, fill = support)) +
  geom_vline(
    xintercept = 80,
    colour = "#666666",
    linetype = "dashed",
    linewidth = 0.35
  ) +
  geom_col(width = 0.64) +
  geom_text(
    aes(x = label_x, label = value_label),
    hjust = 0,
    colour = "#303030",
    size = 2.55,
    family = "Arial"
  ) +
  scale_fill_manual(
    values = c(
      "PP4 >= 80%" = ovary_colour,
      "PP4 < 80%" = background_colour
    )
  ) +
  scale_y_discrete(
    expand = expansion(add = 0.55)
  ) +
  scale_x_continuous(
    breaks = c(0, 25, 50, 75, 100),
    limits = c(0, 115),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(x = "Colocalization PP4 (%)", y = NULL) +
  theme_classic(base_size = 8, base_family = "Arial") +
  theme(
    legend.position = "none",
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(
      size = 8, face = "italic", family = "Arial", colour = "#252525"
    ),
    axis.text.x = element_text(size = 7, colour = "#252525"),
    axis.title.x = element_text(size = 8, margin = margin(t = 4)),
    plot.margin = margin(2, 2, 2, 2, unit = "mm")
  )

output_dir <- dirname(output_prefix)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

plot_width <- 55
plot_height <- 45

save_pdf_windows_safe <- function(plot, destination, width, height) {
  temporary_pdf <- tempfile(pattern = "Figure4E_", fileext = ".pdf")
  on.exit(unlink(temporary_pdf), add = TRUE)

  ggsave(
    temporary_pdf, plot,
    width = width, height = height, units = "mm",
    device = cairo_pdf
  )

  copied <- file.copy(temporary_pdf, destination, overwrite = TRUE)
  final_destination <- destination
  if (!isTRUE(copied)) {
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    final_destination <- sub(
      "\\.pdf$", paste0("_new_", timestamp, ".pdf"), destination,
      ignore.case = TRUE
    )
    copied <- file.copy(temporary_pdf, final_destination, overwrite = FALSE)
    if (!isTRUE(copied)) {
      stop("The PDF was created but could not be copied to the output folder.")
    }
    warning("The existing PDF appears locked. Saved as: ", final_destination)
  }
  final_destination
}

pdf_created <- save_pdf_windows_safe(
  p, paste0(output_prefix, ".pdf"), plot_width, plot_height
)

if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(
    paste0(output_prefix, ".svg"), p,
    width = plot_width, height = plot_height, units = "mm",
    device = svglite::svglite
  )
}

message("Retained genes: ", paste(as.character(rev(levels(plot_data$gene_symbol))), collapse = ", "))
message("Created: ", pdf_created)
