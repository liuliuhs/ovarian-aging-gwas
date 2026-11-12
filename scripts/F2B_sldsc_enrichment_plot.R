# ============================================================
# Figure 3B: SNP-heritability enrichment estimated by S-LDSC
# ============================================================

required_packages <- c("ggplot2", "svglite")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing R package(s): ",
    paste(missing_packages, collapse = ", "),
    ". Install them before running this script."
  )
}

library(ggplot2)

# R may inherit an unavailable UTF-8 locale on this Windows installation.
# Switching to the installed Chinese locale allows graphics devices to write
# to the Chinese output path below.
if (.Platform$OS.type == "windows") {
  suppressWarnings(Sys.setlocale("LC_CTYPE", "Chinese"))
}

# Full paths are used so the script is independent of the working directory.
input_file <- "C:/Users/92975/Desktop/wyw/mvReproAge_baseline.results"
output_dir <- paste0(
  "C:/Users/92975/Desktop/wyw/",
  "\u6587\u7ae0/",
  "\u5375\u5de2gwas\u6587\u7ae0\u56fe\u7247/F2"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

# LDSC *.results files are tab-delimited.
raw <- read.delim(
  input_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_columns <- c(
  "Category",
  "Enrichment",
  "Enrichment_std_error",
  "Enrichment_p"
)

missing_columns <- setdiff(required_columns, names(raw))
if (length(missing_columns) > 0) {
  stop("Missing required column(s): ", paste(missing_columns, collapse = ", "))
}

raw$Enrichment <- as.numeric(raw$Enrichment)
raw$Enrichment_std_error <- as.numeric(raw$Enrichment_std_error)
raw$Enrichment_p <- as.numeric(raw$Enrichment_p)

# The baseline row is not a tested annotation. The input contains 52 tested
# annotations: 28 primary annotations and 24 corresponding +/-500-bp
# extensions. The figure displays the 28 primary annotations only, while the
# Bonferroni denominator remains all 52 tested annotations.
tested <- raw[
  raw$Category != "baseL2_0" & !is.na(raw$Enrichment_p),
  ,
  drop = FALSE
]

n_tested <- nrow(tested)
bonferroni_threshold <- 0.05 / n_tested

plot_data <- tested[
  !grepl("\\.extend\\.500", tested$Category),
  ,
  drop = FALSE
]

# Remove LDSC suffixes while retaining recognizable baseline-model names.
plot_data$Annotation <- sub("L2_0$", "", plot_data$Category)
plot_data$Annotation <- sub("_$", "", plot_data$Annotation)
plot_data$Significant <- plot_data$Enrichment_p < bonferroni_threshold

# Sort from the greatest to the smallest enrichment, as in the current panel.
plot_data <- plot_data[order(plot_data$Enrichment, decreasing = TRUE), ]
plot_data$Display_label <- paste0(
  plot_data$Annotation,
  ifelse(plot_data$Significant, "*", "")
)
plot_data$Display_label <- factor(
  plot_data$Display_label,
  levels = plot_data$Display_label
)

message("Tested annotations (Bonferroni denominator): ", n_tested)
message(
  "Bonferroni threshold: ",
  format(bonferroni_threshold, scientific = TRUE, digits = 5)
)
message("Primary annotations plotted: ", nrow(plot_data))
message("Significant annotations marked with an asterisk: ", sum(plot_data$Significant))

write.table(
  plot_data[, c(
    "Category",
    "Annotation",
    "Enrichment",
    "Enrichment_std_error",
    "Enrichment_p",
    "Significant"
  )],
  file = file.path(output_dir, "F2B_sldsc_enrichment_plot_data.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

figure_b <- ggplot(
  plot_data,
  aes(x = Display_label, y = Enrichment)
) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.55,
    colour = "#5F5F5F"
  ) +
  geom_col(
    width = 0.82,
    fill = "#79A9CC",
    colour = NA
  ) +
  geom_errorbar(
    aes(
      ymin = Enrichment - Enrichment_std_error,
      ymax = Enrichment + Enrichment_std_error
    ),
    width = 0.34,
    linewidth = 0.55,
    colour = "black"
  ) +
  scale_y_continuous(
    limits = c(-7, 20),
    breaks = seq(-5, 20, by = 5),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = NULL,
    y = "Heritability enrichment"
  ) +
  coord_cartesian(ylim = c(-7, 20), clip = "off") +
  theme_classic(base_size = 10, base_family = "Arial") +
  theme(
    axis.title.y = element_text(
      size = 9,
      face = "plain",
      margin = margin(r = 5)
    ),
    axis.text.y = element_text(size = 8, colour = "black"),
    axis.text.x = element_text(
      size = 6.7,
      angle = 52,
      hjust = 1,
      vjust = 1,
      colour = "black"
    ),
    axis.line = element_line(linewidth = 0.55, colour = "black"),
    axis.ticks = element_line(linewidth = 0.45, colour = "black"),
    axis.ticks.length = unit(1.5, "mm"),
    plot.margin = margin(t = 5, r = 5, b = 7, l = 16)
  )

pdf_file <- file.path(output_dir, "F2B_sldsc_enrichment_plot_rebuilt.pdf")
svg_file <- file.path(output_dir, "F2B_sldsc_enrichment_plot_rebuilt.svg")
png_file <- file.path(output_dir, "F2B_sldsc_enrichment_plot_rebuilt.png")

ggsave(
  pdf_file,
  figure_b,
  width = 3.94,
  height = 2.97,
  units = "in",
  device = cairo_pdf
)

ggsave(
  svg_file,
  figure_b,
  width = 3.94,
  height = 2.97,
  units = "in",
  device = svglite::svglite
)

ggsave(
  png_file,
  figure_b,
  width = 3.94,
  height = 2.97,
  units = "in",
  dpi = 600,
  bg = "white"
)

message("Saved PDF: ", pdf_file)
message("Saved SVG: ", svg_file)
message("Saved PNG: ", png_file)
