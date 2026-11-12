#!/usr/bin/env Rscript

# Figure 4D: two MR forest-plot versions.
# 1) All eight Bonferroni-significant proteins from the full cis-pQTL scan;
#    the non-overlapping SF3B4 row is shown in light gray.
# 2) Seven significant proteins overlapping the prespecified 737-gene set.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

if (.Platform$OS.type == "windows") {
  invisible(suppressWarnings(try(Sys.setlocale("LC_CTYPE", "Chinese"), silent = TRUE)))
}

args <- commandArgs(trailingOnly = TRUE)

all_pqtl_file <- if (length(args) >= 1L) args[[1]] else
  "C:/Users/92975/Documents/F4_runtime/cisPQTL_OvarianAging.tsv"
candidate_pqtl_file <- if (length(args) >= 2L) args[[2]] else
  "C:/Users/92975/Documents/F4_runtime/cisPQTL_OvarianAging_737genes.tsv"
output_dir <- if (length(args) >= 3L) args[[3]] else
  "C:/Users/92975/Documents/F4_output"

if (!file.exists(all_pqtl_file)) stop("Full cis-pQTL file not found: ", all_pqtl_file)
if (!file.exists(candidate_pqtl_file)) stop("737-gene cis-pQTL file not found: ", candidate_pqtl_file)

read_tsv_windows_safe <- function(path, temporary_name) {
  read_path <- path
  if (.Platform$OS.type == "windows" && grepl("[^ -~]", path)) {
    read_path <- file.path(tempdir(), temporary_name)
    copied <- file.copy(path, read_path, overwrite = TRUE)
    if (!isTRUE(copied)) stop("Could not copy input to an ASCII temporary path: ", path)
  }
  result <- fread(read_path, sep = "\t", header = TRUE, data.table = TRUE)
  if (!identical(read_path, path)) unlink(read_path)
  result
}

all_results <- read_tsv_windows_safe(all_pqtl_file, "F4D_all_pQTL.tsv")
candidate_results <- read_tsv_windows_safe(candidate_pqtl_file, "F4D_737_pQTL.tsv")

required_all <- c(
  "exposure", "ivw_or_wald.b", "ivw_or_wald.se",
  "ivw_or_wald.p.bonf", "pass_alltests"
)
missing_all <- setdiff(required_all, names(all_results))
if (length(missing_all) > 0L) {
  stop("Missing columns in full cis-pQTL file: ", paste(missing_all, collapse = ", "))
}
if (!"gene_symbol" %in% names(candidate_results)) {
  stop("The 737-gene cis-pQTL file does not contain gene_symbol.")
}

extract_gene_symbol <- function(exposure) {
  pieces <- strsplit(as.character(exposure), "_", fixed = TRUE)
  vapply(pieces, function(x) if (length(x) >= 3L) x[[3L]] else NA_character_, character(1))
}
as_logical_flag <- function(x) {
  toupper(trimws(as.character(x))) %in% c("TRUE", "T", "1")
}

all_results[, gene_symbol := extract_gene_symbol(exposure)]
all_results[, beta := as.numeric(ivw_or_wald.b)]
all_results[, se := as.numeric(ivw_or_wald.se)]
all_results[, bonferroni_p := as.numeric(ivw_or_wald.p.bonf)]
all_results[, pass_flag := as_logical_flag(pass_alltests)]

# Eight proteins significant after correction across the complete pQTL scan.
significant_all8 <- all_results[
  is.finite(beta) & is.finite(se) & is.finite(bonferroni_p) &
    bonferroni_p < 0.05 & pass_flag == TRUE &
    !is.na(gene_symbol) & gene_symbol != ""
]
significant_all8 <- unique(significant_all8, by = "gene_symbol")

candidate_symbols <- unique(toupper(trimws(as.character(candidate_results$gene_symbol))))
significant_all8[, in_737 := toupper(gene_symbol) %chin% candidate_symbols]

if (nrow(significant_all8) != 8L) {
  warning("Expected 8 full-scan significant genes, but found ", nrow(significant_all8), ".")
}
significant_overlap7 <- significant_all8[in_737 == TRUE]
if (nrow(significant_overlap7) != 7L) {
  warning("Expected 7 genes overlapping the 737-gene set, but found ", nrow(significant_overlap7), ".")
}

display_order <- c(
  "VAT1", "TMEM106A", "GDI2", "NRBP1",
  "SF3B4", "ITPKA", "GCKR", "PNKP"
)
positive_colour <- "#D83A3A"
negative_colour <- "#2F6FB2"
normal_label_colour <- "#252525"
nonoverlap_colour <- "#A6A6A6"

make_mr_plot <- function(plot_data, gray_nonoverlap = FALSE) {
  plot_data <- copy(plot_data)
  gene_order <- display_order[display_order %chin% plot_data$gene_symbol]
  gene_order <- c(gene_order, sort(setdiff(plot_data$gene_symbol, gene_order)))

  plot_data[, gene_factor := factor(gene_symbol, levels = rev(gene_order))]
  plot_data[, `:=`(
    ci_lower = beta - 1.96 * se,
    ci_upper = beta + 1.96 * se,
    effect_colour = ifelse(
      gray_nonoverlap & in_737 == FALSE,
      nonoverlap_colour,
      ifelse(beta >= 0, positive_colour, negative_colour)
    ),
    label_colour = ifelse(
      gray_nonoverlap & in_737 == FALSE,
      nonoverlap_colour,
      normal_label_colour
    )
  )]

  ggplot(plot_data, aes(y = gene_factor, x = beta)) +
    geom_vline(
      xintercept = c(-0.50, -0.25, 0.25, 0.50),
      colour = "#9A9A9A", linetype = "dashed", linewidth = 0.3
    ) +
    geom_vline(xintercept = 0, colour = "#303030", linewidth = 0.45) +
    geom_errorbar(
      aes(xmin = ci_lower, xmax = ci_upper, colour = effect_colour),
      width = 0, linewidth = 0.55, orientation = "y"
    ) +
    geom_point(aes(colour = effect_colour), size = 2.2) +
    geom_text(
      aes(x = -0.54, label = gene_symbol, colour = label_colour),
      hjust = 1, size = 2.75, family = "Arial", fontface = "italic"
    ) +
    scale_colour_identity() +
    scale_y_discrete(breaks = NULL, expand = expansion(add = 0.55)) +
    scale_x_continuous(
      breaks = c(-0.50, -0.25, 0, 0.25, 0.50),
      limits = c(-0.55, 0.55), expand = expansion(mult = c(0, 0))
    ) +
    labs(x = "IVW effect (95% CI)", y = NULL) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 8, base_family = "Arial") +
    theme(
      axis.line.y = element_blank(), axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 7, colour = "#252525"),
      axis.title.x = element_text(size = 8, margin = margin(t = 4)),
      plot.margin = margin(2, 2, 2, 15, unit = "mm")
    )
}

plot_all8 <- make_mr_plot(significant_all8, gray_nonoverlap = TRUE)
plot_overlap7 <- make_mr_plot(significant_overlap7, gray_nonoverlap = FALSE)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
all8_pdf <- file.path(output_dir, "Figure4D_MR_all8_SF3B4_gray.pdf")
overlap7_pdf <- file.path(output_dir, "Figure4D_MR_737overlap7.pdf")

# Cairo can fail when it writes directly to a Chinese Windows path. Save to an
# ASCII temporary path first, then copy the completed PDF to its destination.
# If the destination PDF is open in a viewer and therefore locked, keep the
# requested file untouched and write a timestamped alternative instead.
save_pdf_windows_safe <- function(plot, destination, width = 55, height = 45) {
  temporary_pdf <- tempfile(pattern = "Figure4D_", fileext = ".pdf")
  on.exit(unlink(temporary_pdf), add = TRUE)

  ggsave(
    temporary_pdf,
    plot,
    width = width,
    height = height,
    units = "mm",
    device = cairo_pdf
  )

  copied <- file.copy(temporary_pdf, destination, overwrite = TRUE)
  final_destination <- destination

  if (!isTRUE(copied)) {
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    final_destination <- sub(
      "\\.pdf$",
      paste0("_new_", timestamp, ".pdf"),
      destination,
      ignore.case = TRUE
    )
    copied <- file.copy(temporary_pdf, final_destination, overwrite = FALSE)
    if (!isTRUE(copied)) {
      stop(
        "The PDF was created, but it could not be copied to the output folder. ",
        "Close any open PDF viewer and check folder permissions."
      )
    }
    warning(
      "The existing PDF appears to be open or locked. Saved as: ",
      final_destination
    )
  }

  final_destination
}

all8_pdf_created <- save_pdf_windows_safe(plot_all8, all8_pdf)
overlap7_pdf_created <- save_pdf_windows_safe(plot_overlap7, overlap7_pdf)

message("Full-scan significant genes (n = ", nrow(significant_all8), "): ",
        paste(significant_all8$gene_symbol, collapse = ", "))
message("Overlap with 737-gene set (n = ", nrow(significant_overlap7), "): ",
        paste(significant_overlap7$gene_symbol, collapse = ", "))
message("Created: ", all8_pdf_created)
message("Created: ", overlap7_pdf_created)
