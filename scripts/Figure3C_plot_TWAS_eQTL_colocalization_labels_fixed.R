# ============================================================
# TWAS and ovary eQTL colocalization visualization
#
# Input:
#   Supplementary table.xlsx
#   - Table. S16: TWAS.Z, TWAS.P and COLOC.PP4
#   - Table. S22: gene symbols and final eQTLColoc indicator
#
# Output:
#   TWAS_eQTL_colocalization_labels_fixed.pdf
#   TWAS_eQTL_colocalization_labels_fixed.svg
#   TWAS_eQTL_colocalization_labels_fixed.png
# ============================================================

required_packages <- c("readxl", "dplyr", "ggplot2", "svglite")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Please install the required packages first: ",
    paste(missing_packages, collapse = ", "),
    "\nRun: install.packages(c(",
    paste(sprintf("\"%s\"", missing_packages), collapse = ", "),
    "))"
  )
}

library(readxl)
library(dplyr)
library(ggplot2)

# ---------- File locations ----------
input_file <- "C:/Users/92975/Desktop/wyw/文章/Supplementary table.xlsx"
output_dir <- "C:/Users/92975/Desktop/wyw/文章/卵巢gwas文章图片/F3"

if (!file.exists(input_file)) {
  stop("Input workbook was not found: ", input_file)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ---------- Read the two relevant supplementary tables ----------
# The first row of each worksheet is the table title, so skip = 1.
twas_raw <- read_excel(
  input_file,
  sheet = "Table. S16",
  skip = 1,
  na = c("NA", "   NA", "")
)

prioritization <- read_excel(
  input_file,
  sheet = "Table. S22",
  skip = 1,
  na = c("NA", "")
)

required_s16 <- c("ID", "TWAS.Z", "TWAS.P", "COLOC.PP4")
required_s22 <- c("ensID", "symbol", "eQTLColoc")

if (!all(required_s16 %in% names(twas_raw))) {
  stop("Table. S16 is missing: ",
       paste(setdiff(required_s16, names(twas_raw)), collapse = ", "))
}
if (!all(required_s22 %in% names(prioritization))) {
  stop("Table. S22 is missing: ",
       paste(setdiff(required_s22, names(prioritization)), collapse = ", "))
}

# ---------- Prepare and match genes ----------
# Remove Ensembl version suffixes, e.g. ENSG000001234.12 -> ENSG000001234.
gene_key <- prioritization %>%
  transmute(
    ensID_clean = sub("\\..*$", "", as.character(ensID)),
    symbol = as.character(symbol),
    coloc_supported = suppressWarnings(as.numeric(eQTLColoc)) == 1
  ) %>%
  distinct(ensID_clean, .keep_all = TRUE)

n_tested <- nrow(twas_raw)
twas_threshold <- 0.05 / n_tested

plot_data <- twas_raw %>%
  transmute(
    ID = as.character(ID),
    ensID_clean = sub("\\..*$", "", as.character(ID)),
    TWAS_Z = suppressWarnings(as.numeric(TWAS.Z)),
    TWAS_P = suppressWarnings(as.numeric(TWAS.P)),
    PP4 = suppressWarnings(as.numeric(COLOC.PP4))
  ) %>%
  filter(!is.na(TWAS_P), TWAS_P < twas_threshold) %>%
  left_join(gene_key, by = "ensID_clean") %>%
  mutate(
    coloc_supported = if_else(is.na(coloc_supported), FALSE, coloc_supported),
    symbol = if_else(is.na(symbol) | symbol == "", ensID_clean, symbol),
    evidence = if_else(
      coloc_supported,
      "Colocalization-supported",
      "Other significant TWAS genes"
    ),
    neglog10_p = -log10(pmax(TWAS_P, .Machine$double.xmin))
  ) %>%
  filter(!is.na(TWAS_Z), !is.na(PP4))

n_twas <- nrow(plot_data)
n_supported <- sum(plot_data$coloc_supported)

message("Bonferroni threshold: ", signif(twas_threshold, 4))
message("Significant TWAS genes plotted: ", n_twas)
message("Colocalization-supported genes plotted: ", n_supported)

if (n_twas != 128) {
  warning("Expected 128 significant TWAS genes, but found ", n_twas, ".")
}
if (n_supported != 64) {
  warning("Expected 64 colocalization-supported genes, but found ", n_supported, ".")
}

# The manuscript defines colocalization support as PP.H4 > 0.60.
pp4_threshold <- 0.60

# Rank all significant TWAS genes by PP.H4. This avoids the large empty space
# around TWAS Z = 0 that appears when only significant TWAS genes are plotted.
plot_data <- plot_data %>%
  mutate(
    evidence = factor(
      evidence,
      levels = c("Other significant TWAS genes", "Colocalization-supported")
    )
  ) %>%
  arrange(PP4, TWAS_P) %>%
  mutate(gene_rank = row_number())

# These five genes are highlighted in the Results text as the most significant
# TWAS genes among the colocalization-supported set. Fixed label coordinates
# are used so their names remain separated even when the figure is reduced.
genes_to_label <- c("CBY1", "RP3-508I15.9", "FAM227A", "CCDC200", "TMEM106A")
label_data <- plot_data %>%
  filter(symbol %in% genes_to_label) %>%
  mutate(
    label_x = case_when(
      symbol == "TMEM106A"      ~ 136,
      symbol == "CCDC200"       ~ 136,
      symbol == "FAM227A"       ~ 136,
      symbol == "CBY1"          ~ 136,
      symbol == "RP3-508I15.9"  ~ 136,
      TRUE ~ as.numeric(gene_rank)
    ),
    label_y = case_when(
      symbol == "CCDC200"       ~ 0.69,
      symbol == "TMEM106A"      ~ 0.765,
      symbol == "FAM227A"       ~ 0.84,
      symbol == "CBY1"          ~ 0.915,
      symbol == "RP3-508I15.9"  ~ 0.99,
      TRUE ~ PP4
    )
  )

group_labels <- c(
  "Other significant TWAS genes" = "Other significant TWAS genes",
  "Colocalization-supported" = "eQTL colocalization support"
)

p <- ggplot(plot_data, aes(x = gene_rank, y = PP4)) +
  geom_hline(
    yintercept = pp4_threshold,
    linewidth = 0.65,
    linetype = "dashed",
    colour = "#2A9D8F"
  ) +
  geom_point(
    aes(colour = evidence),
    size = 3.1,
    alpha = 0.92,
    stroke = 0.15
  ) +
  geom_segment(
    data = label_data,
    aes(x = gene_rank, y = PP4, xend = label_x - 2, yend = label_y),
    inherit.aes = FALSE,
    linewidth = 0.38,
    colour = "#A85A55"
  ) +
  geom_text(
    data = label_data,
    aes(x = label_x, y = label_y, label = symbol),
    inherit.aes = FALSE,
    size = 3.55,
    colour = "#B2182B",
    fontface = "italic",
    hjust = 0,
    vjust = 0.5
  ) +
  annotate(
    "text",
    x = 2,
    y = pp4_threshold + 0.025,
    hjust = 0,
    vjust = 0,
    label = "PP.H4 > 0.60",
    size = 4.0,
    fontface = "bold",
    colour = "#237A70"
  ) +
  scale_colour_manual(
    values = c(
      "Other significant TWAS genes" = "#BFC3C7",
      "Colocalization-supported" = "#D73027"
    ),
    labels = group_labels,
    name = NULL,
    drop = FALSE
  ) +
  scale_y_continuous(
    limits = c(0, 1.12),
    breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_x_continuous(
    breaks = c(1, 32, 64, 96, 128),
    limits = c(1, 172),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    #title = "eQTL colocalization",
    #subtitle = paste0(
    #   n_twas, " significant TWAS genes  →  ",
    #   n_supported, " genes with ", "PP.H4 > 0.60"
    # ),
    x = "Ranked significant TWAS genes",
    y = "Colocalization probability (PP.H4)"
  ) +
  guides(
    colour = guide_legend(
      override.aes = list(size = 5, alpha = 1),
      order = 1
    )
  ) +
  theme_classic(base_size = 13, base_family = "Arial") +
  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11.5, colour = "#4F4F4F", margin = margin(b = 8)),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 11, colour = "black"),
    legend.position = "top",
    legend.justification = "left",
    legend.box = "horizontal",
    legend.text = element_text(size = 10.5),
    plot.margin = margin(12, 18, 12, 12)
  )

# ---------- Export ----------
pdf_file <- file.path(output_dir, "F3C_TWAS_eQTL_colocalization_labels_fixed.pdf")
svg_file <- file.path(output_dir, "F3C_TWAS_eQTL_colocalization_labels_fixed.svg")

ggsave(pdf_file, p, width = 6, height = 4.5, units = "in", device = cairo_pdf)
ggsave(svg_file, p, width = 6, height = 4.5, units = "in", device = svglite::svglite)

message("Saved files:")
message(pdf_file)
message(svg_file)

print(p)
