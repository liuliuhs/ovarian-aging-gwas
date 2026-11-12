#!/usr/bin/env Rscript

# Figure S4B: TWAS-DEG directional concordance among the 203 prioritized genes.
# Direction is encoded by red and blue circles with a legend.
#
# Usage:
#   Rscript FigureS4B_TWAS_DEG_concordance.R [Supplementary_table.xlsx] [output_dir]

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(readxl)
})

if (.Platform$OS.type == "windows") {
  invisible(suppressWarnings(try(Sys.setlocale("LC_CTYPE", "Chinese"), silent = TRUE)))
}

args <- commandArgs(trailingOnly = TRUE)
workbook_file <- if (length(args) >= 1L) args[[1]] else
  "C:/Users/92975/Documents/F4B_runtime/Supplementary_table_readonly.xlsx"
output_dir <- if (length(args) >= 2L) args[[2]] else
  "C:/Users/92975/Desktop/wyw/文章/卵巢gwas文章图片/S4"

if (!file.exists(workbook_file)) stop("Workbook not found: ", workbook_file)

# readxl can fail on non-ASCII Windows paths, so read from a temporary copy.
workbook_read_path <- workbook_file
if (.Platform$OS.type == "windows" && grepl("[^ -~]", workbook_file)) {
  workbook_read_path <- file.path(tempdir(), "Supplementary_table_S4B.xlsx")
  if (!file.copy(workbook_file, workbook_read_path, overwrite = TRUE)) {
    stop("Could not create a temporary ASCII-path workbook copy.")
  }
}

s22 <- as.data.table(read_excel(workbook_read_path, sheet = "Table. S22", skip = 1))
if (!identical(workbook_read_path, workbook_file)) unlink(workbook_read_path)

required_columns <- c(
  "symbol", "sum", "twas", "deg", "TWAS.direction", "deg.direction"
)
missing_columns <- setdiff(required_columns, names(s22))
if (length(missing_columns) > 0L) {
  stop("Table S22 is missing required columns: ", paste(missing_columns, collapse = ", "))
}

to_flag <- function(x) suppressWarnings(as.numeric(x)) == 1
normalise_direction <- function(x) {
  x <- tolower(trimws(as.character(x)))
  fifelse(
    x %chin% c("up", "positive", "+"), "Up",
    fifelse(x %chin% c("down", "negative", "-"), "Down", NA_character_)
  )
}

s22[, evidence_sum := suppressWarnings(as.numeric(sum))]
hp <- s22[
  is.finite(evidence_sum) & evidence_sum >= 3 &
    !is.na(symbol) & trimws(as.character(symbol)) != ""
]
hp[, symbol_key := toupper(trimws(as.character(symbol)))]
hp <- unique(hp, by = "symbol_key")
if (nrow(hp) != 203L) {
  warning("Expected 203 genes with at least three evidence layers; found ", nrow(hp), ".")
}

shared <- hp[to_flag(twas) & to_flag(deg), .(
  symbol = as.character(symbol),
  evidence_sum,
  TWAS = normalise_direction(TWAS.direction),
  DEG = normalise_direction(deg.direction)
)]
shared <- shared[!is.na(TWAS) & !is.na(DEG)]
shared[, Concordance := fifelse(TWAS == DEG, "Concordant", "Discordant")]
setorder(shared, Concordance, -evidence_sum, symbol)
shared[, y := .N - .I + 1L]

direction_long <- melt(
  shared,
  id.vars = c("symbol", "evidence_sum", "Concordance", "y"),
  measure.vars = c("TWAS", "DEG"),
  variable.name = "Analysis",
  value.name = "Direction"
)
direction_long[, x := fifelse(Analysis == "TWAS", 1.00, 1.75)]

cell_background <- unique(direction_long[, .(x, y)])

p <- ggplot() +
  geom_rect(
    data = cell_background,
    aes(xmin = x - 0.30, xmax = x + 0.30, ymin = y - 0.30, ymax = y + 0.30),
    fill = "#F7F7F7", colour = "#D8D8D8", linewidth = 0.32
  ) +
  geom_point(
    data = direction_long,
    aes(x = x, y = y, fill = Direction),
    shape = 21, size = 4.15, stroke = 0.50, colour = "white"
  ) +
  scale_fill_manual(
    values = c("Up" = "#D83A3A", "Down" = "#3F6FB5"),
    breaks = c("Up", "Down"),
    name = "Direction"
  ) +
  scale_x_continuous(
    position = "top",
    breaks = c(1.00, 1.75),
    labels = c("TWAS", "DEG"),
    limits = c(0.55, 2.20),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    breaks = shared$y,
    labels = shared$symbol,
    limits = c(0.45, nrow(shared) + 0.55),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(x = NULL, y = NULL) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_family = "Arial", base_size = 9) +
  theme(
    panel.grid = element_blank(),
    axis.text.x.top = element_text(
      size = 10.0, face = "bold", colour = "#252525", margin = margin(b = 3)
    ),
    axis.text.y = element_text(
      size = 9.2, face = "italic", colour = "#252525", margin = margin(r = 4)
    ),
    axis.ticks = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 8.2, face = "bold"),
    legend.text = element_text(size = 8.2),
    legend.key.size = grid::unit(4.0, "mm"),
    legend.key.width = grid::unit(5.0, "mm"),
    legend.spacing.x = grid::unit(0.5, "mm"),
    legend.margin = margin(t = 0.5, unit = "mm"),
    legend.background = element_rect(fill = "transparent", colour = NA),
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA),
    plot.margin = margin(1, 2, 1, 2, unit = "mm")
  )

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_prefix <- file.path(output_dir, "FigureS4B_TWAS_DEG_concordance")

save_and_copy <- function(extension, device, ...) {
  temporary_file <- tempfile(pattern = "FigureS4B_", fileext = paste0(".", extension))
  ggsave(
    temporary_file, p,
    width = 75, height = 85, units = "mm",
    device = device, bg = "transparent", ...
  )
  destination <- paste0(output_prefix, ".", extension)
  if (!file.copy(temporary_file, destination, overwrite = TRUE)) {
    unlink(temporary_file)
    stop("Could not write output file: ", destination)
  }
  unlink(temporary_file)
  destination
}

created_pdf <- save_and_copy("pdf", cairo_pdf)
created_svg <- save_and_copy("svg", svglite::svglite)
created_png <- save_and_copy("png", ragg::agg_png, res = 600)

message("Shared genes: ", nrow(shared))
message("Concordant: ", shared[Concordance == "Concordant", .N])
message("Discordant: ", shared[Concordance == "Discordant", .N])
message("Created: ", created_pdf)
message("Created: ", created_svg)
message("Created: ", created_png)
