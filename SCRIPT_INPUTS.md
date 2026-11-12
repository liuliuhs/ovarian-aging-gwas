# Script input and dependency inventory

Automatically extracted references from the supplied scripts. Paths and worksheet names are retained from the originals; they are not bundled data or validated final-manuscript mappings.

## F2B_sldsc_enrichment_plot.R

```r
required_packages <- c("ggplot2", "svglite")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
library(ggplot2)
input_file <- "C:/Users/92975/Desktop/wyw/mvReproAge_baseline.results"
output_dir <- paste0(
# LDSC *.results files are tab-delimited.
  file = file.path(output_dir, "F2B_sldsc_enrichment_plot_data.tsv"),
```

## Figure3D_sixMethods_upsetR_updated.R

```r
  library(readxl)
  library(UpSetR)
input_file <- find_one_file(
output_dir <- dirname(figure3_reference)
output_prefix <- file.path(output_dir,"F3D_sixMethods_upsetR_updated")
  fileext=".xlsx"
  sheet="Table. S22",
```

## Figure3C_plot_TWAS_eQTL_colocalization_labels_fixed.R

```r
#   Supplementary table.xlsx
required_packages <- c("readxl", "dplyr", "ggplot2", "svglite")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
library(readxl)
library(dplyr)
library(ggplot2)
input_file <- "C:/Users/92975/Desktop/wyw/文章/Supplementary table.xlsx"
output_dir <- "C:/Users/92975/Desktop/wyw/文章/卵巢gwas文章图片/F3"
  sheet = "Table. S16",
  sheet = "Table. S22",
```

## Figure4D_MR_forest.R

```r
  library(data.table)
  library(ggplot2)
  "C:/Users/92975/Documents/F4_runtime/cisPQTL_OvarianAging.tsv"
  "C:/Users/92975/Documents/F4_runtime/cisPQTL_OvarianAging_737genes.tsv"
output_dir <- if (length(args) >= 3L) args[[3]] else
all_results <- read_tsv_windows_safe(all_pqtl_file, "F4D_all_pQTL.tsv")
candidate_results <- read_tsv_windows_safe(candidate_pqtl_file, "F4D_737_pQTL.tsv")
```

## Figure4E_pQTL_colocalization.R

```r
  library(data.table)
  library(ggplot2)
input_file <- if (length(args) >= 1L) {
  "C:/Users/92975/Documents/F4_runtime/F4D_pQTL_reference_figure_audit.tsv"
output_prefix <- if (length(args) >= 2L) {
  input_read_path <- file.path(tempdir(), "Figure4E_pQTL_input.tsv")
output_dir <- dirname(output_prefix)
```

## FigureS4B_TWAS_DEG_concordance.R

```r
#   Rscript FigureS4B_TWAS_DEG_concordance.R [Supplementary_table.xlsx] [output_dir]
  library(data.table)
  library(ggplot2)
  library(readxl)
workbook_file <- if (length(args) >= 1L) args[[1]] else
  "C:/Users/92975/Documents/F4B_runtime/Supplementary_table_readonly.xlsx"
output_dir <- if (length(args) >= 2L) args[[2]] else
  workbook_read_path <- file.path(tempdir(), "Supplementary_table_S4B.xlsx")
s22 <- as.data.table(read_excel(workbook_read_path, sheet = "Table. S22", skip = 1))
output_prefix <- file.path(output_dir, "FigureS4B_TWAS_DEG_concordance")
```

## method5_genes.R

```r
wyw <- read.delim("C:/Users/92975/Desktop/wyw/genePrioritization.tsv")
library(readr)
library(dplyr)
```

## vep_missense_variant1.R

```r
library(ggplot2)
library(readxl)
library(dplyr)
library(ggtext)   # Support Markdown and rich-text legend labels.
raw <- read_excel("文章/S5.xlsx", skip = 1)
library(cowplot)
library(grid)
```
