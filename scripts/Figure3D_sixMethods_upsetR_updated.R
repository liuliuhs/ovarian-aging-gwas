#!/usr/bin/env Rscript

# Figure 3D: overlap of genes identified by six prioritization methods
# Data source: Table. S22 in the supplementary workbook

start_time <- Sys.time()

suppressPackageStartupMessages({
  library(readxl)
  library(UpSetR)
})

# ============================================================
# File paths
# ============================================================

# Use an ASCII-only parent path and discover files by their English names.
# This avoids Chinese-path encoding problems in older Windows R locales.
search_root <- "C:/Users/92975/Desktop/wyw"

find_one_file <- function(pattern,file_description) {
  matches <- list.files(
    path=search_root,
    pattern=pattern,
    recursive=TRUE,
    full.names=TRUE,
    ignore.case=TRUE
  )

  matches <- matches[file.exists(matches)]

  if (length(matches) == 0) {
    stop("Could not find ",file_description," under: ",search_root)
  }

  if (length(matches) > 1) {
    stop(
      "More than one ",file_description," was found:\n",
      paste(matches,collapse="\n")
    )
  }

  matches[1]
}

input_file <- find_one_file(
  "^Supplementary_table_lhs_MGI_S20_S21_rebuilt[.]xlsx$",
  "supplementary workbook"
)

figure3_reference <- find_one_file(
  "^Figure3_for_review2[.]pdf$",
  "Figure 3 reference PDF"
)

output_dir <- dirname(figure3_reference)
output_prefix <- file.path(output_dir,"F3D_sixMethods_upsetR_updated")

dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)

message("Input workbook: ",input_file)
message("Output directory: ",output_dir)

# ============================================================
# Read and validate Table. S22
# ============================================================

# readxl/unzip can fail on non-ASCII Windows paths in some R locales.
# Copy the workbook to an ASCII-only temporary path before reading it.
temporary_workbook <- tempfile(
  pattern="Supplementary_table_",
  fileext=".xlsx"
)

if (!file.copy(input_file,temporary_workbook,overwrite=TRUE)) {
  stop("Failed to copy the workbook to a temporary ASCII path.")
}

dat <- read_excel(
  path=temporary_workbook,
  sheet="Table. S22",
  skip=1
)

unlink(temporary_workbook)

required_columns <- c(
  "ensID","symbol","missense","susie","eQTLColoc",
  "finemap","twas","magma"
)

missing_columns <- setdiff(required_columns,colnames(dat))

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns in Table. S22: ",
    paste(missing_columns,collapse=", ")
  )
}

dat <- dat[!is.na(dat$ensID) & nzchar(trimws(as.character(dat$ensID))),]

if (anyDuplicated(dat$ensID) > 0) {
  duplicated_ids <- unique(dat$ensID[duplicated(dat$ensID)])
  stop(
    "Duplicated ensID values were found: ",
    paste(head(duplicated_ids,10),collapse=", ")
  )
}

# Rename the six method columns to the labels displayed in the figure.
upset_data <- data.frame(
  VEP=as.integer(dat$missense),
  SuSiE=as.integer(dat$susie),
  coloc=as.integer(dat$eQTLColoc),
  FINEMAP=as.integer(dat$finemap),
  TWAS=as.integer(dat$twas),
  MAGMA=as.integer(dat$magma),
  check.names=FALSE
)

if (anyNA(upset_data)) {
  stop("Missing values were found in the six method columns.")
}

if (!all(unlist(upset_data,use.names=FALSE) %in% c(0L,1L))) {
  stop("The six method columns must contain only 0 and 1.")
}

# UpSetR displays the supplied set order from bottom to top. Supplying the
# reverse order gives VEP at the top and MAGMA at the bottom, as in Figure 3D.
method_order <- c("MAGMA","TWAS","FINEMAP","coloc","SuSiE","VEP")
set_sizes <- colSums(upset_data[,method_order,drop=FALSE])

message("Number of genes in Table. S22: ",nrow(upset_data))
message(
  "Set sizes: ",
  paste(paste0(names(set_sizes),"=",set_sizes),collapse=", ")
)

# ============================================================
# UpSet plot
# ============================================================

upset_plot <- upset(
  upset_data,
  sets=method_order,
  keep.order=TRUE,
  nsets=length(method_order),
  nintersects=NA,
  order.by=c("degree","freq"),
  decreasing=c(FALSE,TRUE),
  mb.ratio=c(0.70,0.30),
  main.bar.color="#3D3D3D",
  sets.bar.color="#3D3D3D",
  matrix.color="#3D3D3D",
  shade.color="#F2F2F2",
  shade.alpha=0.80,
  matrix.dot.alpha=0.18,
  point.size=2.4,
  line.size=0.7,
  mainbar.y.label="Intersection Size",
  sets.x.label="Set Size",
  set_size.show=FALSE,
  number.angles=0,
  text.scale=c(2.35,2.10,2.35,2.00,2.40,1.90)
)

# The 7 x 7 inch canvas matches the original F3D source PDF.
temporary_pdf <- tempfile(pattern="F3D_upset_",fileext=".pdf")
temporary_png <- tempfile(pattern="F3D_upset_",fileext=".png")

cairo_pdf(
  filename=temporary_pdf,
  width=7,
  height=7,
  family="Arial",
  onefile=TRUE
)

# UpSetR normally starts a new grid page when printed. On a PDF device this
# can leave page 1 blank, so newpage=FALSE draws on the current first page.
print(upset_plot,newpage=FALSE)
dev.off()

png(
  filename=temporary_png,
  width=2100,
  height=2100,
  units="px",
  res=300,
  type="cairo",
  family="Arial"
)

print(upset_plot,newpage=FALSE)
dev.off()

pdf_output <- paste0(output_prefix,".pdf")
png_output <- paste0(output_prefix,".png")

if (!file.copy(temporary_pdf,pdf_output,overwrite=TRUE)) {
  stop("Failed to copy the PDF to the F3 output directory.")
}

if (!file.copy(temporary_png,png_output,overwrite=TRUE)) {
  stop("Failed to copy the PNG to the F3 output directory.")
}

unlink(c(temporary_pdf,temporary_png))

message("PDF saved to: ",pdf_output)
message("PNG saved to: ",png_output)
message(
  "Elapsed time: ",
  round(as.numeric(difftime(Sys.time(),start_time,units="secs")),2),
  " seconds"
)
