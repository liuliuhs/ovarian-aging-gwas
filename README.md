# Ovarian aging GWAS selected scripts

Selected analysis and plotting scripts associated with the manuscript
"Multivariate genetic analysis reveals shared architecture and candidate genes for ovarian aging".

## Scope

This is a partial code collection, not a complete end-to-end analysis pipeline.
Server-side workflows and other plotting scripts are not included.
Script filenames retain historical figure numbers; they do not establish the
panel numbering in the final manuscript. Match panels against final figures
before describing these scripts as reproducing the published figures.

## Contents

| Script | Purpose | Input / important caveat |
| --- | --- | --- |
| `F2B_sldsc_enrichment_plot.R` | S-LDSC annotation enrichment plot | LDSC `.results` file; filename uses an older panel number |
| `Figure3D_sixMethods_upsetR_updated.R` | UpSet plot of six genetic evidence categories | Legacy workbook `Table. S22`; also locates a reference PDF to determine output directory |
| `Figure3C_plot_TWAS_eQTL_colocalization_labels_fixed.R` | TWAS/eQTL colocalization plot | Legacy workbook `Table. S16` and `Table. S22` |
| `Figure4D_MR_forest.R` | Protein MR forest plot | See original input declarations and `SCRIPT_INPUTS.md` |
| `Figure4E_pQTL_colocalization.R` | Protein colocalization plot | Input contains `gene_symbol`, `figure_PP_display`, and a candidate-membership indicator; verify source PP values |
| `FigureS4B_TWAS_DEG_concordance.R` | Directional concordance between TWAS and DEG evidence | Legacy `Table. S22`; supports workbook and output-directory command-line arguments |
| `method5_genes.R` | Counts of genes supported by genetic evidence categories | `genePrioritization.tsv`; requires ggplot2 to be loaded in addition to packages declared in the script |
| `vep_missense_variant1.R` | VEP consequence donut plot | Legacy `S5.xlsx`; displays hard-coded central text `5,180`; returns a plot without saving it automatically |
| `genePrioritization_cleaned.R` | Integration of genetic and functional evidence | Retains original local paths, intermediate workspace dependencies and analysis logic; requires further review before standalone use |

## Requirements and execution

Use R and install the packages required by each script. Dependencies include
(depending on script) ggplot2, readxl, dplyr, data.table, UpSetR, ggtext,
cowplot, svglite and ragg. The analysis script has additional dependencies
listed in its `library()` statements. `SCRIPT_INPUTS.md` inventories declared
inputs and dependencies; it is not a complete environment lockfile.

1. Obtain the exact input tables used by the relevant script.
2. Update its input/output path declarations to your local locations.
3. Verify worksheet names and column names. Do not substitute the current
   supplementary workbook without checking its schema: historical scripts
   may refer to earlier supplementary table numbers.
4. Run scripts individually in an R session with the appropriate packages.
   For example, from the repository root:

```r
source("scripts/F2B_sldsc_enrichment_plot.R")
```

The concordance plot also accepts command-line paths:

```text
Rscript scripts/FigureS4B_TWAS_DEG_concordance.R path/to/legacy_workbook.xlsx outputs
```

The code retains original filtering, thresholds, colors, dimensions, labels,
and output behavior. English comments replace Chinese comments. Obsolete
commented code and selected standalone inspection statements were removed.
Original source files were not overwritten.

## Data

Input datasets are not included. See `data/README.md` and `SCRIPT_INPUTS.md`.
Do not distribute restricted data without the applicable authorization.

## Verification status

Retained executable lines in the eight plotting scripts were checked against
the originals in order. No plots were regenerated, no end-to-end R execution
was performed, and agreement with final figure panels has not been established.
R/package versions must be recorded from the actual analysis environment.

## Recording the original environment

Run the following in the R environment actually used for the analysis:

```r
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
```

## Code availability wording

After public release, adapt the following to the actual repository URL:

Selected R scripts for gene prioritization and figure generation are available
at [repository URL]. The repository documents the scope of the shared scripts,
required inputs and software dependencies.

## License

A license has not been selected. Confirm code ownership and coauthor agreement
before adding a license for public distribution.