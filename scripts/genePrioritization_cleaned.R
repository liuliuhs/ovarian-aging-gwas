library(tidyverse)
library(gprofiler2)
library(biomaRt)
library(rvest)
library(xml2)
library(glue)
library(corrplot)
library(vautils)
library(purrr)
library(broom)
# Stage 1: integrate genetic evidence
# VEP
vep <- read_tsv("~/SEM_analysis/commonfactor_includeAMHE2/genePrior/vep_leadSNP_ld0.8_missense.txt", col_names=FALSE)
vep.ensID <- unique(vep$X7)
# MAGMA
# MHC hg19, chr6:28,477,797-33,448,354
magma <- read_tsv("~/SEM_analysis/commonfactor_includeAMHE2/FUMA_withOvary/snp2gen/magma.genes.out")
magma.excludeMHC <- magma %>% filter(!(CHR == 6 & START >= 28477797 & STOP <= 33448354))
magma.pbonf0.5 <- magma.excludeMHC %>% filter(P < 0.05/nrow(magma.excludeMHC))
magma.ensID <- magma.pbonf0.5$GENE
# TWAS
twas <- read_tsv("~/SEM_analysis/commonfactor_includeAMHE2/twas/ovaryAging_commonfactor_includeAMHE2.twas.dat.excludeMHC.combine")
twas.pbonf0.5 <- twas %>% filter(TWAS.P < 0.05/nrow(twas))
twas.ensID <- twas.pbonf0.5$ID %>% str_remove_all("\\.\\d+")

# coloc
coloc.pp40.6 <- twas %>% filter(TWAS.P < 0.05/nrow(twas) & COLOC.PP4 >= 0.6)
coloc.ensID <- coloc.pp40.6$ID %>% str_remove_all("\\.\\d+")
# fine mapping
finemapping <- read_tsv("~/SEM_analysis/commonfactor_includeAMHE2/FINEMAP/echolocatorR_leadSNP250kb_finemap.tsv")
write_tsv(finemapping %>% filter(SUSIE.PP > 0.95 | FINEMAP.PP > 0.95) %>% distinct(SNP, .keep_all = T), "~/SEM_analysis/commonfactor_includeAMHE2/FINEMAP/filter_finemap.tsv")
finemapping %>% filter(SUSIE.PP > 0.95 | FINEMAP.PP > 0.95) %>% distinct(SNP, .keep_all = T) %>% filter(FINEMAP.PP > 0.95)
finemapping %>% filter(SUSIE.PP > 0.95 | FINEMAP.PP > 0.95) %>% distinct(SNP, .keep_all = T) %>% filter(FINEMAP.PP > 0.95 & SUSIE.PP > 0.95)
susie.pp0.95 <- finemapping %>% filter(SUSIE.PP > 0.95) %>% dplyr::rename(rsid=SNP, chromosome = CHR, position = POS)
susie.pp0.95.gene <-vautils::find_nearest_gene(susie.pp0.95, flanking = 500, build = "hg19", collapse = FALSE)
susie.pp0.95.gene[susie.pp0.95.gene$distance == "intergenic",]$distance = 0
susie.pp0.95.gene$distance <- as.numeric(susie.pp0.95.gene$distance)
susie.pp0.95.gene.sorted <- susie.pp0.95.gene %>% arrange(rsid, abs(distance))
susie.pp0.95.gene.filter <- susie.pp0.95.gene.sorted %>% distinct(rsid, .keep_all = TRUE)
susie.symbol2ensID <- gconvert(query = unique(susie.pp0.95.gene.filter$GENE), organism = "hsapiens", target="ENSG", mthreshold = Inf, filter_na = TRUE)
# LINC00925=ENSG00000255571  C16orf27=ENSG00000260349 FAM208B=ENSG00000108021 ATP5B=ENSG00000110955 LOC101929679=ENSG00000259416 SUV420H2=ENSG00000133247
# C10orf131 = ENSG00000188649 LOC102577424 = LINC01574 = ENSG00000248859
susie.ensID <- c(susie.symbol2ensID$target, "ENSG00000188649", "ENSG00000248859")

finemap.pp0.95 <- finemapping %>% filter(FINEMAP.PP > 0.95) %>% dplyr::rename(rsid=SNP, chromosome = CHR, position = POS)
finemap.pp0.95.gene <-find_nearest_gene(finemap.pp0.95, flanking = 500, build = "hg19", collapse = FALSE)
finemap.pp0.95.gene[finemap.pp0.95.gene$distance == "intergenic",]$distance = 0
finemap.pp0.95.gene$distance <- as.numeric(finemap.pp0.95.gene$distance)
finemap.pp0.95.gene.sorted <- finemap.pp0.95.gene %>% arrange(rsid, abs(distance))
finemap.pp0.95.gene.filter <- finemap.pp0.95.gene.sorted %>% distinct(rsid, .keep_all = TRUE)
finemap.symbol2ensID <- gconvert(query = unique(finemap.pp0.95.gene.filter$GENE), organism = "hsapiens", target="ENSG", mthreshold = Inf, filter_na = TRUE)

# C16orf27=ENSG00000260349 SUV420H2=ENSG00000133247
# LOC101928154=ENSG00000254211 AGPAT9=ENSG00000138678 LOC101929164=ENSG00000261558 FAM175A=ENSG00000163322
# HMP19=ENSG00000170091 PAPD7=ENSG00000112941 LOC101929144=ENSG00000261615 C20orf196=ENSG00000171984
# rs56201163 map to TNFAIP6=ENSG00000123610 not LOC101929319
# LOC101929225=ENSG00000278719 BRE=ENSG00000158019 LOC102577424 = LINC01574 = ENSG00000248859 HIST2H2AA3=ENSG00000288825
finemap.ensID <- c(finemap.symbol2ensID$target, "ENSG00000254211", "ENSG00000138678", "ENSG00000261558", "ENSG00000163322",
                   "ENSG00000170091", "ENSG00000112941", "ENSG00000261615", "ENSG00000171984",
                   "ENSG00000123610", "ENSG00000278719", "ENSG00000158019", "ENSG00000248859", "ENSG00000288825")

# Stage 2: integrate functional evidence
stage1.ensID <- unique(c(vep.ensID, magma.ensID, twas.ensID, coloc.ensID, susie.ensID, finemap.ensID))
stage1.ensID2symbol <- gconvert(query = stage1.ensID, organism = "hsapiens", target="HGNC", mthreshold = Inf, filter_na = TRUE)
results_df <- data.frame(ensid = character(),
                         symbol = character())
for (ensid in setdiff(stage1.ensID, stage1.ensID2symbol$input)) {
    url <- glue("https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?g={ensid}")
    webpage <- read_html(url)
    text_content <- webpage %>%
    html_nodes(xpath = '/html/body/div[1]/div/div[2]/div[2]/div[1]/div[1]/h1') %>%
    html_text()
    splitted_content <- strsplit(text_content, " ")
    if (length(splitted_content) > 0 && length(splitted_content[[1]]) >= 2) {
        symbol <- splitted_content[[1]][2]
    } else {
        symbol <- NA  # Use NA when the gene symbol cannot be extracted.
    }
    temp_df <- data.frame(ensid = ensid,
                         symbol = symbol)
    results_df <- rbind(results_df, temp_df)
}
results_df[results_df$ensid == "ENSG00000275005",]$symbol <- "RP11-234G16.5"
results_df[results_df$ensid == "ENSG00000279583",]$symbol <- "RP11-368N21.5"
results_df[results_df$ensid == "ENSG00000277476",]$symbol <- "RP11-147L13.13"

stage1.genesymbol <- c(stage1.ensID2symbol$target, (results_df %>% filter(!is.na(symbol)))$symbol)
write_lines(stage1.genesymbol, "/home/yww/SEM_analysis/commonfactor_includeAMHE2/genePrior/stage1_geneSymbol.txt")
write_lines(stage1.ensID, "/home/yww/SEM_analysis/commonfactor_includeAMHE2/genePrior/stage1_geneEnsID.txt")
stage1.mgi <- mgi <- read_tsv("~/SEM_analysis/commonfactor_includeAMHE2/MGI/MGIBatchReport.txt")

library(MPO.db)
library(AnnotationDbi)

mgi_category <- c("MP:0005375", "MP:0005386", "MP:0005385", "MP:0005384", "MP:0005382", "MP:0005381", "MP:0005380", "MP:0005379", "MP:0005378",
  "MP:0005377", "MP:0005397", "MP:0005376", "MP:0005387", "MP:0010771", "MP:0005371", "MP:0005370", "MP:0010768", "MP:0005369",
  "MP:0002006", "MP:0003631", "MP:0001186", "MP:0005367", "MP:0005389", "MP:0005388", "MP:0005390", "MP:0005394", "MP:0005391")
anc_table <- toTable(MPOANCESTOR)
anc_list <- AnnotationDbi::as.list(MPOANCESTOR)

# Map MP terms to ancestor phenotype categories.
# Map MP terms to ancestor phenotype categories.
get_MP_category <- function(mp_id, anc_list, mgi_category) {
  vec <- anc_list[[mp_id]]
  if (is.null(vec)) {
    return(NA) # Return NA for an unknown MP term.
  } else if (length(vec) == 1) {
    return(mp_id) # Return the MP term for a single-element ancestor vector.
  }
  else {
    intersect_vec <- intersect(vec, mgi_category)
    if (length(intersect_vec) > 1) {
        return(paste(intersect_vec, collapse = ";"))
    } else if (length(intersect_vec) == 1) {
        return(intersect_vec)
    } else {
        return(NA)
    }
  }
}
# Assign phenotype categories to the MGI report.
stage1.mgi <- stage1.mgi %>%
  mutate(MP_category = map_chr(`MP ID`, get_MP_category, anc_list = anc_list, mgi_category = mgi_category))
stage1.mgi %>% filter(is.na(`MP_category`) & !is.na(`MP ID`) & Term != "no abnormal phenotype detected") %>% print(n = 25)
stage1.mgi <- stage1.mgi %>%
  mutate(MP_category = case_when(
    !is.na(`MP ID`) & is.na(MP_category) & grepl("limb", Term) ~ "MP:0005371",
    !is.na(`MP ID`) & is.na(MP_category) & grepl("embryo", Term) ~ "MP:0005380",
    !is.na(`MP ID`) & is.na(MP_category) & (grepl("testes", Term) | grepl("ovary", Term)) ~ "MP:0005389",
    `MP ID` == "MP:0014237" ~ "MP:0010768",
    `MP ID` == "MP:0014263" ~ "MP:0005397",
    `MP ID` == "MP:0014348" ~ "MP:0005397",
    `MP ID` == "MP:0014316" ~ "MP:0005385",
    `MP ID` == "MP:0014317" ~ "MP:0005385",
    `MP ID` == "MP:0031496" ~ "MP:0005386",
    `MP ID` == "MP:0031500" ~ "MP:0005391",
    `MP ID` == "MP:0031508" ~ "MP:0005384",
    `MP ID` == "MP:0031509" ~ "MP:0005384",
    `MP ID` == "MP:0014286" ~ "MP:0005379",
    `MP ID` == "MP:0031491" ~ "MP:0005390",
    TRUE ~ MP_category
  ))
# Expand semicolon-separated phenotype categories into separate rows.
stage1.mgi_separated <- stage1.mgi %>% filter(!is.na(`MP ID`) & `MP ID` != "MP:0002169") %>%
  separate_rows(MP_category, sep = ";")
stage1.mgi %>% filter(is.na(`MP_category`))

stage1.mgi.counts <- stage1.mgi_separated %>% distinct(Input, MP_category) %>% count(MP_category)
mgi.counts <- read_csv("/home/yww/SEM_analysis/commonfactor_includeAMHE2/MGI/MPO_counts.txt")

stage1.mgi.counts <- merge(stage1.mgi.counts, mgi.counts, by="MP_category")
stage1.mgi.counts <- stage1.mgi.counts %>% mutate(other_gene_counts = all_gene_counts - n)
p_values <- numeric(nrow(stage1.mgi.counts))

# Test enrichment using the original category-count contingency tables.
# Test enrichment using the original category-count contingency tables.
for (i in 1:nrow(stage1.mgi.counts)) {
  current_row <- stage1.mgi.counts[i, ]

  other_rows <- stage1.mgi.counts[-i, ]

  other_sum <- colSums(other_rows[, c('n', 'other_gene_counts')])

  contingency_table <- matrix(c(current_row$n, other_sum[1], current_row$other_gene_counts, other_sum[2]),
                              nrow = 2,
                              byrow = TRUE)

  fisher_test <- fisher.test(contingency_table, alternative = "greater")

  p_values[i] <- fisher_test$p.value
}

# Append P values and prepare enrichment plots.
stage1.mgi.counts$p <- p_values
stage1.mgi.counts <- stage1.mgi.counts %>%
  mutate(neg_log_p = -log10(p))

stage1.mgi.counts <- stage1.mgi.counts %>%
  arrange(neg_log_p) %>%
  mutate(MP_category_name = factor(MP_category_name, levels = unique(MP_category_name)))
write_tsv(stage1.mgi, "/home/yww/SEM_analysis/commonfactor_includeAMHE2/genePrior/mgi_report.tsv")
write_tsv(stage1.mgi.counts %>% arrange(desc(neg_log_p)), "/home/yww/SEM_analysis/commonfactor_includeAMHE2/genePrior/mgi_enrichment.tsv")
ggplot(stage1.mgi.counts, aes(x=MP_category_name, y=neg_log_p)) +
  geom_col() +
  ylab("-log10P") +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  theme(axis.text.x = element_text(hjust = 0.5, size = 18),  # Set x-axis label size.
        axis.text.y = element_text(hjust = 1, size = 18),    # Set y-axis label size.
        axis.title.x = element_text(size=18,face="bold", hjust=0.5),
        axis.title.y=element_text(size=18, face = "bold", hjust=0.5),
        legend.position = "none") +
  coord_flip() +
  scale_y_continuous(expand=c(0,0), limits = c(0,3))
ggsave("/home/yww/SEM_analysis/commonfactor_includeAMHE2/genePrior/MGI_enrichment.png")

stage1.mgi.filter <- stage1.mgi_separated %>% filter(MP_category %in% c("MP:0005384", "MP:0005389")) %>% distinct(Input, .keep_all = TRUE)
stage1.mgi_separated %>% group_by(MP_category, Term) %>% count() %>% arrange(MP_category, desc(n)) %>% filter(MP_category %in% c("MP:0005384"))
stage1.mgi_separated %>% group_by(MP_category, Term) %>% count() %>% arrange(MP_category, desc(n)) %>% filter(MP_category %in% c("MP:0005389"))
stage1.mgi_separated %>% group_by(MP_category, Term) %>% count() %>% arrange(MP_category, desc(n)) %>% filter(MP_category %in% c("MP:0005380"))
stage1.mgi_separated %>% group_by(MP_category, Term) %>% count() %>% arrange(MP_category, desc(n)) %>% filter(MP_category %in% c("MP:0005391"))

embryo.genes <- stage1.mgi_separated %>% filter(MP_category == "MP:0005380") %>% dplyr::select(Input) %>% pull() %>% unique()
repro.genes <- stage1.mgi_separated %>% filter(MP_category == "MP:0005389") %>% dplyr::select(Input) %>% pull() %>% unique()
cell.genes <- stage1.mgi_separated %>% filter(MP_category == "MP:0005384") %>% dplyr::select(Input) %>% pull() %>% unique()
eyes.genes <- stage1.mgi_separated %>% filter(MP_category == "MP:0005391") %>% dplyr::select(Input) %>% pull() %>% unique()
library(VennDiagram)

venn.plot <- venn.diagram(
  x = list(Embryo = embryo.genes, Reproduction = repro.genes, Cell = cell.genes, Vision = eyes.genes),
  category.names = c("Embryo", "Reproduction", "Cell", "Vision/eye"),
   fill = c("blue", "red", "green", "yellow"),  # Set category colors.
  alpha = 0.5,  # Set fill transparency.
  category.cex = 15,
  cex = 2,
  cat.cex = 2,
  fontface = "bold",
  cat.fontface = "bold",
  filename = NULL
)
pdf("/home/yww/SEM_analysis/commonfactor_includeAMHE2/genePrior/venn_diagram.pdf", width = 12.5, height = 12)  # Set the output filename and figure dimensions.
grid.draw(venn.plot)
dev.off()  # Close the PDF graphics device.

stage1.ensID2symbol %>% filter(target %in% stage1.mgi.filter$Input)

# OTUD7B
stage2.mgi.ensID <- (stage1.ensID2symbol %>% filter(target %in% stage1.mgi.filter$Input))$input

# Gene prioritization
# poi gene
poi.symbol <- read_lines("~/SEM_analysis/commonfactor_includeAMHE2/genePrior/poi_gene_list2.txt")
poi.symbol2ensID <- gostres <- gconvert(query = poi.symbol, organism = "hsapiens", target="ENSG", mthreshold = Inf, filter_na = TRUE)
prioritezed.df1 %>% filter(ensID %in% poi.symbol2ensID$target)
# PTHB1=ENSG00000122507 SGOL2=ENSG00000163535 SMC1B=ENSG00000077935 CSB-PGBD3=ENSG00000225830
poi.ensID <- c(poi.symbol2ensID$target, "ENSG00000122507", "ENSG00000163535", "ENSG00000077935", "ENSG00000225830")
stage2.poi.ensID <- poi.ensID[poi.ensID %in% stage1.ensID]
prioritized.gene <- Reduce(intersect, list(vep.ensID, magma.ensID, twas.ensID, coloc.ensID, susie.ensID, finemap.ensID, stage2.mgi.ensID))

#
allEnsIDs <- unique(c(vep.ensID, magma.ensID, twas.ensID, coloc.ensID, susie.ensID, finemap.ensID, stage2.mgi.ensID))
prioritezed.df <- data.frame(ensID = allEnsIDs)
prioritezed.df$vep <- as.integer(prioritezed.df$ensID %in% vep.ensID)
prioritezed.df$magma <- as.integer(prioritezed.df$ensID %in% magma.ensID)
prioritezed.df$twas <- as.integer(prioritezed.df$ensID %in% twas.ensID)
prioritezed.df$coloc <- as.integer(prioritezed.df$ensID %in% coloc.ensID)
prioritezed.df$susie <- as.integer(prioritezed.df$ensID %in% susie.ensID)
prioritezed.df$finemap <- as.integer(prioritezed.df$ensID %in% finemap.ensID)
prioritezed.df$mgi <- as.integer(prioritezed.df$ensID %in% stage2.mgi.ensID)

prioritezed.df$sum <- rowSums(prioritezed.df[, -1])
prioritezed.df <- prioritezed.df %>% arrange(-sum, ensID)
prioritezed.gene <- (prioritezed.df %>% filter(sum == max(sum)))$ensID
gconvert(query = (prioritezed.df %>% filter(sum ==max(sum)))$ensID, organism = "hsapiens", target="ENSG", mthreshold = Inf, filter_na = TRUE)

enrichmentgostres <- gost(query = (prioritezed.df %>% filter(sum >=3 ))$ensID, organism = "hsapiens")
p <- gostplot(enrichmentgostres, capped = FALSE, interactive = FALSE)

prioritezed.df$ensID %in% genes$ensg

prioritezed.df.filterpoi <- prioritezed.df %>% filter(!(poi == 1 & sum == 1))

prioritezed.ens2symbol <- gconvert(query = prioritezed.df$ensID, organism = "hsapiens", target="HGNC", mthreshold = Inf, filter_na = TRUE)

results_df <- data.frame(ensid = character(),
                         symbol = character())
for (ensid in setdiff(prioritezed.df$ensID, prioritezed.ens2symbol$input)) {
    url <- glue("https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?g={ensid}")
    webpage <- read_html(url)
    text_content <- webpage %>%
    html_nodes(xpath = '/html/body/div[1]/div/div[2]/div[2]/div[1]/div[1]/h1') %>%
    html_text()
    splitted_content <- strsplit(text_content, " ")
    if (length(splitted_content) > 0 && length(splitted_content[[1]]) >= 2) {
        symbol <- splitted_content[[1]][2]
    } else {
        symbol <- NA  # Use NA when the gene symbol cannot be extracted.
    }
    temp_df <- data.frame(ensid = ensid,
                         symbol = symbol)
    results_df <- rbind(results_df, temp_df)
}
results_df[results_df$ensid == "ENSG00000275005",]$symbol <- "RP11-234G16.5"
results_df[results_df$ensid == "ENSG00000277476",]$symbol <- "RP11-147L13.13"
results_df[results_df$ensid == "ENSG00000279583",]$symbol <- "RP11-368N21.5"
results.df <- rbind(results_df, prioritezed.ens2symbol %>% dplyr::select(input, target) %>% rename(input="ensid", target="symbol"))
prioritezed.df1 <- left_join(prioritezed.df, results.df, by= c("ensID" = "ensid"))
prioritezed.genesymbol <- c(prioritezed.ens2symbol$target, (results_df %>% filter(!is.na(symbol)))$symbol)
prioritezed.df1 %>% filter(is.na(symbol))

human_lr <- read_tsv("/home/yww/SEM_analysis/commonfactor_includeAMHE2/genePrior/human_lr_pair.txt")
human_lr %>% filter(ligand_ensembl_gene_id %in% prioritezed.df1$ensID & receptor_ensembl_gene_id %in% prioritezed.df1$ensID)
human_lr %>% filter(ligand_gene_symbol %in% prioritezed.df1$symbol & receptor_gene_symbol %in% prioritezed.df1$symbol)
mouse_lr <- read_tsv("/home/yww/SEM_analysis/commonfactor_includeAMHE2/genePrior/mouse_lr_pair.txt")
mouse_lr$ligand_gene_symbol <- toupper(mouse_lr$ligand_gene_symbol)
mouse_lr$receptor_gene_symbol <- toupper(mouse_lr$receptor_gene_symbol)
mouse_lr %>% filter(ligand_gene_symbol %in% prioritezed.df1$symbol & receptor_gene_symbol %in% prioritezed.df1$symbol)

# GTEx tissue expression analysis
Sys.setenv("VROOM_CONNECTION_SIZE"=131072*10)
gtex_tpm <- read_tsv("~/publicdata/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz", skip=2)
gene_expression <- as.data.frame(gtex_tpm %>% dplyr::select(-1,-2))
rownames(gene_expression) <- gtex_tpm$Name
gtex_meta <- read_tsv("~/SEM_analysis/commonfactor_includeAMHE2/genePrior/gtex_sampleInfo.txt")
ovary_tpm <- read_tsv("~/publicdata/gtex/gene_tpm_2017-06-05_v8_ovary.gct.gz", skip=2)

# Identify ovarian and comparison samples.
ovary_samples <- gtex_meta %>%
                 filter(tissue == "ovary") %>%
                 pull(`#sampleId`)

other_tissues <- gtex_meta %>%
                 filter(tissue != "ovary") %>%
                 pull(tissue) %>%
                 unique()

# Compare ovarian expression with each tissue and report progress.
# Compare ovarian expression with each tissue and report progress.
calculate_wilcox <- function(other_tissue, idx, total) {
  remaining <- total - idx

  message(sprintf("Processing: Ovary vs %s (%d/%d). Remaining tissues: %d", other_tissue, idx, total, remaining))

  other_samples <- gtex_meta %>%
                   filter(tissue == other_tissue) %>%
                   pull(`#sampleId`)

  ovary_expression <- gene_expression %>%
                      dplyr::select(all_of(ovary_samples))
  other_expression <- gene_expression %>%
                      dplyr::select(all_of(other_samples))

  pb <- progress_bar$new(
    format = "  [:bar] :percent :current/:total :eta",
    total = nrow(gene_expression),
    clear = FALSE,
    width = 60
  )

  results <- map_dfr(rownames(gene_expression), function(gene) {
    pb$tick()

    ovary_data <- unlist(ovary_expression[gene, ])
    other_data <- unlist(other_expression[gene, ])
    w_test <- wilcox.test(ovary_data, other_data, exact = FALSE)
    ovary_mean <- mean(ovary_data, na.rm = TRUE)
    other_mean <- mean(other_data, na.rm = TRUE)
    log2fc <- log2(ovary_mean / other_mean)
    tibble(
      gene = gene,
      ovary_mean = ovary_mean,
      other_mean = other_mean,
      p_value = w_test$p.value,
      log2fc = log2fc,
      tissue = other_tissue
    )
  })

  message("\n")

  return(results)
}

# Run and combine tissue comparisons.
total_tissues <- length(other_tissues)

# Run comparisons for each non-ovarian tissue.
results_list <- map2(other_tissues, seq_along(other_tissues), function(tissue, idx) {
  calculate_wilcox(tissue, idx, total_tissues)
})

# Combine tissue comparison results.
results_combined <- bind_rows(results_list)
load("~/SEM_analysis/commonfactor_includeAMHE2/genePrior/gtex_res.RData")
# Adjust P values for multiple testing.
results_combined <- results_combined %>% mutate(gene = sub("\\..*$", "", gene)) %>% filter(gene %in% prioritezed.df$ensID)

results_combined <- results_combined %>%
                    mutate(p_adjust = p.adjust(p_value, method = "BH"))
# Identify ovary-enriched genes using the original thresholds.
# Identify ovary-enriched genes using the original thresholds.
selected_genes <- results_combined %>%
  group_by(gene) %>%
# Summarize mean log2 fold changes, positive comparisons and adjusted P values.
  summarize(
    mean_log2fc = mean(log2fc, na.rm = TRUE),
    log2fc_gt0_count = sum(log2fc > 0, na.rm = TRUE),
    min_gt0_p_adjust = min(ifelse(log2fc > 0, p_adjust, 1), na.rm = TRUE)
  ) %>%
# Require mean log2 fold change > 1, positive changes in at least 48 tissues, and minimum adjusted P < 0.05.
  filter(mean_log2fc > 1, log2fc_gt0_count >= 48, min_gt0_p_adjust < 0.05) %>%
  ungroup()
ovarySp.filter <- selected_genes %>% left_join(genePrior, by = c("gene" = "ensID")) %>% dplyr::select(gene, mean_log2fc, log2fc_gt0_count, min_gt0_p_adjust, symbol)
ovarySp.allmerge.filter <- merge(results_combined, ovarySp.filter, by = "gene") %>% arrange(gene, tissue)
write_tsv(ovarySp.allmerge.filter, "~/SEM_analysis/commonfactor_includeAMHE2/genePrior/ovarySp_allmerge_filter.tsv")
prioritezed.df1$ovarySp <- as.integer(prioritezed.df1$ensID %in% selected_genes$gene)

# Integrate age-associated ovarian expression evidence.
load("~/publicdata/ovary_RNA-seq/meta_stat.RData")
p.up.thr <- quantile(permutation.p.up, 0.05)
p.down.thr <- quantile(permutation.p.down, 0.05)
deg.up.filter <- cumulative_results_up %>% arrange(Binom_p_value_up)  %>% filter(Binom_p_value_up < 0.02388964)
deg.down.filter <- cumulative_results_down %>% arrange(Binom_p_value_down)  %>% filter(Binom_p_value_down < 0.01038358)
deg.up.symbol2ensID <- gconvert(query = unique(deg.up.filter$Gene), organism = "mmusculus", target="ENTREZGENE", mthreshold = Inf, filter_na = TRUE)
deg.down.symbol2ensID <- gconvert(query = unique(deg.down.filter$Gene), organism = "mmusculus", target="ENTREZGENE", mthreshold = Inf, filter_na = TRUE)
deg.up.merge <- deg.up.filter %>% left_join(deg.up.symbol2ensID, by = c("Gene" = "input")) %>% select(Gene, Binom_p_value_up, target) %>% mutate(Direction = "up") %>% rename(Binom_p_value = Binom_p_value_up, Symbol = target)
deg.down.merge <- deg.down.filter %>% left_join(deg.down.symbol2ensID, by = c("Gene" = "input")) %>% select(Gene, Binom_p_value_down, target) %>% mutate(Direction = "down") %>% rename(Binom_p_value = Binom_p_value_down, Symbol = target)
deg.merge.filter <- rbind(deg.up.merge, deg.down.merge) %>% filter(toupper(Symbol) %in% (genePrior %>% filter(deg == 1))$symbol)
write_tsv(deg.merge.filter, "~/SEM_analysis/commonfactor_includeAMHE2/genePrior/deg_merge_filter.tsv")
prioritezed.df1 %>% filter(symbol %in% toupper(deg.up.symbol2ensID$target))
prioritezed.df1 %>% filter(symbol %in% toupper(deg.down.symbol2ensID$target))
prioritezed.df1 %>% filter(ensID %in% cumulative_results_down$Gene)
prioritezed.df1$deg <- as.integer(toupper(prioritezed.df1$symbol) %in% c(toupper(deg.up.symbol2ensID$target), toupper(deg.down.symbol2ensID$target) ))

prioritezed.df1 <- prioritezed.df1 %>% mutate(sum = vep+magma+twas+coloc+susie+finemap+mgi+deg+ovarySp) %>% arrange(desc(sum)) %>% dplyr::select(ensID, symbol, vep, magma, twas, coloc, susie, finemap, mgi, deg, ovarySp, sum)

prioritezed.cor <- cor(prioritezed.df1 %>% dplyr::select(vep, magma, twas, coloc, susie, finemap, mouse_expression))
corrplot(prioritezed.cor)
write_tsv(prioritezed.df1, "/home/yww/SEM_analysis/commonfactor_includeAMHE2/genePrior/genePrioritization.tsv")
load("~/SEM_analysis/commonfactor_includeAMHE2/genePrior/genePrior.RData")

knowngene.symbol <- read_lines("~/SEM_analysis/commonfactor_includeAMHE2/genePrior/ovary_related_gene.list")
knowngene.symbol2ensID <- gconvert(query = knowngene.symbol, organism = "hsapiens", target="ENSG", mthreshold = Inf, filter_na = TRUE)
genePrior %>% filter(ensID %in% knowngene.symbol2ensID$target)
genePrior %>% filter(!ensID %in% knowngene.symbol2ensID$target)
genePrior %>% filter(sum >= 3, ensID %in% knowngene.symbol2ensID$target)
genePrior %>% filter(sum >= 3, ensID %in% knowngene.symbol2ensID$target)
genePrior %>% filter(sum >= 3, !ensID %in% knowngene.symbol2ensID$target)

prioritezed.df1 %>% filter(sum >= 3, ensID %in% poi.ensID)
prioritezed.df1 %>% filter(ensID %in% poi.ensID)

prioritezed.df1 %>% filter(ensID %in% selected_genes$gene)
prioritezed.df1 %>% filter(sum >= 3, ensID %in% selected_genes$gene)
prioritezed.df1 %>% mutate(TWAS.direction = case_when())
# Annotate known ovarian genes and evidence directions.
knowngene <- read_lines("~/SEM_analysis/commonfactor_includeAMHE2/genePrior/ovary_related_gene.list")
knowngene.symbol2ensID <- gconvert(query = knowngene, organism = "hsapiens", target="ENSG", mthreshold = Inf, filter_na = TRUE)
twas.fix <- twas %>% mutate(ensID =  str_remove_all(ID, "\\.\\d+")) %>% dplyr::select(ensID, TWAS.Z, TWAS.P, COLOC.PP4)
prioritezed.df1 <- prioritezed.df1 %>% left_join(twas.fix) %>% mutate(TWAS.direction = case_when(TWAS.Z > 0 ~ "up", TWAS.Z < 0 ~ "down", TRUE~"NA"))
prioritezed.df1 %>% filter(symbol %in% toupper(deg.up.symbol2ensID$target))
prioritezed.df1 %>% filter(symbol %in% toupper(deg.down.symbol2ensID$target))
prioritezed.df1 <- prioritezed.df1 %>%
  mutate(deg.direction = case_when(
    symbol %in% toupper(deg.up.symbol2ensID$target) ~ "up",
    symbol %in% toupper(deg.down.symbol2ensID$target) ~ "down",
    TRUE ~ "NA"
  ))
prioritezed.df1 %>% filter(TWAS.direction == deg.direction & TWAS.direction!="NA")
prioritezed.df1 %>% filter(sum >= 3)
# Add protein MR and colocalization evidence.
mr.symbol <- c("TMEM106A", "NRBP1", "SF3B4", "ITPKA", "PNKP", "VAT1", "GDI2", "GCKR")
mr.direction <- c("up", "down", "down", "up", "down", "up", "up", "up")
mr.dat <- data.frame(mr.symbol, mr.direction)
mr.coloc <- c("TMEM106A", "VAT1")
prioritezed.df1 %>% filter(symbol %in% mr.symbol)
prioritezed.df1$MR <- as.integer(prioritezed.df1$symbol %in% mr.symbol)
prioritezed.df1$pQTLColoc <- as.integer(prioritezed.df1$symbol %in% mr.coloc)
prioritezed.df1 <- prioritezed.df1 %>% left_join(mr.dat, by = c("symbol" = "mr.symbol")) %>% rename(eQTLColoc = coloc) %>% mutate(sum = vep+magma+twas+eQTLColoc+susie+finemap+mgi+deg+ovarySp+MR+pQTLColoc) %>% arrange(desc(sum))
prioritezed.df1 <- prioritezed.df1 %>% dplyr::select(colnames(prioritezed.df1)[1:11], MR, pQTLColoc, sum, TWAS.direction, deg.direction, mr.direction)
write_tsv(prioritezed.df1, "/home/yww/SEM_analysis/commonfactor_includeAMHE2/genePrior/genePrioritization.tsv")
prioritezed.df1 %>% filter(sum >= 4)

# Update existing gene-prioritization tables.

genePrior_withscDEG <- read_tsv("~/SEM_analysis/commonfactor_includeAMHE2/genePrior/genePrior_withscDEG.tsv")
genePrior_withscDEG <- genePrior_withscDEG %>% mutate(mgi = ifelse(ensID %in% stage2.mgi.ensID, 1, 0)) %>% mutate(sum = vep + magma + twas + eQTLColoc + susie + finemap + mgi + deg + ovarySp + MR + pQTLColoc + scdeg) %>% arrange(desc(sum))
write_tsv(genePrior_withscDEG, "~/SEM_analysis/commonfactor_includeAMHE2/genePrior/genePrior_withscDEG_fix.tsv")
genePrior_withscDEG <- read_tsv("~/SEM_analysis/commonfactor_includeAMHE2/genePrior/genePrior_withscDEG_fix.tsv")
genePrior <- read_tsv("~/SEM_analysis/commonfactor_includeAMHE2/genePrior/genePrioritization.tsv")
genePrior %>% filter(twas == 1 & deg == 1 & ovarySp == 1)
genePrior <- genePrior %>% mutate(ovarySp = case_when(ensID %in% ovarySp.filter$gene ~ 1, TRUE ~ 0), mgi = case_when(symbol %in% stage1.mgi.filter$Input ~ 1, TRUE ~ 0))
genePrior <- genePrior %>% mutate(sum = vep + magma + twas +eQTLColoc + susie + finemap + mgi + deg + ovarySp + MR + pQTLColoc) %>% arrange(desc(sum))
write_tsv(genePrior, "~/SEM_analysis/commonfactor_includeAMHE2/genePrior/genePrioritization.tsv")
