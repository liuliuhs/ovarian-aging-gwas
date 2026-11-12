library(ggplot2)
library(readxl)
library(dplyr)
library(ggtext)   # Support Markdown and rich-text legend labels.

raw <- read_excel("文章/S5.xlsx", skip = 1)
colnames(raw)[4] <- "Consequence"
counts <- raw %>%
  count(Consequence, name = "count") %>%
  mutate(consequence = case_when(
    Consequence == "intron_variant" ~ "Intron variant",
    Consequence == "upstream_gene_variant" ~ "Upstream gene variant",
    Consequence == "downstream_gene_variant" ~ "Downstream gene variant",
    Consequence == "intron_variant,non_coding_transcript_variant" ~ "Intron + non-coding transcript",
    Consequence == "intron_variant,NMD_transcript_variant" ~ "Intron + NMD transcript",
    Consequence == "regulatory_region_variant" ~ "Regulatory region variant",
    Consequence == "intergenic_variant" ~ "Intergenic variant",
    grepl("missense", Consequence) ~ "Missense variant",
    TRUE ~ "Other"
  )) %>%
  group_by(consequence) %>%
  summarise(count = sum(count)) %>%
  ungroup() %>%
  arrange(desc(count))

total <- sum(counts$count)
counts$pct <- counts$count / total * 100
counts$consequence <- factor(counts$consequence, levels = counts$consequence)
counts$label <- ifelse(counts$pct > 4.5,
                       paste0(round(counts$pct, 1), "%"), "")

colors <- c(
  "Intron variant"                 = "#2a78d6",
  "Upstream gene variant"          = "#1baf7a",
  "Downstream gene variant"        = "#eda100",
  "Intron + non-coding transcript" = "#4a3aa7",
  "Intron + NMD transcript"        = "#e34948",
  "Regulatory region variant"      = "#e87ba4",
  "Intergenic variant"             = "#eb6834",
  "Missense variant"               = "#008300",
  "Other"                          = "#888780"
)

# Generate legend labels with percentages.
legend_labels <- setNames(
  paste0(levels(counts$consequence), "  (",
         round(counts$pct[match(levels(counts$consequence),
                                counts$consequence)], 1), "%)"),
  levels(counts$consequence)
)

# Emphasize the missense-variant legend entry using HTML formatting.
legend_labels["Missense variant"] <- paste0(
  "<b><span style='font-size:11pt'>",
  legend_labels["Missense variant"],
  "</span></b>"
)

p <- ggplot(counts, aes(x = 2, y = count, fill = consequence)) +
  geom_col(color = "white", linewidth = 0.8, width = 1) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    color = "white", size = 3.2, fontface = "bold"
  ) +
  coord_polar(theta = "y", start = 0) +
  xlim(c(0.5, 2.5)) +
  scale_fill_manual(
    values = colors,
    name   = "Consequence",
    breaks = levels(counts$consequence),
    labels = legend_labels
  ) +
  theme_void() +
  theme(
    legend.position  = "right",
    legend.title     = element_text(size = 10, face = "bold", color = "#222"),
    legend.text      = element_markdown(size = 9, color = "#333"),  # Render rich-text legend labels.
    legend.key.size  = unit(0.4, "cm"),
    plot.margin      = margin(10, 10, 10, 10)
  ) +
  guides(fill = guide_legend(ncol = 1))

# Overlay the central annotation using cowplot.
library(cowplot)
library(grid)
center_text <- ggdraw() +
  draw_label("5,180",       x = 0.5, y = 0.54,
             fontface = "bold", size = 18, color = "#1a1a1a")
final <- ggdraw() +
  draw_plot(p) +
  draw_plot(center_text,
            x = 0.04, y = 0.22,   # Align the annotation with the donut center.
            width = 0.5, height = 0.5)
print(final)
