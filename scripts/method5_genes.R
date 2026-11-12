wyw <- read.delim("C:/Users/92975/Desktop/wyw/genePrioritization.tsv")

library(readr)
library(dplyr)

# Prepare gene counts by evidence category.
counts <- data.frame(
  method = factor(c(
    "VEP\n(Missense variants)",
    "Fine-mapping\n(SuSiE & FINEMAP)",
    "MAGMA\n(Gene-based analysis)",
    "TWAS\n(Ovary eQTL)",
    "eQTL\nColocalization"
  ), levels = c(
    "VEP\n(Missense variants)",
    "Fine-mapping\n(SuSiE & FINEMAP)",
    "MAGMA\n(Gene-based analysis)",
    "TWAS\n(Ovary eQTL)",
    "eQTL\nColocalization"
  )),
  count = c(
    nrow(filter(wyw, vep == 1)),
    nrow(filter(wyw, susie == 1 | finemap == 1)),
    nrow(filter(wyw, magma == 1)),
    nrow(filter(wyw, twas == 1)),
    nrow(filter(wyw, eQTLColoc == 1))
  ),
  color = c("#eb6834", "#4a3aa7", "#2a78d6", "#1baf7a", "#e34948")
)

methods_data <- counts

# Draw the gene-count bar plot.
p <- ggplot(methods_data, aes(x = method, y = count, fill = method)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(
    aes(label = count),
    vjust = -0.5,
    fontface = "bold",
    size = 4,
    color = "#222222"
  ) +
  scale_fill_manual(values = setNames(methods_data$color,
                                      methods_data$method)) +
  scale_y_continuous(
    limits = c(0, 690),
    breaks = seq(0, 600, 100),
    expand = c(0, 0)
  ) +
  labs(x = NULL, y = "Number of genes") +
  theme_classic(base_size = 13, base_family = "sans") +
  theme(
    axis.line         = element_line(color = "#bbbbbb", linewidth = 0.5),
    axis.ticks        = element_line(color = "#bbbbbb", linewidth = 0.4),
    axis.ticks.length = unit(0.2, "cm"),

    # Rotate and right-align x-axis labels.
    axis.text.x       = element_text(size = 11, color = "#333333",
                                     angle = 45, hjust = 1, vjust = 1,
                                     lineheight = 1.3),
    axis.text.y       = element_text(size = 11, color = "#555555"),
    axis.title.y      = element_text(size = 12, color = "#333333",
                                     margin = margin(r = 10)),
    panel.grid.major.y = element_line(color = "#ebebeb",
                                      linewidth = 0.4, linetype = "dashed"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.background   = element_rect(fill = "white", color = NA),
    plot.background    = element_rect(fill = "white", color = NA),

    # Reserve margin space for rotated labels.
    plot.margin = margin(20, 20, 20, 20)
  )

# Export the figures.
ggsave("methods_gene_count.pdf", p,
       width = 8, height = 5, device = cairo_pdf)
ggsave("methods_gene_count.png", p,
       width = 8, height = 5, dpi = 300)

print(p)
