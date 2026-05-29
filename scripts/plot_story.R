#!/usr/bin/env Rscript

# plot_story.R — Capsid Integrin Signaling → Nucleolar Stress story
# Plot 1: GSEA pathway timecourse (bubble plot)
# Plot 2: Key gene expression timecourse (line plot)

library(ggplot2)
library(dplyr)
library(tidyr)

RESULTS <- "results/20260525_093658"

# ────────────────────────────────────────────────────────────
# PLOT 1: GSEA pathway timecourse bubble plot
# ────────────────────────────────────────────────────────────

gsea <- read.table(file.path(RESULTS, "pathway/gsea_kegg_signif.tsv"),
                   header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")

# Pathways of interest
target_pathways <- c("Ribosome",
                     "ErbB signaling pathway",
                     "mTOR signaling pathway",
                     "MAPK signaling pathway",
                     "Ras signaling pathway",
                     "Cell cycle",
                     "PI3K-Akt signaling pathway",
                     "TGF-beta signaling pathway",
                     "p53 signaling pathway")

# Contrasts in storyline order
target_contrasts <- c("A549_1h_vs_mock",
                      "E6_1h_vs_mock",
                      "A549_3h_vs_mock",
                      "E6_3h_vs_mock",
                      "A549_3h_vs_1h",
                      "E6_3h_vs_1h")

plot1_data <- gsea %>%
  filter(Description %in% target_pathways, contrast %in% target_contrasts) %>%
  mutate(
    Description = factor(Description, levels = rev(target_pathways)),
    contrast    = factor(contrast, levels = target_contrasts),
    direction   = ifelse(NES > 0, "Up", "Down"),
    absNES      = abs(NES),
    sig         = -log10(p.adjust)
  )

# Group pathways for facet labeling
pathway_groups <- c(
  "Ribosome"                      = "Translation",
  "ErbB signaling pathway"        = "Signaling",
  "mTOR signaling pathway"        = "Signaling",
  "MAPK signaling pathway"        = "Signaling",
  "Ras signaling pathway"         = "Signaling",
  "Cell cycle"                    = "Proliferation",
  "PI3K-Akt signaling pathway"    = "Signaling",
  "TGF-beta signaling pathway"    = "Response",
  "p53 signaling pathway"         = "Response"
)
plot1_data$group <- factor(pathway_groups[as.character(plot1_data$Description)],
                           levels = c("Signaling", "Translation", "Proliferation", "Response"))

# Contrast labels for x-axis
contrast_labels <- c(
  "A549_1h_vs_mock"  = "A549\n1h vs mock",
  "E6_1h_vs_mock"    = "E6\n1h vs mock",
  "A549_3h_vs_mock"  = "A549\n3h vs mock",
  "E6_3h_vs_mock"    = "E6\n3h vs mock",
  "A549_3h_vs_1h"    = "A549\n3h vs 1h",
  "E6_3h_vs_1h"      = "E6\n3h vs 1h"
)

p1 <- ggplot(plot1_data, aes(x = contrast, y = Description)) +
  geom_point(aes(size = absNES, fill = direction, alpha = sig),
             shape = 21, stroke = 0.3) +
  scale_size_area(name = "|NES|", max_size = 12,
                  breaks = c(0.5, 1, 1.5, 2),
                  labels = c("0.5", "1.0", "1.5", "2.0")) +
  scale_fill_manual(name = "Direction", values = c("Up" = "#d73027", "Down" = "#4575b4")) +
  scale_alpha_continuous(name = expression(-log[10](padj)),
                         range = c(0.3, 1),
                         breaks = c(2, 4, 6, 8, 10)) +
  scale_x_discrete(labels = contrast_labels) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y",
             switch = "y") +
  labs(title = "Pathway enrichment timecourse",
       subtitle = "Ad26 capsid → integrin signaling → nucleolar stress cascade",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 0, face = "bold", size = 10),
    strip.placement = "outside",
    legend.position = "bottom",
    legend.box = "vertical",
    legend.margin = margin(t = 4),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40")
  ) +
  guides(
    size  = guide_legend(order = 1, title.position = "top", nrow = 1),
    fill  = guide_legend(order = 2, title.position = "top", nrow = 1),
    alpha = guide_legend(order = 3, title.position = "top", nrow = 1)
  )

ggsave(file.path("plots", "gsea_timeline_bubble.pdf"), p1,
       width = 9, height = 5.5, device = "pdf")

message("Plot 1 saved: plots/gsea_timeline_bubble.pdf")

# ────────────────────────────────────────────────────────────
# PLOT 2: Key gene expression timecourse
# ────────────────────────────────────────────────────────────

# Read VST counts (wide format: gene × samples)
vst <- read.table(file.path(RESULTS, "tables/vst_normalized_counts.tsv"),
                  header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

# Read metadata
meta <- read.table(file.path(RESULTS, "tables/metadata.tsv"),
                   header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Key genes — ENSG → symbol mapping (expanded set)
key_genes <- data.frame(
  gene_id = c(
    # Panel A: mTOR / ribosomal burst (1h)
    "ENSG00000142676", "ENSG00000137154", "ENSG00000187840", "ENSG00000142208",
    "ENSG00000168209", "ENSG00000175634", "ENSG00000115053", "ENSG00000117461",
    # Panel B: p53 / nucleolar stress (3h)
    "ENSG00000105327", "ENSG00000080546", "ENSG00000130766", "ENSG00000149547",
    "ENSG00000135679",
    # Panel C: TGF-beta / ID cascade (3h)
    "ENSG00000092969", "ENSG00000125968", "ENSG00000115738", "ENSG00000117318",
    "ENSG00000172201", "ENSG00000137834", "ENSG00000101665", "ENSG00000183691"
  ),
  symbol  = c(
    "RPL11", "RPS6", "EIF4EBP1", "AKT1",
    "DDIT4", "RPS6KB2", "NCL", "PIK3R3",
    "BBC3", "SESN1", "SESN2", "EI24", "MDM2",
    "TGFB2", "ID1", "ID2", "ID3", "ID4",
    "SMAD6", "SMAD7", "NOG"
  ),
  panel   = c(rep("A — mTOR / ribosomal burst (1h)", 8),
              rep("B — p53 / nucleolar stress (3h)", 5),
              rep("C — TGF-beta / ID cascade (3h)", 8)),
  stringsAsFactors = FALSE
)

# Filter VST to our genes, pivot to long
vst_long <- vst %>%
  filter(gene_id %in% key_genes$gene_id) %>%
  pivot_longer(-gene_id, names_to = "sample_id", values_to = "vst") %>%
  left_join(key_genes, by = "gene_id") %>%
  left_join(meta, by = "sample_id") %>%
  mutate(
    time = factor(time, levels = c("mock", "1h", "3h")),
    cell_line = factor(cell_line, levels = c("A549", "E6")),
    symbol = factor(symbol, levels = key_genes$symbol),
    panel  = factor(panel, levels = unique(key_genes$panel))
  )

# Summary stats: mean ± SEM per group
vst_summary <- vst_long %>%
  group_by(symbol, panel, cell_line, time) %>%
  summarise(
    mean_vst = mean(vst),
    sem_vst  = sd(vst) / sqrt(n()),
    .groups  = "drop"
  )

# Colors
clrs <- c("A549" = "#4575b4", "E6" = "#d73027")

p2 <- ggplot(vst_summary, aes(x = time, y = mean_vst, color = cell_line, group = cell_line)) +
  # Individual replicate points (semi-transparent, smaller)
  geom_point(data = vst_long, aes(y = vst),
             position = position_dodge(width = 0.3),
             size = 1.2, alpha = 0.5) +
  # Mean line
  geom_line(position = position_dodge(width = 0.3), linewidth = 0.8) +
  # Error ribbon
  geom_ribbon(aes(ymin = mean_vst - sem_vst, ymax = mean_vst + sem_vst,
                  fill = cell_line),
              position = position_dodge(width = 0.3),
              alpha = 0.15, color = NA, show.legend = FALSE) +
  # Mean point
  geom_point(position = position_dodge(width = 0.3), size = 2.5, shape = 18) +
  scale_color_manual(name = "Cell line", values = clrs) +
  scale_fill_manual(values = clrs, guide = "none") +
  facet_grid(panel ~ symbol, scales = "free_y", switch = "y") +
  labs(title = "Key gene expression after Ad26 infection",
       subtitle = "VST-normalized counts, n = 3 per group",
       x = "Time post-infection", y = "VST expression") +
  theme_minimal(base_size = 10) +
  theme(
    strip.text.y.left = element_text(angle = 0, hjust = 0, size = 9),
    strip.text.x = element_text(face = "italic", size = 9),
    strip.placement = "outside",
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1, "lines"),
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 8, color = "grey40")
  )

ggsave(file.path("plots", "key_genes_timecourse.pdf"), p2,
       width = 16, height = 8, device = "pdf")

message("Plot 2 saved: plots/key_genes_timecourse.pdf")
message("Done.")
