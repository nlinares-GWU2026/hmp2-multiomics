# =============================================================
# Script 06: Publication-quality visualizations
# =============================================================
# Generates improved versions of all 5 mixOmics plots:
#   Plot 2: Correlation circle — ggrepel labels, direction arrows
#   Plot 3: Microbiome loadings — colored bars, direction labels,
#            descriptive title
#   Plot 4: Metabolomics loadings — same improvements
#   Plot 5: Correlation heatmap — pheatmap with row/col
#            annotations, full species names, clear color legend
# =============================================================

library(ggplot2)
library(ggrepel)
library(dplyr)
library(pheatmap)

cat("=== Loading data and results ===\n")

# Load saved sparse CCA result from script 05
spls_result <- readRDS("results/mixomics/spls_result.rds")

# Reload data to recompute correlations
mgx  <- read.csv("data/processed/microbiome_clr.csv",
                 row.names=1, check.names=FALSE)
mbx  <- read.csv("data/processed/metabolomics_log2.csv",
                 row.names=1, check.names=FALSE)
meta <- read.csv("data/processed/metadata_filtered.csv",
                 check.names=FALSE)
rownames(meta) <- meta$sample_id

# Reconstruct the same X, Y matrices used in script 05
metabolite_var <- apply(mbx, 2, var)
top500         <- names(sort(metabolite_var, decreasing=TRUE))[1:500]
mbx_top500     <- mbx[, top500]
common_samples <- Reduce(intersect, list(
  rownames(mgx), rownames(mbx_top500), meta$sample_id))

X         <- as.matrix(mgx[common_samples, ])
Y         <- as.matrix(mbx_top500[common_samples, ])
colnames(X) <- make.names(colnames(X))
colnames(Y) <- make.names(colnames(Y))
diagnosis <- meta[common_samples, "diagnosis"]

# Helper: clean species names for display
# Converts "s__Alistipes_putredinis" → "Alistipes putredinis"
clean_species <- function(x) {
  x <- gsub("^s__", "", x)   # remove s__ prefix
  x <- gsub("_",   " ", x)   # underscores to spaces
  x
}

dir.create("results/mixomics", showWarnings=FALSE, recursive=TRUE)

# =============================================================
# Plot 2: Correlation circle with ggrepel + direction labels
# =============================================================
cat("\n=== Plot 2: Correlation circle ===\n")

# Correlations between each variable and the latent components
# These are the (x,y) positions of each feature on the circle
cor_X <- cor(X, spls_result$variates$X)
cor_Y <- cor(Y, spls_result$variates$Y)

# Keep only features that were selected (non-zero loadings)
sel_X <- rownames(spls_result$loadings$X)[
  apply(spls_result$loadings$X, 1, function(r) any(r != 0))]
sel_Y <- rownames(spls_result$loadings$Y)[
  apply(spls_result$loadings$Y, 1, function(r) any(r != 0))]

df_species <- data.frame(
  x     = cor_X[sel_X, 1],
  y     = cor_X[sel_X, 2],
  label = clean_species(sel_X),
  type  = "Microbiome species",
  stringsAsFactors = FALSE
)

df_metabolites <- data.frame(
  x     = cor_Y[sel_Y, 1],
  y     = cor_Y[sel_Y, 2],
  label = sel_Y,
  type  = "Metabolite",
  stringsAsFactors = FALSE
)

df_all <- rbind(df_species, df_metabolites)

# Build circle outlines (inner = r=0.5, outer = r=1.0)
theta      <- seq(0, 2*pi, length.out=300)
circ_outer <- data.frame(x=cos(theta), y=sin(theta))
circ_inner <- data.frame(x=0.5*cos(theta), y=0.5*sin(theta))

p2 <- ggplot(df_all, aes(x=x, y=y, color=type)) +
  # Circles
  geom_path(data=circ_outer, aes(x=x, y=y),
            inherit.aes=FALSE, color="gray45", linewidth=0.5) +
  geom_path(data=circ_inner, aes(x=x, y=y),
            inherit.aes=FALSE, color="gray65", linewidth=0.4,
            linetype="dashed") +
  # Axis crosshairs
  geom_hline(yintercept=0, color="gray60", linetype="dashed",
             linewidth=0.4) +
  geom_vline(xintercept=0, color="gray60", linetype="dashed",
             linewidth=0.4) +
  # Feature points
  geom_point(size=2.5, alpha=0.85) +
  # Species labels only — ggrepel prevents overlap
  geom_label_repel(
    data    = filter(df_all, type == "Microbiome species"),
    aes(label = label),
    size         = 3,
    fontface     = "italic",
    max.overlaps = 50,
    box.padding  = 0.5,
    point.padding = 0.2,
    segment.color = "gray50",
    segment.size  = 0.3,
    show.legend  = FALSE
  ) +
  # Component 1 direction arrows below the circle
  annotate("segment",
           x=0.08, xend=0.95, y=-1.12, yend=-1.12,
           arrow=arrow(length=unit(0.2,"cm"), type="closed"),
           color="#C0392B", linewidth=0.9) +
  annotate("text", x=0.52, y=-1.21,
           label="IBD-enriched →",
           color="#C0392B", fontface="bold", size=3.5) +
  annotate("segment",
           x=-0.08, xend=-0.95, y=-1.12, yend=-1.12,
           arrow=arrow(length=unit(0.2,"cm"), type="closed"),
           color="#2980B9", linewidth=0.9) +
  annotate("text", x=-0.52, y=-1.21,
           label="← Health-associated",
           color="#2980B9", fontface="bold", size=3.5) +
  scale_color_manual(
    values = c("Microbiome species"="darkgreen",
               "Metabolite"="darkorange")) +
  coord_fixed(xlim=c(-1.15, 1.15), ylim=c(-1.3, 1.15)) +
  labs(
    title    = "Correlation circle: microbiome-metabolite co-variation",
    subtitle = paste0(
      "Component 1 (x-axis) captures IBD vs health distinction  |  ",
      "Features near outer ring (r = 1.0) have strongest signal\n",
      "Green = microbiome species  |  Orange = metabolite features"),
    x     = "Component 1",
    y     = "Component 2",
    color = ""
  ) +
  theme_bw(base_size=12) +
  theme(
    legend.position  = "bottom",
    plot.title       = element_text(face="bold", size=13),
    plot.subtitle    = element_text(size=10, color="gray30"),
    panel.grid       = element_blank()
  )

ggsave("results/mixomics/02_correlation_circle.png",
       p2, width=11, height=12, dpi=150)
cat("Saved Plot 2\n")

# =============================================================
# Plot 3: Microbiome loadings — colored, labeled, better title
# =============================================================
cat("\n=== Plot 3: Microbiome loadings ===\n")

load_mgx <- data.frame(
  species = rownames(spls_result$loadings$X),
  loading = spls_result$loadings$X[, 1],
  stringsAsFactors = FALSE
) %>%
  filter(loading != 0) %>%
  mutate(
    label     = clean_species(species),
    direction = ifelse(loading > 0,
                       "IBD-enriched",
                       "Health-associated")
  ) %>%
  arrange(loading) %>%
  mutate(label = factor(label, levels=label))

p3 <- ggplot(load_mgx,
             aes(x=loading, y=label, fill=direction)) +
  geom_col(width=0.72, color="white", linewidth=0.3) +
  geom_vline(xintercept=0, color="gray20", linewidth=0.6) +
  # Direction labels at top of plot
  annotate("text", x=0.38, y=Inf, vjust=2,
           label="IBD-enriched →",
           color="#C0392B", fontface="bold", size=4, hjust=1) +
  annotate("text", x=-0.38, y=Inf, vjust=2,
           label="← Health-associated",
           color="#2980B9", fontface="bold", size=4, hjust=0) +
  scale_fill_manual(
    values = c("IBD-enriched"     = "#E8693A",
               "Health-associated" = "#4A90D9")) +
  scale_x_continuous(limits=c(-0.52, 0.52),
                     breaks=seq(-0.4, 0.4, 0.1)) +
  labs(
    title    = "Bacterial species driving IBD vs. healthy gut distinction",
    subtitle = paste0(
      "Sparse CCA Component 1 loadings  |  ",
      "Orange = enriched in IBD  |  Blue = depleted in IBD (health-associated)\n",
      "Bar length = strength of contribution to the IBD–health axis"),
    x    = "Loading weight (Component 1)",
    y    = "",
    fill = ""
  ) +
  theme_bw(base_size=12) +
  theme(
    legend.position = "none",
    plot.title      = element_text(face="bold", size=13),
    plot.subtitle   = element_text(size=10, color="gray30"),
    axis.text.y     = element_text(face="italic", size=11),
    panel.grid.major.y = element_blank(),
    plot.margin = margin(t=35, r=10, b=10, l=10)
  ) +
   coord_cartesian(clip="off")

ggsave("results/mixomics/03_loadings_microbiome_comp1.png",
       p3, width=13, height=8, dpi=150)
cat("Saved Plot 3\n")

# =============================================================
# Plot 4: Metabolomics loadings — top 15 each direction
# =============================================================
cat("\n=== Plot 4: Metabolomics loadings ===\n")

load_mbx <- data.frame(
  metabolite = rownames(spls_result$loadings$Y),
  loading    = spls_result$loadings$Y[, 1],
  stringsAsFactors = FALSE
) %>%
  filter(loading != 0) %>%
  mutate(direction = ifelse(
    loading > 0,
    "Co-elevated with IBD bacteria",
    "Co-depleted with IBD bacteria\n(health-associated)"))

# Show 15 most extreme in each direction for readability
top_pos <- load_mbx %>%
  filter(loading > 0) %>%
  slice_max(order_by=loading, n=15)
top_neg <- load_mbx %>%
  filter(loading < 0) %>%
  slice_min(order_by=loading, n=15)

load_mbx_sub <- bind_rows(top_pos, top_neg) %>%
  arrange(loading) %>%
  mutate(metabolite = factor(metabolite, levels=metabolite))

p4 <- ggplot(load_mbx_sub,
             aes(x=loading, y=metabolite, fill=direction)) +
  geom_col(width=0.72, color="white", linewidth=0.3) +
  geom_vline(xintercept=0, color="gray20", linewidth=0.6) +
  annotate("text", x=max(load_mbx_sub$loading),
         y=Inf, vjust=2,
         label="Co-elevated with IBD bacteria →",
         color="#C0392B", fontface="bold", size=3.5, hjust=1) +
annotate("text", x=min(load_mbx_sub$loading),
         y=Inf, vjust=2,
         label="← Depleted in IBD (health-associated)",
         color="#2980B9", fontface="bold", size=3.5, hjust=0) +
  scale_fill_manual(
    values = c(
      "Co-elevated with IBD bacteria"              = "#E8693A",
      "Co-depleted with IBD bacteria\n(health-associated)" = "#4A90D9")) +
  labs(
    title    = "Metabolite features co-varying with IBD-associated microbiome shifts",
    subtitle = paste0(
      "Top 15 most extreme features in each direction  |  ",
      "Sparse CCA Component 1 loadings\n",
      "HILp = HILIC positive mode (lipids/bile acids)  |  ",
      "HILn = HILIC negative  |  C18n = C18 negative (hydrophobic metabolites)"),
    x    = "Loading weight (Component 1)",
    y    = "Metabolite feature ID",
    fill = ""
  ) +
  theme_bw(base_size=12) +
   theme(
    legend.position = "none",
    plot.title      = element_text(face="bold", size=12),
    plot.subtitle   = element_text(size=10, color="gray30"),
    axis.text.y     = element_text(size=9),
    panel.grid.major.y = element_blank(),
    plot.margin     = margin(t=35, r=10, b=10, l=10)
  ) +
  coord_cartesian(clip="off")

ggsave("results/mixomics/04_loadings_metabolomics_comp1.png",
       p4, width=14, height=10, dpi=150)
cat("Saved Plot 4\n")

# =============================================================
# Plot 5: Heatmap with pheatmap — full names, annotations
# =============================================================
cat("\n=== Plot 5: Correlation heatmap ===\n")

# Get selected features (non-zero loadings on Component 1)
sel_X_names <- rownames(spls_result$loadings$X)[
  spls_result$loadings$X[,1] != 0]
sel_Y_names <- rownames(spls_result$loadings$Y)[
  spls_result$loadings$Y[,1] != 0]

# Compute pairwise Pearson correlations between selected
# species and metabolites across all 388 samples
cor_mat <- cor(X[, sel_X_names], Y[, sel_Y_names])

# Clean species names for row labels
rownames(cor_mat) <- clean_species(rownames(cor_mat))

# Row annotation: IBD-enriched vs health-associated
row_annot <- data.frame(
  Direction = ifelse(
    spls_result$loadings$X[sel_X_names, 1] > 0,
    "IBD-enriched",
    "Health-associated"),
  row.names = clean_species(sel_X_names)
)

# Column annotation: mass spec measurement mode
col_annot <- data.frame(
  Mode = ifelse(grepl("^HILp", colnames(cor_mat)), "HILIC positive",
         ifelse(grepl("^HILn", colnames(cor_mat)), "HILIC negative",
         ifelse(grepl("^C18n", colnames(cor_mat)), "C18 negative",
                "Other"))),
  row.names = colnames(cor_mat)
)

# Annotation color scheme
ann_colors <- list(
  Direction = c(
    "IBD-enriched"      = "#E8693A",
    "Health-associated" = "#4A90D9"),
  Mode = c(
    "HILIC positive" = "#8E44AD",
    "HILIC negative" = "#27AE60",
    "C18 negative"   = "#E67E22",
    "Other"          = "gray70")
)

pheatmap(
  cor_mat,
  color  = colorRampPalette(
    c("#2166AC", "#92C5DE", "white", "#F4A582", "#B2182B"))(100),
  breaks = seq(-0.6, 0.6, length.out=101),
  annotation_row    = row_annot,
  annotation_col    = col_annot,
  annotation_colors = ann_colors,
  annotation_names_row = FALSE,
  annotation_names_col = TRUE,
  show_colnames  = FALSE,    # too many metabolite IDs
  show_rownames  = TRUE,
  fontsize_row   = 10,
  fontsize       = 11,
  cellwidth      = 10,
  main = paste0(
    "Microbiome-metabolite correlations (Component 1)\n",
    "Red = co-elevated  |  Blue = opposing pattern  |  ",
    "Orange rows = IBD-enriched bacteria  |  Blue rows = health-associated bacteria"),
  filename = "results/mixomics/05_cim_heatmap.png",
  width    = 18,
  height   = 9
)

cat("Saved Plot 5\n")
cat("\n=== Script 06 complete! All improved plots saved ===\n")
