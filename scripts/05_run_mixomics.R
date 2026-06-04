# =============================================================
# Script 05: Sparse CCA with mixOmics
# =============================================================
# What this script does:
#   1. Loads the same 388 matched samples from Scripts 03 & 04
#   2. Prepares X (microbiome, 122 species) and
#                Y (metabolomics, top 500 metabolites)
#   3. Runs sparse PLS in canonical mode (= sparse CCA)
#      to find groups of microbes and metabolites that
#      co-vary together across patients
#   4. Produces four key plots:
#      - Sample plot: do patients cluster by diagnosis?
#      - Correlation circle: which features drive each component?
#      - Circos plot: microbe-metabolite correlation network
#      - Loading plot: ranked feature importance
#   5. Saves all results to results/mixomics/
# =============================================================

library(mixOmics)
library(dplyr)

# =============================================================
# Step 1: Load and align data (same 388 samples as Script 04)
# =============================================================
cat("=== Step 1: Loading data ===\n")

mgx <- read.csv("data/processed/microbiome_clr.csv",
                row.names = 1, check.names = FALSE)
mbx <- read.csv("data/processed/metabolomics_log2.csv",
                row.names = 1, check.names = FALSE)
meta <- read.csv("data/processed/metadata_filtered.csv",
                 check.names = FALSE)
rownames(meta) <- meta$sample_id

# Select top 500 most variable metabolites (same as Script 04)
metabolite_var <- apply(mbx, 2, var)
top500 <- names(sort(metabolite_var, decreasing = TRUE))[1:500]
mbx_top500 <- mbx[, top500]

# Find common samples across all three
common_samples <- Reduce(intersect, list(
  rownames(mgx), rownames(mbx_top500), meta$sample_id
))
cat("Common samples:", length(common_samples), "\n")

# Align datasets
X <- as.matrix(mgx[common_samples, ])       # microbiome (samples x species)
Y <- as.matrix(mbx_top500[common_samples,]) # metabolomics (samples x metabolites)
diagnosis <- meta[common_samples, "diagnosis"]

# Clean column names
colnames(X) <- make.names(colnames(X))
colnames(Y) <- make.names(colnames(Y))

cat("X (microbiome):", nrow(X), "samples x", ncol(X), "species\n")
cat("Y (metabolomics):", nrow(Y), "samples x", ncol(Y), "metabolites\n")
cat("Diagnosis table:\n")
print(table(diagnosis))

# =============================================================
# Step 2: Run sparse PLS in canonical mode (sparse CCA)
# =============================================================
# spls() with mode="canonical" maximises the CORRELATION between
# linear combinations of X and Y — this is sparse CCA.
#
# ncomp = number of components to extract
#   We use 2: Component 1 captures the strongest co-variation
#   pattern, Component 2 captures the next strongest, orthogonal
#   to the first.
#
# keepX = how many microbiome species to keep per component
#   We keep 15 per component — forces sparsity so results are
#   interpretable rather than all 122 species contributing weakly
#
# keepY = how many metabolites to keep per component
#   We keep 50 per component — more metabolites because there
#   are many more of them and metabolite space is richer

cat("\n=== Step 2: Running sparse CCA (spls canonical) ===\n")
cat("ncomp=2, keepX=15, keepY=50 per component\n")
cat("This should take 1-3 minutes...\n\n")

spls_result <- spls(
  X     = X,
  Y     = Y,
  ncomp = 2,
  keepX = c(15, 15),   # 15 species selected per component
  keepY = c(50, 50),   # 50 metabolites selected per component
  mode  = "canonical"  # maximise correlation (= sparse CCA)
)

cat("Sparse CCA complete!\n")

# =============================================================
# Step 3: Check explained variance and correlations
# =============================================================
cat("\n=== Step 3: Checking results ===\n")

# Correlation between microbiome and metabolomics components
cor1 <- cor(spls_result$variates$X[,1],
            spls_result$variates$Y[,1])
cor2 <- cor(spls_result$variates$X[,2],
            spls_result$variates$Y[,2])
cat("Component 1 X-Y correlation:", round(cor1, 3), "\n")
cat("Component 2 X-Y correlation:", round(cor2, 3), "\n")

# Which species were selected on Component 1?
loadings_X1 <- spls_result$loadings$X[,1]
selected_species <- sort(loadings_X1[loadings_X1 != 0],
                         decreasing = TRUE)
cat("\nSpecies selected on Component 1 (top positive):\n")
print(head(selected_species, 5))
cat("Species selected on Component 1 (top negative):\n")
print(tail(selected_species, 5))

# =============================================================
# Step 4: Save all plots to results/mixomics/
# =============================================================
cat("\n=== Step 4: Generating plots ===\n")
dir.create("results/mixomics", recursive = TRUE, showWarnings = FALSE)

# --- Plot 1: Sample plot ---
# Shows each patient as a dot, positioned by Component 1 vs 2
# Colored by diagnosis — if CD/UC/nonIBD cluster separately,
# the co-variation we found is biologically meaningful
cat("Saving sample plot...\n")
png("results/mixomics/01_sample_plot.png",
    width = 800, height = 700, res = 120)
plotIndiv(
  spls_result,
  group      = diagnosis,
  ind.names  = FALSE,
  legend     = TRUE,
  title      = "sPLS samples — colored by diagnosis",
  col.per.group = c("CD" = "#E8693A",
                    "nonIBD" = "#4A90D9",
                    "UC" = "#7AC87A"),
  pch        = 16,
  ellipse    = TRUE   # draws 95% confidence ellipses per group
)
dev.off()

# --- Plot 2: Correlation circle ---
# A circle where features are positioned by their loadings
# Features close together are positively correlated
# Features on opposite sides are negatively correlated
# Features near the edge (r ≈ 1) are strongly selected
cat("Saving correlation circle plot...\n")
png("results/mixomics/02_correlation_circle.png",
    width = 900, height = 900, res = 120)
plotVar(
  spls_result,
  var.names  = list(X = TRUE, Y = FALSE), # show species names, hide metabolite IDs
  cex        = c(3, 2),
  col        = c("darkgreen", "darkorange"),
  title      = "Correlation circle — green=microbiome, orange=metabolites"
)
dev.off()

# --- Plot 3: Loading plot for Component 1 (microbiome side) ---
# Ranked bar chart of which species contribute most to Component 1
# Positive = associated with one group of patients
# Negative = associated with the other group
cat("Saving loading plot (microbiome)...\n")
png("results/mixomics/03_loadings_microbiome_comp1.png",
    width = 900, height = 700, res = 120)
plotLoadings(
  spls_result,
  comp       = 1,
  block      = "X",
  title      = "Component 1 — microbiome loadings",
  contrib    = "max",
  method     = "mean",
  col.ties   = "white"
)
dev.off()

# --- Plot 4: Loading plot for Component 1 (metabolomics side) ---
cat("Saving loading plot (metabolomics)...\n")
png("results/mixomics/04_loadings_metabolomics_comp1.png",
    width = 1100, height = 700, res = 120)
plotLoadings(
  spls_result,
  comp       = 1,
  block      = "Y",
  title      = "Component 1 — metabolomics loadings",
  contrib    = "max",
  method     = "mean",
  col.ties   = "white"
)
dev.off()

# --- Plot 5: Clustered Image Map (CIM) ---
# Heatmap showing correlations between selected microbes and metabolites
# Red = positive correlation (move together)
# Blue = negative correlation (move in opposite directions)
# Rows and columns are clustered so correlated features sit together
cat("Saving clustered image map (CIM)...\n")
png("results/mixomics/05_cim_heatmap.png",
    width = 1200, height = 1000, res = 120)
cim(
  spls_result,
  comp       = 1,
  xlab       = "Metabolites",
  ylab       = "Microbiome species",
  title      = "Component 1 — microbiome-metabolite correlations",
  row.names  = TRUE,
  col.names  = FALSE    # too many metabolite IDs to display
)
dev.off()

# =============================================================
# Step 5: Save numerical results
# =============================================================
cat("\n=== Step 5: Saving numerical results ===\n")

# Microbiome loadings for both components
loadings_mgx <- as.data.frame(spls_result$loadings$X)
colnames(loadings_mgx) <- c("comp1", "comp2")
loadings_mgx$species <- rownames(loadings_mgx)
loadings_mgx <- loadings_mgx %>%
  filter(comp1 != 0 | comp2 != 0) %>%  # keep only selected features
  arrange(desc(abs(comp1)))
write.csv(loadings_mgx,
          "results/mixomics/species_loadings.csv",
          row.names = FALSE)

# Metabolomics loadings for both components
loadings_mbx <- as.data.frame(spls_result$loadings$Y)
colnames(loadings_mbx) <- c("comp1", "comp2")
loadings_mbx$metabolite <- rownames(loadings_mbx)
loadings_mbx <- loadings_mbx %>%
  filter(comp1 != 0 | comp2 != 0) %>%
  arrange(desc(abs(comp1)))
write.csv(loadings_mbx,
          "results/mixomics/metabolite_loadings.csv",
          row.names = FALSE)

cat("Saved: results/mixomics/species_loadings.csv\n")
cat("Saved: results/mixomics/metabolite_loadings.csv\n")
cat("Saved: 5 plots to results/mixomics/\n")

cat("\n=== Script 05 complete! ===\n")
cat("Component 1 X-Y correlation:", round(cor1, 3), "\n")
cat("Component 2 X-Y correlation:", round(cor2, 3), "\n")
cat("Selected", nrow(loadings_mgx), "species and",
    nrow(loadings_mbx), "metabolites total\n")
