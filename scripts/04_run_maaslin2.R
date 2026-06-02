# =============================================================
# Script 04: Run MaAsLin2 — IBD-Associated Microbes & Metabolites
# =============================================================
# This script runs MaAsLin2 twice:
#   Run 1: Which microbiome species are associated with IBD diagnosis?
#   Run 2: Which metabolites are associated with IBD diagnosis?
#
# Why split? With 388 samples we can't fit 500+ predictors at once.
# Instead we ask: what differs between CD, UC, and healthy patients?
# Then in Script 05 (mixOmics) we find which dysregulated microbes
# and metabolites co-vary together — that's the actual integration.
# =============================================================

library(Maaslin2)
library(dplyr)

cat("=== Step 1: Loading processed data ===\n")

# Load microbiome CLR matrix (samples x species)
mgx <- read.csv("data/processed/microbiome_clr.csv",
                row.names = 1, check.names = FALSE)
cat("Microbiome:", nrow(mgx), "samples x", ncol(mgx), "species\n")

# Load metabolomics log2 matrix (samples x metabolites)
mbx <- read.csv("data/processed/metabolomics_log2.csv",
                row.names = 1, check.names = FALSE)
cat("Metabolomics:", nrow(mbx), "samples x", ncol(mbx), "metabolites\n")

# Load metadata
meta <- read.csv("data/processed/metadata_filtered.csv",
                 check.names = FALSE)
rownames(meta) <- meta$sample_id
cat("Metadata:", nrow(meta), "samples\n")

# =============================================================
# Step 2: Select top 500 most variable metabolites
# =============================================================
# Variance = how much a metabolite differs across patients
# High variance metabolites carry the most signal
cat("\n=== Step 2: Selecting top 500 most variable metabolites ===\n")
metabolite_var <- apply(mbx, 2, var)
top500 <- names(sort(metabolite_var, decreasing = TRUE))[1:500]
mbx_top500 <- mbx[, top500]
cat("Selected top 500 of", ncol(mbx), "metabolites by variance\n")

# =============================================================
# Step 3: Find samples present in ALL three datasets
# =============================================================
cat("\n=== Step 3: Aligning samples across datasets ===\n")
common_samples <- Reduce(intersect, list(
  rownames(mgx),
  rownames(mbx_top500),
  meta$sample_id
))
cat("Microbiome samples:", nrow(mgx), "\n")
cat("Metabolomics samples:", nrow(mbx_top500), "\n")
cat("Metadata samples:", nrow(meta), "\n")
cat("Samples in ALL three:", length(common_samples), "\n")

# Align all three datasets to the same samples in the same order
mgx_c  <- mgx[common_samples, ]
mbx_c  <- mbx_top500[common_samples, ]
meta_c <- meta[common_samples, ]

# Build a clean metadata dataframe for MaAsLin2
# participant_id lets MaAsLin2 account for repeated measures
# (same patient sampled multiple times over the year)
meta_df <- data.frame(
  diagnosis      = meta_c$diagnosis,
  participant_id = meta_c$`Participant ID`,
  row.names      = common_samples
)

cat("\nFinal diagnosis breakdown:\n")
print(table(meta_df$diagnosis))

# Clean column names — MaAsLin2 doesn't like spaces or special chars
colnames(mgx_c) <- make.names(colnames(mgx_c))
colnames(mbx_c) <- make.names(colnames(mbx_c))

# =============================================================
# Step 4a: MaAsLin2 Run 1 — Microbiome vs Diagnosis
# =============================================================
# Question: which bacterial species are enriched or depleted in
# CD or UC patients compared to healthy (nonIBD) controls?
#
# fixed_effects  = diagnosis (CD / UC / nonIBD — the main predictor)
# random_effects = participant_id (controls for repeated sampling)
# normalization  = NONE (data is already CLR-transformed)
# transform      = NONE (already transformed)
# analysis_method = LM (linear model — appropriate for CLR data)

cat("\n=== Step 4a: Running MaAsLin2 on microbiome ===\n")
cat("Testing 122 species against IBD diagnosis...\n")
cat("This should take 1-3 minutes\n\n")

dir.create("results/maaslin2/microbiome", recursive = TRUE,
           showWarnings = FALSE)

fit_mgx <- Maaslin2(
  input_data     = mgx_c,
  input_metadata = meta_df,
  output         = "results/maaslin2/microbiome",
  fixed_effects  = "diagnosis",
  random_effects = "participant_id",
reference = "diagnosis,nonIBD",
  normalization  = "NONE",
  transform      = "NONE",
  analysis_method = "LM",
  max_significance = 0.25,
  min_prevalence = 0,
  plot_heatmap   = TRUE,
  plot_scatter   = TRUE,
  cores          = 1
)

cat("\nMicrobiome MaAsLin2 complete!\n")

# =============================================================
# Step 4b: MaAsLin2 Run 2 — Metabolomics vs Diagnosis
# =============================================================
# Question: which metabolites are elevated or depleted in
# CD or UC patients compared to healthy controls?
#
# Same structure as Run 1 but now features = metabolites

cat("\n=== Step 4b: Running MaAsLin2 on metabolomics ===\n")
cat("Testing 500 metabolites against IBD diagnosis...\n")
cat("This should take 3-8 minutes\n\n")

dir.create("results/maaslin2/metabolomics", recursive = TRUE,
           showWarnings = FALSE)

fit_mbx <- Maaslin2(
  input_data     = mbx_c,
  input_metadata = meta_df,
  output         = "results/maaslin2/metabolomics",
  fixed_effects  = "diagnosis",
  random_effects = "participant_id",
reference      = "diagnosis,nonIBD",
  normalization  = "NONE",
  transform      = "NONE",
  analysis_method = "LM",
  max_significance = 0.25,
  min_prevalence = 0,
  plot_heatmap   = TRUE,
  plot_scatter   = FALSE,
  cores          = 1
)

cat("\nMetabolomics MaAsLin2 complete!\n")

# =============================================================
# Step 5: Summarize both results
# =============================================================
cat("\n=== Step 5: Summarizing results ===\n")

# Load results tables
res_mgx <- read.table("results/maaslin2/microbiome/all_results.tsv",
                      header = TRUE, sep = "\t")
res_mbx <- read.table("results/maaslin2/metabolomics/all_results.tsv",
                      header = TRUE, sep = "\t")

# Significant associations (q-value < 0.25 is standard for MaAsLin2)
sig_mgx <- res_mgx %>% filter(qval < 0.25) %>% arrange(qval)
sig_mbx <- res_mbx %>% filter(qval < 0.25) %>% arrange(qval)

cat("\n--- Microbiome results ---\n")
cat("Total species-diagnosis tests:", nrow(res_mgx), "\n")
cat("Significant (q < 0.25):", nrow(sig_mgx), "\n")
if (nrow(sig_mgx) > 0) {
  cat("\nTop 10 IBD-associated species:\n")
  print(sig_mgx[1:min(10, nrow(sig_mgx)),
                c("feature", "metadata", "coef", "qval")])
}

cat("\n--- Metabolomics results ---\n")
cat("Total metabolite-diagnosis tests:", nrow(res_mbx), "\n")
cat("Significant (q < 0.25):", nrow(sig_mbx), "\n")
if (nrow(sig_mbx) > 0) {
  cat("\nTop 10 IBD-associated metabolites:\n")
  print(sig_mbx[1:min(10, nrow(sig_mbx)),
                c("feature", "metadata", "coef", "qval")])
}

# Save significant results
write.csv(sig_mgx,
          "results/maaslin2/significant_species.csv",
          row.names = FALSE)
write.csv(sig_mbx,
          "results/maaslin2/significant_metabolites.csv",
          row.names = FALSE)

cat("\nSaved: results/maaslin2/significant_species.csv\n")
cat("Saved: results/maaslin2/significant_metabolites.csv\n")
cat("\n=== Script 04 complete! ===\n")
cat("These results feed directly into Script 05 (mixOmics)\n")
cat("where we find which dysregulated microbes and metabolites\n")
cat("co-vary together across patients.\n")
