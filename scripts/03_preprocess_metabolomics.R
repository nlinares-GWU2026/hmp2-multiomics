# =============================================================
# Script 03: Preprocess Metabolomics Data
# =============================================================
# What this script does:
#   1. Reads the HMP2 metabolomics BIOM file
#   2. Extracts the abundance matrix
#   3. Transposes so samples are rows, metabolites are columns
#   4. Cleans sample IDs to match microbiome
#   5. Keeps only samples already in our microbiome dataset
#   6. Removes metabolites missing in > 80% of samples
#   7. Applies log2(x + 1) transformation
#   8. Saves processed file to data/processed/
# =============================================================

library(biomformat)
library(dplyr)

cat("=== Step 1: Reading metabolomics BIOM file ===\n")
# BIOM is a standard bioinformatics format for feature tables
# It stores both the data matrix and metadata about features/samples
mbx_biom <- read_biom("data/raw/metabolomics/HMP2_metabolomics.biom")

# Extract the raw abundance matrix
# biom_data() pulls out just the numbers as a sparse matrix
mbx_mat <- as.matrix(biom_data(mbx_biom))
cat("Raw dimensions:", nrow(mbx_mat), "metabolites x", ncol(mbx_mat), "samples\n")
cat("First 3 metabolite IDs:\n")
print(rownames(mbx_mat)[1:3])
cat("First 3 sample IDs:\n")
print(colnames(mbx_mat)[1:3])

# =============================================================
# Step 2: Transpose — samples become rows
# =============================================================
cat("\n=== Step 2: Transposing matrix ===\n")
mbx_t <- as.data.frame(t(mbx_mat))
cat("After transpose:", nrow(mbx_t), "samples x", ncol(mbx_t), "metabolites\n")

# =============================================================
# Step 3: Clean sample IDs
# =============================================================
# Check what the sample IDs look like
cat("\n=== Step 3: Checking sample IDs ===\n")
cat("Example metabolomics sample IDs:\n")
print(head(rownames(mbx_t), 5))

# Load the microbiome sample IDs we already processed
# We only want samples that exist in BOTH datasets
mgx_meta <- read.csv("data/processed/metadata_filtered.csv",
                      check.names = FALSE)
cat("\nExample microbiome sample IDs:\n")
print(head(mgx_meta$sample_id, 5))

# =============================================================
# Step 4: Match to microbiome samples
# =============================================================
cat("\n=== Step 4: Matching samples to microbiome dataset ===\n")

# Try direct match first
common_direct <- intersect(rownames(mbx_t), mgx_meta$sample_id)
cat("Direct match:", length(common_direct), "samples\n")

# If direct match is low, try stripping suffixes
mbx_ids_clean <- sub("_[A-Z]+$", "", rownames(mbx_t))
common_clean <- intersect(mbx_ids_clean, mgx_meta$sample_id)
cat("After ID cleaning:", length(common_clean), "samples\n")

# Use whichever matching strategy worked better
if (length(common_direct) >= length(common_clean)) {
  cat("Using direct sample ID matching\n")
  mbx_matched <- mbx_t[rownames(mbx_t) %in% common_direct, ]
  rownames(mbx_matched) <- rownames(mbx_t)[rownames(mbx_t) %in% common_direct]
  common_samples <- common_direct
} else {
  cat("Using cleaned sample ID matching\n")
  rownames(mbx_t) <- mbx_ids_clean
  common_samples <- common_clean
  mbx_matched <- mbx_t[rownames(mbx_t) %in% common_samples, ]
}

cat("Samples matched to microbiome:", nrow(mbx_matched), "\n")

# =============================================================
# Step 5: Filter metabolites missing in > 80% of samples
# =============================================================
# Metabolomics data has lots of zeros/missing values
# Features missing in most samples can't be reliably analyzed
cat("\n=== Step 5: Filtering missing metabolites ===\n")

# A value of 0 in metabolomics means "not detected"
missing_rate <- colSums(mbx_matched == 0) / nrow(mbx_matched)
mbx_filtered <- mbx_matched[, missing_rate <= 0.80]
cat("Metabolites before filter:", ncol(mbx_matched), "\n")
cat("Metabolites after filter:", ncol(mbx_filtered), "\n")

# =============================================================
# Step 6: Log2 transformation
# =============================================================
# Metabolomics intensities span many orders of magnitude
# Log2(x + 1) compresses the range and makes distributions normal
# The +1 pseudocount handles zeros safely
cat("\n=== Step 6: Applying log2(x+1) transformation ===\n")
mbx_log <- log2(mbx_filtered + 1)

# Quick check: what does the data range look like now?
cat("Value range before log:", round(min(mbx_filtered)), "to",
    round(max(mbx_filtered)), "\n")
cat("Value range after log:", round(min(mbx_log), 2), "to",
    round(max(mbx_log), 2), "\n")

# =============================================================
# Step 7: Save output
# =============================================================
cat("\n=== Step 7: Saving processed metabolomics file ===\n")
write.csv(mbx_log,
          "data/processed/metabolomics_log2.csv")
cat("Saved: data/processed/metabolomics_log2.csv\n")
cat("\n=== Preprocessing complete! ===\n")
cat("Final dataset:", nrow(mbx_log), "samples x",
    ncol(mbx_log), "metabolites\n")
