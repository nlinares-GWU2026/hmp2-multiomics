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

# =============================================================
# Script 03: Preprocess Metabolomics Data
# =============================================================
library(dplyr)
library(readr)

cat("=== Step 1: Reading metabolomics TSV file ===\n")
# skip=1 skips the first comment line "# Constructed from biom file"
mbx_raw <- read_tsv(
  "data/raw/metabolomics/HMP2_metabolomics.tsv",
  skip = 1,
  show_col_types = FALSE
)
# Rename the first column from "#OTU ID" to "Feature"
colnames(mbx_raw)[1] <- "Feature"
cat("Raw dimensions:", nrow(mbx_raw), "metabolites x",
    ncol(mbx_raw)-1, "samples\n")
cat("First 3 metabolite IDs:\n")
print(mbx_raw$Feature[1:3])
cat("First 3 sample IDs:\n")
print(colnames(mbx_raw)[2:4])

# =============================================================
# Step 2: Transpose — samples become rows
# =============================================================
cat("\n=== Step 2: Transposing matrix ===\n")
feature_names <- mbx_raw$Feature
mbx_mat <- as.data.frame(t(mbx_raw[, -1]))
colnames(mbx_mat) <- feature_names
cat("After transpose:", nrow(mbx_mat), "samples x",
    ncol(mbx_mat), "metabolites\n")

# =============================================================
# Step 3: Match to microbiome metadata
# =============================================================
# Metabolomics sample IDs are clean (e.g. CSM5FZ3N)
# metadata_filtered.csv sample_id column is also clean
# They should match directly
cat("\n=== Step 3: Matching samples to microbiome dataset ===\n")
meta <- read.csv("data/processed/metadata_filtered.csv",
                 check.names = FALSE)
cat("Microbiome samples available:", nrow(meta), "\n")
cat("Example metabolomics IDs:", head(rownames(mbx_mat), 3), "\n")
cat("Example metadata IDs:", head(meta$sample_id, 3), "\n")

common_samples <- intersect(rownames(mbx_mat), meta$sample_id)
cat("Samples in BOTH datasets:", length(common_samples), "\n")

mbx_matched <- mbx_mat[common_samples, ]

# =============================================================
# Step 4: Filter metabolites missing in > 80% of samples
# =============================================================
cat("\n=== Step 4: Filtering sparse metabolites ===\n")
missing_rate <- colSums(mbx_matched == 0) / nrow(mbx_matched)
mbx_filtered <- mbx_matched[, missing_rate <= 0.80]
cat("Metabolites before filter:", ncol(mbx_matched), "\n")
cat("Metabolites after filter:", ncol(mbx_filtered), "\n")

# =============================================================
# Step 5: Log2 transformation
# =============================================================
# Raw intensities span millions — log2 compresses to a workable range
cat("\n=== Step 5: Applying log2(x+1) transformation ===\n")
mbx_log <- log2(mbx_filtered + 1)
cat("Value range before:", round(min(mbx_filtered)),
    "to", round(max(mbx_filtered)), "\n")
cat("Value range after:", round(min(mbx_log), 2),
    "to", round(max(mbx_log), 2), "\n")

# =============================================================
# Step 6: Save output
# =============================================================
cat("\n=== Step 6: Saving ===\n")
write.csv(mbx_log, "data/processed/metabolomics_log2.csv")
cat("Saved: data/processed/metabolomics_log2.csv\n")
cat("\n=== Preprocessing complete! ===\n")
cat("Final dataset:", nrow(mbx_log), "samples x",
    ncol(mbx_log), "metabolites\n")
