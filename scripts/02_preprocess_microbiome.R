# =============================================================
# Script 02: Preprocess Microbiome Taxonomic Profiles
# =============================================================
# What this script does, step by step:
#   1. Reads the MetaPhlAn3 merged taxonomic profiles
#   2. Filters to species level only (removes kingdom/phylum/etc rows)
#   3. Transposes so samples are rows and species are columns
#   4. Cleans sample IDs to match metadata
#   5. Reads metadata and extracts diagnosis + key variables
#   6. Keeps only samples present in both datasets
#   7. Removes rare species (seen in < 10% of samples)
#   8. Applies CLR transformation (standard for microbiome data)
#   9. Saves processed files to data/processed/
# =============================================================

library(dplyr)
library(readr)

cat("=== Step 1: Reading microbiome data ===\n")
mgx_raw <- read_tsv(
  "data/raw/microbiome/taxonomic_profiles.tsv.gz",
  show_col_types = FALSE
)
# Rename the first column (it has a backslash which causes issues)
colnames(mgx_raw)[1] <- "Feature"
cat("Raw dimensions:", nrow(mgx_raw), "taxa x", ncol(mgx_raw)-1, "samples\n")

# =============================================================
# Step 2: Filter to species level only
# =============================================================
# The file contains taxonomy at ALL levels:
#   k__Archaea
#   k__Archaea|p__Euryarchaeota
#   k__Archaea|p__Euryarchaeota|...|s__Methanobrevibacter_smithii  <- we want this
# We keep rows that contain 's__' (species) but NOT 't__' (strain)

cat("\n=== Step 2: Filtering to species level ===\n")
mgx_species <- mgx_raw %>%
  filter(grepl("s__", Feature)) %>%
  filter(!grepl("t__", Feature))

cat("Species-level rows kept:", nrow(mgx_species), "\n")

# Simplify the species name: keep only the 's__Genus_species' part
mgx_species$Feature <- gsub(".*\\|", "", mgx_species$Feature)

# =============================================================
# Step 3: Transpose — samples become rows, species become columns
# =============================================================
# Right now: rows = species, columns = samples
# We need: rows = samples, columns = species (for statistical analysis)

cat("\n=== Step 3: Transposing matrix ===\n")
species_names <- mgx_species$Feature
mgx_mat <- as.data.frame(t(mgx_species[, -1]))
colnames(mgx_mat) <- species_names

# =============================================================
# Step 4: Clean sample IDs
# =============================================================
# Microbiome sample IDs look like: CSM5FZ3N_profile
# We strip '_profile' to get: CSM5FZ3N

cat("\n=== Step 4: Cleaning sample IDs ===\n")
rownames(mgx_mat) <- gsub("_profile$", "", rownames(mgx_mat))
cat("Example cleaned sample IDs:\n")
print(head(rownames(mgx_mat), 5))

# =============================================================
# Step 5: Read metadata properly using R's CSV parser
# =============================================================
# R handles quoted fields with commas correctly — awk doesn't
# check.names=FALSE preserves original column names with spaces

cat("\n=== Step 5: Reading metadata ===\n")
meta <- read.csv(
  "data/raw/hmp2_metadata.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
cat("Metadata dimensions:", nrow(meta), "rows x", ncol(meta), "columns\n")

# Filter to metagenomics samples only
meta_mgx <- meta %>%
  filter(data_type == "metagenomics")
cat("Metagenomics samples in metadata:", nrow(meta_mgx), "\n")

# Clean the External ID to match microbiome IDs
# Metadata: CSM5FZ3N_P  ->  strip '_P'  ->  CSM5FZ3N
meta_mgx$sample_id <- sub("_[A-Z]+$", "", meta_mgx$`External ID`)

# Check the diagnosis values
cat("\nDiagnosis values found:\n")
print(table(meta_mgx$diagnosis))

# =============================================================
# Step 6: Match samples across microbiome and metadata
# =============================================================
cat("\n=== Step 6: Matching samples ===\n")
common_samples <- intersect(rownames(mgx_mat), meta_mgx$sample_id)
cat("Samples in microbiome file:", nrow(mgx_mat), "\n")
cat("Samples in metadata:", nrow(meta_mgx), "\n")
cat("Samples in BOTH:", length(common_samples), "\n")

# Keep only matched samples
mgx_matched <- mgx_mat[common_samples, ]
meta_matched <- meta_mgx[match(common_samples, meta_mgx$sample_id), ]

# =============================================================
# Step 7: Filter rare species
# =============================================================
# Remove species present in fewer than 10% of samples
# These are too rare to find meaningful associations

cat("\n=== Step 7: Filtering rare species ===\n")
prevalence <- colSums(mgx_matched > 0) / nrow(mgx_matched)
mgx_filtered <- mgx_matched[, prevalence >= 0.10]
cat("Species before filter:", ncol(mgx_matched), "\n")
cat("Species after filter:", ncol(mgx_filtered), "\n")

# =============================================================
# Step 8: CLR transformation
# =============================================================
# Microbiome data is compositional (abundances sum to 100%)
# CLR (Centered Log-Ratio) is the standard transformation
# Formula: CLR(x) = log(x / geometric_mean(all species in sample))
# We add a tiny pseudocount (1e-6) to handle zeros

cat("\n=== Step 8: Applying CLR transformation ===\n")
clr_transform <- function(x) {
  x <- x + 1e-6                    # pseudocount for zeros
  log(x) - mean(log(x))            # log minus log geometric mean
}
mgx_clr <- as.data.frame(t(apply(mgx_filtered, 1, clr_transform)))
cat("CLR-transformed matrix:", nrow(mgx_clr), "samples x", ncol(mgx_clr), "species\n")

# =============================================================
# Step 9: Save outputs
# =============================================================
cat("\n=== Step 9: Saving processed files ===\n")

# Save the CLR-transformed microbiome matrix
write.csv(mgx_clr,
          "data/processed/microbiome_clr.csv")

# Save the matched metadata (just the columns we need)
meta_out <- meta_matched %>%
  select(sample_id, diagnosis, `Participant ID`, week_num)
write.csv(meta_out,
          "data/processed/metadata_filtered.csv",
          row.names = FALSE)

cat("Saved: data/processed/microbiome_clr.csv\n")
cat("Saved: data/processed/metadata_filtered.csv\n")
cat("\n=== Preprocessing complete! ===\n")
cat("Final dataset:", nrow(mgx_clr), "samples,", ncol(mgx_clr), "species\n")
cat("Diagnosis breakdown:\n")
print(table(meta_out$diagnosis))
