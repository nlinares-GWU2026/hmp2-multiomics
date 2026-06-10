#!/bin/bash
# =============================================================
# Script 01: Download HMP2 data from IBDMDB
# =============================================================
# Data source: ibdmdb.org (Inflammatory Bowel Disease Multi-omics Database)
# All files are publicly available without registration
# Run from the project root directory: bash scripts/01_download_data.sh
# =============================================================

set -e  # stop immediately if any download fails

echo "Creating data directories..."
mkdir -p data/raw/microbiome
mkdir -p data/raw/metabolomics

echo "Downloading sample metadata..."
curl -L "https://g-227ca.190ebd.75bc.data.globus.org/ibdmdb/metadata/hmp2_metadata_2018-08-20.csv" \
  -o data/raw/hmp2_metadata.csv

echo "Downloading metagenomics merged taxonomic profiles (MetaPhlAn3)..."
curl -L "https://g-227ca.190ebd.75bc.data.globus.org/ibdmdb/products/HMP2/MGX/2018-05-04/taxonomic_profiles_3.tsv.gz" \
  -o data/raw/microbiome/taxonomic_profiles.tsv.gz

echo "Downloading metabolomics feature table..."
curl -L "https://g-227ca.190ebd.75bc.data.globus.org/ibdmdb/products/HMP2/MBX/HMP2_metabolomics.biom" \
  -o data/raw/metabolomics/HMP2_metabolomics.biom

echo "Converting metabolomics BIOM to TSV format..."
biom convert \
  -i data/raw/metabolomics/HMP2_metabolomics.biom \
  -o data/raw/metabolomics/HMP2_metabolomics.tsv \
  --to-tsv

echo ""
echo "=== All downloads complete ==="
ls -lh data/raw/
ls -lh data/raw/microbiome/
ls -lh data/raw/metabolomics/
