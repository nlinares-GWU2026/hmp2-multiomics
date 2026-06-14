# HMP2 Multi-Omics Integration: Microbiome & Metabolomics in IBD

A bioinformatics portfolio project integrating gut microbiome and metabolomics data from IBD patients to identify disease-associated microbial and chemical signatures using MaAsLin2 and sparase Canonical Correlation Analysis (mixOmics). 

---

## Background

Inflammatory bowel disease (IBD), including Crohn's disease (CD) and ulcerative colitis (UC), is associated with profound disruption of both the gut microbiome and the gut chemical environment. This project asks: *"Which bacteria and metabolites change in IBD, and do those changes track together across patients?"*

To answer this, I used paired metagenomics and metabolomics data from the Human Microbiome Project Phase 2 (HMP2), one of the most comprehensive longitudinal multi-omics IBD datasets publicly available. 

---

## Dataset

**Source:** [IBDMDB — Inflammatory Bowel Disease Multi-omics Database](https://ibdmdb.org)

- **388 matched samples** with paired microbiome + metabolomics data
- **181 Crohn's disease / 102 ulcerative colitis / 105 healthy controls**
- Samples collected longitudinally across 5 hospital sites over ~1 year

---

## Methods 

- **Microbiome preprocessing** (`readr`, `dplyr`) - filtered to 122 species, CLR-transformed
- **Metabolomics preprocessing** (R + Python `biom-format`) - 66,717 features retained, log2-transformed
- **IBD association testing** (`MaAsLin2`) - species and metabolites tested individually against the diagnosis of the individual
- **Multi-omics integration** (`mixOmics` sparse CCA) - co-varying microbe-metabolite groups identified across patients
- **Visualization** (`ggplot2`, `ggrepel`, `pheatmap`) - 5 publication-style figures

## Key Findings

### MaAsLin2 - IBD-Associated Features
- **23 bacterial species** significantly associated with IBD diagnosis (q < 0.25)
- **559 metabolite features** significantly associated with IBD diagnosis (q < 0.25)
- Consistent depletion of short-chain fatty acid producers in IBD patients
- Consistent enrichment of pro-inflammatory *Clostridium* and *Flavonifractor* species

### Sparse CCA - Microbiome-Metabolite Co-variation
- **Component 1 microbiome-metabolome correlation: r = 0.733**
- **Component 2 microbiome-metabolome correlation: r = 0.617**
- Metabolomics separated IBD from healthy patients more cleanly than microbiome alone - consistent with findings in Lloyd-Price et al. (2019)
- IBD-enriched bacteria (*Flavonifractor plautii*, *Ruminococcus gnavus*, *Clostridium symbiosum*) co-vary with a distinct set of elevated metabolites
- Health-associated bacteria (*Alistipes putredinis*, *Alistipes shahii*, *Oscillibacter* spp.) co-vary with a separate set of depleted metabolites

- ---
