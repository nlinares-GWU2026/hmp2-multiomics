# Project Writeup: HMP2 Multi-Omics Integration
## Microbiome and Metabolomics in Inflammatory Bowel Disease

---

## Overview

This project uses publicly available data from the Human Microbiome Project Phase 2 (HMP2) to explore the relationship between gut bacteria and small molecules (metabolites) in patients with inflammatory bowel disease (IBD). The central question is: when the gut microbiome is disrupted in IBD, which specific bacteria are changing, which chemicals in the gut change alongside them, and do these two types of changes track together across patients?

To answer this, I applied two complementary statistical methods: 1) MaAsLin2 for identifying features individually associated with disease, and 2) sparse Canonical Correlation Analysis (sparse CCA) from the mixOmics R package for finding groups of bacteria and metabolites that co-vary together across patients.  
The MaAsLin2 served to answer the question of which specific bacteria go up or down in a certain disease state (with the allowance of adjustment for covariables), while the sparse CCA found relationships between specific microbes and metabolites (finding hidden patterns without becoming overwhelmed with standard statistical methods).

---

## Background 

### Inflammatory Bowel Disease (IBD)

IBD is a chronic condition characterized by persistent inflammation of the gastrointestinal tract. It has two main forms: Crohn's disease (CD), which can affect any part of the digestive tract, and ulcerative colitis (UC), which is limited to the colon. Both conditions are associated with dysbiosis - a disruption in the normal composition of the gut microbiome, but the exact relationship between the microbial shifts and the chemical (metabolite) environment of the gut is not fully understood. 

### The Multi-Omics Approach

The gut microbiome influences health not just through the bacteria that are present but through the metabolites those bacteria produce. Bacteria ferment dietary fiber, modify bile acids, synthesize vitamins, and produce short-chain fatty acids. All of these circulate through the gut and affect immune function and inflammation. By measuring both the microbial composition (metagenomics) and the chemical landscape composition (metabolomics) from the same patient samples, we can investigate which microbe-metabolite relationships are disrupted in disease.

### The HMP2 Dataset

The Human Microbiome Project Phase 2, also known as the Integrative Human Microbiome Project (iHMP), followed approximately 130 participants across five hospital sites over a period of about one year. Participants included patients with Crohn's disease, ulcerative colitis, and healthy controls (nonIBD). Stool samples were collected every two weeks and subjected to multiple assays simultaneously, making it one of the most comprehensive longitudinal multi-omics datasets in IBD research. Data is publicly available through the IBDMDB portal (ibdmdb.org).

---

## Data

**Source:** IBDMDB (Inflammatory Bowel Disease Multi-omics Database)
**Access:** [ibdmdb.org](https://ibdmdb.org) — publicly available without registration

**Files downloaded:**

- **Metagenomics** — `taxonomic_profiles_3.tsv.gz`
  MetaPhlAn3 merged taxonomic profiles; microbial species abundances per sample

- **Metabolomics** — `HMP2_metabolomics.biom`
  LC-MS metabolite feature intensities per sample; converted to TSV using the Python `biom-format` package before processing in R

- **Metadata** — `hmp2_metadata_2018-08-20.csv`
  Sample-level clinical variables, including diagnosis (CD, UC, nonIBD), participant ID, and collection week

---
