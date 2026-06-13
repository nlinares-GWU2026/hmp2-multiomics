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

## Methods 

### Environment

All analysis was performed in R (version 4.3.3) running on Ubuntu 24.04 via WSL2 on Windows. Key packages include: `Maaslin2` (v1.15.1), `mixOmics` (v6.26.0), `dplyr`, `readr`, `ggplot2`, `ggrepel`, `pheatmap`.

### Microbiome Preprocessing (`02_preprocess_microbiome.R`)
The raw taxonomic profiles comprised 932 rows, spanning all taxonomic levels (from kingdom to species). I filtered to species-level rows only (those containing `s__` but not `t__`, the *strain* level), retaining 578 species. Samples were matched to metadata using the External ID field after stripping suffixes, yielding 1,317 matched samples. Species present in fewer than 10%  of samples were removed as too rare to find reliable associations, leaving 122 species. Because microbiome abundance data is compositional (all values sum to 100%), I applied a Centered Log-Ratio (CLR) transformation, which is the standard normalization for this data type. CLR transforms each species abundance by dividing it by the geometric mean (the average used for numbers that multiply together rather than add together) of all species in that sample and then taking the log, making the data approximately normal and suitable for linear models. 

### Metabolomics Preprocessing (`03_preprocess_metabolomics.R`) 

The raw metabolomics data contained 81, 867 chemical features (mass spectrometry peaks) across 546 samples. After matching to samples that also had microbiome data, 388 samples remained. Features missing in more than 80% of samples were removed (indicating that the chemical was detected too rarely to be informative), leaving 66,717 features. A log2(x + 1) transformation was applied to compress the enormous dynamic range of raw mass spectrometry intensities (which spanned from 0 to 1.08 x 10^14) into a workable range of 0 to 46.6.

### Sample Matching for Integration

The final integrated dataset consisted of 388 samples with complete data from both assays: 181 CD, 102 UC, and 105 nonIBD. These 388 samples form the basis for all downstream analysis. 

### MaAsLin2 Association Testing (`04_run_maaslin2.R`)

MaAsLin2 (Multivariate Association Discovery in Population-Scale Meta-Omics Studies) fits a linear mixed model for each feature against each metadata variable. It tests every single biological feature, such as a specific type of bacteria, against every piece of metadata (age, diet, or whether the patient has IBD). Rather than attempting to fit all 122 species and 66, 717 metabolites into a single model (which would lack sufficient statistical power given the sample size of 388), I ran MaAsLin2 in two separate passes: 

**Run 1:** Microbiome species as features, IBD diagnosis as the fixed effect. *This identifies which bacterial species are significantly enriched or depleted in CD and UC patients compared to nonIBD controls.*

**Run 2:** Top 500 most variable metabolites as features, IBD diagnosis as the fixed effect. *This identifies which metabolites are significantly elevated or reduced in IBD patients.* 

In both runs, participant ID was included as a random effect to account for repeated measures (the same patient sampled multiple times over the year). The reference level for diagnosis was set to `nonIBD` to have a baseline of comparison. All data was pre-normalized, so normalization and transformation were set to `NONE` to prevent a second cleaning that could skew results. Significance was determined using a false discovery rate (FDR) of q < 0.25, the conventional benchmark for MaAsLin2 analysis. Simply, out of all the microbes found that are deemed interesting, 25% might actually be insignificant, but the remaining 75% will have real significance. 

### Sparse CCA via mixOmics (`05_run_mixomics.R`)

Sparse Canonical Correlation Analysis (sparse CCA) identifies groups of features from two datasets that covary across samples. Rather than testing one feature at a time, it asks: *which linear combination of bacterial species, when compared to which linear combination of metabolites, produces the maximum correlation across patients?* The "sparse" constraint forces most weights to zero, so instead of all 122 species contributing weakly, only the most informative ones are selected per component. 

I used the `spls()` function from mixOmics with `mode = "canonical"`, which implements sparse CCA. Parameters were set to `ncomp = 2` (two components), `keepX = c(15, 15)` (15 species selected per component), and `keepY = c(50, 50)` (50 metabolites selected per component). The input matrices were the same 388-sample microbiome CLR matrix (X) and the top 500 most variable log2-transformed metabolites (Y). 

---

## Results

### MaAsLin2: IBD-Associated Microbiome Shifts

Of the 244 species-diagnosis association sets (122 species x 2 comparisons), 23 were significant at q < 0.25. The results show a clear biological pattern: IBD patients have a depletion of bacteria known to produce anti-inflammatory short-chain fatty acids, alongside an enrichment of bacteria associated with mucosal inflammation. 

**Species DEPLETED in IBD (health-associated):**

- ***Alistipes putredinis*** - Produces short-chain fatty acids; anti-inflammatory
- ***Alistipes shahii*** - Associated with gut barrier protection
- ***Oscillibacter* sp 57_20** - Produces valerate; anti-inflammatory
- ***Barnesiella intestinihominis*** - Reduced in multiple IBD studies

**Species enriched in IBD:**

- ***Flavonifractor plautii*** - Associated with IBD in multiple cohorts
- ***Ruminococcus gnavus*** - Mucus-degrading; elevated in Crohn's disease
- ***Clostridium symbiosum*** - Pro-inflammatory; associated with dysbiosis
- ***Erysipelatoclostridium ramosum*** - Elevated in gut inflammation
- ***Clostridium clostridioforme*** - Associated with mucosal inflammation

### MaAsLin2: IBD-Associated Metabolomics Shifts

Of the 1,000 metabolite-diagnosis tests (500 metabolites x 2 comparisons), 559 were significant at q < 0.25, a very high proportion that reflects the known sensitivity of metabolomics to gut inflammatory states. The top-ranked metabolites had q-values as low as 8.1 x 10^-10, indicating extremely strong associations and a very low chance that these findings were chance. 

The majority of IBD-associated metabolite features were detected in HILIC positive mode (HILp), a mass spectrometry method particularly sensitive to picking up lipids, bile acids, and amino acid derivatives. This is consistent with published findings showing that bile acid profiles are severely disrupted in IBD, with secondary bile acids (which require microbial processing) being dramatically reduced. 

### Sparse CCA: Co-varying Microbiome-Metabolite Patterns

The sparse CCA identified two components of co-variation between the microbiome and metabolomics datasets: 

- **Component 1:** X-Y correlation = 0.733 (strong)
- **Component 2:** X-Y correlation = 0.617 (moderate)

Component 1 captured the dominant axis of microbe-metabolome co-variation, which aligned strongly with the IBD vs healthy distinction. This component revealed that changes in the gut bacteria (the microbiome) and changes in the gut chemicals (the metabolome) are moving in sync with each other. The sample plot revealed an important asymmetry: patients separated more clearly by diagnosis in the metabolomics space (Block Y) than in the microbiome space (Block X). This suggests that while both data are affected in IBD, the chemical environment of the gut reflects the disease state more consistently than bacterial composition alone. This finding is consistent with the HMP2 paper (Lloyd-Price et al., Nature 2019). 

**Species driving Component 1 (IBD-enriched direction):**
*Clostridium symbiosum*, *Erysipelatoclostridium ramosum*, *Clostridium clostridioforme*, *Flavonifractor plautii*, *Ruminococcus gnavus*

**Species driving Component 1 (health-associated direction):**
*Alistipes putredinis*, *Alistipes shahii*, *Oscillibacter* sp 57_20, *Barnesiella intestinihominis*

The correlation heatmap (Figure 5) revealed a clear reciprocal structure: the IBD-enriched bacteria are positively correlated with one cluster of metabolites and negatively correlated with another. The health-associated bacteria show the exact opposite pattern. This reciprocal signature, where a loss of protective bacteria coincides with a gain of inflammatory ones, and these two groups have opposing metabolite profiles, represents the main multi-omics finding of this analysis.  

The strongest positive metabolite signal co-elevated with IBD-associated bacteria was a HILIC positive mode feature (HILp_QI20096). Consistent with a lipid or modified bile acid class, this specific chemical was elevated during inflammation and increased alongside the proliferation of harmful gut bacteria. The dominant negative metabolite cluster (HILp_QI22835 and related features) represents a class of compounds depleted in IBD and decreases alongside the loss of healthy gut bacteria.

### Cross-Validation Between Methods

A key quality check of any analysis is consistency across methods. In this project, the same microbial species appeared as significant in both MaAsLin2 (tests each microbe one by one against disease) and sparse CCA (looks at how groups of microbes and chemicals change together):

- *Flavonifractor plautii* - significant in both MaAsLin2 and Component 1
- *Alistipes putredinis* - significant in both MaAsLin2 and Component 1
- *Ruminococcus gnavus* - significant in both MaAsLin2 and Component 1

This convergence between two independent statistical methods on the same species provides strong evidence that these findings reflect a genuine biological signal rather than statistical artifacts. 

## Limitations

**Metabolite identity:** The metabolomics features in this dataset are identified by their mass-to-charge ratio and chromatographic retention time (physical traits), not confirmed by official chemical names. Features labeled HILp_QI20096, for example, point towards specific chemical classes but would require further testing (spectral matching or chemical standards) to confirm their identity. This is a standard limitation of untargeted metabolomics. 

**Correlation does not equal causation:** While sparse CCA identifies bacteria and metabolites that co-vary together, it cannot establish whether the bacteria produce the metabolites, whether the metabolites recruit the bacteria, or whether both are driven by a third factor (such as diet or host genetics). Establishing causality would require experimental validation. 

**Repeated measures:** The dataset contains multiple samples per patient taken over time. While participant ID was included as a random effect in MaAsLin2 to account for this, the sparse CCA was run treating all samples as independent. A longitudinal analysis is a logical progression of this project. 

**Sample size for integration:** Only 388 of the 1,638 metagenomics samples had paired metabolomics data. This mismatch reflects the real-world challenge of collecting multiple types of lab data (assays) from the same samples, and limits the statistical power to find meaningful links between them. 

---

## Conclusion

This analysis reproduces and extends core findings from the HMP2 flagship study using a fully scripted, publicly reproducible pipeline. The gut microbiome in IBD patients shows consistent depletion of short-chain fatty acid producers (*Alistipes*, *Oscillibacter*) and enrichment of inflammatory species (*Flavonifractor*, *Ruminococcus gnavus*, *Clostridium* spp.). These microbial shifts co-vary with broad metabolic disruption, particularly in the lipid and bile-acid-related features, and the metabolomics data separates IBD from healthy patients more cleanly than the microbiome data alone. 

The convergence of MaAsLin2 and sparse CCA on the same key species demonstrates the value of using multiple complementary methods in multi-omics research. Together, these results support a model in which IBD-associated microbiome dysbiosis and metabolite dysregulation are tightly coupled, and suggest that integrative analyses capturing both data types simultaneously will be more informative than analyzing either alone.

## References

Lloyd-Price, J., Arze, C., Bhatt, A.S. et al. (2019). Multi-omics of the gut microbial ecosystem in inflammatory bowel diseases. *Nature*, 569, 655–662.
https://doi.org/10.1038/s41586-019-1237-9

Mallick, H., Rahnavard, A., McIver, L.J. et al. (2021). Multivariable association discovery in population-scale meta-omics studies. *PLOS Computational Biology*,
17(11), e1009442.

Rohart, F., Gautier, B., Singh, A., & Lê Cao, K.A. (2017). mixOmics: An R package for 'omics feature selection and multiple data integration. *PLOS Computational
Biology*, 13(11), e1005752.
