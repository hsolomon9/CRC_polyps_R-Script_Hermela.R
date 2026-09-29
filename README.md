# Transcriptional Differences in Colorectal Adenomas by Driver Mutation and Location

Analysis code for my MSc Bioinformatics and Computational Genomics dissertation at Queen's University Belfast (Patrick G. Johnston Centre for Cancer Research):

**"Identifying Transcriptional Differences in Pre-Cancerous Colon Polyps (Adenomas) by Mutation and Location"**

The project compares transcriptomic profiles of colorectal adenomas from the INCISE cohort across driver-mutation groups (APC-only, APC+KRAS, BRAF) and polyp location (left vs. right colon), with a focus on pathway enrichment, stemness programmes and transcription factor activity.

---

## Repository contents

```
.
├── CRC_Adenoma_Analysis_Appendix.R   # Full analysis pipeline (~1,100 lines)
├── README.md
└── figures/                          # Output plots (optional – add if uploading)
```

## Analysis overview

The pipeline in `CRC_Adenoma_Analysis_Appendix.R` runs the following steps:

1. **Cohort definition and summary tables** – assigns samples to mutation and location groups and produces cohort characteristic tables with `gtsummary`. The APC+KRAS group is defined to strictly exclude BRAF-mutant samples (n = 162).
2. **Differential expression** – `DESeq2` comparisons between mutation groups and between left and right colon.
3. **Pathway enrichment (dualGSEA)** – location-based and mutation-based comparisons against Hallmark gene sets, with Bonferroni correction.
4. **Single-sample scoring (ssGSEA)** – per-sample pathway scores.
5. **Stemness signature scoring** – four intestinal stem cell signatures (CBC, RSC, proCSC, revCSC), interpreted as two homeostatic/neoplastic pairs along a proliferative-versus-revival axis.
6. **Transcription factor activity inference** – `decoupleR` using the CollecTRI regulatory network.
7. **Overlap analysis** – Venn diagrams comparing significant results across comparisons.

## Requirements

- R (version X.X or later)
- Key packages:

```r
# CRAN
install.packages(c("tidyverse", "gtsummary", "VennDiagram"))

# Bioconductor
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "GSVA", "decoupleR"))
```

> Update this list to match the `library()` calls at the top of the script, and note how dualGSEA was installed.

## Data availability

The analysis uses transcriptomic and clinical data from the **INCISE** cohort. These data are not publicly available and are **not included in this repository**. Access is subject to the INCISE study's data access procedures.

To run the script, place the expression matrix and sample metadata in a local `data/` folder and update the file paths at the top of `CRC_Adenoma_Analysis_Appendix.R`.

## Usage

```r
source("CRC_Adenoma_Analysis_Appendix.R")
```

Or open the script in RStudio and run it section by section.

## Key findings (summary)

- **Location:** EMT enrichment in left-colon adenomas; immune/inflammatory enrichment in right-colon adenomas.
- **KRAS co-mutation:** E2F targets and G2M checkpoint programmes emerge specifically with KRAS co-mutation.
- **Stemness:** proCSC enrichment in the left colon; revCSC enrichment in APC+KRAS right-colon samples.
- **TF activity:** 335 significant TFs; FOXO3, TGIF1 and LRRFIP1 most enriched in APC+KRAS, and GATA2, LMO2, SP1 and SMAD3 in APC-only. Higher SMAD3 activity in APC-only is consistent with KRAS/MAPK antagonism of TGF-β signalling.

Full interpretation is in the dissertation.

## Acknowledgements

Supervised by [academic supervisor] and [postdoc supervisor] (Dunne Research Team), Patrick G. Johnston Centre for Cancer Research, Queen's University Belfast. Data from the INCISE cohort.

## Author

Hermela — MSc Bioinformatics and Computational Genomics, Queen's University Belfast (2026)
