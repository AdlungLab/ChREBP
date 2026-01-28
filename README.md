# ChREBP

Single-nucleus RNA-sequencing of murine brown adipose tissue under cold challenges

## Citation

If you utilize our data or code, please cite:

Janina Behrens, Tongtong Wang, Christoph Kilian, Anna Worthmann, Mark A. Herman, Joerg Heeren, Lorenz Adlung\*, Ludger Scheja\*: **"Single-nucleus mRNA-sequencing reveals dynamics of lipogenic and thermogenic adipocyte populations in murine brown adipose tissue in response to cold exposure."** *Mol Metab.* 2025 Nov;101:102252. doi: [10.1016/j.molmet.2025.102252](https://doi.org/10.1016/j.molmet.2025.102252)

\*co-corresponding authors

## Overview

This repository contains R scripts to reproduce the main figure panels from our single-nucleus RNA-sequencing (snRNA-seq) study of murine interscapular brown adipose tissue (BAT). We characterized cellular subpopulations and their adaptation to cold exposure using BAT from control and brown adipocyte-specific ChREBP knockout mice maintained at room temperature (22°C) or exposed to acute (24h) or chronic (10 days) cold (6°C).

## Repository Structure

```
ChREBP/
├── data/                           # Supplemental data files
├── AdlungLab_HelperFunctions.R     # Custom helper functions
├── Figure1.R                       # Script for Figure 1 panels
├── Figure2.R                       # Script for Figure 2 panels
├── Figure3.R                       # Script for Figure 3 panels
├── Figure4.R                       # Script for Figure 4 panels
├── LICENSE                         # GPL-3.0 license
└── README.md                       # This file
```

## Requirements

### R packages

The analysis requires the following R packages:

- **Seurat** (v5.1) - Single-cell analysis
- **Harmony** (v1.2.3) - Data integration
- **ggplot2** - Visualization
- **dplyr** - Data manipulation
- **tidyr** - Data tidying
- **patchwork** - Plot composition
- **clusterProfiler** - Pathway analysis
- **progeny** (v1.28.0) - Pathway activity inference
- **liana** (v0.1.14) - Cell-cell communication analysis
- **scVelo** (v0.3.3) - RNA velocity (Python, via reticulate)

### Data files

Before running the scripts, download the processed Seurat object (see Data section below).

## Usage

1. Clone this repository:
```bash
git clone https://github.com/AdlungLab/ChREBP.git
cd ChREBP
```

2. Download the Seurat object from the link below and place it in your working directory.

3. Source the helper functions and run each figure script:
```r
source("AdlungLab_HelperFunctions.R")
source("Figure1.R")  # Generates Figure 1 panels
source("Figure2.R")  # Generates Figure 2 panels
source("Figure3.R")  # Generates Figure 3 panels
source("Figure4.R")  # Generates Figure 4 panels
```

The scripts are designed to be run line-by-line for transparency and reproducibility.

## Data

### Raw sequencing data

Raw snRNA-seq data are available at ArrayExpress:
- [E-MTAB-15819](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-15819)

### Processed Seurat object

The final processed Seurat object containing 36,611 nuclei across 21 cell clusters:
- [ChREBP_SeuratObject_final.rds](https://doi.org/10.25592/uhhfdm.18248)

## Main Findings

- Identification of 21 distinct cell populations in BAT, including 7 brown adipocyte subtypes
- Discovery of lipogenic adipocytes highly expressing ChREBP and DNL enzymes
- Acute cold exposure transiently depletes lipogenic adipocytes, compensated by basal brown adipocytes
- ChREBP is essential for lipogenic adipocyte identity
- Ttc25 identified as a specific marker of lipogenic brown adipocytes and ChREBP target gene
- Wnt-ChREBP axis implicated in maintenance of lipogenic adipocytes

## License

This project is licensed under the GPL-3.0 License - see the [LICENSE](LICENSE) file for details.

## Contact

Please forward any requests or questions to [Dr. Lorenz Adlung](mailto:l.adlung@uke.de)
