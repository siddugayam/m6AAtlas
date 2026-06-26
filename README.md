# m6AAtlas

**Pan-cancer alteration landscape of the m6A epitranscriptome machinery**

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Data: TCGA](https://img.shields.io/badge/Data-TCGA%20PanCancer%20Atlas-green)](https://www.cbioportal.org)
[![R version](https://img.shields.io/badge/R-%3E%3D4.2-blue)](https://www.r-project.org)

---

## Overview

**m6AAtlas** is a reproducible data resource quantifying copy number alteration (CNA) and somatic mutation frequencies for **17 N6-methyladenosine (m6A) epitranscriptome genes** — methyltransferase writers, demethylase erasers and YTH/IGF2BP readers — across **32 TCGA PanCancer Atlas tumour types** (n = 10,195 samples). The m6A machinery is increasingly implicated in tumourigenesis and therapy response.

This repository accompanies the manuscript:

> Gayam PKR, Earny VA, Aranjani JM, Fayaz SM. *Pan-cancer alteration landscape of the m6A epitranscriptome machinery.* [Journal name, 2026]

---

## Key findings

- **ZC3H13** is the most frequently altered gene across cancers (mean **46.2%** any-alteration), and **ZC3H13** is the most frequently *mutated* (mean **2.3%**)
- **Copy number alteration dominates** over somatic mutation across all 17 genes
- **UCS** shows the highest pan-pathway alteration burden (**71.8%**); **THCA** the lowest (**2.7%**)
- In breast cancer (BRCA, n=1,052), **FTO** is the most altered gene (**70.5%**)

---

## Repository structure

```
m6AAtlas/
├── data/                          # All CSV data files
├── notebooks/notebook.ipynb       # Python data retrieval pipeline (cBioPortal API)
├── figures/                       # Generated TIFF figures (300 DPI)
├── m6a_epitranscriptome_figures.R
├── DATA_DICTIONARY.md
├── CITATION.cff
├── requirements.txt
├── README.md
├── LICENSE                        # MIT (code)
└── LICENSE-data.md                # CC BY 4.0 (data)
```

---

## Genes analysed

| Category | Genes |
|---|---|
| Writer | METTL3, METTL14, METTL16, WTAP, VIRMA, ZC3H13, RBM15 |
| Eraser | FTO, ALKBH5 |
| Reader | YTHDF1, YTHDF2, YTHDF3, YTHDC1, YTHDC2, IGF2BP1, IGF2BP2, IGF2BP3 |

---

## Reproducing the figures

### Requirements
- R ≥ 4.2
- Packages: `tidyverse`, `pheatmap`, `RColorBrewer`, `ggpubr`, `scales`, `cowplot`, `viridis`, `ggtext`, `patchwork`

### Steps
```r
source("m6a_epitranscriptome_figures.R")   # Output: 5 TIFF files in figures/
```
All packages are automatically installed if not already present.

---

## Reproducing the data retrieval

The Python notebook (`notebooks/notebook.ipynb`) retrieves data directly from the **cBioPortal public API** and is fully reproducible.

```bash
pip install -r requirements.txt
cd notebooks && jupyter notebook notebook.ipynb
```
> No login or API key required. cBioPortal data is publicly available. Retrieved June 2026 (cBioPortal v7.0.0).

---

## Data description

| File | Description | Rows | Columns |
|---|---|---|---|
| `*_alteration_freqs_long.csv` | Master long-format table | 544 | 23 (see `DATA_DICTIONARY.md`) |
| `*_any_alter_freq_matrix.csv` | Gene × tumour matrix (any alteration) | 17 | 32 |
| `*_any_cna_freq_matrix.csv` | Gene × tumour matrix (CNA only) | 17 | 32 |
| `*_mut_freq_matrix.csv` | Gene × tumour matrix (mutation only) | 17 | 32 |
| `m6a_epitranscriptome_gene_rank_summary.csv` | Gene ranking across 32 studies | 17 | 9 |
| `m6a_epitranscriptome_top5_tumours_per_gene.csv` | Top 5 tumours per gene | 85 | 9 |
| `m6a_epitranscriptome_tumour_burden_ranking.csv` | Tumour burden ranking (17 genes) | 32 | 7 |
| `brca_m6a_epitranscriptome_*.csv` | BRCA deep-dive (n=1,052) | 17 | 15 |
| `tcga_pancan_atlas_sample_set_summary.csv` | Study metadata and sample counts | 32 | 7 |

---

## Tumour type abbreviations

| Code | Cancer type | N |
|---|---|---|
| BRCA | Breast Invasive Carcinoma | 1,052 |
| COADREAD | Colorectal Adenocarcinoma | 532 |
| LGG | Brain Lower Grade Glioma | 511 |
| LUAD | Lung Adenocarcinoma | 511 |
| OV | Ovarian Serous Cystadenocarcinoma | 511 |
| UCEC | Uterine Corpus Endometrial Carcinoma | 511 |
| HNSC | Head and Neck Squamous Cell Carcinoma | 509 |
| PRAD | Prostate Adenocarcinoma | 489 |
| THCA | Thyroid Carcinoma | 487 |
| LUSC | Lung Squamous Cell Carcinoma | 484 |
| *...and 22 more* | | |

Full list in `data/tcga_pancan_atlas_sample_set_summary.csv`

---

## Citation

```bibtex
@article{gayam2026m6aepitranscriptome,
  title   = {Pan-cancer alteration landscape of the m6A epitranscriptome machinery},
  author  = {Gayam, Prasanna Kumar Reddy and Earny, Venkat Abhiram and Aranjani, Jesil Mathew and SM, Fayaz},
  journal = {[Journal name]},
  year    = {2026}
}
```
Please also cite the cBioPortal:
> Cerami E, et al. The cBio Cancer Genomics Portal. *Cancer Discov.* 2012;2(5):401–404.

---

## Data availability

All underlying data were obtained from publicly available TCGA PanCancer Atlas studies via the cBioPortal API. No individual patient data are included.

---

## License

Code is licensed under the **MIT License** ([LICENSE](LICENSE)); data files are licensed under **CC BY 4.0** ([LICENSE-data](LICENSE-data.md)). Every column is documented in [`DATA_DICTIONARY.md`](DATA_DICTIONARY.md).

---

## Contact

**Gayam Prasanna Kumar Reddy**
Manipal Academy of Higher Education, India
siddugayam1994@gmail.com

*Issues and pull requests welcome.*
