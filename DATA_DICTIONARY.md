# Data Dictionary — m6AAtlas

**Dataset:** m6A epitranscriptome (writer/eraser/reader) gene alterations across the TCGA PanCancer Atlas (32 tumour types, n = 10,195).
**Source:** cBioPortal public REST API (v7.0.0, retrieved June 2026). **Genes (17):** METTL3, METTL14, METTL16, WTAP, VIRMA, ZC3H13, RBM15, FTO, ALKBH5, YTHDF1, YTHDF2, YTHDF3, YTHDC1, YTHDC2, IGF2BP1, IGF2BP2, IGF2BP3.

## Conventions (all files)
- **Encoding/delimiter:** UTF-8, comma-delimited, one header row.
- **`*_freq`** values are **fractions in [0, 1]** (×100 for a percentage). **`*_n`** are integer sample counts.
- **`N` (denominator)** = samples in the **intersection of the `_sequenced` and CNA sample lists** for that study (every counted sample has both a mutation and a CNA call).
- **Copy-number classes (cBioPortal discrete / GISTIC2):** `amp` = +2 (high-level amplification), `gain` = +1, `shallow_del` = −1, `deep_del` = −2.
- **`any_cna`** = any of {gain, shallow_del, amp, deep_del} (CN ≠ 0). **`any_alter`** = **union** of mutation and any CNA (a sample with both counts once).

---

## 1. `tcga_pancan_m6a_epitranscriptome_alteration_freqs_long.csv` — master long table (544 rows × 23 cols)
One row per gene × study. Columns: `tumour_acronym`, `studyId`, `studyName`, `geneSymbol`, `entrezGeneId`, `N`, then per-class count+frequency pairs `mut_n`/`mut_freq`, `amp_n`/`amp_freq`, `gain_n`/`gain_freq`, `shallow_del_n`/`shallow_del_freq`, `deep_del_n`/`deep_del_freq`, `any_cna_n`/`any_cna_freq`, `any_alter_n`/`any_alter_freq`, and the union decomposition `mut_only_n`, `cna_only_n`, `mut_and_cna_n`.
**Check:** `any_alter_n == mut_only_n + cna_only_n + mut_and_cna_n` (holds for all rows).

## 2–4. Matrix files (`*_any_alter_freq_matrix.csv`, `*_any_cna_freq_matrix.csv`, `*_mut_freq_matrix.csv`)
17 gene rows × 33 columns = `geneSymbol` index + 32 tumour-type columns; each cell is the frequency for that gene in that tumour (any-alteration / CNA-only / mutation-only respectively).

## 5. `m6a_epitranscriptome_gene_rank_summary.csv` (17 rows × 9 cols)
`rank`, `geneSymbol`, `mean_any_alter_freq`, `mean_mut_freq`, `mean_any_cna_freq`, `mean_amp_freq`, `mean_deep_del_freq`, `max_any_alter_freq`, `n_studies` (=32). Ranked by `mean_any_alter_freq` (desc).

## 6. `m6a_epitranscriptome_top5_tumours_per_gene.csv` (85 rows × 9 cols)
`geneSymbol`, `rank` (1–5), `tumour_acronym`, `N`, `any_alter_freq`, `mut_freq`, `any_cna_freq`, `amp_freq`, `deep_del_freq`.

## 7. `m6a_epitranscriptome_tumour_burden_ranking.csv` (32 rows × 7 cols)
`rank` (1–32), `tumour_acronym`, `mean_any_alter_freq`, `mean_mut_freq`, `mean_any_cna_freq` (means across all 17 genes), `N`, `n_genes` (=17). Ranked by `mean_any_alter_freq` (desc).

## 8. `brca_m6a_epitranscriptome_alteration_frequencies.csv` (17 rows × 15 cols; N = 1052)
Breast-cancer deep-dive: `gene`, `entrez_id`, `mut_n`, `amp_n`, `gain_n`, `shallow_del_n`, `deep_del_n`, `any_cna_n`, `any_alter_n`, `N`, `mut_freq`, `amp_freq`, `deep_del_freq`, `any_cna_freq`, `any_alter_freq`. Uses the same intersection denominator (N = 1052) as the master table — consistent throughout.

## 9. `tcga_pancan_atlas_sample_set_summary.csv` (32 rows × 7 cols)
`studyId`, `cancerTypeId` (OncoTree), `studyName`, `sampleCount_sequenced`, `sampleCount_cna`, `intersection_sequenced_cna_n` (= `N`; sums to 10,195), `ratio_intersection_over_min_seq_cna`.

## 10. `tcga_low_coverage_tumours_for_m6a_epitranscriptome.txt`
Plain-text note listing the seven tumour types with N < 100 (CHOL, DLBC, UCS, KICH, UVM, MESO, ACC).
