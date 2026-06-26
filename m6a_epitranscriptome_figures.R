################################################################################
##   m6AAtlas — Pan-cancer alteration landscape of the m6A epitranscriptome machinery
##   Figure Generation Script (Publication-ready, 300 DPI TIFF)
##   Author : Gayam Prasanna Kumar Reddy, Manipal Academy of Higher Education
##   GitHub : https://github.com/siddugayam/m6AAtlas
##   Data   : TCGA PanCancer Atlas 2018 via cBioPortal API (v7.0.0, June 2026)
##   License: MIT (code) / CC BY 4.0 (data)
################################################################################

## ── 0. Install and load packages ────────────────────────────────────────────
packages <- c("tidyverse","pheatmap","RColorBrewer","ggpubr","scales","cowplot","viridis","ggtext","patchwork")
installed <- rownames(installed.packages())
to_install <- packages[!packages %in% installed]
if (length(to_install) > 0) install.packages(to_install, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages({
  library(tidyverse); library(pheatmap); library(RColorBrewer); library(ggpubr)
  library(scales); library(cowplot); library(viridis); library(ggtext); library(patchwork)
})
theme_set(theme_cowplot(font_size = 11, font_family = "sans"))

## ── 1. Paths ────────────────────────────────────────────────────────────────
data_dir   <- "data"
output_dir <- "figures"
dir.create(output_dir, showWarnings = FALSE)
SLUG  <- "m6a_epitranscriptome"
NGENES <- 17

## ── 2. Load data ────────────────────────────────────────────────────────────
matrix_any <- read_csv(file.path(data_dir, paste0("tcga_pancan_", SLUG, "_any_alter_freq_matrix.csv")), show_col_types = FALSE)
gene_rank  <- read_csv(file.path(data_dir, paste0(SLUG, "_gene_rank_summary.csv")), show_col_types = FALSE)
tumour_rank<- read_csv(file.path(data_dir, paste0(SLUG, "_tumour_burden_ranking.csv")), show_col_types = FALSE)
brca_df    <- read_csv(file.path(data_dir, paste0("brca_", SLUG, "_alteration_frequencies.csv")), show_col_types = FALSE)

## ── 3. Constants ────────────────────────────────────────────────────────────
LOW_N_TUMOURS <- c("CHOL","DLBC","UCS","KICH","UVM","MESO","ACC")
GENE_ORDER    <- c("METTL3", "METTL14", "METTL16", "WTAP", "VIRMA", "ZC3H13", "RBM15", "FTO", "ALKBH5", "YTHDF1", "YTHDF2", "YTHDF3", "YTHDC1", "YTHDC2", "IGF2BP1", "IGF2BP2", "IGF2BP3")
GENE_CATEGORY <- c("METTL3"="Writer", "METTL14"="Writer", "METTL16"="Writer", "WTAP"="Writer", "VIRMA"="Writer", "ZC3H13"="Writer", "RBM15"="Writer", "FTO"="Eraser", "ALKBH5"="Eraser", "YTHDF1"="Reader", "YTHDF2"="Reader", "YTHDF3"="Reader", "YTHDC1"="Reader", "YTHDC2"="Reader", "IGF2BP1"="Reader", "IGF2BP2"="Reader", "IGF2BP3"="Reader")
CAT_COLORS    <- c("Writer"="#534AB7", "Eraser"="#1D9E75", "Reader"="#D85A30")
COL_MUTATION<-"#A32D2D"; COL_CNA<-"#185FA5"; COL_AMP<-"#3B6D11"; COL_DEL<-"#E24B4A"; COL_CNA_OTHER<-"#85B7EB"
HEATMAP_COLS <- colorRampPalette(c("#FFFFFF","#E6F1FB","#85B7EB","#185FA5","#042C53"))(100)

## ── Figure 1: Pan-cancer heatmap ────────────────────────────────────────────
cat("Generating Figure 1 - pan-cancer heatmap...\n")
mat <- matrix_any %>% column_to_rownames("geneSymbol") %>% as.matrix()
mat <- mat[GENE_ORDER[GENE_ORDER %in% rownames(mat)], ]
col_annot <- data.frame(Sample_size = ifelse(colnames(mat) %in% LOW_N_TUMOURS, "N < 100", "N >= 100"), row.names = colnames(mat))
annot_colours <- list(Sample_size = c("N < 100"="#EF9F27","N >= 100"="#B5D4F4"))
row_annot <- data.frame(Category = unname(GENE_CATEGORY[rownames(mat)]), row.names = rownames(mat))
annot_colours$Category <- CAT_COLORS
tiff(file.path(output_dir, "Fig1_pancan_heatmap.tiff"), width = 12, height = 6, units = "in", res = 300, compression = "lzw")
pheatmap(mat, color = HEATMAP_COLS, breaks = seq(0,1,length.out=101),
  cluster_rows = FALSE, cluster_cols = FALSE,
  annotation_col = col_annot, annotation_row = row_annot, annotation_colors = annot_colours,
  display_numbers = TRUE, number_format = "%.2f", number_color = ifelse(mat > 0.55, "white", "black"),
  fontsize = 8, fontsize_number = 5.5, fontsize_row = 9, fontsize_col = 8, angle_col = 45,
  border_color = "grey92", legend_breaks = c(0,0.25,0.50,0.75,1.0),
  legend_labels = c("0%","25%","50%","75%","100%"), main = "")
dev.off(); cat("  --> figures/Fig1_pancan_heatmap.tiff\n")

## ── Figure 2: Gene ranking (stacked bar) ────────────────────────────────────
cat("Generating Figure 2 - gene ranking stacked bar...\n")
fig2_df <- gene_rank %>% arrange(mean_any_alter_freq) %>%
  mutate(geneSymbol = factor(geneSymbol, levels = geneSymbol),
    mut_pct = mean_mut_freq*100, amp_pct = mean_amp_freq*100, del_pct = mean_deep_del_freq*100,
    cna_other = pmax(0,(mean_any_cna_freq - mean_amp_freq - mean_deep_del_freq))*100) %>%
  select(geneSymbol, mut_pct, amp_pct, del_pct, cna_other) %>%
  pivot_longer(-geneSymbol, names_to="class", values_to="freq") %>%
  mutate(class = factor(class, levels=c("del_pct","amp_pct","cna_other","mut_pct"),
    labels=c("Deep deletion","Amplification","Other CNA","Somatic mutation")))
fig2 <- ggplot(fig2_df, aes(x=freq, y=geneSymbol, fill=class)) + geom_col(width=0.72, position="stack") +
  scale_fill_manual(values=c("Somatic mutation"=COL_MUTATION,"Other CNA"=COL_CNA,"Amplification"=COL_AMP,"Deep deletion"=COL_DEL), name="Alteration class") +
  scale_x_continuous(labels=percent_format(scale=1), expand=expansion(mult=c(0,0.02))) +
  labs(x="Mean alteration frequency (%)", y=NULL, caption="n = 32 TCGA PanCancer Atlas studies") +
  theme(legend.position="bottom", legend.key.size=unit(0.35,"cm"), legend.text=element_text(size=9),
    legend.title=element_text(size=9,face="bold"), axis.text.y=element_text(face="italic",size=10),
    axis.text.x=element_text(size=9), axis.title.x=element_text(size=10), plot.caption=element_text(size=8,colour="grey50"),
    panel.grid.major.x=element_line(colour="grey90",linewidth=0.4), panel.grid.major.y=element_blank()) +
  guides(fill=guide_legend(nrow=2, byrow=TRUE))
ggsave(file.path(output_dir,"Fig2_gene_ranking.tiff"), fig2, width=7, height=5.5, dpi=300, device="tiff", compression="lzw")
cat("  --> figures/Fig2_gene_ranking.tiff\n")

## ── Figure 3: Tumour type burden ranking ────────────────────────────────────
cat("Generating Figure 3 - tumour burden ranking...\n")
fig3_df <- tumour_rank %>% arrange(mean_any_alter_freq) %>%
  mutate(tumour_acronym = factor(tumour_acronym, levels=tumour_acronym),
    low_n = tumour_acronym %in% LOW_N_TUMOURS,
    label = if_else(low_n, paste0(as.character(tumour_acronym)," (!)"), as.character(tumour_acronym)),
    label = factor(label, levels=label),
    cna_pct = pmax(0, mean_any_alter_freq - mean_mut_freq)*100, mut_pct = mean_mut_freq*100) %>%
  select(label, low_n, cna_pct, mut_pct) %>%
  pivot_longer(c(cna_pct, mut_pct), names_to="class", values_to="freq") %>%
  mutate(class = factor(class, levels=c("mut_pct","cna_pct"), labels=c("Somatic mutation","CNA")))
fig3 <- ggplot(fig3_df, aes(x=freq, y=label, fill=class)) + geom_col(width=0.75, position="stack") +
  scale_fill_manual(values=c("Somatic mutation"=COL_CNA_OTHER,"CNA"=COL_CNA), name="Alteration class") +
  scale_x_continuous(labels=percent_format(scale=1), expand=expansion(mult=c(0,0.03))) +
  labs(x=paste0("Mean alteration frequency, ",NGENES," genes (%)"), y=NULL, caption="(!) N < 100 samples") +
  theme(legend.position="bottom", legend.key.size=unit(0.35,"cm"), legend.text=element_text(size=9),
    legend.title=element_text(size=9,face="bold"), axis.text.y=element_text(size=8.5), axis.text.x=element_text(size=9),
    axis.title.x=element_text(size=10), plot.caption=element_text(size=8,colour="grey50"),
    panel.grid.major.x=element_line(colour="grey90",linewidth=0.4), panel.grid.major.y=element_blank())
ggsave(file.path(output_dir,"Fig3_tumour_burden_ranking.tiff"), fig3, width=7, height=9, dpi=300, device="tiff", compression="lzw")
cat("  --> figures/Fig3_tumour_burden_ranking.tiff\n")

## ── Figure 4: BRCA deep-dive ────────────────────────────────────────────────
cat("Generating Figure 4 - BRCA deep-dive...\n")
brca_n <- brca_df$N[1]
fig4_df <- brca_df %>% arrange(any_alter_freq) %>%
  mutate(gene = factor(gene, levels=gene), mut_pct = mut_freq*100, amp_pct = amp_freq*100, del_pct = deep_del_freq*100,
    cna_other = pmax(0, any_cna_freq - amp_freq - deep_del_freq)*100) %>%
  select(gene, mut_pct, amp_pct, del_pct, cna_other) %>%
  pivot_longer(-gene, names_to="class", values_to="freq") %>%
  mutate(class = factor(class, levels=c("del_pct","amp_pct","cna_other","mut_pct"),
    labels=c("Deep deletion","Amplification","CNA gain / shallow deletion","Somatic mutation")))
fig4 <- ggplot(fig4_df, aes(x=freq, y=gene, fill=class)) + geom_col(width=0.72, position="stack") +
  scale_fill_manual(values=c("Somatic mutation"=COL_MUTATION,"CNA gain / shallow deletion"=COL_CNA_OTHER,"Amplification"=COL_AMP,"Deep deletion"=COL_DEL), name="Alteration class") +
  scale_x_continuous(labels=percent_format(scale=1), expand=expansion(mult=c(0,0.02))) +
  labs(x="Alteration frequency (%)", y=NULL, caption=paste0("Breast invasive carcinoma (BRCA) - n = ", format(brca_n, big.mark=","), " samples - TCGA PanCancer Atlas 2018")) +
  theme(legend.position="bottom", legend.key.size=unit(0.35,"cm"), legend.text=element_text(size=9),
    legend.title=element_text(size=9,face="bold"), axis.text.y=element_text(face="italic",size=10), axis.text.x=element_text(size=9),
    axis.title.x=element_text(size=10), plot.caption=element_text(size=8,colour="grey50"),
    panel.grid.major.x=element_line(colour="grey90",linewidth=0.4), panel.grid.major.y=element_blank()) +
  guides(fill=guide_legend(nrow=2, byrow=TRUE))
ggsave(file.path(output_dir,"Fig4_BRCA_deepdive.tiff"), fig4, width=7, height=6, dpi=300, device="tiff", compression="lzw")
cat("  --> figures/Fig4_BRCA_deepdive.tiff\n")

## ── Combined panel ──────────────────────────────────────────────────────────
cat("Generating combined panel figure...\n")
combined <- (fig2 | fig4) / (fig3 + theme(plot.margin = unit(c(0,1,0,0),"cm"))) +
  plot_layout(heights = c(1, 1.6)) +
  plot_annotation(title = "m6AAtlas: Pan-cancer alteration landscape of the m6A epitranscriptome machinery",
    subtitle = paste0("TCGA PanCancer Atlas 2018 - ", NGENES, " genes - 32 tumour types - n = 10,195"),
    caption = "Data: cBioPortal public API - github.com/siddugayam/m6AAtlas",
    theme = theme(plot.title=element_text(size=13,face="bold"), plot.subtitle=element_text(size=10,colour="grey40"), plot.caption=element_text(size=8,colour="grey50")))
ggsave(file.path(output_dir,"Fig_combined_panel.tiff"), combined, width=14, height=14, dpi=300, device="tiff", compression="lzw")
cat("  --> figures/Fig_combined_panel.tiff\n")
cat("\n-- All figures saved to:", output_dir, "--\n")
