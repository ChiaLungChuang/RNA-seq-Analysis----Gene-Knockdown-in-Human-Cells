# ============================================================
# Project: A Knockdown RNA-seq Analysis
# Script:  01_A_knockdown_analysis.R
# Purpose: Differential expression visualization for siA
#          vs. siControl RNA-seq. Generates volcano plots and
#          heatmaps at two FC thresholds (log2FC > 1.0 and
#          log2FC > 0.5). Exports processed gene lists.
# Input:   Fulltable_all.xlsx — full DEG table with TPM values
#          (padj and log2FC from upstream DESeq2 analysis)
# Output:  Processed CSVs, volcano PDFs, heatmap PDFs,
#          significant gene list CSVs
# Author:  Chia-Lung Chuang
# Updated: 2026-03
# NOTE:    Scripts only — raw data not included. Underlying
#          dataset is tied to a manuscript in preparation.
# Dependencies: BiocManager packages — install once:
#   BiocManager::install(c("ComplexHeatmap","clusterProfiler",
#                          "org.Hs.eg.db","enrichplot"))
# ============================================================

library(tidyverse)
library(readxl)
library(ggrepel)
library(ComplexHeatmap)
library(RColorBrewer)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(pheatmap)
library(circlize)
library(grid)

# NOTE: Set your working directory to the folder containing Fulltable_all.xlsx
# setwd("path/to/your/data")

# ============================================================
# Shared color scheme
# ============================================================

sig_colors <- c(
  "Upregulated"    = "#D73027",
  "Downregulated"  = "#4575B4",
  "Not Significant"= "gray70"
)

# ============================================================
# 1. Load Data
# ============================================================

data <- read_xlsx("Fulltable_all.xlsx")

cat("Data dimensions:", nrow(data), "genes x", ncol(data), "columns\n")
cat("Comparison: siA (n=3) vs. siControl (n=3)\n\n")

# Confirm knockdown efficiency
A_row <- data %>% filter(Symbol == "A")
if (nrow(A_row) > 0) {
  cat("A knockdown efficiency:",
      round(2^abs(A_row$log2FC_siA_vs_siControl), 2), "fold\n")
  cat("  log2FC:", round(A_row$log2FC_siA_vs_siControl, 3),
      " | padj:", format(A_row$padj_siA_vs_siControl, scientific = TRUE), "\n\n")
}

# ============================================================
# 2. Classify DEGs at Two FC Thresholds
# ============================================================

# Significance criteria: padj < 0.05 AND |log2FC| > fc_threshold
# Two thresholds are compared:
#   fc1.0 → 2-fold change (stringent)
#   fc0.5 → 1.41-fold change (permissive)

process_with_threshold <- function(data, fc_threshold, suffix) {
  data_proc <- data %>%
    mutate(
      significant = case_when(
        padj_siA_vs_siControl < 0.05 &
          log2FC_siA_vs_siControl >  fc_threshold ~ "Upregulated",
        padj_siA_vs_siControl < 0.05 &
          log2FC_siA_vs_siControl < -fc_threshold ~ "Downregulated",
        TRUE ~ "Not Significant"
      ),
      neg_log10_padj = -log10(padj_siA_vs_siControl)
    )
  
  stats <- data_proc %>%
    count(significant) %>%
    pivot_wider(names_from = significant, values_from = n, values_fill = 0)
  
  cat("Threshold log2FC >", fc_threshold, "(", round(2^fc_threshold, 2), "x ):\n")
  cat("  Upregulated:", stats$Upregulated,
      "| Downregulated:", stats$Downregulated,
      "| Not Significant:", stats$`Not Significant`, "\n\n")
  
  write.csv(data_proc,
            paste0("A_knockdown_processed_data_", suffix, ".csv"),
            row.names = FALSE)
  
  list(data = data_proc, stats = stats, fc_threshold = fc_threshold, suffix = suffix)
}

result_fc1   <- process_with_threshold(data, 1.0, "fc1.0")
result_fc0.5 <- process_with_threshold(data, 0.5, "fc0.5")

# ============================================================
# 3. Volcano Plots
# ============================================================

# Two variants per threshold: labeled (top 20 by padj) and clean

create_volcano_plots <- function(result_obj) {
  d   <- result_obj$data
  thr <- result_obj$fc_threshold
  sfx <- result_obj$suffix
  
  top_genes <- d %>%
    filter(significant != "Not Significant") %>%
    arrange(padj_siA_vs_siControl) %>%
    head(20) %>%
    pull(Symbol)
  
  d <- d %>% mutate(label_top20 = ifelse(Symbol %in% top_genes, Symbol, ""))
  
  base_plot <- ggplot(d, aes(x = log2FC_siA_vs_siControl, y = neg_log10_padj)) +
    geom_point(aes(color = significant, alpha = significant), size = 2) +
    scale_color_manual(values = sig_colors) +
    scale_alpha_manual(values = c("Upregulated" = 0.6,
                                  "Downregulated" = 0.6,
                                  "Not Significant" = 0.3)) +
    geom_vline(xintercept = c(-thr, thr), linetype = "dashed", color = "gray40") +
    geom_hline(yintercept = -log10(0.05),  linetype = "dashed", color = "gray40") +
    theme_minimal() +
    theme(legend.position  = "right",
          plot.title        = element_text(size = 14, face = "bold", hjust = 0.5),
          axis.title        = element_text(size = 12, face = "bold")) +
    labs(title   = paste0("A Knockdown vs Control (",
                          round(2^thr, 2), "x fold change)"),
         x       = "log2 Fold Change",
         y       = "-log10(adjusted P-value)",
         color   = "Regulation",
         caption = paste0("Significance: padj < 0.05 & |log2FC| > ", thr))
  
  # Labeled version
  p_labeled <- base_plot +
    geom_text_repel(aes(label = label_top20), max.overlaps = Inf,
                    size = 3, box.padding = 0.5, segment.color = "grey50", force = 2)
  ggsave(paste0("volcano_A_top20_", sfx, ".pdf"),
         p_labeled, width = 10, height = 8, dpi = 600)
  
  # Clean version
  ggsave(paste0("volcano_A_clean_", sfx, ".pdf"),
         base_plot, width = 10, height = 8, dpi = 600)
  
  cat("Saved volcano plots:", sfx, "\n")
}

create_volcano_plots(result_fc1)
create_volcano_plots(result_fc0.5)

# ============================================================
# 4. Heatmaps — Top DEGs (ComplexHeatmap, z-score of log2 TPM)
# ============================================================

# Shows top 50 DEGs (ranked by padj) across all 6 samples.
# Columns are NOT clustered to preserve siA | siControl order.

create_heatmap <- function(data_file, suffix, n_genes = 50) {
  data_proc <- read.csv(data_file)
  
  top_degs <- data_proc %>%
    filter(significant != "Not Significant") %>%
    arrange(padj_siA_vs_siControl) %>%
    head(n_genes)
  
  if (nrow(top_degs) < 2) {
    cat("Not enough DEGs for heatmap:", suffix, "\n")
    return(NULL)
  }
  
  # Extract TPM columns in fixed sample order
  tpm_cols <- c(
    "TPM_siA_1(siA)", "TPM_siA_2(siA)", "TPM_siA_3(siA)",
    "TPM_siControl_1(siControl)", "TPM_siControl_2(siControl)", "TPM_siControl_3(siControl)"
  )
  
  mat <- data %>%
    filter(Symbol %in% top_degs$Symbol) %>%
    select(Symbol, all_of(tpm_cols)) %>%
    { m <- as.matrix(.[, -1]); rownames(m) <- .$Symbol; m }
  
  # Log2 transform then z-score across rows
  mat_scaled <- t(scale(t(log2(mat + 1))))
  colnames(mat_scaled) <- c("siA_1","siA_2","siA_3",
                            "siCtrl_1","siCtrl_2","siCtrl_3")
  
  col_anno <- HeatmapAnnotation(
    Condition = c(rep("siA", 3), rep("siControl", 3)),
    col = list(Condition = c("siA" = "#D73027", "siControl" = "#4575B4")),
    show_annotation_name = FALSE
  )
  
  col_fun <- colorRamp2(c(-2, 0, 2), c("#4575B4", "white", "#D73027"))
  
  pdf(paste0("heatmap_A_", suffix, ".pdf"), width = 8, height = 11)
  ht <- Heatmap(mat_scaled, name = "z-score", col = col_fun,
                top_annotation  = col_anno,
                cluster_columns = FALSE,
                cluster_rows    = TRUE,
                show_row_names  = TRUE,
                row_names_gp    = gpar(fontsize = 8),
                column_title    = paste0("Top DEGs (", suffix, ")"),
                border = TRUE)
  draw(ht)
  dev.off()
  
  cat("Saved heatmap:", suffix, "\n")
}

create_heatmap("A_knockdown_processed_data_fc1.0.csv", "fc1.0")
create_heatmap("A_knockdown_processed_data_fc0.5.csv", "fc0.5")

# ============================================================
# 5. Export Significant Gene Lists
# ============================================================

export_gene_lists <- function(result_obj) {
  d   <- result_obj$data
  sfx <- result_obj$suffix
  
  sig <- d %>%
    filter(significant != "Not Significant") %>%
    select(Symbol, Description,
           log2FC_siA_vs_siControl,
           padj_siA_vs_siControl,
           pvalue_siA_vs_siControl,
           significant) %>%
    arrange(padj_siA_vs_siControl)
  
  write.csv(sig,
            paste0("significant_genes_", sfx, ".csv"), row.names = FALSE)
  write.csv(sig %>% filter(significant == "Upregulated"),
            paste0("upregulated_genes_", sfx, ".csv"),   row.names = FALSE)
  write.csv(sig %>% filter(significant == "Downregulated"),
            paste0("downregulated_genes_", sfx, ".csv"), row.names = FALSE)
  
  cat("Exported gene lists:", sfx,
      "— Up:", sum(sig$significant == "Upregulated"),
      "| Down:", sum(sig$significant == "Downregulated"), "\n")
}

export_gene_lists(result_fc1)
export_gene_lists(result_fc0.5)