# ============================================================
# Project: A Knockdown RNA-seq Analysis
# Script:  02_GO_enrichment.R
# Purpose: GO enrichment analysis (BP, MF, CC) for DEGs from
#          siA vs. siControl RNA-seq, at two FC thresholds.
#          Generates three plot types per ontology:
#          GeneRatio dotplot, Fold Enrichment dotplot, Z-score
#          dotplot. Exports full result tables as CSV.
# Input:   A_knockdown_processed_data_fc1.0.csv
#          A_knockdown_processed_data_fc0.5.csv
#          (output of 01_A_knockdown_analysis.R)
# Output:  GO_GeneRatio_*.pdf, GO_FoldEnrichment_*.pdf,
#          GO_Zscore_*.pdf, GO_table_*.csv
# Author:  Chia-Lung Chuang
# Updated: 2026-03
# NOTE:    Scripts only — raw data not included. Underlying
#          dataset is tied to a manuscript in preparation.
# Statistical approach: ORA via clusterProfiler::enrichGO;
#   BH FDR correction; padj < 0.05; background = all detected
#   genes in the RNA-seq experiment.
# ============================================================

library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# NOTE: Set your working directory to the folder containing the input CSVs
# setwd("path/to/your/data")

# ============================================================
# GO Enrichment — Three Plot Types per Ontology
# ============================================================

# For each direction (ALL / UP / DOWN) and each ontology (BP, MF, CC):
#   Plot 1 — GeneRatio: proportion of query genes in each GO term
#   Plot 2 — Fold Enrichment: GeneRatio / BgRatio (>1 = over-represented)
#   Plot 3 — Z-score: (observed - expected) / sqrt(expected)
#
# Background: all genes detected in the experiment (full table).
# Significance: padj < 0.05, qvalue < 0.2 (BH correction).

run_go_enrichment_enhanced <- function(data_file, suffix, direction = "ALL") {
  data_proc <- read.csv(data_file)
  
  sig_genes <- switch(direction,
                      "UP"   = data_proc %>% filter(significant == "Upregulated")   %>% pull(Symbol),
                      "DOWN" = data_proc %>% filter(significant == "Downregulated") %>% pull(Symbol),
                      data_proc %>% filter(significant != "Not Significant") %>% pull(Symbol)
  )
  direction_label <- switch(direction,
                            "UP" = "Upregulated", "DOWN" = "Downregulated", "All DEGs")
  
  cat("\n→", direction_label, "in", suffix, ":", length(sig_genes), "genes\n")
  if (length(sig_genes) < 3) { cat("  Too few genes — skipping\n"); return(NULL) }
  
  genes_entrez <- bitr(sig_genes, fromType = "SYMBOL",
                       toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  bg_entrez    <- bitr(data_proc %>% pull(Symbol),
                       fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  cat("  Mapped:", nrow(genes_entrez), "genes\n")
  
  ont_names <- c(BP = "Biological Process",
                 MF = "Molecular Function",
                 CC = "Cellular Component")
  
  for (ont in c("BP", "MF", "CC")) {
    cat("    —", ont_names[ont], "...")
    
    ego <- tryCatch(
      enrichGO(gene        = genes_entrez$ENTREZID,
               universe    = bg_entrez$ENTREZID,
               OrgDb       = org.Hs.eg.db,
               ont         = ont,
               pAdjustMethod = "BH",
               pvalueCutoff  = 0.05,
               qvalueCutoff  = 0.2,
               readable      = TRUE),
      error = function(e) NULL
    )
    
    if (is.null(ego) || nrow(as.data.frame(ego)) == 0) { cat(" no terms\n"); next }
    cat("", nrow(as.data.frame(ego)), "terms\n")
    
    write.csv(as.data.frame(ego),
              paste0("GO_table_", ont, "_", direction, "_", suffix, ".csv"),
              row.names = FALSE)
    
    # Compute additional metrics from ratio strings
    ego_df <- as.data.frame(ego) %>%
      rowwise() %>%
      mutate(
        gene_in_term   = as.numeric(strsplit(GeneRatio, "/")[[1]][1]),
        gene_total     = as.numeric(strsplit(GeneRatio, "/")[[1]][2]),
        bg_in_term     = as.numeric(strsplit(BgRatio,   "/")[[1]][1]),
        bg_total       = as.numeric(strsplit(BgRatio,   "/")[[1]][2]),
        GeneRatio_num  = gene_in_term / gene_total,
        BgRatio_num    = bg_in_term   / bg_total,
        FoldEnrichment = GeneRatio_num / BgRatio_num,
        expected       = gene_total * (bg_in_term / bg_total),
        zscore         = (gene_in_term - expected) / sqrt(expected)
      ) %>%
      ungroup() %>%
      arrange(p.adjust) %>%
      head(15)
    
    base_theme <- list(
      theme_bw(),
      theme(plot.title   = element_text(size = 12, face = "bold", hjust = 0.5),
            axis.text    = element_text(size = 9),
            axis.title   = element_text(size = 10, face = "bold"),
            legend.title = element_text(size = 9),
            legend.text  = element_text(size = 8))
    )
    
    # Plot 1: GeneRatio
    p1 <- dotplot(ego, showCategory = 15, font.size = 10,
                  title = paste0(ont_names[ont], " — GeneRatio (",
                                 direction_label, " ", suffix, ")")) +
      base_theme + labs(x = "Gene Ratio")
    ggsave(paste0("GO_GeneRatio_", ont, "_", direction, "_", suffix, ".pdf"),
           p1, width = 10, height = 8, dpi = 600)
    
    # Plot 2: Fold Enrichment
    p2 <- ggplot(ego_df, aes(x = FoldEnrichment,
                             y = reorder(Description, FoldEnrichment))) +
      geom_point(aes(size = Count, color = p.adjust)) +
      scale_color_gradient(low = "#D73027", high = "#4575B4",
                           name = "Adjusted\nP-value") +
      scale_size_continuous(name = "Gene\nCount", range = c(3, 10)) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
      base_theme +
      labs(title = paste0(ont_names[ont], " — Fold Enrichment (",
                          direction_label, " ", suffix, ")"),
           x = "Fold Enrichment", y = "")
    ggsave(paste0("GO_FoldEnrichment_", ont, "_", direction, "_", suffix, ".pdf"),
           p2, width = 10, height = 8, dpi = 600)
    
    # Plot 3: Z-score
    p3 <- ggplot(ego_df, aes(x = zscore,
                             y = reorder(Description, zscore))) +
      geom_point(aes(size = Count, color = p.adjust)) +
      scale_color_gradient(low = "#D73027", high = "#4575B4",
                           name = "Adjusted\nP-value") +
      scale_size_continuous(name = "Gene\nCount", range = c(3, 10)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
      base_theme +
      labs(title = paste0(ont_names[ont], " — Z-score (",
                          direction_label, " ", suffix, ")"),
           x = "Z-score (Enrichment Score)", y = "")
    ggsave(paste0("GO_Zscore_", ont, "_", direction, "_", suffix, ".pdf"),
           p3, width = 10, height = 8, dpi = 600)
    
    cat("      Saved 3 plots for", ont, "\n")
  }
}

# ============================================================
# Run for both thresholds and all directions
# ============================================================

for (sfx in c("fc1.0", "fc0.5")) {
  cat("\n===", sfx, "===\n")
  f <- paste0("A_knockdown_processed_data_", sfx, ".csv")
  run_go_enrichment_enhanced(f, sfx, "ALL")
  run_go_enrichment_enhanced(f, sfx, "UP")
  run_go_enrichment_enhanced(f, sfx, "DOWN")
}