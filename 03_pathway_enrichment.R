# ============================================================
# Project: A Knockdown RNA-seq Analysis
# Script:  03_pathway_enrichment.R
# Purpose: KEGG and Reactome pathway enrichment analysis for
#          DEGs from siA vs. siControl RNA-seq, at two FC
#          thresholds. Generates GeneRatio dotplots, fold
#          enrichment dotplots, barplots, and gene-pathway
#          network plots (KEGG only, if ≥5 pathways).
# Input:   A_knockdown_processed_data_fc1.0.csv
#          A_knockdown_processed_data_fc0.5.csv
#          (output of 01_A_knockdown_analysis.R)
# Output:  KEGG_*.pdf, KEGG_table_*.csv,
#          Reactome_*.pdf, Reactome_table_*.csv
# Author:  Chia-Lung Chuang
# Updated: 2026-03
# NOTE:    Scripts only — raw data not included. Underlying
#          dataset is tied to a manuscript in preparation.
# Statistical approach: ORA via clusterProfiler::enrichKEGG
#   and ReactomePA::enrichPathway; BH FDR correction;
#   padj < 0.05; background = all detected genes.
# ============================================================

library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ReactomePA)

# NOTE: Set your working directory to the folder containing the input CSVs
# setwd("path/to/your/data")

# ============================================================
# Shared helpers
# ============================================================

parse_enrichment_metrics <- function(df) {
  df %>%
    rowwise() %>%
    mutate(
      gene_in_pathway = as.numeric(strsplit(GeneRatio, "/")[[1]][1]),
      gene_total      = as.numeric(strsplit(GeneRatio, "/")[[1]][2]),
      bg_in_pathway   = as.numeric(strsplit(BgRatio,   "/")[[1]][1]),
      bg_total        = as.numeric(strsplit(BgRatio,   "/")[[1]][2]),
      GeneRatio_num   = gene_in_pathway / gene_total,
      BgRatio_num     = bg_in_pathway   / bg_total,
      FoldEnrichment  = GeneRatio_num / BgRatio_num,
      expected        = gene_total * (bg_in_pathway / bg_total),
      zscore          = (gene_in_pathway - expected) / sqrt(expected)
    ) %>%
    ungroup()
}

base_theme <- list(
  theme_bw(),
  theme(plot.title   = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.text    = element_text(size = 9),
        axis.title   = element_text(size = 10, face = "bold"),
        legend.title = element_text(size = 9))
)

get_gene_sets <- function(data_proc, direction) {
  switch(direction,
         "UP"   = list(genes = data_proc %>% filter(significant == "Upregulated")   %>% pull(Symbol),
                       label = "Upregulated"),
         "DOWN" = list(genes = data_proc %>% filter(significant == "Downregulated") %>% pull(Symbol),
                       label = "Downregulated"),
         list(genes = data_proc %>% filter(significant != "Not Significant") %>% pull(Symbol),
              label = "All DEGs")
  )
}

# ============================================================
# KEGG Enrichment
# ============================================================

run_kegg_enrichment <- function(data_file, suffix, direction = "ALL") {
  data_proc <- read.csv(data_file)
  gs        <- get_gene_sets(data_proc, direction)
  
  cat("\n→ KEGG:", gs$label, "in", suffix, ":", length(gs$genes), "genes\n")
  if (length(gs$genes) < 3) { cat("  Too few genes — skipping\n"); return(NULL) }
  
  genes_entrez <- bitr(gs$genes, fromType = "SYMBOL",
                       toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  bg_entrez    <- bitr(data_proc %>% pull(Symbol),
                       fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  cat("  Mapped:", nrow(genes_entrez), "genes\n")
  
  kegg <- tryCatch(
    enrichKEGG(gene         = genes_entrez$ENTREZID,
               universe     = bg_entrez$ENTREZID,
               organism     = "hsa",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.2),
    error = function(e) { cat("  KEGG error:", e$message, "\n"); NULL }
  )
  
  if (is.null(kegg) || nrow(as.data.frame(kegg)) == 0) {
    cat("  No significant KEGG pathways\n"); return(NULL)
  }
  
  n_pathways <- nrow(as.data.frame(kegg))
  cat("  Found", n_pathways, "KEGG pathways\n")
  
  kegg <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  write.csv(as.data.frame(kegg),
            paste0("KEGG_table_", direction, "_", suffix, ".csv"), row.names = FALSE)
  
  top <- parse_enrichment_metrics(as.data.frame(kegg)) %>%
    arrange(p.adjust) %>% head(20)
  
  # GeneRatio dotplot
  p1 <- dotplot(kegg, showCategory = 20, font.size = 10,
                title = paste0("KEGG — GeneRatio (", gs$label, " ", suffix, ")")) +
    base_theme
  ggsave(paste0("KEGG_GeneRatio_", direction, "_", suffix, ".pdf"),
         p1, width = 12, height = 10, dpi = 600)
  
  # Barplot
  p2 <- barplot(kegg, showCategory = 20,
                title = paste0("KEGG — Barplot (", gs$label, " ", suffix, ")")) +
    base_theme
  ggsave(paste0("KEGG_Barplot_", direction, "_", suffix, ".pdf"),
         p2, width = 12, height = 8, dpi = 600)
  
  # Fold Enrichment
  p3 <- ggplot(top, aes(x = FoldEnrichment,
                        y = reorder(Description, FoldEnrichment))) +
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_color_gradient(low = "#D73027", high = "#4575B4",
                         name = "Adjusted\nP-value") +
    scale_size_continuous(name = "Gene\nCount", range = c(3, 10)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    base_theme +
    labs(title = paste0("KEGG — Fold Enrichment (", gs$label, " ", suffix, ")"),
         x = "Fold Enrichment", y = "")
  ggsave(paste0("KEGG_FoldEnrichment_", direction, "_", suffix, ".pdf"),
         p3, width = 12, height = 10, dpi = 600)
  
  # Gene-pathway network (only when ≥5 pathways to avoid sparse graphs)
  if (n_pathways >= 5) {
    tryCatch({
      p4 <- cnetplot(kegg, showCategory = 10,
                     node_label = "all",
                     cex_label_gene = 0.5, cex_label_category = 0.8) +
        labs(title = paste0("KEGG Network (", gs$label, " ", suffix, ")"))
      ggsave(paste0("KEGG_Network_", direction, "_", suffix, ".pdf"),
             p4, width = 14, height = 12, dpi = 600)
    }, error = function(e) cat("  Network plot failed\n"))
  }
  
  cat("  Saved 3-4 KEGG plots\n")
  return(kegg)
}

# ============================================================
# Reactome Enrichment
# ============================================================

run_reactome_enrichment <- function(data_file, suffix, direction = "ALL") {
  data_proc <- read.csv(data_file)
  gs        <- get_gene_sets(data_proc, direction)
  
  cat("\n→ Reactome:", gs$label, "in", suffix, ":", length(gs$genes), "genes\n")
  if (length(gs$genes) < 3) { cat("  Too few genes — skipping\n"); return(NULL) }
  
  genes_entrez <- bitr(gs$genes, fromType = "SYMBOL",
                       toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  cat("  Mapped:", nrow(genes_entrez), "genes\n")
  
  reactome <- tryCatch(
    enrichPathway(gene         = genes_entrez$ENTREZID,
                  organism     = "human",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable     = TRUE),
    error = function(e) { cat("  Reactome error:", e$message, "\n"); NULL }
  )
  
  if (is.null(reactome) || nrow(as.data.frame(reactome)) == 0) {
    cat("  No significant Reactome pathways\n"); return(NULL)
  }
  
  cat("  Found", nrow(as.data.frame(reactome)), "Reactome pathways\n")
  
  write.csv(as.data.frame(reactome),
            paste0("Reactome_table_", direction, "_", suffix, ".csv"), row.names = FALSE)
  
  top <- parse_enrichment_metrics(as.data.frame(reactome)) %>%
    arrange(p.adjust) %>% head(20)
  
  # GeneRatio dotplot
  p1 <- dotplot(reactome, showCategory = 20, font.size = 10,
                title = paste0("Reactome — GeneRatio (", gs$label, " ", suffix, ")")) +
    base_theme
  ggsave(paste0("Reactome_GeneRatio_", direction, "_", suffix, ".pdf"),
         p1, width = 12, height = 10, dpi = 600)
  
  # Fold Enrichment
  p2 <- ggplot(top, aes(x = FoldEnrichment,
                        y = reorder(Description, FoldEnrichment))) +
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_color_gradient(low = "#D73027", high = "#4575B4",
                         name = "Adjusted\nP-value") +
    scale_size_continuous(name = "Gene\nCount", range = c(3, 10)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    base_theme +
    labs(title = paste0("Reactome — Fold Enrichment (", gs$label, " ", suffix, ")"),
         x = "Fold Enrichment", y = "")
  ggsave(paste0("Reactome_FoldEnrichment_", direction, "_", suffix, ".pdf"),
         p2, width = 12, height = 10, dpi = 600)
  
  cat("  Saved 2 Reactome plots\n")
  return(reactome)
}

# ============================================================
# Run all combinations
# ============================================================

for (sfx in c("fc1.0", "fc0.5")) {
  f <- paste0("A_knockdown_processed_data_", sfx, ".csv")
  cat("\n===== KEGG:", sfx, "=====\n")
  for (dir in c("ALL", "UP", "DOWN")) run_kegg_enrichment(f, sfx, dir)
  
  cat("\n===== Reactome:", sfx, "=====\n")
  for (dir in c("ALL", "UP", "DOWN")) run_reactome_enrichment(f, sfx, dir)
}