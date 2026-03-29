# RNA-seq Analysis — Gene Knockdown in Human Cells

R-based differential expression visualization and pathway enrichment pipeline for siRNA-mediated gene knockdown RNA-seq in human cells.

> **Note:** Scripts only — raw data, results, and gene identity are not included. The underlying dataset is tied to a manuscript in preparation.

---

## Project Overview

This repository contains visualization and enrichment scripts for an siRNA knockdown RNA-seq experiment in human cells, comparing knockdown vs. control conditions.

**Experimental design:**
- **Comparison:** siKnockdown vs. siControl
- **Replicates:** 3 per condition
- **Organism:** *Homo sapiens*
- **Upstream pipeline:** DESeq2 (differential expression), TPM values for visualization
- **Significance filter:** padj < 0.05 (BH correction)
- **FC thresholds:** Two compared — log2FC > 1.0 (2-fold, stringent) and log2FC > 0.5 (1.41-fold, permissive)

---

## Repository Structure

```
RNA-seq-knockdown/
│
├── README.md
├── 01_knockdown_analysis.R    # DEG classification, volcano plots, heatmaps
├── 02_GO_enrichment.R         # GO enrichment (BP, MF, CC) — 3 plot types
└── 03_pathway_enrichment.R    # KEGG and Reactome pathway enrichment
```

Scripts are designed to run sequentially. Script 01 generates the processed CSV files that scripts 02 and 03 depend on.

---

## Analytical Workflow

### Script 01 — Differential Expression Visualization

**Steps:**
1. Classify each gene as Upregulated, Downregulated, or Not Significant at two FC thresholds
2. Generate volcano plots at each threshold — labeled variant (top 20 genes by padj) and clean variant
3. Generate heatmaps of top 50 DEGs (log2 TPM → z-score, clustered rows, fixed column order)
4. Export gene lists (all significant, upregulated only, downregulated only)

**Outputs per threshold:**
- Processed data CSV
- Volcano plots (labeled and clean) — PDF
- Heatmap — PDF
- Significant gene lists — CSV

---

### Script 02 — GO Enrichment Analysis

**Method:** Over-representation analysis (ORA) using `clusterProfiler::enrichGO` with the full detected gene list as background. BH FDR correction applied; padj < 0.05 and qvalue < 0.2 required.

**Three plot types per ontology (BP, MF, CC):**
- **GeneRatio** — standard dotplot; proportion of query genes in each GO term
- **Fold Enrichment** — GeneRatio / BgRatio; values > 1 indicate over-representation
- **Z-score** — (observed − expected) / √expected; higher = stronger enrichment signal

Run for All DEGs, Upregulated only, and Downregulated only, at both FC thresholds.

**Outputs per ontology × direction × threshold:**
- GeneRatio dotplot — PDF
- Fold Enrichment dotplot — PDF
- Z-score dotplot — PDF
- Full enrichment results table — CSV

---

### Script 03 — Pathway Enrichment (KEGG and Reactome)

**Method:** ORA using `clusterProfiler::enrichKEGG` (human, organism = "hsa") and `ReactomePA::enrichPathway`. BH FDR correction; padj < 0.05 and qvalue < 0.2.

**KEGG outputs per direction × threshold:**
- GeneRatio dotplot
- Barplot
- Fold Enrichment dotplot
- Gene-pathway network (cnetplot, only generated when ≥5 significant pathways)

**Reactome outputs per direction × threshold:**
- GeneRatio dotplot
- Fold Enrichment dotplot

All analyses run for All DEGs, Upregulated only, and Downregulated only, at both FC thresholds.

---

## Statistical Notes

- Differential expression testing was performed upstream (DESeq2); this pipeline performs downstream visualization only
- All enrichment analyses use **p-value-only filtering** (padj < 0.05) without an additional fold-change threshold on the enrichment result itself
- Background gene set = all genes detected in the RNA-seq experiment, not the human genome

---

## Requirements

**R version:** ≥ 4.0

**CRAN packages:**
```r
install.packages(c(
  "tidyverse", "readxl", "ggrepel",
  "RColorBrewer", "pheatmap"
))
```

**Bioconductor packages** (install once):
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "ComplexHeatmap",
  "clusterProfiler",
  "org.Hs.eg.db",
  "enrichplot",
  "ReactomePA"
))
```

---

## Usage

```r
setwd("path/to/your/data")

source("01_knockdown_analysis.R")    # generates processed CSVs
source("02_GO_enrichment.R")         # requires processed CSVs from 01
source("03_pathway_enrichment.R")    # requires processed CSVs from 01
```

---

## Author

**Chia-Lung Chuang**
Postdoctoral Research Associate
Department of Developmental Neurobiology
St. Jude Children's Research Hospital
