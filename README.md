# DESeq2 Differential Expression Analysis — Shiny App

A browser-based tool for running DESeq2 differential expression analysis on bulk RNA-seq count data, with interactive visualizations and downloadable outputs.

**Live app:** [https://a-sundaresan.shinyapps.io/BulkRNASeq_downstream_processing_after_quantification/](https://a-sundaresan.shinyapps.io/BulkRNASeq_downstream_processing_after_quantification/)

---

## Features

- Upload your own count matrix and sample metadata
- Dynamically select the condition column, control group, and comparison group
- Runs DESeq2 in the browser — no local R installation required
- Interactive results table with download (CSV)
- PCA plot, MA plot, and Volcano plot with download (PDF)

---

## Input File Format

**Count matrix** (`.csv`, `.tsv`, or `.txt`)
- Rows = genes, columns = samples
- First column = gene IDs (used as row names)
- Values should be raw (unnormalized) integer counts

**Metadata** (`.csv`, `.tsv`, or `.txt`)
- Rows = samples, columns = sample attributes
- Row names must match column names in the count matrix
- Must contain at least one column defining the condition/group for each sample

---

## Usage

1. Upload your count matrix and metadata files
2. Select the column in your metadata that defines sample groups
3. Select the control group
4. Select the comparison group
5. Click **Run DESeq2**
6. Explore results across the **DESeq2 Results**, **PCA Plot**, **MA Plot**, and **Volcano Plot** tabs
7. Download any output using the buttons below each panel

---

## Example Datasets

Three example datasets are included in this repository for testing:

| Dataset | Description |
|---------|-------------|
| `GSE183841/` | Bulk RNA-seq with raw counts and TPM |
| `GSE286979/` | Aging bulk RNA-seq (tximport counts) |
| `GSE298096/` | Early infection bulk RNA-seq (raw counts) |

---

## Running Locally

```r
# Install required packages
install.packages(c("shiny", "DT", "shinyjs", "bslib", "ggplot2",
                   "shinycssloaders", "dplyr", "shinyWidgets"))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("DESeq2", "EnhancedVolcano"))

# Launch the app
shiny::runApp("app.R")
```

---

## Dependencies

| Package | Source |
|---------|--------|
| DESeq2 | Bioconductor |
| EnhancedVolcano | Bioconductor |
| shiny, bslib, shinyjs, shinyWidgets, shinycssloaders | CRAN |
| ggplot2, dplyr, DT | CRAN |
