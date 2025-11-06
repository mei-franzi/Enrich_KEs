# Key Event Enrichment Analysis

A Streamlit application for performing enrichment analysis on differentially expressed genes (DEGs) against genes associated with Key Events (KEs) in Adverse Outcome Pathways (AOPs).

## Overview

This tool identifies which Key Events are statistically overrepresented in your dataset using Fisher's exact test with Benjamini-Hochberg correction (FDR < 0.05).

## Features

- Upload your own differential expression results (CSV or Excel)
- Automatic enrichment analysis against AOP Key Events
- Interactive data visualization
- FDR-corrected statistical analysis

## Installation

1. Clone this repository
2. Install required packages:
```bash
pip install streamlit pandas numpy scipy statsmodels openpyxl
```

## Usage

Run the Streamlit app:
```bash
streamlit run app.py
```

## Data Requirements

Your input file should contain the following columns:
- `padj`: Adjusted p-value
- `log2FoldChange`: Log2 fold change
- `human_ensembl_id`: Human Ensembl gene ID

## Data Files

The `data/` directory contains:
- `Genes_to_KEs.txt`: Mapping of genes to Key Events
- `ke_descriptions.csv`: Key Event descriptions and metadata
- `Browder_DEGs.xlsx`: Example dataset


