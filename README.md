# Key Event and Functional Enrichment Analysis

A comprehensive Streamlit application for performing enrichment analysis on differentially expressed genes (DEGs) against Key Events (KEs) in Adverse Outcome Pathways (AOPs) and functional pathways (GO:BP, KEGG).

## Overview

**Key Events** are measurable biological events within Adverse Outcome Pathways (AOPs)—multi-scale models that connect molecular initiating events to adverse health outcomes. This tool maps your DEGs to the AOP Key Event database and identifies statistically overrepresented KEs using Fisher's exact test with Benjamini-Hochberg correction.

## Features

### Data Upload & Filtering
- Upload differential expression results (CSV, TSV, or Excel)
- Flexible column mapping for non-standard column names
- Customizable filtering parameters (padj and log2FC cutoffs)
- Interactive, scrollable data tables with live filtering

### Key Event Enrichment
- Statistical enrichment analysis against AOP Key Events
- Fisher's exact test with FDR correction (Benjamini-Hochberg)
- Detailed gene information for each significant KE
- Interactive heatmaps showing log2FC values for overlapping genes

### Functional Enrichment
Three levels of functional enrichment analysis using g:Profiler:
1. **All DEGs**: Enrichment on your complete filtered DEG list
2. **Per KE**: Enrichment on genes associated with each significant KE
3. **All KE Genes**: Enrichment on the union of all genes from significant KEs

Supported databases:
- GO Biological Processes (GO:BP)
- KEGG Pathways

### Visualization
- Interactive Plotly heatmaps for gene expression patterns
- Bar plots for enrichment results
- Expandable sections for detailed exploration
- Scientific notation formatting for p-values

## Installation

1. Clone this repository:
```bash
git clone https://github.com/mei-franzi/Enrich_KEs.git
cd Enrich_KEs
```

2. Install required packages:
```bash
pip install -r requirements.txt
```

Or install manually:
```bash
pip install streamlit pandas numpy scipy statsmodels openpyxl plotly gprofiler-official matplotlib
```

## Usage

Run the Streamlit app:
```bash
streamlit run app.py
```

Then open your browser to the provided local URL (typically http://localhost:8501)

## Data Requirements

### Required Columns
Your input file should contain columns for:
- Adjusted p-value (e.g., `padj`, `adj_pvalue`)
- Log2 fold change (e.g., `log2FoldChange`, `log2FC`)
- Human Ensembl gene ID (e.g., `human_ensembl_id`, `ensembl_id`)

### Optional Columns
- Gene names/symbols (e.g., `gene`, `gene_name`, `symbol`) - required for functional enrichment

**Note:** The app includes a column mapping feature that allows you to map your custom column names to the expected format.

## Data Files

The `data/` directory contains:
- `Genes_to_KEs.txt`: Mapping of genes to Key Events (tab-separated)
- `ke_descriptions.csv`: Key Event descriptions and metadata
- `Browder_DEGs.xlsx`: Example dataset for testing

## Analysis Workflow

1. **Upload Data**: Upload your DEG file and map column names if needed
2. **Filter DEGs**: Adjust padj and log2FC cutoffs to preview filtered genes
3. **Optional Functional Enrichment**: Run GO:BP/KEGG enrichment on your DEGs
4. **KE Enrichment**: Perform Key Event enrichment analysis with customized parameters
5. **Explore Results**: 
   - View enrichment statistics
   - Examine detailed gene lists and heatmaps for each KE
   - Run functional enrichment on KE-associated genes

## Output

### KE Enrichment Results Table
- KE ID and name
- Number of overlapping DEGs
- Percent of KE covered
- Odds ratio
- p-value and adjusted p-value (FDR)

### Detailed KE Information
For each significant KE:
- Gene list with Ensembl IDs, gene names, log2FC, and padj
- Interactive heatmap of expression changes
- Optional functional enrichment analysis

### Functional Enrichment Results
- Top enriched GO:BP terms and KEGG pathways
- Bar plots showing -log10(p-value)
- Tables with term names, p-values, and gene counts


