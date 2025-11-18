"""
Functional enrichment analysis module using GProfiler
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import importlib.metadata
import sys
import os
from datetime import datetime
from gprofiler import GProfiler


def perform_functional_enrichment(gene_list, organism='hsapiens', sources=['GO:BP', 'KEGG'], 
                                background_genes=None, max_pval=0.05):
    """
    Perform functional enrichment analysis using GProfiler
    
    Parameters:
    -----------
    gene_list : list
        List of gene symbols or IDs to analyze
    organism : str
        Target organism (default: 'hsapiens' for human)
    sources : list
        Databases to query (default: ['GO:BP', 'KEGG'])
    background_genes : list, optional
        Custom background gene set
    max_pval : float
        P-value significance threshold (default: 0.05)
    
    Returns:
    --------
    pandas.DataFrame
        Enrichment results with gene ratios and fold enrichment
    """
    gp = GProfiler(return_dataframe=True)
    
    try:
        result = gp.profile(
            organism=organism,
            query=gene_list,
            sources=sources,
            user_threshold=max_pval,
            significance_threshold_method='fdr',
            background=background_genes,
            domain_scope='custom' if background_genes else 'annotated',
            no_evidences=False
        )
        
        if result.empty:
            return pd.DataFrame()
        
        # Add gene ratio and background ratio
        result['gene_ratio'] = result['intersection_size'] / result['query_size']
        result['bg_ratio'] = result['term_size'] / result['effective_domain_size']
        result['fold_enrichment'] = result['gene_ratio'] / result['bg_ratio']
        
        return result.sort_values('p_value')
        
    except Exception as e:
        print(f"Enrichment analysis failed: {e}")
        return pd.DataFrame()


def filter_enrichment_results(enrichment_df, source_type='KEGG'):
    """
    Filter out technical artifacts and non-informative root terms
    
    Parameters:
    -----------
    enrichment_df : pandas.DataFrame
        Enrichment results to filter
    source_type : str
        Type of enrichment source ('KEGG' or 'GO')
    
    Returns:
    --------
    pandas.DataFrame
        Filtered enrichment results
    """
    if enrichment_df.empty:
        return enrichment_df
    
    # Define technical artifacts to exclude
    if source_type == 'KEGG':
        exclude_terms = ['KEGG root term', 'root']
    elif source_type == 'GO':
        exclude_terms = ['biological_process', 'molecular_function', 'cellular_component', 'root']
    else:
        exclude_terms = ['root']
    
    # Filter out technical artifacts (case-insensitive)
    filtered_df = enrichment_df[
        ~enrichment_df['name'].str.lower().isin([term.lower() for term in exclude_terms])
    ].copy()
    
    return filtered_df


def create_enrichment_barplot(enrichment_df, title, color='skyblue', max_terms=15):
    """
    Create a horizontal bar plot for enrichment results
    
    Parameters:
    -----------
    enrichment_df : pandas.DataFrame
        Filtered enrichment results
    title : str
        Plot title
    color : str
        Bar color (default: 'skyblue')
    max_terms : int
        Maximum number of terms to display (default: 15)
    
    Returns:
    --------
    matplotlib.figure.Figure
        The created figure object
    """
    if enrichment_df.empty:
        return None
    
    # Select top terms
    top_terms = enrichment_df.head(max_terms).copy()
    
    # Truncate long names
    top_terms['short_name'] = top_terms['name'].apply(
        lambda x: x[:50] + '...' if len(x) > 50 else x
    )
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, max(6, len(top_terms) * 0.4)))
    
    # Create horizontal bar plot
    bars = ax.barh(range(len(top_terms)), -np.log10(top_terms['p_value']),
                   color=color, alpha=0.8)
    
    # Add gene count annotations
    for i, (bar, count) in enumerate(zip(bars, top_terms['intersection_size'])):
        width = bar.get_width()
        ax.text(width + 0.1, bar.get_y() + bar.get_height()/2,
                f'{count}', ha='left', va='center', fontsize=9)
    
    ax.set_yticks(range(len(top_terms)))
    ax.set_yticklabels(top_terms['short_name'], fontsize=10)
    ax.set_xlabel('-log₁₀(FDR)', fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.invert_yaxis()
    
    # Add significance line
    ax.axvline(x=-np.log10(0.05), color='red', linestyle='--', alpha=0.7, label='FDR=0.05')
    ax.legend()
    
    plt.tight_layout()
    
    return fig


def convert_intersections_to_gene_names(intersections, ensembl_to_gene_map):
    """
    Convert intersections (Ensembl IDs or gene symbols) to gene names
    
    Parameters:
    -----------
    intersections : str, list, pandas Series, or numpy array
        Gene IDs from enrichment results (can be Ensembl IDs or gene symbols)
    ensembl_to_gene_map : dict
        Dictionary mapping Ensembl IDs to gene names
    
    Returns:
    --------
    str
        Comma-separated string of gene names
    """
    # Handle pandas Series - convert to scalar first
    if isinstance(intersections, pd.Series):
        if intersections.empty:
            return ""
        # Take the first value if it's a Series
        intersections = intersections.iloc[0] if len(intersections) > 0 else None
    
    # Handle numpy arrays - convert to list
    if isinstance(intersections, np.ndarray):
        if intersections.size == 0:
            return ""
        # Flatten and take first element if 1D, or convert to list
        if intersections.ndim == 0:
            intersections = intersections.item()
        elif intersections.ndim == 1:
            intersections = intersections.tolist()
        else:
            # Multi-dimensional array - take first element
            intersections = intersections.flat[0] if intersections.size > 0 else None
    
    # Check for None first (before any pandas operations)
    if intersections is None:
        return ""
    
    # Check for NaN values (scalar check only)
    try:
        # Only check for scalar NaN values
        if isinstance(intersections, (int, float, str)) or (hasattr(intersections, '__len__') and len(intersections) == 1):
            if pd.isna(intersections):
                return ""
    except (TypeError, ValueError):
        # If pd.isna fails, continue
        pass
    
    # Check if empty (only for sequences, not scalars)
    try:
        if hasattr(intersections, '__len__') and len(intersections) == 0:
            return ""
    except (TypeError, AttributeError):
        pass
    
    # Handle both list and string formats
    if isinstance(intersections, str):
        # Check if it's actually empty string
        if not intersections.strip():
            return ""
        gene_ids = [g.strip() for g in intersections.split(',')]
    elif isinstance(intersections, list):
        # Filter out any None/NaN values from list
        gene_ids = [g for g in intersections if g is not None and not (isinstance(g, float) and pd.isna(g))]
    else:
        # Try to convert to string as last resort
        try:
            intersections_str = str(intersections)
            if intersections_str and intersections_str.lower() not in ['none', 'nan', '']:
                gene_ids = [intersections_str.strip()]
            else:
                return ""
        except:
            return ""
    
    gene_names = []
    for gene_id in gene_ids:
        gene_id = str(gene_id).strip()
        # Check if it's an Ensembl ID and convert
        if gene_id in ensembl_to_gene_map:
            gene_names.append(ensembl_to_gene_map[gene_id])
        elif gene_id.startswith('ENS'):
            # Try to find in mapping
            if gene_id in ensembl_to_gene_map:
                gene_names.append(ensembl_to_gene_map[gene_id])
            else:
                gene_names.append(gene_id)  # Keep as is if not found
        else:
            # Already a gene symbol, keep it
            gene_names.append(gene_id)
    
    return ", ".join(gene_names)


def get_version_info():
    """Get version information for all packages used"""
    packages = [
        'pandas', 'numpy', 'matplotlib', 'seaborn', 'scipy', 
        'networkx', 'py4cytoscape', 'gprofiler-official', 
        'requests', 'openpyxl', 'streamlit', 'plotly'
    ]
    
    version_info = {}
    for package in packages:
        try:
            version = importlib.metadata.version(package)
            version_info[package] = version
        except importlib.metadata.PackageNotFoundError:
            version_info[package] = "Not installed"
    
    version_info['Python'] = sys.version.split()[0]
    return version_info
