"""
Key Event Enrichment Analysis Module

This module contains functions for performing statistical enrichment analysis
of differentially expressed genes (DEGs) against Key Events (KEs) in Adverse
Outcome Pathways (AOPs) using Fisher's exact test.

Author: Enrich_KEs Team
"""

import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import plotly.graph_objects as go
from typing import List, Set, Dict, Tuple, Optional


def get_overlapping_genes(degs: Set[str], ke_genes: Set[str]) -> Set[str]:
    """
    Find overlapping genes between DEGs and KE genes.
    
    Parameters
    ----------
    degs : Set[str]
        Set of differentially expressed gene IDs
    ke_genes : Set[str]
        Set of genes associated with a Key Event
    
    Returns
    -------
    Set[str]
        Set of overlapping gene IDs
    """
    return degs & ke_genes


def calculate_contingency_table(
    degs: Set[str],
    ke_genes: Set[str],
    background_genes: Set[str]
) -> Tuple[int, int, int, int]:
    """
    Calculate 2x2 contingency table for Fisher's exact test.
    
    Parameters
    ----------
    degs : Set[str]
        Set of differentially expressed gene IDs
    ke_genes : Set[str]
        Set of genes associated with a Key Event
    background_genes : Set[str]
        Set of all background genes
    
    Returns
    -------
    Tuple[int, int, int, int]
        (a, b, c, d) where:
        a = in DEG and in KE
        b = in DEG, not in KE
        c = not in DEG, in KE
        d = not in DEG and not in KE
    """
    overlap = degs & ke_genes
    
    a = len(overlap)  # in DEG and in KE
    b = len(degs - ke_genes)  # in DEG, not in KE
    c = len(ke_genes - degs)  # not in DEG, in KE
    d = len(background_genes - degs - ke_genes)  # not in DEG and not in KE
    
    return a, b, c, d


def perform_fishers_test(
    degs: Set[str],
    ke_genes: Set[str],
    background_genes: Set[str]
) -> Tuple[float, float, int]:
    """
    Perform Fisher's exact test for gene enrichment.
    
    Parameters
    ----------
    degs : Set[str]
        Set of differentially expressed gene IDs
    ke_genes : Set[str]
        Set of genes associated with a Key Event
    background_genes : Set[str]
        Set of all background genes
    
    Returns
    -------
    Tuple[float, float, int]
        (odds_ratio, p_value, overlap_count)
    """
    a, b, c, d = calculate_contingency_table(degs, ke_genes, background_genes)
    odds_ratio, p_value = fisher_exact([[a, b], [c, d]], alternative="greater")
    
    return odds_ratio, p_value, a


def perform_ke_enrichment(
    degs: Set[str],
    ke_map: pd.DataFrame,
    background_genes: Set[str]
) -> pd.DataFrame:
    """
    Perform KE enrichment analysis on DEGs using Fisher's exact test.
    
    Parameters
    ----------
    degs : Set[str]
        Set of differentially expressed gene IDs
    ke_map : pd.DataFrame
        DataFrame with KE-gene mappings (columns: 'Gene', 'KE', 'ke.name', 'AOP')
    background_genes : Set[str]
        Set of all background genes
    
    Returns
    -------
    pd.DataFrame
        DataFrame with enrichment results for each KE
    """
    results = []
    
    # Perform Fisher's exact test for each KE
    for ke_id, group in ke_map.groupby("KE"):
        ke_genes = set(group["Gene"])
        overlap = get_overlapping_genes(degs, ke_genes)
        
        if len(overlap) == 0:  # Skip KEs with no overlap
            continue
        
        # Perform Fisher's exact test
        odds_ratio, p_value, overlap_count = perform_fishers_test(
            degs, ke_genes, background_genes
        )
        
        # Clean KE name and get AOP IDs
        ke_name = group["ke.name"].iloc[0] if "ke.name" in group.columns else ""
        ke_name = "" if pd.isna(ke_name) or str(ke_name).lower() == 'nan' else ke_name
        aop_ids = ", ".join(sorted(set(group["AOP"].dropna())))
        
        # Only include KEs that have a corresponding KE name
        if ke_name:
            results.append({
                "KE": ke_id,
                "KE name": ke_name,
                "AOP": aop_ids,
                "DEGs in KE": overlap_count,
                "KE size": len(ke_genes),
                "Percent of KE covered": (overlap_count / len(ke_genes) * 100),
                "Overlapping DEGs": ", ".join(sorted(overlap)),
                "Overlapping DEGs List": list(sorted(overlap)),
                "Odds ratio": odds_ratio,
                "p-value": p_value
            })
    
    if not results:
        return pd.DataFrame()
    
    # Create DataFrame
    res_df = pd.DataFrame(results)
    
    return res_df


def apply_fdr_correction(
    enrichment_df: pd.DataFrame,
    alpha: float = 0.05,
    method: str = "fdr_bh"
) -> pd.DataFrame:
    """
    Apply multiple testing correction and filter results.
    
    Parameters
    ----------
    enrichment_df : pd.DataFrame
        DataFrame with enrichment results
    alpha : float, optional
        FDR threshold (default: 0.05)
    method : str, optional
        Multiple testing correction method (default: "fdr_bh")
    
    Returns
    -------
    pd.DataFrame
        DataFrame with FDR-corrected p-values, sorted by adjusted p-value
    """
    if enrichment_df.empty:
        return enrichment_df
    
    # Apply multiple testing correction
    enrichment_df["adjusted p-value"] = multipletests(
        enrichment_df["p-value"], 
        method=method
    )[1]
    
    # Sort by adjusted p-value
    enrichment_df = enrichment_df.sort_values("adjusted p-value")
    
    return enrichment_df


def filter_significant_kes(
    enrichment_df: pd.DataFrame,
    fdr_threshold: float = 0.05
) -> pd.DataFrame:
    """
    Filter for significant KEs based on FDR threshold.
    
    Parameters
    ----------
    enrichment_df : pd.DataFrame
        DataFrame with enrichment results and adjusted p-values
    fdr_threshold : float, optional
        FDR threshold for significance (default: 0.05)
    
    Returns
    -------
    pd.DataFrame
        DataFrame with only significant KEs
    """
    if enrichment_df.empty or "adjusted p-value" not in enrichment_df.columns:
        return pd.DataFrame()
    
    return enrichment_df[enrichment_df["adjusted p-value"] < fdr_threshold].copy()


def create_ke_heatmap(
    gene_data: pd.DataFrame,
    ke_name: str,
    colorscale: str = "RdBu_r"
) -> go.Figure:
    """
    Create an interactive Plotly heatmap for gene expression in a KE.
    
    Parameters
    ----------
    gene_data : pd.DataFrame
        DataFrame with columns: 'Gene Name', 'log2FC'
    ke_name : str
        Name of the Key Event
    colorscale : str, optional
        Plotly colorscale name (default: "RdBu_r")
    
    Returns
    -------
    go.Figure
        Plotly figure object with heatmap
    """
    if gene_data.empty or "log2FC" not in gene_data.columns:
        return None
    
    # Sort by log2FC
    heatmap_data = gene_data.sort_values("log2FC", ascending=False)
    
    # Prepare data for heatmap
    log2fc_values = heatmap_data["log2FC"].values.reshape(1, -1)
    
    # Get gene labels
    if "Gene Name" in heatmap_data.columns:
        gene_labels = heatmap_data["Gene Name"].tolist()
    else:
        gene_labels = [f"Gene {i+1}" for i in range(len(heatmap_data))]
    
    # Create heatmap
    fig = go.Figure(data=go.Heatmap(
        z=log2fc_values,
        x=gene_labels,
        y=["log2FC"],
        colorscale=colorscale,
        text=np.round(log2fc_values, 2),
        texttemplate="%{text}",
        textfont={"size": 14},
        colorbar=dict(title="log2FC"),
        hovertemplate="Gene: %{x}<br>log2FC: %{z:.2f}<extra></extra>"
    ))
    
    # Calculate size
    n_genes = len(gene_labels)
    cell_size = 60
    
    # Update layout
    fig.update_layout(
        title=ke_name,
        height=cell_size + 100,
        width=max(n_genes * cell_size, 300),
        margin=dict(l=5, r=5, t=50, b=100),
        xaxis=dict(
            side="bottom",
            tickfont=dict(size=14),
            tickangle=-45
        ),
        yaxis=dict(
            tickfont=dict(size=14)
        ),
        font=dict(size=14)
    )
    
    return fig


def format_ke_results_for_display(enrichment_df: pd.DataFrame) -> pd.DataFrame:
    """
    Format KE enrichment results for display in Streamlit.
    
    Parameters
    ----------
    enrichment_df : pd.DataFrame
        DataFrame with enrichment results
    
    Returns
    -------
    pd.DataFrame
        Formatted DataFrame suitable for display
    """
    if enrichment_df.empty:
        return enrichment_df
    
    df_display = enrichment_df.copy()
    
    # Format p-values in scientific notation
    df_display["p-value"] = df_display["p-value"].apply(lambda x: f"{x:.2e}")
    df_display["adjusted p-value"] = df_display["adjusted p-value"].apply(lambda x: f"{x:.2e}")
    
    # Format other numeric columns
    df_display["Percent of KE covered"] = df_display["Percent of KE covered"].apply(lambda x: f"{x:.1f}%")
    df_display["Odds ratio"] = df_display["Odds ratio"].apply(lambda x: f"{x:.2f}")
    
    # Remove the internal list column from display
    if "Overlapping DEGs List" in df_display.columns:
        df_display = df_display.drop(columns=["Overlapping DEGs List"])
    
    return df_display

