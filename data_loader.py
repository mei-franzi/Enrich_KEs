"""
Data Loading Module

This module handles loading and preprocessing of data files including
DEG files, KE mappings, and KE descriptions.

Author: Enrich_KEs Team
"""

import pandas as pd
import os
from typing import Optional, Tuple, Dict, List
import streamlit as st
from utils import get_file_extension, validate_ensembl_ids


def get_excel_sheet_names(file_source) -> Optional[List[str]]:
    """
    Get list of sheet names from an Excel file.
    
    Parameters
    ----------
    file_source : UploadedFile or str
        Either a Streamlit uploaded file object or a file path
    
    Returns
    -------
    Optional[List[str]]
        List of sheet names, or None if reading fails
    """
    try:
        excel_file = pd.ExcelFile(file_source)
        return excel_file.sheet_names
    except Exception as e:
        st.error(f"Error reading Excel sheets: {str(e)}")
        return None


def load_excel_sheet(file_source, sheet_name: str = 0) -> Optional[pd.DataFrame]:
    """
    Load a specific sheet from an Excel file.
    
    Parameters
    ----------
    file_source : UploadedFile or str
        Either a Streamlit uploaded file object or a file path
    sheet_name : str or int, optional
        Name or index of the sheet to load (default: 0 for first sheet)
    
    Returns
    -------
    Optional[pd.DataFrame]
        Loaded DataFrame, or None if loading fails
    """
    try:
        df = pd.read_excel(file_source, sheet_name=sheet_name)
        return df
    except Exception as e:
        st.error(f"Error loading Excel sheet '{sheet_name}': {str(e)}")
        return None


def load_deg_file(uploaded_file, sheet_name: Optional[str] = None) -> Optional[pd.DataFrame]:
    """
    Load a DEG file from various formats (CSV, TSV, Excel).
    
    Parameters
    ----------
    uploaded_file : UploadedFile
        Streamlit uploaded file object
    sheet_name : Optional[str], optional
        For Excel files, name of the sheet to load (default: None, loads first sheet)
    
    Returns
    -------
    Optional[pd.DataFrame]
        Loaded DataFrame, or None if loading fails
    """
    if uploaded_file is None:
        return None
    
    try:
        file_extension = get_file_extension(uploaded_file.name)
        
        if file_extension == 'csv':
            # Try different separators
            try:
                df = pd.read_csv(uploaded_file, sep=',')
            except:
                uploaded_file.seek(0)
                try:
                    df = pd.read_csv(uploaded_file, sep=';')
                except:
                    uploaded_file.seek(0)
                    df = pd.read_csv(uploaded_file, sep='\t')
        
        elif file_extension == 'tsv':
            df = pd.read_csv(uploaded_file, sep='\t')
        
        elif file_extension in ['xlsx', 'xls']:
            # Use sheet_name if provided, otherwise load first sheet (0)
            sheet = sheet_name if sheet_name is not None else 0
            df = load_excel_sheet(uploaded_file, sheet_name=sheet)
            if df is None:
                return None
        
        else:
            st.error(f"Unsupported file format: {file_extension}")
            return None
        
        return df
    
    except Exception as e:
        st.error(f"Error loading file: {str(e)}")
        return None


def load_deg_from_path(filepath: str, sheet_name: Optional[str] = None) -> Optional[pd.DataFrame]:
    """
    Load a DEG file from a file path (for example data).
    
    Parameters
    ----------
    filepath : str
        Path to the DEG file
    sheet_name : Optional[str], optional
        For Excel files, name of the sheet to load (default: None, loads first sheet)
    
    Returns
    -------
    Optional[pd.DataFrame]
        Loaded DataFrame, or None if loading fails
    """
    if not os.path.exists(filepath):
        return None
    
    try:
        file_extension = get_file_extension(filepath)
        
        if file_extension == 'csv':
            # Try different separators (semicolon first for European CSVs)
            try:
                df = pd.read_csv(filepath, sep=';')
                # Check if it parsed correctly (more than 1 column)
                if len(df.columns) == 1:
                    raise ValueError("Only 1 column detected, trying different separator")
            except:
                try:
                    df = pd.read_csv(filepath, sep=',')
                except:
                    df = pd.read_csv(filepath, sep='\t')
        
        elif file_extension == 'tsv':
            df = pd.read_csv(filepath, sep='\t')
        
        elif file_extension in ['xlsx', 'xls']:
            # Use sheet_name if provided, otherwise load first sheet (0)
            sheet = sheet_name if sheet_name is not None else 0
            df = load_excel_sheet(filepath, sheet_name=sheet)
            if df is None:
                return None
        
        else:
            return None
        
        return df
    
    except Exception as e:
        return None


@st.cache_data
def load_ke_mapping(filepath: str) -> Optional[pd.DataFrame]:
    """
    Load KE-to-gene mapping file.
    
    Parameters
    ----------
    filepath : str
        Path to the KE mapping file
    
    Returns
    -------
    Optional[pd.DataFrame]
        Loaded DataFrame with columns: Gene, KE, ke.name, AOP
    """
    if not os.path.exists(filepath):
        return None
    
    try:
        df = pd.read_csv(filepath, sep="\t")
        
        # Validate required columns
        required_cols = ["Gene", "KE"]
        missing_cols = [col for col in required_cols if col not in df.columns]
        
        if missing_cols:
            st.error(f"KE mapping file missing required columns: {missing_cols}")
            return None
        
        # Drop rows with missing Gene or KE
        df = df.dropna(subset=["Gene", "KE"])
        
        return df
    
    except Exception as e:
        st.error(f"Error loading KE mapping file: {str(e)}")
        return None


@st.cache_data
def load_ke_descriptions(filepath: str) -> Optional[pd.DataFrame]:
    """
    Load KE descriptions file.
    
    Parameters
    ----------
    filepath : str
        Path to the KE descriptions file
    
    Returns
    -------
    Optional[pd.DataFrame]
        Loaded DataFrame with KE descriptions
    """
    if not os.path.exists(filepath):
        return None
    
    try:
        df = pd.read_csv(filepath)
        
        # Validate that 'KE' column exists
        if 'KE' not in df.columns:
            st.error("KE descriptions file missing 'KE' column")
            return None
        
        return df
    
    except Exception as e:
        st.error(f"Error loading KE descriptions file: {str(e)}")
        return None


@st.cache_data(show_spinner="Loading KE data...")
def prepare_ke_data(ke_map_path: str, ke_desc_path: str) -> Tuple[Optional[pd.DataFrame], Optional[set]]:
    """
    Load and prepare KE mapping and description data.
    
    Parameters
    ----------
    ke_map_path : str
        Path to KE mapping file
    ke_desc_path : str
        Path to KE descriptions file
    
    Returns
    -------
    Tuple[Optional[pd.DataFrame], Optional[set]]
        (ke_map_with_descriptions, background_genes) or (None, None) if loading fails
    """
    # Load KE mapping
    ke_map = load_ke_mapping(ke_map_path)
    if ke_map is None:
        return None, None
    
    # Load KE descriptions
    ke_desc = load_ke_descriptions(ke_desc_path)
    if ke_desc is None:
        # KE descriptions are optional, can proceed without them
        ke_desc = pd.DataFrame({'KE': ke_map['KE'].unique()})
    
    # Merge KE map with descriptions
    ke_map_merged = ke_map.merge(ke_desc, on="KE", how="left")
    
    # Get background genes
    background_genes = set(ke_map["Gene"].unique())
    
    return ke_map_merged, background_genes


def apply_column_mapping(df: pd.DataFrame, mapping: Dict[str, str]) -> pd.DataFrame:
    """
    Apply column name mapping to DataFrame.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to rename columns
    mapping : Dict[str, str]
        Dictionary mapping old column names to new column names
    
    Returns
    -------
    pd.DataFrame
        DataFrame with renamed columns
    """
    # Only rename columns that exist in the DataFrame
    valid_mapping = {old: new for old, new in mapping.items() if old in df.columns}
    
    if valid_mapping:
        df = df.rename(columns=valid_mapping)
    
    return df


def filter_degs(
    df: pd.DataFrame,
    padj_cutoff: float,
    log2fc_cutoff: float,
    padj_col: str = "padj",
    log2fc_col: str = "log2FoldChange",
    ensembl_col: str = "human_ensembl_id"
) -> pd.DataFrame:
    """
    Filter DEGs based on padj and log2FC cutoffs.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with DEG data
    padj_cutoff : float
        Adjusted p-value threshold
    log2fc_cutoff : float
        Absolute log2 fold change threshold
    padj_col : str, optional
        Name of adjusted p-value column (default: "padj")
    log2fc_col : str, optional
        Name of log2 fold change column (default: "log2FoldChange")
    ensembl_col : str, optional
        Name of Ensembl ID column (default: "human_ensembl_id")
    
    Returns
    -------
    pd.DataFrame
        Filtered DataFrame
    """
    # Check if required columns exist
    required_cols = [padj_col, log2fc_col, ensembl_col]
    missing_cols = [col for col in required_cols if col not in df.columns]
    
    if missing_cols:
        st.warning(f"Cannot filter: missing columns {missing_cols}")
        return df
    
    # Apply filters
    filtered_df = df[
        (df[padj_col] < padj_cutoff) & 
        (df[log2fc_col].abs() > log2fc_cutoff) & 
        (df[ensembl_col].notna()) & 
        (validate_ensembl_ids(df[ensembl_col]))
    ].copy()
    
    return filtered_df


def validate_deg_columns(df: pd.DataFrame) -> Dict[str, bool]:
    """
    Check which standard DEG columns are present in the DataFrame.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to validate
    
    Returns
    -------
    Dict[str, bool]
        Dictionary indicating presence of each standard column
    """
    standard_cols = {
        'padj': False,
        'log2FoldChange': False,
        'human_ensembl_id': False,
        'pvalue': False,
        'gene_name': False
    }
    
    # Check for padj column
    padj_aliases = ['padj', 'adj_pvalue', 'FDR', 'qvalue', 'adjusted_pvalue']
    if any(col in df.columns for col in padj_aliases):
        standard_cols['padj'] = True
    
    # Check for log2FoldChange column
    log2fc_aliases = ['log2FoldChange', 'log2FC', 'logFC', 'log2fc']
    if any(col in df.columns for col in log2fc_aliases):
        standard_cols['log2FoldChange'] = True
    
    # Check for Ensembl ID column
    ensembl_aliases = ['human_ensembl_id', 'ensembl_id', 'ensembl', 'gene_id']
    if any(col in df.columns for col in ensembl_aliases):
        standard_cols['human_ensembl_id'] = True
    
    # Check for pvalue column
    pvalue_aliases = ['pvalue', 'p_value', 'pval', 'PValue']
    if any(col in df.columns for col in pvalue_aliases):
        standard_cols['pvalue'] = True
    
    # Check for gene name column
    gene_aliases = ['gene', 'gene_name', 'symbol', 'gene_symbol', 'Gene']
    if any(col in df.columns for col in gene_aliases):
        standard_cols['gene_name'] = True
    
    return standard_cols


def get_column_options(df: pd.DataFrame) -> list:
    """
    Get list of column names from DataFrame for selection widgets.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame
    
    Returns
    -------
    list
        List of column names
    """
    return list(df.columns)


def preview_data(df: pd.DataFrame, n_rows: int = 5) -> pd.DataFrame:
    """
    Get a preview of the DataFrame.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to preview
    n_rows : int, optional
        Number of rows to show (default: 5)
    
    Returns
    -------
    pd.DataFrame
        Preview DataFrame
    """
    return df.head(n_rows)

