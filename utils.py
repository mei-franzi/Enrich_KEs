"""
Utility Functions Module

This module contains shared utility functions used across the application
for data validation, formatting, and file operations.

Author: Enrich_KEs Team
"""

import pandas as pd
import numpy as np
from typing import List, Optional, Dict, Any
import os


def format_scientific_notation(value: float, decimals: int = 2) -> str:
    """
    Format a number in scientific notation.
    
    Parameters
    ----------
    value : float
        Number to format
    decimals : int, optional
        Number of decimal places (default: 2)
    
    Returns
    -------
    str
        Formatted string in scientific notation
    
    Examples
    --------
    >>> format_scientific_notation(0.000001)
    '1.00e-06'
    """
    if pd.isna(value):
        return "NA"
    return f"{value:.{decimals}e}"


def format_percentage(value: float, decimals: int = 1) -> str:
    """
    Format a number as a percentage.
    
    Parameters
    ----------
    value : float
        Number to format (as decimal, e.g., 0.5 for 50%)
    decimals : int, optional
        Number of decimal places (default: 1)
    
    Returns
    -------
    str
        Formatted percentage string
    
    Examples
    --------
    >>> format_percentage(50.0)
    '50.0%'
    """
    if pd.isna(value):
        return "NA"
    return f"{value:.{decimals}f}%"


def format_pvalue_column(df: pd.DataFrame, column: str) -> pd.DataFrame:
    """
    Format a p-value column in scientific notation.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing the column
    column : str
        Name of the column to format
    
    Returns
    -------
    pd.DataFrame
        DataFrame with formatted column
    """
    if column in df.columns:
        df[column] = df[column].apply(lambda x: format_scientific_notation(x) if pd.notna(x) else x)
    return df


def validate_deg_file(df: pd.DataFrame, required_cols: List[str]) -> Dict[str, Any]:
    """
    Validate that a DEG file has required columns.
    
    Parameters
    ----------
    df : pd.DataFrame
        DEG DataFrame to validate
    required_cols : List[str]
        List of required column names
    
    Returns
    -------
    Dict[str, Any]
        Dictionary with 'valid' (bool) and 'missing' (list) keys
    """
    missing_cols = [col for col in required_cols if col not in df.columns]
    
    return {
        'valid': len(missing_cols) == 0,
        'missing': missing_cols
    }


def find_column_by_aliases(df: pd.DataFrame, aliases: List[str]) -> Optional[str]:
    """
    Find a column in DataFrame by checking multiple possible names.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to search
    aliases : List[str]
        List of possible column names to check
    
    Returns
    -------
    Optional[str]
        Name of the found column, or None if not found
    
    Examples
    --------
    >>> df = pd.DataFrame({'gene_name': [1, 2, 3]})
    >>> find_column_by_aliases(df, ['gene', 'gene_name', 'symbol'])
    'gene_name'
    """
    for alias in aliases:
        if alias in df.columns:
            return alias
    return None


def get_gene_name_column(df: pd.DataFrame) -> Optional[str]:
    """
    Find the gene name column in a DataFrame.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to search
    
    Returns
    -------
    Optional[str]
        Name of the gene name column, or None if not found
    """
    gene_aliases = ['gene', 'gene_name', 'symbol', 'gene_symbol', 'Gene', 'Gene Name']
    return find_column_by_aliases(df, gene_aliases)


def check_file_exists(filepath: str) -> bool:
    """
    Check if a file exists at the given path.
    
    Parameters
    ----------
    filepath : str
        Path to the file
    
    Returns
    -------
    bool
        True if file exists, False otherwise
    """
    return os.path.exists(filepath) and os.path.isfile(filepath)


def validate_ensembl_ids(ensembl_ids: pd.Series) -> pd.Series:
    """
    Validate that Ensembl IDs start with 'ENS'.
    
    Parameters
    ----------
    ensembl_ids : pd.Series
        Series of Ensembl IDs
    
    Returns
    -------
    pd.Series
        Boolean series indicating valid IDs
    """
    return ensembl_ids.astype(str).str.startswith("ENS")


def clean_string_column(series: pd.Series) -> pd.Series:
    """
    Clean a string column by removing NaN strings and converting to proper NA.
    
    Parameters
    ----------
    series : pd.Series
        Series to clean
    
    Returns
    -------
    pd.Series
        Cleaned series
    """
    # Replace string 'nan' with actual NA
    series = series.replace(['nan', 'NaN', 'NA', 'None', ''], np.nan)
    return series


def get_file_extension(filename: str) -> str:
    """
    Get the file extension from a filename.
    
    Parameters
    ----------
    filename : str
        Filename or path
    
    Returns
    -------
    str
        File extension (lowercase, without dot)
    
    Examples
    --------
    >>> get_file_extension('data.csv')
    'csv'
    >>> get_file_extension('results.xlsx')
    'xlsx'
    """
    return os.path.splitext(filename)[1].lower().strip('.')


def format_number(value: float, decimals: int = 2) -> str:
    """
    Format a number with specified decimal places.
    
    Parameters
    ----------
    value : float
        Number to format
    decimals : int, optional
        Number of decimal places (default: 2)
    
    Returns
    -------
    str
        Formatted string
    """
    if pd.isna(value):
        return "NA"
    return f"{value:.{decimals}f}"


def create_gene_id_mapping(df: pd.DataFrame, ensembl_col: str, gene_name_col: Optional[str]) -> Dict[str, str]:
    """
    Create a mapping from Ensembl IDs to gene names.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing gene information
    ensembl_col : str
        Name of Ensembl ID column
    gene_name_col : Optional[str]
        Name of gene name column (can be None)
    
    Returns
    -------
    Dict[str, str]
        Dictionary mapping Ensembl IDs to gene names
    """
    if gene_name_col is None or gene_name_col not in df.columns:
        return {}
    
    mapping = {}
    for _, row in df.iterrows():
        ensembl_id = row[ensembl_col]
        gene_name = row[gene_name_col]
        if pd.notna(ensembl_id) and pd.notna(gene_name):
            mapping[ensembl_id] = gene_name
    
    return mapping


def truncate_string(text: str, max_length: int = 50, suffix: str = "...") -> str:
    """
    Truncate a string to a maximum length.
    
    Parameters
    ----------
    text : str
        String to truncate
    max_length : int, optional
        Maximum length (default: 50)
    suffix : str, optional
        Suffix to add when truncating (default: "...")
    
    Returns
    -------
    str
        Truncated string
    """
    if len(text) <= max_length:
        return text
    return text[:max_length - len(suffix)] + suffix


def get_version_info() -> Dict[str, str]:
    """
    Get version information for key packages.
    
    Returns
    -------
    Dict[str, str]
        Dictionary with package names as keys and versions as values
    """
    import importlib.metadata
    import sys
    
    packages = ['streamlit', 'pandas', 'numpy', 'scipy', 'plotly', 'gprofiler-official']
    versions = {'python': sys.version.split()[0]}
    
    for package in packages:
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = 'not installed'
    
    return versions

