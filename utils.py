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
import io
import matplotlib.pyplot as plt
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter, A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak, Image
from reportlab.lib.enums import TA_CENTER, TA_LEFT, TA_RIGHT
from reportlab.pdfgen import canvas
from datetime import datetime


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


def create_ke_heatmap_figure(gene_names: List[str], log2fc_values: List[float], 
                              ke_name: str, ke_id: str, aop: Optional[str] = None) -> plt.Figure:
    """
    Create a matplotlib figure for a Key Event heatmap.
    
    Parameters
    ----------
    gene_names : List[str]
        List of gene names
    log2fc_values : List[float]
        List of log2 fold change values
    ke_name : str
        Name of the Key Event
    ke_id : str
        ID of the Key Event
    aop : Optional[str]
        AOP identifier (optional, will be included in title if provided)
    
    Returns
    -------
    plt.Figure
        Matplotlib figure object
    """
    n_genes = len(gene_names)
    if n_genes > 35:
        fig_height = max(n_genes * 0.15, 10)
    else:
        fig_height = max(n_genes * (4/18), 3)
    
    fig, ax = plt.subplots(figsize=(9, fig_height))
    colors_list = ['#FF6B6B' if x > 0 else '#4ECDC4' for x in log2fc_values]
    
    ax.barh(gene_names, log2fc_values, color=colors_list, alpha=1)
    ax.set_xlabel('log2 Fold Change', fontsize=12)
    
    # Match app title format: include AOP if available
    if aop:
        ax.set_title(f'{ke_name} ({ke_id}, {aop})', fontsize=13)
    else:
        ax.set_title(f'{ke_name} ({ke_id})', fontsize=13)
    
    ax.axvline(x=0, color='black', linestyle='-', linewidth=0.5)
    ax.grid(axis='x', alpha=0.3)
    
    plt.tight_layout()
    return fig


def generate_ke_pdf(ke_data_list: List[Dict[str, Any]], analysis_name: str = "KE Enrichment Analysis") -> bytes:
    """
    Generate a PDF document containing all Key Event heatmaps, info boxes, and gene tables.
    
    Parameters
    ----------
    ke_data_list : List[Dict[str, Any]]
        List of dictionaries, each containing:
        - 'ke_id': str
        - 'ke_name': str
        - 'ke_row': pd.Series with KE statistics
        - 'gene_details': List[Dict] with gene information
        - 'gene_names': List[str]
        - 'log2fc_values': List[float]
    analysis_name : str
        Name of the analysis (default: "KE Enrichment Analysis")
    
    Returns
    -------
    bytes
        PDF file as bytes
    """
    buffer = io.BytesIO()
    doc = SimpleDocTemplate(buffer, pagesize=letter,
                           rightMargin=72, leftMargin=72,
                           topMargin=72, bottomMargin=18)
    
    # Container for the 'Flowable' objects
    elements = []
    
    # Define styles
    styles = getSampleStyleSheet()
    title_style = ParagraphStyle(
        'CustomTitle',
        parent=styles['Heading1'],
        fontSize=18,
        textColor=colors.HexColor('#1f77b4'),
        spaceAfter=30,
        alignment=TA_CENTER
    )
    
    heading_style = ParagraphStyle(
        'CustomHeading',
        parent=styles['Heading2'],
        fontSize=14,
        textColor=colors.HexColor('#2c3e50'),
        spaceAfter=12,
        spaceBefore=12
    )
    
    ke_title_style = ParagraphStyle(
        'KETitle',
        parent=styles['Heading2'],
        fontSize=14,
        textColor=colors.HexColor('#4ECDC4'),
        spaceAfter=10,
        spaceBefore=20
    )
    
    normal_style = styles['Normal']
    normal_style.fontSize = 10
    
    # Title
    elements.append(Paragraph(analysis_name, title_style))
    elements.append(Spacer(1, 0.2*inch))
    
    # Date
    date_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    elements.append(Paragraph(f"Generated on: {date_str}", normal_style))
    elements.append(Spacer(1, 0.3*inch))
    
    # Summary
    elements.append(Paragraph(f"Total Key Events: {len(ke_data_list)}", heading_style))
    elements.append(Spacer(1, 0.2*inch))
    
    # Process each Key Event
    for idx, ke_data in enumerate(ke_data_list, 1):
        ke_id = ke_data['ke_id']
        ke_name = ke_data['ke_name']
        ke_row = ke_data['ke_row']
        gene_names = ke_data['gene_names']
        log2fc_values = ke_data['log2fc_values']
        gene_details = ke_data['gene_details']
        
        # Key Event Title
        elements.append(Paragraph(f"{idx}. {ke_name} ({ke_id})", ke_title_style))
        elements.append(Spacer(1, 0.1*inch))
        
        # Key Event Information Box
        info_data = [
            ['KE ID:', ke_id],
            ['KE Name:', ke_name],
            ['AOP:', str(ke_row.get('AOP', 'N/A'))],
            ['DEGs in KE:', str(ke_row.get('DEGs in KE', 'N/A'))],
            ['KE Size:', str(ke_row.get('KE size', 'N/A'))],
            ['Percent Covered:', f"{ke_row.get('Percent of KE covered', 0):.1f}%"],
            ['Adjusted p-value:', f"{ke_row.get('adjusted p-value', 0):.2e}"],
            ['Odds Ratio:', f"{ke_row.get('Odds ratio', 0):.2f}"]
        ]
        
        info_table = Table(info_data, colWidths=[2*inch, 4*inch])
        info_table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (0, -1), colors.HexColor('#f0f2f6')),
            ('TEXTCOLOR', (0, 0), (-1, -1), colors.black),
            ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
            ('FONTNAME', (0, 0), (0, -1), 'Helvetica-Bold'),
            ('FONTNAME', (1, 0), (1, -1), 'Helvetica'),
            ('FONTSIZE', (0, 0), (-1, -1), 10),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 8),
            ('TOPPADDING', (0, 0), (-1, -1), 8),
            ('LEFTPADDING', (0, 0), (-1, -1), 12),
            ('RIGHTPADDING', (0, 0), (-1, -1), 12),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.grey),
        ]))
        elements.append(info_table)
        elements.append(Spacer(1, 0.2*inch))
        
        # Heatmap Figure
        try:
            # Get AOP from ke_row if available
            aop = ke_row.get('AOP', None)
            fig = create_ke_heatmap_figure(gene_names, log2fc_values, ke_name, ke_id, aop=aop)
            img_buffer = io.BytesIO()
            
            # Save at high DPI to match app quality (300 DPI for print quality)
            # Use the same figure size as displayed in app (9 inches width)
            fig.savefig(img_buffer, format='png', dpi=300, bbox_inches='tight', 
                       facecolor='white', edgecolor='none')
            img_buffer.seek(0)
            
            # Calculate dimensions matching the app display
            # App uses figsize=(9, fig_height), so we preserve aspect ratio
            n_genes = len(gene_names)
            if n_genes > 35:
                fig_height = max(n_genes * 0.15, 10)
            else:
                fig_height = max(n_genes * (4/18), 3)
            
            # Use full width available (letter size is 8.5 inches, minus margins = ~7 inches)
            # Scale to fit page width while preserving aspect ratio
            pdf_width = 7*inch  # Full width minus margins
            aspect_ratio = fig_height / 9.0  # height/width ratio
            pdf_height = pdf_width * aspect_ratio
            
            # But limit max height to prevent overly tall images (max 9 inches)
            max_height = 9*inch
            if pdf_height > max_height:
                pdf_height = max_height
                pdf_width = pdf_height / aspect_ratio
            
            # Add image to PDF with high quality
            img = Image(img_buffer, width=pdf_width, height=pdf_height)
            elements.append(img)
            elements.append(Spacer(1, 0.2*inch))
            plt.close(fig)
        except Exception as e:
            elements.append(Paragraph(f"Error generating heatmap: {str(e)}", normal_style))
        
        # Gene Table
        elements.append(Paragraph("DEGs in Key Event:", heading_style))
        elements.append(Spacer(1, 0.1*inch))
        
        # Prepare table data
        table_data = [['Gene Name', 'Ensembl ID', 'log2FC', 'padj']]
        for gene in gene_details:
            gene_name = gene.get('Gene Name', gene.get('Ensembl ID', 'N/A'))
            ensembl_id = gene.get('Ensembl ID', 'N/A')
            log2fc = f"{gene.get('log2FoldChange', 0):.3f}"
            padj = format_scientific_notation(gene.get('padj', 0))
            table_data.append([gene_name, ensembl_id, log2fc, padj])
        
        # Create table
        gene_table = Table(table_data, colWidths=[1.5*inch, 2*inch, 1*inch, 1.5*inch])
        gene_table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#4ECDC4')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, 0), 11),
            ('FONTSIZE', (0, 1), (-1, -1), 9),
            ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
            ('TOPPADDING', (0, 0), (-1, 0), 12),
            ('BACKGROUND', (0, 1), (-1, -1), colors.white),
            ('TEXTCOLOR', (0, 1), (-1, -1), colors.black),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.grey),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#f9f9f9')]),
        ]))
        elements.append(gene_table)
        
        # Add page break if not last item
        if idx < len(ke_data_list):
            elements.append(PageBreak())
    
    # Build PDF
    doc.build(elements)
    buffer.seek(0)
    return buffer.getvalue()

