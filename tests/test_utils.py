"""
Unit tests for utils module
"""

import pytest
import pandas as pd
import numpy as np
from utils import (
    format_scientific_notation,
    format_percentage,
    find_column_by_aliases,
    get_gene_name_column,
    validate_ensembl_ids,
    get_file_extension,
    format_number
)


class TestFormattingFunctions:
    """Test formatting utility functions."""
    
    def test_format_scientific_notation(self):
        """Test scientific notation formatting."""
        assert format_scientific_notation(0.00001) == "1.00e-05"
        assert format_scientific_notation(1.5e-10) == "1.50e-10"
        assert format_scientific_notation(0.05) == "5.00e-02"
    
    def test_format_scientific_notation_decimals(self):
        """Test scientific notation with different decimal places."""
        assert format_scientific_notation(0.00001, decimals=3) == "1.000e-05"
        assert format_scientific_notation(0.00001, decimals=1) == "1.0e-05"
    
    def test_format_scientific_notation_na(self):
        """Test formatting with NA values."""
        assert format_scientific_notation(pd.NA) == "NA"
        assert format_scientific_notation(np.nan) == "NA"
    
    def test_format_percentage(self):
        """Test percentage formatting."""
        assert format_percentage(50.0) == "50.0%"
        assert format_percentage(33.333) == "33.3%"
        assert format_percentage(100.0) == "100.0%"
    
    def test_format_percentage_decimals(self):
        """Test percentage formatting with different decimals."""
        assert format_percentage(33.333, decimals=2) == "33.33%"
        assert format_percentage(33.333, decimals=0) == "33%"
    
    def test_format_number(self):
        """Test number formatting."""
        assert format_number(3.14159) == "3.14"
        assert format_number(3.14159, decimals=3) == "3.142"
        assert format_number(100.0, decimals=0) == "100"


class TestColumnFunctions:
    """Test column detection and manipulation functions."""
    
    def test_find_column_by_aliases(self):
        """Test finding column by multiple possible names."""
        df = pd.DataFrame({'gene_name': [1, 2, 3], 'other': [4, 5, 6]})
        
        # Should find existing column
        assert find_column_by_aliases(df, ['gene', 'gene_name', 'symbol']) == 'gene_name'
        
        # Should return None if not found
        assert find_column_by_aliases(df, ['missing1', 'missing2']) is None
    
    def test_get_gene_name_column(self):
        """Test finding gene name column."""
        df1 = pd.DataFrame({'gene': [1], 'value': [2]})
        assert get_gene_name_column(df1) == 'gene'
        
        df2 = pd.DataFrame({'gene_name': [1], 'value': [2]})
        assert get_gene_name_column(df2) == 'gene_name'
        
        df3 = pd.DataFrame({'symbol': [1], 'value': [2]})
        assert get_gene_name_column(df3) == 'symbol'
        
        df4 = pd.DataFrame({'value': [1], 'data': [2]})
        assert get_gene_name_column(df4) is None


class TestValidationFunctions:
    """Test validation functions."""
    
    def test_validate_ensembl_ids(self):
        """Test Ensembl ID validation."""
        ids = pd.Series(['ENSG00000139618', 'ENST00000380152', 'INVALID123', 'ENSG00000141510'])
        result = validate_ensembl_ids(ids)
        
        assert result[0] == True  # ENSG00000139618
        assert result[1] == True  # ENST00000380152
        assert result[2] == False  # INVALID123
        assert result[3] == True  # ENSG00000141510
    
    def test_get_file_extension(self):
        """Test file extension extraction."""
        assert get_file_extension('data.csv') == 'csv'
        assert get_file_extension('results.xlsx') == 'xlsx'
        assert get_file_extension('/path/to/file.tsv') == 'tsv'
        assert get_file_extension('file.txt') == 'txt'


class TestVersionInfo:
    """Test version information function."""
    
    def test_get_version_info(self):
        """Test getting version information."""
        from utils import get_version_info
        
        versions = get_version_info()
        
        # Should have python version
        assert 'python' in versions
        assert isinstance(versions['python'], str)
        
        # Should have key packages
        assert 'pandas' in versions
        assert 'numpy' in versions
        
        # Versions should be strings
        for package, version in versions.items():
            assert isinstance(version, str)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])

