"""
Unit tests for ke_enrichment module
"""

import pytest
import pandas as pd
import numpy as np
from ke_enrichment import (
    get_overlapping_genes,
    calculate_contingency_table,
    perform_fishers_test,
    filter_significant_kes,
    format_ke_results_for_display
)


class TestOverlappingGenes:
    """Test gene overlap functions."""
    
    def test_get_overlapping_genes(self):
        """Test finding overlapping genes."""
        degs = {'GENE1', 'GENE2', 'GENE3', 'GENE4'}
        ke_genes = {'GENE2', 'GENE3', 'GENE5', 'GENE6'}
        
        overlap = get_overlapping_genes(degs, ke_genes)
        
        assert overlap == {'GENE2', 'GENE3'}
        assert len(overlap) == 2
    
    def test_get_overlapping_genes_no_overlap(self):
        """Test with no overlapping genes."""
        degs = {'GENE1', 'GENE2'}
        ke_genes = {'GENE3', 'GENE4'}
        
        overlap = get_overlapping_genes(degs, ke_genes)
        
        assert overlap == set()
        assert len(overlap) == 0
    
    def test_get_overlapping_genes_all_overlap(self):
        """Test when all genes overlap."""
        degs = {'GENE1', 'GENE2', 'GENE3'}
        ke_genes = {'GENE1', 'GENE2', 'GENE3'}
        
        overlap = get_overlapping_genes(degs, ke_genes)
        
        assert overlap == degs
        assert len(overlap) == 3


class TestContingencyTable:
    """Test contingency table calculation."""
    
    def test_calculate_contingency_table(self):
        """Test 2x2 contingency table calculation."""
        degs = {'GENE1', 'GENE2', 'GENE3'}
        ke_genes = {'GENE2', 'GENE3', 'GENE4'}
        background_genes = {'GENE1', 'GENE2', 'GENE3', 'GENE4', 'GENE5', 'GENE6'}
        
        a, b, c, d = calculate_contingency_table(degs, ke_genes, background_genes)
        
        assert a == 2  # GENE2, GENE3 (in DEG and in KE)
        assert b == 1  # GENE1 (in DEG, not in KE)
        assert c == 1  # GENE4 (not in DEG, in KE)
        assert d == 2  # GENE5, GENE6 (not in DEG, not in KE)
        
        # Total should equal background
        assert a + b + c + d == len(background_genes)
    
    def test_calculate_contingency_table_no_overlap(self):
        """Test contingency table with no overlap."""
        degs = {'GENE1'}
        ke_genes = {'GENE2'}
        background_genes = {'GENE1', 'GENE2', 'GENE3'}
        
        a, b, c, d = calculate_contingency_table(degs, ke_genes, background_genes)
        
        assert a == 0  # No overlap
        assert b == 1  # GENE1
        assert c == 1  # GENE2
        assert d == 1  # GENE3


class TestFishersTest:
    """Test Fisher's exact test."""
    
    def test_perform_fishers_test(self):
        """Test Fisher's exact test execution."""
        degs = {'GENE1', 'GENE2', 'GENE3'}
        ke_genes = {'GENE2', 'GENE3', 'GENE4'}
        background_genes = {'GENE1', 'GENE2', 'GENE3', 'GENE4', 'GENE5', 'GENE6'}
        
        odds_ratio, p_value, overlap_count = perform_fishers_test(degs, ke_genes, background_genes)
        
        # Check types
        assert isinstance(odds_ratio, float)
        assert isinstance(p_value, float)
        assert isinstance(overlap_count, int)
        
        # Check ranges
        assert odds_ratio >= 0
        assert 0 <= p_value <= 1
        assert overlap_count == 2
    
    def test_perform_fishers_test_no_overlap(self):
        """Test Fisher's test with no overlap."""
        degs = {'GENE1'}
        ke_genes = {'GENE2'}
        background_genes = {'GENE1', 'GENE2', 'GENE3', 'GENE4'}
        
        odds_ratio, p_value, overlap_count = perform_fishers_test(degs, ke_genes, background_genes)
        
        assert overlap_count == 0
        assert p_value > 0.5  # Should not be significant


class TestFilteringAndFormatting:
    """Test result filtering and formatting."""
    
    def test_filter_significant_kes(self):
        """Test filtering for significant KEs."""
        df = pd.DataFrame({
            'KE': ['KE1', 'KE2', 'KE3'],
            'adjusted p-value': [0.001, 0.049, 0.051],
            'p-value': [0.0001, 0.005, 0.006]
        })
        
        significant = filter_significant_kes(df, fdr_threshold=0.05)
        
        assert len(significant) == 2
        assert 'KE1' in significant['KE'].values
        assert 'KE2' in significant['KE'].values
        assert 'KE3' not in significant['KE'].values
    
    def test_filter_significant_kes_custom_threshold(self):
        """Test filtering with custom threshold."""
        df = pd.DataFrame({
            'KE': ['KE1', 'KE2', 'KE3'],
            'adjusted p-value': [0.001, 0.005, 0.02],
            'p-value': [0.0001, 0.0005, 0.002]
        })
        
        significant = filter_significant_kes(df, fdr_threshold=0.01)
        
        assert len(significant) == 2
        assert 'KE3' not in significant['KE'].values
    
    def test_filter_significant_kes_empty(self):
        """Test filtering with empty DataFrame."""
        df = pd.DataFrame()
        
        result = filter_significant_kes(df)
        
        assert result.empty
    
    def test_format_ke_results_for_display(self):
        """Test formatting results for display."""
        df = pd.DataFrame({
            'KE': ['KE1'],
            'p-value': [0.00001],
            'adjusted p-value': [0.0001],
            'Percent of KE covered': [45.5],
            'Odds ratio': [2.345],
            'Overlapping DEGs List': [['GENE1', 'GENE2']]
        })
        
        formatted = format_ke_results_for_display(df)
        
        # Check p-value formatting
        assert formatted['p-value'].iloc[0] == '1.00e-05'
        assert formatted['adjusted p-value'].iloc[0] == '1.00e-04'
        
        # Check percentage formatting
        assert formatted['Percent of KE covered'].iloc[0] == '45.5%'
        
        # Check odds ratio formatting
        assert formatted['Odds ratio'].iloc[0] == '2.35'
        
        # Check that list column is removed
        assert 'Overlapping DEGs List' not in formatted.columns


class TestKEEnrichmentIntegration:
    """Integration tests for KE enrichment workflow."""
    
    def test_enrichment_workflow(self):
        """Test complete enrichment workflow."""
        # Create test data
        ke_map = pd.DataFrame({
            'Gene': ['GENE1', 'GENE2', 'GENE3', 'GENE4', 'GENE5'],
            'KE': ['KE1', 'KE1', 'KE2', 'KE2', 'KE2'],
            'ke.name': ['Event 1', 'Event 1', 'Event 2', 'Event 2', 'Event 2'],
            'AOP': ['AOP1', 'AOP1', 'AOP2', 'AOP2', 'AOP2']
        })
        
        degs = {'GENE1', 'GENE2', 'GENE3'}
        background_genes = {'GENE1', 'GENE2', 'GENE3', 'GENE4', 'GENE5', 'GENE6'}
        
        # Import the main function
        from ke_enrichment import perform_ke_enrichment, apply_fdr_correction
        
        # Perform enrichment
        results = perform_ke_enrichment(degs, ke_map, background_genes)
        
        # Should have results
        assert not results.empty
        assert 'KE' in results.columns
        assert 'p-value' in results.columns
        
        # Apply FDR correction
        results = apply_fdr_correction(results)
        
        # Should have adjusted p-values
        assert 'adjusted p-value' in results.columns


if __name__ == '__main__':
    pytest.main([__file__, '-v'])

