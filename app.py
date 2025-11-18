# Libraries 
import pandas as pd
import os
import streamlit as st
import matplotlib.pyplot as plt
from datetime import datetime

# Project modules
from enrichment import perform_functional_enrichment, filter_enrichment_results, create_enrichment_barplot, convert_intersections_to_gene_names, wrap_gene_names
from ke_enrichment import (
    perform_ke_enrichment,
    apply_fdr_correction,
    filter_significant_kes,
    format_ke_results_for_display
)
from data_loader import load_deg_file, load_deg_from_path, prepare_ke_data, apply_column_mapping, filter_degs, get_excel_sheet_names
from utils import format_scientific_notation, get_gene_name_column, generate_ke_pdf, generate_ke_html_report

st.set_page_config(layout="wide", page_title="KE & Functional Enrichment")

# =============================================================================
# SIDEBAR
# =============================================================================
with st.sidebar:
    st.title("🧬 DEG Analyzer")
    st.markdown("Analysis Management")
    st.markdown("---")
    
    # Initialize number of analysis tabs
    if 'num_analyses' not in st.session_state:
        st.session_state.num_analyses = 1
    
    # Tab management buttons
    if st.button("➕ Add Analysis", help="Add new analysis tab", use_container_width=True):
        if st.session_state.num_analyses < 10:
            st.session_state.num_analyses += 1
            st.rerun()
    
    if st.button("➖ Remove Analysis", help="Remove last analysis tab", use_container_width=True):
        if st.session_state.num_analyses > 1:
            st.session_state.num_analyses -= 1
            st.rerun()
    
    if st.button("🗑️ Clear All", help="Clear all analyses", use_container_width=True):
        st.session_state.num_analyses = 1
        # Clear all analysis-specific session state
        keys_to_delete = [k for k in st.session_state.keys() if k.startswith('tab')]
        for k in keys_to_delete:
            del st.session_state[k]
        st.rerun()
    
    st.markdown("---")
    st.caption(f"Current analyses: {st.session_state.num_analyses}")

# =============================================================================
# MAIN CONTENT
# =============================================================================

st.title("🧬 Enrich DEGs for toxicological Key Events")
st.markdown("Differential Expression Gene Analysis with KE Enrichment")

# File paths 
data_dir = "data"
ke_map_path = os.path.join(data_dir, "Genes_to_KEs.txt")
ke_desc_path = os.path.join(data_dir, "ke_descriptions.csv")

# Check if reference files exist
if not os.path.exists(ke_map_path):
    st.error(f"Missing file: {ke_map_path}")
    st.info("Please create a 'data' folder in your project directory and add the required reference files:")
    st.write("- Genes_to_KEs.txt")
    st.write("- ke_descriptions.csv")
    st.stop()
if not os.path.exists(ke_desc_path):
    st.error(f"Missing file: {ke_desc_path}")
    st.stop()

# Load KE mapping and descriptions
ke_map, background_genes = prepare_ke_data(ke_map_path, ke_desc_path)
if ke_map is None:
    st.error("Failed to load KE data. Please check the data files.")
    st.stop()

# Create tabs with custom names
tab_names = []
for i in range(1, st.session_state.num_analyses + 1):
    key_prefix = f"tab{i}"
    # Get custom name from session state, or use default
    if f"{key_prefix}_name" in st.session_state and st.session_state[f"{key_prefix}_name"]:
        tab_names.append(st.session_state[f"{key_prefix}_name"])
    else:
        tab_names.append(f"Analysis {i}")

tabs = st.tabs(tab_names)

# =============================================================================
# ANALYSIS WORKFLOW (repeated in each tab with unique session state keys)
# =============================================================================

for tab_idx in range(st.session_state.num_analyses):
    with tabs[tab_idx]:
        analysis_num = tab_idx + 1
        
        # Unique key prefix for this tab
        key_prefix = f"tab{analysis_num}"
        
        # Analysis name input
        analysis_name = st.text_input(
            "Name this analysis",
            value=st.session_state.get(f"{key_prefix}_name", ""),
            placeholder=f"Analysis {analysis_num}",
            key=f"{key_prefix}_name_input",
            help="Name this analysis (will appear in tab title)"
        )
        
        # Update session state if name changed
        if analysis_name != st.session_state.get(f"{key_prefix}_name", ""):
            st.session_state[f"{key_prefix}_name"] = analysis_name
            st.rerun()

        #st.markdown("---")
        
        # =============================================================================
        # 1. UPLOAD DEGs
        # =============================================================================
        #st.markdown("#### Upload DEGs")

        uploaded_file = st.file_uploader(
            "**Upload differential expression results**",
            type=["csv", "tsv", "xlsx", "xls"],
            help="File should contain: adjusted p-values, log2 fold change, and human Ensembl IDs. You can map your column names after upload.",
            key=f"{key_prefix}_uploader"
        )

        # Example data buttons
        st.write("Or use example data:")

        example_data_path = os.path.join(data_dir, "GSE255602_DEGs.csv")
        browder_data_path = os.path.join(data_dir, "Browder_DEGs.xlsx")

        use_example = False
        use_browder = False

        if os.path.exists(example_data_path):
            if st.button("GSE255602_DEGs.csv", key=f"{key_prefix}_gse"):
                use_example = True
                st.session_state[f"{key_prefix}_example_type"] = 'gse'

        if os.path.exists(browder_data_path):
            if st.button("Browder_DEGs.xlsx", key=f"{key_prefix}_browder"):
                use_browder = True
                st.session_state[f"{key_prefix}_example_type"] = 'browder'

        # Determine which file source to use
        file_source = None
        file_name = None
        is_excel = False

        if uploaded_file is not None:
            file_source = uploaded_file
            file_name = uploaded_file.name
            is_excel = file_name.endswith(('.xlsx', '.xls'))
        elif use_example or (f"{key_prefix}_example_type" in st.session_state and st.session_state[f"{key_prefix}_example_type"] == 'gse'):
            file_source = example_data_path
            file_name = "GSE255602_DEGs.csv"
            is_excel = False
        elif use_browder or (f"{key_prefix}_example_type" in st.session_state and st.session_state[f"{key_prefix}_example_type"] == 'browder'):
            file_source = browder_data_path
            file_name = "Browder_DEGs.xlsx"
            is_excel = True

        # Sheet selection for Excel files
        selected_sheet = None
        if file_source and is_excel:
            # Get sheet names
            sheet_names = get_excel_sheet_names(file_source)
            
            if sheet_names:
                selected_sheet = st.selectbox(
                    "Choose which sheet to analyze",
                    options=sheet_names,
                    help="Select the sheet containing your DEG data",
                    key=f"{key_prefix}_sheet"
                )
            else:
                st.error("Could not read Excel sheets")
                file_source = None

        # Load data
        if file_source is not None:
            try:
                    # Load the file with sheet selection if applicable
                    if isinstance(file_source, str):  # File path (example data)
                        deg_file_data = load_deg_from_path(file_source, sheet_name=selected_sheet)
                    else:  # Uploaded file
                        deg_file_data = load_deg_file(file_source, sheet_name=selected_sheet)
        
                    if deg_file_data is None:
                        st.error("Failed to load file. Please check the format.")
                        st.stop()

                    # Store dataset label for reporting
                    dataset_label = file_name if file_name else "Uploaded dataset"
                    st.session_state[f"{key_prefix}_dataset_label"] = dataset_label
                    sheet_label = selected_sheet if selected_sheet else "Not specified"
                    st.session_state[f"{key_prefix}_dataset_sheet"] = sheet_label
        
                    # =============================================================================
                    # 2. COLUMN MAPPING
                    # =============================================================================
                    #st.markdown("---")
                    with st.expander("Specify which columns in your file correspond to the required data", expanded=True):
                        col_map1, col_map2 = st.columns(2)
        
                        with col_map1:
                            padj_col = st.selectbox(
                                "Adjusted p-value column",
                                options=deg_file_data.columns.tolist(),
                                index=deg_file_data.columns.tolist().index('padj') if 'padj' in deg_file_data.columns else 0,
                                help="Column containing adjusted p-values",
                                key=f"{key_prefix}_padj_col"
                            )
        
                        with col_map2:
                            log2fc_col = st.selectbox(
                                "Log2 Fold Change column",
                                options=deg_file_data.columns.tolist(),
                                index=deg_file_data.columns.tolist().index('log2FoldChange') if 'log2FoldChange' in deg_file_data.columns else 0,
                                help="Column containing log2 fold change values",
                                key=f"{key_prefix}_log2fc_col"
                            )
        
                        col_map3, col_map4 = st.columns(2)
        
                        with col_map3:
                            ensembl_col = st.selectbox(
                                "Human Ensembl ID column",
                                options=deg_file_data.columns.tolist(),
                                index=deg_file_data.columns.tolist().index('human_ensembl_id') if 'human_ensembl_id' in deg_file_data.columns else 0,
                                help="Column containing human Ensembl gene IDs",
                                key=f"{key_prefix}_ensembl_col"
                            )
        
                        with col_map4:
                            # Try to find gene name column with common aliases
                            gene_col_idx = 0
                            for idx, col in enumerate(deg_file_data.columns.tolist()):
                                if col.lower() in ['gene', 'gene_name', 'symbol', 'gene_symbol', 'genesymbol']:
                                    gene_col_idx = idx
                                    break
            
                            gene_name_col = st.selectbox(
                                "Gene name/symbol column",
                                options=deg_file_data.columns.tolist(),
                                index=gene_col_idx,
                                help="Column containing gene names or symbols (required for functional enrichment)",
                                key=f"{key_prefix}_gene_col"
                            )
        
                        # Rename columns to standard names
                        rename_dict = {
                            padj_col: 'padj',
                            log2fc_col: 'log2FoldChange',
                            ensembl_col: 'human_ensembl_id',
                            gene_name_col: 'gene'
                        }
                        deg_file_data = apply_column_mapping(deg_file_data, rename_dict)
        
                    # =============================================================================
                    # 3. VIEW AND FILTER DEGs
                    # =============================================================================
                    with st.expander("View Differentially expressed Genes Table", expanded=False):
                        col_filter1, col_filter2, col_filter3 = st.columns(3)
        
                        with col_filter1:
                            padj_cutoff = st.number_input(
                                "Adjusted p-value cutoff",
                                min_value=0.0,
                                max_value=1.0,
                                value=0.05,
                                step=0.01,
                                format="%.3f",
                                key=f"{key_prefix}_padj_cutoff"
                            )
        
                        with col_filter2:
                            log2fc_cutoff = st.number_input(
                                "Log2 Fold Change cutoff",
                                min_value=0.0,
                                max_value=10.0,
                                value=0.1,
                                step=0.1,
                                format="%.2f",
                                key=f"{key_prefix}_log2fc_cutoff"
                            )
        
                        with col_filter3:
                            direction_options = ["All DEGs", "Up-regulated only", "Down-regulated only"]
                            direction_selection = st.selectbox(
                                "Direction",
                                options=direction_options,
                                help="Filter genes by regulation direction",
                                key=f"{key_prefix}_direction"
                            )
                        
                        # Apply filters
                        filtered_df = filter_degs(deg_file_data, padj_cutoff, log2fc_cutoff)
                        
                        # Apply directional filter
                        if direction_selection == "Up-regulated only":
                            filtered_df = filtered_df[filtered_df['log2FoldChange'] > 0]
                        elif direction_selection == "Down-regulated only":
                            filtered_df = filtered_df[filtered_df['log2FoldChange'] < 0]
                        
                        st.write(f"**Filtered DEGs: {len(filtered_df)} genes** (padj < {padj_cutoff:.3g}, |log2FC| > {log2fc_cutoff:.3g}, {direction_selection})")
                        
                        # Format preview dataframe
                        preview_df = filtered_df.copy()
                        if 'padj' in preview_df.columns:
                            preview_df['padj'] = preview_df['padj'].apply(format_scientific_notation)
                        if 'pvalue' in preview_df.columns:
                            preview_df['pvalue'] = preview_df['pvalue'].apply(format_scientific_notation)
                        
                        st.dataframe(preview_df, use_container_width=True, height=400)
        
                    # =============================================================================
                    # 4. FUNCTIONAL ENRICHMENT
                    # =============================================================================
                    #st.markdown("---")
                    
                    # Simple button (not primary/red)
                    run_functional = st.button("Run Functional Enrichment Analysis", key=f"{key_prefix}_enrich_btn")
        
                    if run_functional:
                        with st.spinner("Running functional enrichment analysis..."):
                            # Check if gene names are available
                            gene_col = get_gene_name_column(filtered_df)
                            if gene_col:
                                gene_list = filtered_df[gene_col].dropna().unique().tolist()
                            else:
                                st.warning("⚠️ Gene names not found in data. Enrichment analysis requires gene symbols.")
                                gene_list = []
                
                            if gene_list:
                                # Create mapping from Ensembl IDs to gene names for intersection conversion
                                ensembl_to_gene = {}
                                if 'human_ensembl_id' in filtered_df.columns and 'gene' in filtered_df.columns:
                                    for _, row in filtered_df.iterrows():
                                        ensembl_id = row.get('human_ensembl_id')
                                        gene_name = row.get('gene')
                                        if pd.notna(ensembl_id) and pd.notna(gene_name):
                                            ensembl_to_gene[str(ensembl_id)] = str(gene_name)
                                
                                # Run GO:BP enrichment
                                gobp_results = perform_functional_enrichment(gene_list, sources=['GO:BP'])
                                gobp_filtered = filter_enrichment_results(gobp_results, 'GO')
                    
                                # Run KEGG enrichment
                                kegg_results = perform_functional_enrichment(gene_list, sources=['KEGG'])
                                kegg_filtered = filter_enrichment_results(kegg_results, 'KEGG')
                    
                                # Display results in an expander
                                with st.expander("Functional Enrichment Results", expanded=True):
                                    col1, col2 = st.columns(2)
                        
                                    with col1:
                                        st.markdown("### GO Biological Processes")
                                        if not gobp_filtered.empty:
                                            st.write(f"Found {len(gobp_filtered)} significant terms")
                                            # Create and display plot
                                            fig_gobp = create_enrichment_barplot(gobp_filtered, "GO:BP Enrichment", color='skyblue', max_terms=15)
                                            if fig_gobp:
                                                st.pyplot(fig_gobp)
                                            # Display table with term_id and intersections
                                            # Work with head(20) from the start to ensure consistent indexing
                                            gobp_head = gobp_filtered.head(20).copy()
                                            
                                            # Check what columns gprofiler actually returned
                                            available_cols = list(gobp_head.columns)
                                            
                                            # Find term_id column - check all possible names
                                            term_id_col = None
                                            for col in available_cols:
                                                if col.lower() in ['native', 'native_id', 'term_id', 'native_term_id', 'native_term']:
                                                    term_id_col = col
                                                    break
                                            
                                            # If not found by name, check if any column contains term IDs (GO: or hsa:)
                                            if not term_id_col and not gobp_head.empty:
                                                for col in available_cols:
                                                    try:
                                                        sample_val = str(gobp_head[col].iloc[0])
                                                        if sample_val.startswith('GO:') or sample_val.startswith('hsa:'):
                                                            term_id_col = col
                                                            break
                                                    except:
                                                        continue
                                            
                                            # Find intersections column - check all possible names
                                            intersections_col = None
                                            for col in available_cols:
                                                if col.lower() in ['intersections', 'intersection', 'evidence', 'intersection_gene_names', 'intersection_genes']:
                                                    intersections_col = col
                                                    break
                                            # Also check in original results
                                            if not intersections_col:
                                                for col in gobp_results.columns:
                                                    if col.lower() in ['intersections', 'intersection', 'evidence', 'intersection_gene_names', 'intersection_genes']:
                                                        intersections_col = col
                                                        break
                                            
                                            # Build display dataframe with ALL 5 columns from the start
                                            n_rows = len(gobp_head)
                                            
                                            # Create Term ID column
                                            if term_id_col and term_id_col in gobp_head.columns:
                                                term_id_data = gobp_head[term_id_col].reset_index(drop=True).tolist()
                                            else:
                                                term_id_data = ['N/A'] * n_rows
                                            
                                            # Create Genes column
                                            if intersections_col:
                                                source_df = gobp_head if intersections_col in gobp_head.columns else gobp_results.head(20)
                                                if intersections_col in source_df.columns:
                                                    source_df_reset = source_df.reset_index(drop=True)
                                                    genes_data = source_df_reset[intersections_col].apply(
                                                        lambda x: convert_intersections_to_gene_names(x, ensembl_to_gene)
                                                    ).tolist()
                                                    # Apply text wrapping to genes
                                                    genes_data = [wrap_gene_names(g, genes_per_line=8) for g in genes_data]
                                                else:
                                                    genes_data = ['N/A'] * n_rows
                                            else:
                                                genes_data = ['N/A'] * n_rows
                                            
                                            # Get base columns
                                            name_data = gobp_head['name'].reset_index(drop=True).tolist()
                                            pvalue_data = gobp_head['p_value'].reset_index(drop=True).apply(format_scientific_notation).tolist()
                                            intersection_size_data = gobp_head['intersection_size'].reset_index(drop=True).tolist()
                                            
                                            # Build final dataframe with all 5 columns
                                            display_df = pd.DataFrame({
                                                'Term ID': term_id_data,
                                                'Term Name': name_data,
                                                'p-value': pvalue_data,
                                                'Genes in Term': intersection_size_data,
                                                'Genes': genes_data
                                            })
                                            
                                            # Display with column configuration for text wrapping
                                            st.dataframe(
                                                display_df,
                                                use_container_width=True,
                                                column_config={
                                                    "Genes": st.column_config.TextColumn(
                                                        "Genes",
                                                        help="Genes in the enriched term",
                                                        width="large"
                                                    )
                                                },
                                                hide_index=True
                                            )
                                        else:
                                            st.info("No significant GO:BP terms found")
                        
                                    with col2:
                                        st.markdown("### KEGG Pathways")
                                        if not kegg_filtered.empty:
                                            st.write(f"Found {len(kegg_filtered)} significant pathways")
                                            # Create and display plot
                                            fig_kegg = create_enrichment_barplot(kegg_filtered, "KEGG Enrichment", color='lightcoral', max_terms=15)
                                            if fig_kegg:
                                                st.pyplot(fig_kegg)
                                            # Display table with term_id and intersections
                                            # Work with head(20) from the start to ensure consistent indexing
                                            kegg_head = kegg_filtered.head(20).copy()
                                            
                                            # Check what columns gprofiler actually returned
                                            available_cols = list(kegg_head.columns)
                                            
                                            # Find term_id column - check all possible names
                                            term_id_col = None
                                            for col in available_cols:
                                                if col.lower() in ['native', 'native_id', 'term_id', 'native_term_id', 'native_term']:
                                                    term_id_col = col
                                                    break
                                            
                                            # If not found by name, check if any column contains term IDs (GO: or hsa:)
                                            if not term_id_col and not kegg_head.empty:
                                                for col in available_cols:
                                                    try:
                                                        sample_val = str(kegg_head[col].iloc[0])
                                                        if sample_val.startswith('GO:') or sample_val.startswith('hsa:'):
                                                            term_id_col = col
                                                            break
                                                    except:
                                                        continue
                                            
                                            # Find intersections column - check all possible names
                                            intersections_col = None
                                            for col in available_cols:
                                                if col.lower() in ['intersections', 'intersection', 'evidence', 'intersection_gene_names', 'intersection_genes']:
                                                    intersections_col = col
                                                    break
                                            # Also check in original results
                                            if not intersections_col:
                                                for col in kegg_results.columns:
                                                    if col.lower() in ['intersections', 'intersection', 'evidence', 'intersection_gene_names', 'intersection_genes']:
                                                        intersections_col = col
                                                        break
                                            
                                            # Build display dataframe with ALL 5 columns from the start
                                            n_rows = len(kegg_head)
                                            
                                            # Create Term ID column
                                            if term_id_col and term_id_col in kegg_head.columns:
                                                term_id_data = kegg_head[term_id_col].reset_index(drop=True).tolist()
                                            else:
                                                term_id_data = ['N/A'] * n_rows
                                            
                                            # Create Genes column
                                            if intersections_col:
                                                source_df = kegg_head if intersections_col in kegg_head.columns else kegg_results.head(20)
                                                if intersections_col in source_df.columns:
                                                    source_df_reset = source_df.reset_index(drop=True)
                                                    genes_data = source_df_reset[intersections_col].apply(
                                                        lambda x: convert_intersections_to_gene_names(x, ensembl_to_gene)
                                                    ).tolist()
                                                    # Apply text wrapping to genes
                                                    genes_data = [wrap_gene_names(g, genes_per_line=8) for g in genes_data]
                                                else:
                                                    genes_data = ['N/A'] * n_rows
                                            else:
                                                genes_data = ['N/A'] * n_rows
                                            
                                            # Get base columns
                                            name_data = kegg_head['name'].reset_index(drop=True).tolist()
                                            pvalue_data = kegg_head['p_value'].reset_index(drop=True).apply(format_scientific_notation).tolist()
                                            intersection_size_data = kegg_head['intersection_size'].reset_index(drop=True).tolist()
                                            
                                            # Build final dataframe with all 5 columns
                                            display_df = pd.DataFrame({
                                                'Term ID': term_id_data,
                                                'Term Name': name_data,
                                                'p-value': pvalue_data,
                                                'Genes in Term': intersection_size_data,
                                                'Genes': genes_data
                                            })
                                            
                                            # Display with column configuration for text wrapping
                                            st.dataframe(
                                                display_df,
                                                use_container_width=True,
                                                column_config={
                                                    "Genes": st.column_config.TextColumn(
                                                        "Genes",
                                                        help="Genes in the enriched term",
                                                        width="large"
                                                    )
                                                },
                                                hide_index=True
                                            )
                                        else:
                                            st.info("No significant KEGG pathways found")
        
                    # =============================================================================
                    # 4. RUN KE ENRICHMENT
                    # =============================================================================
                    st.markdown("---")
                    st.subheader("Enrich for Key Events")
        
                    # Get unique human ENSGs for enrichment
                    degs = set(filtered_df["human_ensembl_id"].dropna())
        
                    if len(degs) == 0:
                        st.warning("⚠️ No DEGs found with the current filters. Try relaxing the cutoff values.")
                    else:
                        if st.button("Run KE Enrichment", key=f"{key_prefix}_run_ke_enrichment", type="primary"):
                            with st.spinner("Running KE enrichment analysis..."):
                                # Perform KE enrichment
                                res_df = perform_ke_enrichment(degs, ke_map, background_genes)
                    
                                if not res_df.empty:
                                    # Apply FDR correction
                                    res_df = apply_fdr_correction(res_df, alpha=0.05, method="fdr_bh")
                        
                                    # Filter for significant results
                                    significant_df = filter_significant_kes(res_df, fdr_threshold=0.05)
                        
                                    # Store results in TAB-SPECIFIC session state
                                    st.session_state[f"{key_prefix}_ke_results"] = res_df
                                    st.session_state[f"{key_prefix}_significant_kes"] = significant_df
                                    st.session_state[f"{key_prefix}_filtered_degs"] = filtered_df
                                else:
                                    st.warning("No KE enrichment results found.")
            
                        # Display results if they exist (check TAB-SPECIFIC session state)
                        if f"{key_prefix}_significant_kes" in st.session_state and st.session_state[f"{key_prefix}_significant_kes"] is not None and not st.session_state[f"{key_prefix}_significant_kes"].empty:
                            significant_df = st.session_state[f"{key_prefix}_significant_kes"]
                            filtered_df = st.session_state[f"{key_prefix}_filtered_degs"]
                
                            #st.markdown("---")
                            st.subheader("Results")
                            st.success(f"Found {len(significant_df)} significant KEs (FDR < 0.05)")
                
                            # Format and display results table
                            df_display_clean = format_ke_results_for_display(significant_df)
                            st.dataframe(df_display_clean, use_container_width=True)
                            
                            # Collect data for PDF generation
                            pdf_ke_data_list = []
                            summary_table_data = []
                
                            # Collect gene data for visualizations
                            all_ke_gene_data = {}
                
                            for idx, row in significant_df.iterrows():
                                ke_name = row["KE name"]
                                ke_id = row["KE"]
                                overlapping_ensembl = row["Overlapping DEGs List"]
                    
                                gene_details = []
                                for ensembl_id in overlapping_ensembl:
                                    gene_row = filtered_df[filtered_df["human_ensembl_id"] == ensembl_id]
                                    if not gene_row.empty:
                                        gene_info = {
                                            "Ensembl ID": ensembl_id,
                                            "log2FoldChange": gene_row.iloc[0]["log2FoldChange"],
                                            "padj": gene_row.iloc[0]["padj"]
                                        }
                                        if "gene" in gene_row.columns and pd.notna(gene_row.iloc[0]["gene"]):
                                            gene_info["Gene Name"] = gene_row.iloc[0]["gene"]
                                        else:
                                            gene_info["Gene Name"] = ensembl_id
                            
                                        gene_details.append(gene_info)
                    
                                if gene_details:
                                    all_ke_gene_data[ke_id] = {
                                        'ke_name': ke_name,
                                        'gene_details': gene_details
                                    }
                                    
                                    # Prepare data for PDF
                                    viz_data = pd.DataFrame(gene_details)
                                    viz_data = viz_data.sort_values("log2FoldChange", ascending=False)
                                    gene_names = viz_data['Gene Name'].tolist() if 'Gene Name' in viz_data.columns else [f"Gene {i}" for i in range(len(viz_data))]
                                    log2fc_values = viz_data['log2FoldChange'].tolist()
                                    
                                    pdf_ke_data_list.append({
                                        'ke_id': ke_id,
                                        'ke_name': ke_name,
                                        'ke_row': row.to_dict(),  # Convert Series to dict for PDF generation
                                        'gene_details': gene_details,
                                        'gene_names': gene_names,
                                        'log2fc_values': log2fc_values
                                    })

                                    summary_table_data.append({
                                        'KE': ke_id,
                                        'KE name': ke_name,
                                        'DEGs in KE': row.get('DEGs in KE', 0),
                                        'Percent covered': f"{row.get('Percent of KE covered', 0):.1f}%",
                                        'Odds Ratio': f"{row.get('Odds ratio', 0):.2f}",
                                        'adjusted p-value': f"{row.get('adjusted p-value', 0):.2e}"
                                    })
                            
                            # Store PDF data in session state
                            st.session_state[f"{key_prefix}_pdf_data"] = pdf_ke_data_list
                            st.session_state[f"{key_prefix}_summary_table"] = summary_table_data
                            
                            # Add download buttons
                            if pdf_ke_data_list:
                                st.markdown("---")
                                col_download1, col_download2, col_download3 = st.columns([2, 1, 1])
                                with col_download1:
                                    st.markdown("**Download Results**")
                                with col_download2:
                                    pdf_filename = f"KE_Enrichment_{analysis_name or f'Analysis_{analysis_num}'}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.pdf"
                                    
                                    # Generate PDF when button is clicked
                                    if st.button("📥 Generate PDF", key=f"{key_prefix}_download_pdf", use_container_width=True):
                                        with st.spinner("Generating PDF..."):
                                            try:
                                                pdf_bytes = generate_ke_pdf(
                                                    pdf_ke_data_list,
                                                    analysis_name=analysis_name or f"Analysis {analysis_num}",
                                                    dataset_name=st.session_state.get(f"{key_prefix}_dataset_label", "Not specified"),
                                                    sheet_name=st.session_state.get(f"{key_prefix}_dataset_sheet", "Not specified"),
                                                    summary_table=summary_table_data,
                                                    fdr_threshold=0.05
                                                )
                                                st.session_state[f"{key_prefix}_pdf_bytes"] = pdf_bytes
                                                st.session_state[f"{key_prefix}_pdf_filename"] = pdf_filename
                                                st.success("PDF generated successfully!")
                                            except Exception as e:
                                                st.error(f"Error generating PDF: {str(e)}")
                                    
                                    # Show download button if PDF is ready
                                    if f"{key_prefix}_pdf_bytes" in st.session_state:
                                        st.download_button(
                                            label="⬇️ Download PDF",
                                            data=st.session_state[f"{key_prefix}_pdf_bytes"],
                                            file_name=st.session_state[f"{key_prefix}_pdf_filename"],
                                            mime="application/pdf",
                                            key=f"{key_prefix}_pdf_download_btn"
                                        )
                                
                                with col_download3:
                                    html_filename = f"KE_Enrichment_{analysis_name or f'Analysis_{analysis_num}'}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.html"
                                    
                                    # Generate HTML when button is clicked
                                    if st.button("🌐 Generate HTML", key=f"{key_prefix}_download_html", use_container_width=True):
                                        with st.spinner("Generating HTML report..."):
                                            try:
                                                # Get functional enrichment data from session state if available
                                                fe_data = st.session_state.get(f"{key_prefix}_functional_enrichment", {})
                                                
                                                html_content = generate_ke_html_report(
                                                    pdf_ke_data_list,
                                                    analysis_name=analysis_name or f"Analysis {analysis_num}",
                                                    dataset_name=st.session_state.get(f"{key_prefix}_dataset_label", "Not specified"),
                                                    sheet_name=st.session_state.get(f"{key_prefix}_dataset_sheet", "Not specified"),
                                                    summary_table=summary_table_data,
                                                    fdr_threshold=0.05,
                                                    functional_enrichment_data=fe_data
                                                )
                                                st.session_state[f"{key_prefix}_html_content"] = html_content
                                                st.session_state[f"{key_prefix}_html_filename"] = html_filename
                                                st.success("HTML report generated successfully!")
                                            except Exception as e:
                                                st.error(f"Error generating HTML: {str(e)}")
                                    
                                    # Show download button if HTML is ready
                                    if f"{key_prefix}_html_content" in st.session_state:
                                        st.download_button(
                                            label="⬇️ Download HTML",
                                            data=st.session_state[f"{key_prefix}_html_content"],
                                            file_name=st.session_state[f"{key_prefix}_html_filename"],
                                            mime="text/html",
                                            key=f"{key_prefix}_html_download_btn"
                                        )
                                st.markdown("---")
                
                            # Display visualizations
                            #st.markdown("---")
                            #st.subheader("View Detailed Results")
                
                            for ke_id, data in all_ke_gene_data.items():
                                ke_name = data['ke_name']
                                gene_details = data['gene_details']
                                ke_row = significant_df[significant_df['KE'] == ke_id].iloc[0]
                    
                                viz_data = pd.DataFrame(gene_details)
                                viz_data = viz_data.sort_values("log2FoldChange", ascending=False)
                    
                                st.markdown("---")
                                st.markdown(f"##### {ke_name} ({ke_id}, {ke_row['AOP']})")
                    
                                gene_names = viz_data['Gene Name'].tolist() if 'Gene Name' in viz_data.columns else [f"Gene {i}" for i in range(len(viz_data))]
                                log2fc_values = viz_data['log2FoldChange'].tolist()
                    
                                # Calculate height
                                n_genes = len(gene_names)
                                if n_genes > 35:
                                    fig_height = max(n_genes * 0.15, 10)
                                else:
                                    fig_height = max(n_genes * (4/18), 3)
                    
                                col_left, col_right = st.columns([1, 1])
                    
                                with col_left:
                                    st.markdown("<br>", unsafe_allow_html=True)
                        
                                    fig, ax = plt.subplots(figsize=(9, fig_height))
                                    colors = ['#FF6B6B' if x > 0 else '#4ECDC4' for x in log2fc_values]
                        
                                    ax.barh(gene_names, log2fc_values, color=colors, alpha=1)
                                    ax.set_xlabel('log2 Fold Change', fontsize=12)
                                    ax.set_title(f'{ke_name} ({ke_id})', fontsize=13)
                                    ax.axvline(x=0, color='black', linestyle='-', linewidth=0.5)
                                    ax.grid(axis='x', alpha=0.3)
                        
                                    plt.tight_layout()
                                    st.pyplot(fig)
                                    plt.close()
                    
                                with col_right:
                                    st.markdown("##### Key Event Information")
                        
                                    info_html = f"""
                                    <div style="background-color: #f0f2f6; padding: 20px; border-radius: 10px; border-left: 5px solid #4ECDC4;">
                                        <p style="margin: 5px 0;"><strong>KE ID:</strong> {ke_id}</p>
                                        <p style="margin: 5px 0;"><strong>KE Name:</strong> {ke_name}</p>
                                        <p style="margin: 5px 0;"><strong>AOP:</strong> {ke_row['AOP']}</p>
                                        <p style="margin: 5px 0;"><strong>DEGs in KE:</strong> {ke_row['DEGs in KE']}</p>
                                        <p style="margin: 5px 0;"><strong>KE Size:</strong> {ke_row['KE size']}</p>
                                        <p style="margin: 5px 0;"><strong>Percent Covered:</strong> {ke_row['Percent of KE covered']:.1f}%</p>
                                        <p style="margin: 5px 0;"><strong>Adjusted p-value:</strong> {ke_row['adjusted p-value']:.2e}</p>
                                        <p style="margin: 5px 0;"><strong>Odds Ratio:</strong> {ke_row['Odds ratio']:.2f}</p>
                                    </div>
                                    """
                                    st.markdown(info_html, unsafe_allow_html=True)
                    
                                # Gene details expander
                                with st.expander(f"View DEGs in KE: {ke_name}", expanded=False):
                                    # Format gene details for display
                                    gene_display_df = pd.DataFrame(gene_details)
                                    if 'Gene Name' in gene_display_df.columns:
                                        gene_display_df = gene_display_df[['Gene Name', 'Ensembl ID', 'log2FoldChange', 'padj']]
                                        gene_display_df['log2FoldChange'] = gene_display_df['log2FoldChange'].apply(lambda x: f"{x:.3f}")
                                        gene_display_df['padj'] = gene_display_df['padj'].apply(format_scientific_notation)
                                        gene_display_df = gene_display_df.sort_values('log2FoldChange', key=lambda x: x.astype(float).abs(), ascending=False)
                                    st.dataframe(gene_display_df, use_container_width=True, hide_index=True)
                    
                                # Functional enrichment button
                                if st.button(f"Run Functional Enrichment", key=f"{key_prefix}_enrich_ke_{ke_id}"):
                                    with st.spinner("Running functional enrichment..."):
                                        ke_genes_df = pd.DataFrame(gene_details)
                                        gene_col = get_gene_name_column(ke_genes_df)
                                        ke_gene_list = []
                                        if gene_col:
                                            ke_gene_list = ke_genes_df[gene_col].dropna().tolist()
                            
                                        if ke_gene_list:
                                            # Create mapping from Ensembl IDs to gene names for intersection conversion
                                            ke_ensembl_to_gene = {}
                                            for gene_info in gene_details:
                                                ensembl_id = gene_info.get('Ensembl ID')
                                                gene_name = gene_info.get('Gene Name')
                                                if ensembl_id and gene_name:
                                                    ke_ensembl_to_gene[str(ensembl_id)] = str(gene_name)
                                            
                                            gobp_ke = perform_functional_enrichment(ke_gene_list, sources=['GO:BP'])
                                            kegg_ke = perform_functional_enrichment(ke_gene_list, sources=['KEGG'])
                                
                                            gobp_ke_filtered = filter_enrichment_results(gobp_ke, 'GO')
                                            kegg_ke_filtered = filter_enrichment_results(kegg_ke, 'KEGG')
                                            
                                            # Store raw filtered data for barplot generation
                                            if f"{key_prefix}_functional_enrichment" not in st.session_state:
                                                st.session_state[f"{key_prefix}_functional_enrichment"] = {}
                                            if ke_id not in st.session_state[f"{key_prefix}_functional_enrichment"]:
                                                st.session_state[f"{key_prefix}_functional_enrichment"][ke_id] = {}
                                            st.session_state[f"{key_prefix}_functional_enrichment"][ke_id]['GO:BP_raw'] = gobp_ke_filtered
                                            st.session_state[f"{key_prefix}_functional_enrichment"][ke_id]['KEGG_raw'] = kegg_ke_filtered
                                
                                            with st.expander("Enrichment Results", expanded=True):
                                                col_gobp, col_kegg = st.columns(2)
                                    
                                                with col_gobp:
                                                    st.markdown("**GO:BP**")
                                                    if not gobp_ke_filtered.empty:
                                                        fig_gobp_ke = create_enrichment_barplot(gobp_ke_filtered, f"GO:BP - DEGs in {ke_name}", color='skyblue', max_terms=10)
                                                        if fig_gobp_ke:
                                                            st.pyplot(fig_gobp_ke)
                                                        
                                                        # Build display dataframe with ALL 5 columns
                                                        gobp_ke_head = gobp_ke_filtered.head(10).copy()
                                                        n_rows = len(gobp_ke_head)
                                                        
                                                        # Find term_id and intersections columns
                                                        ke_term_id_col = None
                                                        for col in gobp_ke_head.columns:
                                                            if col.lower() in ['native', 'native_id', 'term_id', 'native_term_id', 'native_term']:
                                                                ke_term_id_col = col
                                                                break
                                                        if not ke_term_id_col and not gobp_ke_head.empty:
                                                            for col in gobp_ke_head.columns:
                                                                try:
                                                                    sample_val = str(gobp_ke_head[col].iloc[0])
                                                                    if sample_val.startswith('GO:') or sample_val.startswith('hsa:'):
                                                                        ke_term_id_col = col
                                                                        break
                                                                except:
                                                                    continue
                                                        
                                                        ke_intersections_col = None
                                                        for col in gobp_ke_head.columns:
                                                            if col.lower() in ['intersections', 'intersection', 'evidence', 'intersection_gene_names', 'intersection_genes']:
                                                                ke_intersections_col = col
                                                                break
                                                        if not ke_intersections_col:
                                                            for col in gobp_ke.columns:
                                                                if col.lower() in ['intersections', 'intersection', 'evidence', 'intersection_gene_names', 'intersection_genes']:
                                                                    ke_intersections_col = col
                                                                    break
                                                        
                                                        # Create all 5 columns
                                                        if ke_term_id_col and ke_term_id_col in gobp_ke_head.columns:
                                                            term_id_data = gobp_ke_head[ke_term_id_col].reset_index(drop=True).tolist()
                                                        else:
                                                            term_id_data = ['N/A'] * n_rows
                                                        
                                                        if ke_intersections_col:
                                                            source_df = gobp_ke_head if ke_intersections_col in gobp_ke_head.columns else gobp_ke.head(10)
                                                            if ke_intersections_col in source_df.columns:
                                                                source_df_reset = source_df.reset_index(drop=True)
                                                                genes_data = source_df_reset[ke_intersections_col].apply(
                                                                    lambda x: convert_intersections_to_gene_names(x, ke_ensembl_to_gene)
                                                                ).tolist()
                                                                # Apply text wrapping to genes
                                                                genes_data = [wrap_gene_names(g, genes_per_line=8) for g in genes_data]
                                                            else:
                                                                genes_data = ['N/A'] * n_rows
                                                        else:
                                                            genes_data = ['N/A'] * n_rows
                                                        
                                                        name_data = gobp_ke_head['name'].reset_index(drop=True).tolist()
                                                        pvalue_data = gobp_ke_head['p_value'].reset_index(drop=True).apply(format_scientific_notation).tolist()
                                                        intersection_size_data = gobp_ke_head['intersection_size'].reset_index(drop=True).tolist()
                                                        
                                                        display_gobp = pd.DataFrame({
                                                            'Term ID': term_id_data,
                                                            'Term Name': name_data,
                                                            'p-value': pvalue_data,
                                                            'Genes in Term': intersection_size_data,
                                                            'Genes': genes_data
                                                        })
                                                        
                                                        # Store for HTML report (remove wrapping for storage)
                                                        if f"{key_prefix}_functional_enrichment" not in st.session_state:
                                                            st.session_state[f"{key_prefix}_functional_enrichment"] = {}
                                                        if ke_id not in st.session_state[f"{key_prefix}_functional_enrichment"]:
                                                            st.session_state[f"{key_prefix}_functional_enrichment"][ke_id] = {}
                                                        
                                                        # Store without wrapping for HTML (HTML will handle wrapping)
                                                        gobp_for_html = gobp_ke_head.copy()
                                                        if ke_intersections_col and ke_intersections_col in gobp_for_html.columns:
                                                            gobp_for_html['Genes'] = gobp_for_html[ke_intersections_col].apply(
                                                                lambda x: convert_intersections_to_gene_names(x, ke_ensembl_to_gene)
                                                            )
                                                        else:
                                                            gobp_for_html['Genes'] = 'N/A'
                                                        
                                                        if ke_term_id_col and ke_term_id_col in gobp_for_html.columns:
                                                            gobp_for_html['Term ID'] = gobp_for_html[ke_term_id_col]
                                                        else:
                                                            gobp_for_html['Term ID'] = 'N/A'
                                                        
                                                        gobp_html_df = pd.DataFrame({
                                                            'Term ID': gobp_for_html.get('Term ID', 'N/A'),
                                                            'Term Name': gobp_for_html['name'],
                                                            'p-value': gobp_for_html['p_value'].apply(format_scientific_notation),
                                                            'Genes in Term': gobp_for_html['intersection_size'],
                                                            'Genes': gobp_for_html['Genes']
                                                        }).head(10)
                                                        
                                                        st.session_state[f"{key_prefix}_functional_enrichment"][ke_id]['GO:BP'] = gobp_html_df
                                                        
                                                        # Display with column configuration for text wrapping
                                                        st.dataframe(
                                                            display_gobp,
                                                            use_container_width=True,
                                                            column_config={
                                                                "Genes": st.column_config.TextColumn(
                                                                    "Genes",
                                                                    help="Genes in the enriched term",
                                                                    width="large"
                                                                )
                                                            },
                                                            hide_index=True
                                                        )
                                                    else:
                                                        st.info("No significant terms")
                                    
                                                with col_kegg:
                                                    st.markdown("**KEGG**")
                                                    if not kegg_ke_filtered.empty:
                                                        fig_kegg_ke = create_enrichment_barplot(kegg_ke_filtered, f"KEGG - DEGs in {ke_name}", color='lightcoral', max_terms=10)
                                                        if fig_kegg_ke:
                                                            st.pyplot(fig_kegg_ke)
                                                        
                                                        # Build display dataframe with ALL 5 columns
                                                        kegg_ke_head = kegg_ke_filtered.head(10).copy()
                                                        n_rows = len(kegg_ke_head)
                                                        
                                                        # Find term_id and intersections columns
                                                        ke_term_id_col = None
                                                        for col in kegg_ke_head.columns:
                                                            if col.lower() in ['native', 'native_id', 'term_id', 'native_term_id', 'native_term']:
                                                                ke_term_id_col = col
                                                                break
                                                        if not ke_term_id_col and not kegg_ke_head.empty:
                                                            for col in kegg_ke_head.columns:
                                                                try:
                                                                    sample_val = str(kegg_ke_head[col].iloc[0])
                                                                    if sample_val.startswith('GO:') or sample_val.startswith('hsa:'):
                                                                        ke_term_id_col = col
                                                                        break
                                                                except:
                                                                    continue
                                                        
                                                        ke_intersections_col = None
                                                        for col in kegg_ke_head.columns:
                                                            if col.lower() in ['intersections', 'intersection', 'evidence', 'intersection_gene_names', 'intersection_genes']:
                                                                ke_intersections_col = col
                                                                break
                                                        if not ke_intersections_col:
                                                            for col in kegg_ke.columns:
                                                                if col.lower() in ['intersections', 'intersection', 'evidence', 'intersection_gene_names', 'intersection_genes']:
                                                                    ke_intersections_col = col
                                                                    break
                                                        
                                                        # Create all 5 columns
                                                        if ke_term_id_col and ke_term_id_col in kegg_ke_head.columns:
                                                            term_id_data = kegg_ke_head[ke_term_id_col].reset_index(drop=True).tolist()
                                                        else:
                                                            term_id_data = ['N/A'] * n_rows
                                                        
                                                        if ke_intersections_col:
                                                            source_df = kegg_ke_head if ke_intersections_col in kegg_ke_head.columns else kegg_ke.head(10)
                                                            if ke_intersections_col in source_df.columns:
                                                                source_df_reset = source_df.reset_index(drop=True)
                                                                genes_data = source_df_reset[ke_intersections_col].apply(
                                                                    lambda x: convert_intersections_to_gene_names(x, ke_ensembl_to_gene)
                                                                ).tolist()
                                                                # Apply text wrapping to genes
                                                                genes_data = [wrap_gene_names(g, genes_per_line=8) for g in genes_data]
                                                            else:
                                                                genes_data = ['N/A'] * n_rows
                                                        else:
                                                            genes_data = ['N/A'] * n_rows
                                                        
                                                        name_data = kegg_ke_head['name'].reset_index(drop=True).tolist()
                                                        pvalue_data = kegg_ke_head['p_value'].reset_index(drop=True).apply(format_scientific_notation).tolist()
                                                        intersection_size_data = kegg_ke_head['intersection_size'].reset_index(drop=True).tolist()
                                                        
                                                        display_kegg = pd.DataFrame({
                                                            'Term ID': term_id_data,
                                                            'Term Name': name_data,
                                                            'p-value': pvalue_data,
                                                            'Genes in Term': intersection_size_data,
                                                            'Genes': genes_data
                                                        })
                                                        
                                                        # Store for HTML report
                                                        if f"{key_prefix}_functional_enrichment" not in st.session_state:
                                                            st.session_state[f"{key_prefix}_functional_enrichment"] = {}
                                                        if ke_id not in st.session_state[f"{key_prefix}_functional_enrichment"]:
                                                            st.session_state[f"{key_prefix}_functional_enrichment"][ke_id] = {}
                                                        
                                                        # Store without wrapping for HTML
                                                        kegg_for_html = kegg_ke_head.copy()
                                                        if ke_intersections_col and ke_intersections_col in kegg_for_html.columns:
                                                            kegg_for_html['Genes'] = kegg_for_html[ke_intersections_col].apply(
                                                                lambda x: convert_intersections_to_gene_names(x, ke_ensembl_to_gene)
                                                            )
                                                        else:
                                                            kegg_for_html['Genes'] = 'N/A'
                                                        
                                                        if ke_term_id_col and ke_term_id_col in kegg_for_html.columns:
                                                            kegg_for_html['Term ID'] = kegg_for_html[ke_term_id_col]
                                                        else:
                                                            kegg_for_html['Term ID'] = 'N/A'
                                                        
                                                        kegg_html_df = pd.DataFrame({
                                                            'Term ID': kegg_for_html.get('Term ID', 'N/A'),
                                                            'Term Name': kegg_for_html['name'],
                                                            'p-value': kegg_for_html['p_value'].apply(format_scientific_notation),
                                                            'Genes in Term': kegg_for_html['intersection_size'],
                                                            'Genes': kegg_for_html['Genes']
                                                        }).head(10)
                                                        
                                                        st.session_state[f"{key_prefix}_functional_enrichment"][ke_id]['KEGG'] = kegg_html_df
                                                        
                                                        # Display with column configuration for text wrapping
                                                        st.dataframe(
                                                            display_kegg,
                                                            use_container_width=True,
                                                            column_config={
                                                                "Genes": st.column_config.TextColumn(
                                                                    "Genes",
                                                                    help="Genes in the enriched term",
                                                                    width="large"
                                                                )
                                                            },
                                                            hide_index=True
                                                        )
                                                    else:
                                                        st.info("No significant pathways")
                                        else:
                                            st.warning("⚠️ Gene names not available for enrichment")
            
            except Exception as e:
                st.error(f"❌ Error reading file: {e}")
                st.info("Make sure your file is in the correct format (CSV or Excel)")
        else:
            st.info("ℹ️ Upload a file or select example data to get started.")

