# Libraries 
import pandas as pd
import os
import streamlit as st
import matplotlib.pyplot as plt

# Project modules
from enrichment import perform_functional_enrichment, filter_enrichment_results, create_enrichment_barplot
from ke_enrichment import (
    perform_ke_enrichment,
    apply_fdr_correction,
    filter_significant_kes,
    format_ke_results_for_display
)
from data_loader import load_deg_file, load_deg_from_path, prepare_ke_data, apply_column_mapping, filter_degs, get_excel_sheet_names
from utils import format_scientific_notation, get_gene_name_column

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

st.title("🧬 Key Event Enrichment Analysis")
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
                                            # Display table
                                            display_df = gobp_filtered[['name', 'p_value', 'intersection_size', 'term_size']].head(20)
                                            display_df['p_value'] = display_df['p_value'].apply(format_scientific_notation)
                                            st.dataframe(display_df, use_container_width=True)
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
                                            # Display table
                                            display_df = kegg_filtered[['name', 'p_value', 'intersection_size', 'term_size']].head(20)
                                            display_df['p_value'] = display_df['p_value'].apply(format_scientific_notation)
                                            st.dataframe(display_df, use_container_width=True)
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
                        
                                    # Store results in session state
                                    st.session_state.ke_results = res_df
                                    st.session_state.significant_kes = significant_df
                                    st.session_state.filtered_degs = filtered_df
                                else:
                                    st.warning("No KE enrichment results found.")
            
                        # Display results if they exist
                        if 'significant_kes' in st.session_state and st.session_state.significant_kes is not None and not st.session_state.significant_kes.empty:
                            significant_df = st.session_state.significant_kes
                            filtered_df = st.session_state.filtered_degs
                
                            #st.markdown("---")
                            st.subheader("Results")
                            st.success(f"Found {len(significant_df)} significant KEs (FDR < 0.05)")
                
                            # Format and display results table
                            df_display_clean = format_ke_results_for_display(significant_df)
                            st.dataframe(df_display_clean, use_container_width=True)
                
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
                                st.markdown(f"##### {ke_name} ({ke_id})")
                    
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
                                            gobp_ke = perform_functional_enrichment(ke_gene_list, sources=['GO:BP'])
                                            kegg_ke = perform_functional_enrichment(ke_gene_list, sources=['KEGG'])
                                
                                            gobp_ke_filtered = filter_enrichment_results(gobp_ke, 'GO')
                                            kegg_ke_filtered = filter_enrichment_results(kegg_ke, 'KEGG')
                                
                                            with st.expander("Enrichment Results", expanded=True):
                                                col_gobp, col_kegg = st.columns(2)
                                    
                                                with col_gobp:
                                                    st.markdown("**GO:BP**")
                                                    if not gobp_ke_filtered.empty:
                                                        fig_gobp_ke = create_enrichment_barplot(gobp_ke_filtered, f"GO:BP - {ke_name}", color='skyblue', max_terms=10)
                                                        if fig_gobp_ke:
                                                            st.pyplot(fig_gobp_ke)
                                                        display_gobp = gobp_ke_filtered[['name', 'p_value', 'intersection_size']].head(10)
                                                        display_gobp['p_value'] = display_gobp['p_value'].apply(format_scientific_notation)
                                                        st.dataframe(display_gobp, use_container_width=True, hide_index=True)
                                                    else:
                                                        st.info("No significant terms")
                                    
                                                with col_kegg:
                                                    st.markdown("**KEGG**")
                                                    if not kegg_ke_filtered.empty:
                                                        fig_kegg_ke = create_enrichment_barplot(kegg_ke_filtered, f"KEGG - {ke_name}", color='lightcoral', max_terms=10)
                                                        if fig_kegg_ke:
                                                            st.pyplot(fig_kegg_ke)
                                                        display_kegg = kegg_ke_filtered[['name', 'p_value', 'intersection_size']].head(10)
                                                        display_kegg['p_value'] = display_kegg['p_value'].apply(format_scientific_notation)
                                                        st.dataframe(display_kegg, use_container_width=True, hide_index=True)
                                                    else:
                                                        st.info("No significant pathways")
                                        else:
                                            st.warning("⚠️ Gene names not available for enrichment")
            
            except Exception as e:
                st.error(f"❌ Error reading file: {e}")
                st.info("Make sure your file is in the correct format (CSV or Excel)")
        else:
            st.info("ℹ️ Upload a file or select example data to get started.")

