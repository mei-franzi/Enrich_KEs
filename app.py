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
    format_ke_results_for_display,
    create_ke_heatmap
)
from data_loader import load_deg_file, load_deg_from_path, prepare_ke_data, apply_column_mapping, filter_degs, get_excel_sheet_names
from utils import format_scientific_notation, get_gene_name_column

st.set_page_config(layout="wide", page_title="KE & Functional Enrichment")

# =============================================================================
# SIDEBAR
# =============================================================================
with st.sidebar:
    st.title("")


    # Global experiment name
    experiment_name = st.text_input(
        "Experiment Name",
        value="My Experiment",
      #  help="Name your analysis"
    )
    
    # DEG Filtering Parameters (Global)
    st.subheader("DEG Filter Thresholds")
    padj_cutoff = st.number_input(
        "Adjusted p-value cutoff",
        min_value=0.0,
        max_value=1.0,
        value=0.05,
        step=0.01,
        format="%.3f",
        help="can be p-value, padj, FDR, q-value, ..."
    )
    
    log2fc_cutoff = st.number_input(
        "Log2 Fold Change cutoff",
        min_value=0.0,
        max_value=10.0,
        value=0.1,
        step=0.1,
        format="%.2f",
      #  help="Minimum absolute log2 fold change"
    )
    
    st.markdown("---")
    st.caption("Data files expected in `data/` folder")

# =============================================================================
# MAIN CONTENT
# =============================================================================

st.title("DEGs Analyzer")

# 1: KE Enrichment
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

# Initialize session state for data persistence across tabs
if 'deg_file_data' not in st.session_state:
    st.session_state.deg_file_data = None
if 'filtered_df' not in st.session_state:
    st.session_state.filtered_df = None

# =============================================================================
# TABS
# =============================================================================

tab1, tab2, tab3 = st.tabs(["Data Upload & Mapping", "Functional Enrichment", "KE Enrichment"])

# =============================================================================
# TAB 1: Data Upload & Column Mapping
# =============================================================================
with tab1:
    st.subheader("Upload DEGs")
    
    # File upload section
   # st.subheader("")

    uploaded_file = st.file_uploader(
        "Upload differential expression results",
        type=["csv", "tsv", "xlsx", "xls"],
        help="File should contain: adjusted p-values, log2 fold change, and human Ensembl IDs. You can map your column names after upload."
    )
    
    # Example data buttons
    st.write("**Or use example data:**")
    col_ex1, col_ex2 = st.columns(2)
    
    example_data_path = os.path.join(data_dir, "GSE255602_DEGs.csv")
    browder_data_path = os.path.join(data_dir, "Browder_DEGs.xlsx")
    
    use_example = False
    use_browder = False
    
    with col_ex1:
        if os.path.exists(example_data_path):
            if st.button("Use GSE255602_DEGs.csv"):
                use_example = True
                st.session_state.example_type = 'gse'
    
    with col_ex2:
        if os.path.exists(browder_data_path):
            if st.button("Use Browder_DEGs.xlsx"):
                use_browder = True
                st.session_state.example_type = 'browder'
    
    # Determine which file source to use
    file_source = None
    file_name = None
    is_excel = False
    
    if uploaded_file is not None:
        file_source = uploaded_file
        file_name = uploaded_file.name
        is_excel = file_name.endswith(('.xlsx', '.xls'))
    elif use_example or (hasattr(st.session_state, 'example_type') and st.session_state.example_type == 'gse'):
        file_source = example_data_path
        file_name = "GSE255602_DEGs.csv"
        is_excel = False
    elif use_browder or (hasattr(st.session_state, 'example_type') and st.session_state.example_type == 'browder'):
        file_source = browder_data_path
        file_name = "Browder_DEGs.xlsx"
        is_excel = True
    
    # Sheet selection for Excel files
    selected_sheet = None
    if file_source and is_excel:
        st.subheader("Select Excel Sheet")
        
        # Get sheet names
        sheet_names = get_excel_sheet_names(file_source)
        
        if sheet_names:
            selected_sheet = st.selectbox(
                "Choose which sheet to analyze",
                options=sheet_names,
                help="Select the sheet containing your DEG data"
            )
            st.info(f"Selected sheet: **{selected_sheet}**")
        else:
            st.error("Could not read Excel sheets")
            file_source = None
    
    # Load data
    if file_source is not None:
        st.success(f"Using data: {file_name}")
        
        try:
            # Load the file with sheet selection if applicable
            if isinstance(file_source, str):  # File path (example data)
                deg_file_data = load_deg_from_path(file_source, sheet_name=selected_sheet)
            else:  # Uploaded file
                deg_file_data = load_deg_file(file_source, sheet_name=selected_sheet)
            
            if deg_file_data is None:
                st.error("Failed to load file. Please check the format.")
                st.stop()
            
            # Column mapping section
            st.markdown("---")
            st.subheader("2. Column Mapping")
            st.write("Specify which columns in your file correspond to the required data:")
            
            col_map1, col_map2 = st.columns(2)
            
            with col_map1:
                padj_col = st.selectbox(
                    "Adjusted p-value column",
                    options=deg_file_data.columns.tolist(),
                    index=deg_file_data.columns.tolist().index('padj') if 'padj' in deg_file_data.columns else 0,
                    help="Column containing adjusted p-values"
                )
            
            with col_map2:
                log2fc_col = st.selectbox(
                    "Log2 Fold Change column",
                    options=deg_file_data.columns.tolist(),
                    index=deg_file_data.columns.tolist().index('log2FoldChange') if 'log2FoldChange' in deg_file_data.columns else 0,
                    help="Column containing log2 fold change values"
                )
            
            col_map3, col_map4 = st.columns(2)
            
            with col_map3:
                ensembl_col = st.selectbox(
                    "Human Ensembl ID column",
                    options=deg_file_data.columns.tolist(),
                    index=deg_file_data.columns.tolist().index('human_ensembl_id') if 'human_ensembl_id' in deg_file_data.columns else 0,
                    help="Column containing human Ensembl gene IDs"
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
                    help="Column containing gene names or symbols (required for functional enrichment)"
                )
            
            # Rename columns to standard names
            rename_dict = {
                padj_col: 'padj',
                log2fc_col: 'log2FoldChange',
                ensembl_col: 'human_ensembl_id',
                gene_name_col: 'gene'
            }
            deg_file_data = apply_column_mapping(deg_file_data, rename_dict)
            
            # Store in session state
            st.session_state.deg_file_data = deg_file_data
            
            # Display DEGs table
            st.markdown("---")
            st.subheader("3. Differentially Expressed Genes Table")
            
            st.write(f"**Total genes loaded:** {deg_file_data.shape[0]} rows, {deg_file_data.shape[1]} columns")
            
            # Apply filters
            if 'padj' in deg_file_data.columns and 'log2FoldChange' in deg_file_data.columns:
                filtered_df = deg_file_data[
                    (deg_file_data["padj"] < padj_cutoff) & 
                    (deg_file_data["log2FoldChange"].abs() > log2fc_cutoff)
                ]
                st.session_state.filtered_df = filtered_df
                
                st.write(f"**Filtered DEGs:** {filtered_df.shape[0]} genes (padj < {padj_cutoff:.3g}, |log2FC| > {log2fc_cutoff:.3g})")
                
                # Format preview dataframe
                preview_df = filtered_df.copy()
                if 'padj' in preview_df.columns:
                    preview_df['padj'] = preview_df['padj'].apply(format_scientific_notation)
                if 'pvalue' in preview_df.columns:
                    preview_df['pvalue'] = preview_df['pvalue'].apply(format_scientific_notation)
                
                st.dataframe(preview_df, use_container_width=True, height=400)
            else:
                st.warning("Required columns not found for filtering")
                
        except Exception as e:
            st.error(f"❌ Error reading file: {e}")
            st.info("Make sure your file is in the correct format (CSV or Excel)")
    else:
        st.info("ℹ️ Upload a file to get started.")

# =============================================================================
# TAB 2: Functional Enrichment
# =============================================================================
with tab2:
    st.header("Functional Enrichment Analysis")
    
    if st.session_state.filtered_df is not None and not st.session_state.filtered_df.empty:
        filtered_df = st.session_state.filtered_df
        
        st.write(f"**Analyzing {len(filtered_df)} filtered DEGs** (padj < {padj_cutoff:.3g}, |log2FC| > {log2fc_cutoff:.3g})")
        
        st.markdown("---")
        
        if st.button("Run Functional Enrichment Analysis", key="enrich_all_degs", type="primary"):
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
                    
                    # Display results with plots
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
    else:
        st.info("ℹ️ Please upload and configure data in the 'Data Upload & Mapping' tab first.")

# =============================================================================
# TAB 3: KE Enrichment
# =============================================================================
with tab3:
    st.header(f"KE-Enrichment Analysis: {experiment_name}")
    
    if st.session_state.deg_file_data is not None:
        deg_file_data = st.session_state.deg_file_data
        
        # Input parameters section
        st.subheader("Adjust Input Parameters")
        
        col_param1, col_param2 = st.columns(2)
        
        with col_param1:
            ke_padj_cutoff = st.number_input(
                "Adjusted p-value cutoff",
                min_value=0.0,
                max_value=1.0,
                value=padj_cutoff,
                step=0.01,
                format="%.3f",
                help="Threshold for statistical significance",
                key="ke_padj_cutoff"
            )
        
        with col_param2:
            ke_log2fc_cutoff = st.number_input(
                "Log2 Fold Change cutoff",
                min_value=0.0,
                max_value=10.0,
                value=log2fc_cutoff,
                step=0.1,
                format="%.2f",
                help="Minimum absolute log2 fold change",
                key="ke_log2fc_cutoff"
            )
        
        # Filter DEGs based on KE parameters
        try:
            df = deg_file_data.copy()
            df = filter_degs(df, ke_padj_cutoff, ke_log2fc_cutoff)
            
            # Show DEG table used for input in an expander
            with st.expander(f"View DEGs for KE Enrichment ({len(df)} genes)", expanded=False):
                st.write(f"**Filtered with:** padj < {ke_padj_cutoff:.3g}, |log2FC| > {ke_log2fc_cutoff:.3g}")
                
                # Format for display
                preview_ke_df = df.copy()
                if 'padj' in preview_ke_df.columns:
                    preview_ke_df['padj'] = preview_ke_df['padj'].apply(format_scientific_notation)
                if 'pvalue' in preview_ke_df.columns:
                    preview_ke_df['pvalue'] = preview_ke_df['pvalue'].apply(format_scientific_notation)
                
                st.dataframe(preview_ke_df, use_container_width=True, height=300)
            
            st.markdown("---")
            
            # Get unique human ENSGs for enrichment
            degs = set(df["human_ensembl_id"].dropna())
            
            if len(degs) == 0:
                st.warning("⚠️ No DEGs found with the current filters. Try relaxing the cutoff values.")
            else:
                st.info(f"Running KE enrichment on {len(degs)} unique human Ensembl IDs...")
                
                # Perform KE enrichment analysis
                res_df = perform_ke_enrichment(degs, ke_map, background_genes)
                
                if not res_df.empty:
                    # Apply FDR correction
                    res_df = apply_fdr_correction(res_df, alpha=0.05, method="fdr_bh")
                    
                    # Filter for significant results
                    significant_df = filter_significant_kes(res_df, fdr_threshold=0.05)
                    
                    # Display results
                    st.markdown("---")
                    st.subheader("Results")
                    
                    if not significant_df.empty:
                        st.success(f"Found {len(significant_df)} significant KEs (FDR < 0.05)")
                        
                        # Format the dataframe for display
                        df_display_clean = format_ke_results_for_display(significant_df)
                        
                        st.dataframe(df_display_clean, use_container_width=True)
                        
                        # Collect gene data for all KEs
                        all_ke_gene_data = {}
                        
                        for idx, row in significant_df.iterrows():
                            ke_name = row["KE name"]
                            ke_id = row["KE"]
                            overlapping_ensembl = row["Overlapping DEGs List"]
                            
                            # Get gene details from the original DEG data
                            gene_details = []
                            for ensembl_id in overlapping_ensembl:
                                gene_row = df[df["human_ensembl_id"] == ensembl_id]
                                if not gene_row.empty:
                                    gene_info = {
                                        "Ensembl ID": ensembl_id,
                                        "log2FoldChange": gene_row.iloc[0]["log2FoldChange"],
                                        "padj": gene_row.iloc[0]["padj"]
                                    }
                                    # Check if there's a gene name column (now standardized as 'gene')
                                    if "gene" in gene_row.columns and pd.notna(gene_row.iloc[0]["gene"]):
                                        gene_info["Gene Name"] = gene_row.iloc[0]["gene"]
                                    elif "gene_name" in gene_row.columns and pd.notna(gene_row.iloc[0]["gene_name"]):
                                        gene_info["Gene Name"] = gene_row.iloc[0]["gene_name"]
                                    elif "symbol" in gene_row.columns and pd.notna(gene_row.iloc[0]["symbol"]):
                                        gene_info["Gene Name"] = gene_row.iloc[0]["symbol"]
                                    else:
                                        # Fallback to Ensembl ID if no gene name available
                                        gene_info["Gene Name"] = ensembl_id
                                    
                                    gene_details.append(gene_info)
                            
                            # Store gene data for this KE
                            if gene_details:
                                all_ke_gene_data[ke_id] = {
                                    'ke_name': ke_name,
                                    'gene_details': gene_details
                                }
                        
                        # Display gene expression visualizations
                        st.markdown("---")
                        st.subheader("Gene Expression Visualizations for Significant KEs")
                        
                        for ke_id, data in all_ke_gene_data.items():
                            ke_name = data['ke_name']
                            gene_details = data['gene_details']
                            
                            # Prepare data for visualization
                            viz_data = pd.DataFrame(gene_details)
                            viz_data = viz_data.sort_values("log2FoldChange", ascending=False)
                            
                            # Create matplotlib bar chart (more reliable than plotly in tabs)
                            st.markdown(f"**{ke_name}** ({ke_id})")
                            
                            # Get gene names and log2FC values
                            gene_names = viz_data['Gene Name'].tolist() if 'Gene Name' in viz_data.columns else [f"Gene {i}" for i in range(len(viz_data))]
                            log2fc_values = viz_data['log2FoldChange'].tolist()
                            
                            # Calculate height based on number of genes
                            n_genes = len(gene_names)
                            if n_genes > 35:
                                # For large gene sets, use more compact spacing
                                fig_height = max(n_genes * 0.15, 10)  # ~0.15 inches per gene
                            else:
                                # For smaller gene sets, use comfortable spacing (18 genes = 4 inches)
                                fig_height = max(n_genes * (4/18), 3)
                            
                            # Use columns to position plot on the left half
                            col_left, col_right = st.columns([1, 1])
                            
                            with col_left:
                                fig, ax = plt.subplots(figsize=(9, fig_height))
                                
                                # Create bar colors based on up/down regulation
                                colors = ['red' if x > 0 else 'blue' for x in log2fc_values]
                                
                                # Create horizontal bar chart
                                ax.barh(gene_names, log2fc_values, color=colors, alpha=0.7)
                                ax.set_xlabel('log2 Fold Change', fontsize=12)
                                ax.set_title(f'{ke_name}', fontsize=12, fontweight='bold')
                                ax.axvline(x=0, color='black', linestyle='-', linewidth=0.5)
                                ax.grid(axis='x', alpha=0.3)
                                
                                # Adjust layout
                                plt.tight_layout()
                                st.pyplot(fig)
                                plt.close()
                        
                        # NOW display gene detail tables in expanders
                        st.markdown("---")
                        st.subheader("Detailed Gene Information for Each KE")
                        
                        for idx, row in significant_df.iterrows():
                            ke_name = row["KE name"]
                            ke_id = row["KE"]
                            overlapping_ensembl = row["Overlapping DEGs List"]
                            
                            # Get the gene details we collected earlier
                            if ke_id in all_ke_gene_data:
                                gene_details = all_ke_gene_data[ke_id]['gene_details']
                                
                                with st.expander(f"{ke_id} {ke_name} - {len(overlapping_ensembl)} genes"):
                                    # Heading with KE name
                                    st.markdown(f"### {ke_name}")
                                    
                                    # Create and display gene table
                                    gene_df = pd.DataFrame(gene_details)
                                    gene_df = gene_df.rename(columns={"log2FoldChange": "log2FC"})
                                    gene_df["log2FC"] = gene_df["log2FC"].apply(lambda x: f"{x:.2f}" if pd.notna(x) else x)
                                    gene_df["padj"] = gene_df["padj"].apply(lambda x: f"{x:.2e}" if pd.notna(x) else x)
                                    
                                    if "Gene Name" in gene_df.columns:
                                        gene_df = gene_df[["Ensembl ID", "Gene Name", "log2FC", "padj"]]
                                    else:
                                        gene_df = gene_df[["Ensembl ID", "log2FC", "padj"]]
                                    
                                    gene_df["log2FC_numeric"] = gene_df["log2FC"].astype(float)
                                    gene_df = gene_df.sort_values("log2FC_numeric", key=abs, ascending=False)
                                    gene_df = gene_df.drop(columns=["log2FC_numeric"])
                                    
                                    st.dataframe(gene_df, hide_index=True)
                                    
                                    # Functional enrichment on this KE's genes
                                    st.markdown("---")
                                    
                                    if st.button(f"Run Functional Enrichment", key=f"enrich_ke_{ke_id}"):
                                        with st.spinner("Running functional enrichment..."):
                                            # Recreate dataframe for functional enrichment
                                            ke_genes_df = pd.DataFrame(gene_details)
                                            gene_col = get_gene_name_column(ke_genes_df)
                                            ke_gene_list = []
                                            if gene_col:
                                                ke_gene_list = ke_genes_df[gene_col].dropna().tolist()
                                            
                                            if ke_gene_list:
                                                # Run enrichment
                                                gobp_ke = perform_functional_enrichment(ke_gene_list, sources=['GO:BP'])
                                                kegg_ke = perform_functional_enrichment(ke_gene_list, sources=['KEGG'])
                                                
                                                gobp_ke_filtered = filter_enrichment_results(gobp_ke, 'GO')
                                                kegg_ke_filtered = filter_enrichment_results(kegg_ke, 'KEGG')
                                                
                                                # Display in expandable section
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

                        # Functional enrichment on union of all KE genes
                        st.markdown("---")
                        st.subheader("Run Functional Enrichment on KE-associated DEGs")
                        st.write(f"Analyzing the union of all genes from {len(significant_df)} significant KEs")
                        
                        if st.button("Run Enrichment on All KE Genes", key="enrich_all_ke_genes"):
                            with st.spinner("Running functional enrichment on combined KE genes..."):
                                # Collect all unique genes from all significant KEs
                                all_ke_genes = set()
                                for idx, row in significant_df.iterrows():
                                    overlapping_ensembl = row["Overlapping DEGs List"]
                                    for ensembl_id in overlapping_ensembl:
                                        gene_row = df[df["human_ensembl_id"] == ensembl_id]
                                        if not gene_row.empty:
                                            if "gene" in gene_row.columns:
                                                gene_name = gene_row.iloc[0]["gene"]
                                            elif "gene_name" in gene_row.columns:
                                                gene_name = gene_row.iloc[0]["gene_name"]
                                            elif "symbol" in gene_row.columns:
                                                gene_name = gene_row.iloc[0]["symbol"]
                                            else:
                                                continue
                                            if pd.notna(gene_name):
                                                all_ke_genes.add(gene_name)
                                
                                all_ke_gene_list = list(all_ke_genes)
                                st.write(f"Analyzing {len(all_ke_gene_list)} unique genes")
                                
                                if all_ke_gene_list:
                                    # Run enrichment
                                    gobp_all = perform_functional_enrichment(all_ke_gene_list, sources=['GO:BP'])
                                    kegg_all = perform_functional_enrichment(all_ke_gene_list, sources=['KEGG'])
                                    
                                    gobp_all_filtered = filter_enrichment_results(gobp_all, 'GO')
                                    kegg_all_filtered = filter_enrichment_results(kegg_all, 'KEGG')
                                    
                                    # Display results with plots
                                    col1_all, col2_all = st.columns(2)
                                    
                                    with col1_all:
                                        st.markdown("**GO Biological Processes**")
                                        if not gobp_all_filtered.empty:
                                            st.write(f"Found {len(gobp_all_filtered)} significant terms")
                                            fig_gobp_all = create_enrichment_barplot(gobp_all_filtered, "GO:BP - All KE Genes", color='skyblue', max_terms=15)
                                            if fig_gobp_all:
                                                st.pyplot(fig_gobp_all)
                                            display_df = gobp_all_filtered[['name', 'p_value', 'intersection_size', 'term_size']].head(20)
                                            display_df['p_value'] = display_df['p_value'].apply(format_scientific_notation)
                                            st.dataframe(display_df, use_container_width=True)
                                        else:
                                            st.info("No significant GO:BP terms found")
                                    
                                    with col2_all:
                                        st.markdown("**KEGG Pathways**")
                                        if not kegg_all_filtered.empty:
                                            st.write(f"Found {len(kegg_all_filtered)} significant pathways")
                                            fig_kegg_all = create_enrichment_barplot(kegg_all_filtered, "KEGG - All KE Genes", color='lightcoral', max_terms=15)
                                            if fig_kegg_all:
                                                st.pyplot(fig_kegg_all)
                                            display_df = kegg_all_filtered[['name', 'p_value', 'intersection_size', 'term_size']].head(20)
                                            display_df['p_value'] = display_df['p_value'].apply(format_scientific_notation)
                                            st.dataframe(display_df, use_container_width=True)
                                        else:
                                            st.info("No significant KEGG pathways found")
                                else:
                                    st.warning("⚠️ Gene names not available for enrichment")
                    else:
                        st.info("ℹ️ No significant enrichment found (FDR < 0.05)")
                else:
                    st.warning("⚠️ No enrichment results - no overlapping genes found")
                    
        except Exception as e:
            st.error(f"❌ Error processing data: {e}")
    else:
        st.info("ℹ️ Please upload and configure data in the 'Data Upload & Mapping' tab first.")