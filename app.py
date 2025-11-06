# libraries 
import pandas as pd  # data manipulation and analysis
import os  # operating system interface for file paths
import numpy as np # numerical operations and arrays
import scipy
import statsmodels
from scipy.stats import fisher_exact  # statistical tests
from statsmodels.stats.multitest import multipletests 
import importlib.metadata  # package version information
import streamlit as st
import plotly.graph_objects as go
import matplotlib.pyplot as plt
from enrichment import perform_functional_enrichment, filter_enrichment_results, create_enrichment_barplot

st.set_page_config(layout="wide")
st.title("Key Event and Functional Enrichment of DEGs")

# Description in a box
st.markdown("""
<div style="border: 2px solid white; padding: 20px; border-radius: 10px;">
    <h3 style="margin-top: 0;">Perform enrichment analyses on your DEGs against Key Events and Functional Pathways</h3>
    <p><strong>Key Events</strong> are measurable biological events within Adverse Outcome Pathways (AOPs)—multi-scale models that connect molecular initiating events to adverse health outcomes.</p>
    <p>By mapping your DEGs to the AOP Key Event database, KEs that are statistically overrepresented in your dataset using Fisher's exact test and corrected for false discovery with Benjamini-Hochberg correction.</p>
</div>
""", unsafe_allow_html=True)

# File paths - using relative paths for easier deployment
# Create a 'data' folder in your project directory and place your files there
data_dir = "data"

ke_map_path = os.path.join(data_dir, "Genes_to_KEs.txt")
ke_desc_path = os.path.join(data_dir, "ke_descriptions.csv")

# Check if reference files exist
if not os.path.exists(ke_map_path):
    st.error(f"⚠️ Missing file: {ke_map_path}")
    st.info("Please create a 'data' folder in your project directory and add the required reference files:")
    st.write("- Genes_to_KEs.txt")
    st.write("- ke_descriptions.csv")
    st.stop()
if not os.path.exists(ke_desc_path):
    st.error(f"⚠️ Missing file: {ke_desc_path}")
    st.stop()

# =============================================================================
# KEY EVENT ENRICHMENT 
# =============================================================================

# Key Events identification through Fisher's exact test enrichment analysis 
# of differentially expressed genes
# apply Multiple testing correction using the Benjamini-Hochberg FDR method, with FDR < 0.05. 

st.markdown("---")
st.subheader("Upload Your Data")

uploaded_file = st.file_uploader(
    "Upload your differential expression results (csv, tsv or excel)",
    type=["csv", "tsv", "xlsx", "xls"],
    help="File should contain: adjusted p-values, log2 fold change, and human Ensembl IDs. You can map your column names after upload."
)

if uploaded_file is not None:
    st.success("File upload successful")
    
    try:
        # Load the uploaded file
        if uploaded_file.name.endswith('.csv'):
            # Try to detect the delimiter (comma or semicolon)
            deg_file_data = pd.read_csv(uploaded_file, sep=None, engine='python')
        else:
            deg_file_data = pd.read_excel(uploaded_file)
        
        # Column mapping section
        st.markdown("---")
        st.subheader("Column Mapping")
        st.write("Specify which columns in your file correspond to the required data:")
        
        col_map1, col_map2, col_map3 = st.columns(3)
        
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
        
        with col_map3:
            ensembl_col = st.selectbox(
                "Human Ensembl ID column",
                options=deg_file_data.columns.tolist(),
                index=deg_file_data.columns.tolist().index('human_ensembl_id') if 'human_ensembl_id' in deg_file_data.columns else 0,
                help="Column containing human Ensembl gene IDs"
            )
        
        # Rename columns to standard names for internal processing
        # Only rename if the column names are different from the standard names
        rename_dict = {}
        if padj_col != 'padj':
            rename_dict[padj_col] = 'padj'
        if log2fc_col != 'log2FoldChange':
            rename_dict[log2fc_col] = 'log2FoldChange'
        if ensembl_col != 'human_ensembl_id':
            rename_dict[ensembl_col] = 'human_ensembl_id'
        
        if rename_dict:
            deg_file_data = deg_file_data.rename(columns=rename_dict)
        
    except Exception as e:
        st.error(f"❌ Error reading file: {e}")
        st.info("Make sure your file is in the correct format (CSV or Excel)")
        deg_file_data = None
    
    # Analysis parameters
    st.markdown("---")
    st.subheader("DEG Filtering Parameters for KE Enrichment")
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        experiment_name = st.text_input(
            "Experiment Name",
            value="My Experiment",
            help="Name your analysis"
        )
    
    with col2:
        padj_cutoff = st.number_input(
            "Adjusted p-value cutoff (DEGs)",
            min_value=0.0,
            max_value=1.0,
            value=0.05,
            step=0.01,
            format="%.3f",
            help="Threshold for statistical significance (default: 0.05)"
        )
    
    with col3:
        log2fc_cutoff = st.number_input(
            "Log2 Fold Change cutoff",
            min_value=0.0,
            max_value=10.0,
            value=0.1,
            step=0.1,
            format="%.2f",
            help="Minimum absolute log2 fold change (default: 0.1)"
        )
    
    # Apply filtering to preview dataframe
    st.markdown("---")
    st.write(f"**Loaded file:** {deg_file_data.shape[0]} DEGs, {deg_file_data.shape[1]} columns")
    
    # Filter the data based on user parameters
    filtered_df = deg_file_data.copy()
    
    # Apply filters if the required columns exist
    if 'padj' in filtered_df.columns and 'log2FoldChange' in filtered_df.columns:
        filtered_df = filtered_df[
            (filtered_df["padj"] < padj_cutoff) & 
            (filtered_df["log2FoldChange"].abs() > log2fc_cutoff)
        ]
        st.write(f"**Filtered DEGs:** {filtered_df.shape[0]} genes (padj < {padj_cutoff:.3g}, |log2FC| > {log2fc_cutoff:.3g})")
    
    # Format preview dataframe
    preview_df = filtered_df.copy()
    
    # Format p-value columns in scientific notation if they exist
    if 'padj' in preview_df.columns:
        preview_df['padj'] = preview_df['padj'].apply(lambda x: f"{x:.2e}" if pd.notna(x) else x)
    if 'pvalue' in preview_df.columns:
        preview_df['pvalue'] = preview_df['pvalue'].apply(lambda x: f"{x:.2e}" if pd.notna(x) else x)
    
    # Display filtered dataframe
    st.dataframe(preview_df, use_container_width=True, height=400)
    
    # ENRICHMENT #1: Functional enrichment on all filtered DEGs
    st.markdown("---")
    st.subheader("Optional: Functional Enrichment on your DEGs")
    
    if st.button("Run Enrichment", key="enrich_all_degs"):
        with st.spinner("Running functional enrichment analysis..."):
            # Check if gene names are available
            if "gene" in filtered_df.columns or "gene_name" in filtered_df.columns or "symbol" in filtered_df.columns:
                gene_col = "gene" if "gene" in filtered_df.columns else ("gene_name" if "gene_name" in filtered_df.columns else "symbol")
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
                    st.markdown("**GO Biological Processes**")
                    if not gobp_filtered.empty:
                        st.write(f"Found {len(gobp_filtered)} significant terms")
                        # Create and display plot
                        fig_gobp = create_enrichment_barplot(gobp_filtered, "GO:BP Enrichment", color='skyblue', max_terms=15)
                        if fig_gobp:
                            st.pyplot(fig_gobp)
                        # Display table
                        display_df = gobp_filtered[['name', 'p_value', 'intersection_size', 'term_size']].head(20)
                        display_df['p_value'] = display_df['p_value'].apply(lambda x: f"{x:.2e}")
                        st.dataframe(display_df, use_container_width=True)
                    else:
                        st.info("No significant GO:BP terms found")
                
                with col2:
                    st.markdown("**KEGG Pathways**")
                    if not kegg_filtered.empty:
                        st.write(f"Found {len(kegg_filtered)} significant pathways")
                        # Create and display plot
                        fig_kegg = create_enrichment_barplot(kegg_filtered, "KEGG Enrichment", color='lightcoral', max_terms=15)
                        if fig_kegg:
                            st.pyplot(fig_kegg)
                        # Display table
                        display_df = kegg_filtered[['name', 'p_value', 'intersection_size', 'term_size']].head(20)
                        display_df['p_value'] = display_df['p_value'].apply(lambda x: f"{x:.2e}")
                        st.dataframe(display_df, use_container_width=True)
                    else:
                        st.info("No significant KEGG pathways found")

else:
    st.info("ℹ️ Upload a file to get started.")
    deg_file_data = None
    experiment_name = "Default Data"
    padj_cutoff = 0.05
    log2fc_cutoff = 0.1


# Load KE mapping and descriptions
ke_map = pd.read_csv(ke_map_path, sep="\t").dropna(subset=["Gene", "KE"])
ke_desc = pd.read_csv(ke_desc_path)
ke_map = ke_map.merge(ke_desc, on="KE", how="left")
background_genes = set(ke_map["Gene"].unique())

# Run enrichment analysis only if file is uploaded
if deg_file_data is not None:
    # Use uploaded file
    st.markdown("---")
    st.subheader(f"KE-Enrichment Analysis: {experiment_name}")
    
    # Second set of filtering parameters for easy access
    st.markdown("**Adjust filtering parameters**")
    
    # Use a narrower column for the inputs
    col_narrow, col_spacer = st.columns([1, 2])
    
    with col_narrow:
        padj_cutoff = st.number_input(
            "Adjusted p-value (KE results)",
            min_value=0.0,
            max_value=1.0,
            value=padj_cutoff,
            step=0.01,
            format="%.3f",
            help="Threshold for statistical significance (default: 0.05)",
            key="padj_cutoff_2"
        )
        
        log2fc_cutoff = st.number_input(
            "Log2 Fold Change cutoff of DEGs",
            min_value=0.0,
            max_value=10.0,
            value=log2fc_cutoff,
            step=0.1,
            format="%.2f",
            help="Minimum absolute log2 fold change (default: 0.1)",
            key="log2fc_cutoff_2"
        )
        
    all_results = {}
    
    try:
        # Load DEGs with human orthologs
        df = deg_file_data.copy()
        
        # Filter for significant DEGs with valid human Ensembl IDs that start with "ENS"
        df = df[
            (df["padj"] < padj_cutoff) & 
            (df["log2FoldChange"].abs() > log2fc_cutoff) & 
            (df["human_ensembl_id"].notna()) & 
            (df["human_ensembl_id"].astype(str).str.startswith("ENS"))
        ].copy()
        
        # Get unique human ENSGs for enrichment
        degs = set(df["human_ensembl_id"].dropna())
        
        results = []
        
        # Perform Fisher's exact test for each KE
        for ke_id, group in ke_map.groupby("KE"):
            ke_genes = set(group["Gene"])
            overlap = degs & ke_genes
            
            if len(overlap) == 0:  # Skip KEs with no overlap
                continue
                
            # Build contingency table and perform test
            a = len(overlap)  # in DEG and in KE
            b = len(degs - ke_genes)  # in DEG, not in KE
            c = len(ke_genes - degs)  # not in DEG, in KE
            d = len(background_genes - degs - ke_genes)  # not in DEG and not in KE
            
            odds_ratio, p = fisher_exact([[a, b], [c, d]], alternative="greater")
            
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
                    "DEGs in KE": a,
                    "KE size": len(ke_genes),
                    "Percent of KE covered": (a/len(ke_genes)*100),
                    "Overlapping DEGs": ", ".join(sorted(overlap)),
                    "Overlapping DEGs List": list(sorted(overlap)),  # Keep as list for later use
                    "Odds ratio": odds_ratio,
                    "p-value": p
                })
        
        if results:
            # Apply multiple testing correction and sort
            res_df = pd.DataFrame(results)
            res_df["adjusted p-value"] = multipletests(res_df["p-value"], method="fdr_bh")[1]
            res_df = res_df.sort_values("adjusted p-value")
            
            # Store significant results (FDR < 0.05)
            significant_df = res_df[res_df["adjusted p-value"] < 0.05].copy()
            
            # Display results
            st.markdown("---")
            st.subheader("Results")
            
            if not significant_df.empty:
                st.success(f"Found {len(significant_df)} significant KEs")
                
                # Format the dataframe for display
                df_display = significant_df.copy()
                
                # Format p-values in scientific notation
                df_display["p-value"] = df_display["p-value"].apply(lambda x: f"{x:.2e}")
                df_display["adjusted p-value"] = df_display["adjusted p-value"].apply(lambda x: f"{x:.2e}")
                
                # Format other numeric columns
                df_display["Percent of KE covered"] = df_display["Percent of KE covered"].apply(lambda x: f"{x:.1f}%")
                df_display["Odds ratio"] = df_display["Odds ratio"].apply(lambda x: f"{x:.2f}")
                
                # Remove the internal list column from display
                df_display_clean = df_display.drop(columns=["Overlapping DEGs List"])
                
                st.dataframe(df_display_clean, use_container_width=True)
                
                # Add expandable sections for each KE showing gene details
                st.markdown("---")
                st.subheader("Detailed Gene Information for Each KE")
                
                for idx, row in significant_df.iterrows():
                    ke_name = row["KE name"]
                    ke_id = row["KE"]
                    overlapping_ensembl = row["Overlapping DEGs List"]
                    
                    with st.expander(f"{ke_id} {ke_name} - {len(overlapping_ensembl)} genes"):
                        # Heading with KE name
                        st.markdown(f"### {ke_name}")
                        
                        # Get gene details from the original DEG data
                        gene_details = []
                        for ensembl_id in overlapping_ensembl:
                            # Find this gene in the original filtered DEG data
                            gene_row = df[df["human_ensembl_id"] == ensembl_id]
                            if not gene_row.empty:
                                gene_info = {
                                    "Ensembl ID": ensembl_id,
                                    "log2FoldChange": gene_row.iloc[0]["log2FoldChange"],
                                    "padj": gene_row.iloc[0]["padj"]
                                }
                                # Check if there's a gene name column
                                if "gene" in gene_row.columns:
                                    gene_info["Gene Name"] = gene_row.iloc[0]["gene"]
                                elif "gene_name" in gene_row.columns:
                                    gene_info["Gene Name"] = gene_row.iloc[0]["gene_name"]
                                elif "symbol" in gene_row.columns:
                                    gene_info["Gene Name"] = gene_row.iloc[0]["symbol"]
                                
                                gene_details.append(gene_info)
                        
                        if gene_details:
                            gene_df = pd.DataFrame(gene_details)
                            # Rename for display
                            gene_df = gene_df.rename(columns={"log2FoldChange": "log2FC"})
                            # Format log2FC to 2 decimal places
                            gene_df["log2FC"] = gene_df["log2FC"].apply(lambda x: f"{x:.2f}" if pd.notna(x) else x)
                            # Format padj in scientific notation
                            gene_df["padj"] = gene_df["padj"].apply(lambda x: f"{x:.2e}" if pd.notna(x) else x)
                            # Reorder columns: Ensembl ID, Gene Name, log2FC, padj
                            if "Gene Name" in gene_df.columns:
                                gene_df = gene_df[["Ensembl ID", "Gene Name", "log2FC", "padj"]]
                            else:
                                gene_df = gene_df[["Ensembl ID", "log2FC", "padj"]]
                            # Sort by absolute log2FC (convert back to float for sorting)
                            gene_df["log2FC_numeric"] = gene_df["log2FC"].astype(float)
                            gene_df = gene_df.sort_values("log2FC_numeric", key=abs, ascending=False)
                            gene_df = gene_df.drop(columns=["log2FC_numeric"])
                            
                            # Display table
                            st.dataframe(gene_df, hide_index=True)
                            
                            # Prepare data for heatmap - sorted by log2FC (not absolute)
                            heatmap_data = pd.DataFrame(gene_details)
                            heatmap_data = heatmap_data.sort_values("log2FoldChange", ascending=False)
                            
                            # Get gene labels (use Gene Name if available, else Ensembl ID)
                            if "Gene Name" in heatmap_data.columns:
                                gene_labels = heatmap_data["Gene Name"].tolist()
                            else:
                                gene_labels = heatmap_data["Ensembl ID"].tolist()
                            
                            # Get log2FC values
                            log2fc_values = heatmap_data["log2FoldChange"].values
                            
                            # Create horizontal Plotly heatmap (genes on x-axis)
                            max_abs = max(abs(log2fc_values.min()), abs(log2fc_values.max()))
                            
                            fig = go.Figure(data=go.Heatmap(
                                z=log2fc_values.reshape(1, -1),  # Reshape for horizontal layout
                                x=gene_labels,  # Genes on x-axis
                                y=["log2FC"],  # Single row
                                colorscale="RdBu_r",  # Red for high, Blue for low
                                zmid=0,
                                zmin=-max_abs,
                                zmax=max_abs,
                                text=[[f"{val:.2f}" for val in log2fc_values]],
                                texttemplate="%{text}",
                                textfont={"size": 14},  # Cell values size 14
                                hovertemplate="Gene: %{x}<br>log2FC: %{z:.2f}<extra></extra>",
                                showscale=False
                            ))
                            
                            # Update layout for horizontal display with squared cells
                            n_genes = len(gene_labels)
                            cell_size = 60  # Squared cells (60px x 60px)
                            fig.update_layout(
                                height=cell_size + 100,  # Cell height + margin for labels
                                width=max(n_genes * cell_size, 300),  # 60px per gene
                                margin=dict(l=5, r=5, t=5, b=100),
                                xaxis=dict(side="bottom", tickfont=dict(size=14), tickangle=-45),  # Bigger x-axis labels
                                yaxis=dict(tickfont=dict(size=14)),
                                font=dict(size=11)
                            )
                            
                            # Display the heatmap in Streamlit with unique key
                            st.plotly_chart(fig, use_container_width=True, key=f"heatmap_{ke_id}")
                            
                            # ENRICHMENT #2: Functional enrichment on this KE's genes
                            st.markdown("---")
                            
                            if st.button(f"Run Functional Enrichment ({ke_name})", key=f"enrich_ke_{ke_id}"):
                                with st.spinner("Running functional enrichment..."):
                                    # Get gene names for this KE
                                    ke_gene_list = []
                                    if "Gene Name" in heatmap_data.columns:
                                        ke_gene_list = heatmap_data["Gene Name"].dropna().tolist()
                                    
                                    if ke_gene_list:
                                        # Run enrichment
                                        gobp_ke = perform_functional_enrichment(ke_gene_list, sources=['GO:BP'])
                                        kegg_ke = perform_functional_enrichment(ke_gene_list, sources=['KEGG'])
                                        
                                        gobp_ke_filtered = filter_enrichment_results(gobp_ke, 'GO')
                                        kegg_ke_filtered = filter_enrichment_results(kegg_ke, 'KEGG')
                                        
                                        # Display in columns with plots
                                        col_gobp, col_kegg = st.columns(2)
                                        
                                        with col_gobp:
                                            st.markdown("**GO:BP**")
                                            if not gobp_ke_filtered.empty:
                                                # Create and display plot
                                                fig_gobp_ke = create_enrichment_barplot(gobp_ke_filtered, f"GO:BP - {ke_name}", color='skyblue', max_terms=10)
                                                if fig_gobp_ke:
                                                    st.pyplot(fig_gobp_ke)
                                                # Display table
                                                display_gobp = gobp_ke_filtered[['name', 'p_value', 'intersection_size']].head(10)
                                                display_gobp['p_value'] = display_gobp['p_value'].apply(lambda x: f"{x:.2e}")
                                                st.dataframe(display_gobp, use_container_width=True, hide_index=True)
                                            else:
                                                st.info("No significant terms")
                                        
                                        with col_kegg:
                                            st.markdown("**KEGG**")
                                            if not kegg_ke_filtered.empty:
                                                # Create and display plot
                                                fig_kegg_ke = create_enrichment_barplot(kegg_ke_filtered, f"KEGG - {ke_name}", color='lightcoral', max_terms=10)
                                                if fig_kegg_ke:
                                                    st.pyplot(fig_kegg_ke)
                                                # Display table
                                                display_kegg = kegg_ke_filtered[['name', 'p_value', 'intersection_size']].head(10)
                                                display_kegg['p_value'] = display_kegg['p_value'].apply(lambda x: f"{x:.2e}")
                                                st.dataframe(display_kegg, use_container_width=True, hide_index=True)
                                            else:
                                                st.info("No significant pathways")
                                    else:
                                        st.warning("⚠️ Gene names not available for enrichment")
                # ENRICHMENT #3: Functional enrichment on union of all KE genes
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
                                    # Create and display plot
                                    fig_gobp_all = create_enrichment_barplot(gobp_all_filtered, "GO:BP - All KE Genes", color='skyblue', max_terms=15)
                                    if fig_gobp_all:
                                        st.pyplot(fig_gobp_all)
                                    # Display table
                                    display_df = gobp_all_filtered[['name', 'p_value', 'intersection_size', 'term_size']].head(20)
                                    display_df['p_value'] = display_df['p_value'].apply(lambda x: f"{x:.2e}")
                                    st.dataframe(display_df, use_container_width=True)
                                else:
                                    st.info("No significant GO:BP terms found")
                            
                            with col2_all:
                                st.markdown("**KEGG Pathways**")
                                if not kegg_all_filtered.empty:
                                    st.write(f"Found {len(kegg_all_filtered)} significant pathways")
                                    # Create and display plot
                                    fig_kegg_all = create_enrichment_barplot(kegg_all_filtered, "KEGG - All KE Genes", color='lightcoral', max_terms=15)
                                    if fig_kegg_all:
                                        st.pyplot(fig_kegg_all)
                                    # Display table
                                    display_df = kegg_all_filtered[['name', 'p_value', 'intersection_size', 'term_size']].head(20)
                                    display_df['p_value'] = display_df['p_value'].apply(lambda x: f"{x:.2e}")
                                    st.dataframe(display_df, use_container_width=True)
                                else:
                                    st.info("No significant KEGG pathways found")
                        else:
                            st.warning("⚠️ Gene names not available for enrichment")
            else:
                st.info(f"ℹ️ No significant enrichment found (FDR < 0.05)")
        else:
            st.warning(f"⚠️ No enrichment results - no overlapping genes found")
            
    except Exception as e:
        st.error(f"❌ Error processing data: {e}")