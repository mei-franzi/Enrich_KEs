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

st.set_page_config(layout="wide")
st.title("Key Event Enrichment")

# Description in a box
st.markdown("""
<div style="border: 2px solid white; padding: 20px; border-radius: 10px;">
    <h3 style="margin-top: 0;">Perform KE enrichment analysis with your DEGs</h3>
    <p><strong>Key Events</strong> are measurable biological events within Adverse Outcome Pathways (AOPs)—multi-scale models that connect molecular initiating events to adverse health outcomes.</p>
    <p>By mapping your DEGs to the AOP Key Event database, KEs that are statistically overrepresented in your dataset using Fisher's exact test and corrected for false discovery with Benjamini-Hochberg correction.</p>
</div>
""", unsafe_allow_html=True)

# File paths - using relative paths for easier deployment
# Create a 'data' folder in your project directory and place your files there
data_dir = "data"

ke_map_path = os.path.join(data_dir, "Genes_to_KEs.txt")
ke_desc_path = os.path.join(data_dir, "ke_descriptions.csv")
deg_file = os.path.join(data_dir, "Browder_DEGs.xlsx")

# Check if files exist, provide helpful message if not
if not os.path.exists(ke_map_path):
    st.error(f"⚠️ Missing file: {ke_map_path}")
    st.info("Please create a 'data' folder in your project directory and add the required files:")
    st.write("- Genes_to_KEs.txt")
    st.write("- ke_descriptions.csv")
    st.write("- differential_expression_results_with_human_orthologs.xlsx")
    st.stop()
if not os.path.exists(ke_desc_path):
    st.error(f"⚠️ Missing file: {ke_desc_path}")
    st.stop()
if not os.path.exists(deg_file):
    st.error(f"⚠️ Missing file: {deg_file}")
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
    "Upload your differential expression results (CSV, tsv or excel). Required columns: padj, log2FoldChange, human_ensembl_id",
    type=["csv", "tsv", "xlsx", "xls"],
    help="File should contain columns: padj, log2FoldChange, human_ensembl_id"
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
    except Exception as e:
        st.error(f"❌ Error reading file: {e}")
        st.info("Make sure your file has the correct format with columns: padj, log2FoldChange, human_ensembl_id")
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
    st.write(f"**Loaded file:** {deg_file_data.shape[0]} total rows, {deg_file_data.shape[1]} columns")
    
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
else:
    st.info("ℹ️ Upload a file to get started, or view the example results below.")
    deg_file_data = None
    experiment_name = "Default Data"
    padj_cutoff = 0.05
    log2fc_cutoff = 0.1

st.markdown("---")
st.write("Loading reference data...")

# Load KE mapping and descriptions
ke_map = pd.read_csv(ke_map_path, sep="\t").dropna(subset=["Gene", "KE"])
ke_desc = pd.read_csv(ke_desc_path)
ke_map = ke_map.merge(ke_desc, on="KE", how="left")
background_genes = set(ke_map["Gene"].unique())

# Determine which data to use
if deg_file_data is not None:
    # Use uploaded file
    st.subheader(f"🔍 Running Enrichment Analysis: {experiment_name}")
    
    # Second set of filtering parameters for easy access
    st.markdown("**Adjust filtering parameters for enrichment analysis:**")
    col1b, col2b, col3b = st.columns(3)
    
    with col1b:
        experiment_name = st.text_input(
            "Experiment Name",
            value=experiment_name,
            help="Name your analysis",
            key="experiment_name_2"
        )
    
    with col2b:
        padj_cutoff = st.number_input(
            "Adjusted p-value cutoff (KE results)",
            min_value=0.0,
            max_value=1.0,
            value=padj_cutoff,
            step=0.01,
            format="%.3f",
            help="Threshold for statistical significance (default: 0.05)",
            key="padj_cutoff_2"
        )
    
    with col3b:
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
    
    st.markdown("---")
    df = deg_file_data.copy()
    sheet_names = [experiment_name]
else:
    # Use default file
    st.subheader("🔍 Running Enrichment Analysis on Default Data")
    sheet_names = ["Kidney_old_young"]

all_results = {}

for sheet in sheet_names:
    try:
        # Load DEGs with human orthologs
        if deg_file_data is not None:
            df = deg_file_data.copy()
        else:
            df = pd.read_excel(deg_file, sheet_name=sheet)
        
        # Filter for significant DEGs with valid human Ensembl IDs that start with "ENS"
        df = df[
            (df["padj"] < padj_cutoff) & 
            (df["log2FoldChange"].abs() > log2fc_cutoff) & 
            (df["human_ensembl_id"].notna()) & 
            (df["human_ensembl_id"].astype(str).str.startswith("ENS"))
        ].copy()
        
        # Get unique human ENSGs for enrichment
        degs = set(df["human_ensembl_id"].dropna())
        st.info(f"📊 Filtered DEGs: {len(degs)} genes (padj < {padj_cutoff:.3g}, |log2FC| > {log2fc_cutoff:.3g})")
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
            if not significant_df.empty:
                all_results[sheet] = significant_df
                st.success(f"✅ {sheet}: Found {len(significant_df)} significant KEs")
            else:
                st.info(f"ℹ️ {sheet}: No significant enrichment")
        else:
            st.warning(f"⚠️ {sheet}: No enrichment results")
            
    except Exception as e:
        st.error(f"❌ {sheet}: Error processing sheet - {e}")

# Display results
if all_results:
    st.markdown("---")
    st.subheader("Results")
    for sheet, df in all_results.items():
        st.write(f"**{sheet}**")
        
        # Format the dataframe for display
        df_display = df.copy()
        
        # Format p-values in scientific notation
        df_display["p-value"] = df_display["p-value"].apply(lambda x: f"{x:.2e}")
        df_display["adjusted p-value"] = df_display["adjusted p-value"].apply(lambda x: f"{x:.2e}")
        
        # Format other numeric columns
        df_display["Percent of KE covered"] = df_display["Percent of KE covered"].apply(lambda x: f"{x:.1f}%")
        df_display["Odds ratio"] = df_display["Odds ratio"].apply(lambda x: f"{x:.2f}")
        
        st.dataframe(df_display, use_container_width=True)
else:
    st.warning("No results to display. Please check your data files.")