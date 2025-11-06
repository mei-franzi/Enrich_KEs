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
st.title("📊 Key Event Enrichment")

st.markdown("""
This tool performs enrichment analysis on your differentially expressed genes (DEGs) against genes associated with Key Events (KEs). 

**Key Events** are measurable biological events within Adverse Outcome Pathways (AOPs)—multi-scale models that connect molecular initiating events to adverse health outcomes. By mapping your DEGs to the AOP Key Event database, this analysis identifies which KEs are statistically overrepresented in your dataset using Fisher's exact test with Benjamini-Hochberg correction (FDR < 0.05).

This approach embeds the AOP framework into molecular data interpretation, enabling novel AOP-based approaches in biomedical research and supporting the development of new approach methodologies (NAMs) that reduce reliance on animal experimentation.
""")

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
st.subheader("📁 Upload Your Data")

uploaded_file = st.file_uploader(
    "Upload your differential expression results (CSV or Excel format)",
    type=["csv", "xlsx", "xls"],
    help="File should contain columns: padj, log2FoldChange, human_ensembl_id"
)

if uploaded_file is not None:
    st.success("✅ File uploaded successfully!")
    
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
    
    st.write(f"Loaded file shape: {deg_file_data.shape[0]} rows, {deg_file_data.shape[1]} columns")
    st.write("First few rows of your data:")
    st.dataframe(deg_file_data.head(), use_container_width=True)
else:
    st.info("ℹ️ Upload a file to get started, or use the default data below.")
    deg_file_data = None

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
    st.subheader("🔍 Running Enrichment Analysis on Your Data")
    df = deg_file_data.copy()
    sheet_names = ["Uploaded Data"]
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
            (df["padj"] < 0.05) & 
            (df["log2FoldChange"].abs() > 0.1) & 
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
        st.dataframe(df, use_container_width=True)
else:
    st.warning("No results to display. Please check your data files.")