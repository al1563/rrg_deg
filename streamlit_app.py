import streamlit as st
import pandas as pd
import plotly.express as px
import numpy as np

# Set the page configuration for a wider layout
st.set_page_config(layout="wide")

@st.cache_data # Use Streamlit's caching to load data only once
def load_data():
    """Loads all DGE and GSEA files and merges them."""
    
    # Function to safely read a CSV file
    def read_data(path):
        try:
            return pd.read_csv(path, low_memory=False)
        except FileNotFoundError:
            st.warning(f"File not found: {path}")
            return None

    # Load all DGE files
    dge_files = {
        "cd45pos_rrg": "data/cd45pos_rrgcell_degs_wilcoxon.csv",
        "cd45pos_jeff": "data/cd45pos_jeffcell_degs_wilcoxon.csv",
        "cd45neg_rrg": "data/cd45neg_rrgcell_degs_wilcoxon.csv"
    }
    dge_list = []
    for name, path in dge_files.items():
        df = read_data(path)
        if df is not None:
            df['dataset'] = name
            dge_list.append(df)
    all_dge = pd.concat(dge_list, ignore_index=True).rename(columns={'cell': 'cell_type'})

    # Load all GSEA files
    gsea_files = {
        "cd45pos_rrg": "data/cd45pos_rrgcell_gsea.csv",
        "cd45pos_jeff": "data/cd45pos_jeffcell_gsea.csv",
        "cd45neg_rrg": "data/cd45neg_rrgcell_gsea.csv"
    }
    gsea_list = []
    for name, path in gsea_files.items():
        df = read_data(path)
        if df is not None:
            df['dataset'] = name
            gsea_list.append(df)
    all_gsea = pd.concat(gsea_list, ignore_index=True).rename(columns={'Term': 'pathway', 'FDR q-val': 'padj'})
    
    return all_dge, all_gsea

all_dge, all_gsea = load_data()

st.title("Interactive Single-Cell Analysis Explorer")

# --- Sidebar for Controls ---
with st.sidebar:
    st.header("Data Selection")

    main_group = st.selectbox(
        "Select Main Group:",
        options=[("CD45 Positive", "cd45pos"), ("CD45 Negative", "cd45neg")],
        format_func=lambda x: x[0]
    )[1]

    # Conditional UI for annotation type
    if main_group == "cd45pos":
        annotation_type = st.selectbox(
            "Select Annotation:",
            options=[("RRG Cell", "rrg"), ("Jeff Cell", "jeff")],
            format_func=lambda x: x[0]
        )[1]
        selected_dataset = f"cd45pos_{annotation_type}"
    else:
        selected_dataset = "cd45neg_rrg"

    # Filter data based on selected dataset
    dge_subset = all_dge[all_dge['dataset'] == selected_dataset]
    
    # Dynamic selectbox for cell type
    cell_type = st.selectbox(
        "Select Cell Type:",
        options=sorted(dge_subset['cell_type'].unique())
    )
    
    # Further subset data for the chosen cell type
    dge_cell_type_subset = dge_subset[dge_subset['cell_type'] == cell_type]

    # Dynamic selectboxes for comparisons
    comp1 = st.selectbox(
        "Select Comparison Group 1:",
        options=sorted(dge_cell_type_subset['comp1'].unique())
    )
    
    comp2 = st.selectbox(
        "Select Comparison Group 2 (Reference):",
        options=sorted(dge_cell_type_subset[dge_cell_type_subset['comp1'] == comp1]['comp2'].unique())
    )

    st.markdown("---")
    st.header("Volcano Plot Controls")

    logfc_cutoff = st.slider(
        "Log2 Fold Change Cutoff:",
        min_value=0.0, max_value=4.0, value=0.5, step=0.1
    )
    
    pval_cutoff_log = st.slider(
        "Adjusted P-value Cutoff (-log10):",
        min_value=0.0, max_value=50.0, value=-np.log10(0.05), step=0.5
    )
    pval_cutoff = 10**-pval_cutoff_log

    st.markdown("---")
    st.header("GSEA Plot Controls")
    
    gsea_pathway_count = st.slider(
        "Number of Pathways to Display:",
        min_value=5, max_value=50, value=20, step=1
    )

# --- Main Panel with Tabs ---
tab1, tab2 = st.tabs(["DGE Volcano View", "GSEA Pathway View"])

# --- Tab 1: DGE Volcano Plot and Table ---
with tab1:
    # Filter data based on all selections
    filtered_dge = dge_cell_type_subset[
        (dge_cell_type_subset['comp1'] == comp1) & 
        (dge_cell_type_subset['comp2'] == comp2)
    ]
    
    if filtered_dge.empty:
        st.warning("No data available for the current selection.")
    else:
        # Define significance based on cutoffs
        filtered_dge['significant'] = (
            (filtered_dge['p_val_adj'] < pval_cutoff) & 
            (abs(filtered_dge['avg_log2FC']) > logfc_cutoff)
        ).map({True: 'Significant', False: 'Not Significant'})


        st.subheader("Volcano Plot")
        fig_volcano = px.scatter(
            filtered_dge,
            x="avg_log2FC",
            y=-np.log10(filtered_dge['p_val_adj']+min(filtered_dge['p_val_adj'][filtered_dge['p_val_adj'] > 0])),
            color="significant",
            color_discrete_map={"Significant": "#d62728", "Not Significant": "grey"},
            hover_name="gene",
            hover_data={"avg_log2FC": ':.3f', "p_val_adj": ':.3e'},
            labels={"y": "-log10(Adjusted P-value)", "x": "Log2 Fold Change"}
        )
        fig_volcano.add_vline(x=logfc_cutoff, line_dash="dash", line_color="grey")
        fig_volcano.add_vline(x=-logfc_cutoff, line_dash="dash", line_color="grey")
        fig_volcano.add_hline(y=pval_cutoff_log, line_dash="dash", line_color="grey")
        fig_volcano.update_layout(showlegend=False, height=600)
        st.plotly_chart(fig_volcano, width=True)

        st.subheader("DGE Results Table")
        st.dataframe(
            filtered_dge[['gene', 'avg_log2FC', 'p_val_adj']].round(3),
            width=True
        )

# --- Tab 2: GSEA Pathway View and Table ---
with tab2:
    # Filter GSEA data
    filtered_gsea = all_gsea[
        (all_gsea['dataset'] == selected_dataset) &
        (all_gsea['cell_type'] == cell_type) &
        (all_gsea['comp1'] == comp1) &
        (all_gsea['comp2'] == comp2)
    ].sort_values('padj').head(gsea_pathway_count)

    if filtered_gsea.empty:
        st.warning("No GSEA data available for the current selection.")
    else:
        st.subheader("GSEA Pathway Enrichment")
        fig_gsea = px.scatter(
            filtered_gsea,
            x="NES",
            y="pathway",
            color="NES",
            color_continuous_scale="RdBu_r",
            color_continuous_midpoint = 0,
            hover_name="pathway",
            hover_data={"NES": ':.3f', "padj": ':.3e'},
            labels={"y": "Pathway", "x": "Normalized Enrichment Score (NES)"}
        )
        fig_gsea.update_layout(
            yaxis={'categoryorder':'total ascending'},
            height=600
        )
        st.plotly_chart(fig_gsea, width=True)

        st.subheader("GSEA Results Table")
        st.dataframe(
            filtered_gsea[['path_name', 'reference', 'NES', 'padj', 'Lead_genes','Tag %','Gene %']].round(3),
            width=True
        )
