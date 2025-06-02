import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
from scipy.cluster.hierarchy import linkage, dendrogram
from streamlit_extras.switch_page_button import switch_page
from utils import load_sequence, extract_s_gene, align_sequences, find_mutations

# CONFIG
S_START, S_END = 21563, 25384

variant_paths = {
    "Wuhan": "./wuhan.fasta",
    "Delta": "./delta.fasta",
    "Gamma": "./gamma.fasta",
    "Omicron": "./omicron.fasta"
}

variant_colors = {
    "Wuhan": "#1f77b4",
    "Delta": "#ff7f0e",
    "Gamma": "#2ca02c",
    "Omicron": "#d62728"
}

# STYLE
st.markdown("""
<style>
    [data-testid="stSidebarNav"] {
        display: none;
    }
    .main-header {
        background: linear-gradient(90deg, #667eea 0%, #764ba2 100%);
        padding: 2rem;
        border-radius: 10px;
        margin-bottom: 2rem;
        text-align: center;
        color: white;
    }
    
    .main-header h1 {
        margin: 0;
        font-size: 2.5rem;
        font-weight: 700;
    }
    
    .main-header p {
        margin: 0.5rem 0 0 0;
        font-size: 1.2rem;
        opacity: 0.9;
    }
    
    .stButton > button {
        background: linear-gradient(90deg, #667eea 0%, #764ba2 100%);
        color: white;
        border: none;
        border-radius: 8px;
        padding: 0.5rem 1rem;
        font-weight: 600;
        transition: all 0.3s ease;
    }
    
    .stButton > button:hover {
        box-shadow: 0 4px 8px rgba(0,0,0,0.2);
        transform: translateY(-2px);
    }
    
    .sidebar .stSelectbox > div > div {
        background-color: #f1f3f4;
    }
</style>
""", unsafe_allow_html=True)

# Header
st.markdown("""
<div class="main-header">
    <h1>🧬 Evolution from Wuhan to Variants</h1>
    <p>This section visualizes how SARS-CoV-2 has mutated from the original Wuhan strain into the Delta, Gamma, and Omicron variants.</p>
</div>
""", unsafe_allow_html=True)

# Sidebar configuration
with st.sidebar:
    st.title("Mutation Evolution")

    if st.button("🧬 See Variant Classifier"):
        st.switch_page("./app.py")

    st.markdown("---")
    st.markdown("### 🔧 Configuration")
    st.markdown("---")
    
    # Reference genome selection
    st.markdown("#### Reference Genomes")
    st.info("📁 Specify paths to your reference FASTA files")
    
    ref_paths = {}
    variant_colors = {
        "Wuhan": "#1f77b4",
        "Delta": "#ff7f0e", 
        "Gamma": "#2ca02c",
        "Omicron": "#d62728"
    }
    
    for name in ["Wuhan", "Delta", "Gamma", "Omicron"]:
        with st.expander(f"{name} Reference", expanded=False):
            ref_paths[name] = st.text_input(
                f"Path to {name} FASTA:",
                value=f"{name.lower()}.fasta",
                key=f"ref_{name.lower()}",
                help=f"Full path to the {name} reference genome FASTA file"
            )
    
    st.markdown("---")


with st.expander("🧬 What is a Mutation? Learn more"):
    st.markdown("""
    A **mutation** is a change in the nucleotide sequence of the genome. In viruses like SARS-CoV-2,
    mutations can affect how the virus spreads and responds to immunity.

    ### Common Types of Mutations:

    - **Substitution**  
      Replacement of one nucleotide by another (e.g., A → G). This can change amino acids and protein function.
      
    - **Insertion/Addition**  
      Addition of one or more nucleotides into the genome, potentially changing the reading frame.
      
    - **Deletion**  
      Removal of one or more nucleotides, which can disrupt gene coding or protein structure.

    ### Why Mutations Matter
    These changes lead to new variants, which may have different transmissibility or immune evasion abilities.
    """)
    st.image("mutation_types_illustration.png", caption="Types of mutations in genomes")

# --- LOAD & PREPARE SEQUENCES ---
seqs = {}
s_genes = {}
for var, path in variant_paths.items():
    seq = load_sequence(path)
    if seq is None:
        st.error(f"Failed to load sequence for {var} from {path}")
        st.stop()
    seqs[var] = seq
    s_genes[var] = extract_s_gene(seq, S_START, S_END)

# --- FIND MUTATIONS PER VARIANT COMPARED TO WUHAN ---
wuhan_s = s_genes["Wuhan"]

mutation_records = []
for var in ["Delta", "Gamma", "Omicron"]:
    aligned = align_sequences(wuhan_s, s_genes[var])
    muts = find_mutations(aligned[0], aligned[1], wuhan_s)
    for m in muts:
        m["variant"] = var
    mutation_records.extend(muts)

df_mut = pd.DataFrame(mutation_records)

if df_mut.empty:
    st.warning("No mutation data detected.")
    st.stop()

# --- 1. TABEL MUTASI TERPISAH PER VARIAN ---
st.markdown("### 🧬 Detailed Mutation Tables per Variant")

st.markdown("""##### 📋 Explanation
The tables below list all mutations detected in each variant (Delta, Gamma, Omicron) compared to the Wuhan reference sequence. 
These include mutation positions, types such as substitutions, deletions, or insertions, and details of nucleotide or protein changes.
""")

for var in ["Delta", "Gamma", "Omicron"]:
    st.markdown(f"#### Mutations in {var}")
    df_var = df_mut[df_mut["variant"] == var].reset_index(drop=True)
    if df_var.empty:
        st.write(f"No mutation data found for {var}.")
    else:
        st.dataframe(df_var)
st.markdown("---")

# --- 2. STEPWISE MUTATION EVOLUTION TIMELINE ---
st.markdown("### 📈 Stepwise Mutation Evolution Timeline")

st.markdown("""##### 🕒 Explanation 
This timeline illustrates the sequential accumulation of mutations from the original Wuhan strain to subsequent variants (Delta, Gamma, Omicron).
Each point marks a mutation that appeared in the respective variant, highlighting the chronological development of variant mutations.
""")

mutation_positions = {
    var: sorted(df_mut[df_mut["variant"] == var]["position_ref"].unique())
    for var in ["Delta", "Gamma", "Omicron"]
}

fig_stepwise = go.Figure()
order_map = {"Delta": 3, "Gamma": 2, "Omicron": 1}
for var in ["Delta", "Gamma", "Omicron"]:
    y_val = order_map[var]
    fig_stepwise.add_trace(go.Scatter(
        x=mutation_positions[var],
        y=[y_val] * len(mutation_positions[var]),
        mode="markers+lines",
        name=var,
        line=dict(shape="hv", color=variant_colors[var]),
        marker=dict(size=9),
        text=[f"Pos {pos}" for pos in mutation_positions[var]],
        hoverinfo="text+name"
    ))

fig_stepwise.update_layout(
    yaxis=dict(
        tickvals=list(order_map.values()),
        ticktext=list(order_map.keys()),
        autorange="reversed",
        title="Variant"
    ),
    xaxis_title="Genome Position",
    title="Mutations Accumulated from Wuhan (Stepwise)",
    height=400,
    plot_bgcolor='white',
    hovermode="closest"
)
st.plotly_chart(fig_stepwise, use_container_width=True)
st.markdown("---")

# --- 3. DENDROGRAM: HIERARKI JARAK MUTASI VARIAN ---
st.markdown("### 🌳 Dendrogram: Variant Mutation Distance")

st.markdown("""##### 🌳 Explanation 
This dendrogram visualizes the evolutionary relationship among variants based on their mutation patterns.
Variants positioned closer in the dendrogram share more similar genetic mutation profiles.
""")

# Membuat matriks biner mutasi per posisi per varian
positions = sorted(df_mut["position_ref"].unique())
variant_list = ["Wuhan", "Delta", "Gamma", "Omicron"]
binary_matrix = []

for var in variant_list:
    muts = set(df_mut[df_mut["variant"] == var]["position_ref"]) if var != "Wuhan" else set()
    binary_matrix.append([1 if pos in muts else 0 for pos in positions])

binary_matrix = np.array(binary_matrix)

# Menghitung linkage (jarak) dengan metode ward
linked = linkage(binary_matrix, method='ward')

fig_dendro = go.Figure()

dendro = dendrogram(linked, labels=variant_list, orientation='left', no_plot=True)
for i, dcoord in enumerate(dendro['dcoord']):
    x = dcoord
    y = dendro['icoord'][i]
    fig_dendro.add_trace(go.Scatter(x=x, y=y, mode='lines', line=dict(color='black')))

fig_dendro.update_layout(
    # yaxis=dict(tickmode='array', tickvals=np.arange(len(dendro['ivl'])), ticktext=dendro['ivl'], autorange='reversed'),
    yaxis=dict(tickmode='array', tickvals=[5,15,25,35], ticktext=dendro['ivl'], autorange='reversed'),
    xaxis_title="Distance",
    title="Hierarchical Clustering of Variants by Mutation Patterns",
    height=400,
    plot_bgcolor='white'
)
st.plotly_chart(fig_dendro, use_container_width=True)