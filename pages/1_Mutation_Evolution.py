import streamlit as st
import pandas as pd
import plotly.express as px
from utils import load_sequence, extract_s_gene, align_sequences, find_mutations

st.title("🧬 Evolution from Wuhan to Variants")
st.markdown("""
This section visualizes how SARS-CoV-2 has mutated from the original Wuhan strain into the Delta, Gamma, and Omicron variants.
""")

# Load Wuhan and variant sequences
variant_paths = {
    "Wuhan": "wuhan.fasta",
    "Delta": "delta.fasta",
    "Gamma": "gamma.fasta",
    "Omicron": "omicron.fasta"
}
S_START, S_END = 21563, 25384

ref_seq = load_sequence(variant_paths["Wuhan"])
ref_s = extract_s_gene(ref_seq, start=S_START, end=S_END)

mutation_summary = []

for var, path in variant_paths.items():
    if var == "Wuhan":
        continue
    var_seq = load_sequence(path)
    var_s = extract_s_gene(var_seq, start=S_START, end=S_END)
    aligned = align_sequences(ref_s, var_s)
    muts = find_mutations(aligned[0], aligned[1], ref_s)
    for m in muts:
        m["variant"] = var
    mutation_summary.extend(muts)

df_mut = pd.DataFrame(mutation_summary)

if not df_mut.empty:
    st.markdown("### 🔍 Mutation Distribution by Variant")
    fig = px.scatter(
        df_mut[df_mut["type"] == "Substitution"],
        x="position_ref", y="variant",
        color="variant",
        title="Substitution Positions Across Variants",
        labels={"position_ref": "Reference Position"},
        height=500
    )
    st.plotly_chart(fig, use_container_width=True)

    st.markdown("### 🧬 Mutation Table")
    st.dataframe(df_mut)
else:
    st.warning("No mutation data found.")
