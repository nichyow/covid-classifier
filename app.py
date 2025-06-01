import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from utils import (
    load_sequence,
    extract_s_gene,
    align_sequences,
    find_mutations,
    count_matching_mutations,
    classify_variant
)

# Page configuration
st.set_page_config(
    page_title="SARS-CoV-2 Variant Classifier",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for better styling
st.markdown("""
<style>
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
    
    .info-card {
        background: #f8f9fa;
        padding: 1.5rem;
        border-radius: 10px;
        border-left: 4px solid #667eea;
        margin: 1rem 0;
    }
    
    .success-card {
        background: #d4edda;
        padding: 1.5rem;
        border-radius: 10px;
        border-left: 4px solid #28a745;
        margin: 1rem 0;
    }
    
    .warning-card {
        background: #fff3cd;
        padding: 1.5rem;
        border-radius: 10px;
        border-left: 4px solid #ffc107;
        margin: 1rem 0;
    }
    
    .error-card {
        background: #f8d7da;
        padding: 1.5rem;
        border-radius: 10px;
        border-left: 4px solid #dc3545;
        margin: 1rem 0;
    }
    
    .metric-container {
        background: white;
        padding: 1rem;
        border-radius: 8px;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        text-align: center;
        margin: 0.5rem 0;
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
    <h1>🧬 SARS-CoV-2 Variant Classifier</h1>
    <p>Advanced genomic analysis for COVID-19 variant identification</p>
</div>
""", unsafe_allow_html=True)

# Define S gene coordinates
S_START = 21563
S_END = 25384

# Sidebar configuration
with st.sidebar:
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
    
    # Analysis parameters
    # st.markdown("#### Analysis Parameters")
    # wuhan_threshold = st.slider(
    #     "Wuhan-like threshold",
    #     min_value=1,
    #     max_value=10,
    #     value=3,
    #     help="Minimum mutations required to classify as non-Wuhan variant"
    # )

# Main content area
col1, col2 = st.columns([2, 1])

with col1:
    st.markdown("### 📂 Upload Patient Sample")
    
    # File uploader with enhanced styling
    uploaded = st.file_uploader(
        "Choose your FASTA file",
        type=["fasta", "fa"],
        key="patient_upload",
        help="Upload either a full genome or Spike gene sequence in FASTA format"
    )
    
    # Information about file requirements
    with st.expander("ℹ️ File Requirements & Information"):
        st.markdown("""
        **Accepted formats:** `.fasta`, `.fa`
        
        **Input types:**
        - Full SARS-CoV-2 genome sequence
        - Spike (S) gene sequence only
        
        **S-gene coordinates:** Position {start:,} - {end:,} (1-based)
        
        **Analysis process:**
        1. Extract Spike gene from full genome (if applicable)
        2. Align with Wuhan reference
        3. Identify mutations
        4. Compare with variant signatures
        5. Classification based on mutation matches
        """.format(start=S_START, end=S_END))

with col2:
    st.markdown("### 📊 System Status")
    
    # Load and process reference sequences
    ref_s = {}
    ref_status = {}
    
    for name, path in ref_paths.items():
        try:
            seq = load_sequence(path)
            if seq:
                sseq = extract_s_gene(seq, start=S_START, end=S_END)
                ref_s[name] = sseq
                ref_status[name] = "✅ Loaded"
            else:
                ref_status[name] = "❌ Failed"
        except Exception as e:
            ref_status[name] = "❌ Error"
    
    # Display reference status
    for name, status in ref_status.items():
        color = variant_colors.get(name, "#666666")
        if "✅" in status:
            st.success(f"{name}: {status}")
        else:
            st.error(f"{name}: {status}")

# Build variant profiles
variant_profiles = {}
if "Wuhan" in ref_s:
    wuhan_s = ref_s["Wuhan"]
    for var, sseq in ref_s.items():
        if var == "Wuhan":
            continue
        try:
            aln = align_sequences(wuhan_s, sseq)
            if aln:
                aligned_w, aligned_v, *_ = aln
                variant_profiles[var] = find_mutations(aligned_w, aligned_v, wuhan_s)
        except Exception as e:
            st.warning(f"Failed to process {var}: {e}")

# Analysis section
if uploaded:
    st.markdown("---")
    st.markdown("### 🔬 Analysis Results")
    
    try:
        # Load patient sequence
        seq = load_sequence(uploaded)
        if seq is None:
            raise ValueError("Failed to load sequence from uploaded file")
        
        # Extract S-gene or use sequence as-is
        if len(seq) >= S_END:
            patient_s = extract_s_gene(seq, start=S_START, end=S_END)
        else:
            patient_s = seq
        # Perform alignment
        if "Wuhan" not in ref_s:
            st.error("❌ Wuhan reference not available for comparison")
        else:
            aln = align_sequences(ref_s["Wuhan"], patient_s)
            if not aln:
                st.error("❌ Alignment failed between patient and reference sequences")
            else:
                aw, ap, score, *_ = aln
                
                # Display alignment metrics
                col1, col2, col3 = st.columns(3)
                with col1:
                    st.metric("Alignment Score", f"{score:.1f}")
                with col2:
                    st.metric("Patient Sequence Length", len(patient_s))
                with col3:
                    st.metric("Reference Length", len(ref_s["Wuhan"]))
                
                # Find mutations
                patient_mutations = find_mutations(aw, ap, ref_s["Wuhan"])
                
                # Calculate matches
                match_counts = {}
                for var, profile in variant_profiles.items():
                    match_counts[var] = count_matching_mutations(patient_mutations, profile)
                
                # Visualizations
                if match_counts:
                    col1, col2 = st.columns([2, 1])
                    
                    with col1:
                        st.markdown("#### 📈 Variant Signature Matches")
                        
                        # Create enhanced bar chart with Plotly
                        df = pd.DataFrame.from_dict(match_counts, orient='index', columns=['Matches'])
                        df['Variant'] = df.index
                        df['Color'] = [variant_colors.get(var, "#666666") for var in df.index]
                        
                        fig = px.bar(
                            df, 
                            x='Variant', 
                            y='Matches',
                            color='Variant',
                            color_discrete_map=variant_colors,
                            title="Mutation Matches by Variant"
                        )
                        fig.update_layout(
                            showlegend=False,
                            plot_bgcolor='rgba(0,0,0,0)',
                            paper_bgcolor='rgba(0,0,0,0)',
                        )
                        st.plotly_chart(fig, use_container_width=True)
                    
                    with col2:
                        st.markdown("#### 🎯 Classification")
                        
                        # Perform classification
                        identified_var, match_score = classify_variant(
                            patient_mutations, variant_profiles
                        )
                        
                        # Display result with appropriate styling
                        if "Wuhan-like" in identified_var:
                            st.markdown(f"""
                            <div class="warning-card">
                                <h4>🟡 Classification Result</h4>
                                <p><strong>{identified_var}</strong></p>
                                <p>Mutations matched: {match_score}</p>
                            </div>
                            """, unsafe_allow_html=True)
                        else:
                            st.markdown(f"""
                            <div class="success-card">
                                <h4>🟢 Classification Result</h4>
                                <p><strong>{identified_var}</strong></p>
                                <p>Mutations matched: {match_score}</p>
                            </div>
                            """, unsafe_allow_html=True)
                
                # Mutation details
                if patient_mutations:
                    with st.expander("🧬 Detailed Mutation Analysis"):
                        st.markdown("#### Identified Mutations")
                        
                        mutation_df = pd.DataFrame(patient_mutations)
                        if not mutation_df.empty:
                            st.dataframe(mutation_df, use_container_width=True)
                        
                        # Mutation type breakdown
                        mutation_types = {}
                        for mut in patient_mutations:
                            mut_type = mut.get('type', 'Unknown')
                            mutation_types[mut_type] = mutation_types.get(mut_type, 0) + 1
                        
                        if mutation_types:
                            st.markdown("#### Mutation Type Distribution")
                            type_df = pd.DataFrame.from_dict(mutation_types, orient='index', columns=['Count'])
                            st.bar_chart(type_df)
    
    except Exception as e:
        st.markdown(f"""
        <div class="error-card">
            <h4>❌ Analysis Error</h4>
            <p>{str(e)}</p>
        </div>
        """, unsafe_allow_html=True)

# Footer
st.markdown("---")
st.markdown("""
<div style="text-align: center; color: #666; padding: 2rem;">
    <p>🧬 SARS-CoV-2 Variant Classifier | Powered by Bioinformatics & Streamlit</p>
</div>
""", unsafe_allow_html=True)