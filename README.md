<h1 align="center">🧪 SARS-CoV-2 Variant Classifier 🧪</h1>
<div align="center"> A web-based application to identify SARS-CoV-2 variants (Delta, Gamma, and Omicron) by analyzing DNA/RNA sequences of the spike protein (gen S). </div>

## Table of Contents 📚
* [General Information](#general-information-🧬)
* [How to Run](#how-to-run-ℹ️)
* [Contributors](#contributors-👩‍💻)

## General Information 🧬
This project provides a user-friendly platform to identify SARS-CoV-2 variants through genomic sequence analysis. Using global alignment algorithms, the system detects mutations from user-submitted FASTA files and classifies the variant based on predefined mutational profiles.

**Key Features:**
* Alignment with the reference Wuhan strain using Needleman-Wunsch algorithm.
* Detection of key mutations (substitution, deletion, insertion).
* Real-time classification and mutation visualization.
* No technical expertise required — just upload your sequence!

**File Requirements & Information**

**Accepted formats:** `.fasta`, `.fa`  
**Input types:**
- Full SARS-CoV-2 genome sequence
- Spike (S) gene sequence only  

**S-gene coordinates:** Position **21,563 – 25,384** (1-based)

**Analysis process:**
1. Extract Spike gene from full genome (if applicable)  
2. Align with Wuhan reference  
3. Identify mutations  
4. Compare with variant signatures  
5. Classify based on mutation matches


## How to Run ℹ️
The application is live and can be accessed here:
```
https://covid-classifier-kds.streamlit.app/
```
1. Prepare your DNA/RNA sequence file in `.fasta` format,  
   **or use one of the example files available in this repository**:
   - `delta1.fasta`, `delta2.fasta`
   - `gamma1.fasta`
   - `omicron1.fasta`, `omicron2.fasta`  
   > These sequences are taken from **NCBI GenBank**.

2. Upload the file using the provided form.

3. View classification results and mutation insights instantly.

## Contributors 👩‍💻
### Kelompok 1 - K05

| NIM      | Nama                      |
|----------|---------------------------|
| 18222058 | Matthew Nicholas Gunawan     |
| 18222074 | Kayla Dyara               |
| 18222078 | Monica Angela Hartono     |