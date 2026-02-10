# Computational Methods Documentation
## SARS-CoV-2 Transcriptomics Analysis Pipeline

### 1. Data Acquisition
**Dataset:** GSE147507  
**Source:** NCBI Gene Expression Omnibus (GEO)  
**Files:**
- Raw read counts: GSE147507_RawReadCounts_Human.tsv
- Sample metadata: GSE147507-GPL18573_series_matrix.txt

**Experimental Design:**
- Cell types: A549-ACE2, Calu-3 (lung epithelial cells)
- Conditions: Mock (n=9), SARS-CoV-2 infected (n=8), Ruxolitinib-treated (n=3)
- Sequencing platform: Illumina NextSeq 500

### 2. Quality Control & Preprocessing

**Library Size Filtering:**
```python
threshold = 1,000,000 reads
samples_retained = samples[library_size >= threshold]
```

**Gene Filtering:**
```python
min_count = 10
min_samples = 3
genes_retained = genes[(counts >= min_count).sum(axis=1) >= min_samples]
```

**Result:** 365 genes × 20 samples

### 3. Normalization

**Method:** DESeq2-style median-of-ratios normalization

**Algorithm:**
1. Calculate geometric mean per gene across samples
2. For each sample, compute ratio of counts to geometric mean
3. Take median ratio per sample (size factor)
4. Divide raw counts by size factors

**Transformation:** log2(normalized_counts + 1)

### 4. Differential Expression Analysis

**Statistical Model:** Welch's t-test (unequal variances)

**Comparison:** SARS-CoV-2 infected vs. Mock controls

**Multiple Testing Correction:** Benjamini-Hochberg FDR

**Significance Criteria:**
- |log2 Fold Change| ≥ 1.5
- FDR-adjusted p-value < 0.05

**Results:**
- Total genes tested: 13803
- Significant DEGs: 365
- Upregulated: 331
- Downregulated: 34

### 5. Functional Enrichment Analysis

**Tool:** Enrichr API (Ma'ayan Lab, Icahn School of Medicine)

**Gene Set Libraries:**
- GO Biological Process 2023
- GO Molecular Function 2023
- GO Cellular Component 2023
- KEGG Pathways 2021 Human

**Statistical Test:** Hypergeometric distribution with FDR correction

**Results:**
- Biological Process: 205 terms
- KEGG Pathways: 67 terms

### 6. Protein Annotation

**Source:** UniProt REST API

**Query:** Gene symbol → Protein function, recommended name

**Coverage:** Top 40 DEGs annotated

### 7. Network Analysis

**Method:** Gene co-expression network construction

**Correlation Metric:** Pearson correlation coefficient

**Edge Criteria:**
- |correlation| ≥ 0.7
- p-value < 0.05

**Network Properties:**
- Nodes: 80 genes
- Edges: 1,886 interactions
- Density: 0.597

**Centrality Measures:**
- Degree centrality: Number of connections / (N-1)
- Betweenness centrality: Fraction of shortest paths through node

**Visualization:** NetworkX spring layout algorithm

### 8. Software & Dependencies

**Python Version:** 3.13.2

**Core Packages:**
- pandas 2.2.3 (data manipulation)
- numpy 2.2.1 (numerical computing)
- scipy 1.15.1 (statistical tests)
- statsmodels 0.14.4 (FDR correction)
- networkx 3.4.2 (network analysis)
- matplotlib 3.10.0 (visualization)
- seaborn 0.13.2 (statistical plotting)
- requests 2.32.3 (API calls)

**Bioinformatics Packages:**
- scanpy 1.10.4 (single-cell analysis framework)
- biopython 1.84 (sequence analysis)

### 9. Computational Environment

**Operating System:** Windows 11  
**Processor:** Intel Core i5  
**GPU:** NVIDIA RTX 3050 (not utilized for this analysis)  
**RAM:** Sufficient for in-memory processing  
**Storage:** Local SSD

### 10. Reproducibility

**Version Control:** Git  
**Repository Structure:** Organized scripts, data, results  
**Execution:** Sequential phase-gate workflow  
**Runtime:** ~10-15 minutes total (excluding API calls)

---

**Documentation Generated:** 2026-02-10 07:44:10
