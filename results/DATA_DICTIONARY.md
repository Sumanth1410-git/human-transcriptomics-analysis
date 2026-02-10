# Data Dictionary

## Key Result Files

### Differential Expression
| File | Description | Columns | Rows |
|------|-------------|---------|------|
| `deg_full_results.csv` | All tested genes | Gene, BaseMean, Log2FoldChange, PValue, FDR, Significance | 13,803 |
| `deg_significant_only.csv` | Significant DEGs only | Same as above | 365 |
| `deg_top40_genes.csv` | Top 40 DEGs by |log2FC| | Same as above | 40 |
| `deg_top40_annotated.csv` | Top 40 with UniProt annotations | +UniProt_ID, Protein_Name, Function | 40 |

### Functional Enrichment
| File | Description | Columns |
|------|-------------|---------|
| `go_enrichment_Biological_Process.csv` | Enriched BP terms | Term, P_value, Adjusted_P_value, Genes |
| `go_enrichment_Molecular_Function.csv` | Enriched MF terms | Same |
| `go_enrichment_Cellular_Component.csv` | Enriched CC terms | Same |
| `go_enrichment_KEGG_Pathways.csv` | Enriched pathways | Same |

### Network Analysis
| File | Description | Columns | Rows |
|------|-------------|---------|------|
| `network_hub_genes.csv` | Network centrality metrics | Gene, Degree_Centrality, Betweenness_Centrality, Log2FC, FDR | 80 |

### Processed Data
| File | Description |
|------|-------------|
| `counts_filtered_raw.csv` | QC-filtered raw counts |
| `counts_normalized.csv` | DESeq2-normalized counts |
| `counts_log2_transformed.csv` | Log2-transformed normalized counts |
| `metadata_covid_vs_mock.csv` | Sample annotations |

## Column Definitions

**Gene:** HGNC gene symbol  
**Log2FoldChange:** log2(Infected / Mock), positive = upregulated in infection  
**FDR:** False Discovery Rate (Benjamini-Hochberg adjusted p-value)  
**Significance:** Classification (Upregulated, Downregulated, Not Significant)  
**BaseMean:** Mean normalized expression across all samples  
**Degree_Centrality:** (# connections) / (N-1), range [0, 1]  
**Betweenness_Centrality:** Fraction of shortest paths through gene  
**Combined_Score:** Enrichr metric = log(p-value) × z-score  

---

**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
