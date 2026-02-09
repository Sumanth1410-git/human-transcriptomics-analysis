"""
Functional Annotation & Gene Ontology Enrichment
Phase 6: Map DEGs to Biological Processes

Methods:
1. GO Enrichment Analysis (Biological Process, Molecular Function, Cellular Component)
2. Over-representation test (hypergeometric distribution)
3. FDR correction for multiple testing
4. UniProt protein function annotation

Input: deg_significant_only.csv (365 DEGs)
Output: Enriched GO terms, annotated gene tables
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests
import requests
import time
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

# Set plotting style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (14, 8)

print("="*70)
print("PHASE 6: Functional Annotation & GO Enrichment")
print("="*70 + "\n")

# ============================================================================
# 1. LOAD SIGNIFICANT DEGs
# ============================================================================
print("[1/7] Loading significant DEGs...")

deg_sig = pd.read_csv('results/tables/deg_significant_only.csv')
all_genes = pd.read_csv('results/tables/deg_full_results.csv')

print(f"✓ Significant DEGs: {len(deg_sig)} genes")
print(f"  - Upregulated: {(deg_sig['Significance'] == 'Upregulated').sum()}")
print(f"  - Downregulated: {(deg_sig['Significance'] == 'Downregulated').sum()}")
print(f"✓ Background genes: {len(all_genes):,} genes")
print()

# Extract gene lists
deg_genes = set(deg_sig['Gene'].tolist())
background_genes = set(all_genes['Gene'].tolist())

# ============================================================================
# 2. GENE ONTOLOGY ENRICHMENT via Enrichr API
# ============================================================================
print("[2/7] Performing GO enrichment analysis via Enrichr...")
print("(This may take 1-2 minutes...)\n")

def enrichr_analysis(gene_list, gene_set_library='GO_Biological_Process_2023'):
    """
    Perform enrichment analysis using Enrichr API
    
    Parameters:
    -----------
    gene_list : list
        List of gene symbols
    gene_set_library : str
        Enrichr library name
        
    Returns:
    --------
    pd.DataFrame : Enrichment results
    """
    # Submit gene list
    ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/addList'
    genes_str = '\n'.join(gene_list)
    payload = {
        'list': (None, genes_str),
    }
    
    try:
        response = requests.post(ENRICHR_URL, files=payload)
        if not response.ok:
            print(f"  ✗ Error submitting genes to Enrichr")
            return pd.DataFrame()
        
        data = response.json()
        user_list_id = data['userListId']
        
        # Retrieve enrichment results
        ENRICHR_ENRICH_URL = f'https://maayanlab.cloud/Enrichr/enrich'
        query_string = f'?userListId={user_list_id}&backgroundType={gene_set_library}'
        
        response = requests.get(ENRICHR_ENRICH_URL + query_string)
        if not response.ok:
            print(f"  ✗ Error retrieving enrichment results")
            return pd.DataFrame()
        
        data = response.json()
        
        # Parse results
        if gene_set_library in data:
            results = []
            for entry in data[gene_set_library]:
                results.append({
                    'Term': entry[1],
                    'P_value': entry[2],
                    'Adjusted_P_value': entry[6],
                    'Odds_Ratio': entry[3],
                    'Combined_Score': entry[4],
                    'Genes': ';'.join(entry[5])
                })
            return pd.DataFrame(results)
        else:
            return pd.DataFrame()
            
    except Exception as e:
        print(f"  ✗ Enrichr API error: {e}")
        return pd.DataFrame()

# Perform enrichment for different GO categories
go_categories = {
    'Biological_Process': 'GO_Biological_Process_2023',
    'Molecular_Function': 'GO_Molecular_Function_2023',
    'Cellular_Component': 'GO_Cellular_Component_2023',
    'KEGG_Pathways': 'KEGG_2021_Human'
}

enrichment_results = {}

for category_name, library in go_categories.items():
    print(f"  Analyzing: {category_name}...")
    time.sleep(1)  # Rate limiting
    results = enrichr_analysis(list(deg_genes), library)
    
    if not results.empty:
        # Filter significant terms (Adjusted P < 0.05)
        results_sig = results[results['Adjusted_P_value'] < 0.05].copy()
        enrichment_results[category_name] = results_sig
        print(f"    ✓ {len(results_sig)} significant terms found")
    else:
        enrichment_results[category_name] = pd.DataFrame()
        print(f"    ⚠ No results returned")

print()

# ============================================================================
# 3. SAVE ENRICHMENT RESULTS
# ============================================================================
print("[3/7] Saving GO enrichment results...")

for category, results in enrichment_results.items():
    if not results.empty:
        filename = f"results/tables/go_enrichment_{category}.csv"
        results.to_csv(filename, index=False)
        print(f"  ✓ {filename} ({len(results)} terms)")

print()

# ============================================================================
# 4. VISUALIZE TOP ENRICHED TERMS
# ============================================================================
print("[4/7] Visualizing top enriched GO terms...")

for category, results in enrichment_results.items():
    if results.empty or len(results) == 0:
        print(f"  ⚠ Skipping {category} (no significant terms)")
        continue
    
    # Take top 15 terms by combined score
    top_terms = results.nlargest(15, 'Combined_Score')
    
    # Create horizontal bar plot
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    
    # Clean term names (remove GO ID)
    top_terms['Term_Clean'] = top_terms['Term'].str.split('\\(GO:').str[0].str.strip()
    
    # Plot
    bars = ax.barh(range(len(top_terms)), top_terms['Combined_Score'].values, 
                   color='steelblue', edgecolor='black')
    
    # Color by significance
    colors = plt.cm.RdYlGn_r(np.log10(top_terms['Adjusted_P_value'].values + 1e-10))
    for bar, color in zip(bars, colors):
        bar.set_color(color)
    
    ax.set_yticks(range(len(top_terms)))
    ax.set_yticklabels(top_terms['Term_Clean'].values, fontsize=10)
    ax.set_xlabel('Combined Score', fontsize=12, fontweight='bold')
    ax.set_title(f'Top 15 Enriched Terms: {category.replace("_", " ")}', 
                 fontsize=14, fontweight='bold')
    ax.invert_yaxis()
    ax.grid(axis='x', alpha=0.3)
    
    # Add colorbar for p-values
    sm = plt.cm.ScalarMappable(cmap='RdYlGn_r', 
                                norm=plt.Normalize(vmin=np.log10(top_terms['Adjusted_P_value'].min()),
                                                  vmax=np.log10(top_terms['Adjusted_P_value'].max())))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label('log10(Adjusted P-value)', fontsize=10)
    
    plt.tight_layout()
    filename = f"results/figures/07_go_enrichment_{category}.png"
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    print(f"  ✓ {filename}")
    plt.close()

print()

# ============================================================================
# 5. UNIPROT ANNOTATION (Top 40 DEGs)
# ============================================================================
print("[5/7] Retrieving UniProt protein annotations...")
print("(Fetching data for top 40 DEGs...)\n")

top_degs = pd.read_csv('results/tables/deg_top40_genes.csv')

def get_uniprot_annotation(gene_symbol):
    """
    Retrieve protein function from UniProt API
    
    Parameters:
    -----------
    gene_symbol : str
        Human gene symbol
        
    Returns:
    --------
    dict : Protein annotation
    """
    try:
        # Query UniProt for human gene
        url = f"https://rest.uniprot.org/uniprotkb/search?query=gene:{gene_symbol}+AND+organism_id:9606&fields=accession,gene_names,protein_name,cc_function&format=json"
        
        response = requests.get(url, timeout=10)
        
        if response.ok:
            data = response.json()
            if 'results' in data and len(data['results']) > 0:
                entry = data['results'][0]
                
                # Extract function (first entry if multiple)
                function = "N/A"
                if 'comments' in entry:
                    for comment in entry['comments']:
                        if comment['commentType'] == 'FUNCTION':
                            function = comment['texts'][0]['value'][:200] + "..."
                            break
                
                return {
                    'UniProt_ID': entry.get('primaryAccession', 'N/A'),
                    'Protein_Name': entry.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', 'N/A'),
                    'Function': function
                }
        
        return {'UniProt_ID': 'N/A', 'Protein_Name': 'N/A', 'Function': 'N/A'}
        
    except Exception as e:
        return {'UniProt_ID': 'Error', 'Protein_Name': 'Error', 'Function': str(e)}

# Annotate top DEGs
annotations = []
for idx, row in top_degs.iterrows():
    gene = row['Gene']
    print(f"  Fetching: {gene}...", end=' ')
    
    annotation = get_uniprot_annotation(gene)
    annotations.append({
        'Gene': gene,
        'Log2FoldChange': row['Log2FoldChange'],
        'FDR': row['FDR'],
        'Significance': row['Significance'],
        **annotation
    })
    
    print("✓")
    time.sleep(0.5)  # Rate limiting

annotations_df = pd.DataFrame(annotations)
annotations_df.to_csv('results/tables/deg_top40_annotated.csv', index=False)

print(f"\n✓ Annotations saved: results/tables/deg_top40_annotated.csv")
print()

# Display top upregulated annotated genes
print("Top 10 Upregulated Genes with Annotations:")
print(annotations_df[annotations_df['Significance'] == 'Upregulated']
      .head(10)[['Gene', 'Log2FoldChange', 'Protein_Name', 'Function']]
      .to_string(index=False))
print()

# ============================================================================
# 6. CREATE SUMMARY REPORT
# ============================================================================
print("[6/7] Generating comprehensive summary report...")

# Count enriched terms per category
summary_stats = []
for category, results in enrichment_results.items():
    summary_stats.append({
        'Category': category,
        'Significant_Terms': len(results),
        'Top_Term': results.iloc[0]['Term'] if not results.empty else 'N/A',
        'Top_P_value': results.iloc[0]['Adjusted_P_value'] if not results.empty else np.nan
    })

summary_df = pd.DataFrame(summary_stats)
summary_df.to_csv('results/tables/enrichment_summary.csv', index=False)

print("✓ Summary saved: results/tables/enrichment_summary.csv")
print()

# ============================================================================
# 7. DISPLAY KEY BIOLOGICAL PROCESSES
# ============================================================================
print("[7/7] Extracting key biological insights...")

if 'Biological_Process' in enrichment_results and not enrichment_results['Biological_Process'].empty:
    bp_results = enrichment_results['Biological_Process']
    
    print("\nTop 10 Enriched Biological Processes (SARS-CoV-2 Response):")
    print("="*70)
    for idx, row in bp_results.head(10).iterrows():
        term_clean = row['Term'].split('(GO:')[0].strip()
        print(f"{idx+1}. {term_clean}")
        print(f"   P-value: {row['Adjusted_P_value']:.2e} | Genes: {row['Genes'][:80]}...")
        print()

# ============================================================================
# FINAL SUMMARY
# ============================================================================
print("="*70)
print("PHASE 6 FUNCTIONAL ANNOTATION SUMMARY")
print("="*70)
print(f"Total DEGs analyzed: {len(deg_sig)}")
print(f"\nEnrichment Results:")
for category, results in enrichment_results.items():
    print(f"  {category}: {len(results)} significant terms")
print(f"\nTop annotated DEGs: {len(annotations_df)}")
print(f"\nKey finding: SARS-CoV-2 infection activates:")
if 'Biological_Process' in enrichment_results and not enrichment_results['Biological_Process'].empty:
    top_process = enrichment_results['Biological_Process'].iloc[0]['Term'].split('(GO:')[0].strip()
    print(f"  - {top_process}")
print(f"\nNext Step: Phase 7 - Network Analysis & Visualization")
print("="*70)

# Save final summary
with open('results/tables/phase6_annotation_summary.txt', 'w', encoding='utf-8') as f:
    f.write("PHASE 6 FUNCTIONAL ANNOTATION SUMMARY\n")
    f.write("="*50 + "\n")
    f.write(f"DEGs analyzed: {len(deg_sig)}\n")
    for category, results in enrichment_results.items():
        f.write(f"{category}: {len(results)} terms\n")
    f.write(f"Top DEGs annotated: {len(annotations_df)}\n")

print("\n✓ Summary saved: results/tables/phase6_annotation_summary.txt")