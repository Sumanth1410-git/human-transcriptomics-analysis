"""
Differential Expression Analysis (DEG)
SARS-CoV-2 Infected vs Mock Controls

Statistical Method:
- Welch's t-test (log2-transformed normalized counts)
- FDR correction (Benjamini-Hochberg)
- Fold change calculation

Significance Criteria:
- |log2FC| >= 1.5
- FDR < 0.05
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings('ignore')

# Set plotting style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (14, 8)

print("="*70)
print("PHASE 5: Differential Expression Analysis")
print("="*70 + "\n")

# ============================================================================
# 1. LOAD PROCESSED DATA
# ============================================================================
print("[1/6] Loading preprocessed data...")

# Load log2-transformed normalized counts
counts_log2 = pd.read_csv('data/processed/counts_log2_transformed.csv', index_col=0)
metadata = pd.read_csv('data/processed/metadata_covid_vs_mock.csv')

print(f"✓ Loaded: {counts_log2.shape[0]:,} genes × {counts_log2.shape[1]} samples")
print(f"Sample distribution:")
print(metadata['Condition'].value_counts())
print()

# ============================================================================
# 2. DEFINE COMPARISON GROUPS
# ============================================================================
print("[2/6] Defining comparison groups...")

# Primary comparison: SARS-CoV-2 infected vs Mock
# Exclude Ruxolitinib-treated samples for initial analysis
mock_samples = metadata[metadata['Condition'] == 'Mock']['Sample_ID'].tolist()
infected_samples = metadata[metadata['Condition'] == 'Infected']['Sample_ID'].tolist()

print(f"Mock controls: {len(mock_samples)} samples")
print(f"SARS-CoV-2 infected: {len(infected_samples)} samples")
print()

# Subset data
mock_expr = counts_log2[mock_samples]
infected_expr = counts_log2[infected_samples]

# ============================================================================
# 3. CALCULATE FOLD CHANGES & STATISTICS
# ============================================================================
print("[3/6] Computing differential expression statistics...")

deg_results = []

for gene in counts_log2.index:
    # Get expression values for each group
    mock_vals = mock_expr.loc[gene].values
    infected_vals = infected_expr.loc[gene].values
    
    # Mean expression (log2 scale)
    mean_mock = np.mean(mock_vals)
    mean_infected = np.mean(infected_vals)
    
    # Log2 Fold Change (Infected / Mock)
    log2fc = mean_infected - mean_mock  # Already in log2 scale
    
    # Welch's t-test (unequal variances)
    statistic, pval = stats.ttest_ind(infected_vals, mock_vals, equal_var=False)
    
    # Base mean (for filtering)
    base_mean = np.mean([mean_mock, mean_infected])
    
    deg_results.append({
        'Gene': gene,
        'BaseMean': base_mean,
        'Log2FoldChange': log2fc,
        'Mean_Mock': mean_mock,
        'Mean_Infected': mean_infected,
        'Statistic': statistic,
        'PValue': pval
    })

# Create DataFrame
deg_df = pd.DataFrame(deg_results)

# FDR correction (Benjamini-Hochberg)
_, padj, _, _ = multipletests(deg_df['PValue'], method='fdr_bh')
deg_df['FDR'] = padj

print(f"✓ Differential expression computed for {len(deg_df):,} genes")
print(f"Raw p-values range: [{deg_df['PValue'].min():.2e}, {deg_df['PValue'].max():.2e}]")
print(f"FDR-adjusted p-values range: [{deg_df['FDR'].min():.2e}, {deg_df['FDR'].max():.2e}]")
print()

# ============================================================================
# 4. CLASSIFY DIFFERENTIALLY EXPRESSED GENES
# ============================================================================
print("[4/6] Identifying significant DEGs...")

# Significance thresholds
log2fc_threshold = 1.5
fdr_threshold = 0.05

# Classify genes
deg_df['Significance'] = 'Not Significant'
deg_df.loc[(deg_df['Log2FoldChange'] >= log2fc_threshold) & 
           (deg_df['FDR'] < fdr_threshold), 'Significance'] = 'Upregulated'
deg_df.loc[(deg_df['Log2FoldChange'] <= -log2fc_threshold) & 
           (deg_df['FDR'] < fdr_threshold), 'Significance'] = 'Downregulated'

# Summary statistics
n_up = (deg_df['Significance'] == 'Upregulated').sum()
n_down = (deg_df['Significance'] == 'Downregulated').sum()
n_total_sig = n_up + n_down

print(f"Significance criteria:")
print(f"  - |log2FC| >= {log2fc_threshold}")
print(f"  - FDR < {fdr_threshold}")
print(f"\nDEG Summary:")
print(f"  Upregulated: {n_up:,} genes ({100*n_up/len(deg_df):.2f}%)")
print(f"  Downregulated: {n_down:,} genes ({100*n_down/len(deg_df):.2f}%)")
print(f"  Total significant: {n_total_sig:,} genes ({100*n_total_sig/len(deg_df):.2f}%)")
print(f"  Not significant: {(deg_df['Significance'] == 'Not Significant').sum():,} genes")
print()

# ============================================================================
# 5. EXTRACT TOP DEGS
# ============================================================================
print("[5/6] Extracting top differentially expressed genes...")

# Sort by absolute fold change
deg_df_sorted = deg_df.sort_values('Log2FoldChange', ascending=False)

# Top upregulated
top_up = deg_df_sorted[deg_df_sorted['Significance'] == 'Upregulated'].head(20)
print("\nTop 20 Upregulated Genes (SARS-CoV-2 Infection):")
print(top_up[['Gene', 'Log2FoldChange', 'FDR', 'Mean_Mock', 'Mean_Infected']].to_string(index=False))

# Top downregulated
top_down = deg_df_sorted[deg_df_sorted['Significance'] == 'Downregulated'].tail(20).iloc[::-1]
print("\nTop 20 Downregulated Genes (SARS-CoV-2 Infection):")
print(top_down[['Gene', 'Log2FoldChange', 'FDR', 'Mean_Mock', 'Mean_Infected']].to_string(index=False))
print()

# ============================================================================
# 6. VISUALIZATION: VOLCANO PLOT
# ============================================================================
print("[6/6] Generating volcano plot...")

fig, ax = plt.subplots(1, 1, figsize=(12, 8))

# Plot all genes
colors = {'Upregulated': 'red', 'Downregulated': 'blue', 'Not Significant': 'gray'}
for sig_type in ['Not Significant', 'Downregulated', 'Upregulated']:
    subset = deg_df[deg_df['Significance'] == sig_type]
    ax.scatter(subset['Log2FoldChange'], 
               -np.log10(subset['FDR']),
               c=colors[sig_type],
               alpha=0.6 if sig_type == 'Not Significant' else 0.8,
               s=20 if sig_type == 'Not Significant' else 40,
               label=f"{sig_type} (n={len(subset):,})",
               edgecolors='none')

# Add threshold lines
ax.axhline(-np.log10(fdr_threshold), color='black', linestyle='--', 
           linewidth=1.5, label=f'FDR = {fdr_threshold}')
ax.axvline(log2fc_threshold, color='black', linestyle='--', linewidth=1.5)
ax.axvline(-log2fc_threshold, color='black', linestyle='--', linewidth=1.5)

# Annotate top genes
top_genes_to_label = pd.concat([
    deg_df_sorted[deg_df_sorted['Significance'] == 'Upregulated'].head(5),
    deg_df_sorted[deg_df_sorted['Significance'] == 'Downregulated'].tail(5)
])

for _, row in top_genes_to_label.iterrows():
    ax.annotate(row['Gene'], 
                xy=(row['Log2FoldChange'], -np.log10(row['FDR'])),
                xytext=(5, 5), textcoords='offset points',
                fontsize=9, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.7))

# Labels and formatting
ax.set_xlabel('Log2 Fold Change (Infected / Mock)', fontsize=14, fontweight='bold')
ax.set_ylabel('-Log10(FDR)', fontsize=14, fontweight='bold')
ax.set_title('Volcano Plot: SARS-CoV-2 Infected vs Mock\nDifferential Gene Expression',
             fontsize=16, fontweight='bold')
ax.legend(loc='upper right', frameon=True, shadow=True)
ax.grid(True, alpha=0.3)
ax.set_xlim([-12, 12])

plt.tight_layout()
plt.savefig('results/figures/06_volcano_plot.png', dpi=300, bbox_inches='tight')
print("✓ Figure saved: results/figures/06_volcano_plot.png")

# ============================================================================
# 7. SAVE DEG RESULTS
# ============================================================================
print("\n[7/6] Saving DEG results...")

# Save full results
deg_df_sorted.to_csv('results/tables/deg_full_results.csv', index=False)

# Save significant genes only
deg_sig = deg_df_sorted[deg_df_sorted['Significance'] != 'Not Significant']
deg_sig.to_csv('results/tables/deg_significant_only.csv', index=False)

# Save top genes for annotation
top_deg = pd.concat([top_up, top_down])
top_deg.to_csv('results/tables/deg_top40_genes.csv', index=False)

print("✓ Saved results:")
print("  1. results/tables/deg_full_results.csv (all genes)")
print("  2. results/tables/deg_significant_only.csv (significant DEGs only)")
print("  3. results/tables/deg_top40_genes.csv (top 40 genes for annotation)")

# ============================================================================
# SUMMARY REPORT
# ============================================================================
print("\n" + "="*70)
print("PHASE 5 DEG ANALYSIS SUMMARY")
print("="*70)
print(f"Total genes tested: {len(deg_df):,}")
print(f"Significantly upregulated: {n_up:,} ({100*n_up/len(deg_df):.2f}%)")
print(f"Significantly downregulated: {n_down:,} ({100*n_down/len(deg_df):.2f}%)")
print(f"Total DEGs: {n_total_sig:,} ({100*n_total_sig/len(deg_df):.2f}%)")
print(f"\nTop upregulated gene: {top_up.iloc[0]['Gene']} (log2FC = {top_up.iloc[0]['Log2FoldChange']:.2f})")
print(f"Top downregulated gene: {top_down.iloc[0]['Gene']} (log2FC = {top_down.iloc[0]['Log2FoldChange']:.2f})")
print(f"\nNext Step: Phase 6 - Functional Annotation (GO Enrichment)")
print("="*70)

# Save summary with UTF-8 encoding (fixes Windows encoding error)
with open('results/tables/phase5_deg_summary.txt', 'w', encoding='utf-8') as f:
    f.write("PHASE 5 DEG ANALYSIS SUMMARY\n")
    f.write("="*50 + "\n")
    f.write(f"Total genes tested: {len(deg_df):,}\n")
    f.write(f"Upregulated: {n_up:,}\n")
    f.write(f"Downregulated: {n_down:,}\n")
    f.write(f"Total DEGs: {n_total_sig:,}\n")
    f.write(f"Thresholds: |log2FC| >= {log2fc_threshold}, FDR < {fdr_threshold}\n")

print("\n✓ Summary saved: results/tables/phase5_deg_summary.txt")