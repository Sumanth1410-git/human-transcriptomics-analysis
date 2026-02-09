"""
Data Preprocessing & Normalization Pipeline
GSE147507: SARS-CoV-2 Transcriptomics

Steps:
1. Filter low-expression genes
2. Subset SARS-CoV-2 vs Mock samples
3. Normalize counts (DESeq2-style median-of-ratios)
4. Log2-transform with pseudocount
5. Save processed data for DEG analysis
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings('ignore')

# Set plotting style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (14, 6)

print("="*70)
print("PHASE 4: Data Preprocessing & Normalization")
print("="*70 + "\n")

# ============================================================================
# 1. LOAD RAW DATA
# ============================================================================
print("[1/7] Loading raw count matrix...")

counts_df = pd.read_csv('data/raw/covid19_raw_counts.tsv', sep='\t', index_col=0)
print(f"✓ Original data: {counts_df.shape[0]} genes × {counts_df.shape[1]} samples\n")

# ============================================================================
# 2. PARSE SAMPLE METADATA FROM COLUMN NAMES
# ============================================================================
print("[2/7] Extracting sample metadata from column names...")

# Parse sample names: "Series1_NHBE_Mock_1" → Extract cell type, condition, replicate
metadata_list = []

for sample_name in counts_df.columns:
    parts = sample_name.split('_')
    
    # Extract components
    series = parts[0]
    cell_type = parts[1]
    
    # Determine condition (handle variable naming)
    if 'Mock' in sample_name:
        condition = 'Mock'
        virus = 'None'
    elif 'SARS-CoV-2' in sample_name:
        condition = 'Infected'
        virus = 'SARS-CoV-2'
        # Check for drug treatment (Ruxolitinib)
        if 'Rux' in sample_name:
            condition = 'Infected_Rux'
    elif 'RSV' in sample_name:
        condition = 'Infected'
        virus = 'RSV'
    elif 'IAV' in sample_name:
        condition = 'Infected'
        virus = 'IAV'
    elif 'HPIV3' in sample_name:
        condition = 'Infected'
        virus = 'HPIV3'
    elif 'IFNB' in sample_name:
        condition = 'IFNB_treated'
        virus = 'None'
    elif 'COVID19Lung' in sample_name:
        condition = 'COVID19_patient'
        virus = 'SARS-CoV-2'
        cell_type = 'LungBiopsy'
    elif 'HealthyLung' in sample_name:
        condition = 'Healthy_patient'
        virus = 'None'
        cell_type = 'LungBiopsy'
    else:
        condition = 'Unknown'
        virus = 'Unknown'
    
    metadata_list.append({
        'Sample_ID': sample_name,
        'Series': series,
        'Cell_Type': cell_type,
        'Condition': condition,
        'Virus': virus,
        'Library_Size': counts_df[sample_name].sum()
    })

metadata_df = pd.DataFrame(metadata_list)
print(metadata_df.groupby(['Condition', 'Virus']).size())
print(f"\n✓ Metadata extracted for {len(metadata_df)} samples")

# Save metadata
metadata_df.to_csv('data/processed/sample_metadata_annotated.csv', index=False)
print("✓ Saved: data/processed/sample_metadata_annotated.csv\n")

# ============================================================================
# 3. FILTER LOW-QUALITY SAMPLES
# ============================================================================
print("[3/7] Filtering low-quality samples...")

# Remove samples with library size < 1 million (outliers from QC)
low_lib_threshold = 1e6
filtered_samples = metadata_df[metadata_df['Library_Size'] >= low_lib_threshold]['Sample_ID'].tolist()

removed_samples = set(counts_df.columns) - set(filtered_samples)
print(f"Samples removed (library size < {low_lib_threshold:.0e}):")
for sample in removed_samples:
    lib_size = counts_df[sample].sum()
    print(f"  - {sample}: {lib_size:.2e} reads")

counts_filtered = counts_df[filtered_samples].copy()
metadata_filtered = metadata_df[metadata_df['Sample_ID'].isin(filtered_samples)].copy()

print(f"\n✓ Retained: {counts_filtered.shape[1]} samples\n")

# ============================================================================
# 4. SUBSET: SARS-CoV-2 vs MOCK (Primary Comparison)
# ============================================================================
print("[4/7] Subsetting SARS-CoV-2 infected vs Mock controls...")

# Focus on A549-ACE2 and Calu-3 cell lines (most relevant for COVID-19)
# These cell types support productive SARS-CoV-2 infection

subset_conditions = metadata_filtered[
    ((metadata_filtered['Virus'] == 'SARS-CoV-2') | 
     (metadata_filtered['Condition'] == 'Mock')) &
    (metadata_filtered['Cell_Type'].isin(['A549-ACE2', 'Calu3']))
]['Sample_ID'].tolist()

counts_subset = counts_filtered[subset_conditions].copy()
metadata_subset = metadata_filtered[metadata_filtered['Sample_ID'].isin(subset_conditions)].copy()

print(f"Selected cell types: A549-ACE2, Calu-3")
print(f"Sample distribution:")
print(metadata_subset.groupby(['Cell_Type', 'Condition']).size())
print(f"\n✓ Subset: {counts_subset.shape[1]} samples for primary analysis\n")

# Save subset metadata
metadata_subset.to_csv('data/processed/metadata_covid_vs_mock.csv', index=False)

# ============================================================================
# 5. FILTER LOW-EXPRESSION GENES
# ============================================================================
print("[5/7] Filtering low-expression genes...")

# Criterion: Gene must have ≥10 counts in ≥3 samples (standard RNA-seq threshold)
min_count = 10
min_samples = 3

genes_pass_filter = (counts_subset >= min_count).sum(axis=1) >= min_samples
counts_filtered_genes = counts_subset[genes_pass_filter].copy()

print(f"Filtering criteria: ≥{min_count} counts in ≥{min_samples} samples")
print(f"Genes before filtering: {counts_subset.shape[0]:,}")
print(f"Genes after filtering: {counts_filtered_genes.shape[0]:,}")
print(f"Genes removed: {counts_subset.shape[0] - counts_filtered_genes.shape[0]:,} "
      f"({100*(1 - counts_filtered_genes.shape[0]/counts_subset.shape[0]):.1f}%)\n")

# ============================================================================
# 6. NORMALIZATION: DESeq2-style Median-of-Ratios
# ============================================================================
print("[6/7] Normalizing counts (DESeq2 median-of-ratios method)...")

def deseq2_normalization(counts):
    """
    DESeq2-style normalization
    
    Rationale:
    - Robust to outlier genes (uses median, not mean)
    - Accounts for compositional bias (library size + RNA composition)
    - Preserves count nature for statistical modeling
    
    Steps:
    1. Calculate geometric mean for each gene across samples
    2. For each sample, compute ratio of counts to geometric mean
    3. Take median ratio per sample (size factor)
    4. Divide counts by size factors
    """
    # Step 1: Geometric mean per gene (exclude zeros)
    log_counts = np.log(counts)
    log_counts[np.isinf(log_counts)] = np.nan
    geometric_means = np.exp(np.nanmean(log_counts, axis=1))
    
    # Step 2: Ratio to geometric mean
    ratios = counts.div(geometric_means, axis=0)
    
    # Step 3: Median ratio per sample (size factor)
    size_factors = ratios.median(axis=0)
    
    # Step 4: Normalize
    normalized = counts.div(size_factors, axis=1)
    
    return normalized, size_factors

counts_normalized, size_factors = deseq2_normalization(counts_filtered_genes)

print(f"✓ Normalization complete")
print(f"Size factors (sample-specific scaling):")
print(size_factors.describe())

# Visualize size factors
fig, ax = plt.subplots(1, 1, figsize=(12, 5))
ax.bar(range(len(size_factors)), size_factors.values, color='teal')
ax.axhline(1.0, color='red', linestyle='--', label='No scaling')
ax.set_xlabel('Sample Index')
ax.set_ylabel('Size Factor')
ax.set_title('DESeq2 Normalization Size Factors')
ax.legend()
plt.tight_layout()
plt.savefig('results/figures/03_size_factors.png', dpi=300, bbox_inches='tight')
print("✓ Figure saved: results/figures/03_size_factors.png\n")

# ============================================================================
# 7. LOG2 TRANSFORMATION + PSEUDOCOUNT
# ============================================================================
print("[7/7] Log2-transforming normalized counts...")

# Add pseudocount to avoid log(0)
pseudocount = 1
counts_log2 = np.log2(counts_normalized + pseudocount)

print(f"✓ Log2 transformation complete (pseudocount = {pseudocount})")
print(f"Log2-transformed data range: [{counts_log2.min().min():.2f}, {counts_log2.max().max():.2f}]\n")

# ============================================================================
# 8. EXPLORATORY DATA ANALYSIS: PCA
# ============================================================================
print("[8/7 - BONUS] Principal Component Analysis (Batch Effect Detection)...")

# Standardize for PCA
scaler = StandardScaler()
counts_scaled = scaler.fit_transform(counts_log2.T)  # Transpose: samples × genes

# PCA
pca = PCA(n_components=10)
pca_coords = pca.fit_transform(counts_scaled)

# Create PCA DataFrame
pca_df = pd.DataFrame(
    pca_coords[:, :5],  # First 5 PCs
    columns=[f'PC{i+1}' for i in range(5)],
    index=counts_log2.columns
)
pca_df = pca_df.merge(metadata_subset[['Sample_ID', 'Cell_Type', 'Condition']], 
                       left_index=True, right_on='Sample_ID')

# Plot PC1 vs PC2
fig, axes = plt.subplots(1, 2, figsize=(16, 6))

# Color by condition
for condition in pca_df['Condition'].unique():
    subset = pca_df[pca_df['Condition'] == condition]
    axes[0].scatter(subset['PC1'], subset['PC2'], 
                    label=condition, s=100, alpha=0.7, edgecolors='black')
axes[0].set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% variance)')
axes[0].set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% variance)')
axes[0].set_title('PCA: Colored by Condition')
axes[0].legend()
axes[0].grid(True, alpha=0.3)

# Color by cell type
for cell_type in pca_df['Cell_Type'].unique():
    subset = pca_df[pca_df['Cell_Type'] == cell_type]
    axes[1].scatter(subset['PC1'], subset['PC2'], 
                    label=cell_type, s=100, alpha=0.7, edgecolors='black', marker='s')
axes[1].set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% variance)')
axes[1].set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% variance)')
axes[1].set_title('PCA: Colored by Cell Type')
axes[1].legend()
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('results/figures/04_pca_batch_detection.png', dpi=300, bbox_inches='tight')
print("✓ Figure saved: results/figures/04_pca_batch_detection.png")

# Scree plot (variance explained)
fig, ax = plt.subplots(1, 1, figsize=(10, 5))
variance_pct = pca.explained_variance_ratio_ * 100
ax.bar(range(1, 11), variance_pct, color='steelblue', edgecolor='black')
ax.set_xlabel('Principal Component')
ax.set_ylabel('Variance Explained (%)')
ax.set_title('PCA Scree Plot')
ax.set_xticks(range(1, 11))
plt.tight_layout()
plt.savefig('results/figures/05_pca_scree_plot.png', dpi=300, bbox_inches='tight')
print("✓ Figure saved: results/figures/05_pca_scree_plot.png\n")

# ============================================================================
# 9. SAVE PROCESSED DATA
# ============================================================================
print("[9/7] Saving processed datasets...")

# Save all processed versions
counts_filtered_genes.to_csv('data/processed/counts_filtered_raw.csv')
counts_normalized.to_csv('data/processed/counts_normalized.csv')
counts_log2.to_csv('data/processed/counts_log2_transformed.csv')

print("✓ Saved processed data files:")
print("  1. data/processed/counts_filtered_raw.csv (filtered raw counts)")
print("  2. data/processed/counts_normalized.csv (DESeq2 normalized)")
print("  3. data/processed/counts_log2_transformed.csv (log2-transformed)")
print("  4. data/processed/metadata_covid_vs_mock.csv (sample annotations)\n")

# ============================================================================
# SUMMARY REPORT
# ============================================================================
print("="*70)
print("PHASE 4 PREPROCESSING SUMMARY")
print("="*70)
print(f"Original dataset: {counts_df.shape[0]:,} genes × {counts_df.shape[1]} samples")
print(f"After quality filtering: {counts_filtered_genes.shape[0]:,} genes × {counts_subset.shape[1]} samples")
print(f"\nSample distribution (SARS-CoV-2 vs Mock):")
print(metadata_subset.groupby('Condition').size())
print(f"\nCell types included: {', '.join(metadata_subset['Cell_Type'].unique())}")
print(f"\nNormalization method: DESeq2 median-of-ratios")
print(f"Transformation: log2(normalized_counts + 1)")
print(f"\nData ready for Phase 5: Differential Expression Analysis")
print("="*70)

# Save summary
with open('results/tables/phase4_preprocessing_summary.txt', 'w') as f:
    f.write("PHASE 4 PREPROCESSING SUMMARY\n")
    f.write("="*50 + "\n")
    f.write(f"Filtered genes: {counts_filtered_genes.shape[0]:,}\n")
    f.write(f"Samples analyzed: {counts_subset.shape[1]}\n")
    f.write(f"Conditions: SARS-CoV-2 infected vs Mock\n")
    f.write(f"Cell types: A549-ACE2, Calu-3\n")
    f.write(f"Normalization: DESeq2 median-of-ratios\n")

print("\n✓ Summary saved: results/tables/phase4_preprocessing_summary.txt")