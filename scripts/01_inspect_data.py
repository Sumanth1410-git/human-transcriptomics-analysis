"""
Data Inspection & Initial QC
GSE147507: SARS-CoV-2 Transcriptomics
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Set plotting style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (12, 6)

print("="*70)
print("PHASE 3: Data Inspection & Quality Control")
print("Dataset: GSE147507 - SARS-CoV-2 Infection Transcriptomics")
print("="*70 + "\n")

# ============================================================================
# 1. LOAD RAW COUNT MATRIX
# ============================================================================
print("[1/5] Loading raw count matrix...")

try:
    counts_df = pd.read_csv(
        'data/raw/covid19_raw_counts.tsv',
        sep='\t',
        index_col=0  # First column is gene IDs
    )
    print(f"✓ Count matrix loaded successfully")
    print(f"  Shape: {counts_df.shape[0]} genes × {counts_df.shape[1]} samples\n")
    
except FileNotFoundError:
    print("✗ Error: covid19_raw_counts.tsv not found in data/raw/")
    print("  Please verify file extraction")
    exit(1)

# ============================================================================
# 2. BASIC DATA STRUCTURE
# ============================================================================
print("[2/5] Data Structure Overview:")
print("-" * 70)
print(counts_df.head(10))
print("\nColumn names (sample IDs):")
print(counts_df.columns.tolist())
print(f"\nData types:\n{counts_df.dtypes.value_counts()}")

# Check for missing values
missing = counts_df.isnull().sum().sum()
print(f"\nMissing values: {missing}")

if missing > 0:
    print("⚠ Warning: Missing values detected - will handle in preprocessing")

# ============================================================================
# 3. LIBRARY SIZE DISTRIBUTION (Total Reads per Sample)
# ============================================================================
print("\n[3/5] Library Size Analysis (Total Counts per Sample):")
print("-" * 70)

library_sizes = counts_df.sum(axis=0)
print(library_sizes.describe())

# Visualize
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Bar plot
axes[0].bar(range(len(library_sizes)), library_sizes.values, color='steelblue')
axes[0].set_xlabel('Sample Index')
axes[0].set_ylabel('Total Read Count')
axes[0].set_title('Library Sizes Across Samples')
axes[0].axhline(library_sizes.mean(), color='red', linestyle='--', 
                label=f'Mean: {library_sizes.mean():.2e}')
axes[0].legend()
axes[0].tick_params(axis='x', rotation=90)

# Box plot
axes[1].boxplot(library_sizes.values)
axes[1].set_ylabel('Total Read Count')
axes[1].set_title('Library Size Distribution')
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('results/figures/01_library_sizes.png', dpi=300, bbox_inches='tight')
print("✓ Figure saved: results/figures/01_library_sizes.png")

# ============================================================================
# 4. GENE EXPRESSION DISTRIBUTION
# ============================================================================
print("\n[4/5] Gene Expression Distribution:")
print("-" * 70)

# Count genes with zero expression across all samples
zero_genes = (counts_df.sum(axis=1) == 0).sum()
print(f"Genes with zero counts across all samples: {zero_genes}")

# Mean expression per gene
mean_expression = counts_df.mean(axis=1)
print(f"\nMean expression per gene (summary):")
print(mean_expression.describe())

# Visualize
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Histogram of mean expression (log scale)
axes[0].hist(np.log10(mean_expression[mean_expression > 0] + 1), 
             bins=50, color='coral', edgecolor='black')
axes[0].set_xlabel('log10(Mean Expression + 1)')
axes[0].set_ylabel('Number of Genes')
axes[0].set_title('Distribution of Gene Expression Levels')
axes[0].grid(True, alpha=0.3)

# Count distribution per sample (first 5 samples)
sample_subset = counts_df.iloc[:, :5]
sample_subset_log = np.log10(sample_subset + 1)
sample_subset_log.boxplot(ax=axes[1])
axes[1].set_ylabel('log10(Count + 1)')
axes[1].set_title('Expression Distribution (First 5 Samples)')
axes[1].tick_params(axis='x', rotation=45)

plt.tight_layout()
plt.savefig('results/figures/02_expression_distribution.png', dpi=300, bbox_inches='tight')
print("✓ Figure saved: results/figures/02_expression_distribution.png")

# ============================================================================
# 5. PARSE SAMPLE METADATA
# ============================================================================
print("\n[5/5] Parsing Sample Metadata:")
print("-" * 70)

try:
    # Series matrix files have a specific format - need custom parsing
    with open('data/metadata/covid19_sample_metadata.txt', 'r') as f:
        lines = f.readlines()
    
    # Find sample characteristic lines
    metadata_dict = {}
    for line in lines:
        if line.startswith('!Sample_'):
            parts = line.strip().split('\t')
            key = parts[0].replace('!Sample_', '')
            values = parts[1:]
            metadata_dict[key] = values
    
    # Create metadata DataFrame
    metadata_df = pd.DataFrame(metadata_dict)
    
    # Display key metadata fields
    print("Available metadata fields:")
    print(metadata_df.columns.tolist())
    
    if 'characteristics_ch1' in metadata_df.columns:
        print("\nSample characteristics (first 5):")
        print(metadata_df[['geo_accession', 'characteristics_ch1']].head())
    
    # Save cleaned metadata
    metadata_df.to_csv('data/processed/sample_metadata_cleaned.csv', index=False)
    print("\n✓ Metadata saved: data/processed/sample_metadata_cleaned.csv")
    
except Exception as e:
    print(f"⚠ Warning: Could not parse metadata file automatically")
    print(f"  Error: {e}")
    print("  We'll manually curate metadata in next phase")

# ============================================================================
# SUMMARY REPORT
# ============================================================================
print("\n" + "="*70)
print("PHASE 3 SUMMARY")
print("="*70)
print(f"✓ Count matrix dimensions: {counts_df.shape}")
print(f"✓ Genes: {counts_df.shape[0]:,}")
print(f"✓ Samples: {counts_df.shape[1]}")
print(f"✓ Missing values: {missing}")
print(f"✓ Zero-count genes: {zero_genes:,} ({100*zero_genes/counts_df.shape[0]:.2f}%)")
print(f"✓ Mean library size: {library_sizes.mean():.2e} ± {library_sizes.std():.2e}")
print("\nQC Figures Generated:")
print("  1. results/figures/01_library_sizes.png")
print("  2. results/figures/02_expression_distribution.png")
print("\nNext Step: Review QC plots, then proceed to Phase 4 (Normalization)")
print("="*70)

# Save summary to file
with open('results/tables/phase3_qc_summary.txt', 'w') as f:
    f.write(f"GSE147507 Quality Control Summary\n")
    f.write(f"{'='*50}\n")
    f.write(f"Genes: {counts_df.shape[0]:,}\n")
    f.write(f"Samples: {counts_df.shape[1]}\n")
    f.write(f"Missing values: {missing}\n")
    f.write(f"Zero-count genes: {zero_genes:,}\n")
    f.write(f"Mean library size: {library_sizes.mean():.2e}\n")

print("\n✓ Summary saved: results/tables/phase3_qc_summary.txt")