"""
Gene Network Analysis & Visualization
Phase 7: Build gene interaction networks with fallback mechanisms

Methods:
1. STRING database API (primary)
2. Gene co-expression network (fallback)
3. Literature-based core COVID-19 interactions (fallback)
4. Hub gene identification via centrality
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from scipy.stats import pearsonr
import requests
import warnings
warnings.filterwarnings('ignore')

# Set plotting style
sns.set_style("white")
plt.rcParams['figure.figsize'] = (16, 12)

print("="*70)
print("PHASE 7: Gene Network Analysis")
print("="*70 + "\n")

# ============================================================================
# 1. LOAD DATA
# ============================================================================
print("[1/7] Loading data...")

deg_sig = pd.read_csv('results/tables/deg_significant_only.csv')
counts_log2 = pd.read_csv('data/processed/counts_log2_transformed.csv', index_col=0)

# Focus on top DEGs for clearer visualization
top_degs = deg_sig.nlargest(80, 'Log2FoldChange', keep='all')
deg_genes = top_degs['Gene'].tolist()

print(f"✓ Network analysis for top {len(deg_genes)} DEGs")
print(f"  - Upregulated: {(top_degs['Significance'] == 'Upregulated').sum()}")
print(f"  - Downregulated: {(top_degs['Significance'] == 'Downregulated').sum()}")
print()

# ============================================================================
# 2. ATTEMPT STRING DATABASE QUERY
# ============================================================================
print("[2/7] Attempting STRING database query...")

def get_string_interactions(genes, species=9606, score_threshold=500):
    """
    Retrieve protein interactions from STRING with robust error handling
    """
    try:
        string_api_url = "https://string-db.org/api/tsv/network"
        
        params = {
            "identifiers": "%0d".join(genes),
            "species": species,
            "caller_identity": "transcriptomics_analysis"
        }
        
        print(f"  Querying STRING for {len(genes)} genes...")
        response = requests.post(string_api_url, data=params, timeout=30)
        
        if response.ok and len(response.text) > 100:
            lines = response.text.strip().split('\n')
            
            interactions = []
            for line in lines[1:]:  # Skip header
                parts = line.split('\t')
                if len(parts) >= 6:
                    try:
                        score = float(parts[5])
                        if score >= score_threshold:
                            interactions.append({
                                'node1': parts[2],
                                'node2': parts[3],
                                'score': score
                            })
                    except (ValueError, IndexError):
                        continue
            
            if len(interactions) > 0:
                print(f"  ✓ STRING API success: {len(interactions)} interactions")
                return pd.DataFrame(interactions)
        
        print("  ⚠ STRING API returned no valid interactions")
        return pd.DataFrame()
        
    except Exception as e:
        print(f"  ✗ STRING API error: {e}")
        return pd.DataFrame()

interactions_df = get_string_interactions(deg_genes, score_threshold=500)

# ============================================================================
# 3. FALLBACK: GENE CO-EXPRESSION NETWORK
# ============================================================================
if interactions_df.empty:
    print("\n[3/7] Building gene co-expression network (STRING fallback)...")
    
    # Calculate pairwise correlations
    gene_subset = [g for g in deg_genes if g in counts_log2.index]
    expr_matrix = counts_log2.loc[gene_subset].T  # Samples × Genes
    
    interactions = []
    correlation_threshold = 0.7
    
    print(f"  Computing correlations for {len(gene_subset)} genes...")
    
    for i, gene1 in enumerate(gene_subset):
        for gene2 in gene_subset[i+1:]:
            try:
                corr, pval = pearsonr(expr_matrix[gene1], expr_matrix[gene2])
                if abs(corr) >= correlation_threshold and pval < 0.05:
                    interactions.append({
                        'node1': gene1,
                        'node2': gene2,
                        'score': abs(corr) * 1000  # Scale to STRING-like scores
                    })
            except:
                continue
    
    if len(interactions) > 0:
        interactions_df = pd.DataFrame(interactions)
        print(f"  ✓ Co-expression network: {len(interactions_df)} edges (|r| ≥ {correlation_threshold})")
    else:
        print("  ⚠ Low correlation - using literature-based core network")
        
        # FALLBACK 2: Known COVID-19 gene interactions
        core_interactions = [
            ('IFNB1', 'STAT1'), ('IFNB1', 'STAT2'), ('IFNB1', 'IRF9'),
            ('STAT1', 'STAT2'), ('STAT1', 'IRF9'), ('STAT2', 'IRF9'),
            ('TNF', 'NFKB1'), ('TNF', 'NFKBIA'), ('TNF', 'TNFAIP3'),
            ('IL6', 'STAT3'), ('IL6', 'JAK2'), ('IL6', 'SOCS3'),
            ('IFIT1', 'IFIT2'), ('IFIT1', 'IFIT3'), ('IFIT2', 'IFIT3'),
            ('CXCL2', 'CXCL3'), ('CXCL2', 'CXCL8'), ('CCL20', 'CCL5'),
            ('EGR1', 'JUN'), ('EGR1', 'FOS'), ('JUN', 'FOS'),
            ('IFNL1', 'IFNL2'), ('IFNL1', 'IFNL3'), ('IFNL2', 'IFNL3'),
            ('NFKB1', 'NFKBIA'), ('NFKB1', 'RELA')
        ]
        
        interactions_df = pd.DataFrame([
            {'node1': g1, 'node2': g2, 'score': 700}
            for g1, g2 in core_interactions
            if g1 in deg_genes and g2 in deg_genes
        ])
        
        print(f"  ✓ Literature network: {len(interactions_df)} core interactions")
else:
    print("\n[3/7] Using STRING database interactions")

print()

# ============================================================================
# 4. BUILD NETWORK GRAPH
# ============================================================================
print("[4/7] Constructing network graph...")

G = nx.Graph()

if not interactions_df.empty:
    # Add weighted edges
    for _, row in interactions_df.iterrows():
        G.add_edge(row['node1'], row['node2'], weight=row['score']/1000)
    
    # Keep only DEGs in network
    deg_genes_set = set(deg_genes)
    nodes_to_remove = [n for n in G.nodes() if n not in deg_genes_set]
    G.remove_nodes_from(nodes_to_remove)
    
    # Keep largest connected component
    if len(G.nodes()) > 1:
        components = list(nx.connected_components(G))
        if len(components) > 1:
            largest_cc = max(components, key=len)
            G = G.subgraph(largest_cc).copy()
            print(f"  (Keeping largest component: {len(largest_cc)} nodes)")

print(f"✓ Network constructed:")
print(f"  - Nodes (genes): {G.number_of_nodes()}")
print(f"  - Edges (interactions): {G.number_of_edges()}")

if G.number_of_nodes() > 0:
    print(f"  - Density: {nx.density(G):.3f}")
    print(f"  - Avg degree: {sum(dict(G.degree()).values())/G.number_of_nodes():.1f}")
print()

# ============================================================================
# 5. CENTRALITY ANALYSIS (if network exists)
# ============================================================================
if G.number_of_nodes() > 0:
    print("[5/7] Computing network centrality...")
    
    degree_centrality = nx.degree_centrality(G)
    betweenness_centrality = nx.betweenness_centrality(G)
    
    centrality_df = pd.DataFrame({
        'Gene': list(degree_centrality.keys()),
        'Degree_Centrality': list(degree_centrality.values()),
        'Betweenness_Centrality': list(betweenness_centrality.values())
    })
    
    # Merge with DEG data
    centrality_df = centrality_df.merge(
        top_degs[['Gene', 'Log2FoldChange', 'FDR', 'Significance']], 
        on='Gene', how='left'
    )
    centrality_df = centrality_df.sort_values('Degree_Centrality', ascending=False)
    
    print(f"Top 15 Hub Genes (by Degree Centrality):")
    print(centrality_df.head(15)[['Gene', 'Degree_Centrality', 'Log2FoldChange']].to_string(index=False))
    print()
    
    # Save
    centrality_df.to_csv('results/tables/network_hub_genes.csv', index=False)
    print("✓ Hub genes saved: results/tables/network_hub_genes.csv\n")
    
else:
    print("[5/7] Skipping centrality (empty network)\n")
    centrality_df = pd.DataFrame(columns=['Gene', 'Degree_Centrality', 'Log2FoldChange'])

# ============================================================================
# 6. NETWORK VISUALIZATION
# ============================================================================
print("[6/7] Generating network visualization...")

if G.number_of_nodes() >= 3:
    fig, ax = plt.subplots(1, 1, figsize=(18, 16))
    
    # Layout
    if G.number_of_nodes() < 50:
        pos = nx.spring_layout(G, k=2, iterations=100, seed=42)
    else:
        pos = nx.spring_layout(G, k=1.5, iterations=50, seed=42)
    
    # Node properties
    if len(centrality_df) > 0:
        node_sizes = [degree_centrality.get(node, 0.01) * 6000 for node in G.nodes()]
    else:
        node_sizes = [500 for _ in G.nodes()]
    
    gene_to_fc = dict(zip(top_degs['Gene'], top_degs['Log2FoldChange']))
    node_colors = [gene_to_fc.get(node, 0) for node in G.nodes()]
    
    # Draw edges
    nx.draw_networkx_edges(G, pos, alpha=0.25, width=1, edge_color='gray', ax=ax)
    
    # Draw nodes
    nodes = nx.draw_networkx_nodes(
        G, pos, node_size=node_sizes, node_color=node_colors,
        cmap='RdBu_r', vmin=-3, vmax=6, alpha=0.85,
        edgecolors='black', linewidths=1.5, ax=ax
    )
    
    # Label top hubs
    if len(centrality_df) > 0:
        top_hubs = centrality_df.head(min(20, len(centrality_df)))['Gene'].tolist()
        labels = {node: node if node in top_hubs else '' for node in G.nodes()}
        nx.draw_networkx_labels(G, pos, labels, font_size=9, font_weight='bold', ax=ax)
    
    # Colorbar
    cbar = plt.colorbar(nodes, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('Log2 Fold Change', fontsize=14, fontweight='bold')
    
    ax.set_title(
        f'Protein-Protein Interaction Network\nSARS-CoV-2 DEGs ({G.number_of_nodes()} genes, {G.number_of_edges()} interactions)',
        fontsize=16, fontweight='bold', pad=20
    )
    ax.axis('off')
    
    plt.tight_layout()
    plt.savefig('results/figures/08_network_visualization.png', dpi=300, bbox_inches='tight')
    print("✓ Figure saved: results/figures/08_network_visualization.png")
    plt.close()
else:
    print("⚠ Network too small for visualization (skipping)")

# ============================================================================
# 7. HUB SUBNETWORK
# ============================================================================
print("\n[7/7] Creating hub gene subnetwork...")

if len(centrality_df) >= 15 and G.number_of_nodes() >= 15:
    n_hubs = min(30, len(centrality_df))
    top_hubs_list = centrality_df.head(n_hubs)['Gene'].tolist()
    G_sub = G.subgraph(top_hubs_list).copy()
    
    if G_sub.number_of_nodes() >= 3:
        fig, ax = plt.subplots(1, 1, figsize=(14, 14))
        pos_sub = nx.spring_layout(G_sub, k=2, iterations=100, seed=42)
        
        node_sizes_sub = [degree_centrality.get(node, 0.01) * 8000 for node in G_sub.nodes()]
        node_colors_sub = [gene_to_fc.get(node, 0) for node in G_sub.nodes()]
        
        nx.draw_networkx_edges(G_sub, pos_sub, alpha=0.3, width=2.5, edge_color='gray', ax=ax)
        nodes_sub = nx.draw_networkx_nodes(
            G_sub, pos_sub, node_size=node_sizes_sub, node_color=node_colors_sub,
            cmap='RdBu_r', vmin=-3, vmax=6, alpha=0.9,
            edgecolors='black', linewidths=2, ax=ax
        )
        nx.draw_networkx_labels(G_sub, pos_sub, font_size=10, font_weight='bold', ax=ax)
        
        cbar = plt.colorbar(nodes_sub, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label('Log2 Fold Change', fontsize=12, fontweight='bold')
        
        ax.set_title(
            f'Hub Gene Subnetwork (Top {n_hubs} Central Genes)\nSARS-CoV-2 Response',
            fontsize=16, fontweight='bold', pad=20
        )
        ax.axis('off')
        
        plt.tight_layout()
        plt.savefig('results/figures/09_hub_subnetwork.png', dpi=300, bbox_inches='tight')
        print("✓ Figure saved: results/figures/09_hub_subnetwork.png")
    else:
        print("⚠ Hub subnetwork too small (skipping)")
else:
    print("⚠ Not enough hub genes for subnetwork (skipping)")

# ============================================================================
# SUMMARY
# ============================================================================
print("\n" + "="*70)
print("PHASE 7 NETWORK ANALYSIS SUMMARY")
print("="*70)

if G.number_of_nodes() > 0:
    print(f"Network nodes: {G.number_of_nodes()} genes")
    print(f"Network edges: {G.number_of_edges()} interactions")
    print(f"Network density: {nx.density(G):.3f}")
    
    if len(centrality_df) > 0:
        top_hub = centrality_df.iloc[0]
        print(f"\nTop hub gene: {top_hub['Gene']}")
        print(f"  - Centrality: {top_hub['Degree_Centrality']:.3f}")
        print(f"  - Log2FC: {top_hub['Log2FoldChange']:.2f}")
        
        print(f"\nHub genes (top 10): {', '.join(centrality_df.head(10)['Gene'].tolist())}")
else:
    print("Network: Empty (API failure - manual review recommended)")
    print("\nKey DEGs identified:")
    print(f"  {', '.join(deg_genes[:15])}")

print(f"\nNext Step: Phase 8 - Documentation & Final Report")
print("="*70)

# Save summary
with open('results/tables/phase7_network_summary.txt', 'w', encoding='utf-8') as f:
    f.write("PHASE 7 NETWORK ANALYSIS SUMMARY\n")
    f.write("="*50 + "\n")
    
    if G.number_of_nodes() > 0:
        f.write(f"Network nodes: {G.number_of_nodes()}\n")
        f.write(f"Network edges: {G.number_of_edges()}\n")
        f.write(f"Density: {nx.density(G):.3f}\n")
        
        if len(centrality_df) > 0:
            f.write(f"Top hub: {centrality_df.iloc[0]['Gene']}\n")
            f.write(f"Top 10 hubs: {', '.join(centrality_df.head(10)['Gene'].tolist())}\n")
    else:
        f.write("Network construction failed - using fallback gene list\n")
        f.write(f"Top DEGs: {', '.join(deg_genes[:20])}\n")

print("\n✓ Summary saved: results/tables/phase7_network_summary.txt")