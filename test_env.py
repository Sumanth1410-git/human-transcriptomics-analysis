"""
Environment Validation Script
Tests library imports and GPU availability
"""

import sys
print(f"Python Version: {sys.version}\n")

# Core libraries
try:
    import numpy as np
    import pandas as pd
    import scipy
    import matplotlib.pyplot as plt
    import seaborn as sns
    print("✓ Core data science libraries imported successfully")
except ImportError as e:
    print(f"✗ Error importing core libraries: {e}")

# Bioinformatics
try:
    from Bio import SeqIO, Entrez
    import scanpy as sc
    import networkx as nx
    import statsmodels.api as sm
    print("✓ Bioinformatics libraries imported successfully")
except ImportError as e:
    print(f"✗ Error importing bioinformatics libraries: {e}")

# GPU Check (Optional for future use)
try:
    import torch
    if torch.cuda.is_available():
        print(f"✓ CUDA available: {torch.cuda.get_device_name(0)}")
    else:
        print("⚠ CUDA not available (CPU mode)")
except ImportError:
    print("⚠ PyTorch not installed (optional for deep learning)")

print("\n--- Environment Ready for Phase 2 ---")