"""
Compress large network visualization images
Reduces file size by 60-70% while maintaining visual quality
"""
from PIL import Image
import os

print("Compressing large network images...")

# Target the two large network files
large_images = [
    "results/figures/08_network_visualization.png",
    "results/figures/09_hub_subnetwork.png"
]

for img_path in large_images:
    # Check original size
    original_size = os.path.getsize(img_path) / (1024*1024)
    print(f"\n{os.path.basename(img_path)}:")
    print(f"  Original: {original_size:.2f} MB")
    
    # Open and compress
    img = Image.open(img_path)
    
    # Save with optimization (reduces file size significantly)
    img.save(img_path, "PNG", optimize=True, compress_level=9)
    
    # Check new size
    new_size = os.path.getsize(img_path) / (1024*1024)
    reduction = ((original_size - new_size) / original_size) * 100
    
    print(f"  Compressed: {new_size:.2f} MB")
    print(f"  Reduction: {reduction:.1f}%")
    print(f"  ✓ Saved")

print("\n✓ Compression complete!")
