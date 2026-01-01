import visualize_results as vis
import sys

print("Checking for visualize_marker_genes_violin...")
if hasattr(vis, 'visualize_marker_genes_violin'):
    print("SUCCESS: Function found!")
else:
    print("FAILURE: Function NOT found.")
    # Print dir(vis) to see what's there
    print(dir(vis))
