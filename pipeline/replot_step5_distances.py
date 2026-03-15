"""Replot step5_expansion_distances.png with µm units and peak annotation.

Usage:
    python pipeline/replot_step5_distances.py
"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial import cKDTree

# --- Paths (from config.yaml) ---
BASE = "/data/project/lyrak/masld-xenium/xenium-output/Xenium_V1_FFPE_Human_Brain_Alzheimers_With_Addon_outs"
STEP5_DIR = os.path.join(BASE, "step5_optimal_expansion")
CSV_PATH = os.path.join(STEP5_DIR, "human_alzheimers_step5_expanded_transcripts.csv")
OUT_PATH = os.path.join(STEP5_DIR, "human_alzheimers_step5_expansion_distances.png")

# Xenium: 4.70588 pixels per µm
PX_PER_UM = 4.70588
SUBSAMPLE_FRAC = 0.01

print(f"Loading expanded transcripts from: {CSV_PATH}")
print(f"  (this may take a moment for ~16M rows...)")

# Only load columns we need
cols = ['cell_id', 'x_location', 'y_location', 'domain']
df = pd.read_csv(CSV_PATH, usecols=cols, low_memory=False)
print(f"  Loaded {len(df):,} transcripts")

# Separate assigned vs unassigned (original cell_id)
if df['cell_id'].dtype in [np.int64, np.int32, np.float64, np.float32]:
    mask_unassigned = df['cell_id'] == -1
else:
    mask_unassigned = df['cell_id'].astype(str).str.upper().isin(['-1', 'UNASSIGNED', 'NAN', ''])

assigned = df[~mask_unassigned & df['domain'].notna()]
unassigned = df[mask_unassigned]
print(f"  Assigned: {len(assigned):,}  |  Unassigned: {len(unassigned):,}")

# Build KDTree from subsample of assigned reads
n_sample = max(int(len(assigned) * SUBSAMPLE_FRAC), 1000)
n_sample = min(n_sample, len(assigned))
anchors = assigned.sample(n=n_sample, random_state=42)
print(f"  Building KDTree with {n_sample:,} anchors...")

tree = cKDTree(anchors[['x_location', 'y_location']].values)

# Query unassigned reads
coords = unassigned[['x_location', 'y_location']].values
print(f"  Querying {len(coords):,} unassigned reads...")

chunk_size = 1_000_000
all_dists = []
n_chunks = int(np.ceil(len(coords) / chunk_size))
for i in range(n_chunks):
    start = i * chunk_size
    end = min((i + 1) * chunk_size, len(coords))
    dists, _ = tree.query(coords[start:end], k=1)
    all_dists.append(dists)
    print(f"    chunk {i+1}/{n_chunks} done")

concat_distances = np.concatenate(all_dists)

# Convert to µm
distances_um = concat_distances / PX_PER_UM
print(f"\n  Distance stats (µm):")
print(f"    min={distances_um.min():.1f}, median={np.median(distances_um):.1f}, "
      f"mean={distances_um.mean():.1f}, max={distances_um.max():.1f}")

# Clip to 99th percentile for readable x-axis
p99 = np.percentile(distances_um, 99)
xlim_max = min(p99 * 1.5, 50)  # cap at 50 µm
distances_plot = distances_um[distances_um <= xlim_max]
n_clipped = len(distances_um) - len(distances_plot)
print(f"  Plot x-limit: {xlim_max:.1f} µm  (clipped {n_clipped:,} reads > {xlim_max:.1f} µm)")

# --- Plot ---
fig, ax = plt.subplots(figsize=(8, 5))
counts, bin_edges, patches = ax.hist(distances_plot, bins=50, color='orange', alpha=0.7, edgecolor='white')
sns.kdeplot(distances_plot, color='darkorange', linewidth=1.5, ax=ax)

# Peak annotation
peak_idx = int(np.argmax(counts))
peak_count = int(counts[peak_idx])
peak_dist = (bin_edges[peak_idx] + bin_edges[peak_idx + 1]) / 2

ax.annotate(
    f'Peak: {peak_dist:.1f} µm, n={peak_count:,}',
    xy=(peak_dist, peak_count),
    xytext=(peak_dist + xlim_max * 0.25, peak_count * 0.75),
    arrowprops=dict(arrowstyle='->', color='black', lw=1.2),
    fontsize=10, fontweight='bold',
    bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='gray', alpha=0.9)
)

# Median line
med = np.median(distances_um)
ax.axvline(med, color='red', linestyle='--', linewidth=1, alpha=0.7)
ax.text(med + 0.5, ax.get_ylim()[1] * 0.95, f'Median: {med:.1f} µm',
        fontsize=9, color='red', va='top')

# Reference lines from our dataset's actual turnover analysis
OPT_TXT = os.path.join(STEP5_DIR, "human_alzheimers_step5_optimal_expansion.txt")
ref = {}
if os.path.exists(OPT_TXT):
    for line in open(OPT_TXT):
        k, v = line.strip().split('\t')
        try:
            ref[k] = float(v)
        except ValueError:
            ref[k] = v

opt_exp = ref.get('optimal_expansion', 4.10)
turnover = ref.get('median_turnover', 8.90)
nuc_size = ref.get('median_nuclei_size', 4.80)
ymax = ax.get_ylim()[1]

ax.axvline(opt_exp, color='blue', linestyle=':', linewidth=1.2, alpha=0.7)
ax.text(opt_exp + 0.3, ymax * 0.80, f'Optimal exp.\n({opt_exp:.1f} µm)',
        fontsize=8, color='blue', va='top')

ax.axvline(turnover, color='green', linestyle=':', linewidth=1.2, alpha=0.7)
ax.text(turnover + 0.3, ymax * 0.65, f'Turnover dist.\n({turnover:.1f} µm)',
        fontsize=8, color='green', va='top')

ax.axvline(15.0, color='gray', linestyle=':', linewidth=1, alpha=0.5)
ax.text(15.0 + 0.3, ymax * 0.50, 'Xenium default\n(15 µm)',
        fontsize=8, color='gray', va='top')

ax.set_title("Distance to Nearest Domain-Anchor (Unassigned Reads)", fontsize=13)
ax.set_xlabel("Distance (µm)", fontsize=11)
ax.set_ylabel("Count", fontsize=11)
ax.set_xlim(0, xlim_max)

fig.tight_layout()
fig.savefig(OUT_PATH, dpi=300, bbox_inches='tight')
plt.close(fig)
print(f"\n  Saved: {OUT_PATH}")
