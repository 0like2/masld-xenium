# pipeline/xenium_eda.py
# ============================================================
# Xenium Data Exploratory Data Analysis (EDA)
# ============================================================
# Standalone script to inspect raw Xenium data and pipeline outputs.
# Provides: file inventory, data head/info, cell_id formats,
# QC control vs raw distributions, coordinate ranges, DAPI metadata,
# negative control analysis, and RAM usage estimation.
#
# Usage:
#   python pipeline/xenium_eda.py                       # default config
#   python pipeline/xenium_eda.py --config path/to/config.yaml
#   python pipeline/xenium_eda.py --skip-plots           # text-only mode
#   python pipeline/xenium_eda.py --section all          # run all sections
#   python pipeline/xenium_eda.py --section files,qc     # specific sections

import os
import sys
import argparse
import json
import logging
import time
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

logging.basicConfig(level=logging.INFO, format='%(asctime)s [EDA] %(message)s')
logger = logging.getLogger(__name__)

# Optional imports (plots)
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import seaborn as sns
    HAS_PLOTS = True
except ImportError:
    HAS_PLOTS = False

try:
    import scanpy as sc
    HAS_SCANPY = True
except ImportError:
    HAS_SCANPY = False

try:
    import tifffile as tf
    HAS_TIFF = True
except ImportError:
    HAS_TIFF = False

try:
    import psutil
    HAS_PSUTIL = True
except ImportError:
    HAS_PSUTIL = False


# ============================================================
# Utility helpers
# ============================================================

def _sizeof_fmt(num_bytes):
    """Human-readable file size."""
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if abs(num_bytes) < 1024.0:
            return f"{num_bytes:.1f} {unit}"
        num_bytes /= 1024.0
    return f"{num_bytes:.1f} PB"


def _separator(title):
    w = 70
    print("\n" + "=" * w)
    print(f"  {title}")
    print("=" * w)


# ============================================================
# Section 1: File Inventory
# ============================================================

def section_files(input_path, output_dir, sample_name, **kwargs):
    """List all input files with sizes and identify key data files."""
    _separator("1. FILE INVENTORY")

    print(f"\n[Input Directory] {input_path}")
    if not os.path.isdir(input_path):
        print("  ERROR: Input directory does not exist!")
        return

    total_size = 0
    file_list = []
    for f in sorted(os.listdir(input_path)):
        fp = os.path.join(input_path, f)
        if os.path.isfile(fp):
            sz = os.path.getsize(fp)
            total_size += sz
            file_list.append((f, sz))
        elif os.path.isdir(fp):
            dir_size = sum(
                os.path.getsize(os.path.join(dp, fn))
                for dp, _, fns in os.walk(fp) for fn in fns
            )
            total_size += dir_size
            file_list.append((f + "/", dir_size))

    print(f"\n{'File/Dir':<55} {'Size':>12}")
    print("-" * 67)
    for name, sz in file_list:
        print(f"  {name:<53} {_sizeof_fmt(sz):>12}")
    print("-" * 67)
    print(f"  {'TOTAL':<53} {_sizeof_fmt(total_size):>12}")

    # Key files check
    print("\n[Key Files Status]")
    key_files = [
        ("transcripts.parquet", "Transcript data (fast)"),
        ("transcripts.csv", "Transcript data (CSV)"),
        ("cells.csv", "Cell metadata"),
        ("cell_feature_matrix/", "Count matrix"),
        ("morphology_focus.ome.tif", "DAPI image (focus)"),
        ("morphology_mip.ome.tif", "DAPI image (MIP)"),
        ("morphology.ome.tif", "DAPI image (full stack)"),
        ("gene_panel.json", "Gene panel definition"),
        ("nucleus_boundaries.parquet", "Nucleus boundaries"),
        ("cell_boundaries.parquet", "Cell boundaries"),
        ("experiment.xenium", "Experiment metadata"),
        ("metrics_summary.csv", "10X QC metrics"),
        ("analysis_summary.html", "10X analysis report"),
    ]
    for fname, desc in key_files:
        fp = os.path.join(input_path, fname)
        exists = os.path.exists(fp)
        status = "OK" if exists else "MISSING"
        print(f"  [{status:>7}] {fname:<40} - {desc}")

    # Pipeline output check
    output_sample_dir = os.path.join(output_dir, sample_name)
    if os.path.isdir(output_sample_dir):
        print(f"\n[Pipeline Output] {output_sample_dir}")
        for d in sorted(os.listdir(output_sample_dir)):
            dp = os.path.join(output_sample_dir, d)
            if os.path.isdir(dp):
                n_files = sum(1 for _ in os.listdir(dp) if os.path.isfile(os.path.join(dp, _)))
                dir_sz = sum(
                    os.path.getsize(os.path.join(dp, fn))
                    for fn in os.listdir(dp) if os.path.isfile(os.path.join(dp, fn))
                )
                print(f"  {d + '/':<40} {n_files:>4} files  {_sizeof_fmt(dir_sz):>10}")


# ============================================================
# Section 2: Transcripts Data Inspection
# ============================================================

def section_transcripts(input_path, **kwargs):
    """Inspect transcripts.csv/parquet: head, dtypes, shape, coordinate ranges."""
    _separator("2. TRANSCRIPTS DATA")

    parquet_path = os.path.join(input_path, "transcripts.parquet")
    csv_path = os.path.join(input_path, "transcripts.csv")

    if os.path.exists(parquet_path):
        print(f"\n[Loading] {parquet_path} (parquet - fast)")
        spots = pd.read_parquet(parquet_path)
    elif os.path.exists(csv_path):
        print(f"\n[Loading] {csv_path} (CSV - slower)")
        spots = pd.read_csv(csv_path, nrows=500000)
        print(f"  NOTE: Only loaded first 500K rows for preview. Full file has more rows.")
    else:
        print("  ERROR: No transcripts file found!")
        return

    # Decode bytes columns if needed (parquet sometimes stores as bytes)
    for col in spots.columns:
        if spots[col].dtype == object:
            try:
                sample = spots[col].dropna().iloc[0]
                if isinstance(sample, bytes):
                    spots[col] = spots[col].str.decode('utf-8')
            except (IndexError, AttributeError):
                pass

    print(f"\n[Shape] {spots.shape[0]:,} rows x {spots.shape[1]} columns")

    print(f"\n[Columns & Dtypes]")
    for col in spots.columns:
        n_null = spots[col].isnull().sum()
        n_unique = spots[col].nunique()
        print(f"  {col:<25} dtype={str(spots[col].dtype):<12} nulls={n_null:>8,}  unique={n_unique:>10,}")

    print(f"\n[Head - First 5 Rows]")
    print(spots.head().to_string(max_colwidth=30))

    print(f"\n[Tail - Last 5 Rows]")
    print(spots.tail().to_string(max_colwidth=30))

    # Cell ID analysis
    if 'cell_id' in spots.columns:
        print(f"\n[Cell ID Analysis]")
        cid = spots['cell_id']
        print(f"  dtype: {cid.dtype}")
        print(f"  unique values: {cid.nunique():,}")
        print(f"  sample values: {list(cid.dropna().unique()[:10])}")

        # Check for UNASSIGNED
        if cid.dtype == object or cid.dtype.name == 'category':
            unassigned = (cid == 'UNASSIGNED').sum()
            print(f"  UNASSIGNED transcripts: {unassigned:,} ({unassigned/len(spots)*100:.1f}%)")
        else:
            zero_cells = (cid == 0).sum()
            neg_cells = (cid < 0).sum()
            print(f"  cell_id == 0 (unassigned): {zero_cells:,} ({zero_cells/len(spots)*100:.1f}%)")
            if neg_cells > 0:
                print(f"  cell_id < 0: {neg_cells:,}")

        # Float-string issue check
        if cid.dtype == np.float64:
            print(f"  WARNING: cell_id is float64! (e.g., 1.0 instead of '1')")
            print(f"  This causes dtype mismatch in P2R merging. Pipeline _normalize_cell_ids() handles this.")
        elif cid.dtype == object:
            has_float = cid.dropna().str.contains(r'\.0$', regex=True).any()
            if has_float:
                print(f"  WARNING: Some cell_ids have '.0' suffix (e.g., '1.0')")

    # Coordinate ranges
    print(f"\n[Coordinate Ranges]")
    for col in ['x_location', 'y_location', 'z_location']:
        if col in spots.columns:
            vals = spots[col].dropna()
            print(f"  {col}: min={vals.min():.2f}, max={vals.max():.2f}, "
                  f"range={vals.max()-vals.min():.2f} um")

    # Feature/Gene analysis
    gene_col = 'feature_name' if 'feature_name' in spots.columns else 'gene'
    if gene_col in spots.columns:
        genes = spots[gene_col]
        print(f"\n[Gene/Feature Analysis]")
        print(f"  Total unique features: {genes.nunique()}")

        # Control probes
        ctrl_mask = genes.str.contains('NegControl|BLANK|antisense', case=False, na=False)
        n_ctrl = ctrl_mask.sum()
        print(f"  Control probe transcripts: {n_ctrl:,} ({n_ctrl/len(spots)*100:.2f}%)")

        ctrl_types = genes[ctrl_mask].unique()
        if len(ctrl_types) <= 50:
            print(f"  Control probe types ({len(ctrl_types)}):")
            for ct in sorted(ctrl_types)[:20]:
                n = (genes == ct).sum()
                print(f"    {ct}: {n:,}")
            if len(ctrl_types) > 20:
                print(f"    ... and {len(ctrl_types)-20} more")

        # Top genes by count
        gene_counts = genes.value_counts()
        print(f"\n  Top 20 genes by transcript count:")
        for gene, cnt in gene_counts.head(20).items():
            print(f"    {gene:<30} {cnt:>10,} ({cnt/len(spots)*100:.1f}%)")

    # QV (Quality Value) analysis
    if 'qv' in spots.columns:
        qv = spots['qv']
        print(f"\n[Quality Value (QV)]")
        print(f"  min={qv.min():.1f}, max={qv.max():.1f}, "
              f"median={qv.median():.1f}, mean={qv.mean():.1f}")
        print(f"  QV >= 20: {(qv >= 20).sum():,} ({(qv>=20).sum()/len(spots)*100:.1f}%)")
        print(f"  QV < 20:  {(qv < 20).sum():,} ({(qv<20).sum()/len(spots)*100:.1f}%)")

    # Nucleus overlap
    if 'overlaps_nucleus' in spots.columns:
        ov = spots['overlaps_nucleus']
        print(f"\n[Nucleus Overlap]")
        print(f"  overlaps_nucleus=1 (nuclear): {(ov==1).sum():,} ({(ov==1).sum()/len(spots)*100:.1f}%)")
        print(f"  overlaps_nucleus=0 (other):   {(ov==0).sum():,} ({(ov==0).sum()/len(spots)*100:.1f}%)")

    return spots


# ============================================================
# Section 3: Cells Metadata
# ============================================================

def section_cells(input_path, **kwargs):
    """Inspect cells.csv: head, dtypes, spatial coordinate ranges."""
    _separator("3. CELLS METADATA")

    cells_path = os.path.join(input_path, "cells.csv")
    if not os.path.exists(cells_path):
        cells_gz = cells_path + ".gz"
        if os.path.exists(cells_gz):
            print(f"[Loading] {cells_gz}")
            cells = pd.read_csv(cells_gz, compression='gzip')
        else:
            print("  No cells.csv found!")
            return
    else:
        print(f"[Loading] {cells_path}")
        cells = pd.read_csv(cells_path)

    print(f"\n[Shape] {cells.shape[0]:,} cells x {cells.shape[1]} columns")

    print(f"\n[Columns & Dtypes]")
    for col in cells.columns:
        print(f"  {col:<30} dtype={str(cells[col].dtype):<12}")

    print(f"\n[Head - First 5 Rows]")
    print(cells.head().to_string(max_colwidth=25))

    # Cell ID format
    if 'cell_id' in cells.columns:
        cid = cells['cell_id']
        print(f"\n[Cell ID in cells.csv]")
        print(f"  dtype: {cid.dtype}")
        print(f"  range: {cid.min()} ~ {cid.max()}")
        print(f"  sample: {list(cid.head(5))}")

    # Spatial coordinate ranges
    print(f"\n[Spatial Coordinates]")
    for col in ['x_centroid', 'y_centroid']:
        if col in cells.columns:
            v = cells[col]
            print(f"  {col}: min={v.min():.2f}, max={v.max():.2f}, range={v.max()-v.min():.2f} um")

    # Cell area
    if 'cell_area' in cells.columns:
        a = cells['cell_area']
        print(f"\n[Cell Area]")
        print(f"  min={a.min():.1f}, max={a.max():.1f}, median={a.median():.1f}, mean={a.mean():.1f}")

    # Nucleus area
    if 'nucleus_area' in cells.columns:
        na = cells['nucleus_area']
        print(f"\n[Nucleus Area]")
        print(f"  min={na.min():.1f}, max={na.max():.1f}, median={na.median():.1f}, mean={na.mean():.1f}")

    # Transcript counts per cell
    if 'transcript_counts' in cells.columns:
        tc = cells['transcript_counts']
        print(f"\n[Transcript Counts per Cell]")
        print(f"  min={tc.min()}, max={tc.max()}, median={tc.median():.0f}, mean={tc.mean():.1f}")
        print(f"  cells with 0 transcripts: {(tc==0).sum()}")

    return cells


# ============================================================
# Section 4: Gene Panel
# ============================================================

def section_gene_panel(input_path, **kwargs):
    """Inspect gene_panel.json: number of genes, control probes, descriptors."""
    _separator("4. GENE PANEL")

    json_path = os.path.join(input_path, "gene_panel.json")
    if not os.path.exists(json_path):
        print("  gene_panel.json not found!")
        return

    with open(json_path) as f:
        data = json.load(f)

    targets = data.get('payload', {}).get('targets', [])
    print(f"\n[Panel Info]")
    print(f"  Total targets (probes): {len(targets)}")

    # Categorize
    genes = []
    controls = []
    descriptors = {}

    for t in targets:
        name = t['type']['data']['name']
        desc = t['type'].get('descriptor', 'unknown')
        descriptors[desc] = descriptors.get(desc, 0) + 1

        if any(kw in name.lower() for kw in ['negcontrol', 'blank', 'antisense']):
            controls.append(name)
        else:
            genes.append(name)

    print(f"  Biological genes: {len(genes)}")
    print(f"  Control probes: {len(controls)}")

    print(f"\n[Probe Descriptors]")
    for desc, cnt in sorted(descriptors.items(), key=lambda x: -x[1]):
        print(f"  {desc}: {cnt}")

    if controls:
        print(f"\n[Control Probes] ({len(controls)} total)")
        for c in sorted(controls)[:30]:
            print(f"  {c}")
        if len(controls) > 30:
            print(f"  ... and {len(controls) - 30} more")

    # Addon panel detection
    addon_genes = [t for t in targets if t['type'].get('descriptor', '') == 'addOn']
    if addon_genes:
        print(f"\n[Add-on Panel Genes] ({len(addon_genes)})")
        for g in addon_genes[:20]:
            print(f"  {g['type']['data']['name']}")

    return genes, controls


# ============================================================
# Section 5: QC - Control vs Raw Distributions (교수님 질문 ③)
# ============================================================

def section_qc_control(input_path, output_dir, sample_name, sample_tag, skip_plots=False, **kwargs):
    """Compare total_counts_control vs total_counts_raw distributions.

    Loads Step 0 h5ad (which has pre-computed QC columns) and visualizes
    the distributions that the professor asked about.
    """
    _separator("5. QC: CONTROL vs RAW DISTRIBUTIONS")

    if not HAS_SCANPY:
        print("  scanpy not available; skipping h5ad-based QC analysis.")
        return

    # Try to load Step 0 output
    step0_path = os.path.join(output_dir, sample_name, "step0_formatting", f"{sample_tag}.h5ad")
    if not os.path.exists(step0_path):
        print(f"  Step 0 output not found: {step0_path}")
        print("  Run the pipeline first, or provide the h5ad path.")
        # Fallback: compute from raw data
        return _qc_from_raw(input_path, output_dir, sample_name, sample_tag, skip_plots)

    print(f"[Loading] {step0_path}")
    adata = sc.read_h5ad(step0_path)
    print(f"  AnnData shape: {adata.shape}")
    print(f"  obs columns: {list(adata.obs.columns)}")

    # Check for QC columns
    has_raw = 'total_counts_raw' in adata.obs.columns
    has_ctrl = 'total_counts_control' in adata.obs.columns
    has_pct = 'pct_counts_control' in adata.obs.columns

    if has_raw and has_ctrl:
        raw = adata.obs['total_counts_raw']
        ctrl = adata.obs['total_counts_control']
        pct = adata.obs['pct_counts_control'] if has_pct else None

        print(f"\n[total_counts_raw]  (total reads per cell, before control removal)")
        print(f"  min={raw.min():.0f}, max={raw.max():.0f}, median={raw.median():.0f}, mean={raw.mean():.1f}")
        print(f"  cells with 0: {(raw==0).sum()}")

        print(f"\n[total_counts_control]  (control probe reads per cell)")
        print(f"  min={ctrl.min():.0f}, max={ctrl.max():.0f}, median={ctrl.median():.0f}, mean={ctrl.mean():.2f}")
        print(f"  cells with 0 control reads: {(ctrl==0).sum()} ({(ctrl==0).sum()/len(ctrl)*100:.1f}%)")
        print(f"  cells with >=1 control read: {(ctrl>=1).sum()} ({(ctrl>=1).sum()/len(ctrl)*100:.1f}%)")
        print(f"  cells with >=5 control reads: {(ctrl>=5).sum()} ({(ctrl>=5).sum()/len(ctrl)*100:.1f}%)")

        global_ctrl_pct = ctrl.sum() / raw.sum() * 100
        print(f"\n[Global Control %]")
        print(f"  Total control reads: {ctrl.sum():,.0f}")
        print(f"  Total raw reads: {raw.sum():,.0f}")
        print(f"  Global control %: {global_ctrl_pct:.3f}%")

        if pct is not None:
            print(f"\n[pct_counts_control per cell]")
            print(f"  min={pct.min():.3f}%, max={pct.max():.1f}%, median={pct.median():.3f}%, mean={pct.mean():.3f}%")

        # Estimated FDR (per 10x method)
        n_control_probes = ctrl.sum()  # total detected NCP transcripts
        # Count control probe types from gene panel
        gene_panel_path = os.path.join(input_path, 'gene_panel.json')
        n_ncp_types = 20  # default
        n_target_genes = adata.n_vars
        if os.path.exists(gene_panel_path):
            with open(gene_panel_path) as f:
                gp = json.load(f)
            targets = gp.get('payload', {}).get('targets', [])
            ncp_names = [t['type']['data']['name'] for t in targets
                         if 'negcontrol' in t['type']['data']['name'].lower()]
            if ncp_names:
                n_ncp_types = len(ncp_names)

        total_detected = raw.sum()
        ncp_per_probe = n_control_probes / n_ncp_types if n_ncp_types > 0 else 0
        fdr = (ncp_per_probe * n_target_genes / total_detected * 100) if total_detected > 0 else 0

        print(f"\n[Estimated FDR (10x method)]")
        print(f"  NCP probe types: {n_ncp_types}")
        print(f"  Target genes: {n_target_genes}")
        print(f"  NCP detections per probe: {ncp_per_probe:.1f}")
        print(f"  Estimated FDR: {fdr:.3f}%")

        # --- PLOTS ---
        if not skip_plots and HAS_PLOTS:
            plot_dir = os.path.join(output_dir, sample_name, "eda_plots")
            os.makedirs(plot_dir, exist_ok=True)
            _plot_qc_distributions(raw, ctrl, pct, plot_dir, sample_tag)
    else:
        print("  QC columns (total_counts_raw/total_counts_control) not found in adata.obs")
        print(f"  Available columns: {list(adata.obs.columns)}")

    return adata


def _qc_from_raw(input_path, output_dir, sample_name, sample_tag, skip_plots):
    """Fallback: compute QC from raw cell_feature_matrix when Step 0 not available."""
    print("\n[Fallback] Computing control stats from raw cell_feature_matrix...")

    cfm_path = os.path.join(input_path, 'cell_feature_matrix')
    features_path = os.path.join(cfm_path, 'features.tsv')
    if not os.path.exists(features_path):
        features_path += '.gz'
    if not os.path.exists(features_path):
        print("  Cannot compute: features.tsv not found")
        return

    features = pd.read_csv(features_path, header=None, sep='\t')
    if features.shape[1] >= 2:
        gene_ids = features.iloc[:, 0]
        ctrl_mask = (
            gene_ids.str.contains('NegControlProbe_', case=False) |
            gene_ids.str.contains('NegControlCodeword_', case=False) |
            gene_ids.str.contains('antisense_', case=False) |
            gene_ids.str.contains('BLANK', case=False)
        )
        print(f"  Total features: {len(gene_ids)}")
        print(f"  Control probes: {ctrl_mask.sum()}")
        print(f"  Biological genes: {(~ctrl_mask).sum()}")


def _plot_qc_distributions(raw, ctrl, pct, plot_dir, sample_tag):
    """Generate QC distribution plots for professor's review."""
    sns.set_style("white")

    # Plot 1: Histogram comparison - control vs raw (log scale)
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # 1a: Raw counts distribution
    axes[0].hist(raw.values, bins=100, color='steelblue', alpha=0.8, edgecolor='none')
    axes[0].set_xlabel('Total counts (raw)')
    axes[0].set_ylabel('Number of cells')
    axes[0].set_title('Raw Counts per Cell')
    axes[0].axvline(raw.median(), color='red', linestyle='--', label=f'median={raw.median():.0f}')
    axes[0].legend()

    # 1b: Control counts distribution
    axes[1].hist(ctrl.values, bins=np.arange(-0.5, ctrl.max()+1.5, 1),
                 color='salmon', alpha=0.8, edgecolor='none')
    axes[1].set_xlabel('Total counts (control probes)')
    axes[1].set_ylabel('Number of cells')
    axes[1].set_title('Control Probe Counts per Cell')
    axes[1].axvline(ctrl.median(), color='red', linestyle='--', label=f'median={ctrl.median():.0f}')
    axes[1].legend()
    # Annotate zero-count fraction
    zero_pct = (ctrl == 0).sum() / len(ctrl) * 100
    axes[1].text(0.95, 0.95, f'{zero_pct:.1f}% = 0',
                 transform=axes[1].transAxes, ha='right', va='top',
                 fontsize=11, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    # 1c: Overlay on log scale
    axes[2].hist(raw.values, bins=100, alpha=0.6, label='Raw', color='steelblue')
    axes[2].hist(ctrl.values, bins=100, alpha=0.6, label='Control', color='salmon')
    axes[2].set_yscale('log')
    axes[2].set_xlabel('Counts per cell')
    axes[2].set_ylabel('Number of cells (log)')
    axes[2].set_title('Raw vs Control (log scale)')
    axes[2].legend()

    fig.suptitle(f'{sample_tag} - QC: Raw vs Control Counts Distribution', fontsize=14)
    fig.tight_layout()
    path = os.path.join(plot_dir, f"{sample_tag}_eda_qc_raw_vs_control_hist.png")
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Plot saved: {path}")

    # Plot 2: Violin plot comparison
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    data_violin = pd.DataFrame({
        'Raw Counts': raw.values,
        'Control Counts': ctrl.values,
    })
    # 2a: Violin - raw
    axes[0].violinplot([raw.values], positions=[0], showmedians=True)
    axes[0].set_xticks([0])
    axes[0].set_xticklabels(['Raw'])
    axes[0].set_ylabel('Counts per cell')
    axes[0].set_title('Raw Counts Distribution')

    # 2b: Violin - control (zoomed)
    axes[1].violinplot([ctrl.values], positions=[0], showmedians=True)
    axes[1].set_xticks([0])
    axes[1].set_xticklabels(['Control'])
    axes[1].set_ylabel('Counts per cell')
    axes[1].set_title('Control Counts Distribution')

    fig.suptitle(f'{sample_tag} - QC: Violin Plots', fontsize=14)
    fig.tight_layout()
    path = os.path.join(plot_dir, f"{sample_tag}_eda_qc_violin.png")
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Plot saved: {path}")

    # Plot 3: pct_counts_control per cell
    if pct is not None:
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        axes[0].hist(pct.values, bins=100, color='darkorange', alpha=0.8, edgecolor='none')
        axes[0].set_xlabel('% Control Reads per Cell')
        axes[0].set_ylabel('Number of Cells')
        axes[0].set_title('Control Read Percentage Distribution')
        axes[0].axvline(pct.median(), color='red', linestyle='--',
                        label=f'median={pct.median():.3f}%')
        axes[0].legend()

        # Zoomed: cells with >1% control
        high_ctrl = pct[pct > 1.0]
        if len(high_ctrl) > 0:
            axes[1].hist(high_ctrl.values, bins=50, color='red', alpha=0.7)
            axes[1].set_xlabel('% Control Reads per Cell')
            axes[1].set_ylabel('Number of Cells')
            axes[1].set_title(f'Cells with >1% Control ({len(high_ctrl)} cells)')
        else:
            axes[1].text(0.5, 0.5, 'No cells with >1%\ncontrol reads',
                         transform=axes[1].transAxes, ha='center', va='center', fontsize=14)
            axes[1].set_title('Cells with >1% Control')

        fig.suptitle(f'{sample_tag} - QC: Control Read Percentage', fontsize=14)
        fig.tight_layout()
        path = os.path.join(plot_dir, f"{sample_tag}_eda_qc_pct_control.png")
        fig.savefig(path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f"  Plot saved: {path}")

    # Plot 4: Scatter - raw vs control per cell
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(raw.values, ctrl.values, s=0.5, alpha=0.3, rasterized=True)
    ax.set_xlabel('Total Raw Counts')
    ax.set_ylabel('Total Control Counts')
    ax.set_title(f'{sample_tag} - Raw vs Control per Cell')
    ax.axhline(y=ctrl.median(), color='red', linestyle='--', alpha=0.5, label=f'ctrl median={ctrl.median():.0f}')
    ax.legend()
    fig.tight_layout()
    path = os.path.join(plot_dir, f"{sample_tag}_eda_qc_raw_vs_control_scatter.png")
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Plot saved: {path}")


# ============================================================
# Section 6: DAPI Image Metadata & Coordinate Alignment (교수님 질문 ①)
# ============================================================

def section_dapi_alignment(input_path, **kwargs):
    """Inspect DAPI image dimensions and compare with transcript coordinate ranges."""
    _separator("6. DAPI-TRANSCRIPT COORDINATE ALIGNMENT")

    if not HAS_TIFF:
        print("  tifffile not available; skipping DAPI analysis.")
        return

    # Pixel size from Xenium (standard)
    pixel_size_um = 0.2125  # um per pixel

    # Load DAPI metadata (without loading full image into memory)
    dapi_files = [
        ("morphology_focus.ome.tif", "DAPI Focus"),
        ("morphology_mip.ome.tif", "DAPI MIP (Max Intensity Projection)"),
    ]

    for fname, desc in dapi_files:
        fpath = os.path.join(input_path, fname)
        if not os.path.exists(fpath):
            continue

        print(f"\n[{desc}] {fname}")
        try:
            with tf.TiffFile(fpath) as tif:
                # Read metadata without loading image
                page0 = tif.pages[0]
                shape = page0.shape
                dtype = page0.dtype
                print(f"  Shape: {shape}")
                print(f"  Dtype: {dtype}")
                print(f"  Pages/frames: {len(tif.pages)}")

                # Compute physical dimensions
                if len(shape) >= 2:
                    h, w = shape[-2], shape[-1]
                    phys_h = h * pixel_size_um
                    phys_w = w * pixel_size_um
                    print(f"  Pixel dimensions: {h:,} x {w:,}")
                    print(f"  Physical size: {phys_h:.1f} x {phys_w:.1f} um "
                          f"({phys_h/1000:.2f} x {phys_w/1000:.2f} mm)")
                    print(f"  Pixel size: {pixel_size_um} um/px")

                # OME metadata
                if tif.ome_metadata:
                    print(f"  OME metadata: present ({len(tif.ome_metadata)} chars)")
        except Exception as e:
            print(f"  Error reading: {e}")

    # Compare with transcript coordinate ranges
    print(f"\n[Coordinate Alignment Check]")
    parquet_path = os.path.join(input_path, "transcripts.parquet")
    if os.path.exists(parquet_path):
        # Read only coordinate columns for speed
        spots = pd.read_parquet(parquet_path, columns=['x_location', 'y_location'])
        x_range = (spots['x_location'].min(), spots['x_location'].max())
        y_range = (spots['y_location'].min(), spots['y_location'].max())

        print(f"  Transcript X range: {x_range[0]:.2f} ~ {x_range[1]:.2f} um")
        print(f"  Transcript Y range: {y_range[0]:.2f} ~ {y_range[1]:.2f} um")

        # Convert to pixel coords (like pipeline Step 3 does)
        um_per_pixel_inv = 4.70588  # px/um (from config)
        x_px_range = (x_range[0] * um_per_pixel_inv, x_range[1] * um_per_pixel_inv)
        y_px_range = (y_range[0] * um_per_pixel_inv, y_range[1] * um_per_pixel_inv)
        print(f"  Transcript X (pixels): {x_px_range[0]:.0f} ~ {x_px_range[1]:.0f}")
        print(f"  Transcript Y (pixels): {y_px_range[0]:.0f} ~ {y_px_range[1]:.0f}")

        # Compare with DAPI dimensions
        focus_path = os.path.join(input_path, "morphology_focus.ome.tif")
        if os.path.exists(focus_path):
            with tf.TiffFile(focus_path) as tif:
                page0 = tif.pages[0]
                img_h, img_w = page0.shape[-2], page0.shape[-1]

            print(f"\n  DAPI image: {img_h} x {img_w} pixels")
            print(f"  DAPI physical: {img_h*pixel_size_um:.1f} x {img_w*pixel_size_um:.1f} um")

            # Check overlap
            x_in_range = x_px_range[1] <= img_w
            y_in_range = y_px_range[1] <= img_h
            print(f"\n  Transcripts fit within DAPI X: {'YES' if x_in_range else 'NO'} "
                  f"(max_tx_px={x_px_range[1]:.0f}, img_w={img_w})")
            print(f"  Transcripts fit within DAPI Y: {'YES' if y_in_range else 'NO'} "
                  f"(max_tx_px={y_px_range[1]:.0f}, img_h={img_h})")

            if not x_in_range or not y_in_range:
                print("  WARNING: Some transcripts fall outside DAPI image bounds!")
                print("  This may indicate coordinate system mismatch or different pixel scale.")
            else:
                coverage_x = (x_px_range[1] - x_px_range[0]) / img_w * 100
                coverage_y = (y_px_range[1] - y_px_range[0]) / img_h * 100
                print(f"  Coverage: X={coverage_x:.1f}%, Y={coverage_y:.1f}%")

    print(f"\n[Alignment Summary]")
    print(f"  Xenium uses the SAME coordinate system for DAPI images and transcripts.")
    print(f"  Pixel size: {pixel_size_um} um/pixel (standard Xenium).")
    print(f"  Known limitation: DAPI stains nuclei only. Many transcripts are")
    print(f"  cytoplasmic and fall OUTSIDE nuclear boundaries - this is biological,")
    print(f"  not a coordinate alignment error.")
    print(f"  Our pipeline (Step 3) uses Cellpose + expand_labels to capture")
    print(f"  cytoplasmic transcripts via nuclear expansion.")


# ============================================================
# Section 7: 10X Metrics Summary
# ============================================================

def section_10x_metrics(input_path, **kwargs):
    """Read and display the 10X Xenium metrics_summary.csv."""
    _separator("7. 10X XENIUM METRICS SUMMARY")

    metrics_path = os.path.join(input_path, "metrics_summary.csv")
    if not os.path.exists(metrics_path):
        print("  metrics_summary.csv not found!")
        return

    metrics = pd.read_csv(metrics_path)
    print(f"\n[Metrics from Xenium Instrument]")
    for _, row in metrics.iterrows():
        for col in metrics.columns:
            val = row[col]
            print(f"  {col}: {val}")
    print()


# ============================================================
# Section 8: Pipeline h5ad Inspection
# ============================================================

def section_h5ad_inspect(output_dir, sample_name, sample_tag, **kwargs):
    """Inspect all pipeline h5ad outputs: shape, layers, obs columns."""
    _separator("8. PIPELINE H5AD OUTPUT INSPECTION")

    if not HAS_SCANPY:
        print("  scanpy not available.")
        return

    sample_dir = os.path.join(output_dir, sample_name)
    if not os.path.isdir(sample_dir):
        print(f"  Pipeline output not found: {sample_dir}")
        return

    step_dirs = [
        ("step0_formatting", f"{sample_tag}.h5ad"),
        ("step1_exploration", f"{sample_tag}_step1_exploration.h5ad"),
        ("step2_segmentation_free", f"{sample_tag}_step2_points2regions.h5ad"),
        ("step3_resegmentation", f"{sample_tag}_step3_resegmented.h5ad"),
    ]

    for step_name, h5ad_file in step_dirs:
        h5ad_path = os.path.join(sample_dir, step_name, h5ad_file)
        if not os.path.exists(h5ad_path):
            print(f"\n  [{step_name}] NOT FOUND: {h5ad_file}")
            continue

        try:
            adata = sc.read_h5ad(h5ad_path)
            file_size = os.path.getsize(h5ad_path)

            print(f"\n  [{step_name}] {h5ad_file} ({_sizeof_fmt(file_size)})")
            print(f"    Shape: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
            print(f"    X dtype: {adata.X.dtype}, sparse: {hasattr(adata.X, 'nnz')}")

            if adata.layers:
                print(f"    Layers: {list(adata.layers.keys())}")
            if list(adata.obsm.keys()):
                print(f"    obsm: {list(adata.obsm.keys())}")
            if list(adata.uns.keys()):
                print(f"    uns keys: {list(adata.uns.keys())[:10]}")

            print(f"    obs columns ({len(adata.obs.columns)}):")
            for col in adata.obs.columns:
                dtype = adata.obs[col].dtype
                sample_val = adata.obs[col].dropna().iloc[0] if len(adata.obs[col].dropna()) > 0 else "N/A"
                if isinstance(sample_val, (float, np.floating)):
                    sample_val = f"{sample_val:.4f}"
                print(f"      {col:<30} dtype={str(dtype):<15} sample={sample_val}")

            print(f"    var columns ({len(adata.var.columns)}):")
            for col in adata.var.columns:
                print(f"      {col}")

            # obs_names (cell IDs) format
            print(f"    obs_names (cell IDs):")
            print(f"      dtype: {adata.obs_names.dtype}")
            print(f"      first 5: {list(adata.obs_names[:5])}")
            print(f"      last 5: {list(adata.obs_names[-5:])}")

            del adata
        except Exception as e:
            print(f"\n  [{step_name}] ERROR reading {h5ad_file}: {e}")


# ============================================================
# Section 9: RAM Usage Estimation (교수님 질문 ④)
# ============================================================

def section_ram_estimation(input_path, **kwargs):
    """Estimate peak RAM usage for each pipeline step."""
    _separator("9. RAM USAGE ESTIMATION")

    # Get actual data sizes
    parquet_path = os.path.join(input_path, "transcripts.parquet")
    csv_path = os.path.join(input_path, "transcripts.csv")

    if os.path.exists(csv_path):
        csv_size = os.path.getsize(csv_path)
    else:
        csv_size = 0

    if os.path.exists(parquet_path):
        parquet_size = os.path.getsize(parquet_path)
    else:
        parquet_size = 0

    # Count transcripts for estimation
    if os.path.exists(parquet_path):
        n_transcripts = len(pd.read_parquet(parquet_path, columns=['x_location']))
    else:
        n_transcripts = 16_300_000  # approximate

    # Cells
    cells_path = os.path.join(input_path, "cells.csv")
    if os.path.exists(cells_path):
        n_cells = len(pd.read_csv(cells_path, usecols=[0]))
    else:
        n_cells = 120_000  # approximate

    # DAPI image size
    focus_path = os.path.join(input_path, "morphology_focus.ome.tif")
    if os.path.exists(focus_path) and HAS_TIFF:
        with tf.TiffFile(focus_path) as tif:
            page0 = tif.pages[0]
            dapi_h, dapi_w = page0.shape[-2], page0.shape[-1]
            dapi_bytes_per_px = page0.dtype.itemsize
    else:
        dapi_h, dapi_w = 40000, 37000
        dapi_bytes_per_px = 2  # uint16

    dapi_size_gb = dapi_h * dapi_w * dapi_bytes_per_px / 1e9

    print(f"\n[Current Dataset Statistics]")
    print(f"  Transcripts: {n_transcripts:,}")
    print(f"  Cells (10X): {n_cells:,}")
    print(f"  DAPI image: {dapi_h:,} x {dapi_w:,} ({dapi_size_gb:.1f} GB)")
    print(f"  transcripts.csv: {_sizeof_fmt(csv_size)}")
    print(f"  transcripts.parquet: {_sizeof_fmt(parquet_size)}")

    # Current server
    if HAS_PSUTIL:
        vm = psutil.virtual_memory()
        print(f"\n[Server Memory]")
        print(f"  Total RAM: {vm.total/1e9:.0f} GB")
        print(f"  Available: {vm.available/1e9:.0f} GB")
        print(f"  Used: {vm.used/1e9:.0f} GB ({vm.percent}%)")

    print(f"\n{'Step':<35} {'Peak RAM (est.)':>15} {'Notes'}")
    print("-" * 80)

    # Step 0: Formatting
    # Loads: matrix.mtx (dense), cells.csv, transcripts DataFrame
    step0_ram = (n_cells * 313 * 8 / 1e9  # dense matrix (float64)
                 + csv_size / 1e9 * 3      # DataFrame overhead ~3x CSV
                 + 0.5)                     # misc
    print(f"  {'Step 0: Formatting':<33} {step0_ram:>12.1f} GB  "
          f"dense matrix + transcripts DataFrame")

    # Step 1: Exploration
    step1_ram = (n_cells * 313 * 8 / 1e9 * 2  # adata + copy for operations
                 + 1.0)                          # scanpy intermediates
    print(f"  {'Step 1: Exploration':<33} {step1_ram:>12.1f} GB  "
          f"adata + scanpy PCA/UMAP/neighbors")

    # Step 2: Points2Regions
    step2_ram = (csv_size / 1e9 * 3   # transcripts DataFrame
                 + 2.0                  # P2R binning grids
                 + 1.0)                 # adata
    print(f"  {'Step 2: Segmentation-Free (P2R)':<33} {step2_ram:>12.1f} GB  "
          f"transcript DF + P2R grids")

    # Step 3: Resegmentation
    # Peak: DAPI image + Cellpose masks (int32) + transcripts DF
    masks_size_gb = dapi_h * dapi_w * 4 / 1e9  # int32 masks
    step3_ram_gpu = dapi_size_gb  # DAPI loaded to memory
    step3_ram_cpu = (dapi_size_gb                       # DAPI image
                     + masks_size_gb * 2                # nuclei + expanded masks
                     + csv_size / 1e9 * 3               # transcripts DF
                     + n_cells * 313 * 8 / 1e9          # output AnnData
                     + masks_size_gb * 0.5)             # expand_labels working copy
    print(f"  {'Step 3: Resegmentation (Cellpose)':<33} {step3_ram_cpu:>12.1f} GB  "
          f"DAPI({dapi_size_gb:.1f}G) + masks({masks_size_gb:.1f}G x2) + transcripts")
    print(f"  {'  + GPU VRAM (per tile)':<33} {'~4-8':>12} GB  "
          f"Cellpose flow dynamics (tiled to fit 24GB A5000)")

    # Step 4: Techniques Comparison
    step4_ram = 3.0  # lightweight comparisons
    print(f"  {'Step 4: Techniques Comparison':<33} {step4_ram:>12.1f} GB  "
          f"adata comparisons + scRNA ref subset")

    # Step 5: Optimal Expansion
    step5_ram = (csv_size / 1e9 * 3 + 2.0)
    print(f"  {'Step 5: Optimal Expansion':<33} {step5_ram:>12.1f} GB  "
          f"transcripts + KDTree + correlation")

    # Step 6: Benchmark (PEAK)
    # Loads: multiple adata objects + Baysor output + raw layers
    step6_ram = (n_cells * 313 * 8 / 1e9 * 4  # 4 methods x adata
                 + csv_size / 1e9 * 3           # transcripts for Baysor
                 + 3.0)                          # preprocessing intermediates
    print(f"  {'Step 6: Benchmark':<33} {step6_ram:>12.1f} GB  "
          f"4 adata objects + Baysor transcripts")

    # Step 7: Simulation
    step7_ram = 4.0
    print(f"  {'Step 7: Simulation':<33} {step7_ram:>12.1f} GB  "
          f"CellxGene download + simulation grid")

    # Overall peak
    peak = max(step0_ram, step1_ram, step2_ram, step3_ram_cpu, step4_ram, step5_ram, step6_ram, step7_ram)
    print("-" * 80)
    print(f"  {'ESTIMATED PEAK (this dataset)':<33} {peak:>12.1f} GB")

    # 5K Panel estimation
    print(f"\n[5K Panel Estimation]")
    # 5K panel: ~5000 genes instead of ~313
    scale_genes = 5000 / 313
    # More transcripts: roughly 5-10x more (Bilous 2025 notes more transcripts)
    scale_tx = 7  # conservative estimate

    n_tx_5k = int(n_transcripts * scale_tx)
    csv_5k = csv_size * scale_tx

    print(f"  Gene count scaling: {313} -> 5000 ({scale_genes:.0f}x)")
    print(f"  Transcript count estimate: {n_transcripts:,} -> {n_tx_5k:,} ({scale_tx}x)")

    step0_5k = n_cells * 5000 * 8 / 1e9 + csv_5k / 1e9 * 3 + 0.5
    step3_5k = (dapi_size_gb + masks_size_gb * 2 + csv_5k / 1e9 * 3
                + n_cells * 5000 * 8 / 1e9 + masks_size_gb * 0.5)
    step6_5k = n_cells * 5000 * 8 / 1e9 * 4 + csv_5k / 1e9 * 3 + 5.0

    print(f"\n{'Step':<35} {'313-gene':>12} {'5K panel':>12}")
    print("-" * 65)
    print(f"  {'Step 0: Formatting':<33} {step0_ram:>10.1f} GB {step0_5k:>10.1f} GB")
    print(f"  {'Step 3: Resegmentation':<33} {step3_ram_cpu:>10.1f} GB {step3_5k:>10.1f} GB")
    print(f"  {'Step 6: Benchmark':<33} {step6_ram:>10.1f} GB {step6_5k:>10.1f} GB")
    peak_5k = max(step0_5k, step3_5k, step6_5k)
    print("-" * 65)
    print(f"  {'ESTIMATED PEAK (5K panel)':<33} {peak:>10.1f} GB {peak_5k:>10.1f} GB")

    if HAS_PSUTIL:
        vm = psutil.virtual_memory()
        avail_gb = vm.available / 1e9
        print(f"\n[Verdict]")
        print(f"  Server available RAM: {avail_gb:.0f} GB")
        print(f"  Current dataset peak: {peak:.1f} GB -> {'OK' if peak < avail_gb * 0.8 else 'TIGHT'}")
        print(f"  5K panel peak: {peak_5k:.1f} GB -> {'OK' if peak_5k < avail_gb * 0.8 else 'TIGHT'}")
        if peak_5k < avail_gb * 0.5:
            print(f"  Conclusion: Server has ample headroom for 5K panel ({avail_gb:.0f}GB >> {peak_5k:.0f}GB)")

    # Baysor-specific notes
    print(f"\n[Baysor RAM Notes]")
    print(f"  Baysor runs as an external Julia process.")
    print(f"  For full-tissue (16M transcripts): expect 30-50 GB RAM, 24-72h runtime.")
    print(f"  For crop mode (2000x2000 um): ~5-10 GB RAM, 5-10h runtime.")
    print(f"  Config currently uses crop mode (baysor.crop.enabled=true).")
    print(f"  For 5K panel full-tissue ({n_tx_5k/1e6:.0f}M tx): could need 200+ GB RAM.")
    print(f"  Recommendation: Always use crop mode or tiling for 5K panel Baysor runs.")


# ============================================================
# Section 10: Experiment Metadata
# ============================================================

def section_experiment(input_path, **kwargs):
    """Read experiment.xenium metadata."""
    _separator("10. EXPERIMENT METADATA")

    exp_path = os.path.join(input_path, "experiment.xenium")
    if not os.path.exists(exp_path):
        print("  experiment.xenium not found!")
        return

    with open(exp_path) as f:
        try:
            data = json.load(f)
            print(json.dumps(data, indent=2))
        except json.JSONDecodeError:
            f.seek(0)
            print(f.read())


# ============================================================
# Main
# ============================================================

SECTIONS = {
    'files': section_files,
    'transcripts': section_transcripts,
    'cells': section_cells,
    'gene_panel': section_gene_panel,
    'qc': section_qc_control,
    'dapi': section_dapi_alignment,
    'metrics': section_10x_metrics,
    'h5ad': section_h5ad_inspect,
    'ram': section_ram_estimation,
    'experiment': section_experiment,
}


def main():
    parser = argparse.ArgumentParser(description="Xenium Data EDA")
    parser.add_argument("--config", default=os.path.join(os.path.dirname(__file__), "config.yaml"),
                        help="Path to pipeline config.yaml")
    parser.add_argument("--skip-plots", action="store_true",
                        help="Skip generating plot files")
    parser.add_argument("--section", default="all",
                        help="Comma-separated sections to run (or 'all'). "
                             f"Available: {','.join(SECTIONS.keys())}")
    parser.add_argument("--output-report", default=None,
                        help="Save text output to file")
    args = parser.parse_args()

    # Load config
    if not os.path.exists(args.config):
        print(f"Config not found: {args.config}")
        return

    with open(args.config) as f:
        config = yaml.safe_load(f)

    input_path = config['input_path']
    output_dir = config.get('output_dir', 'xenium-output')
    sample_tag = config.get('sample_tag', 'human_alzheimers')

    # Derive sample name (same logic as pipeline_main.py)
    sample_name = os.path.basename(input_path.rstrip(os.sep))

    print("=" * 70)
    print("  XENIUM DATA - EXPLORATORY DATA ANALYSIS (EDA)")
    print("=" * 70)
    print(f"  Input: {input_path}")
    print(f"  Output: {output_dir}")
    print(f"  Sample: {sample_tag}")
    print(f"  Date: {time.strftime('%Y-%m-%d %H:%M:%S')}")

    # Redirect output if requested
    if args.output_report:
        import io
        import contextlib
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            _run_sections(args, input_path, output_dir, sample_name, sample_tag)
        report = buf.getvalue()
        print(report)
        with open(args.output_report, 'w') as f:
            f.write(report)
        print(f"\nReport saved to: {args.output_report}")
    else:
        _run_sections(args, input_path, output_dir, sample_name, sample_tag)


def _run_sections(args, input_path, output_dir, sample_name, sample_tag):
    """Run selected EDA sections."""
    if args.section == 'all':
        sections_to_run = list(SECTIONS.keys())
    else:
        sections_to_run = [s.strip() for s in args.section.split(',')]

    common_kwargs = dict(
        input_path=input_path,
        output_dir=output_dir,
        sample_name=sample_name,
        sample_tag=sample_tag,
        skip_plots=args.skip_plots,
    )

    for section_name in sections_to_run:
        if section_name not in SECTIONS:
            print(f"\n  WARNING: Unknown section '{section_name}'. "
                  f"Available: {', '.join(SECTIONS.keys())}")
            continue
        try:
            SECTIONS[section_name](**common_kwargs)
        except Exception as e:
            print(f"\n  ERROR in section '{section_name}': {e}")
            import traceback
            traceback.print_exc()

    _separator("EDA COMPLETE")


if __name__ == '__main__':
    main()
