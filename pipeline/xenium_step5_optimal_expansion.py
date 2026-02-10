# Step 5: Optimal Expansion (Ref: notebooks/4_optimal_expansion/4_1)
# Assigns unassigned transcripts to nearest annotated cell domain via KDTree,
# then computes correlation-based turnover to find optimal expansion radius.
#
# Flow:
#   5-1. Load original transcripts (Step 0)
#   5-2. Load annotated cells (Step 1/2/4)
#   5-3. Map domain assignments to reads
#   5-4. KDTree nearest-domain expansion
#   5-5. Distance threshold filter (optional)
#   5-6. Spatial map + distance histogram
#   5-7. Save expanded transcripts
#   5-8. Turnover analysis
#     5-8a. Nuclear vs background expression profiles
#     5-8b. Distance-bin correlation curves
#     5-8c. Turnover distance detection
#     5-8d. Optimal expansion = turnover − nuclei_size
#     5-8e. Summary barplot + CSVs

import os
import logging
import pandas as pd
import numpy as np
import scanpy as sc
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree, ConvexHull
import random as rd
import math

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


def dist_nuc(reads_ctdsub):
    """Median distance to nuclei edges via ConvexHull vertices."""
    allds = []
    for g, n in reads_ctdsub.groupby('cell_id'):
        try:
            if len(n) < 3:
                continue
            hull = ConvexHull(np.array(n.loc[:, ['x_location', 'y_location']]))
            if 'distance' in n.columns:
                allds.append(np.mean(n.iloc[hull.vertices]['distance']))
        except Exception:
            pass
    if len(allds) > 0:
        median_dist = np.median(allds)
    else:
        median_dist = np.nan
    return median_dist


def distance_calc(x1, y1, x2, y2):
    """Euclidean distance between two points."""
    return math.sqrt(((x1 - x2) ** 2) + ((y1 - y2) ** 2))


def hex_to_rgb(value):
    """Hex color string to RGB tuple."""
    value = value.lstrip('#')
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))


logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


# --- Correlation-based Turnover / Crossover Analysis ---

def calculate_turnover(reads_assigned, reads_not_assigned, output_dir, sample_tag,
                       min_reads_per_domain=5000, diff_threshold=0.1):
    """
    Correlation-based turnover/crossover for each cell-type x domain.

    Algorithm per (cell-type, domain) pair with enough reads:
      1. Build distance x gene crosstab from assigned reads.
      2. Get nuclear-expression profile (overlaps_nucleus == 1).
      3. Get background profile from unassigned reads.
      4. Correlate expression at each distance bin with nuclear and background.
      5. Turnover distance = first bin where (corr_nuc - corr_back) < diff_threshold.
      6. Compute nuclei_size and cell_size via dist_nuc.

    Returns (turnover_summary, per_celltype_df, optimal_expansion_value).
    """
    print("\n    [Turnover] Starting correlation-based turnover analysis ...")

    # Ensure overlaps_nucleus exists; use distance-based proxy if missing
    has_overlap_col = 'overlaps_nucleus' in reads_assigned.columns
    if not has_overlap_col:
        print("    [Turnover] 'overlaps_nucleus' column missing -- using distance-based proxy.")
        med_dist = reads_assigned['distance'].median()
        reads_assigned = reads_assigned.copy()
        reads_assigned['overlaps_nucleus'] = (reads_assigned['distance'] < med_dist).astype(int)
        print(f"    [Turnover] Proxy threshold (median distance): {med_dist:.2f}")

    # Background expression from unassigned reads
    reads_not_assigned_clean = reads_not_assigned.dropna(subset=['domain', 'feature_name'])
    if len(reads_not_assigned_clean) == 0:
        print("    [Turnover] No usable background reads -- aborting turnover calculation.")
        return None, None, np.nan

    # --- 5-8a. Nuclear vs background expression profiles ---
    background_express = pd.crosstab(
        reads_not_assigned_clean['domain'],
        reads_not_assigned_clean['feature_name']
    )

    unique_celltypes = reads_assigned['initial_annotation'].dropna().unique()
    unique_domains = reads_assigned['domain'].dropna().unique()

    turnover_summ = pd.DataFrame(
        columns=unique_celltypes,
        index=unique_domains,
        dtype=float
    )

    meand_celltype = []
    nuclimall = []
    cellmall = []
    ct_labels = []

    figures_dir = os.path.join(output_dir, "crossover_plots")
    os.makedirs(figures_dir, exist_ok=True)

    # --- Main loop: cell type x domain ---
    for celltype in unique_celltypes:
        if pd.isna(celltype):
            continue
        print(f"    [Turnover] Processing cell type: {celltype}")

        reads_ct = reads_assigned[reads_assigned['initial_annotation'] == celltype]

        domain_counts = reads_ct.groupby('domain').size()
        good_domains = domain_counts[domain_counts > min_reads_per_domain].index

        if len(good_domains) == 0:
            print(f"      - No domains with > {min_reads_per_domain} reads. Skipping.")
            continue

        tdistanceall = []
        nuclim = []
        cellm = []

        for dom in good_domains:
            reads_ctd = reads_ct[reads_ct['domain'] == dom]

            # 1. Distance x gene crosstab
            expression_distances = pd.crosstab(
                reads_ctd['distance'], reads_ctd['feature_name']
            )

            # 2. Nuclear expression profile
            reads_ctd_nucl = pd.crosstab(
                reads_ctd['overlaps_nucleus'], reads_ctd['feature_name']
            )
            if 1 not in reads_ctd_nucl.index:
                continue

            # 3. Background expression for this domain
            if dom not in background_express.index:
                continue
            bck_sub = background_express.copy()

            # Align columns across all three matrices
            common_genes = (
                set(reads_ctd_nucl.columns) &
                set(bck_sub.columns) &
                set(expression_distances.columns)
            )
            if len(common_genes) < 5:
                continue

            common_genes = sorted(common_genes)
            reads_ctd_nucl = reads_ctd_nucl[common_genes]
            bck_sub = bck_sub[common_genes]
            expression_distances = expression_distances[common_genes]

            nuc_vec = reads_ctd_nucl.loc[1, :].values.astype(float)
            bck_vec = bck_sub.loc[dom, :].values.astype(float)

            # --- 5-8b. Distance-bin correlation curves ---
            corr_nuc_list = []
            corr_bck_list = []
            for dist in expression_distances.index:
                dist_vec = expression_distances.loc[dist, :].values.astype(float)
                if np.std(dist_vec) == 0 or np.std(nuc_vec) == 0:
                    corr_nuc_list.append(np.nan)
                else:
                    corr_nuc_list.append(np.corrcoef(dist_vec, nuc_vec)[0, 1])
                if np.std(dist_vec) == 0 or np.std(bck_vec) == 0:
                    corr_bck_list.append(np.nan)
                else:
                    corr_bck_list.append(np.corrcoef(dist_vec, bck_vec)[0, 1])

            summary = pd.DataFrame({
                'corr_nuc': corr_nuc_list,
                'corr_back': corr_bck_list
            }, index=expression_distances.index)
            summary.index.name = 'distance'
            summary['corr_nuc'] = summary['corr_nuc'].astype(float)
            summary['corr_back'] = summary['corr_back'].astype(float)
            summary['diff'] = summary['corr_nuc'] - summary['corr_back']

            # Per-domain crossover plot
            try:
                fig, ax = plt.subplots(figsize=(6, 4))
                summary_plot = summary.reset_index()
                sns.lineplot(data=summary_plot, x='distance', y='corr_nuc', label='corr_nuc', ax=ax)
                sns.lineplot(data=summary_plot, x='distance', y='corr_back', label='corr_back', ax=ax)
                dom_str = str(dom).replace('/', '_')
                ct_str = str(celltype).replace('/', '_')
                ax.set_title(f"Domain {dom}, celltype {celltype}")
                ax.set_ylabel("Correlation")
                ax.set_xlabel("Distance (rounded)")
                fig_path = os.path.join(figures_dir, f"{dom_str}_{ct_str}.png")
                fig.savefig(fig_path, dpi=150, bbox_inches='tight')
                plt.close(fig)
            except Exception as e:
                logger.debug(f"      Plot failed for {dom}/{celltype}: {e}")
                plt.close('all')

            # --- 5-8c. Turnover distance detection ---
            below_threshold = summary.loc[summary['diff'] < diff_threshold, :]
            try:
                tdistance = np.nanmin(below_threshold.index)
                if np.isnan(tdistance):
                    tdistance = np.nanmax(summary.index)
            except (ValueError, TypeError):
                tdistance = np.nanmax(summary.index)

            turnover_summ.loc[dom, celltype] = tdistance
            tdistanceall.append(tdistance)

            # 6. nuclei_size and cell_size
            reads_ctdsub_nuc = reads_ctd[reads_ctd['overlaps_nucleus'] == 1]
            cellm.append(dist_nuc(reads_ctd))
            nuclim.append(dist_nuc(reads_ctdsub_nuc))

        ct_labels.append(celltype)
        meand_celltype.append(np.nanmean(tdistanceall) if len(tdistanceall) > 0 else np.nan)
        nuclimall.append(np.nanmean(nuclim) if len(nuclim) > 0 else np.nan)
        cellmall.append(np.nanmean(cellm) if len(cellm) > 0 else np.nan)

    # --- Per-celltype summary ---
    per_celltype = pd.DataFrame({
        'celltype': ct_labels,
        'turnover': meand_celltype,
        'nuclei_size': nuclimall,
        'cell_size': cellmall
    })

    # --- 5-8d. Optimal expansion = turnover − nuclei_size ---
    optimal_expansion = np.nanmean(meand_celltype) - np.nanmean(nuclimall)

    print(f"\n    [Turnover] Mean turnover distance : {np.nanmean(meand_celltype):.3f}")
    print(f"    [Turnover] Mean nuclei size       : {np.nanmean(nuclimall):.3f}")
    print(f"    [Turnover] >>> Optimal expansion   : {optimal_expansion:.3f}")

    # --- 5-8e. Summary barplot + CSVs ---
    try:
        per_ct_sorted = per_celltype.dropna(subset=['turnover']).sort_values('turnover')
        if len(per_ct_sorted) > 0:
            fig, ax = plt.subplots(figsize=(10, max(4, len(per_ct_sorted) * 0.5)))
            sns.barplot(
                data=per_ct_sorted, y='celltype', x='turnover',
                palette='tab20', alpha=0.85, ax=ax
            )
            sns.scatterplot(
                data=per_ct_sorted, y='celltype', x='nuclei_size',
                color='#D83066', edgecolor=None, alpha=0.7, s=80, label='nuclei_size', ax=ax
            )
            sns.scatterplot(
                data=per_ct_sorted, y='celltype', x='cell_size',
                color='black', edgecolor='gray', s=80, label='cell_size', ax=ax
            )
            ax.set_xlabel("Distance")
            ax.set_title("Turnover per cell type")
            ax.legend(loc='lower right')
            summary_plot_path = os.path.join(output_dir, f"{sample_tag}_step5_turnover_barplot.png")
            fig.savefig(summary_plot_path, dpi=300, bbox_inches='tight')
            plt.close(fig)
            print(f"    [Turnover] Saved barplot: {summary_plot_path}")
    except Exception as e:
        logger.warning(f"    [Turnover] Summary barplot failed: {e}")
        plt.close('all')

    # --- Save CSVs ---
    csv_summ = os.path.join(output_dir, f"{sample_tag}_step5_turnover_summary.csv")
    turnover_summ.to_csv(csv_summ)
    print(f"    [Turnover] Saved turnover summary matrix: {csv_summ}")

    csv_ct = os.path.join(output_dir, f"{sample_tag}_step5_turnover_per_celltype.csv")
    per_celltype.to_csv(csv_ct, index=False)
    print(f"    [Turnover] Saved per-celltype table:      {csv_ct}")

    txt_opt = os.path.join(output_dir, f"{sample_tag}_step5_optimal_expansion.txt")
    with open(txt_opt, 'w') as fh:
        fh.write(f"optimal_expansion\t{optimal_expansion}\n")
        fh.write(f"mean_turnover\t{np.nanmean(meand_celltype)}\n")
        fh.write(f"mean_nuclei_size\t{np.nanmean(nuclimall)}\n")
    print(f"    [Turnover] Saved optimal expansion value:  {txt_opt}")

    return turnover_summ, per_celltype, optimal_expansion


# --- Main entry point ---

def run_step5(config):
    """Step 5 entry point: domain expansion + turnover analysis."""
    print("\n" + "=" * 60)
    print("[Step 5] Optimal Expansion (Ref: 4_1 Notebook)")
    print("=" * 60)

    output_dir = config["output_dir"]
    sample_tag = config["sample_tag"]

    exp_config = config.get("optimal_expansion", {})
    run_exp = exp_config.get("run_expansion", True)
    subsample_frac = exp_config.get("subsample_fraction", 0.01)
    dist_threshold = exp_config.get("distance_threshold", None)
    min_reads_per_domain = exp_config.get("min_reads_per_domain", 5000)
    diff_threshold = exp_config.get("diff_threshold", 0.1)

    if not run_exp:
        print("    - [Info] 'run_expansion' is False in config. Skipping Step 5.")
        return

    # --- 5-1. Load original transcripts (Step 0) ---
    print("\n[Step 5-1] Loading Data...")

    step0_file = os.path.join(output_dir, f"{sample_tag}.h5ad")
    if not os.path.exists(step0_file):
        parent_dir = os.path.dirname(output_dir)
        step0_file_alt = os.path.join(parent_dir, "step0_formatting", f"{sample_tag}.h5ad")
        if os.path.exists(step0_file_alt):
            step0_file = step0_file_alt

    if not os.path.exists(step0_file):
        logging.error(f"    - Step 0 output not found at {step0_file} or in sibling dir. Cannot run expansion.")
        return

    print(f"    - Loading Step 0 Data from: {step0_file}")
    adata_step0 = sc.read_h5ad(step0_file)
    if 'spots' not in adata_step0.uns:
        logging.error("    - 'spots' dataframe missing in Step 0 adata.uns.")
        return

    reads_original = adata_step0.uns['spots'].copy()
    print(f"    - Loaded {len(reads_original)} original reads.")

    # --- 5-2. Load annotated cells (Step 1/2/4) ---
    # Load annotated cells (from Step 2, 1, or 4)
    input_adata_file = None
    if 'previous_step_adata_path' in config and config['previous_step_adata_path']:
        input_adata_file = config['previous_step_adata_path']

    if not input_adata_file or not os.path.exists(input_adata_file):
        print("    - [Info] 'previous_step_adata_path' missing or invalid. Searching...")
        parent_dir = os.path.dirname(output_dir)
        possible_paths = [
            os.path.join(parent_dir, "step4_resegmentation", f"{sample_tag}_resegmented.h5ad"),
            os.path.join(parent_dir, "step1_exploration", f"{sample_tag}_step1_exploration.h5ad"),
            os.path.join(parent_dir, "step2_segmentation_free", f"{sample_tag}_step2_points2regions.h5ad")
        ]
        for p in possible_paths:
            if os.path.exists(p):
                input_adata_file = p
                break

    print(f"    - Loading Annotated Data from: {input_adata_file}")
    if not input_adata_file or not os.path.exists(input_adata_file):
        logging.error("    - Annotated adata not found. Run Step 1, 2, or 4 first.")
        return

    adata_annotated = sc.read_h5ad(input_adata_file)
    print(f"    - Loaded annotated cells: {adata_annotated.shape}")

    # --- 5-3. Map domain assignments to reads ---
    print("\n[Step 5-2] Identifying Domains & Unassigned Reads...")

    domain_key = None
    priority_keys = ['spatial_annotation', 'Class', 'leiden', 'cluster', 'graph_clusters']
    for key in priority_keys:
        if key in adata_annotated.obs.columns:
            domain_key = key
            break

    if not domain_key:
        logging.error("    - No suitable domain/cluster key found in annotated adata. Cannot assign domains.")
        return
    print(f"    - Using '{domain_key}' as domain source.")

    if 'cell_id' in adata_annotated.obs.columns:
        annotated_ids = adata_annotated.obs['cell_id']
    else:
        annotated_ids = adata_annotated.obs.index

    domain_map = dict(zip(annotated_ids, adata_annotated.obs[domain_key]))

    print("    - Mapping existing domains to reads...")
    if 'cell_id' not in reads_original.columns:
        if reads_original.index.name == 'cell_id':
            reads_original = reads_original.reset_index()
        else:
            logging.error("    - 'cell_id' missing in reads. Cannot link to cells.")
            return

    reads_original['domain'] = reads_original['cell_id'].map(domain_map)

    ct_map = dict(zip(annotated_ids, adata_annotated.obs[domain_key]))
    reads_original['initial_annotation'] = reads_original['cell_id'].map(ct_map)

    nancells = reads_original[reads_original['domain'].isna()]
    annotatedcells = reads_original[~reads_original['domain'].isna()]

    n_assigned = len(annotatedcells)
    n_unassigned = len(nancells)
    print(f"    - Assigned Reads (with domain): {n_assigned}")
    print(f"    - Unassigned Reads (to be expanded): {n_unassigned}")

    if n_unassigned == 0:
        print("    - No unassigned reads found. Expansion not needed.")
        return

    # --- 5-4. KDTree nearest-domain expansion ---
    print(f"\n[Step 5-3] Building cKDTree (Subsample Fraction: {subsample_frac})...")

    n_sample = int(np.round(n_assigned * subsample_frac))
    if n_sample < 100:
        n_sample = min(n_assigned, 1000)

    print(f"    - Subsampling {n_sample} reads as anchors...")
    annotated_sub = annotatedcells.sample(n=n_sample, random_state=42)

    coords1 = annotated_sub[['x_location', 'y_location']].values

    print("    - Constructing KDTree...")
    tree = cKDTree(coords1)

    # Query unassigned reads against KDTree
    print(f"\n[Step 5-4] Assigning {n_unassigned} reads to nearest domains...")

    coords2 = nancells[['x_location', 'y_location']].values

    chunk_size = 1000000
    n_chunks = int(np.ceil(len(coords2) / chunk_size))

    all_indices = []
    all_distances = []

    print(f"    - Processing in {n_chunks} chunks...")
    for i in range(n_chunks):
        start = i * chunk_size
        end = min((i + 1) * chunk_size, len(coords2))
        chunk = coords2[start:end]
        dists, idxs = tree.query(chunk, k=1)
        all_indices.append(idxs)
        all_distances.append(dists)

    concat_indices = np.concatenate(all_indices)
    concat_distances = np.concatenate(all_distances)

    # --- 5-5. Distance threshold filter (optional) ---
    print("\n[Step 5-5] Finalizing Assignment...")

    anchor_domains = annotated_sub['domain'].values
    assigned_domains = anchor_domains[concat_indices]

    if dist_threshold:
        print(f"    - Applying distance threshold: {dist_threshold}")
        mask_too_far = concat_distances > dist_threshold
        params_assigned = assigned_domains.copy()
        params_assigned[mask_too_far] = np.nan
        assigned_domains = params_assigned
        n_filtered = np.sum(mask_too_far)
        print(f"    - {n_filtered} reads exceeded distance threshold.")

    anchor_annotations = annotated_sub['initial_annotation'].values
    assigned_annotations = anchor_annotations[concat_indices]

    nancells = nancells.copy()
    nancells.loc[:, 'domain'] = assigned_domains
    nancells.loc[:, 'initial_annotation'] = assigned_annotations
    reads_original.loc[nancells.index, 'domain'] = assigned_domains
    reads_original.loc[nancells.index, 'initial_annotation'] = assigned_annotations

    print("    - Merging results complete.")

    # --- 5-6. Spatial map + distance histogram ---
    print("\n[Step 5-6] Generating Visualization...")
    plot_df = reads_original.sample(n=min(len(reads_original), 100000), random_state=42)

    plt.figure(figsize=(10, 10))
    plot_df_clean = plot_df.dropna(subset=['domain'])
    sns.scatterplot(
        data=plot_df_clean,
        x='x_location',
        y='y_location',
        hue='domain',
        s=1,
        linewidth=0,
        palette='tab20',
        legend=False
    )
    plt.title(f"Optimal Expansion Result (Subsample {subsample_frac})")
    plt.axis('equal')

    plot_file = os.path.join(output_dir, f"{sample_tag}_step5_expansion_map.png")
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"    - Saved Expansion Map to: {plot_file}")

    plt.figure(figsize=(6, 4))
    sns.histplot(concat_distances, bins=50, kde=True, color='orange')
    plt.title("Distance to Nearest Domain-Anchor")
    plt.xlabel("Distance (pixels/units)")
    plt.ylabel("Count")
    dist_file = os.path.join(output_dir, f"{sample_tag}_step5_expansion_distances.png")
    plt.savefig(dist_file)
    plt.close()
    print(f"    - Saved Distance Histogram to: {dist_file}")

    # --- 5-7. Save expanded transcripts ---
    output_file = os.path.join(output_dir, f"{sample_tag}_step5_expanded_transcripts.csv")
    print(f"    - Saving Expanded Transcripts to: {output_file}")
    reads_original.to_csv(output_file, index=False)

    # --- 5-8. Turnover analysis ---
    print("\n[Step 5-8] Running Turnover/Crossover Analysis...")

    # Separate reads by cell assignment:
    #   reads_not_assigned = cell_id == -1 or 'UNASSIGNED' (background)
    #   reads_assigned = valid cell_id (segmented cells)
    print("    - Preparing reads for turnover analysis...")

    if reads_original['cell_id'].dtype in [np.int64, np.int32, np.float64, np.float32, int, float]:
        mask_not_assigned = reads_original['cell_id'] == -1
    else:
        mask_not_assigned = (
            reads_original['cell_id'].astype(str).str.upper().isin(['-1', 'UNASSIGNED', 'NAN', ''])
        )

    reads_not_assigned = reads_original[mask_not_assigned].copy()
    reads_assigned_all = reads_original[~mask_not_assigned].copy()

    print(f"    - Reads assigned to cells:   {len(reads_assigned_all)}")
    print(f"    - Reads not assigned (background): {len(reads_not_assigned)}")

    # Compute distance from each assigned read to its cell centroid
    if 'x_centroid' in adata_annotated.obs.columns and 'y_centroid' in adata_annotated.obs.columns:
        cx_map = dict(zip(annotated_ids, adata_annotated.obs['x_centroid']))
        cy_map = dict(zip(annotated_ids, adata_annotated.obs['y_centroid']))
    elif 'x_centroid' in adata_step0.obs.columns and 'y_centroid' in adata_step0.obs.columns:
        if 'cell_id' in adata_step0.obs.columns:
            cx_map = dict(zip(adata_step0.obs['cell_id'], adata_step0.obs['x_centroid']))
            cy_map = dict(zip(adata_step0.obs['cell_id'], adata_step0.obs['y_centroid']))
        else:
            cx_map = dict(zip(adata_step0.obs.index, adata_step0.obs['x_centroid']))
            cy_map = dict(zip(adata_step0.obs.index, adata_step0.obs['y_centroid']))
    else:
        print("    - [Warning] Centroid data (x_centroid, y_centroid) missing. Skipping Turnover Analysis.")
        print("\n=== Step 5 Optimal Expansion Complete ===")
        return

    reads_assigned_val = reads_assigned_all[reads_assigned_all['cell_id'].isin(cx_map.keys())].copy()
    print(f"    - Reads with centroid info: {len(reads_assigned_val)}")

    if len(reads_assigned_val) == 0:
        print("    - [Warning] No reads matched to cells with centroids. Skipping Turnover.")
        print("\n=== Step 5 Optimal Expansion Complete ===")
        return

    reads_assigned_val['x_cell'] = reads_assigned_val['cell_id'].map(cx_map).astype(float)
    reads_assigned_val['y_cell'] = reads_assigned_val['cell_id'].map(cy_map).astype(float)
    reads_assigned_val['distance'] = np.sqrt(
        (reads_assigned_val['x_location'] - reads_assigned_val['x_cell']) ** 2 +
        (reads_assigned_val['y_location'] - reads_assigned_val['y_cell']) ** 2
    )
    reads_assigned_val['distance'] = reads_assigned_val['distance'].round(0)

    if 'initial_annotation' not in reads_assigned_val.columns or reads_assigned_val['initial_annotation'].isna().all():
        reads_assigned_val['initial_annotation'] = reads_assigned_val['cell_id'].map(
            dict(zip(annotated_ids, adata_annotated.obs[domain_key]))
        )

    if 'feature_name' not in reads_assigned_val.columns:
        print("    - [Warning] 'feature_name' missing in reads. Skipping turnover.")
        print("\n=== Step 5 Optimal Expansion Complete ===")
        return

    turnover_summ, per_celltype, optimal_expansion = calculate_turnover(
        reads_assigned=reads_assigned_val,
        reads_not_assigned=reads_not_assigned,
        output_dir=output_dir,
        sample_tag=sample_tag,
        min_reads_per_domain=min_reads_per_domain,
        diff_threshold=diff_threshold
    )

    if per_celltype is not None:
        print(f"\n    Turnover per-celltype summary:")
        print(per_celltype.to_string(index=False))

    print("\n=== Step 5 Optimal Expansion Complete ===")


if __name__ == "__main__":
    pass
