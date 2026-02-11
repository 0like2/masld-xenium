# Step 6: Segmentation Benchmark -- Notebook vs Pipeline Review

**Pipeline file:** `pipeline/xenium_step6_segmentation_benchmark.py` (622 lines)
**Notebook:** `notebooks/5_segmentation_benchmark/5_1_Compare_Clustering on_different_segmentations.ipynb`
**Notebook metrics:** `notebooks/5_segmentation_benchmark/metrics.py` (352 lines)
**Notebook util:** `notebooks/5_segmentation_benchmark/util.py` (217 lines)
**Pipeline metrics:** `pipeline/benchmark_utils/metrics.py` (352 lines)
**Pipeline util:** `pipeline/benchmark_utils/util.py` (217 lines)
**Shared library:** `xb/preprocessing.py` (`preprocess_adata()`)

**Reviewer:** Claude Opus 4.6
**Date:** 2025-02-11

---

## Executive Summary

Step 6 benchmarks multiple segmentation methods (nuclei-only, Cellpose, expansion, Baysor) by
concatenating their outputs into a single AnnData object, preprocessing jointly, transferring
cell type annotations from a reference, computing quality metrics, and generating comparison
visualizations. The pipeline faithfully reproduces the core flow (load, concat, preprocess, annotate,
plot) and adds several improvements (Baysor integration with dry-run, configurable parameters,
robust column detection). However, several **critical and medium-severity** gaps exist in the
metrics computation and annotation transfer that reduce the pipeline's ability to quantitatively
rank segmentation methods -- the central purpose of this step.

---

## Sub-step Analysis

### 6-1. Baysor Execution

**Notebook Intent:**
The notebook loads a pre-computed Baysor h5ad file (`counts_custom-0_baysor-82_area-0_majority-0.h5ad`)
rather than running Baysor itself. Baysor was executed externally before the notebook was run.

**Paper Methodology:**
Baysor v0.7.x is used with prior segmentation masks and configurable scale/min-molecules parameters.

**Pipeline Implementation (lines 34-174):**
- `prep_xenium_data_for_baysor()` (lines 36-78): Converts Xenium transcripts to Baysor CSV format
  (gene, x, y). Handles both CSV and Parquet inputs.
- `run_baysor()` (lines 83-174): Full Baysor invocation including:
  - Config-driven parameters: `scale`, `min_molecules_per_cell`, `prior_segmentation_confidence`
  - Prior segmentation TIF support (line 122-135)
  - `dry_run` mode producing empty placeholder CSV (lines 140-148)
  - Binary existence check via `shutil.which()` (line 151)
  - Comprehensive error handling with specific messages

**Differences:**
- Pipeline **improves** over notebook by automating Baysor execution end-to-end.
- Dry-run mode is well-designed for CI/testing.
- No z-coordinate usage (only x, y) -- notebook also does not use z.

**Severity:** OK -- Well implemented, exceeds notebook scope.

---

### 6-2. Load Segmentation Results

**Notebook Intent:**
Loads two pre-computed h5ad files: nuclei-only (Xenium default segmentation with nuclear reads)
and Baysor output. Subsets nuclei data to a region of interest matching Baysor's ROI. Adds
`segmentation` and `sample` obs columns.

**Paper Methodology:**
Compares nuclei-only, Cellpose, boundary-expansion, and Baysor segmentations on the same tissue ROI.

**Pipeline Implementation (lines 505-546):**
- Loads up to four methods: `input_nuclei`, `input_advanced` (Cellpose), `input_expansion`, Baysor
- Nuclei and Cellpose loaded via `sc.read()` for h5ad files (lines 514-524)
- Expansion loaded via `load_transcripts_as_adata()` (lines 526-533) which aggregates CSV transcripts
- Baysor loaded via `load_transcripts_as_adata()` (lines 535-542)
- `load_transcripts_as_adata()` (lines 179-229): Robust CSV loading with multiple column-name
  conventions (`cell_id`/`cell`, `feature_name`/`gene`), unassigned transcript filtering, and
  `pd.crosstab()` aggregation

**Differences:**
1. Notebook subsets nuclei data to a spatial ROI (`x_interval`, `y_interval`) to match Baysor's region.
   Pipeline does not perform spatial subsetting -- assumes inputs are already aligned. If they are not,
   cell counts will be incomparable.
2. Notebook rescales coordinates by pixel-to-um factor (`*4.70588`). Pipeline does not handle
   coordinate unit conversion.
3. Pipeline supports 4 methods vs notebook's 2 -- this is an improvement.
4. `load_transcripts_as_adata()` does not store raw spots in `adata.uns['spots']`, which means
   `proportion_of_assigned_reads()` will never be callable (it requires `adata.uns['spots']`).
   The pipeline checks for `'spots' in subset.uns` at line 602, so it gracefully skips, but
   the metric is effectively unreachable for CSV-loaded methods.

**Severity:** MEDIUM -- Spatial ROI subsetting gap can lead to mismatched comparisons. Missing
`uns['spots']` makes `proportion_of_assigned_reads` unreachable.

---

### 6-3. Concatenate and Preprocess

**Notebook Intent:**
Concatenates nuclei and Baysor AnnData objects, computes `total_counts` and `expressed_genes`,
then calls `preprocess_adata()` with explicit parameters:
```python
clustering_params = {
    'normalization_target_sum': 100,
    'min_counts_x_cell': 40,
    'min_genes_x_cell': 15,
    'scale': False,
    'clustering_alg': 'louvain',
    'resolutions': [1.0],
    'n_neighbors': 15,
    'umap_min_dist': 0.1,
    'n_pcs': 0
}
```

**Paper Methodology:**
Standard scanpy preprocessing: filter, normalize, log1p, HVG selection, PCA, neighbors, clustering,
UMAP. Louvain algorithm specified (but implementation falls back to Leiden -- see below).

**Pipeline Implementation (lines 297-340, `preprocess_benchmark()`):**
- Config-driven parameters via `config["benchmark"]["preprocessing"]` (lines 300-309)
- Defaults match notebook: `target_sum=100`, `min_counts=40`, `min_genes=15`, `n_neighbors=15`,
  `n_pcs=0`, `umap_min_dist=0.1`, `resolution=1.0`, `scale=False`
- Stores `layers['raw']` before and after filtering (lines 317, 325)
- Runs: `filter_cells` -> `normalize_total` -> `log1p` -> optional `scale` -> `pca` -> `neighbors`
  -> `umap` -> `leiden` (lines 320-338)

**Differences:**

| Aspect | Notebook (`preprocess_adata`) | Pipeline (`preprocess_benchmark`) |
|--------|------------------------------|----------------------------------|
| Clustering algorithm label | `'louvain'` in config, key_added=`'louvain_1.0'` | Hardcoded `'leiden'`, key_added=`'leiden'` |
| Actual algorithm used | **Leiden** (line 142 of xb/preprocessing.py: `sc.tl.leiden` even when `clustering_alg=='louvain'`) | **Leiden** (line 338) |
| HVG selection | **YES** -- `sc.pp.highly_variable_genes(min_mean=0.3, max_mean=7, min_disp=-0.5)` (line 129 of xb/preprocessing.py) | **NO** -- omitted entirely |
| Marker gene analysis | `sc.tl.rank_genes_groups()` + dotplot + DEG plots | Not included |
| Resolution loop | Supports list of resolutions | Single resolution only |
| Spatial plots per sample | Per-sample spatial cluster maps | Handled separately in 6-6b |

Key finding: **Both notebook and pipeline actually use Leiden**, despite the notebook labeling it
"louvain". The pipeline comment at line 337 is correct. However, three material differences remain:

1. **HVG selection is missing in the pipeline.** The notebook always runs HVG filtering
   (`min_mean=0.3, max_mean=7, min_disp=-0.5`) which reduces genes to informative features before
   PCA. Without HVG, PCA will be computed on all genes including noise, potentially degrading
   cluster quality.

2. **Clustering key naming differs.** Notebook stores clusters as `louvain_1.0`, pipeline stores as
   `leiden`. Downstream annotation transfer uses `cluster_key='leiden'` (line 569), which is
   consistent within the pipeline, but any external code expecting `louvain_*` keys will break.

3. **Marker gene ranking is absent.** The notebook runs `rank_genes_groups` with Wilcoxon test and
   generates dotplots/DEG lists. This is useful for interpreting clusters but is not strictly required
   for benchmarking.

**Severity:** MEDIUM -- HVG selection gap may affect cluster quality and downstream comparisons.

---

### 6-4. Annotation Transfer

**Notebook Intent:**
Two-phase annotation transfer:
1. **Direct cell ID mapping:** For nuclei cells that exist in the reference, map `cell_id` to `Class`
   directly via dictionary lookup (`id2class`).
2. **Cluster consensus:** Build a crosstab of `louvain_1.0` clusters vs `Class` labels, pick the
   dominant class per cluster (`idxmax`), then assign that class to ALL cells in the cluster
   (including Baysor cells that had no direct mapping). This ensures Baysor cells get annotations
   through cluster-level consensus.

**Paper Methodology:**
Transfer known cell type annotations from a previously annotated reference (Figure 1 dataset) to
the benchmark dataset so that cell type proportions can be compared across segmentation methods.

**Pipeline Implementation (lines 234-292, `annotate_by_majority_voting()`):**
- Uses **kNN majority voting** in PCA space (sklearn `NearestNeighbors`, lines 262-265)
- For each target cell, finds k nearest reference cells and assigns the most common label (lines 267-275)
- Also does per-cluster consensus: majority vote within each Leiden cluster (lines 278-289)
- Stores `celltype_majority` (per-cell) and `celltype_cluster` (per-cluster) in obs

**Differences:**

| Aspect | Notebook | Pipeline |
|--------|----------|----------|
| Method | Direct ID mapping + cluster crosstab consensus | kNN majority voting in PCA space |
| Requires shared cell IDs | YES (nuclei cells must exist in reference) | NO (works in embedding space) |
| Handles unseen cells | Only via cluster-level consensus | Per-cell kNN voting covers all cells |
| Confidence score | Notebook `util.py` `run_majority_voting()` computes `ct_majority_cert` (proportion of reads supporting winner label) | **NOT computed** -- no confidence/certainty metric |
| Certainty filtering | `generate_adata()` applies `ct_certainty_threshold=0.7` to set low-confidence cells to "Unknown" | Not available |
| Gene intersection check | Not explicit | Explicit shared gene check with error handling (lines 252-254) |

The pipeline's kNN approach is arguably **more general and robust** than the notebook's direct ID
mapping, since it does not require cell IDs to be shared between the reference and target. However,
the absence of a confidence score is a meaningful gap: without it, there is no way to identify or
filter unreliable annotations, and all annotations carry equal weight in downstream analyses.

**Severity:** MEDIUM -- Missing confidence scores prevent filtering unreliable annotations.

---

### 6-5. Benchmark Metrics

**Notebook Intent / Paper Methodology:**
The notebook's `metrics.py` defines a comprehensive metric suite for quantitative segmentation
comparison:

| Metric | Function | Purpose |
|--------|----------|---------|
| Proportion of assigned reads | `proportion_of_assigned_reads()` | Transcript capture efficiency |
| Number of cells | `number_of_cells()` | Segmentation granularity |
| 5th percentile reads/cell | `percentile_5th_reads_cells()` | Low-quality cell detection |
| Median reads/cell | `median_reads_cells()` | Central tendency of transcript capture |
| 5th percentile genes/cell | `percentile_5th_genes_cells()` | Gene diversity floor |
| Median genes/cell | `median_genes_cells()` | Central tendency of gene diversity |
| Negative marker purity (cells) | `negative_marker_purity_cells()` | Read leakage between cells (cell-level) |
| Negative marker purity (reads) | `negative_marker_purity_reads()` | Read leakage between cells (read-level) |
| Adjusted Rand Index | `rand_idx()` | Agreement between segmentation assignments |

These metrics are central to the paper's Table 1 and supplementary figures comparing segmentation
methods.

**Pipeline Implementation (lines 584-608):**
```python
for seg_method in adata.obs['segmentation'].unique():
    subset = adata[adata.obs['segmentation'] == seg_method]
    n_cells = subset.shape[0]
    results[f'{seg_method}_n_cells'] = n_cells

    if 'raw' in subset.layers:
        results[f'{seg_method}_median_reads'] = metrics.median_reads_cells(subset)
        results[f'{seg_method}_median_genes'] = metrics.median_genes_cells(subset)

    if 'spots' in subset.uns:
        results[f'{seg_method}_assigned_prop'] = metrics.proportion_of_assigned_reads(subset)
```

**Differences:**

| Metric | Available in `benchmark_utils/metrics.py` | Called in `run_step6()` | Gap |
|--------|-------------------------------------------|------------------------|-----|
| `number_of_cells` | YES | YES (inline `subset.shape[0]`) | OK |
| `median_reads_cells` | YES | YES (line 594) | OK |
| `median_genes_cells` | YES | YES (line 595) | OK |
| `proportion_of_assigned_reads` | YES | Conditionally (line 603) | Unreachable* |
| `percentile_5th_reads_cells` | YES | **NOT CALLED** | MISSING |
| `percentile_5th_genes_cells` | YES | **NOT CALLED** | MISSING |
| `negative_marker_purity_cells` | YES | **NOT CALLED** | MISSING |
| `negative_marker_purity_reads` | YES | **NOT CALLED** | MISSING |
| `rand_idx` | YES | **NOT CALLED** | MISSING |
| Silhouette score | NO | NO | MISSING |
| Calinski-Harabasz score | NO | NO | MISSING |
| Davies-Bouldin score | NO | NO | MISSING |

*`proportion_of_assigned_reads` requires `adata.uns['spots']` which is never populated by the
pipeline's data loading path.

**Critical missing metrics:**

1. **Negative marker purity (cells and reads):** These measure read leakage between neighboring cells,
   which is the primary indicator of segmentation boundary quality. The functions exist in
   `pipeline/benchmark_utils/metrics.py` lines 109-326 but are never invoked. They require:
   - A scRNA-seq reference (`adata_sc`) with matching cell types
   - `layers['raw']` on both spatial and reference data
   - `celltype` annotation in `.obs`
   The pipeline has the annotation transfer step (6-4) that could provide cell types, and the
   reference AnnData is already loaded at line 559. The infrastructure exists but the call is missing.

2. **Adjusted Rand Index:** Measures agreement between different segmentation methods at the
   transcript level. The `rand_idx()` function (lines 329-351) takes a DataFrame where each column
   is a segmentation method's cell assignments per transcript. The pipeline does not construct this
   DataFrame or call the function.

3. **5th percentile metrics:** `percentile_5th_reads_cells` and `percentile_5th_genes_cells` are
   trivial to add -- they use the same `layers['raw']` data already available.

4. **Clustering quality scores (Silhouette, Calinski-Harabasz, Davies-Bouldin):** Standard
   sklearn metrics for evaluating cluster separation quality. Neither notebook nor pipeline compute
   these, but they are commonly used in the literature and would strengthen the quantitative
   comparison.

**Severity:** CRITICAL -- The pipeline computes only 3 of 9 available metrics. The most
scientifically important metrics (negative marker purity, Rand index) that distinguish segmentation
quality are implemented in the utility module but never called. Without these, the pipeline cannot
produce the paper's key quantitative comparison table.

---

### 6-6a. UMAP per Method

**Notebook Intent:**
UMAP colored by `segmentation` method (nuclei vs Baysor) with publication-quality settings (dpi=500,
PDF output, custom palette `['#7BB542','#B84E9D']`).

**Pipeline Implementation (lines 346-359, `_save_umap()`):**
- Two-panel figure: segmentation method (left), Leiden clusters (right)
- PNG output at dpi=150
- Uses scanpy default palette

**Differences:**
- Pipeline adds Leiden cluster coloring (second panel) -- an improvement.
- Pipeline uses PNG (dpi=150) vs notebook's PDF (dpi=500). Lower resolution and non-vector format.
- Pipeline does not support custom color palette via config.
- Pipeline colors by `segmentation` (matching notebook).

**Severity:** LOW -- Functional, minor quality/format differences.

---

### 6-6b. Spatial Scatter

**Notebook Intent:**
Spatial scatter plot using `xb.plotting.map_of_clusters()` colored by cell type (`Class`),
one panel per sample/segmentation.

**Pipeline Implementation (lines 362-404, `_save_spatial_map()`):**
- Multi-panel spatial scatter, one panel per segmentation method
- Auto-detects coordinate columns (`x_centroid`/`y_centroid`, `x_location`/`y_location`, `x`/`y`)
- Colored uniformly per method (not by cell type)
- Inverts y-axis for correct image orientation

**Differences:**
- Notebook colors by **cell type**, pipeline colors by **method** (all dots same color per panel).
  This means the pipeline spatial plot shows where cells are but not what type they are -- less
  informative for spatial comparison.
- Pipeline does not use `xb.plotting.map_of_clusters()` -- reimplements with matplotlib.
- Pipeline handles missing coordinates gracefully.

**Severity:** LOW -- Functional but less informative coloring.

---

### 6-6c. Cell Type Barplot

**Notebook Intent:**
Bar plot of cell type frequencies (`Class`) across segmentation methods using `pd.crosstab` and
pandas `.plot(kind='bar')`.

**Pipeline Implementation (lines 407-447, `_save_celltype_barplot()`):**
- Grouped bar plot of cell type proportions per segmentation method
- Uses matplotlib grouped bars with proper offset calculation
- Normalizes to proportions (notebook uses raw counts)

**Differences:**
- Pipeline uses **proportions** (line 423: `count / totals`), notebook uses **raw counts**.
  Proportions are more appropriate for comparing methods with different total cell counts.
- Pipeline uses `celltype_majority` key (from kNN), notebook uses `Class` (from cluster consensus).
- Pipeline handles missing `celltype_majority` gracefully.

**Severity:** OK -- Pipeline approach (proportions) is actually more correct than notebook (raw counts).

---

### 6-6d. Counts Violin

**Notebook Intent:**
Violin plot of `total_counts` per cell grouped by segmentation method using `seaborn.violinplot`.

**Pipeline Implementation (lines 450-486, `_save_counts_violin()`):**
- Computes total counts from `layers['raw']` (handles sparse matrices)
- Uses matplotlib `violinplot` with means and medians shown
- Proper handling of missing raw layer

**Differences:**
- Pipeline uses matplotlib instead of seaborn -- slightly different aesthetics.
- Pipeline shows both means and medians (more informative).
- Pipeline handles sparse raw layer explicitly (line 458-461).

**Severity:** OK -- Functionally equivalent with minor visual differences.

---

## Additional Findings

### Notebook metrics.py vs Pipeline metrics.py -- NOT identical

Despite the project memory stating these are "same as notebook", there is one important difference:

| Location | Line 206 / 320 | Effect |
|----------|----------------|--------|
| `notebooks/.../metrics.py` | `.cip(0, None)` | **TYPO** -- `cip` is not a valid pandas method; would raise `AttributeError` at runtime |
| `pipeline/benchmark_utils/metrics.py` | `.clip(0, None)` | **FIXED** -- correct method name |

The pipeline version has already fixed this bug. This means the notebook's `negative_marker_purity_*`
functions with `pipeline_output=False` would crash in the notebook but work correctly in the pipeline
(if they were called).

### Notebook util.py vs Pipeline util.py -- NOT identical

| Location | Difference |
|----------|------------|
| `notebooks/.../util.py` line 176 | `from descartes import PolygonPatch` (deprecated library, unused import) |
| `pipeline/benchmark_utils/util.py` line 176 (absent) | Import removed |

The pipeline correctly removed the unused `descartes` dependency.

### HVG Selection Gap Detail

The notebook's `preprocess_adata()` (xb/preprocessing.py line 129) always runs:
```python
sc.pp.highly_variable_genes(adata, min_mean=0.3, max_mean=7, min_disp=-0.5)
```

Although the notebook does not subset to HVG genes (no `adata = adata[:, adata.var.highly_variable]`
line follows), the HVG annotation is stored in `adata.var.highly_variable` which scanpy's PCA
implementation uses automatically when present. The pipeline's `preprocess_benchmark()` never calls
`sc.pp.highly_variable_genes()`, so PCA runs on all genes. With spatial transcriptomics panels
(typically 300-500 genes), HVG filtering may have minimal impact, but for larger panels it could
matter.

### Annotation Transfer Method Comparison

The notebook uses a simple but effective two-step method:
1. Map shared cell IDs directly from reference
2. Use cluster-level crosstab consensus for unmapped cells

The pipeline uses a more general kNN approach:
1. Fit kNN on reference PCA embeddings
2. Query with target PCA embeddings
3. Majority vote among k nearest neighbors

The kNN approach is more robust when cell IDs do not overlap (e.g., Baysor creates entirely new
cell identities), but it depends on the PCA embedding quality and gene overlap between reference
and target. The pipeline correctly checks for shared genes (line 252) and PCA existence (line 241).

---

## Prioritized Fix List

### CRITICAL

**C1. Call negative marker purity metrics (lines 584-608)**
- Add calls to `metrics.negative_marker_purity_cells()` and `metrics.negative_marker_purity_reads()`
- Requires scRNA-seq reference with `layers['raw']` and `obs['celltype']`
- Could use the already-loaded `adata_ref` (line 559) if it has a raw layer
- Add config key `benchmark.scrna_reference` for scRNA-seq reference path
- Estimated effort: ~30 lines

**C2. Call Rand Index metric (lines 584-608)**
- Construct a per-transcript assignment DataFrame with one column per segmentation method
- Requires loading raw transcript CSVs and aligning cell assignments
- Call `metrics.rand_idx()` on the DataFrame
- Estimated effort: ~50 lines (need transcript-level alignment)

**C3. Call 5th percentile metrics (lines 584-608)**
- Trivial addition alongside existing `median_reads_cells` / `median_genes_cells` calls:
  ```python
  results[f'{seg_method}_p5_reads'] = metrics.percentile_5th_reads_cells(subset)
  results[f'{seg_method}_p5_genes'] = metrics.percentile_5th_genes_cells(subset)
  ```
- Estimated effort: ~4 lines

### MEDIUM

**M1. Add HVG selection to `preprocess_benchmark()` (line 333)**
- Add HVG step matching notebook parameters:
  ```python
  sc.pp.highly_variable_genes(adata, min_mean=0.3, max_mean=7, min_disp=-0.5)
  ```
- Make configurable via `preprocessing.hvg` boolean and `preprocessing.hvg_params`
- Estimated effort: ~10 lines

**M2. Add confidence scores to `annotate_by_majority_voting()` (line 275)**
- Compute `votes_for_winner / k` as certainty for each cell:
  ```python
  max_count = np.max(counts)
  cell_labels.append(values[np.argmax(counts)])
  cell_confidences.append(max_count / n_neighbors)
  ```
- Store as `adata_target.obs['celltype_confidence']`
- Optionally filter low-confidence cells to "Unknown" via config threshold
- Estimated effort: ~15 lines

**M3. Make clustering algorithm configurable (line 337-338)**
- Add `preprocessing.clustering_alg` config key (default `'leiden'`)
- Support both `sc.tl.leiden` and `sc.tl.louvain` based on config
- Estimated effort: ~10 lines

**M4. Populate `uns['spots']` for CSV-loaded methods (lines 179-229)**
- In `load_transcripts_as_adata()`, store the raw transcript DataFrame in `adata.uns['spots']`
  before filtering unassigned reads
- This enables `proportion_of_assigned_reads()` to work
- Estimated effort: ~5 lines

**M5. Spatial ROI subsetting (lines 514-518)**
- Add optional `benchmark.roi` config with `x_interval` and `y_interval`
- Subset all loaded AnnData objects to the same spatial region before concatenation
- Prevents comparing cells from different tissue areas
- Estimated effort: ~20 lines

### LOW

**L1. Add clustering quality scores**
- Compute Silhouette, Calinski-Harabasz, Davies-Bouldin scores per method using
  `sklearn.metrics` on PCA embeddings + Leiden labels
- Not in the original notebook but would strengthen the benchmark
- Estimated effort: ~20 lines

**L2. Area normalization support (util.py lines 122-217)**
- `normalize_by_area()` and `calculate_alpha_area()` are available in `benchmark_utils/util.py`
  but never called
- Useful for Baysor which produces variable cell sizes
- Estimated effort: ~15 lines + optional dependency (`alphashape`)

**L3. Marker gene analysis**
- Add `sc.tl.rank_genes_groups()` per segmentation method to identify method-specific markers
- Present in notebook's `preprocess_adata()` but not in pipeline
- Estimated effort: ~15 lines

**L4. Spatial plot colored by cell type**
- Current spatial plot colors by method (uniform per panel)
- Add option to color by `celltype_majority` for more informative comparison
- Estimated effort: ~20 lines

**L5. Publication-quality output format**
- Support PDF output and configurable DPI for publication figures
- Notebook uses dpi=500/700 PDF; pipeline uses dpi=150 PNG
- Estimated effort: ~10 lines

---

## Summary Table

| Sub-step | Description | Status | Severity |
|----------|-------------|--------|----------|
| 6-1 | Baysor execution | WELL IMPLEMENTED | OK |
| 6-2 | Load segmentation results | Functional with gaps | MEDIUM |
| 6-3 | Concatenate & preprocess | Missing HVG | MEDIUM |
| 6-4 | Annotation transfer | Missing confidence scores | MEDIUM |
| 6-5 | Benchmark metrics | **5 of 8 metrics not called** | **CRITICAL** |
| 6-6a | UMAP per method | Functional | LOW |
| 6-6b | Spatial scatter | Less informative coloring | LOW |
| 6-6c | Cell type barplot | Improved (proportions) | OK |
| 6-6d | Counts violin | Equivalent | OK |

**Bottom line:** The pipeline's structural implementation (Baysor integration, data loading,
preprocessing, visualization) is solid and in several cases improves upon the notebook. The
critical gap is in metrics computation: the `benchmark_utils/metrics.py` module contains all
necessary functions, but `run_step6()` only calls 3 of them. Adding calls to the existing
functions (especially negative marker purity and Rand index) is the highest-priority fix and
would bring the pipeline to full feature parity with the notebook's analytical capabilities.
