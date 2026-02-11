# Step 4: Techniques Comparison -- Notebook vs Pipeline Review

**Reviewed:** 2026-02-11
**Pipeline file:** `pipeline/xenium_step4_techniques_comparison.py` (665 lines)
**Notebooks reviewed:**
- `3_3_efficiency_between_methods.ipynb`
- `3_4_negative_marker_purity_for_specificity.ipynb`
- `3_5_Computing_positivity_after_preprocessing_for_all_ST_techs.ipynb`
- `3_6_Diffussion_on_resegmented_data.ipynb`
- `3_7_Xenium_vs_Visium_comparison.ipynb`

---

## Overview

Step 4 performs four main analyses (efficiency, specificity, positivity, diffusion) that
compare resegmented data against original segmentation and, in the notebook context,
against multiple spatial transcriptomics (ST) technologies (CosMx, HybrISS, MERFISH,
Resolve Biosciences, Vizgen, Xenium) using scRNAseq as the ground-truth reference.

The pipeline collapses this multi-technology comparison into a two-condition comparison
(Resegmented vs Original), which is appropriate for a single-sample pipeline but loses
the cross-technology benchmarking that is central to the paper's narrative.

---

## Sub-step 4-2: Efficiency Analysis (Notebook 3_3)

### Notebook Intent
Quantify detection efficiency of each ST technology by computing per-gene expression
ratios relative to scRNAseq. The notebook operates across 6 ST technologies + scRNAseq,
filters to genes present in >= 3 datasets, and analyzes three brain regions separately
(Cortex, Hippocampus, Thalamus). Key outputs:
- NxN pairwise scatter plots of per-gene mean expression between ALL method pairs
- Per-gene ST/scRNAseq ratio boxplots, one per region
- Probe-count vs efficiency correlation analysis (Xenium-specific)
- Region-specific ROI spatial plots confirming annotation consistency

### Paper Methodology
`median_calculator()` from `xb/comparing.py` computes per-gene median expression ratios
for each technology vs scRNAseq. `combine_med()` aggregates these into a tidy DataFrame.
Per-gene ratios are plotted as boxplots with strip overlay per method, sorted by median.
The pairwise scatter (NxN grid) shows log-log per-gene mean expression for each pair of
methods with identity line.

### Pipeline Implementation
- **Lines 195-254:** `analyze_efficiency()` -- Transcripts/genes per cell histograms and
  summary statistics CSV. Compares Resegmented vs Original using QC metrics.
- **Lines 256-301:** Expression ratio (ST/scRNAseq) with CPM normalization. Computes
  `st_mean / sc_mean` per gene for common genes. Saves ratio CSV and histogram.
- **Lines 303-340:** Region-based efficiency breakdown using auto-detected region column.
  Generates boxplots and CSV of per-region stats.

### Differences

| Aspect | Notebook | Pipeline | Severity |
|--------|----------|----------|----------|
| Number of methods compared | 6 ST + scRNAseq | 2 (Reseg vs Original) | Design choice |
| Pairwise NxN scatter grid | Full method-vs-method grid with log-log axes | Not implemented | MEDIUM |
| Per-region analysis | Explicit Cortex/Hippocampus/Thalamus with curated cluster lists | Auto-detected region column, single-pass | LOW |
| Probe-count correlation | Xenium probe count vs efficiency ratio correlation | Not implemented | LOW |
| Normalization | `median_calculator()` using raw counts, 100x scaling | CPM (1e6) normalization | MEDIUM |
| Multi-region ratio CSV export | Separate CSV per region (cortex, hippocampus, thalamus) | Single combined region CSV | LOW |

### Severity: MEDIUM
The missing pairwise scatter plots (NxN grid) are a key paper figure. However, this
requires loading all 6+ technologies simultaneously, which is outside the single-sample
pipeline design. The CPM normalization difference vs the notebook's 100x scaling is a
methodological divergence that could affect numerical ratio values, though the relative
ordering should be preserved.

---

## Sub-step 4-3: Specificity Analysis (Notebook 3_4)

### Notebook Intent
Quantify specificity of each ST technology using Negative Marker Purity (NMP/NCP).
Genes expected to NOT be co-expressed (based on scRNAseq reference) are checked for
spurious co-expression in ST data. The notebook:
- Filters genes to those with efficiency ratio < 10 (removes highly inflated genes)
- Computes NMP per gene per method using `negative_marker_purity_coexpression()` from
  `xb/calculating.py`
- Returns per-gene, per-celltype purity breakdowns (6 return values in notebook mode)
- Produces the critical **efficiency vs specificity scatter plot** (ratio on x-axis,
  purity on y-axis)

### Paper Methodology
`negative_marker_purity_coexpression()` with `pipeline_output=False` returns 6 values:
`(negative_marker_purity, purity_per_gene, purity_per_celltype, lowvals_sc, lowvals_sp,
commongenes)`. The notebook uses this to build a per-gene purity DataFrame, merges it
with the efficiency ratios from 3_3, and plots the combined scatter.

### Pipeline Implementation
- **Lines 349-367:** NMP score calculation. Calls the inlined
  `_negative_marker_purity_coexpression()` with `pipeline_output=True`, which returns
  only the scalar NMP score. Writes scores to a text file.
- **Lines 369-397:** Gene-gene correlation heatmap for top 50 expressed genes. Uses
  `np.corrcoef()` on the expression matrix, not the co-expression calculation.

### Differences

| Aspect | Notebook | Pipeline | Severity |
|--------|----------|----------|----------|
| NMP return detail | Per-gene and per-celltype purity breakdown (6 return values) | Scalar NMP only (`pipeline_output=True`) | HIGH |
| Gene filtering | Filters to genes with efficiency ratio < 10 | No pre-filtering | MEDIUM |
| Efficiency vs Specificity scatter | Key paper figure combining ratio and purity per gene | Not implemented | MEDIUM |
| Per-method NMP boxplot | Per-gene purity boxplot+strip per method, sorted by median | Not implemented (scalar only) | MEDIUM |
| Gene-gene correlation | Uses co-expression (fraction of positive cells) | Uses Pearson correlation (`np.corrcoef`) | LOW |
| Inlined function signature | `negative_marker_purity_coexpression()` returns 6 values | `_negative_marker_purity_coexpression()` returns 1 or 3 values (missing `lowvals_sc`, `lowvals_sp`, `commongenes`) | HIGH |

### Severity: HIGH
The inlined function at lines 47-97 is a simplified version of `xb.calculating.negative_marker_purity_coexpression()`. The original notebook version returns 6 values
including `lowvals_sc`, `lowvals_sp`, and `commongenes`, which are needed to build the
per-gene purity breakdown. The pipeline's inlined version returns either a scalar or 3
values (NMP score, purity_per_gene, purity_per_celltype), losing the intermediate data
needed for the efficiency-vs-specificity scatter plot. Even with `pipeline_output=False`,
the pipeline cannot reproduce the paper figure.

---

## Sub-step 4-4: Positivity Analysis (Notebook 3_5)

### Notebook Intent
Preprocess all ST datasets with identical parameters, cluster, and compare gene detection
rates (positivity = fraction of cells with non-zero expression) across technologies.
The notebook:
- Applies identical preprocessing to all 6 ST datasets:
  `filter_cells(min_counts=10, min_genes=3)`, `normalize_total(target_sum=None)`,
  `log1p`, `neighbors(n_neighbors=8, n_pcs=0)`, Leiden at 3 resolutions (0.6, 1.4, 2.2),
  UMAP(`min_dist=0.1`)
- Computes per-gene detection rates across all methods
- Identifies optimal cluster per gene (highest mean expression)
- Plots violin plots for marker genes in their best cluster, across methods

### Pipeline Implementation
- **Lines 400-419:** Gene positivity (fraction of positive cells) per dataset. Saves
  top 50 genes CSV per dataset.
- **Lines 423-430:** Positivity distribution histogram comparing datasets.
- **Lines 432-441:** Preprocessing: `normalize_total()` (default target_sum=None is
  actually used by pipeline without arguments -- but the subsequent `log1p`, `pca`,
  `neighbors`, Leiden at resolution=1.0, UMAP are applied.
- **Lines 443-495:** Violin plots for top 10 positive genes by cluster. UMAP colored by
  top 4 genes. Optimal cluster per gene CSV.

### Differences

| Aspect | Notebook | Pipeline | Severity |
|--------|----------|----------|----------|
| Preprocessing parameters | `n_neighbors=8, n_pcs=0`, 3 Leiden resolutions (0.6, 1.4, 2.2), `min_dist=0.1` | Default `n_neighbors`, PCA-based, single Leiden at 1.0, default `min_dist` | MEDIUM |
| Cell filtering | `min_counts=10, min_genes=3` before clustering | No explicit filtering before clustering | LOW |
| Number of methods | 6 ST technologies | 1-2 (Reseg + optional Original) | Design choice |
| Cross-method violin | Same gene in its best cluster compared across all methods | Only within resegmented data | Design choice |
| Genes/cell and Counts/cell boxplots | Boxplots per method, log2 scale | Histogram overlay (KDE) | LOW |
| Layers | Saves `raw` layer for positivity from unprocessed counts | Uses QC metrics (`n_cells_by_counts`) | LOW |

### Severity: MEDIUM
The preprocessing parameter differences (especially `n_pcs=0` meaning no PCA in notebook
vs PCA in pipeline, and `n_neighbors=8` vs default 15) could lead to substantially
different clustering results. The single Leiden resolution (1.0) vs three resolutions
(0.6, 1.4, 2.2) limits exploratory analysis. However, for a single-sample pipeline run,
these differences are acceptable with appropriate documentation.

---

## Sub-step 4-5: Diffusion Analysis (Notebook 3_6)

### Notebook Intent
Measure transcript-to-centroid distance distributions for each technology to assess
signal diffusion/leakage. The notebook:
- Loads transcript-level data for all 6 technologies
- Computes Euclidean distance from each transcript to its assigned cell centroid
- Converts pixel distances to micrometers using technology-specific factors
- Plots complementary CDF (all genes, all methods overlaid)
- Plots per-gene complementary CDF in a 3x3 subplot grid for top 9 common genes
- Builds a gene-by-method heatmap of mean distances for top 25 genes
- Computes proportion of reads assigned to cells per method (stacked bar)

### Paper Methodology
Distance is computed as `sqrt((x - closest_cell_y)^2 + (y - closest_cell_x)^2)` -- note
the notebook has swapped x/y centroid columns (`closest_cell_y` for x-coordinate), which
is a known coordinate convention issue. Per-gene average distances are stored in a
`results` DataFrame for heatmap generation. The top 25 genes are selected by cross-method
presence (genes detected in most methods).

### Pipeline Implementation
- **Lines 498-567:** `analyze_diffusion()` with inner `calculate_distances()`. Computes
  per-transcript distance to centroid, handles multiple column naming conventions,
  aligns index types (str vs int).
- **Lines 569-581:** Saves global distance statistics (`describe()`) per dataset.
- **Lines 588-601:** Complementary CDF plot with 50k subsample for performance.
- **Lines 602-631:** Per-gene ECDF subplots for top 9 genes by transcript count.
- **Lines 633-660:** Gene x Method median distance heatmap for top 30 genes. Saves
  pivot table as CSV (`diffusion_gene_method_median_distances.csv`).

### Differences

| Aspect | Notebook | Pipeline | Severity |
|--------|----------|----------|----------|
| Per-gene summary CSV | `results` DataFrame with mean distance per gene per method, used for heatmap | Heatmap pivot CSV exists (line 659) but only for top 30 genes; no comprehensive per-gene summary with mean/median/std | CRITICAL |
| Number of top genes in heatmap | Top 25 (by cross-method presence) | Top 30 (by transcript count in combined data) | LOW |
| Distance statistic | Mean (for heatmap) | Median (for heatmap, line 639) | MEDIUM |
| Subsampling | 10% random sample at load time | 50k sample for CDF plot only; full data for heatmap | LOW |
| Assigned reads proportion | Stacked barplot of in-cell vs out-of-cell | Not implemented | LOW |
| Coordinate swap | Notebook swaps x/y centroids (convention issue) | Pipeline uses consistent `x_centroid`/`y_centroid` mapping | BUG FIX (pipeline is correct) |

### Severity: CRITICAL
The pipeline generates the heatmap pivot CSV (`diffusion_gene_method_median_distances.csv`)
at line 659 for the top 30 genes, which partially addresses the data export need.
However, a comprehensive per-gene diffusion summary CSV containing mean, median, and std
distance per gene per method for ALL genes (not just top 30) is missing. This full
dataset is needed for downstream analysis and reproducing paper figures that correlate
diffusion distances with other metrics (e.g., efficiency ratios). The pipeline should
export a complete `diffusion_per_gene_summary.csv` with columns:
`[gene, method, mean_distance_um, median_distance_um, std_distance_um, n_transcripts]`.

---

## Sub-step 4-6: Xenium vs Visium Comparison (Notebook 3_7)

### Notebook Intent
Compare Xenium single-cell resolution data against Visium (Fresh Frozen) spot-based data.
The notebook:
- Loads both Xenium and Visium FF datasets
- Assigns unannotated Xenium cells to nearest annotated region
- Computes tissue area using alpha shapes (`xb.calculating.alphashape_fun()`)
- Normalizes total counts by area (counts per cm^2) for fair comparison
- Generates per-gene scatter plots (Visium counts/cm^2 vs Xenium counts/cm^2) for
  Cortex, Hippocampus, and Thalamus separately
- Reports aggregate Visium/Xenium total count ratio (~13.7x in cortex and thalamus,
  ~12.9x in hippocampus)

### Pipeline Implementation
**Not implemented.** The pipeline has no Xenium vs Visium comparison.

### Differences

| Aspect | Notebook | Pipeline | Severity |
|--------|----------|----------|----------|
| Entire sub-step | Full area-normalized cross-technology comparison | Missing | MEDIUM |
| Alpha shape area computation | Uses `xb.calculating.alphashape_fun()` | N/A | N/A |
| Region-specific area ratios | Cortex, Hippocampus, Thalamus | N/A | N/A |

### Severity: MEDIUM
This analysis is specific to the original benchmarking study where both Xenium and Visium
data exist for the same tissue. For a general-purpose pipeline, this is optional --
it requires a paired Visium dataset that may not always be available. However, if the
pipeline aims to fully reproduce the paper, this sub-step should be added as an optional
analysis gated behind a config flag (e.g., `comparison.visium_path`).

---

## Cross-cutting: Inlined xb Functions (Code Quality)

### Inlined Code

| Function | Lines | Original Location | Differences |
|----------|-------|-------------------|-------------|
| `_coexpression_calculation()` | 36-44 | `xb/calculating.py:302` (`coexpression_calculation`) | Identical logic |
| `_negative_marker_purity_coexpression()` | 47-97 | `xb/calculating.py:196` (`negative_marker_purity_coexpression`) | Simplified: returns 1 or 3 values vs original's 6. Missing `lowvals_sc`, `lowvals_sp`, `commongenes` return. Missing `min_number_cells` filtering logic from original. |
| `PIXEL_TO_UM` dict | 99-107 | Hardcoded in notebook 3_6 cells | Equivalent values |

### Maintenance Concerns
- ~70 lines of duplicated code between pipeline and `xb/` library
- The inlined `_negative_marker_purity_coexpression()` diverges from the original
  by removing intermediate return values, creating a functional gap
- If `xb/calculating.py` is updated, the pipeline copy will not receive fixes
- The `_coexpression_calculation()` function uses a Python loop over columns with
  `tqdm`, which is very slow for large gene panels (O(n_genes^2))

### Recommendation
Import from `xb.calculating` directly instead of inlining. If self-containment is
required, add a fallback import pattern:
```python
try:
    from xb.calculating import negative_marker_purity_coexpression, coexpression_calculation
except ImportError:
    # Inline fallback
    ...
```

---

## Prioritized Fix List

### CRITICAL

1. **Export comprehensive per-gene diffusion summary CSV** (Sub-step 4-5)
   - File: `xenium_step4_techniques_comparison.py`, after line 586
   - Action: After building `df_all`, compute and export a full per-gene summary:
     ```python
     gene_summary = df_all.groupby(['Gene', 'Dataset'])['Distance_um'].agg(
         ['mean', 'median', 'std', 'count']
     ).reset_index()
     gene_summary.columns = ['gene', 'method', 'mean_distance_um',
                              'median_distance_um', 'std_distance_um', 'n_transcripts']
     gene_summary.to_csv(os.path.join(output_dir, 'diffusion_per_gene_summary.csv'), index=False)
     ```
   - Impact: Enables downstream cross-metric correlation analysis and paper figure
     reproduction

### HIGH

2. **Restore full NMP return values in inlined function** (Sub-step 4-3)
   - File: `xenium_step4_techniques_comparison.py`, lines 47-97
   - Action: Align `_negative_marker_purity_coexpression()` with the original
     `xb.calculating.negative_marker_purity_coexpression()` to return all 6 values
     when `pipeline_output=False`. At minimum, add `lowvals_sc`, `lowvals_sp`, and
     `commongenes` to the non-pipeline return path.
   - Impact: Enables per-gene purity analysis needed for efficiency-vs-specificity plot

3. **Add per-gene NMP breakdown and export** (Sub-step 4-3)
   - File: `xenium_step4_techniques_comparison.py`, lines 349-367
   - Action: Call `_negative_marker_purity_coexpression()` with `pipeline_output=False`
     to get per-gene purity scores. Save as CSV.
   - Impact: Provides granular specificity data for individual genes

### MEDIUM

4. **Add efficiency-vs-specificity scatter plot** (Sub-step 4-3)
   - File: `xenium_step4_techniques_comparison.py`, after line 397
   - Action: If both efficiency ratio CSV and per-gene NMP CSV are available, merge on
     gene name and generate a scatter plot (x=ratio, y=purity). This is a key paper
     figure.
   - Impact: Reproduces Figure showing the tradeoff between detection and specificity

5. **Align preprocessing parameters with notebook** (Sub-step 4-4)
   - File: `xenium_step4_techniques_comparison.py`, lines 434-441
   - Action: Make Leiden resolution, n_neighbors, n_pcs, and min_dist configurable
     via `comp_config`. Default to notebook values: `n_neighbors=8, n_pcs=0,
     resolution=2.2, min_dist=0.1`.
   - Impact: Clustering reproducibility

6. **Use mean instead of median for diffusion heatmap** (Sub-step 4-5)
   - File: `xenium_step4_techniques_comparison.py`, line 639
   - Action: Change `.median()` to `.mean()` to match notebook, or make it configurable.
   - Impact: Numerical consistency with paper figures

7. **Add gene filtering for NMP (ratio < 10)** (Sub-step 4-3)
   - File: `xenium_step4_techniques_comparison.py`, before line 360
   - Action: If efficiency ratio CSV is available, filter to genes with ratio < 10
     before computing NMP. This removes inflated genes that distort purity.
   - Impact: Methodological accuracy

### LOW

8. **Add Xenium vs Visium optional comparison** (Sub-step 4-6)
   - File: `xenium_step4_techniques_comparison.py`
   - Action: Add optional `analyze_visium_comparison()` gated behind
     `comp_config.get('visium_path')`. Requires alpha shape area computation.
   - Impact: Paper figure reproduction for paired Visium+Xenium datasets

9. **Replace inlined functions with imports** (Code quality)
   - File: `xenium_step4_techniques_comparison.py`, lines 34-107
   - Action: Use `from xb.calculating import ...` with fallback pattern
   - Impact: Maintenance, consistency with upstream fixes

10. **Add assigned-reads proportion barplot** (Sub-step 4-5)
    - File: `xenium_step4_techniques_comparison.py`, in `analyze_diffusion()`
    - Action: Compute fraction of transcripts with valid cell assignment per dataset.
      Plot stacked barplot.
    - Impact: Minor diagnostic visualization

11. **Make cell filtering explicit before clustering** (Sub-step 4-4)
    - File: `xenium_step4_techniques_comparison.py`, line 435
    - Action: Add `sc.pp.filter_cells(ad_proc, min_counts=10)` and
      `sc.pp.filter_cells(ad_proc, min_genes=3)` before normalization.
    - Impact: Alignment with notebook methodology

---

## Summary Table

| Sub-step | Notebook | Pipeline Status | Max Severity |
|----------|----------|----------------|--------------|
| 4-2 Efficiency | 3_3 | Implemented (simplified) | MEDIUM |
| 4-3 Specificity | 3_4 | Partially implemented | HIGH |
| 4-4 Positivity | 3_5 | Well implemented | MEDIUM |
| 4-5 Diffusion | 3_6 | Mostly implemented | CRITICAL |
| 4-6 Xenium vs Visium | 3_7 | Not implemented | MEDIUM |

**Total issues:** 11 (1 CRITICAL, 2 HIGH, 4 MEDIUM, 4 LOW)
