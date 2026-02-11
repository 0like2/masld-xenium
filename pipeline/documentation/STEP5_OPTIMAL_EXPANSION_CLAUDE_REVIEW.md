# Step 5 -- Optimal Expansion: Notebook vs Pipeline Review

**Pipeline file:** `pipeline/xenium_step5_optimal_expansion.py` (612 lines)
**Notebook:** `notebooks/4_optimal_expansion/4_1_Optimal_expansion_multisection.ipynb`
**Library dependency:** `xb/calculating.py` (`dist_nuc`, `dispersion`, `distance_calc`)
**Review date:** 2026-02-11

---

## Overview

Step 5 determines the **optimal expansion radius** from the nucleus border for each
cell type, using a correlation-based turnover analysis.  The algorithm compares
the expression profile at increasing distances from the cell centroid to (a) the
nuclear expression profile and (b) the domain background profile.  Where the
nuclear correlation drops below the background correlation (diff < threshold),
the "turnover distance" is declared.  Subtracting the mean nuclei size yields the
optimal expansion from the nucleus boundary.

The pipeline reimplements notebook `4_1_Optimal_expansion_multisection.ipynb`
across eight logical substeps (5-1 through 5-8e).

---

## 5-1. Load Original Transcripts (Step 0)

### Notebook Intent
Load the raw Xenium `.h5ad` file and extract the reads table from `adata.uns['spots']`.
The notebook reads `ms_brain_multisection1.h5ad` directly.

### Paper Methodology
Xenium output transcripts are the starting point; every read carries `x_location`,
`y_location`, `cell_id`, and `feature_name`.

### Pipeline Implementation
- Lines 331--352 of `xenium_step5_optimal_expansion.py`.
- Searches for the Step 0 output file at `{output_dir}/{sample_tag}.h5ad`.
- Falls back to `{parent_dir}/step0_formatting/{sample_tag}.h5ad` if the
  primary path does not exist.
- Extracts `adata_step0.uns['spots']` and copies it into `reads_original`.

### Differences
| Aspect | Notebook | Pipeline |
|--------|----------|----------|
| File path | Hard-coded relative path | Config-driven with fallback |
| Missing spots | Would crash | Explicit error + early return |

### Severity: LOW
Functionally equivalent.  The pipeline adds robustness for missing data.

---

## 5-2. Load Annotated Cells (Domain Source)

### Notebook Intent
Load a separately processed AnnData (`adata_multisection_nuclei_r1_with_annotations.h5ad`)
containing cell-level annotations: `Class` (cell type) and `spatial_annotation` (domain).

### Paper Methodology
Cell-type and domain labels come from upstream annotation (clustering + spatial
domain labeling).

### Pipeline Implementation
- Lines 354--379.
- First tries `config['previous_step_adata_path']`.
- If missing, falls back to searching sibling step directories:
  `step4_resegmentation`, `step1_exploration`, `step2_segmentation_free`.
- Loads `adata_annotated` and reports shape.

### Differences
| Aspect | Notebook | Pipeline |
|--------|----------|----------|
| Input path | Single hard-coded file | Multi-path fallback |
| Fallback naming | N/A | Uses `step4_resegmentation` -- actual directory is `step3_resegmentation` |

### Severity: MEDIUM
The fallback path `step4_resegmentation` does not match the pipeline's actual
output directory name `step3_resegmentation`.  If `previous_step_adata_path` is
unset AND the step3 output directory has the correct name, the fallback will
silently miss it.

**Recommendation:** Change `step4_resegmentation` to `step3_resegmentation` in
the fallback search list (line 364).

---

## 5-3. Map Domain Assignments to Reads

### Notebook Intent
Map `cell_id -> spatial_annotation` (domain) and `cell_id -> Class`
(initial_annotation) onto reads.  Split reads into `annotatedcells` (have domain)
and `nancells` (domain is NaN).

```python
# Notebook (cells 6, 8)
splocdic = dict(zip(adata1.obs['cell_id'], adata1.obs['spatial_annotation']))
reads_original['domain'] = reads_original['cell_id'].map(splocdic)
```

### Paper Methodology
Each transcript is associated with its originating cell's domain to enable
domain-specific background comparison.

### Pipeline Implementation
- Lines 381--426.
- Searches `adata_annotated.obs` for a domain key using priority list:
  `['spatial_annotation', 'Class', 'leiden', 'cluster', 'graph_clusters']`.
- Maps both `domain` and `initial_annotation` from the same key onto reads.
- Separates `annotatedcells` / `nancells`.

### Differences
| Aspect | Notebook | Pipeline |
|--------|----------|----------|
| Domain key | Always `spatial_annotation` | Priority-based search |
| Annotation key | Always `Class` | Same key as domain (line 413) |
| `region_annotation` | N/A | Not in priority list |

### Severity: MEDIUM
Two concerns:

1. **`initial_annotation` uses the same key as `domain`** (line 413 duplicates
   the domain mapping into `initial_annotation`).  The notebook uses `Class` for
   cell-type annotation and `spatial_annotation` for domain -- these are
   different columns.  In the pipeline, if `spatial_annotation` is found first as
   the domain key, then `initial_annotation` will also be set to
   `spatial_annotation` values rather than cell-type labels.

2. **`region_annotation`** (used in `xb.calculating.domainassign()`) is absent
   from the priority list.

**Recommendation:** Separate the domain key search from the cell-type annotation
key search.  Use a distinct priority list for cell-type annotation
(e.g., `['Class', 'initial_annotation', 'celltype', 'cell_type']`).

---

## 5-4. KDTree Nearest-Domain Expansion

### Notebook Intent
Subsample 1% of annotated reads as anchors, build a cKDTree, and query every
unassigned read to find the nearest anchor's domain.

```python
# Notebook (cells 10, 12)
sub = rd.sample(list(annotatedcells.index), int(np.round(int(annotatedcells.shape[0]) * 0.01)))
annotatedsub = annotatedcells.loc[sub, :]
tree1 = cKDTree(coords1)
_, index = tree1.query(coord)
```

### Paper Methodology
Unassigned transcripts are assigned to the nearest segmented cell's domain via
spatial proximity, enabling domain-level background expression estimation.

### Pipeline Implementation
- Lines 428--489.
- Subsamples using `subsample_fraction` (default 0.01, configurable).
- Applies minimum anchor count guard: `n_sample = max(n_sample, min(n_assigned, 1000))`.
- Queries in chunks of 1,000,000 for memory efficiency.
- Optionally applies `distance_threshold` to reject reads too far from any anchor.
- Transfers both `domain` and `initial_annotation` from the anchor to the
  unassigned read.

### Differences
| Aspect | Notebook | Pipeline |
|--------|----------|----------|
| Subsampling | `random.sample`, no seed | `DataFrame.sample(random_state=42)` -- reproducible |
| Chunked query | No (iterates one point at a time) | Yes, 1M-point chunks -- much faster |
| Distance threshold | None | Optional config parameter |
| Min anchor count | None | Floor of 1000 |

### Severity: LOW
The pipeline is a strict improvement: reproducible seeding, vectorized querying,
and optional guardrails.  Core algorithm is identical.

---

## 5-5. Distance Threshold Filter (Optional)

### Notebook Intent
No distance threshold is applied in the notebook.

### Paper Methodology
The paper does not mention a hard distance cutoff; the turnover analysis itself
serves as the quality criterion.

### Pipeline Implementation
- Lines 466--479.
- If `config['optimal_expansion']['distance_threshold']` is set, reads assigned
  beyond that distance are set to `NaN` domain.
- Default is `None` (no filtering).

### Differences
This is a pipeline-only addition (operational guardrail).

### Severity: NONE
Does not affect results when disabled (default).

---

## 5-6. Spatial Map + Distance Histogram

### Notebook Intent
The notebook saves a scatter plot of reads colored by domain but does not
generate a dedicated distance histogram.

### Pipeline Implementation
- Lines 492--524.
- Generates two plots:
  1. Spatial scatter of 100k sampled reads colored by domain.
  2. Histogram of distances from unassigned reads to their nearest anchor.
- Both saved as PNG.

### Differences
| Aspect | Notebook | Pipeline |
|--------|----------|----------|
| Spatial plot | Interactive `plt.show()` | Saved PNG |
| Distance histogram | Not present | Added for QC |
| Subsampling for plot | Full dataset | 100k sample |

### Severity: NONE
Visualization-only; does not affect computed outputs.

---

## 5-7. Save Expanded Transcripts

### Notebook Intent
Saves the fully domain-assigned reads table to CSV.

```python
# Notebook (cell 14)
reads_original.to_csv('.../reads_multisection_r1_with_domain.csv')
```

### Pipeline Implementation
- Lines 526--529.
- Saves to `{output_dir}/{sample_tag}_step5_expanded_transcripts.csv`.

### Differences
Output naming convention only.

### Severity: NONE

---

## 5-8. Turnover / Crossover Analysis

This is the core algorithm of Step 5.  It is implemented in
`calculate_turnover()` (lines 73--306) and invoked from `run_step5()` at
lines 531--601.

---

### 5-8 (pre). Distance Column Computation

### Notebook Intent
Compute per-read distance from transcript position to its cell centroid:
```python
# Notebook (cell 23)
dist2 = (reads_assigned['x_location'] - reads_assigned['x_cell'])**2 + \
        (reads_assigned['y_location'] - reads_assigned['y_cell'])**2
reads_assigned['distance'] = np.sqrt(dist2)
reads_assigned['distance'] = reads_assigned['distance'].round(0)
```
This is equivalent to `xb.calculating.dispersion()`, which computes the same
Euclidean distance and stores it in the `distance` column.

### Pipeline Implementation
- Lines 552--582.
- Maps `x_centroid` / `y_centroid` from `adata_annotated.obs` (or
  `adata_step0.obs` as fallback) onto reads via `cell_id`.
- Computes Euclidean distance and rounds to integer.
- Guards against missing centroid data (line 563: early return with warning).

### Differences
| Aspect | Notebook | Pipeline |
|--------|----------|----------|
| Centroid source | `adata1.obs` directly | Two-level fallback (annotated, then step0) |
| Rounding | `round(0)` | `round(0)` -- identical |
| Missing centroid | Would crash | Warns and skips turnover |

### Severity: LOW
Functionally identical.  Pipeline adds resilience.

---

### 5-8a. Nuclear vs Background Expression Profiles

### Notebook Intent
- **Background:** `pd.crosstab(reads_not_assigned['domain'], reads_not_assigned['feature_name'])`
  gives a domain x gene count matrix from unassigned reads.
- **Nuclear:** `pd.crosstab(reads_ctd['overlaps_nucleus'], reads_ctd['feature_name'])`
  gives a 2-row (0/1) x gene matrix; row `1` = nuclear reads expression profile.

```python
# Notebook (cell 25-26)
background_express = pd.crosstab(reads_not_assigned['domain'], reads_not_assigned['feature_name'])
reads_ctd_nucl = pd.crosstab(reads_ctd['overlaps_nucleus'], reads_ctd['feature_name'])
```

### Paper Methodology
The nuclear expression profile serves as the "true" cell-type signature.  The
background profile represents ambient/leaked reads.  Comparing them at increasing
distances reveals where cell-specific signal fades into noise.

### Pipeline Implementation
- Lines 99--109 (background crosstab).
- Lines 156--158 (nuclear crosstab per domain).
- Lines 90--97: if `overlaps_nucleus` is missing, creates a proxy using the
  median distance as a threshold (reads closer than median = "nuclear").

### Differences
| Aspect | Notebook | Pipeline |
|--------|----------|----------|
| `overlaps_nucleus` | Always present | Proxy if missing (median distance) |
| Background reads | `cell_id == -1` | Flexible dtype handling (int/str) |

### Severity: LOW
The `overlaps_nucleus` proxy is a reasonable fallback.  The notebook assumes the
column exists because it uses 10X Xenium output directly, which includes it.
Pipeline datasets from other upstream steps may lack it.

**Recommendation:** Log a warning when the proxy is used, so users know the
nuclear profile is approximate.  (Already done at lines 93--97.)

---

### 5-8b. Distance-Bin Correlation Curves

### Notebook Intent
For each (cell-type, domain) pair with >5000 reads:
1. Build a distance x gene crosstab (`expression_distances`).
2. At each distance bin, compute Pearson correlation with the nuclear profile
   (`corr_nuc`) and the background profile (`corr_bck`).
3. Store in a summary DataFrame with columns `corr_nuc`, `corr_back`, `diff`.

```python
# Notebook (cell 26, inner loop)
for dist in expression_distances.index:
    corr_bck.append(np.corrcoef(expression_distances.loc[dist,:], bck_sub.loc[dom,:])[0,1])
    corr_nuc.append(np.corrcoef(expression_distances.loc[dist,:], reads_ctd_nucl.loc[1,:])[0,1])
```

### Paper Methodology
Pearson correlation at each 1-um distance bin captures how expression profile
transitions from cell-type-specific (nuclear) to ambient (background).

### Pipeline Implementation
- Lines 150--206.
- Identical algorithm structure: distance x gene crosstab, per-bin correlation.
- **Improvement:** pipeline aligns columns across all three matrices using set
  intersection (`common_genes`) before computing correlations (lines 168--179).
  The notebook uses `isin()` filtering and `sort_values()` which achieves the
  same goal less robustly.
- **Improvement:** pipeline guards against zero-variance vectors (`np.std == 0`)
  to avoid NaN correlations (lines 189--196).  The notebook does not.

### Differences
| Aspect | Notebook | Pipeline |
|--------|----------|----------|
| Column alignment | `isin()` + `sort_values()` | Set intersection + sorted list |
| Zero-variance guard | None | `np.std == 0` check, appends NaN |
| Min genes threshold | None | `len(common_genes) < 5` skips pair |

### Severity: NONE
Pipeline is more robust.  Core correlation algorithm is identical.

---

### 5-8c. Turnover Distance Detection

### Notebook Intent
Find the **minimum distance** where `diff = corr_nuc - corr_back < 0.1`:
```python
# Notebook (cell 26)
try:
    tdistance = np.min(summary.loc[summary['diff'] < 0.1, :].index)
except:
    tdistance = np.max(summary.index)
if str(np.min(summary.loc[summary['diff'] < 0.1, :].index)) == 'nan':
    tdistance = np.max(summary.index)
```
The notebook uses a bare `except` and a string comparison for NaN detection.

### Paper Methodology
The crossover point where background correlation exceeds nuclear correlation
marks the boundary of meaningful cell-type signal.

### Pipeline Implementation
- Lines 225--235.
- Uses `np.nanmin()` instead of `np.min()` for robustness.
- Falls back to `np.nanmax(summary.index)` if no bin meets the threshold.
- Configurable `diff_threshold` (default 0.1, matching notebook).

### Differences
| Aspect | Notebook | Pipeline |
|--------|----------|----------|
| Threshold | Hard-coded 0.1 | Configurable (`diff_threshold`) |
| NaN handling | String comparison (`'nan'`) | `np.isnan()` check |
| Fallback | `np.max` | `np.nanmax` |

### Severity: NONE
Pipeline is strictly more robust.  Same algorithm, better edge-case handling.

---

### 5-8d. Nuclei Size and Cell Size via `dist_nuc()`

### Notebook Intent
After computing the turnover distance for a (cell-type, domain) pair, compute
`nuclei_size` and `cell_size` using `dist_nuc()`:
```python
# Notebook (cell 26)
reads_ctdsub = reads_ctd[reads_ctd['overlaps_nucleus'] == 1]
cellm.append(dist_nuc(reads_ctd))        # all reads -> "cell size"
nuclim.append(dist_nuc(reads_ctdsub))    # nuclear reads only -> "nuclei size"
```

### Original `dist_nuc()` (xb/calculating.py, lines 23--44)
```python
def dist_nuc(reads_ctdsub):
    allds = []
    for g, n in reads_ctdsub.groupby('cell_id'):
        try:
            hull = ConvexHull(np.array(n.loc[:, ['x_location', 'y_location']]))
            allds.append(np.mean(n.iloc[hull.vertices]['distance']))
        except:
            print()
    median_dist = np.median(allds)
    return median_dist
```

**What this function actually computes:**
1. Group reads by `cell_id`.
2. For each cell, compute the ConvexHull of its read positions.
3. Take the mean of the `distance` column (= distance to centroid) at the hull
   vertices.  This is the mean centroid-distance of the outermost reads.
4. Return the median across all cells.

So `dist_nuc()` computes the **median (across cells) of the mean hull-vertex
distance-to-centroid**.  When called on nuclear reads only, this estimates the
nuclear radius; when called on all reads, it estimates the cell radius.

### Pipeline `dist_nuc()` (lines 36--52)
```python
def dist_nuc(reads_ctdsub):
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
```

### Differences
| Aspect | Original (`xb`) | Pipeline |
|--------|-----------------|----------|
| `len(n) < 3` guard | None (ConvexHull would raise) | Explicit skip |
| `distance` column check | None (assumes present) | `if 'distance' in n.columns` |
| Exception handling | Bare `except: print()` | `except Exception: pass` |
| Empty `allds` | `np.median([])` = warning + NaN | Explicit `np.nan` return |

### Severity: MEDIUM

The pipeline version is functionally faithful but has a subtle issue:

- The `if 'distance' in n.columns` guard means that if the `distance` column is
  missing, the hull is still computed but **nothing is appended to `allds`**.
  The function silently returns `NaN` without any warning.  This could mask data
  preparation bugs.

- The `distance` column **must** be computed before calling `calculate_turnover()`.
  The pipeline does this at lines 578--582, so under normal execution the column
  exists.  However, if `calculate_turnover()` is ever called from a different
  context without the distance column, results will be silently empty.

**Recommendation:** Add a warning log when `distance` column is missing in
`dist_nuc()` to make debugging easier.

---

### 5-8d (cont). Optimal Expansion Formula

### Notebook Intent
```python
# Notebook (cells 28, 30, 32)
optimal_expansion = np.nanmean(meand_celltype) - np.nanmean(nuclimall)
# Example output: 10.71 - 5.06 = 5.65
```

### Paper Methodology
Optimal expansion from the nucleus border = mean turnover distance (from
centroid) minus mean nuclei radius (estimated via `dist_nuc` on nuclear reads).

### Pipeline Implementation
- Line 256: `optimal_expansion = np.nanmean(meand_celltype) - np.nanmean(nuclimall)`
- Identical formula.

### Differences
| Aspect | Notebook | Pipeline |
|--------|----------|----------|
| Formula | `mean(turnover) - mean(nuclei_size)` | Identical |
| Aggregation | `np.mean` (no NaN handling) | `np.nanmean` (NaN-safe) |

### Severity: NONE
Exact match.  Pipeline version is more NaN-robust.

---

### 5-8e. Summary Barplot + CSV Outputs

### Notebook Intent
- Create a horizontal barplot showing turnover per cell type, with scatter
  markers for `nuclei_size` (pink) and `cell_size` (black).
- Uses cell-type-specific colors from `adata.uns['Class_colors']`.
- Also creates a filtered version showing only cell types with scores from >5
  domains.

```python
# Notebook (cell 36)
sns.barplot(data=tfmerge, y=tfmerge.cluster, x='score', palette=tf.colors, ...)
sns.scatterplot(data=tfmerge, y=tfmerge.cluster, x='nuclei_size', color='#D83066', ...)
sns.scatterplot(data=tfmerge, y=tfmerge.cluster, x='cell_size', ...)
```

### Paper Methodology
Figure 3-style visualization showing per-cell-type turnover distances with
nuclear and cell boundary markers.

### Pipeline Implementation
- Lines 262--304.
- Generates a barplot with turnover bars + scatter markers for `nuclei_size` and
  `cell_size`.
- Uses `tab20` palette instead of custom `Class_colors`.
- Saves three output files:
  1. `{sample_tag}_step5_turnover_summary.csv` -- domain x cell-type turnover matrix.
  2. `{sample_tag}_step5_turnover_per_celltype.csv` -- per-cell-type summary.
  3. `{sample_tag}_step5_optimal_expansion.txt` -- scalar result with
     `optimal_expansion`, `mean_turnover`, `mean_nuclei_size`.

### Differences
| Aspect | Notebook | Pipeline |
|--------|----------|----------|
| Color palette | `Class_colors` from adata | `tab20` generic |
| Domain-count filter | Shows subset with >5 domains | Not implemented |
| Output format | PDF figures, no CSV | PNG figure + 3 CSV/TXT files |
| Per-domain scores | Merged into barplot (`tfmerge`) | Not merged; simpler per-celltype view |

### Severity: LOW
The pipeline simplifies the visualization (no per-domain overlay, generic palette)
but captures all quantitative outputs.  The >5 domain filter is a presentation
choice and does not affect the computed optimal expansion value.

---

## Cross-Cutting Concerns

### A. `initial_annotation` vs `domain` Key Confusion (Lines 401--414)

In the notebook, `initial_annotation` comes from `Class` (cell type) and `domain`
comes from `spatial_annotation`.  These are semantically distinct columns.

In the pipeline, both are mapped from the **same** domain key (lines 401 and 413):
```python
domain_map = dict(zip(annotated_ids, adata_annotated.obs[domain_key]))
ct_map = dict(zip(annotated_ids, adata_annotated.obs[domain_key]))  # same key!
```

This means if `spatial_annotation` is found, both `domain` and
`initial_annotation` will contain spatial domain labels instead of cell-type
labels.  The turnover analysis iterates over `initial_annotation` as cell types
(line 129), which would then be iterating over domains -- collapsing the
cell-type x domain matrix into a domain x domain analysis.

**Severity: HIGH** -- This is the most significant discrepancy.  The analysis
loses its cell-type resolution when `spatial_annotation` is found first.

**Recommendation:** Use two separate key searches:
```python
# Domain key
domain_key = find_key(adata_annotated.obs, ['spatial_annotation', 'region_annotation', 'domain'])
# Cell-type key
celltype_key = find_key(adata_annotated.obs, ['Class', 'initial_annotation', 'celltype', 'cell_type', 'leiden'])
```

### B. Reads Split Logic (Lines 539--547)

The notebook splits reads simply:
```python
reads_not_assigned = reads_original[reads_original['cell_id'] == -1]
reads_assigned = reads_original[reads_original['cell_id'] != -1]
```

The pipeline handles mixed dtypes (int vs string `cell_id` values), which is
correct for robustness.  However, the split happens **after** the KDTree
expansion (line 534), meaning the "not assigned" reads here refer to reads with
`cell_id == -1` (never belonging to any cell), not "reads without domain."  This
is consistent with the notebook.

**Severity: NONE** -- Correct behavior.

### C. Pixel-to-Micron Conversion

Neither the notebook nor the pipeline explicitly converts pixel coordinates to
micrometers.  The `distance` column and all size metrics are in the native
coordinate system of the Xenium data.  For Xenium, the conversion factor is
~4.70588 pixels/um (per MEMORY.md).  The paper reports optimal expansion in
micrometers (~5.6 um), and the notebook's raw output (5.65) appears to already be
in appropriate units, suggesting the Xenium spatial coordinates used in this
dataset are in micrometers.

**Severity: NONE** -- Consistent between notebook and pipeline.

---

## Summary Table

| Substep | Description | Match Quality | Severity |
|---------|-------------|---------------|----------|
| 5-1 | Load original transcripts | Faithful | LOW |
| 5-2 | Load annotated cells | Mostly faithful | MEDIUM (fallback path) |
| 5-3 | Map domain/annotation to reads | **Divergent** | **HIGH** (key confusion) |
| 5-4 | KDTree expansion | Improved | LOW |
| 5-5 | Distance threshold | Pipeline-only | NONE |
| 5-6 | Visualization | Simplified | NONE |
| 5-7 | Save expanded transcripts | Faithful | NONE |
| 5-8 (pre) | Distance column | Faithful | LOW |
| 5-8a | Nuclear/background profiles | Faithful + proxy | LOW |
| 5-8b | Distance-bin correlations | Improved | NONE |
| 5-8c | Turnover detection | Improved | NONE |
| 5-8d | `dist_nuc()` + formula | Faithful | MEDIUM (silent NaN) |
| 5-8e | Summary barplot + CSVs | Simplified | LOW |

---

## Prioritized Fix List

### Priority 1 -- HIGH: `initial_annotation` / `domain` Key Separation (Lines 401--414)

**Problem:** Both `domain` and `initial_annotation` are mapped from the same
`domain_key`, causing the turnover analysis to lose cell-type resolution when
`spatial_annotation` is found.

**Fix:** Introduce a separate cell-type key search:
```python
celltype_key = None
celltype_priority = ['Class', 'initial_annotation', 'celltype', 'cell_type']
for key in celltype_priority:
    if key in adata_annotated.obs.columns:
        celltype_key = key
        break
if celltype_key is None:
    celltype_key = domain_key  # fallback to domain if no cell type found
    logger.warning("No cell-type column found; using domain key as cell-type proxy.")
ct_map = dict(zip(annotated_ids, adata_annotated.obs[celltype_key]))
```

**Impact:** Restores the intended cell-type x domain analysis structure.

---

### Priority 2 -- MEDIUM: Fallback Path Name Mismatch (Line 364)

**Problem:** Fallback path uses `step4_resegmentation` but the pipeline step is
named `step3_resegmentation`.

**Fix:** Change line 364:
```python
os.path.join(parent_dir, "step3_resegmentation", f"{sample_tag}_resegmented.h5ad"),
```

**Impact:** Enables fallback to find resegmentation output when
`previous_step_adata_path` is not set.

---

### Priority 3 -- MEDIUM: `dist_nuc()` Silent NaN on Missing Distance Column (Lines 44--47)

**Problem:** If `distance` column is absent, `dist_nuc()` returns `NaN` without
any diagnostic output, making data preparation bugs hard to detect.

**Fix:** Add a warning:
```python
if 'distance' not in n.columns:
    logger.warning("dist_nuc(): 'distance' column missing in group %s", g)
    continue
```

**Impact:** Improves debuggability without changing behavior.

---

### Priority 4 -- LOW: Domain Key Priority List (Line 385)

**Problem:** `region_annotation` (used by `xb.calculating.domainassign()`) is not
in the priority search list.

**Fix:** Add `'region_annotation'` to the priority list:
```python
priority_keys = ['spatial_annotation', 'region_annotation', 'Class', 'leiden', 'cluster', 'graph_clusters']
```

**Impact:** Better compatibility with datasets using `domainassign()` output.

---

### Priority 5 -- LOW: Visualization Color Palette (Line 269)

**Problem:** Pipeline uses generic `tab20` instead of dataset-specific
`Class_colors`.

**Fix (optional):** Check for `adata_annotated.uns.get(f'{celltype_key}_colors')`
and use it if available:
```python
palette = adata_annotated.uns.get(f'{celltype_key}_colors', 'tab20')
```

**Impact:** Visual consistency with notebook figures.

---

## Overall Assessment

Step 5 is **well-implemented** and closely follows the notebook's core algorithm.
The correlation-based turnover analysis, KDTree expansion, and optimal expansion
formula are all faithfully reproduced with added robustness (NaN guards, chunked
processing, configurable thresholds, proxy for `overlaps_nucleus`).

The single HIGH-severity issue -- using the same key for both domain and cell-type
annotation -- is a structural bug that could produce incorrect turnover matrices
when `spatial_annotation` exists in the input data.  This should be fixed before
production use.

All other issues are MEDIUM or LOW severity and relate to operational resilience
(fallback paths, diagnostic logging) rather than algorithmic correctness.
