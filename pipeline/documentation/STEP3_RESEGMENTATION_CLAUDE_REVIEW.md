# Step 3: Resegmentation -- Notebook vs Pipeline Implementation Review

**Pipeline file:** `pipeline/xenium_step3_resegmentation.py` (241 lines)
**Notebook:** `notebooks/3_techniques_comparison/3_1_resegmentation_notebooks/batch_segmentation_cellpose-Xenium.ipynb`
**xb library:** `xb/calculating.py`
**Config section:** `resegmentation` in `pipeline/config.yaml` (lines 74-86)

---

## Section 3-1: Device Detection (CUDA/MPS/CPU)

### Notebook Intent
The notebook uses `models.Cellpose(gpu=False, model_type='nuclei')` with a hardcoded `gpu=False` flag. No device auto-detection is performed; all segmentation runs on CPU.

### Paper Methodology
The paper does not specify a particular compute device requirement. Cellpose segmentation should produce identical masks regardless of device.

### Pipeline Implementation
Lines 31-41: `_detect_device()` auto-detects CUDA > MPS > CPU and returns a `torch.device` object.
Line 79: `device, use_gpu = _detect_device()`
Line 81: `model = models.CellposeModel(gpu=use_gpu, model_type='nuclei', device=device)`

### Specific Differences
1. The pipeline uses `models.CellposeModel` (the newer API), while the notebook uses `models.Cellpose` (legacy API). Both produce equivalent results for nuclei segmentation.
2. The pipeline enables GPU acceleration when available, which is an improvement over the notebook's CPU-only approach.
3. No functional divergence in segmentation output.

### Severity Rating
**NONE** -- This is a valid improvement. No action required.

---

## Section 3-2: Load DAPI Image

### Notebook Intent
Loads DAPI image from a per-dataset directory using `tifffile.imread()`. Handles both `.tif` and `.tiff` extensions via try/except fallback. For large images (Xenium), the notebook also demonstrates a tile-based approach (cells `a46c4ea3` through `13340b8e`), segmenting in 2000x2000 tiles and recomposing into a full mask.

### Paper Methodology
Standard DAPI fluorescence image loading as the input for nuclear segmentation.

### Pipeline Implementation
Lines 59-71: Loads DAPI image from `dapi_image_path` using `tifffile.TiffFile` context manager, then `tf.imread()`. Prints shape.

### Specific Differences
1. The pipeline opens a `TiffFile` context manager (line 65) but then calls `tf.imread()` separately (line 66), making the context manager redundant. Minor inefficiency, not a bug.
2. The notebook's tile-based recomposition approach (for very large Xenium images) is not implemented in the pipeline. Instead, the pipeline relies on Cellpose's internal `tile=True` parameter (line 86) when the image exceeds 2000px. This is functionally equivalent and simpler.

### Severity Rating
**LOW** -- The redundant `TiffFile` context manager is a cosmetic issue. Cellpose internal tiling is an acceptable substitute for manual tiling.

---

## Section 3-3: Cellpose Nuclei Segmentation

### Notebook Intent
Runs `model.eval(dapi_image, diameter=None, channels=[0, 0])` to segment nuclei from DAPI. Returns `masks, flows, styles, diams`. Then applies `label()` to convert Cellpose masks to labeled integer masks:
```python
labeled_nuclei = label(masks)
```

### Paper Methodology
Cellpose with the `nuclei` model is applied to DAPI images. The `diameter=None` setting triggers Cellpose's auto-diameter estimation. This is a standard approach for nuclear segmentation in spatial transcriptomics.

### Pipeline Implementation
Lines 75-91: Loads Cellpose parameters from config, auto-detects device, runs `model.eval()` with configurable `diameter` and `channels=[0,0]`. Uses `tile=True` for large images. Cleans up DAPI image and intermediate arrays with `del` and `gc.collect()`.

### Specific Differences
1. **CRITICAL: Missing `label()` call.** The notebook applies `label(masks)` to convert the Cellpose output into properly labeled connected components. The pipeline uses raw Cellpose `masks` directly (line 96, 118, 136). In practice, Cellpose `masks` are already integer-labeled (each cell gets a unique integer), so `label()` has minimal effect unless mask IDs are non-contiguous. However, the notebook's subsequent `expand_labels()` depends on `labeled_nuclei`, not raw masks. This omission is part of the larger expand_labels gap described in Section 3-5.
2. Memory management (`del`, `gc.collect()`) is an improvement in the pipeline.

### Severity Rating
**LOW** -- The `label()` call is mostly redundant for Cellpose output, but its absence contributes to the missing `expand_labels` pipeline described below.

---

## Section 3-4: Save Masks to TIF

### Notebook Intent
Saves the full mask (after tile recomposition or direct segmentation) to a TIF file:
```python
tf.imwrite(path_unsegmented_datasets + ds + "segmented_mask.tif", full_mask)
```
Also visualizes the expanded nuclei using `label2rgb`:
```python
color1 = label2rgb(expanded_nuclei, bg_label=0)
plt.imshow(color1)
```

### Paper Methodology
Intermediate mask storage for reproducibility and downstream use.

### Pipeline Implementation
Lines 93-96: Saves Cellpose masks to a TIF file with a `{sample_tag}_step3_resegmented_masks.tif` naming convention.

### Specific Differences
1. The pipeline does not save the expanded mask (because it never creates one -- see Section 3-5).
2. No visualization of the mask is generated. This is acceptable for a batch pipeline, but diagnostic QC plots would be beneficial.

### Severity Rating
**LOW** -- Saving only nuclei masks is consistent with the pipeline's current (incomplete) implementation. Once expand_labels is added, the expanded mask should also be saved or at least used downstream.

---

## Section 3-5: Map Transcripts to Cell Masks

### Notebook Intent
This is the core assignment step. The notebook performs **dual assignment** using both nuclei masks and expanded masks:

```python
labeled_nuclei = label(masks)
expanded_nuclei = expand_labels(labeled_nuclei, distance=400)

# Build centroid dictionaries from nuclei (not expanded) regionprops
centroid_dictx = {}
centroid_dicty = {}
for nucleus_props in regionprops(labeled_nuclei):
    centroid_dictx[nucleus_props.label] = nucleus_props.centroid[0]
    centroid_dicty[nucleus_props.label] = nucleus_props.centroid[1]

# Dual assignment: each transcript gets both labels
for ind in tqdm(read_positions.index):
    try:
        in_cell.append(labeled_nuclei[int(...), int(...)])
        closest_cell.append(expanded_nuclei[int(...), int(...)])
    except:
        in_cell.append(0)
        closest_cell.append(0)

read_positions['in_cell'] = in_cell
read_positions['closest_cell'] = closest_cell
```

Key aspects:
- `in_cell`: transcript falls directly within a nucleus mask (nuclear transcripts)
- `closest_cell`: transcript falls within the expanded mask (cytoplasmic + nuclear transcripts)
- `expand_labels(labeled_nuclei, distance=400)` grows each nucleus region by 400 pixels in all directions, capturing cytoplasmic transcripts that lie outside the nucleus but belong to the same cell
- Both columns are stored in the transcript DataFrame

### Paper Methodology
The paper explicitly discusses that nuclear-only segmentation misses cytoplasmic transcripts. The expansion approach (similar to 10x Genomics' default cell expansion) captures transcripts in the perinuclear/cytoplasmic space. The dual `in_cell`/`closest_cell` tracking is essential for:
1. Measuring transcript assignment efficiency (what fraction of transcripts fall in nuclei vs expanded cells)
2. Diffusion analysis in Step 4 (comparing nuclear vs cytoplasmic transcript profiles)
3. Computing distance-to-centroid distributions that characterize read dispersion

### Pipeline Implementation
Lines 98-138:
```python
# Line 107: um_per_pixel loaded but not used for expansion
um_per_pixel = reseg_config.get('um_per_pixel', 1.0)

# Lines 118-134: Direct coordinate-to-mask lookup
max_h, max_w = masks.shape
y_coords = df_transcripts[y_col].values.astype(int)
x_coords = df_transcripts[x_col].values.astype(int)
valid_mask = (y_coords >= 0) & (y_coords < max_h) & (x_coords >= 0) & (x_coords < max_w)

# Line 136: Single assignment to nuclei masks only
cell_labels = masks[y_coords, x_coords]
df_valid['cell_id_reseg'] = cell_labels
df_assigned = df_valid[df_valid['cell_id_reseg'] > 0].copy()
```

### Specific Differences

1. **CRITICAL: `expand_labels` is imported (line 23) but NEVER CALLED.** The masks go directly from Cellpose output to transcript assignment at line 136. There is no expanded mask. The `distance=400` expansion that captures cytoplasmic transcripts is completely absent.

2. **CRITICAL: No dual assignment.** The pipeline produces only `cell_id_reseg` (equivalent to `in_cell` in the notebook). There is no `closest_cell` column. Transcripts outside nuclei but within the cell cytoplasm are assigned `cell_id_reseg = 0` and then discarded at line 138.

3. **Impact on transcript counts:** Without expansion, only transcripts that fall directly on nuclear pixels are assigned. For a typical spatial transcriptomics dataset, nuclear-only assignment captures roughly 30-50% of transcripts. With expansion (distance=400), assignment can reach 70-90%. This means the pipeline may be losing 40-60% of assignable transcripts.

4. **Impact on downstream steps:**
   - Step 4 (Techniques Comparison): Requires `in_cell` vs `closest_cell` distinction for diffusion analysis. The pipeline's single `cell_id_reseg` column cannot support this.
   - Step 5 (Optimal Expansion): Computes the optimal expansion radius, but Step 3 never actually applies any expansion, creating a logical disconnect.
   - Step 6 (Benchmark): The Cellpose resegmentation input will have artificially low transcript counts, skewing benchmark comparisons.

5. **Vectorized vs loop-based assignment:** The pipeline uses vectorized NumPy indexing (`masks[y_coords, x_coords]`), which is orders of magnitude faster than the notebook's per-transcript Python loop. This is a significant performance improvement, but it does not compensate for the missing expansion.

### Severity Rating
**CRITICAL** -- This is the most impactful difference between the notebook and pipeline. The missing `expand_labels` step fundamentally changes the biological interpretation of the resegmentation output. Transcripts in the cytoplasm are lost, cell transcript counts are artificially low, and the dual `in_cell`/`closest_cell` tracking needed for downstream analysis is absent.

### Recommended Fix
```python
# After Cellpose segmentation (line 89), add:
labeled_nuclei = label(masks)
expansion_distance = reseg_config.get('expansion_distance', 400)
expanded_masks = expand_labels(labeled_nuclei, distance=expansion_distance)

# At transcript assignment (replace lines 136-138):
cell_labels_nuclei = labeled_nuclei[y_coords, x_coords]
cell_labels_expanded = expanded_masks[y_coords, x_coords]

df_valid['in_cell'] = cell_labels_nuclei
df_valid['closest_cell'] = cell_labels_expanded
df_valid['cell_id_reseg'] = cell_labels_expanded  # Primary assignment uses expanded

df_assigned = df_valid[df_valid['cell_id_reseg'] > 0].copy()
```

---

## Section 3-6: Build Cell x Gene Matrix (AnnData)

### Notebook Intent
Builds a cell-by-gene count matrix using `pd.crosstab` on the `in_cell` column (nuclear assignment), excluding cell_id=0 (unassigned). Creates an AnnData object:
```python
cellxgene = pd.crosstab(read_positions['in_cell'], read_positions['gene'])
cellxgene = cellxgene.loc[~cellxgene.index.isin([0]), :]
adata = sc.AnnData(cellxgene)
```

Note: The notebook uses `in_cell` (nuclear mask) for the count matrix, not `closest_cell` (expanded mask). This is intentional -- the AnnData represents nuclear expression, while the transcript-level `closest_cell` is used separately for diffusion analysis.

### Paper Methodology
Cell-by-gene matrices from nuclei segmentation serve as the baseline expression profiles. The expanded assignment is used for separate analyses (diffusion, efficiency).

### Pipeline Implementation
Lines 142-148:
```python
gene_col = 'feature_name' if 'feature_name' in df_assigned.columns else 'gene'
cell_gene_matrix = pd.crosstab(df_assigned['cell_id_reseg'], df_assigned[gene_col])
adata = sc.AnnData(cell_gene_matrix)
```

### Specific Differences
1. The pipeline uses `cell_id_reseg` (which equals `in_cell` since no expansion exists). In isolation this matches the notebook's AnnData construction. However, once the expand_labels fix is applied, the pipeline should explicitly use `in_cell` for the AnnData and store `closest_cell` in the transcript CSV.
2. The pipeline dynamically detects the gene column name (`feature_name` vs `gene`), which is a useful generalization.
3. Cell ID 0 is excluded implicitly because `df_assigned` already filters out `cell_id_reseg == 0` at line 138.

### Severity Rating
**LOW** -- The AnnData construction logic is correct given the current (nuclei-only) pipeline. After the expand_labels fix, the column used for crosstab should be explicitly documented.

---

## Section 3-7: Extract Cell Centroids (regionprops)

### Notebook Intent
Computes centroids from `labeled_nuclei` (not expanded masks) using `regionprops`:
```python
centroid_dictx = {}
centroid_dicty = {}
for nucleus_props in regionprops(labeled_nuclei):
    centroid_dictx[nucleus_props.label] = nucleus_props.centroid[0]
    centroid_dicty[nucleus_props.label] = nucleus_props.centroid[1]
```
These centroids are stored both in AnnData `.obs` and used for distance-to-centroid calculations.

### Paper Methodology
Cell centroids are fundamental spatial coordinates used throughout downstream analysis (domain assignment, spatial autocorrelation, visualization).

### Pipeline Implementation
Lines 150-158:
```python
props = regionprops(masks)
centroid_dict = {p.label: p.centroid for p in props}

centroids = [centroid_dict.get(idx, (np.nan, np.nan)) for idx in adata.obs.index.astype(int)]
centroids = np.array(centroids)
adata.obs['y_centroid'] = centroids[:, 0]
adata.obs['x_centroid'] = centroids[:, 1]
```

### Specific Differences
1. **MEDIUM: Missing cell area and perimeter.** The notebook's `regionprops` call provides access to `.area` and `.perimeter` for each cell, but neither the notebook nor the pipeline explicitly extracts these. However, the notebook framework makes it easy to add (since `nucleus_props` is iterated), while the pipeline's tuple-only extraction (`p.centroid`) would need modification.

2. The pipeline uses raw `masks` for regionprops instead of `labeled_nuclei`. As noted in Section 3-3, this is functionally equivalent for Cellpose output.

3. Both store centroids as `x_centroid` and `y_centroid` in `adata.obs`, maintaining column name compatibility.

### Severity Rating
**MEDIUM** -- Centroids are correctly computed. The missing area/perimeter data is a gap if morphological analysis is needed downstream, but it is not used in the current notebook either.

---

## Section 3-8: Distance to Centroid & Transcript Output

### Notebook Intent
After dual assignment, the notebook computes per-transcript distance to the centroid of the assigned (expanded) cell:
```python
read_positions['closest_cell_x'] = read_positions['closest_cell'].map(centroid_dictx)
read_positions['closest_cell_y'] = read_positions['closest_cell'].map(centroid_dicty)
read_positions['distance_to_centroid'] = np.sqrt(
    (read_positions['closest_cell_y'] - read_positions['y'])**2 +
    (read_positions['closest_cell_x'] - read_positions['x'])**2
)
```
This `distance_to_centroid` column is essential for:
- Diffusion analysis in Step 4 (how transcript density varies with distance from cell center)
- Identifying read dispersion patterns (nuclear vs cytoplasmic enrichment)
- Quality metrics comparing segmentation methods

The notebook saves both the AnnData and the full transcript table with all assignment columns:
```python
adata.write(path + '/adata.h5ad')
read_positions.to_csv(path + '/transcripts_with_cell_assignment.csv')
```

### Paper Methodology
Distance-to-centroid distributions are a key analytical tool in the paper. They reveal how transcript assignment behaves as a function of distance from the nuclear center, which is the basis for the optimal expansion analysis in Step 5.

### Pipeline Implementation
Lines 160-167:
```python
adata_out_path = os.path.join(output_dir, f"{config['sample_tag']}_step3_resegmented.h5ad")
transcripts_out_path = os.path.join(output_dir, f"{config['sample_tag']}_step3_transcripts_resegmented.csv")
adata.write(adata_out_path)
df_assigned.to_csv(transcripts_out_path)
```

### Specific Differences
1. **MEDIUM: No `distance_to_centroid` computation.** The pipeline does not compute transcript-to-centroid distances in Step 3. This is partially mitigated by `xb.calculating.dispersion()` being available for use in Step 4, but the notebook computes it here as part of the resegmentation output.

2. **CRITICAL (related to 3-5): No `in_cell`, `closest_cell`, `closest_cell_x`, `closest_cell_y` columns.** The saved transcript CSV only contains `cell_id_reseg`. Downstream steps expecting dual assignment columns will fail or produce incorrect results.

3. The pipeline uses descriptive filenames with `sample_tag` prefix, which is an improvement over the notebook's generic naming.

### Severity Rating
**MEDIUM** for distance_to_centroid alone (since Step 4 can handle it), **CRITICAL** when combined with the missing dual assignment from Section 3-5.

---

## Section 3-9: Domain Assignment (Optional, JSON Polygons)

### Notebook Intent
The notebook does not include domain assignment in the `batch_segmentation_cellpose-Xenium.ipynb` file. However, the `xb/calculating.py` library provides `domainassign()` (line 165), which assigns cells to predefined polygon regions using Shapely point-in-polygon tests:
```python
def domainassign(plsin, adatadom):
    adatadom.obs['region_annotation'] = 'None'
    for sel in plsin['region_annotation'].unique():
        plsub = plsin[plsin['region_annotation'] == sel]
        if plsub.shape[0] > 2:
            coord = np.array(plsub[['y','x']]).tolist()
            coord.append(coord[0])
            poli = Polygon(coord)
            for n in adatadom.obs.index:
                pnt = Point(adatadom.obs.loc[n,'y_centroid'], adatadom.obs.loc[n,'x_centroid'])
                if pnt.within(poli):
                    adatadom.obs.loc[n,'region_annotation'] = sel
```

### Paper Methodology
Domain assignment maps cells to tissue regions (e.g., lobular zones in liver, cortical layers in brain). This spatial annotation is used in Step 5 for computing domain-specific turnover and optimal expansion.

### Pipeline Implementation
Lines 173-238: A complete domain assignment implementation reading from a JSON file with polygon coordinates. Includes:
- Bounding-box pre-filter for performance (lines 210-217)
- Polygon closure check (lines 204-205)
- Point-in-polygon test using Shapely (lines 219-226)
- Re-saving AnnData with `region_annotation` column (line 231)

### Specific Differences
1. **Input format divergence:** The notebook's `domainassign()` expects a DataFrame with `region_annotation`, `y`, `x` columns. The pipeline expects a JSON file with `name` and `coordinates` fields. This is not necessarily wrong -- it may reflect the actual data format used in the project -- but it means the pipeline cannot use the `xb.domainassign()` function directly.

2. **Coordinate convention:** The notebook `domainassign()` uses `(y, x)` ordering for Point construction: `Point(adatadom.obs.loc[n,'y_centroid'], adatadom.obs.loc[n,'x_centroid'])`. The pipeline also uses `(y, x)` ordering (line 222): `pnt = Point(y, x)`. This is consistent.

3. **Bounding-box optimization:** The pipeline adds a bounding-box pre-filter (lines 210-217) before expensive point-in-polygon checks. The notebook iterates over all cells for every polygon. This is a performance improvement.

4. **Visualization:** The notebook's `domainassign()` plots the polygons and cell scatter. The pipeline does not generate any visualization.

### Severity Rating
**LOW** -- Domain assignment is functionally correct. The JSON vs DataFrame input format difference is a design choice. The bounding-box optimization is a valid improvement.

---

## CosMx x/y Coordinate Swap (Cross-cutting)

### Notebook Intent
The notebook handles platform-specific coordinate conventions explicitly:
```python
if ds == 'CosMx':
    read_positions['x'] = read_positions['y_global_px']
    read_positions['y'] = read_positions['x_global_px']
    # Note: x and y are SWAPPED for CosMx
```
For Xenium, coordinates are scaled by the pixel-to-um factor:
```python
if ds in ['Xenium']:
    read_positions['x'] = read_positions['x'] * 4.70588
    read_positions['y'] = read_positions['y'] * 4.70588
```

### Pipeline Implementation
Lines 109-114: Generic coordinate column detection:
```python
if 'x_global_px' in df_transcripts.columns:
    x_col, y_col = 'x_global_px', 'y_global_px'
elif 'global_x' in df_transcripts.columns:
    x_col, y_col = 'global_x', 'global_y'
else:
    x_col, y_col = 'x_location', 'y_location'
```

### Specific Differences
1. **LOW: No CosMx x/y swap.** The pipeline maps `x_global_px` to `x_col` and `y_global_px` to `y_col`, but the notebook swaps them (`x = y_global_px`, `y = x_global_px`). If the pipeline is ever used with CosMx data, transcript-to-mask assignment will be transposed.
2. **No pixel scaling.** The pipeline does not apply the Xenium 4.70588 scaling factor or any platform-specific coordinate transformation. The `um_per_pixel` config value (line 107) is loaded but never used. If transcripts are in micrometer coordinates and the DAPI image is in pixels, the coordinates will not align.
3. For Xenium data specifically, coordinates may already be in pixel space (from Step 0 formatting), so this may not be an issue in practice. But the lack of explicit coordinate transformation is fragile.

### Severity Rating
**LOW** (for Xenium-only use) to **MEDIUM** (if multi-platform support is intended). The unused `um_per_pixel` variable suggests this was planned but not implemented.

---

## Final Judgment

### Summary of Findings

| Section | Issue | Severity |
|---------|-------|----------|
| 3-1 | Device detection | NONE (improvement) |
| 3-2 | DAPI loading | LOW |
| 3-3 | Missing `label()` call | LOW |
| 3-4 | No expanded mask saved | LOW |
| 3-5 | **`expand_labels` imported but never called** | **CRITICAL** |
| 3-5 | **No dual `in_cell`/`closest_cell` assignment** | **CRITICAL** |
| 3-6 | AnnData construction | LOW |
| 3-7 | Missing area/perimeter from regionprops | MEDIUM |
| 3-8 | No `distance_to_centroid` in transcript output | MEDIUM |
| 3-8 | Missing assignment columns in CSV | CRITICAL (linked to 3-5) |
| 3-9 | Domain assignment JSON vs DataFrame format | LOW |
| Cross | CosMx x/y swap not handled | LOW-MEDIUM |
| Cross | `um_per_pixel` loaded but unused | LOW |

### Overall Assessment

The pipeline captures the high-level flow of the notebook (Cellpose segmentation, transcript assignment, AnnData creation, domain assignment) but **misses the most scientifically important step: label expansion**. The `expand_labels(labeled_nuclei, distance=400)` call is the mechanism by which cytoplasmic transcripts are captured, and without it, the resegmentation output is equivalent to nuclei-only segmentation. This defeats the purpose of the resegmentation step, which is to provide an improved cell-level assignment that captures both nuclear and cytoplasmic transcripts.

The dual `in_cell`/`closest_cell` tracking is also absent, which means downstream diffusion analysis (Step 4) cannot distinguish nuclear from cytoplasmic transcripts -- a distinction central to the paper's methodology.

### Prioritized Fix List

**Priority 1 (CRITICAL) -- Add label expansion and dual assignment:**
1. After Cellpose segmentation, apply `label()` and then `expand_labels(labeled_nuclei, distance=400)`.
2. Make `expansion_distance` configurable in `config.yaml` under `resegmentation`.
3. Assign each transcript to both `in_cell` (nuclei mask) and `closest_cell` (expanded mask).
4. Set `cell_id_reseg = closest_cell` as the primary assignment for downstream use.
5. Save the expanded mask alongside the nuclei mask.

**Priority 2 (MEDIUM) -- Add distance-to-centroid computation:**
1. After dual assignment, compute `distance_to_centroid` using expanded cell centroids (matching the notebook formula).
2. Store `closest_cell_x`, `closest_cell_y`, and `distance_to_centroid` in the transcript DataFrame.
3. Save all columns in the output CSV.

**Priority 3 (MEDIUM) -- Extract morphological metrics from regionprops:**
1. Extract `area` and `perimeter` from regionprops and store in `adata.obs`.
2. This supports downstream morphological comparisons across segmentation methods in Step 6.

**Priority 4 (LOW) -- Coordinate handling robustness:**
1. Add CosMx x/y swap logic (or a configurable coordinate mapping).
2. Apply `um_per_pixel` scaling when coordinates are not already in pixel space.
3. Add a config option for `technology` to drive platform-specific coordinate handling (mirroring the `comparison.technology` field already in config).

**Priority 5 (LOW) -- QC visualization:**
1. Add an optional diagnostic plot showing nuclei masks overlaid with expanded masks and transcript positions.
2. Print summary statistics: number of transcripts assigned to nuclei vs expanded cells vs unassigned.

### Estimated Lines of Code for Fixes
- Priority 1: ~25 lines (expansion + dual assignment logic)
- Priority 2: ~10 lines (distance computation + column storage)
- Priority 3: ~5 lines (regionprops extraction)
- Priority 4: ~15 lines (coordinate handling)
- Priority 5: ~20 lines (visualization)
- **Total: ~75 lines of new/modified code**
