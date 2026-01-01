# Ground Truth ARI Implementation Plan

## Goal
Implement a "Ground Truth" Adjusted Rand Index (ARI) calculation module ("Step 7") for the Xenium pipeline. This module will validate the clustering performance of the Xenium analysis by comparing it against detailed cell-type annotations derived from an external single-cell RNA sequencing (scRNA-seq) reference atlas.

## User Review Required
> [!IMPORTANT]
> **Data Download Size**: Downloading scRNA-seq reference data (e.g., from Cellxgene Census) can involve large files. Ensure sufficient disk space and internet bandwidth.
> **Methodology**: We will use `scanpy.tl.ingest` for label transfer as a robust and standard baseline. More advanced methods (like scVI) could be added later but `ingest` is sufficient for a first implementation.

## Proposed Changes

### [NEW] Restoring Missing Visualizations (Steps 1 & 2)
The current pipeline uses a lightweight preprocessing version that omits critical QC and biological validation plots present in the original research code (`xb/preprocessing.py`).

#### [MODIFY] `Xenium_benchmarking/visualize_results.py`
- **Implement QC Plots**:
    - Add `plot_cell_counts` equivalent: Histograms/Violin plots for `n_genes` and `n_counts`.
- **Implement Dimensionality Reduction Plots**:
    - Add UMAP calculation (`sc.tl.umap`) and plotting (`sc.pl.umap`). Currently only PCA is plotted.
- **Implement Biological Validation**:
    - Add Differential Expression Gene (DEG) analysis (`sc.tl.rank_genes_groups`).
    - Add Dotplots (`sc.pl.rank_genes_groups_dotplot`) for cluster markers.
- **Implement Segmentation-Free Distribution Plots (Step 2)**:
    - Instead of just CSV stats, generate **Boxplots/Violin plots** for gene-to-nuclei distances (Top 20 spatially variable genes).

### [NEW] `Xenium_benchmarking/step7_ground_truth.py`
New module to handle the entire Ground Truth workflow.
- **Function `download_reference_data(tissue_type)`**:
    - Uses `cellxgene_census` to query and download Human Brain or Breast scRNA-seq data.
    - Filters for relevant disease states if applicable (or healthy reference).
    - Saves the reference `AnnData` to `data/reference/`.
- **Function `transfer_labels(adata_xenium, adata_reference)`**:
    - Preprocesses reference data (normalize, log1p, find highly variable genes).
    - Aligns Xenium data to Reference features (gene intersection).
    - Trains a classifier (or uses `scanpy.tl.ingest`) to predict labels for Xenium cells based on Reference.
    - Adds `ground_truth_celltype` column to `adata_xenium.obs`.
- **Function `calculate_ground_truth_ari(adata_xenium)`**:
    - Calculates ARI between `adata_xenium.obs['leiden_1_4']` (from Step 1) and `adata_xenium.obs['ground_truth_celltype']`.
    - Returns the ARI score.

### [MODIFY] `Xenium_benchmarking/xenium_pipeline_main.py`
- Import `step7_ground_truth`.
- Add `step7_ground_truth_analysis()` method to `XeniumPipeline` class.
- Call `download_reference_data`, `transfer_labels`, and `calculate_ground_truth_ari` sequentially.
- Log and save the results (ARI score) to `output/`.

### [MODIFY] `Xenium_benchmarking/config.yaml` (if exists, or via kwargs)
- Add parameters:
    - `reference_tissue`: "brain" or "breast" (inferred from sample tag or explicit).
    - `download_dir`: `data/reference/`.

### [NEW] `Xenium_benchmarking/step7_spatial_domains.py` (New Step)
**Goal**: Match `7_1_SpaGCN_domains.ipynb`.
- **Implement**: `run_spagcn(adata)`
    - Auto-search for `l` parameter (contribution of neighbors).
    - Auto-search for `res` parameter to match target cluster number or stability.
    - Save results to `obs['spatial_domain']`.

### [NEW] `Xenium_benchmarking/step8_svf.py` (New Step)
**Goal**: Match `8_1_batch_processing_SpatialDE_SVF.ipynb`.
- **Implement**: `run_spatialde(adata)`
    - Run SpatialDE significance test.
    - Identify spatially variable genes.
    - Save results to CSV/AnnData.

### [MODIFY] `Xenium_benchmarking/xenium_pipeline_main.py`
- **Step 1**: **Logic Verified** (calls library). No changes needed to logic, only visualization (see above).
- **Step 4**: Fix hardcoded expansion (see above).
- **Step 6**: Fix sampling bug (see above).
- **Step 7/8**: Integrate new modules as optional steps.

### Phase 3: Implementing Missing Modules (Current)
Based on the audit, the following modules and visualizations are being added to reach full parity with the notebooks.

#### [NEW] [step2_overlaps.py](file:///data1/project/20rak/masld_xenium/Xenium_benchmarking/step2_overlaps.py)
- **Goal**: Analyze cell overlaps and z-axis signal coherence.
- **Logic**: Histogram of overlaps, Z-distribution stats.
- **Dependencies**: `ovrlpy` (Need to verify availability).

#### [NEW] [step2_ssam.py](file:///data1/project/20rak/masld_xenium/Xenium_benchmarking/step2_ssam.py)
- **Goal**: Segmentation-free analysis using KDE (SSAM equivalent).
- **Logic**: KDE density estimation, Vector field calculation (if feasible).
- **Dependencies**: `ssam` (Need to verify availability).

#### [MODIFY] [visualize_results.py](file:///data1/project/20rak/masld_xenium/Xenium_benchmarking/visualize_results.py)
- **Fix**: Implement missing plotting functions.
- `visualize_step6_simulation`: Simulation ARI Heatmaps.
- `visualize_step8_svf`: Spatial Gene Expression Maps for top SVGs.

#### [MODIFY] [xenium_pipeline_main.py](file:///data1/project/20rak/masld_xenium/Xenium_benchmarking/xenium_pipeline_main.py)
- **Integration**: Add calls to `step2_overlaps` and `step2_ssam` within `run_full_pipeline`.

#### [NEW] [h_breast_2 Dataset]
- **Source**: Zenodo Record 11121221.
- **Status**: Downloading (background process).
- **Goal**: Expand benchmarking to a second breast tissue sample.

- **Structure Scores (Notebook 1_7)**: Identified as missing but separate from core benchmarking scope.
- **Comparison Plots (Series 3/5)**: Out of scope for single-sample pipeline.

## Verification Plan
1. **Dependency Check**: Verify imports for `ovrlpy` and `ssam`.
2. **Visualization Test**: Run `visualize_results.py` manually on existing `human_brain` output to generate new plots immediately.
3. **Module Test**: Run new steps in isolation on `human_brain`.
4. **Data Verification**: Confirm `h_breast_2.h5ad` integrity once download completes.
- **Visual Check**:
    - Check logs for successful download and processing.
    - Verify `output/` contains the ARI score.
    - (Optional) Plot UMAP of Xenium data colored by `ground_truth_celltype` vs `leiden_1_4` to visually compare.
