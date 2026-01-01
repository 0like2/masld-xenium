# Xenium Pipeline Detailed Analysis

## Step 0: Formatting & Ground Truth Integration

### 1. Visualization
- **Current Status**: `visualize_results.py` focuses on Step 1 (PCA, Clustering), Step 2 (Gene Distances), Step 4 (Expansion), and Step 6 (Simulation).
- **Step 0 Specific**: There is no dedicated visualization for Step 0 (raw formatted data).
- **Recommendation**: Add a basic spatial plot for Step 0 to verify data loading and spatial coordinates before preprocessing. Visualizing 'total_counts' spatially could be useful quality check.

### 2. scRNA-seq Integration
- **Implemented**: Yes, via `step7_ground_truth` module integrated into `step0_format_xenium`.
- **Functionality**: Downloads reference from Cellxgene Census, performs `scanpy.tl.ingest` for label transfer.
- **Data Flow**: `format_xenium` -> `download_reference` -> `transfer_labels` -> `adata.obs['ground_truth_celltype']`.

### 3. Code Comparison (Notebook vs Pipeline)
- **Notebook**: `0_0_Formatting xenium to anndata .ipynb`
    - Uses `xb.formatting.format_xenium_adata` (2022), `format_xenium_adata_2023`, `format_xenium_adata_mid_2023`.
    - Also calls `xf.format_background`.
    - Handles multiple datasets manually.
- **Pipeline**: `xenium_pipeline_main.py`
    - Uses `format_xenium_adata_mid_2023` by default (hardcoded for newer format).
    - Includes `xf.format_background`.
    - Adds Ground Truth generation (new feature).
- **Discrepancy**: Pipeline assumes "mid_2023" format (2024 standards). If user has older data, it might fail. The config has a `xenium_format` section but the code in `step0` currently calls `format_xenium_adata_mid_2023` directly.

### 4. Missing Items / Improvements
- **Format Selection**: The pipeline should respect the `xenium_format` config to choose the correct formatting function (`format_xenium_adata` vs `format_xenium_adata_2023` vs `format_xenium_adata_mid_2023`).
- **Validation**: No explicit check if the output `h5ad` matches expected structure (e.g., check for `x_centroid` columns).
- **Ground Truth**: The reference download might be slow; we should ensure it's cached (implemented in `step7` module, verified).

---



## Step 1: Preprocessing

### 1. Visualization
- **Current Status**: `visualize_results.py` provides:
    - PCA Variance Ratio (`sc.pl.pca_variance_ratio`).
    - PCA 2D Projection (`sc.pl.pca`).
    - Spatial Cluster Map (`sns.scatterplot` on x/y centroids).
- **Original Code**: `xb.preprocessing.preprocess_adata` (not used in pipeline) includes:
    - `hvg.pdf`: Highly Variable Genes plot.
    - `umap_*.pdf`: UMAP plots (standard for scRNA-seq).
    - `deg_*.pdf` & `deg_dotplot_*.pdf`: Differential Expression analysis (identifies cluster markers).
    - `plot_cell_counts`: Histograms of genes/counts per cell (Vital QC).
- **Gap Analysis**: The pipeline uses `main_preprocessing` (simulation-focused, light) instead of `preprocess_adata` (exploration-focused, rich). We are missing **Critical QC Plots** (filtering effects), **UMAP** (better cluster visualization than PCA), and **Biology Checks** (DEG/Markers).

### 2. Logic Comparison
- **Notebook**: `1_6_Batch_preprocessing_real_Xenium_datasets.ipynb`.
    - Primarily runs `allcombs` or `main_preprocessing`.
    - Parameters: `mincounts=10`, `mingenes=3`, `neigh=15`, `target_sum=100`, `log=True`.
- **Pipeline Implementation**: `step1_preprocess` calls `xp.main_preprocessing`.
    - **Enforced Clustering**: The most significant logic is the `while` loop in `main_preprocessing` that iteratively adjusts Leiden resolution (±0.05/0.1) to force the number of clusters to match `total_clusters` (default 30). This is somewhat arbitrary and might over/under-cluster heterogeneous/homogeneous tissues.
    - **Parameters**: Defaults align with notebook (`target_sum=100`, `scale=False`).

### 3. Missing Items / Improvements
- **QC Dashboard**: Add a pre/post filtering QC plot (Violin plots of `n_counts`, `n_genes`) to `visualize_results.py`.
- **UMAP**: Replace or augment PCA visualization with UMAP calculation and plotting (`sc.pp.neighbors` + `sc.tl.umap` + `sc.pl.umap`).
- **Cluster Annotation**: Computing DEGs (`sc.tl.rank_genes_groups`) and saving a top-marker CSV would allow users to verify biological identity of clusters.
- **Configurable Clusters**: The `total_clusters=30` target should be made configurable in `config.yaml` to allow flexibility for different tissues.

---


## Step 2: Segmentation-Free Analysis

### 1. Visualization
- **Current Status**: `visualize_results.py` includes `visualize_step2_gene_distances`.
    - It plots "Average Distance to Nucleus by Gene" (Barplot).
    - **Original Code**: `2_1_batch_processing_distance_to_nuclei_across_samples.ipynb` generates:
        - Boxplots of distance per high-dispersion gene.
    - **Note**: The current pipeline visualization uses a Barplot of the *mean* distance, whereas the notebook uses boxplots which show distribution. Boxplots would be richer.

### 2. Logic Comparison
- **Notebook**: Uses `xb.calculating.dispersion` to calculate Euclidean distance between each read and its assigned cell centroid.
- **Pipeline Implementation**: `step2_segmentation_free_analysis` directly calls `xb.calculating.dispersion`.
    - Logic is identical: computes distance and overlaps with nucleus.
    - Saves results to a CSV file or returns DataFrame.


### 3. Missing Items / Improvements
- **Distribution Plots (Critical)**: The notebook generates a boxplot for *each of the top variable genes*, showing the distribution of distances for all reads of that gene.
    - **Why it matters**: A simple mean (current pipeline) hides the spread. A gene might have a high mean because of outliers, or bimodal distribution. Boxplots clearly show if a gene is tightly packed in the nucleus vs. diffused in the cytoplasm.
    - **Action**: Implement `sns.boxplot` or `sns.violinplot` in `visualize_results.py` for the top 20 most spatially variable genes.
- **QC Metrics**: The notebook outputs "Read specific dispersion metrics". We should ensure our CSV output includes not just mean, but also variance/dispersion indices to help users select "Nuclear Markers" systematically.

### 4. Unimplemented Modules (Scope Check)
The `notebooks/2_segmentation_free_analysis/` directory contains other significant analyses NOT present in the current pipeline:
- **`2_3_brain_cell_overlaps.ipynb`**: Analyzes signal coherence and overlapping cells using `ovrlpy`. This is advanced QC/signal analysis.
- **`2_4_brain_ssam.ipynb`**: Implements **SSAM** (Segmentation-Free Sparse Patch-Based Analysis) using Kernel Density Estimation (KDE) and integrates with scRNA-seq data.
    - **Current Status**: Completely missing from `xenium_pipeline_main.py`.
    - **Recommendation**: Decide if SSAM is required for this phase. If so, it requires a new step (e.g., Step 2B) and external dependencies (`ssam`, scRNA-seq reference).

---


## Step 4: Optimal Expansion

### 1. Visualization
- **Current Status**: `visualize_results.py` includes `visualize_step4_optimal_expansion`.
    - It plots `Optimal Cell Expansion Parameters` as a bar plot.
- **Original Code**: `4_1_Optimal_expansion_multisection.ipynb` generates:
    - Dot plots (Scatter) of "Cell Type Purity" vs "Signal Capture" for different expansion distances.
    - This allows choosing the trade-off point (maximum purity and capture).
- **Comparison**: The pipeline visualization shows the *result* (barplot of values), but not the *trade-off optimization curve* which is critical for understanding *why* a specific distance was chosen.

### 2. Logic Comparison
- **Notebook**:
    - Calculates "Signal Capture" (reads assigned to correct cell type / total reads).
    - Calculates "Purity" (reads assigned to correct cell type / total reads assigned to that cell).
    - Iterates varying expansion distances (e.g., 5, 10, 15, ... 100 µm/pixels?).
- **Pipeline Implementation**: `step4_optimal_expansion` **hardcodes values** based on literature (`10.71` and `5.65` µm).
    - **CRITICAL GAP**: It does **NOT** iterate distances or perform optimization on the current dataset.
    - It merely calculates `mean_nuclear_distance` for validation but does not use it to select parameters.

### 3. Missing Items / Improvements
- **Optimization Plot**: We should implement the Purity vs Capture trade-off plot. It's much more informative than a simple bar chart of final values.
- **Configurability**: Allow users to set the range of expansion distances to test in `config.yaml`.

---

## Step 6: Preprocessing Simulation & Optimization

### 1. Visualization
- **Current Status**: `visualize_results.py` includes `visualize_step6_simulation`.
    - Plots "Stability (ARI)" vs "Quality (Silhouette)".
    - Annotated "Best ARI" and "Best Silhouette" points.
    - **Status**: **Excellent**. Matches the notebook's intent well. No changes needed for visualization.

### 2. Logic Comparison & Critical Bug
- **Notebook**:
    - Generates ~1000 simulated datasets (subsampling + noise).
    - Checks 618 preprocessing combinations.
- **Pipeline Implementation**: `step6_preprocessing_simulation` handles the logic but has a **Critical Performance Bug**.
    - **Code Issue**: It correctly creates a subsample `adata_sample = self.adata[:maxcell, :].copy()`.
    - **BUT**: It calls `allcombs(self.adata)` (passing the *full* 100% dataset) instead of `adata_sample`.
    - **Impact**: On large datasets (e.g., millions of cells), running 618 full clustering jobs will cause the pipeline to hang or crash indefinitely.

### 3. Missing Items / Improvements
- **Fix Sampling**: Change the function call to use `adata_sample`.
- **Report Generation**: Currently logs to stdout. It should output a `step6_recommendation.json` or text file stating: "Recommended Params: Norm=True, Log=True..." so downstream steps can read it programmatically.

---

## Step 7: Spatial Domains (New Finding)
**Notebooks**: `7_domain_exploration/7_1_SpaGCN_domains.ipynb` (and others like Banksy, DeepST).
- **Status**: **Completely Missing**. The current pipeline stops after Step 6.
- **Why it matters**: Identifying spatial domains (tissue architecture) is a key output of spatial transcriptomics. The original study extensively compares methods here.
- **Action**: Implement a new `step7_spatial_domains` module using **SpaGCN** (Python-based, used in `7_1`).
    - **Logic**: It requires an iterative search for hyper-parameters (`l` and `res`). We should implement this automated search loop as seen in the notebook logs.

## Step 8: Spatially Variable Features (New Finding)
**Notebooks**: `8_SVF_identification/8_1_batch_processing_SpatialDE_SVF.ipynb`.
- **Status**: **Completely Missing**.
- **Action**: Implement a new `step8_svf` module using **SpatialDE** (Python).
    - **Note**: SVF detection is computationally expensive. Should be an optional step.



