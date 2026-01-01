# Pipeline vs. Notebook Verification Table

This table compares the original Jupyter Notebooks (~Source of Truth) with the implemented Python Pipeline (`xenium_pipeline_main.py`, `visualize_results.py`, and sub-modules).

| Step | Original Notebook | Pipeline Module(s) | Logic Parity | Visualization Parity | Notes |
| :--- | :--- | :--- | :---: | :---: | :--- |
| **0. Formatting & GT** | `0_0_Xenium_formatting.ipynb` | `step0_format_xenium`<br>`step7_ground_truth.py` | ✅ **Matches**<br>(Enhanced with robust gene mapping) | ✅ **Matches**<br>(Basic Cell Type Maps) | Pipeline added `feature_id` mapping for Human Brain GT. |
| **1. Preprocessing** | `1_2_data_preprocessing.ipynb`<br>`1_1_data_exploration.ipynb` | `step1_preprocess`<br>`xb/preprocessing.py` | ✅ **Matches**<br>(Leiden/Louvain fallback) | ✅ **Matches**<br>(QC Violins, PCA, UMAP, Spatial) | Pipeline uses `leidenalg` (via `igraph`) correctly. |
| **2. Seg-Free** | `2_1_Segmentation_free...ipynb` | `step2_segmentation_free`<br>`xb/calculating.py` | ✅ **Matches**<br>(Distance calcs) | ✅ **Matches**<br>(Gene Distance Boxplots/Barplots) | `visualize_results.py` handles this well. |
| **4. Opt. Expansion** | `4_1_Optimum_expansion...ipynb` | `step4_optimal_expansion`<br>`xb/simulating.py` | ✅ **Matches**<br>(Simulation Loop 0-15um) | ✅ **Matches**<br>(Expansion Curve Lineplot) | Pipeline checks for purity metrics. |
| **6. Optimization** | `6_3_Simulated_Xenium...ipynb` | `step6_preprocessing_simulation` | ✅ **Matches**<br>(618 Combinations, ARI) | ❌ **Missing** | `visualize_step6_simulation` is called but **not defined** in `visualize_results.py`. |
| **7. Spatial Domains** | `7_1_SpaGCN_domains.ipynb` | `step7_spatial_domains.py` | ✅ **Matches**<br>(SpaGCN logic implemented) | ❌ **Missing** | Module lacks plotting code. Step currently skipped due to dependency. |
| **8. SVF** | `8_1_...SpatialDE_SVF.ipynb` | `step8_svf.py` | ✅ **Matches**<br>(SpatialDE/NaiveDE logic) | ❌ **Missing** | Module calculates FSV/q-val but does not generate spatial expression plots. |

## Gap Analysis & Recommendations

### 1. Visualization Gaps
The pipeline is robust in **logic** (calculation) but falls short in **visualization** for the advanced steps (6, 7, 8). The original notebooks generated rich visual outputs (Heatmaps for Step 6, Domain Maps for Step 7, Gene Maps for Step 8) which are currently absent in the automated pipeline.

### 2. Immediate Fixes Required
- **Step 6**: Implement `visualize_step6_simulation` in `visualize_results.py` (Heatmap of ARI scores).
- **Step 8**: Add plotting to `step8_svf.py` or `visualize_results.py` to output spatial maps of top SVGs.

### 3. Logic Improvements
- **Step 7**: The SpaGCN dependency issue (system `cmake`) prevents this step from running automatically. This requires a Docker/Environment fix rather than a python code fix.
