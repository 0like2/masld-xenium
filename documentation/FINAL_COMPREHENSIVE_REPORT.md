# Xenium Pipeline Master Report: Comprehensive Analysis & Status
**Date:** 2026-01-01
**Version:** 2.0 (Master Compilation)

## 1. Executive Summary
This document serves as the single source of truth for the Xenium Benchmarking Pipeline. It consolidates the technical architecture, visualization capabilities, gap analysis, and future roadmap.
We have successfully refactored a collection of 20+ disparate Jupyter Notebooks into a **robust, automated Python pipeline** that matches publication-grade standards.

### Key Achievements
*   **Automation:** Steps 0-6 (Preprocessing to Optimization) are fully automated via `run_all.sh`.
*   **Visualization:** Implemented critical publication figures (Images 3e, 3f, 5c).
*   **Robustness:** Added "Ground Truth" logic from external scRNA-seq references and ROI optimization simulations.

---

## 2. Technical Deep Dive (Pipeline Logic)

### 2.1 Ground Truth Generation (Step 0)
*   **Objective**: Create a "Gold Standard" cell type label for every cell to benchmark accuracy.
*   **Method**: `step7_ground_truth.py` (integrated into Step 0).
    1.  **Reference**: Download curated scRNA-seq data from Cellxgene Census.
    2.  **Mapping**: Resolve gene symbol discrepancies (Ensembl vs Symbol) using intersection logic.
    3.  **Transfer**: Use `scanpy.tl.ingest` to transfer labels from Reference (scRNA-seq) to Target (Xenium) based on shared expression profiles.

### 2.2 Segmentation-Free Analysis (Step 2)
*   **Objective**: Analyze transcript distribution without relying on imperfect cell boundaries.
*   **Method**:
    *   **Distance Analysis**: Computes Euclidean distance from every transcript to the nearest nucleus centroid.
    *   **Overlap/SSAM**: (Partially implemented) Uses Z-axis signal coherence and Kernel Density Estimation to assess cell packing.

### 2.3 ROI Optimization (Step 4)
*   **Objective**: Find the optimal "expansion distance" to define a cell's boundary.
*   **Logic**:
    *   **Simulation**: Iteratively expands the boundary (0µm to 15µm in 0.5µm steps).
    *   **Trade-off**: Measures **Capture Efficiency** (more reads included) vs **Purity** (nuclear fraction).
    *   **Result**: The distance where purity drops significantly (Elbow Point) is selected as optimal (e.g., ~10µm for Brain).

---

## 3. Visualization Capabilities Audit

We have implemented a `visualize_results.py` module and a regeneration script `run_vis_all.py` to produce consistent plots.

### A. Implemented Visualizations (Status: ✅ Ready)

| Category | Plot Type | Description | File Pattern |
| :--- | :--- | :--- | :--- |
| **QC & Preprocessing** | Violin Plots | QC Metrics (Counts, Genes) | `*_step1_qc_violin.png` |
| | PCA Variance | Variance ratio per PC | `*_step1_pca_variance.png` |
| | UMAP | Leiden Clustering Map | `*_step1_umap_plot.png` |
| | Spatial Scatter | Cell Types on XY Coords | `*_spatial_clusters.png` |
| **Segmentation-Free** | Box/Bar Plots | Gene Distance to Nucleus | `*_step2_gene_dist_boxplot.png` |
| **Optimization** | Dual-Axis Line | Purity vs Capture Optimization (Image 3f) | `*_step4_optimization_curve.png` |
| **Advanced** | **Violin Plot** | Marker Gene Expression (Image 3e) | `*_marker_genes_violin.png` |
| | **Seg. Overlay** | Transcripts + Centroids/Boundaries (Image 5c) | `*_segmentation_overlay_roi.png` |

### B. Missing Visualizations (Status: ❌ Todo)

| Category | Plot Type | Reason for Omission | Priority |
| :--- | :--- | :--- | :--- |
| **Spatial Stats** | Neighborhood Heatmap (Image 5e) | Requires `squidpy` module integration. | High |
| | Spatial Co-occurrence | Requires `squidpy`. | High |
| **ARI Heatmap** | Heatmap Image | Data exists in CSV, PNG generation pending. | Medium |
| **Raw Data** | DAPI / Raw Transcripts | Too heavy for standard pipeline; skipped for speed. | Low |

---

## 4. Pipeline Step-by-Step Analysis

| Step | Module | Status | Visualization Support | Gap / Improvement Area |
| :--- | :--- | :--- | :--- | :--- |
| **Step 0** | `step0_formatting` | ✅ | N/A | Format validation check (ensure `x_centroid` exists). |
| **Step 1** | `step1_preprocess` | ✅ | **Excellent** (QC, UMAP, PCA) | Add cluster annotation (Auto-annotate based on markers). |
| **Step 2** | `step2_seg_free` | ✅ | **Good** (Distance Plots) | Implement full SSAM (Density Map) if needed. |
| **Step 4** | `step4_expansion` | ✅ | **Excellent** (Opt. Curve) | Allow user-configurable expansion ranges in `config.yaml`. |
| **Step 6** | `step6_optimization` | ✅ | **Partial** (CSV Only) | **CRITICAL**: Fixed performance bug (sampling), need Heatmap PNG. |
| **Step 7** | `step7_domains` | ⚠️ | None | **Blocked**: `SpaGCN`/`cmake` build failure. |
| **Step 8** | `step8_svf` | ✅ | Spatial Maps | SVG detection implies heavy compute; currently optional. |

---

## 5. Dataset Status Report

**1. Human Brain (`human_brain`)**
*   **Processing**: Complete.
*   **Issues Resolved**: Fixed missing `leiden` key for Violin Plot (auto-detects `leiden_1_4`).
*   **Visualizations**: All standard plots + Marker Violin + Optimization Curve + Segmentation Overlay (Fallback Mode).

**2. Human Breast (`h_breast_1`)**
*   **Processing**: Complete.
*   **Performance**: Large dataset (~7GB) processed successfully. UMAP visualization takes time but completes.

---

## 6. How to Use

### Run Full Pipeline
```bash
# Process all configured datasets from Step 0 to 6
bash run_all.sh
```

### Regenerate Visualizations Only
```bash
# Detects output files and redraws all plots (including new features)
python run_vis_all.py
```

### Generate Specific Plot (Violin)
```bash
# Targeted script for marker genes
python run_violin_only.py
```

---

## 7. Future Roadmap (Recommended)

1.  **Implement Step 6 Visualization**: Create a heatmap PNG generator for the ARI/Silhouette scores CSV.
2.  **Integrate Squidpy**: Add `step9_spatial_stats` to compute and plot Neighborhood Enrichment and Co-occurrence (High Scientific Value).
3.  **Refine Segmentation Overlay**: Add logic to find and use raw polygon files (`.parquet`) if available, for "True" segmentation overlay (currently using Centroids).
