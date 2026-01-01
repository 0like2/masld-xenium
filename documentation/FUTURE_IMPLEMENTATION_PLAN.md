# Future Implementation Plan

`notebook_visualization.md`의 분석 결과를 바탕으로, 향후 파이프라인에 추가 구현할 수 있는 기능들을 중요도 순으로 정리한 문서입니다.

## 1. High Priority (즉시 도입 권장)
분석 결과의 **생물학적 해석**을 돕기 위해 시급한 시각화 기능입니다.

### [1-1] Marker Gene Heatmap / Dotplot
*   **Original Code**: `sns.heatmap(gene_sorted)` (Notebook 1-5)
*   **Need**: 현재 UMAP 클러스터가 숫자로만 표시되어, 실제 해당 클러스터가 어떤 세포 타입(e.g., Neuron, Astrocyte)인지 알기 어렵습니다.
*   **Implementation Plan**:
    1.  `visualize_results.py`에 `visualize_marker_genes(adata)` 함수 추가.
    2.  `scanpy.tl.rank_genes_groups`를 실행하여 클러스터별 상위 유전자(Top 5) 추출.
    3.  `scanpy.pl.dotplot` 또는 `scanpy.pl.heatmap`을 사용하여 시각화 및 저장.
*   **Value**: 연구자가 클러스터링 결과를 한눈에 검증할 수 있게 됨.

---

## 2. Medium Priority (기능 확장)
파이프라인의 **분석 깊이(Depth)**를 더하기 위한 기능들입니다.

### [2-1] User-defined Gene UMAP
*   **Original Code**: `sc.pl.umap(color='Gfap')`
*   **Need**: 현재 Step 8(SVF)에서 통계적으로 유의한 유전자만 자동으로 보여주지만, 연구자는 'Gfap'이나 'Olig2' 같은 **"관심 유전자(Known Marker)"**의 분포를 보고 싶어 합니다.
*   **Implementation Plan**:
    1.  `config.yaml`에 `genes_to_plot: ['Gfap', 'Olig2']` 리스트 추가.
    2.  `visualize_results.py`에서 해당 리스트를 읽어 `sc.pl.umap` 및 `sc.pl.spatial`을 그리는 로직 추가.

### [2-2] Spatial Regions (Step 7 재가동)
*   **Original Code**: `sc.pl.spatial(color='region_annotation')`
*   **Need**: 세포를 단순 클러스터링하는 것을 넘어, "Cortex Layer 1", "Striatum" 같은 **공간적 영역(Spatial Domain)**을 자동으로 나누는 기능입니다.
*   **Implementation Plan**:
    1.  현재 빌드 문제(`cmake` dependency)로 막혀 있는 Step 7 (`step7_spatial_domains.py`)의 `SpaGCN` 의존성을 해결하거나 `Squidpy`로 대체.
    2.  영역 분석이 가능해지면 `region_level` 시각화를 활성화.

---

## 3. Low Priority (탐색적 분석 / 고비용)
자동화 벤치마킹보다는 **개별 연구(Deep Dive)**에 적합한 기능들입니다.

### [3-1] Advanced SSAM (Vector Field)
*   **Original Code**: `ds.plot_celltypes_map(...)`
*   **Status**: 현재 계산 부하 문제로 단순 밀도 맵(Density Map)으로 대체됨.
*   **Plan**: 전체 이미지를 처리하지 않고, 특정 **ROI(관심 영역)**만 잘라서 정밀 분석하는 "Interactive Mode" 모듈을 별도로 개발.

### [3-2] DAPI / Cellpose Re-segmentation
*   **Original Code**: Notebook 3-2 전반
*   **Status**: 벤치마킹 범위를 벗어남 (Raw Image Processing).
*   **Plan**: Xenium 제공 Segmentation 결과가 만족스럽지 않을 경우를 대비해, **Step 0-B (Alternative Segmentation)** 모듈로 분리하여 개발 고려. (매우 큰 작업 소요)

---

## 4. Paper-Grade Visualizations (Derived from Uploaded Images)
These are advanced visualizations identified from the provided paper figures, essential for high-quality benchmarking.

### High Priority: Segmentation Comparison (Images 5c, 5d, 5h)
*   **ARI Heatmap (Image 5d)**:
    *   **Description**: Heatmap showing Adjusted Rand Index (ARI) similarity between different segmentation methods (Xenium vs Cellpose vs Baysor).
    *   **Implementation**: Requires running multiple segmentations (Series 3 & 5) and computing pairwise ARI.
*   **Cell Count Bar Charts (Image 5h)**:
    *   **Description**: Grouped bar chart comparing number of cells detected per cell type across methods.
    *   **Implementation**: Simple Pandas bar plot after merging `adata.obs` from different methods.
*   **Segmentation Mask Overlays (Image 5c)**:
    *   **Description**: Zoomed-in spatial crops showing cell boundaries from different methods side-by-side or overlaid.
    *   **Implementation**: Use `squidpy` or `matplotlib` with polygon patches.

### Medium Priority: Advanced Quality Control (Images 2f, 3e, 3f)
*   **Nuclear vs Cytoplasmic Gene Boxplots (Image 2f)**:
    *   **Description**: Boxplots showing "Distance to Centroid" for top nuclear vs cytoplasmic genes.
    *   **Status**: Step 2 generic analysis exists, but needs this specific "Top 5 vs Bottom 5" boxplot visualization.
*   **Marker Gene Violin Plots (Image 3e)**:
    *   **Description**: Stacked violin plots for key marker genes (Sox10, Pvalb, etc.) to assess cell type specificity.
    *   **Implementation**: `sc.pl.stacked_violin`.
*   **Expansion Optimization Curve (Image 3f)**:
    *   **Description**: "Proportion of reads" vs "Distance" curve.
    *   **Status**: Step 4 computes this data (`step4_optimal_expansion.csv`), but we need a dedicated plotting function for it.

### Low Priority: 3D/Z-Axis Analysis (Image 2e)
*   **Z-Axis Cross Sections (Image 2e)**:
    *   **Description**: "Side view" (X-Z or Y-Z) of the tissue to show cell layers.
    *   **Implementation**: select specific Y slice, plot X vs Z scatter.

---

## 5. Summary

| 우선순위 | 기능명 | 예상 소요 시간 | 가치 |
|---|---|---|---|
| 🔥 **High** | **Marker Gene Dotplot** | 1~2 시간 | 해석력 극대화 |
| ☁️ **Medium** | **Custom Gene Plot** | 1 시간 | 사용자 편의 증대 |
| ☁️ **Medium** | **Fix Step 7 (Regions)** | 4~8 시간 | 공간 미세 구조 분석 가능 |
| 📉 **Low** | **Full SSAM / Raw Image** | 수 일(Days) | 연구 목적 확장 시 고려 |
