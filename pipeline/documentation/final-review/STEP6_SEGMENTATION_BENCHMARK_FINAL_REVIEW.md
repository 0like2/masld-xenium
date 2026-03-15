
## 1. Overview

### 1.1 What This Step Does
Step 6은 여러 segmentation 방법을 **하나의 통일된 벤치마크 프레임워크**에서 비교한다:
1. **Nuclei** (Xenium 기본 segmentation)
2. **Cellpose** (Step 3 재분할)
3. **Expansion** (Step 5 최적 확장)
4. **Baysor** (transcript 기반 segmentation, 선택적)

모든 방법의 결과를 concatenate하여 동일 전처리를 적용한 후, UMAP, cell type annotation, 정량적 메트릭으로 비교한다.

### 1.2 Paper Context
논문의 "Baysor and Cellpose outperform standard Xenium segmentation" 섹션 (p.817-818):
- "We benchmarked the segmentation methods Baysor, MESMER, Watershed, Cellpose and Clustermap against the segmentations provided by 10x Genomics"
- "We found that Baysor in combination with Cellpose segmentation outperformed the rest of the strategies"
- NMP (Negative Marker Purity)와 proportion of assigned reads를 주요 비교 메트릭으로 사용
- Fig. 3에서 segmentation 비교의 핵심 결과 시각화

### 1.3 논문 Figure 매칭

> **주의**: 논문 Fig. 3의 segmentation 비교는 **8가지 방법/변형** (Xenium cell, Xenium nuc, MESMER, Clustermap, Cellpose, Watershed, Baysor conf=0, Baysor conf=0.8)을 사용한다. 파이프라인은 이 중 4가지 (Nuclei, Cellpose, Expansion, Baysor)를 구현하며, 추가로 Expansion (논문에서는 별도 Figure 3a-b에서 다룸)을 포함한다.

| 논문 Figure | 논문 원문 캡션 | 논문 내용 | Step 6 구현 | 노트북 출처 |
|------------|-------------|---------|------------|------------|
| **Fig. 3a** | "Mouse brain region with reads overlaid on DAPI staining, colored by distance to the nearest cell centroid" | DAPI 위 reads (거리별 색상) + oligodendrocyte PCC line plot (nuclear vs background signature) | Step 5 결과 활용 (Step 6 아님) | `notebooks/4_optimal_expansion/4_1_Optimal_expansion_multisection.ipynb` |
| **Fig. 3b** | "Predicted optimal expansion by cell type" | Cell type별 optimal expansion 거리 (cell edge/nuclei edge), Default Xenium 15µm 기준선 포함 | Step 5 결과 활용 (Step 6 아님) | `notebooks/4_optimal_expansion/4_1_Optimal_expansion_multisection.ipynb` |
| **Fig. 3c** | "Comparison of cells identified with different segmentation algorithms in an ROI (160 × 160 µm), using DAPI background" | **8가지 방법**: Xenium cell, Xenium nuc, MESMER, Clustermap, Cellpose, Watershed, Baysor conf=0, Baysor conf=0.8 | `_save_dapi_segmentation_comparison()` (파이프라인은 4가지) | `notebooks/5_segmentation_benchmark/run_segmentation.py` + 개별 Cellpose 노트북 |
| **Fig. 3d** | "Adjusted rand index (ARI) comparison of segmentation outputs (52 top performers) when applied to mouse brain section 2" | **52개 top performer** 간 ARI heatmap (315 config에서 선별). Staining-based vs read-based vs mixed 클러스터링 | `_compute_rand_index()` (파이프라인은 4가지 방법 간) | `notebooks/5_segmentation_benchmark/metrics.py` → `rand_idx()` |
| **Fig. 3e** | "Scatter plot of reads assigned versus negative marker purity for segmentation strategies applied to mouse brain section 2" | X축=Proportion of assigned reads, Y축=NMP. 7가지 방법 × expansion(0-14.9µm) × prior conf(0-0.99) 조합. Baysor=●, Xenium=✱, Cellpose=◆ 등 | `_save_reads_vs_nmp_scatter()` | `notebooks/3_techniques_comparison/3_4_negative_marker_purity_for_specificity.ipynb` |
| **Fig. 3f** | "UMAP from coprocessed cells using Baysor and Xenium's nuclear segmentation in a mouse brain ROI" | **Baysor BA2 P0.8 + Nuclei만** 2가지 co-processed. Cell type 색상 | `_save_umap()` (파이프라인은 모든 방법 포함, 3-panel) | `notebooks/5_segmentation_benchmark/5_1_Compare_Clustering on_different_segmentations.ipynb` → `UMAP_between_segmentations.pdf` |
| **Fig. 3g** | "Violin plot comparing cell counts segmented by Baysor versus Xenium nuclear segmentation methods" | **Baysor BA2 P0.8 vs Nuclei만** 2가지. Counts per cell 분포 | `_save_counts_violin()` (파이프라인은 모든 방법 포함) | `notebooks/5_segmentation_benchmark/5_1_Compare_Clustering on_different_segmentations.ipynb` → `counts_between_segmentations.pdf` |
| **Fig. 3h** | "Bar plot of cell counts per population using different segmentation strategies" | **Baysor vs Nuclei** 2가지의 cell type별 **절대 세포 수** (proportion 아님) | `_save_celltype_barplot()` (파이프라인은 모든 방법 포함) | `notebooks/5_segmentation_benchmark/5_1_Compare_Clustering on_different_segmentations.ipynb` → `frequencies_between_segmentations.pdf` |
| **ExtData Fig. 5a** | "Localization of regions of interest represented in Extended Data Fig. 5b and Fig. 3c" | 마우스 뇌 조직 위 ROI 위치 표시 | `_save_dapi_zoomed_roi()` | `notebooks/5_segmentation_benchmark/run_segmentation.py` |
| **ExtData Fig. 5b** | "Regions of interest representing the cells identified using different segmentation algorithms" | DAPI 배경 + 방법별 color-specific mask (160×160µm ROI). 8가지 방법 side-by-side | `_save_dapi_segmentation_comparison()` | `notebooks/5_segmentation_benchmark/run_segmentation.py` |
| **ExtData Fig. 5c** | "Heat map representing the segmentation metrics of all segmentation strategies described in Fig. 3d" | 52 top performer × 7 메트릭 (assigned_prop, n_cells, median_reads, p5_reads, median_genes, p5_genes, NMP) | `benchmark_metrics.csv` | `notebooks/5_segmentation_benchmark/metrics.py` |
| **ExtData Fig. 5d** | "Adjusted rand index (ARI) between the different outputs produced by combinations of segmentation algorithms... when applied to human breast sections" | Human breast 데이터에서의 ARI heatmap | `_compute_rand_index()` (데이터셋만 다름) | `notebooks/5_segmentation_benchmark/metrics.py` |
| **ExtData Fig. 5e** | "Scatter plot representing the number of reads assigned (x-axis) and the negative marker purity (y-axis) of different assessed segmentation strategies in human breast tumor samples" | Human breast에서의 NMP vs Assigned reads scatter | `_save_reads_vs_nmp_scatter()` (데이터셋만 다름) | `notebooks/3_techniques_comparison/3_4_negative_marker_purity_for_specificity.ipynb` |

**파이프라인 vs 논문 주요 차이점**:
- 논문 Fig. 3c/d/e: **8가지 방법 × 다양한 expansion/confidence 조합 (315 configs → 52 top performers)** → 파이프라인: **4가지 방법 (Nuclei, Cellpose, Expansion, Baysor)**
- 논문 Fig. 3f/g/h: **Baysor BA2 P0.8 vs Nuclei만 (2가지)** → 파이프라인: **모든 방법 포함 (확장된 비교)**
- 논문에는 MESMER, Clustermap, Watershed, Binning 포함 → 파이프라인에는 미구현

---

## 2. Pipeline Process Flow

```
Step 0: Nuclei .h5ad          Step 3: Cellpose .h5ad
Step 5: Expanded CSV          Baysor: 별도 실행
         ↓                              ↓
    [Step 6: Segmentation Benchmark]
         ↓
    ┌───── 6-1. Baysor Execution (선택) ────────────────────┐
    │  Transcript CSV → Baysor format                        │
    │  Optional: ROI crop / tiled processing                 │
    │  baysor run → segmentation.csv                         │
    │  → cell×gene matrix 구축                               │
    └────────────────────────────────────────────────────────┘
         ↓
    ┌───── 6-2. Load Segmentation Results ──────────────────┐
    │  각 방법별 AnnData 로드/생성                           │
    │  Nuclei: Step 0 h5ad                                   │
    │  Cellpose: Step 3 h5ad                                 │
    │  Expansion: Step 5 CSV → AnnData                       │
    │  Baysor: segmentation.csv → AnnData                    │
    └────────────────────────────────────────────────────────┘
         ↓
    ┌───── 6-3. Concatenate & Preprocess ───────────────────┐
    │  ad.concat([nuclei, cellpose, expansion, baysor])      │
    │  normalize → log1p → (HVG) → PCA → neighbors          │
    │  → UMAP → Leiden                                       │
    └────────────────────────────────────────────────────────┘
         ↓
    ┌───── 6-4. Annotation Transfer ────────────────────────┐
    │  Primary: Cluster-level crosstab majority vote         │
    │  Fallback: kNN in aligned PCA space (overlap < 1%)     │
    │  + Per-cluster consensus + Confidence score            │
    └────────────────────────────────────────────────────────┘
         ↓
    ┌───── 6-4b. Benchmark Mode 결정 ────────────────────────┐
    │  Dual mode: Baysor crop + subset_other_methods=true    │
    │    → A. Full tissue (Baysor 제외)                      │
    │    → B. Crop region (모든 방법 ROI 제한)               │
    │  Single mode: 모든 방법 전체 데이터                     │
    └────────────────────────────────────────────────────────┘
         ↓
    ┌───── 6-5. Benchmark Metrics ──────────────────────────┐
    │  기본: n_cells, median_reads/genes, p5_reads/genes     │
    │  조건부: assigned_prop, NMP_cells, NMP_reads           │
    │  sklearn: silhouette, calinski_harabasz, davies_bouldin│
    │  Rand Index: 방법 간 transcript 할당 일치도            │
    └────────────────────────────────────────────────────────┘
         ↓
    ┌───── 6-5b. Marker Gene Analysis ─────────────────────┐
    │  rank_genes_groups per segmentation method             │
    │  rank_genes_groups per cell type                       │
    │  DEG dotplot                                           │
    └────────────────────────────────────────────────────────┘
         ↓
    ┌───── 6-6. Visualizations ─────────────────────────────┐
    │  6-6a. UMAP per method (2~3-panel: method/cluster/type)│
    │  6-6b. Spatial scatter (방법별)                        │
    │  6-6c. Cell type barplot                               │
    │  6-6d. Counts violin                                   │
    │  6-6e. Spatial cell type map                           │
    │  6-6f. Reads vs NMP scatter (ExtData Fig. 5e)          │
    │  6-6g. DAPI segmentation comparison                    │
    │  6-6h. DAPI zoomed ROI                                 │
    │  6-6i. Cell boundary overlay                           │
    └────────────────────────────────────────────────────────┘
         ↓
├── <output_dir>/ (full_tissue 또는 single mode 출력)
│   ├── benchmark_combined.h5ad
│   ├── benchmark_metrics.csv
│   ├── rand_index_matrix.csv
│   ├── rand_index_heatmap.png
│   ├── umap_benchmark.png
│   ├── spatial_segmentation_map.png
│   ├── spatial_celltype_map.png
│   ├── celltype_frequency_barplot.png
│   ├── counts_violin_by_method.png
│   ├── reads_vs_nmp_scatter.png
│   ├── marker_genes_by_segmentation.png
│   ├── marker_genes_by_celltype.png
│   ├── deg_dotplot_by_segmentation.png
│   ├── dapi_segmentation_comparison.png
│   ├── dapi_zoomed_roi_comparison.png
│   ├── dapi_nucleus_boundaries.png
│   ├── qc_counts_genes_histogram.png
│   ├── hvg_selection_plot.png
│   ├── pca_scree_plot.png
│   ├── crop_comparison/ (Dual mode 시 crop 벤치마크 결과)
│   │   └── (위와 동일 구조)
│   └── baysor_run/ (Baysor 실행 결과)
│       └── crop/ (또는 tile_XXXX/ 또는 prep/)
└── step6_done.txt
```

---

## 3. Sub-step 상세 분석

### 3.1 Baysor Execution (`run_baysor()`, lines 484-569)

**Baysor란**: Baysor (Petukhov et al., Nat. Biotechnol. 2021)는 transcript 위치와 밀도를 기반으로 세포를 정의하는 알고리즘이다. DAPI 기반 방법과 달리, read density와 composition으로 세포 경계를 결정한다.

**실행 모드**:

| 모드 | 설명 | 용도 |
|------|------|------|
| **Single** | 전체 데이터셋에 Baysor 실행 | 소규모 데이터 |
| **Crop** | 밀도가 높은 ROI만 실행 | 대규모 데이터 (기본, ~5-10h) |
| **Tiled** | 전체를 타일로 분할하여 병렬 실행 | 미래 확장용 |

**Crop 모드 상세**:
```python
def find_densest_region(transcripts_df, size_um=2000):
    """2000×2000 µm 영역 중 가장 transcript 밀도가 높은 영역 탐색"""
    # Grid search로 최적 ROI 탐색
    # → crop 좌표 반환 (x_min, x_max, y_min, y_max)
```

**Baysor 명령** (현재 설정 기반):
```bash
baysor run \
    -s 15 \
    -m 10 \
    --prior-segmentation-confidence 0.8 \
    -o baysor_output/ \
    transcripts.csv \
    [prior_masks.tif]  # 선택: Cellpose masks를 prior로 사용
```

**설정값**:
```yaml
baysor:
  enabled: true
  dry_run: false                    # true=시뮬레이션만
  executable_path: "/path/to/baysor"
  params:
    scale: 15                       # 예상 세포 반경 (µm). 뇌 조직: ~15 µm
    min_molecules_per_cell: 10      # 최소 transcript/cell
    prior_segmentation_confidence: 0.8  # Prior mask 신뢰도
    save_polygons: true             # 세포 경계 polygon 저장
  crop:
    enabled: true
    size_um: 2000                   # ROI 크기 (µm)
    coords: null                    # null=자동 밀도 기반, 또는 [x_min, x_max, y_min, y_max]
    subset_other_methods: true      # 다른 방법도 같은 ROI로 제한
```

**scale 파라미터 참고**: 노트북에서 사용된 파일명 `baysor-82`의 82는 **pixel 단위** (~82px × 0.2125 µm/px ≈ 17 µm). 파이프라인은 µm 좌표를 직접 전달하므로 scale=15가 동등한 값이다. scale=50은 과도한 over-merging을 유발 (2.2M transcripts에서 295개 세포만 생성).

### 3.2 Load Segmentation Results (`load_transcripts_as_adata()`, lines 574-649)

각 방법별 데이터 로딩:

| 방법 | 소스 | 로딩 방식 |
|------|------|---------|
| Nuclei | Step 0 `*.h5ad` | 직접 로드 |
| Cellpose | Step 3 `resegmented.h5ad` | 직접 로드 |
| Expansion | Step 5 `expanded_transcripts.csv` | CSV → groupby → AnnData 생성 |
| Baysor | `baysor_run/segmentation.csv` | CSV → groupby → AnnData 생성 |

**CSV → AnnData 변환**:
```python
def load_transcripts_as_adata(csv_path, ...):
    df = pd.read_csv(csv_path)
    # 컬럼 자동 해석 (우선순위: closest_cell > in_cell > cell_col > cell_id > cell)
    # 유전자 컬럼: gene_col > feature_name > gene
    # cell_id > 0인 할당된 transcript만 사용 (숫자/문자열 모두 처리)
    df_assigned = df[df[actual_cell_col] > 0]  # 또는 != 'unassigned'/'noise'
    # Cell×Gene count matrix (crosstab 사용)
    counts = pd.crosstab(df_assigned[actual_cell_col], df_assigned[actual_gene_col])
    adata = ad.AnnData(counts)
    # Centroids 계산 (좌표 자동 감지: x_location > x_global_px > x)
    centroids = df_assigned.groupby(actual_cell_col)[[xc, yc]].mean()
    adata.obs[['x_centroid','y_centroid']] = centroids
    # Raw spots 저장 (proportion_of_assigned_reads + Rand Index용)
    adata.uns['spots'] = df[keep_cols].copy()
    return adata
```

### 3.3 Concatenate & Preprocess (`preprocess_benchmark()`, lines 831-890)

**통일된 전처리**:
```python
# 1. 모든 방법을 concatenate
combined = ad.concat(
    [nuclei, cellpose, expansion, baysor],
    keys=['nuclei', 'cellpose', 'expansion', 'baysor'],
    axis=0
)
combined.obs['segmentation'] = combined.obs.index.get_level_values(0)

# 2. Raw layer 저장
combined.layers['raw'] = combined.X.copy()

# 3. 전처리 파이프라인
sc.pp.filter_cells(combined, min_counts=40)
sc.pp.filter_cells(combined, min_genes=15)
combined.raw = combined  # scanpy .raw 속성 보존
combined.layers['raw'] = combined.X.copy()
sc.pp.normalize_total(combined, target_sum=100)
sc.pp.log1p(combined)

# 4. HVG (설정에 따라)
if config['benchmark']['preprocessing'].get('hvg', True):
    sc.pp.highly_variable_genes(combined,
        min_mean=0.3, max_mean=7, min_disp=-0.5)

# 5. Scaling (코드 기본값: False, config에서 true로 오버라이드)
if config['benchmark']['preprocessing'].get('scale', False):
    sc.pp.scale(combined)

# 6. PCA → Neighbors → UMAP → Leiden
sc.tl.pca(combined)
sc.pp.neighbors(combined, n_neighbors=15, n_pcs=0)  # 코드 기본값 15, config에서 16으로 설정
sc.tl.umap(combined, min_dist=0.1)
sc.tl.leiden(combined, resolution=1.0)

# 7. Auto-reduce resolution (클러스터 >40개면 해상도 절반으로 재클러스터링)
if n_clusters > 40:
    sc.tl.leiden(combined, resolution=resolution * 0.5)
```

**논문의 Best-Practice Preprocessing**:
```
(1) Library-size normalization (target_sum=100)
(2) Log transformation
(3) Scaling
(4) All principal components
(5) k-NN graph (k=16)
(6) Louvain clustering
```
파이프라인은 이 best-practice를 따르되:
- Louvain 대신 Leiden 사용 (실질적으로 동일한 알고리즘의 개선판)
- 코드 기본값 `scale=False`, `n_neighbors=15`이지만 config.yaml에서 `scale=true`, `n_neighbors=16`으로 오버라이드
- Leiden 클러스터 >40개 시 자동으로 해상도 절반 재클러스터링

### 3.4 Annotation Transfer (`annotate_by_majority_voting()`, lines 654-826)

**이중 전략 (Dual Strategy)**:

**전략 1 (Primary): Cluster-level Crosstab Majority Vote** (overlap_ratio > 1% 시):
```python
# 1. Target과 Reference 사이의 overlapping cells 매칭
#    (obs_names 또는 cell_id 컬럼으로 매칭 시도)
shared_cells = adata_target.obs_names.intersection(adata_reference.obs_names)

# 2. Crosstab: cluster × cell_type counts
overlap_labels = ref_labels.loc[shared_cells]
overlap_clusters = adata_target.obs.loc[shared_cells, cluster_key]
crosstab = pd.crosstab(overlap_clusters, overlap_labels)

# 3. 각 cluster에서 최빈 cell type 할당
cluster_map = {cluster: crosstab.loc[cluster].idxmax() for cluster in crosstab.index}
adata_target.obs['celltype_cluster'] = adata_target.obs[cluster_key].map(cluster_map)
adata_target.obs['celltype_majority'] = adata_target.obs['celltype_cluster']

# 4. Confidence score 계산 (클러스터 내 winner 비율)
adata_target.obs['celltype_confidence'] = ...  # crosstab.loc[cid].max() / total
```

**전략 2 (Fallback): kNN in Aligned PCA Space** (overlap < 1% 시):
```python
from sklearn.neighbors import NearestNeighbors
# Shared genes만으로 Reference PCA 구축 → Target 투영
# Reference statistics (mean/std)로 Target scaling
# kNN 검색 후 majority voting
# → celltype_majority + celltype_confidence + celltype_cluster (per-cluster consensus)
```

**Confidence Score**: 두 전략 모두 `celltype_confidence`를 계산함 (crosstab: winner 비율, kNN: k개 이웃 중 winner 비율). 단, 논문의 0.7 threshold filtering은 미적용.

### 3.5 Benchmark Metrics (`_run_benchmark_pipeline()`, lines 1752-2067)

**호출되는 메트릭** (모든 segmentation method별):
| 메트릭 | 함수 | 호출 라인 | 조건 | 결과 키 |
|--------|------|---------|------|---------|
| `n_cells` | inline `shape[0]` | 1879 | 무조건 | `{method}_n_cells` |
| `median_reads_cells` | `metrics.median_reads_cells()` | 1882 | `layers['raw']` 존재 시 | `{method}_median_reads` |
| `median_genes_cells` | `metrics.median_genes_cells()` | 1883 | `layers['raw']` 존재 시 | `{method}_median_genes` |
| `percentile_5th_reads_cells` | `metrics.percentile_5th_reads_cells()` | 1884 | `layers['raw']` 존재 시 | `{method}_p5_reads` |
| `percentile_5th_genes_cells` | `metrics.percentile_5th_genes_cells()` | 1885 | `layers['raw']` 존재 시 | `{method}_p5_genes` |
| `proportion_of_assigned_reads` | `metrics.proportion_of_assigned_reads()` | 1893 | `uns['spots']` 존재 시 | `{method}_assigned_prop` |
| `negative_marker_purity_cells` | `metrics.negative_marker_purity_cells()` | 1954 | scRNAseq reference 제공 시 | `{method}_nmp_cells` |
| `negative_marker_purity_reads` | `metrics.negative_marker_purity_reads()` | 1962 | scRNAseq reference 제공 시 | `{method}_nmp_reads` |
| `rand_idx` | `metrics.rand_idx()` | 1644 | `_compute_rand_index()` 내 | `rand_index_matrix.csv` |

**sklearn 클러스터링 품질 메트릭** (추가):
| 메트릭 | 함수 | 호출 라인 | 결과 키 |
|--------|------|---------|---------|
| `silhouette_score` | `sklearn.metrics.silhouette_score()` | 1902 | `{method}_silhouette` |
| `calinski_harabasz_score` | `sklearn.metrics.calinski_harabasz_score()` | 1903 | `{method}_calinski_harabasz` |
| `davies_bouldin_score` | `sklearn.metrics.davies_bouldin_score()` | 1904 | `{method}_davies_bouldin` |

**참고**: 모든 metrics.py 함수(8/8)가 호출됨. NMP 메트릭은 scRNAseq reference가 필요하며, proportion_of_assigned_reads는 `adata.uns['spots']`가 필요. 모든 메트릭은 `adata.layers['raw']`를 기대함.

### 3.6 Marker Gene Analysis (lines 1987-2034)

벤치마크 메트릭 이후, 방법별/세포유형별 marker gene 분석:

```python
# 1. Segmentation 방법별 DEG
sc.tl.rank_genes_groups(adata, 'segmentation')
# → marker_genes_by_segmentation.png

# 2. Cell type별 DEG (annotation 있을 경우)
sc.tl.rank_genes_groups(adata, 'celltype')
# → marker_genes_by_celltype.png

# 3. DEG dotplot (segmentation별)
sc.pl.dotplot(adata, var_names=top_genes, groupby='segmentation')
# → deg_dotplot_by_segmentation.png
```

### 3.7 Dual Benchmark Mode (`run_step6()`, lines 2072-2208)

Baysor crop 모드 사용 시, 공정한 비교를 위해 **두 개의 벤치마크**를 실행:

| 모드 | 포함 방법 | ROI | 출력 디렉토리 |
|------|---------|-----|-------------|
| **Full tissue** | Nuclei, Cellpose, Expansion (Baysor 제외) | 전체 조직 | `<output_dir>/` (루트) |
| **Crop comparison** | Nuclei, Cellpose, Expansion, Baysor (모두 crop) | 2000×2000µm ROI | `<output_dir>/crop_comparison/` |

Baysor crop이 비활성이면 **Single mode**로 모든 방법을 `<output_dir>/` 루트에 직접 출력 (서브디렉토리 없음).

### 3.8 Visualizations

**6-6a. UMAP per Method** (`_save_umap()`, lines 975-1005):
- 3-panel: segmentation method / Leiden cluster / cell type 색상
- 모든 방법이 동일 UMAP 공간에 투영
- → 논문 Fig. 3f 관련(확장) — 논문은 Baysor+Nuclei 2가지만, 파이프라인은 4가지
- 노트북: `5_1_Compare_Clustering on_different_segmentations.ipynb` → `UMAP_between_segmentations.pdf`

**6-6b. Spatial Scatter** (`_save_spatial_map()`, lines 1256-1305):
- 방법별 패널: 세포 위치 scatter plot
- 좌표 자동 감지 (x_centroid/y_centroid, x_location/y_location, x/y)
- Y축 반전 (이미지 좌표계)
- → 논문 직접 대응 없음 (파이프라인 자체 확장)

**6-6c. Cell Type Barplot** (`_save_celltype_barplot()`, lines 1309-1348):
- 방법별 cell type proportion 비교
- → 논문 Fig. 3h 관련(확장) — 논문은 Baysor vs Nuclei 2가지의 **절대 세포 수**, 파이프라인은 4가지의 proportion
- 노트북: `5_1_Compare_Clustering on_different_segmentations.ipynb` → `frequencies_between_segmentations.pdf`

**6-6d. Counts Violin** (`_save_counts_violin()`, lines 1352-1387):
- 방법별 total counts per cell 분포
- → 논문 Fig. 3g 관련(확장) — 논문은 Baysor vs Nuclei 2가지만, 파이프라인은 4가지
- 노트북: `5_1_Compare_Clustering on_different_segmentations.ipynb` → `counts_between_segmentations.pdf`

**6-6e. Spatial Cell Type Map** (`_save_spatial_celltype_map()`, lines 1529-1584):
- 방법별 패널: cell type 색상으로 공간 scatter
- 일관된 색상 매핑 (tab20)
- → 논문 직접 대응 없음 (파이프라인 자체 확장)

**6-6f. Reads vs NMP Scatter** (`_save_reads_vs_nmp_scatter()`, lines 1391-1525):
- → 논문 Fig. 3e / ExtData Fig. 5e 재현
- 논문 Fig. 3e: 7방법×expansion(0-14.9µm)×prior confidence(0-0.99) 조합, Baysor=●, Xenium=✱, Cellpose=◆ 등
- 파이프라인: 4가지 방법, 각각 단일 점
- X축: proportion of assigned reads (또는 median reads), Y축: NMP
- 방법별 고유 마커/색상: nuclei(○,#00BFFF), cellpose(△,#FF6347), expansion(□,#32CD32), baysor(◇,#FFD700)
- NMP-cells / NMP-reads 2-panel subplot
- NMP=0.8 기준선 표시
- 노트북: `3_4_negative_marker_purity_for_specificity.ipynb`

**6-6g. DAPI Comparison** (`_save_dapi_segmentation_comparison()`, lines 1036-1093):
- DAPI 이미지 위에 각 방법의 cell centroids overlay
- 방법별 색상: nuclei=#00BFFF, cellpose=#FF6347, expansion=#32CD32, baysor=#FFD700
- → 논문 Fig. 3c / ExtData Fig. 5b 관련(축소) — 논문 8가지 방법 → 파이프라인 4가지
- 노트북: `run_segmentation.py`

**6-6h. DAPI Zoomed ROI** (`_save_dapi_zoomed_roi()`, lines 1096-1179):
- Dense region 자동 선택 (median cell position 기반, roi_fraction=0.2)
- 첫 패널: DAPI only, 이후 패널: 방법별 ROI overlay
- → 논문 ExtData Fig. 5a-b 관련
- 노트북: `run_segmentation.py`

**6-6i. Boundary Overlay** (`_save_boundary_overlay()`, lines 1182-1252):
- nucleus_boundaries.parquet에서 polygon 로드
- 최대 500개 cell boundary 샘플링하여 DAPI 위에 표시
- → 논문 직접 대응 없음 (파이프라인 자체 확장)

---

## 4. Notebook vs Pipeline 구현 비교

### 4.1 원본 Notebook
| 노트북 | 위치 | Pipeline 함수 | 상태 |
|--------|------|--------------|------|
| `5_1_Compare_Clustering on_different_segmentations.ipynb` | `notebooks/5_segmentation_benchmark/` | `run_step6()` 전체 (UMAP, violin, barplot) | 구현됨 |
| `3_4_negative_marker_purity_for_specificity.ipynb` | `notebooks/3_techniques_comparison/` | `_save_reads_vs_nmp_scatter()` (NMP 분석) | 구현됨 |
| `3_6_Diffussion_on_resegmented_data.ipynb` | `notebooks/3_techniques_comparison/` | `_compute_rand_index()` (read 할당 비교) | 구현됨 |
| `run_segmentation.py` | `notebooks/5_segmentation_benchmark/` | `_save_dapi_segmentation_comparison()` (ROI 비교) | 구현됨 |
| `metrics.py` | `notebooks/5_segmentation_benchmark/` | `benchmark_metrics.csv` (모든 정량 메트릭) | 구현됨 |
| 개별 Cellpose 노트북 (`batch_segmentation_cellpose-*.ipynb`) | `notebooks/3_techniques_comparison/3_1_resegmentation_notebooks/` | Cellpose segmentation 결과 생성 | Step 3에서 구현 |

### 4.2 상세 비교

| 기능 | Notebook | Pipeline | 심각도 |
|------|----------|----------|--------|
| Baysor 실행 | 외부 실행 후 결과 로드 | 파이프라인 내 자동 실행 (crop/tile/full) | 개선 |
| 데이터 로딩 | 수동 경로 지정 | config 기반 자동 | 개선 |
| Concatenation | `ad.concat()` | `ad.concat()` | 동일 |
| 전처리 | HVG 포함 | HVG 설정 가능 (best-practice: scale+HVG) | OK |
| 클러스터링 | Louvain (실제 Leiden) | Leiden | 동일 |
| Annotation transfer | ID 매핑 + cluster consensus | Crosstab majority vote (primary) + kNN (fallback) + cluster consensus | 다름 (더 일반적) |
| NMP metrics | 호출됨 | 호출됨 (scRNAseq reference 필요) | OK |
| Rand Index | 호출됨 | 호출됨 (`_compute_rand_index()` 내) | OK |
| 5th percentile | 호출됨 | 호출됨 | OK |
| Clustering quality | - | sklearn: silhouette, calinski_harabasz, davies_bouldin | 개선 |
| Marker gene analysis | - | `rank_genes_groups` per method/celltype + dotplot | 개선 |
| Dual benchmark | - | Full tissue + Crop region 분리 비교 | 개선 |
| Confidence score | 계산됨 (0.7 threshold로 필터링) | 계산됨 (`celltype_confidence`), **threshold 필터링 미적용** | LOW |
| UMAP | 방법별 색상 | 방법별 + cluster + cell type (3-panel) | OK (개선) |
| Spatial plot | Cell type 색상 | 방법별 + cell type 별도 패널 | OK |
| Publication quality | PDF, dpi=500 | PNG, dpi=150 | LOW |

### 4.3 남은 차이점

1. **Confidence score threshold**: 파이프라인에서 `celltype_confidence`를 계산하지만 (crosstab: winner 비율, kNN: neighbor winner 비율), 노트북의 0.7 threshold 기반 불확실 세포 필터링은 미적용.
2. **Publication quality**: 노트북은 PDF/dpi=500, 파이프라인은 PNG/dpi=150-200 (미관상 차이만, 분석에 영향 없음).
3. **`number_of_cells()` 함수**: metrics.py에 정의되어 있으나 직접 `shape[0]`으로 계산 (기능적으로 동일).
4. **코드 기본값 vs Config**: `scale`(코드 기본값=False, config=True), `n_neighbors`(코드 기본값=15, config=16). Config 설정이 논문 best-practice에 맞춰져 있음.

---

## 5. 전체 시각화 및 데이터 출력 목록

### 5.1 QC 시각화 (전처리 단계)
| # | 파일명 | 유형 | 분석법 |
|---|-------|------|--------|
| 1 | `qc_counts_genes_histogram.png` | Histogram | Cell counts/genes 분포 + QC threshold 라인 |
| 2 | `hvg_selection_plot.png` | Scatter | HVG 선택 결과 |
| 3 | `pca_scree_plot.png` | Line | PCA variance explained (cumulative) |

### 5.2 메인 벤치마크 시각화
| # | 파일명 | 유형 | 논문 Figure | 분석법 |
|---|-------|------|-----------|--------|
| 4 | `umap_benchmark.png` | UMAP (3-panel) | Fig 3f | 방법/cluster/cell type 비교 |
| 5 | `spatial_segmentation_map.png` | Spatial | - | 방법별 세포 공간 분포 |
| 6 | `celltype_frequency_barplot.png` | Barplot | Fig 3h | 방법별 cell type 구성 |
| 7 | `counts_violin_by_method.png` | Violin | Fig 3g | 방법별 counts 분포 |
| 8 | `spatial_celltype_map.png` | Spatial | - | Cell type별 공간 분포 |

### 5.2b NMP 시각화
| # | 파일명 | 유형 | 논문 Figure | 분석법 |
|---|-------|------|-----------|--------|
| 9 | `reads_vs_nmp_scatter.png` | Scatter (2-panel) | Ext Fig 5e | 방법별 Assigned reads vs NMP |

### 5.3 DAPI Morphology 시각화
| # | 파일명 | 유형 | 논문 Figure | 분석법 |
|---|-------|------|-----------|--------|
| 9 | `dapi_segmentation_comparison.png` | Multi-panel | Fig 3c | 방법별 segmentation 비교 |
| 10 | `dapi_zoomed_roi_comparison.png` | Zoomed | Ext Fig 5a-b | ROI 확대 비교 |
| 11 | `dapi_nucleus_boundaries.png` | Overlay | - | DAPI + 핵 경계 오버레이 |

### 5.4 Marker Gene 분석 시각화
| # | 파일명 | 유형 | 분석법 |
|---|-------|------|--------|
| 12 | `marker_genes_by_segmentation.png` | Rank plot | 방법별 marker genes |
| 13 | `marker_genes_by_celltype.png` | Rank plot | 세포유형별 marker genes |
| 14 | `deg_dotplot_by_segmentation.png` | Dotplot | 방법별 DEG 시각화 |

### 5.5 정량적 데이터 출력
| # | 파일명 | 유형 | 논문 Figure | 내용 |
|---|-------|------|-----------|------|
| 15 | `benchmark_metrics.csv` | CSV | Ext Fig 5c | 방법별 11+ 메트릭 요약 |
| 16 | `rand_index_matrix.csv` | CSV | Fig 3d | 방법 간 ARI 매트릭스 |
| 17 | `rand_index_heatmap.png` | Heatmap | Fig 3d | ARI 매트릭스 히트맵 시각화 (YlOrRd colormap) |
| 18 | `reads_vs_nmp_scatter.png` | Scatter | Ext Fig 5e | NMP vs Assigned reads (방법별 마커/색상) |
| 19 | `benchmark_combined.h5ad` | HDF5 | - | 전체 concatenated AnnData (캐시) |

---

## 6. Input / Output 상세

### 6.1 Input
| 파일 | 형식 | 설명 |
|------|------|------|
| Step 0 `*.h5ad` | AnnData | Nuclei segmentation (기준) |
| Step 3 `resegmented.h5ad` | AnnData | Cellpose segmentation |
| Step 5 `expanded_transcripts.csv` | CSV | Expansion 데이터 |
| Xenium `transcripts.csv` | CSV | Baysor 입력용 |
| Step 3 `masks.tif` (선택) | TIFF | Baysor prior segmentation |
| scRNAseq reference (선택) | h5ad | Annotation transfer, NMP |

### 6.2 Output
| 파일 | 형식 | 설명 |
|------|------|------|
| `<output_dir>/` | Directory | 모든 벤치마크 결과 (full_tissue 또는 single) |
| `crop_comparison/` | Directory | Crop 벤치마크 결과 (Dual mode 시) |
| `baysor_run/` | Directory | Baysor 실행 결과 (crop/tile/prep) |
| `benchmark_metrics.csv` | CSV | 방법별 정량 메트릭 (11+ metrics) |
| `benchmark_combined.h5ad` | HDF5 | Concatenated AnnData (캐시) |
| 16개+ 시각화 파일 | PNG | 비교 플롯 (UMAP, spatial, DAPI, NMP 등) |
| `step6_done.txt` | Marker | 완료 표시 |

---

## 7. 관련 설정값 정리

```yaml
baysor:
  enabled: true
  dry_run: false
  executable_path: "/path/to/baysor"
  params:
    scale: 15                       # 예상 세포 반경 (µm). 뇌: ~15 µm
    min_molecules_per_cell: 10
    prior_segmentation_confidence: 0.8
    save_polygons: true
  crop:
    enabled: true
    size_um: 2000
    coords: null                    # null=자동 밀도 기반
    subset_other_methods: true      # Dual benchmark mode 활성화
  tiling:
    enabled: false                  # 미래 확장용

benchmark:
  preprocessing:
    target_sum: 100             # Library-size normalization target (코드 기본값: 100)
    min_counts: 40              # 최소 counts/cell (QC) (코드 기본값: 40)
    min_genes: 15               # 최소 genes/cell (QC) (코드 기본값: 15)
    n_neighbors: 16             # k-NN (config 설정값, 코드 기본값: 15)
    n_pcs: 0                    # 0 = all PCs (코드 기본값: 0)
    umap_min_dist: 0.1
    resolution: 1.0             # Leiden 해상도 (코드 기본값: 1.0)
    scale: true                 # Scaling (config 설정값, 코드 기본값: false)
    hvg: true                   # HVG 선택 (min_mean=0.3, max_mean=7, min_disp=-0.5)
  reference_adata: null          # 비교용 reference h5ad (annotation transfer용)
  ref_celltype_key: "celltype"  # Reference cell type 컬럼
  sc_reference_path: null        # scRNAseq (NMP 메트릭용)
```

### 7.1 파라미터 조정 가이드

| 파라미터 | 현재값 | 범위/옵션 | 영향도 | 조정 가이드 |
|----------|--------|----------|--------|------------|
| **Baysor 파라미터** | | | | |
| `baysor.enabled` | `true` | true/false | - | Baysor 실행 여부. false면 nuclei/cellpose/expansion만 비교 |
| `baysor.dry_run` | `false` | true/false | - | true=시뮬레이션만(실제 Baysor 미실행). 설정 검증용 |
| `baysor.params.scale` | `15` | 5-30 µm | **TOP 2** | **핵심 파라미터**. 예상 세포 반경(µm). 뇌 조직=15µm. 과소(5)=over-segmentation, 과대(50)=over-merging(295개 세포만 생성). 조직 유형별 세포 크기에 맞춤 |
| `baysor.params.min_molecules_per_cell` | `10` | 5-50 | MEDIUM | 세포당 최소 transcript. 높이면 noise 세포 제거, 낮으면 작은 세포 보존 |
| `baysor.params.prior_segmentation_confidence` | `0.8` | 0.0-0.99 | **HIGH** | Prior mask(Cellpose) 신뢰도. 0=무시(순수 density), 0.8=강한 prior(논문 최적), 0.99=prior 거의 그대로 |
| `baysor.crop.enabled` | `true` | true/false | MEDIUM | Crop 모드. true=밀도 높은 ROI만 분석(5-10h), false=전체 조직(3일+) |
| `baysor.crop.size_um` | `2000` | 500-5000 µm | MEDIUM | Crop ROI 크기. 2000=4mm² (전체의 ~4%). 작으면 빠르지만 대표성 감소 |
| `baysor.crop.coords` | `null` | null/[x_min,x_max,y_min,y_max] | LOW | null=자동 밀도 감지, 좌표=특정 ROI 지정 |
| `baysor.crop.subset_other_methods` | `true` | true/false | MEDIUM | true=Dual benchmark(crop+full 두 번), false=전체만. 공정한 비교를 위해 true 권장 |
| **벤치마크 전처리** | | | | |
| `benchmark.preprocessing.target_sum` | `100` | 100/1000/10000 | MEDIUM | Normalization target. 논문 best-practice: 100. scRNAseq 기본(10000)보다 Xenium에 최적 |
| `benchmark.preprocessing.min_counts` | `40` | 10-100 | MEDIUM | 벤치마크 QC. Step 0보다 엄격(40 vs 10). 낮추면 더 많은 세포 포함 |
| `benchmark.preprocessing.min_genes` | `15` | 3-30 | MEDIUM | 벤치마크 QC. Step 0보다 엄격(15 vs 3) |
| `benchmark.preprocessing.n_neighbors` | `16` | 10-50 | MEDIUM | k-NN 그래프. 논문: 16이 최적. 코드 기본값: 15. config에서 16으로 설정 |
| `benchmark.preprocessing.scale` | `true` | true/false | **HIGH** | 유전자 간 분산 표준화. 논문 best-practice: true. 코드 기본값: false. config에서 true로 설정 |
| `benchmark.preprocessing.hvg` | `true` | true/false | MEDIUM | HVG 선택. Xenium 패널은 이미 curated이므로 영향 제한적. true가 약간 더 나은 결과 |
| `benchmark.preprocessing.n_pcs` | `0` | 0/15-50 | LOW | PCA 컴포넌트. 0=전체 사용(논문 권장). Xenium은 유전자 수가 적어 차원 축소 불필요 |
| `benchmark.preprocessing.resolution` | `1.0` | 0.5-2.0 | MEDIUM | Leiden 해상도. 세포 유형 수에 따라 조정 |
| **Reference 설정** | | | | |
| `benchmark.reference_adata` | `null` | 파일 경로/null | MEDIUM | Annotation transfer용 reference. null이면 Leiden 클러스터만 사용 |
| `benchmark.sc_reference_path` | `null` | 파일 경로/null | MEDIUM | NMP 메트릭용 scRNAseq. null이면 NMP 계산 건너뜀 |

**튜닝 팁**:
- `baysor.params.scale`이 가장 중요: **반드시 조직의 세포 크기에 맞춰야 함**. 뇌=15µm, 간=20µm. 잘못 설정하면 극단적 결과
- `prior_segmentation_confidence=0.8`은 논문에서 검증된 최적값. Cellpose segmentation이 불안정하면 낮춤 (0.5), 정확하면 높임 (0.9)
- 전처리 `scale=true`와 `target_sum=100`은 **변경하지 않는 것을 권장** (논문 best-practice)
- Baysor crop 모드가 기본. 전체 조직 Baysor가 필요하면 `crop.enabled: false`로 설정하되 실행 시간(3일+) 주의

---

## 8. 핵심 로직의 의미

### 8.1 Baysor scale = 15 (µm)
Baysor의 `-s` 파라미터로, **예상 세포 반경** (µm 단위):
- **scale=15**: 뇌 조직의 평균 세포 반경 (~15 µm). 파이프라인은 µm 좌표를 직접 전달.
- **이전 값 scale=50**: 과도한 over-merging 유발 (2.2M transcripts → 295개 세포만 생성)
- **노트북 참조**: `baysor-82` 파일명은 pixel 단위 (82px × 0.2125 µm/px ≈ 17 µm)

### 8.2 prior_segmentation_confidence = 0.8
Baysor에서 prior segmentation (Cellpose masks)을 얼마나 신뢰할지:
- **0.0**: Prior 무시 (순수 read-density 기반)
- **0.8**: Prior를 강하게 참고하되 read density로 조정
- **0.99**: Prior 거의 그대로 사용
- 논문: "Baysor combined with Xenium's nuclei segmentation (BA2 P0.8)" → 최고 성능

### 8.3 Crop ROI Size = 2000 µm
- 2000 × 2000 µm = 4 mm² 영역
- 전체 Xenium 이미지 (~1 cm²)의 약 4%
- Baysor 전체 실행: 3일+ → Crop: 5-10시간

### 8.4 subset_other_methods
Baysor를 crop 모드로 실행할 때, 다른 방법 (nuclei, cellpose, expansion)도 같은 ROI로 제한하여 **공정한 비교**를 보장.

### 8.5 Benchmark Preprocessing 파라미터
논문 Fig. 4b-c에서 식별된 best-practice:
- **target_sum=100**: Library-size normalization (SCTransform보다 나음)
- **scale=true**: 유전자 간 분산 표준화 (클러스터링 품질 향상)
- **n_neighbors=16**: 16 이상에서 안정적 그래프
- **n_pcs=0 (all)**: 모든 PC 사용 (Xenium 패널은 유전자 수가 적으므로)

---

## 10. 시각화 상세 분석 가이드

### 10.1 UMAP Benchmark (`umap_benchmark.png`) → 논문 Fig. 3f 관련 (확장)

**논문 위치**: Fig. 3f (p.818) — "UMAP from coprocessed cells using Baysor and Xenium's nuclear segmentation in a mouse brain ROI"
**노트북 출처**: `notebooks/5_segmentation_benchmark/5_1_Compare_Clustering on_different_segmentations.ipynb` → `UMAP_between_segmentations.pdf`

> **논문 vs 파이프라인 차이**: 논문 Fig. 3f는 **Baysor BA2 P0.8 + Nuclei 2가지만** co-process한 단일 UMAP (cell type 색상). 파이프라인은 **모든 방법 (Nuclei, Cellpose, Expansion, Baysor)을 함께** co-process하여 3-panel UMAP (방법별/클러스터별/cell type별)으로 확장.

**무엇을 봐야 하는가**:
```
    ┌── Panel 1: Method 색상 ──┐  ┌── Panel 2: Cluster 색상 ──┐  ┌── Panel 3: Cell type 색상 ──┐
    │  ○○○   ●●●   □□□        │  │  111   222   333            │  │  AAA   BBB   CCC             │
    │ ○○○○  ●●●●  □□□□        │  │ 1111  2222  3333            │  │ AAAA  BBBB  CCCC             │
    │  ○○○   ●●●   □□□        │  │  111   222   333            │  │  AAA   BBB   CCC             │
    │                          │  │                             │  │                               │
    │  ▲▲▲                     │  │  444   555                  │  │  DDD   EEE                    │
    │ ▲▲▲▲                     │  │ 4444  5555                  │  │ DDDD  EEEE                    │
    └──────────────────────────┘  └─────────────────────────────┘  └───────────────────────────────┘

    Panel 1: ○=Nuclei ●=Cellpose □=Expansion ▲=Baysor
    Panel 2: 1-5 = Leiden clusters
    Panel 3: A-E = Annotated cell types (e.g., Oligo, Neuron, Astro, ...)

    핵심 확인 사항:
    ① 모든 방법이 같은 클러스터에 혼재 (intermingled) → 방법 간 일관성
    ② 특정 방법이 별도 클러스터 형성 → 해당 방법의 systematic bias
    ③ 클러스터 분리도 → 전처리 품질
    ④ Cell type 색상 패널에서 클러스터와 cell type의 일치도
```

**해석법**:
- **잘 혼재된 UMAP (Panel 1)**: 모든 segmentation 방법의 세포가 같은 UMAP 영역에 분포 → 방법에 관계없이 동일한 세포 유형을 포착
- **분리된 영역**: 특정 방법의 세포만 있는 영역 → batch effect 또는 방법별 systematic bias
  - Expansion이 별도로 몰려있으면: 과도한 확장으로 인한 transcript 혼합 → 새로운 발현 패턴
  - Baysor만 별도: Baysor의 세포 정의가 다른 방법과 근본적으로 다름
- **논문 결과** (Fig. 3f): Baysor BA2 P0.8로 정의된 세포가 Nuclei와 co-process 시 동일 UMAP 공간에서 잘 분리된 클러스터 형성 → "the identified cellular populations were the same across both segmentation strategies" (p.817)

---

### 10.2 DAPI Segmentation Comparison (`dapi_segmentation_comparison.png`) → 논문 Fig. 3c / ExtData Fig. 5b ★핵심★

**논문 위치**: Fig. 3c (p.818) — "Comparison of cells identified with different segmentation algorithms in an ROI (160 × 160 µm), using DAPI background"
**노트북 출처**: `notebooks/5_segmentation_benchmark/run_segmentation.py` + 개별 Cellpose 노트북

> **논문 vs 파이프라인 차이**: 논문 Fig. 3c는 **8가지 방법/변형** (Xenium cell, Xenium nuc, MESMER, Clustermap, Cellpose, Watershed, Baysor conf=0, Baysor conf=0.8)의 2×4 grid 비교. 파이프라인은 구현된 **4가지 방법** (Nuclei, Cellpose, Expansion, Baysor)만 비교. MESMER, Clustermap, Watershed는 파이프라인에 미구현.

**무엇을 봐야 하는가**:
```
    논문 Fig. 3c (8 panels, 2×4 grid):
    ┌─ Xenium(cell)─┐  ┌─ Xenium(nuc) ─┐  ┌─── MESMER ────┐  ┌── Clustermap ─┐
    │  큰 mask       │  │  작은 핵 mask  │  │  중간 mask     │  │  read 기반    │
    └────────────────┘  └────────────────┘  └────────────────┘  └────────────────┘
    ┌── Cellpose ───┐  ┌── Watershed ──┐  ┌─ Baysor(c=0) ─┐  ┌─ Baysor(c=0.8)┐
    │  핵 정밀 탐지  │  │  경계 기반    │  │  순수 density  │  │  Prior+density │
    └────────────────┘  └────────────────┘  └────────────────┘  └────────────────┘

    파이프라인 출력 (4 panels):
    ┌── Nuclei ──┐  ┌─ Cellpose ─┐  ┌─ Expansion ┐  ┌── Baysor ──┐
    │ ╭─╮  ╭─╮  │  │ ╭──╮ ╭──╮  │  │ ╭────╮╭───╮│  │ ╭──╮╭──╮  │
    │ │●│  │●│  │  │ │● │ │● │  │  │ │●   ││ ● ││  │ │● ││● │  │
    │ ╰─╯  ╰─╯  │  │ ╰──╯ ╰──╯  │  │ ╰────╯╰───╯│  │ ╰──╯╰──╯  │
    │  ← 미할당  │  │  ← 핵 기반 │  │  ← 과도확장 │  │  ← 적응적  │
    └────────────┘  └─────────────┘  └─────────────┘  └────────────┘

    각 패널 = 같은 ROI에서의 다른 segmentation 방법
    DAPI (회색 배경) + 세포별 color-specific mask

    핵심 비교:
    ① Nuclei (=Xenium nuc): 핵만 탐지 → 작은 mask, 많은 미할당 reads
    ② Cellpose: nuclei 정밀 탐지 → Xenium default보다 정확한 핵 segmentation
    ③ Expansion: Step 5 최적 확장 → 가장 큰 coverage (논문의 1-14.9µm expansion에 해당)
    ④ Baysor: read density 기반 → 불규칙하지만 생물학적으로 적응적
```

**해석법 (각 방법별 강점/약점)**:

| 방법 | 강점 | 약점 | 시각적 특징 |
|------|------|------|-----------|
| **Nuclei** | 높은 순도, 정확한 경계 | reads 손실 많음 | 작은 타이트한 mask |
| **Cellpose** | nuclei 정밀 탐지 | 확장이 균일 (생물학적으로 부자연스러움) | 균일한 ring 형태 |
| **Expansion** | 높은 coverage | misassignment 위험 | 매우 큰, 겹치는 mask |
| **Baysor** | 적응적 경계 | 실행 시간 김, 의존성 | 불규칙한 자연스러운 경계 |

- **논문 결과**: "Baysor + Cellpose nuclei (BA2 P0.8)" 조합이 최고 성능

---

### 10.3 Cell Type Frequency Barplot (`celltype_frequency_barplot.png`) → 논문 Fig. 3h 관련 (확장)

**논문 위치**: Fig. 3h (p.818) — "Bar plot of cell counts per population using different segmentation strategies"
**노트북 출처**: `notebooks/5_segmentation_benchmark/5_1_Compare_Clustering on_different_segmentations.ipynb` → `frequencies_between_segmentations.pdf`

> **논문 vs 파이프라인 차이**: 논문 Fig. 3h는 **Baysor BA2 P0.8 vs Nuclei 2가지만** 비교하며, Y축은 **절대 세포 수** (No. of cells). 파이프라인은 **모든 방법 (Nuclei, Cellpose, Expansion, Baysor) 포함**하며, Y축은 cell type proportion.

**무엇을 봐야 하는가**:
```
    논문 Fig. 3h (2가지만):
    No. of cells
    3000 │  ████░░░░
    2000 │  ████░░░░  ████░░░░
    1000 │  ████░░░░  ████░░░░  ████░░░░  ...
       0 ├──Astro──Oligo──Neuron──Micro── ...──→
         ████ = Baysor (파랑)   ░░░░ = Nuclei (주황)

    파이프라인 출력 (4가지):
    Cell Count (%)
    40 ├────────────────────────────
    30 │  ████
    20 │  ████  ████  ████
    10 │  ████  ████  ████  ████  ████
     0 ├──Oligo─Neuron─Astro─Micro─OPC──Endo──→
       ████ = Nuclei    ░░░░ = Cellpose
       ▓▓▓▓ = Expansion ╠╣╠╣ = Baysor

    핵심 비교:
    ① 세포 유형 비율이 방법 간 일관적인지
    ② 특정 방법에서 사라지거나 과대 표현되는 세포 유형
    ③ scRNAseq reference 비율과의 일치도
```

**해석법**:
- **일관된 비율**: 모든 방법에서 비슷한 세포 유형 구성 → robust segmentation
- **Expansion에서 특정 유형 증가**: 과도한 확장으로 작은 세포(microglia)가 큰 세포(neuron)에 흡수 → microglia 감소, neuron 증가
- **Baysor에서 더 많은 세포**: Baysor가 nuclei 외 추가 세포를 탐지
- **논문 결과** (p.817): "cells defined by Baysor had a higher count per cell... the identified cellular populations were the same across both segmentation strategies (Fig. 3f-h), with mostly mild differences in cell-type abundance"

---

### 10.4 Counts Violin (`counts_violin_by_method.png`) → 논문 Fig. 3g 관련 (확장)

**논문 위치**: Fig. 3g (p.818) — "Violin plot comparing cell counts segmented by Baysor versus Xenium nuclear segmentation methods"
**노트북 출처**: `notebooks/5_segmentation_benchmark/5_1_Compare_Clustering on_different_segmentations.ipynb` → `counts_between_segmentations.pdf`

> **논문 vs 파이프라인 차이**: 논문 Fig. 3g는 **Baysor BA2 P0.8 vs Nuclei 2가지만** violin 비교. 파이프라인은 **모든 방법 포함**.

**무엇을 봐야 하는가**:
```
    논문 Fig. 3g (2가지):             파이프라인 (4가지):
    Counts/cell                       Counts/cell
    1000 ┤                            500 ┤
         │  ╭─╮                           │              ╭─╮
     500 ┤  │ │  ╭─╮                  400 ┤              │ │
         │  │━│  │ │                      │     ╭─╮      │ │
         │  │ │  │━│                  300 ┤     │ │  ╭─╮ │ │  ╭─╮
         │  │ │  │ │                      │  ╭─╮│ │  │ │ │ │  │ │
         │  ╰─╯  ╰─╯                 200 ┤  │ ││━│  │━│ │━│  │━│
       0 ├─Baysor──Nuclei──→              │  │━││ │  │ │ │ │  │ │
                                      100 ┤  │ │╰─╯  ╰─╯ ╰─╯  ╰─╯
    Baysor > Nuclei (높은 counts)        0 ├──Nuc─CP──Exp──Bay──→

    ① Median (가로선): 방법별 중앙 reads/cell
    ② Violin 폭: 넓으면 분포가 산포, 좁으면 균일
    ③ 상위 꼬리: 매우 많은 reads를 가진 세포 (큰 세포 또는 doublet)
    ④ 하위 꼬리: 저품질 세포 (QC 경계)
```

**해석법**:
- **Nuclei**: 가장 낮은 median → 핵 내 reads만 포함
- **Expansion**: 가장 높은 median → 가장 많은 reads 할당 (but misassignment 위험)
- **Cellpose**: Nuclei보다 약간 높음 → 확장된 mask로 추가 reads
- **Baysor**: Cellpose+Baysor 조합이 높은 median + 낮은 편차
- **논문 결과** (p.817): "cells defined by Baysor had a higher count per cell" — Baysor가 cytoplasmic reads를 효과적으로 할당하므로 Nuclei보다 counts 높음
- **이상적**: 높은 median + 좁은 분포 = 일관된 품질의 세포 탐지

---

### 10.5 Rand Index Heatmap (`rand_index_heatmap.png`) → 논문 Fig. 3d ★핵심★

**논문 위치**: Fig. 3d (p.818) — "Adjusted rand index (ARI) comparison of segmentation outputs (52 top performers) when applied to one of the mouse brain samples profiled (mouse brain section 2)"
**노트북 출처**: `notebooks/5_segmentation_benchmark/metrics.py` → `rand_idx()`

> **논문 vs 파이프라인 차이**: 논문 Fig. 3d는 **315개 config에서 선별한 52 top performers** 간의 대규모 ARI heatmap. Staining-based (DAPI), read-based (Baysor), mixed (Baysor+prior) 방법들이 자연스럽게 클러스터링됨. 파이프라인은 **4가지 방법 간** 소규모 ARI 매트릭스. ExtData Fig. 5d는 human breast 데이터에서 동일 분석.

**무엇을 봐야 하는가**:
```
    논문 Fig. 3d (52×52 대규모):              파이프라인 (4×4 소규모):
    ┌─────────────────────┐                           Nuc  CP   Exp  Bay
    │ BA2        ████████ │  ← Baysor 계열     Nuc   1.0  0.72 0.45 0.68
    │ BA2 P0.8   ████████ │    자체적으로      CP    0.72 1.0  0.58 0.82
    │ BA3        ████████ │    높은 ARI        Exp   0.45 0.58 1.0  0.55
    │ ...        ........ │                    Bay   0.68 0.82 0.55 1.0
    │ Xenium nuc ████████ │  ← DAPI 계열
    │ CPn        ████████ │    자체적으로
    │ Mesmer     ████████ │    높은 ARI
    │ ...        ........ │
    │ CM r20 s1  ████████ │  ← Expansion 계열
    │ CM r30 s5  ████████ │
    └─────────────────────┘
    색상: 어두움 = 높은 ARI (유사)
          밝음 = 낮은 ARI (상이)

    핵심 패턴 (논문):
    ① DAPI-based 방법끼리 높은 ARI (같은 staining 기반)
    ② Baysor 계열끼리 높은 ARI (같은 read-density 기반)
    ③ DAPI vs Baysor: 중간 ARI (다른 접근)
    ④ Expansion 크기에 따라 ARI 변화
```

**해석법**:
- **ARI > 0.8**: 두 방법이 거의 동일한 transcript 할당
- **ARI 0.5-0.8**: 대부분 일치하지만 일부 차이 존재
- **ARI < 0.5**: 상당히 다른 segmentation → 방법 선택이 결과에 큰 영향
- **같은 카테고리 내 높은 ARI**: DAPI-based끼리, Baysor 계열끼리 자연 클러스터링
- **Expansion이 낮은 ARI**: 과도한 확장으로 인해 다른 방법과 크게 다른 할당 패턴

**논문 결과** (p.817):
- "We next identified groups of strategies that performed similarly... Staining-based strategies using DAPI generated similar outputs, with cell expansion being the force driving their differences"
- "Baysor-based, Clustermap-based and binning strategies clustered according to method"

---

### 10.6 Benchmark Metrics Table (`benchmark_metrics.csv`) → 논문 ExtData Fig. 5c

**논문 위치**: ExtData Fig. 5c (p.25/830) — "Heat map representing the segmentation metrics of all segmentation strategies described in Fig. 3d"
**노트북 출처**: `notebooks/5_segmentation_benchmark/metrics.py` (median_reads_cells, median_genes_cells, percentile_5th_*, proportion_of_assigned_reads, negative_marker_purity_*) + `notebooks/3_techniques_comparison/3_4_negative_marker_purity_for_specificity.ipynb` (NMP 분석)

> **논문 vs 파이프라인 차이**: 논문 ExtData Fig. 5c는 **52 top performers × 7 메트릭** (assigned_prop, n_cells, median_reads, p5_reads, median_genes, p5_genes, NMP) heatmap. 파이프라인은 **4가지 방법 × 11+ 메트릭** (기본 7 + sklearn clustering quality 3 + Rand Index) CSV 테이블.

**확인해야 할 메트릭**:
```
    Method      n_cells  median_reads  median_genes  NMP    assigned%
    Nuclei      45,230   142           38            0.92   76.8%
    Cellpose    48,100   168           42            0.88   82.3%
    Expansion   52,500   215           48            0.75   95.1%
    Baysor      47,800   185           44            0.90   85.6%

    ★ NMP vs Assigned Reads 트레이드오프 분석 ★
    높은 NMP + 높은 Assigned% = 최고의 segmentation
    높은 NMP + 낮은 Assigned% = 보수적 (reads 손실)
    낮은 NMP + 높은 Assigned% = 과도 확장 (reads 오염)
```

**해석법**:
- **n_cells**: 탐지된 세포 수. 너무 많으면 over-segmentation, 너무 적으면 under-segmentation
- **median_reads**: 높을수록 정보 풍부. 하지만 NMP와 함께 봐야 함
- **NMP > 0.85 + Assigned > 80%**: 이상적인 segmentation
- **NMP < 0.7**: segmentation 경계 불량 → 이웃 세포 reads 혼입
- **Assigned > 95%**: 과도한 확장 경고
- 논문 Fig. 3e에서 NMP (y축) vs Assigned reads (x축) scatter plot으로 트레이드오프 시각화

---

### 10.7 Spatial Cell Type Map (`spatial_celltype_map.png`)

**논문 위치**: 직접 대응 없음 (파이프라인 자체 확장 — 논문에서는 Fig. 3f UMAP으로 cell type 분포를 보여주지만, 공간 scatter map은 별도로 제공하지 않음)
**노트북 출처**: `notebooks/5_segmentation_benchmark/5_1_Compare_Clustering on_different_segmentations.ipynb` → `xp.map_of_clusters()` (유사 기능)

**무엇을 봐야 하는가**:
```
    ┌── Nuclei ──────┐  ┌── Cellpose ────┐  ┌── Baysor ──────┐
    │ ●●●  ○○○  ▲▲▲  │  │ ●●●  ○○○  ▲▲▲  │  │ ●●●  ○○○  ▲▲▲  │
    │ ●●   ○○   ▲▲   │  │ ●●●  ○○○  ▲▲▲  │  │ ●●●  ○○○  ▲▲   │
    │      □□□       │  │  □□□□□□       │  │  □□□□□□        │
    │ □□□□□□□        │  │ □□□□□□□       │  │ □□□□□□□□       │
    └────────────────┘  └────────────────┘  └────────────────┘

    각 색상 = 세포 유형 (tab20 색상 매핑)
    패널별로 같은 ROI에서의 다른 segmentation

    비교 포인트:
    ① 조직 구조가 모든 방법에서 보존되는지
    ② 특정 방법에서만 나타나는/사라지는 세포 유형
    ③ 경계 영역에서의 세포 유형 할당 차이
```

**해석법**:
- 모든 방법에서 동일한 공간 패턴 → robust cell typing
- 경계 영역에서 Expansion이 다른 유형을 할당 → misassignment
- Baysor가 추가로 탐지한 세포들의 공간 분포 확인

---

## 9. Summary

Step 6은 **통합 벤치마크**의 핵심 단계로, 여러 segmentation 방법을 동일 프레임워크에서 비교한다.

| 컴포넌트 | 구현 상태 | 심각도 |
|---------|---------|--------|
| Baysor execution | 완전 (crop/tile/full 모드, scale=15µm) | OK |
| Data loading | 완전 (4개 방법 자동 로딩) | OK |
| Dual benchmark mode | 완전 (Full tissue + Crop comparison) | OK |
| Concatenation | 완전 | OK |
| Preprocessing | 완전 (config: scale=true, HVG=true, n_neighbors=16; 코드 기본값: scale=false, n_neighbors=15) | OK |
| Annotation transfer | 구현 (crosstab majority vote + kNN fallback + cluster consensus) | OK (confidence 계산됨, threshold 필터링만 미적용) |
| Benchmark metrics (8/8) | **완전** (median, p5, assigned_prop, NMP_cells, NMP_reads, rand_idx) | OK |
| Clustering quality metrics | 완전 (silhouette, calinski_harabasz, davies_bouldin) | OK |
| Marker gene analysis | 완전 (per method + per celltype + dotplot) | OK |
| UMAP visualization | 완전 (3-panel) | OK |
| Spatial visualization | 완전 (method + celltype) | OK |
| Cell type barplot | 완전 | OK |
| Counts violin | 완전 | OK |
| DAPI comparison | 완전 (overlay + zoomed ROI + boundary) | OK |
| Rand Index | 완전 (`_compute_rand_index()`) | OK |
| NMP | 완전 (scRNAseq reference 조건부) | OK |
| Reads vs NMP scatter | 완전 (ExtData Fig. 5e 재현, 방법별 마커/색상) | OK |

**구현 완성도**: **HIGH** - 모든 핵심 메트릭(8/8) 호출됨. Dual benchmark 모드, sklearn clustering quality metrics, marker gene analysis 등 노트북 이상의 기능 제공.

**남은 개선 사항**:
1. Annotation confidence threshold 필터링 (노트북의 0.7 threshold) — LOW (confidence score는 이미 계산됨, 필터링만 미적용)
2. Publication quality 출력 옵션 (PDF/dpi=500) — LOW

---

## 11. 논문 Figure 직접 대응 및 시각화 정상 판별 종합 가이드

### 11.1 파이프라인 출력 → 논문 Figure 매핑 종합표

> **범례**: "관련(확장)" = 논문과 같은 분석이지만 파이프라인이 더 많은 방법/패널을 포함하여 확장한 것

| # | 파이프라인 출력 파일 | 논문 Figure | 논문 캡션 (원문) | 논문↔파이프라인 차이 | 노트북 출처 |
|---|---|---|---|---|---|
| 1 | `umap_benchmark.png` | **Fig. 3f** 관련(확장) | "UMAP from coprocessed cells using Baysor and Xenium's nuclear segmentation in a mouse brain ROI" | 논문: Baysor+Nuclei 2가지, 단일 UMAP → 파이프라인: 4가지 방법, 3-panel | `5_1_Compare_Clustering on_different_segmentations.ipynb` |
| 2 | `dapi_segmentation_comparison.png` | **Fig. 3c** / **ExtData Fig. 5b** 관련(축소) | "Comparison of cells identified with different segmentation algorithms in an ROI (160 × 160 µm), using DAPI background" | 논문: 8가지 방법 2×4 grid → 파이프라인: 4가지 방법 | `run_segmentation.py` |
| 3 | `dapi_zoomed_roi_comparison.png` | **ExtData Fig. 5a-b** 관련 | "Localization of regions of interest" + "Regions of interest representing the cells identified using different segmentation algorithms" | 논문: ROI 위치 + 8가지 비교 → 파이프라인: auto-dense ROI + 4가지 비교 | `run_segmentation.py` |
| 4 | `rand_index_heatmap.png` | **Fig. 3d** 관련(축소) | "Adjusted rand index (ARI) comparison of segmentation outputs (52 top performers) when applied to mouse brain section 2" | 논문: 52×52 heatmap → 파이프라인: 4×4 matrix | `metrics.py` → `rand_idx()` |
| 5 | `celltype_frequency_barplot.png` | **Fig. 3h** 관련(확장) | "Bar plot of cell counts per population using different segmentation strategies" | 논문: Baysor vs Nuclei 2가지, 절대수 → 파이프라인: 4가지, proportion | `5_1_Compare_Clustering on_different_segmentations.ipynb` |
| 6 | `counts_violin_by_method.png` | **Fig. 3g** 관련(확장) | "Violin plot comparing cell counts segmented by Baysor versus Xenium nuclear segmentation methods" | 논문: Baysor vs Nuclei 2가지 → 파이프라인: 4가지 | `5_1_Compare_Clustering on_different_segmentations.ipynb` |
| 7 | `reads_vs_nmp_scatter.png` | **Fig. 3e** / **ExtData Fig. 5e** | "Scatter plot of reads assigned versus negative marker purity for segmentation strategies" | 논문: 7방법×expansion×confidence 조합 → 파이프라인: 4가지 방법 | `3_4_negative_marker_purity_for_specificity.ipynb` |
| 8 | `spatial_segmentation_map.png` | 직접 대응 없음 | - | 파이프라인 자체 확장 (공간 scatter) | `5_1_Compare_Clustering on_different_segmentations.ipynb` |
| 9 | `spatial_celltype_map.png` | 직접 대응 없음 | - | 파이프라인 자체 확장 (cell type 공간 맵) | `5_1_Compare_Clustering on_different_segmentations.ipynb` |
| 10 | `marker_genes_by_segmentation.png` | 직접 대응 없음 | - | 파이프라인 자체 확장 | - |
| 11 | `marker_genes_by_celltype.png` | 직접 대응 없음 | - | 파이프라인 자체 확장 | - |
| 12 | `deg_dotplot_by_segmentation.png` | 직접 대응 없음 | - | 파이프라인 자체 확장 | - |
| 13 | `benchmark_metrics.csv` | **ExtData Fig. 5c** 관련(확장) | "Heat map representing the segmentation metrics of all segmentation strategies described in Fig. 3d" | 논문: 52 top performers × 7 메트릭 heatmap → 파이프라인: 4방법 × 11+ 메트릭 CSV | `metrics.py` |
| 14 | `rand_index_matrix.csv` | **Fig. 3d** 데이터 | (위와 동일) | 논문: 52×52 → 파이프라인: 4×4 | `metrics.py` → `rand_idx()` |
| 15 | `benchmark_combined.h5ad` | - | - | 캐시용 | - |
| 16 | `dapi_nucleus_boundaries.png` | - | - | 파이프라인 자체 확장 | - |
| 17 | `qc_counts_genes_histogram.png` | - | - | 파이프라인 자체 QC | - |
| 18 | `hvg_selection_plot.png` | - | - | 파이프라인 자체 QC | - |
| 19 | `pca_scree_plot.png` | - | - | 파이프라인 자체 QC | - |
| 20 | `crop_comparison/` | **ExtData Fig. 5a-b** 관련 | "Localization of regions of interest" | 논문: 8방법 ROI 비교 → 파이프라인: 4방법 Dual mode | `run_segmentation.py` |

### 11.2 정상 결과 판별 체크리스트

#### 시각화 1: `umap_benchmark.png` (UMAP) → Fig. 3f 관련(확장)
- **논문 Fig. 3f**: Baysor BA2 P0.8 + Nuclei **2가지만** co-processed, cell type 색상 단일 UMAP
- **파이프라인**: 4가지 방법 모두 co-processed, 3-panel (방법/클러스터/cell type)
- [ ] **클러스터 분리**: 주요 세포 유형이 잘 분리됨
- [ ] **방법 간 intermingling**: 같은 cell type의 세포가 방법에 관계없이 같은 클러스터에 혼재
- [ ] **방법 간 세포 유형 보존**: 같은 세포 유형이 일관되게 존재
- **왜 정상인가**: 논문 (p.817): "the identified cellular populations were the same across both segmentation strategies (Fig. 3f–h), with mostly mild differences in cell-type abundance"
- **읽는법**: Panel 1에서 방법별 색상이 고르게 섞여 있으면 OK. Panel 3에서 클러스터와 cell type이 일치하면 좋은 전처리.
- **비정상 신호**: 특정 방법의 세포만 별도 클러스터 형성 → batch effect 또는 systematic bias

#### 시각화 2: `dapi_segmentation_comparison.png` (DAPI ROI) → Fig. 3c / ExtData Fig. 5b 관련(축소)
- **논문 Fig. 3c**: **8가지** 방법 (Xenium cell/nuc, MESMER, Clustermap, Cellpose, Watershed, Baysor×2) 160×160µm ROI
- **파이프라인**: **4가지** 방법 (Nuclei, Cellpose, Expansion, Baysor)
- [ ] **Cellpose mask**: DAPI 핵과 정확히 일치
- [ ] **Baysor mask**: 핵 + cytoplasm 영역 포함 (불규칙 경계)
- [ ] **Nuclei**: 핵만 탐지 → 작은 타이트한 mask
- [ ] **Expansion**: 가장 넓은 coverage (과도확장 주의)
- **왜 정상인가**: 논문 Fig. 3c 캡션 — "Comparison of cells identified with different segmentation algorithms in an ROI (160 × 160 µm), using DAPI background". 각 방법의 mask 형태가 다른 것이 정상.
- **읽는법**: DAPI 배경 위 color-specific mask. Nuclei < Cellpose < Expansion 순으로 mask 크기 증가. Baysor는 불규칙 형태.
- **비정상 신호**: 모든 방법에서 mask가 없는 영역 → DAPI 품질 문제

#### 시각화 3: `rand_index_heatmap.png` (ARI Heatmap) → Fig. 3d 관련(축소)
- **논문 Fig. 3d**: **52 top performers** (315 config에서 선별) 대규모 ARI heatmap
- **파이프라인**: **4가지 방법** 간 소규모 4×4 ARI matrix
- [ ] **대각선**: 1.0 (자기 자신과의 비교)
- [ ] **DAPI-based 방법 간**: 비교적 높은 ARI (같은 staining 기반)
- [ ] **Cellpose vs Baysor**: 중간 ARI (다른 접근법)
- **왜 정상인가**: 논문 (p.817): "Staining-based strategies using DAPI generated similar outputs... Baysor-based, Clustermap-based and binning strategies clustered according to method". DAPI vs read-based는 다른 ARI가 정상.
- **읽는법**: 행/열=segmentation 방법. 색상=ARI (0-1). 높은 ARI=유사한 분할. 같은 카테고리 내 높은 ARI가 기대됨.
- **비정상 신호**: 모든 ARI < 0.3 → 방법들이 완전히 다른 결과; ARI = 1.0 (비대각선) → 두 방법이 동일 (설정 오류)

#### 시각화 4: `celltype_frequency_barplot.png` (Cell Type Barplot) → Fig. 3h 관련(확장)
- **논문 Fig. 3h**: Baysor vs Nuclei **2가지**, Y축=**절대 세포 수** (No. of cells)
- **파이프라인**: 4가지 방법, Y축=cell type proportion
- [ ] **주요 세포 유형 비율 보존**: 방법 간 비슷한 비율
- [ ] **작은 세포 유형**: 방법에 따라 약간 다를 수 있음 (정상)
- **왜 정상인가**: 논문 (p.817): "the identified cellular populations were the same across both segmentation strategies... with mostly mild differences in cell-type abundance"
- **읽는법**: X축=방법, Y축=proportion. 색상=cell type. 비율이 방법 간 유사하면 일관된 결과.

#### 시각화 5: `counts_violin_by_method.png` (Counts Violin) → Fig. 3g 관련(확장)
- **논문 Fig. 3g**: Baysor BA2 P0.8 vs Nuclei **2가지만** violin
- **파이프라인**: 4가지 방법 모두 포함
- [ ] **Baysor**: 가장 높은 reads/cell
- [ ] **Nuclei**: 가장 낮은 reads/cell (핵 내만)
- [ ] **분포 폭**: 적당 (너무 넓지 않음)
- **왜 정상인가**: 논문 (p.817): "cells defined by Baysor had a higher count per cell". Baysor가 cytoplasmic reads를 추가 할당하므로 nuclei보다 높은 counts가 정상.
- **읽는법**: X축=방법, Y축=counts/cell. Violin 폭=분포. 높은 median + 좁은 분포 = 일관된 품질.

#### 시각화 6: `reads_vs_nmp_scatter.png` (NMP vs Assigned Reads) → Fig. 3e / ExtData Fig. 5e
- **논문 Fig. 3e**: 7방법 × expansion(0-14.9µm) × prior confidence(0-0.99) 조합. 방법별 마커: Baysor=●, Xenium=✱, Cellpose=◆, MESMER=✦, Clustermap=◆, Watershed=▲, Binning=✱. Expansion 크기별 색상, prior confidence별 크기
- **파이프라인**: 4가지 방법, 각각 단일 점
- [ ] **NMP > 0.8**: 논문에서 "high specificity (NCP > 0.8)" 기준
- [ ] **Baysor BA2 P0.8**: 높은 NMP + 높은 Assigned% (우상단) → 최적
- [ ] **Nuclei**: 높은 NMP + 낮은 Assigned% (보수적)
- [ ] **과도 Expansion**: 높은 Assigned% + 낮은 NMP (우하단) → 경고
- **왜 정상인가**: 논문 (p.817): "We defined the optimal segmentation strategy as the one maximizing the proportion of reads assigned to cells while maintaining specific expression patterns, quantified by negative marker purity (NMP)"
- **읽는법**: X축=Assigned reads 비율, Y축=NMP. 우상단=최적. 좌상단=보수적(reads 손실). 우하단=과도확장(reads 오염).

### 11.3 논문 원문 캡션 인용 (정확한 인용)

> **참고**: 아래는 논문 PDF p.818 Figure 3 caption과 p.25 Extended Data Fig. 5 caption에서 **직접 발췌**한 것이다.

| Figure | 논문 원문 캡션 (직접 발췌) | 위치 |
|---|---|---|
| **Fig. 3a** | "Mouse brain region with reads overlaid on DAPI staining, colored by distance to the nearest cell centroid (left). The line plot (right) shows the PCC of oligodendrocytes to the nuclear signature (blue) and to the background signature (orange), depending on the distance to the cell centroid." | p.818 caption |
| **Fig. 3b** | "Bar plot showing the distance in micrometers of the intersection between nuclei and the domain-specific regions, as in a, across cell types. Error bars represent the 95% confidence intervals. Mean nuclei and cell radius are also shown." | p.818 caption |
| **Fig. 3c** | "Comparison of cells identified with different segmentation algorithms in an ROI (160 × 160 µm), using DAPI background." | p.818 caption |
| **Fig. 3d** | "Adjusted rand index (ARI) comparison of segmentation outputs (52 top performers) when applied to one of the mouse brain samples profiled (mouse brain section 2)." | p.818 caption |
| **Fig. 3e** | "Scatter plot of reads assigned versus negative marker purity for segmentation strategies applied to mouse brain section 2." | p.818 caption |
| **Fig. 3f** | "UMAP from coprocessed cells using Baysor and Xenium's nuclear segmentation in a mouse brain ROI." | p.818 caption |
| **Fig. 3g** | "Violin plot comparing cell counts segmented by Baysor versus Xenium nuclear segmentation methods." | p.818 caption |
| **Fig. 3h** | "Bar plot of cell counts per population using different segmentation strategies." | p.818 caption |
| **ExtData Fig. 5a** | "Localization of regions of interest represented in Extended Data Fig. 5b and Fig. 3c." | p.25 caption |
| **ExtData Fig. 5b** | "Regions of interest representing the cells identified using different segmentation algorithms in a region of interest outlined in Extended Data Fig 5B. DAPI background is represented as a background and individual isolated color-specific masks represent individual cells." | p.25 caption |
| **ExtData Fig. 5c** | "Heat map representing the segmentation metrics of all segmentation strategies described in Fig. 3d." | p.25 caption |
| **ExtData Fig. 5d** | "Adjusted rand index (ARI) between the different outputs produced by combinations of segmentation algorithms, hyperparameters and expansions when applied to human breast sections." | p.25 caption |
| **ExtData Fig. 5e** | "Scatter plot representing the number of reads assigned (x-axis) and the negative marker purity (y-axis) of different assessed segmentation strategies in human breast tumor samples." | p.25 caption |

**논문 본문 핵심 인용 (segmentation 관련)**:

| 주제 | 논문 원문 | 위치 |
|---|---|---|
| Segmentation 비교 | "We benchmarked the segmentation methods Baysor, MESMER, Watershed, Cellpose and Clustermap against the segmentations provided by 10x Genomics" | p.817 |
| NMP 정의 | "We defined the optimal segmentation strategy as the one maximizing the proportion of reads assigned to cells while maintaining specific expression patterns, quantified by negative marker purity (NMP)" | p.817 |
| 최적 전략 | "Baysor combined with Xenium's nuclei segmentation (BA2 P0.8), represent the best segmentation strategy (Fig. 3e and Extended Data Fig. 5e)" | p.817 |
| Co-processed 결과 | "the identified cellular populations were the same across both segmentation strategies (Fig. 3f–h), with mostly mild differences in cell-type abundance" | p.817 |
| Baysor counts | "cells defined by Baysor had a higher count per cell" | p.817 |
| ARI 클러스터링 | "Staining-based strategies using DAPI generated similar outputs, with cell expansion being the force driving their differences" | p.817 |
| 최적 파이프라인 결론 | "the optimal algorithm involves two steps: first, identifying nuclei using Cellpose and second, assigning reads to individual cells using Baysor" | p.821 |
| NCP 기준 | "all the different SRT technologies presented a mean high specificity (NCP > 0.8)" | p.816-817 |
