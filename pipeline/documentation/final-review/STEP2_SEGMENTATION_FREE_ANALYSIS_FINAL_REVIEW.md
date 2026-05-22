# Step 2: Segmentation-Free Analysis - Final Review

## 1. Overview

### 1.1 What This Step Does
Step 2는 **세포 경계를 정의하지 않고** 공간 transcript 패턴을 분석하는 segmentation-free 접근법을 구현한다. Points2Regions를 사용한 spatial domain clustering, 핵 경계 기반 거리 분석(EDT), centroid 기반 거리 분석, ovrlpy를 통한 3D coherence 분석, SSAM 기반 de novo cell typing을 수행한다.

### 1.2 Paper Context
논문의 "Xenium retains key 3D and subcellular cell information" 섹션 (p.814)에서 segmentation-free 분석의 의의를 설명한다:
- Points2Regions를 사용하여 44개 cell-type-specific 클러스터를 식별 (Extended Data Fig. 3a)
- 핵 경계 기반으로 nuclear/cytoplasmic 클러스터 분류
- SSAM de novo 모드로 세포 유형 맵 생성 (Extended Data Fig. 2)
- 3D z-axis coherence 분석으로 overlapping cells 탐지

### 1.3 논문 Figure 매칭
| 논문 Figure | 설명 | Step 2 구현 |
|------------|------|------------|
| **Fig. 1e** | 3D coherence map (ovrlpy) | `run_ovrlpy_analysis()` |
| **Fig. 1f** | Distance to centroid boxplot (nuclear vs cyto) | `plot_centroid_distance_analysis()` |
| **Fig. 1i** | Points2Regions clusters (Oligo) | `classify_p2r_subcellular()` + `plot_p2r_heatmap()` |
| **Fig. 1j** | Distance to nuclear edge boxplot | `plot_p2r_boundary_distance()` |
| **Fig. 1k** | Top genes per subcellular cluster | `plot_p2r_top_genes()` |
| **Extended Data Fig. 2a** | SSAM cell type map (whole brain) | `run_ssam_analysis()` |
| **Extended Data Fig. 2b** | SSAM coherence map | `run_ssam_analysis()` |
| **Extended Data Fig. 2c** | SSAM UMAP | `_ssam_plots()` |
| **Extended Data Fig. 3a** | P2R spatial map (whole tissue) | `plot_p2r_heatmap()` |
| **Extended Data Fig. 3b** | P2R cluster × cell type confusion matrix | `plot_p2r_heatmap()` |
| **Extended Data Fig. 3c** | Distance to nuclear edge per P2R cluster | `plot_p2r_boundary_distance()` |
| **Extended Data Fig. 3d** | Top genes per subcellular cluster (astrocytes) | `plot_p2r_top_genes()` |

---

## 2. Pipeline Process Flow

```
Step 1 Output
├── {sample_tag}_step1_exploration.h5ad
├── {sample_tag}_transcripts.parquet
         ↓
    [Step 2: Segmentation-Free Analysis]
         ↓
    ┌───── 2-0. Celltype Pre-population ───────────────────────┐
    │  spots['cell_id'] → adata.obs[leiden/Class] 매핑          │
    │  _normalize_cell_ids()로 float→int 문자열 정규화          │
    │  → spots['celltype'] 컬럼 생성 (downstream P2R 분석용)    │
    └──────────────────────────────────────────────────────────┘
         ↓
    ┌───── 2-1. Points2Regions ──────────────────────────────┐
    │  Transcript → spatial bins → Gaussian blur → KMeans    │
    │  → domain clusters (nuclear/cytoplasmic 분류)           │
    │  → P2R 색상 할당 (HLS palette, nuclear proximity 기반) │
    └────────────────────────────────────────────────────────┘
         ↓
    ┌───── 2-2. ovrlpy (선택) ──────────────────────────────┐
    │  z-axis coherence (VSI) → overlapping cell 탐지       │
    │  → per-cell incoherence 점수 계산                      │
    └────────────────────────────────────────────────────────┘
         ↓
    ┌───── 2-3. Distance Metrics ────────────────────────────┐
    │  Centroid distance: transcript → cell centroid          │
    │  Boundary distance: signed EDT (neg=inside nucleus)    │
    └────────────────────────────────────────────────────────┘
         ↓
    ┌───── 2-4. Visualizations ─────────────────────────────┐
    │  Centroid analysis, P2R heatmap/top genes,             │
    │  P2R boundary distance, distance histograms,           │
    │  boundary spatial overlay (전체+확대 ROI)               │
    └────────────────────────────────────────────────────────┘
         ↓
    ┌───── 2-5. SSAM (선택) ─────────────────────────────────┐
    │  KDE-based gene expression field → cell type map       │
    └────────────────────────────────────────────────────────┘
         ↓
├── {sample_tag}_step2_points2regions.h5ad
├── points2regions/ 분석 결과 디렉토리
```

---

## 3. Sub-step 상세 분석

### 3.1 Points2Regions Clustering

**개념**: Points2Regions (Andersson et al., Cytometry A, 2024)는 공간 transcript 데이터를 bin 단위로 집계하여 molecular pattern을 클러스터링하는 segmentation-free 도구이다.

**처리 과정**:
```
1. Transcript 좌표를 공간 bin으로 집계
   - bin_width: 1.0 µm (config 설정 가능)

2. 각 bin의 유전자 조성 벡터 생성
   - bin별 gene × count 매트릭스

3. Gaussian filter로 인접 bin 블러링
   - sigma: 3.0 µm (σ가 클수록 더 넓은 도메인)

4. 저발현 bin 필터링
   - min_genes_per_bin: 15 (최소 15개 유전자 발현 필요)

5. Mini-batch KMeans 클러스터링
   - n_clusters: [50, 100, 200, 500] (다중 해상도)

6. Nuclear/Cytoplasmic 분류
   - overlaps_nucleus > 0.5 → nuclear cluster
   - overlaps_nucleus ≤ 0.5 → cytoplasmic cluster
```

**설정값**:
```yaml
segmentation_free:
  points2regions:
    sigma: 3.0              # Gaussian blur σ (µm)
    n_clusters: [50, 100, 200, 500]  # 클러스터 수 목록
    bin_width: 1.0           # 공간 bin 너비 (µm)
    min_genes_per_bin: 15    # 최소 유전자/bin
```

**출력**:
- `{sample_tag}_step2_p2r_classification.csv` - P2R 클러스터별 nuclear/cyto 분류
- `{sample_tag}_step2_points2regions_k{n}_bins.h5ad` - n개 클러스터의 bin-level AnnData

### 3.2 P2R Subcellular Classification (`classify_p2r_subcellular()`, lines 82-154)

**핵심 로직**:
```python
# 각 P2R 클러스터 내 transcript의 overlaps_nucleus 비율 계산
for cluster_id in unique_clusters:
    cluster_transcripts = transcripts[transcripts['p2r_cluster'] == cluster_id]
    nucleus_fraction = cluster_transcripts['overlaps_nucleus'].mean()

    if nucleus_fraction > 0.5:
        classification[cluster_id] = 'nuclear'
    else:
        classification[cluster_id] = 'cytoplasmic'
```

**의미**:
- Nuclear cluster: 대다수 transcript가 DAPI 핵 mask 내부에 위치 → nuclear-enriched genes
- Cytoplasmic cluster: 대다수 transcript가 핵 외부 → cytoplasmic/membrane genes
- 논문 Extended Data Fig. 3b에서 이 분류를 confusion matrix로 시각화

### 3.2b Cell ID Normalization (`_normalize_cell_ids()`, lines 157-161)

Xenium 데이터에서 cell_id의 dtype이 일관되지 않는 문제를 해결하는 헬퍼 함수.

**문제**: spots의 `cell_id`가 float (1.0, 2.0 등)인 경우, `.astype(str)`이 `"1.0"`, `"2.0"`을 생성하여 adata.obs의 `"1"`, `"2"`와 매핑이 실패한다.

```python
def _normalize_cell_ids(series):
    """'1.0' → '1', '2' → '2' 등 float-string을 int-string으로 정규화"""
    s = series.astype(str)
    return s.str.replace(r'\.0$', '', regex=True)
```

### 3.2c Cell Type Resolution (`_resolve_celltype_column()`, lines 164-226)

Spots에서 cell type 컬럼을 찾거나 생성하는 3단계 fallback 로직:

1. **직접 컬럼 검색**: spots에 이미 `Class`, `celltype`, `cell_type`, `leiden_*` 컬럼이 있으면 반환
2. **cell_id 매핑**: `_normalize_cell_ids()`로 정규화한 후 spots['cell_id'] → adata.obs 매핑
3. **공간 최근접 매핑**: cell_id 매핑 실패 시, cKDTree로 50µm 이내 최근접 세포의 cell type 할당

**노트북 대응**: P2R `compute_distance.ipynb`에서 직접 매핑하던 패턴:
```python
cell_id_to_class = adata_processed.obs.set_index("cell_id")["Class"].to_dict()
adata.uns["spots"]["Class"] = adata.uns["spots"]["cell_id"].map(cell_id_to_class)
```
을 `_normalize_cell_ids()`를 통해 dtype 안전하게 구현.

### 3.2d Celltype Pre-population in `run_step2()` (lines 1574-1594)

P2R 분석 시작 전, spots에 celltype 컬럼을 미리 매핑하는 단계:

```python
spots_cid = _normalize_cell_ids(spots['cell_id'])
obs_cid = _normalize_cell_ids(adata.obs.index.to_series())
cid_to_ct = dict(zip(obs_cid, adata.obs[ct_col]))
spots['celltype'] = spots_cid.map(cid_to_ct).fillna('Background')
```

**이유**: 이전에는 `_resolve_celltype_column()`이 downstream P2R 함수에서만 호출되어, dtype 불일치로 매핑이 실패할 수 있었다. 이제 `run_step2()` 초반에 미리 매핑하여 일관된 celltype 라벨을 보장한다.

### 3.3 Distance to Centroid (`calculate_distance_to_centroid()`, lines 697-754)

**계산식**:
```
d_centroid(transcript_i) = sqrt((x_i - x_cell)² + (y_i - y_cell)²)
```
여기서 `(x_cell, y_cell)`은 해당 transcript가 할당된 세포의 centroid 좌표.

**의미**: Centroid distance는 transcript가 세포 중심에서 얼마나 떨어져 있는지를 측정한다. 핵 유전자(e.g., Opalin)는 짧은 거리, 세포질 유전자(e.g., MAG)는 긴 거리를 보인다.

**출력**: `spots['dist_to_centroid']` 컬럼 + `adata.uns['gene_distance_stats']` (유전자별 평균/중앙값/표준편차)

### 3.4 Distance to Boundary - Signed EDT (`calculate_distance_to_boundary()`, lines 763-817)

**개념**: Euclidean Distance Transform (EDT)을 사용하여 각 transcript의 가장 가까운 핵 경계까지 거리를 계산한다.

**핵심 특징 - Signed Convention**:
```
d_boundary < 0 → 핵 내부 (inside nucleus)
d_boundary = 0 → 핵 경계 위 (on boundary)
d_boundary > 0 → 핵 외부 (outside nucleus)
```

**처리 과정**:
1. Nuclear segmentation mask → binary image
2. EDT 계산: 각 pixel → nearest boundary pixel 거리
3. 핵 내부 pixel에 음수 부호 부여 (signed convention)
4. Transcript 좌표를 pixel 단위로 변환
5. EDT map에서 각 transcript 위치의 거리값 조회

**논문과의 관계**: Fig. 1j와 Extended Data Fig. 3c에서 P2R 클러스터별 nuclear edge까지 거리를 boxplot으로 시각화.

### 3.5 Centroid Distance Visualization (`plot_centroid_distance_analysis()`, lines 812-964)

**시각화 목록**:

| Plot | 설명 | 논문 Figure |
|------|------|------------|
| `centroid_stripplot.png` | 극단값 유전자의 centroid distance (상위/하위) | Fig. 1f |
| `centroid_boxplot.png` | 거리 카테고리별 boxplot (Core/Intermediate/Peripheral/Outer) | - |
| `centroid_ecdf.png` | 전체 + 상위 유전자별 ECDF | Fig. 2f |
| `centroid_welchs_heatmap.png` | Domain 간 거리 분포 차이 (Welch's t-test) | - |

**거리 카테고리 정의**:
- Core (0-25%): 핵 중심 근처
- Intermediate (25-50%): 핵 내부-경계 사이
- Peripheral (50-75%): 경계 근처
- Outer (75-100%): 세포 외곽

### 3.6 P2R Boundary Distance (`plot_p2r_boundary_distance()`, lines 972-1038)

각 P2R 클러스터의 transcript들이 핵 경계로부터 어떤 거리에 분포하는지 boxplot으로 시각화.
- X축: 거리 to nuclear edge (µm)
- Y축: P2R 클러스터 (nuclear/cyto 분류 라벨)
- 빨간 점선: 핵 경계 (distance = 0)

### 3.7 P2R Top Genes (`plot_p2r_top_genes()`, lines 340-440)

각 세포 유형의 subcellular cluster (nuclear/cytoplasmic)별 상위 발현 유전자를 stacked bar chart로 시각화.
- Y축: 유전자별 상대 발현 비율
- 패널: 각 subcellular cluster (e.g., Astrocytes_nuclear, Astrocytes_cyto_38)
- 논문 Extended Data Fig. 3d에 해당

### 3.7b P2R Mean Distance Scatter (`plot_p2r_mean_distance_scatter()`, lines 1059-1123)

각 P2R 클러스터의 평균 signed distance를 scatter plot으로 시각화. X축=클러스터 인덱스, Y축=평균 signed distance, 색상=HLS palette. 극단값 클러스터에 라벨 표시.

### 3.7c Distance Histograms (`_plot_distance_histograms()`, lines 1742-1784)

Boundary distance (signed + absolute)와 centroid distance의 요약 히스토그램. Log scale 사용, boundary plot에서 x=0 (핵 경계)을 빨간 점선으로 표시.

**출력**:
- `{sample_tag}_step2_distance_boundary_dist.png` - 2패널: signed (좌) + absolute (우)
- `{sample_tag}_step2_distance_centroid_dist.png` - centroid distance 분포

### 3.7d Boundary Spatial Overlay (`plot_boundary_spatial_overlay()`, lines 1360-1540)

핵 경계 polygon을 transcript 공간 맵 위에 오버레이하는 시각화.

**출력**:
- `{sample_tag}_step2_boundary_overview.png` - 2패널 전체 조직 뷰 (100K transcript + 1K boundaries 서브샘플)
- `{sample_tag}_step2_boundary_zoomed_roi.png` - 3패널 500µm ROI 확대 (transcripts, boundaries, centroids)

### 3.8 ovrlpy Analysis (`run_ovrlpy_analysis()`, lines 462-597)

**개념**: ovrlpy는 tissue의 z-axis coherence를 분석하여 3D 공간에서 신호가 겹치는 (overlapping) 세포를 탐지한다.

**처리 과정**:
1. Transcript (x, y, z) 좌표로 KDE 생성
2. Tissue를 top/bottom 반으로 분할 (z 중앙값 기준)
3. 각 위치에서 top/bottom 유전자 발현 signature 비교 (cosine similarity)
4. Coherence map 생성: cosine similarity < threshold → overlapping region

**설정값**:
```yaml
segmentation_free:
  overlaps:
    um_per_pixel: 2.0    # KDE 해상도
    radius: 50.0          # 분석 반경 (µm)
    bw: 1                 # KDE bandwidth
```

**시각화**:
- `coherence_map.png` - 전체 조직의 z-coherence map (논문 Fig. 1e)
- `ovrlpy_pseudocells.png` - Pseudo-cell 기반 coherence
- `ovrlpy_tissue.png` - 조직 레벨 coherence
- `ovrlpy_doublets.csv` - 검출된 overlapping cell 목록

**분석법**:
- Coherence 0에 가까움 = top/bottom 신호 불일치 → overlapping cells 가능성
- Coherence 0.8 이상 = top/bottom 일관됨 → 정상 세포

### 3.8b Per-Cell Incoherence (`_compute_and_save_cell_incoherence()`, lines 600-689)

ovrlpy coherence map에서 각 세포의 incoherence 점수를 계산. 핵 경계 polygon 내부 영역의 integrity map 값을 추출하여 세포별 평균/최대/중앙값 incoherence를 산출한다.

**출력**:
- `{sample_tag}_step2_cell_incoherence.csv` - 세포별 incoherence 메트릭
- `{sample_tag}_step2_incoherence_hist.png` - Incoherence 분포 히스토그램 (low-coherence threshold 표시)

---

## 4. SSAM Analysis 상세 검토

### 4.1 SSAM 이론적 배경

**SSAM** (Spot-based Spatial cell-type Analysis by Multidimensional mRNA density estimation, Park et al., Nature Communications 12, 3545, 2021)은 세포 segmentation 없이 공간 전사체 데이터에서 세포 유형 맵을 생성하는 알고리즘이다.

**SSAM 논문 4단계 알고리즘**:

```
Step 1: Spatial mRNA Density 계산 (KDE)
  ┌─────────────────────────────────────────────────────────────┐
  │ 각 유전자 i에 대해 KDE 밀도 추정:                            │
  │                                                              │
  │   σ̂_h(x) = (1/Nh) × Σ K((x - x_j) / h)                   │
  │                                                              │
  │ 여기서:                                                      │
  │   K = Gaussian kernel: K(x) = (1/√2π) × exp(-0.5|x|²)     │
  │   h = bandwidth (세포 직경에 비례, ~10 µm → FWHM ≈ 2h)     │
  │   N = 해당 유전자의 총 mRNA 수                               │
  │                                                              │
  │ Gene expression at position x:                               │
  │   E_i(x) = σ̂_i(x) × N_i                                   │
  │                                                              │
  │ 결과: n_genes 차원의 벡터 필드 (각 pixel = gene expression   │
  │       vector). 예: 33 유전자 × 3380×2080 pixel grid          │
  │       = 7,030,400 개의 33차원 벡터                           │
  └─────────────────────────────────────────────────────────────┘
                    ↓
Step 2: Cell-type Signature 식별
  ┌─────────────────────────────────────────────────────────────┐
  │ 2a. L1 norm 계산 → local maxima 탐지 (= pseudo-cell 위치) │
  │ 2b. Local maxima의 벡터를 embedding (tSNE/UMAP)            │
  │ 2c. 클러스터링 → cell-type cluster 식별                     │
  │ 2d. 클러스터 centroids = representative signature vectors   │
  │                                                              │
  │ 또는 Guided Mode:                                           │
  │ 2e. 외부 reference (scRNA-seq)에서 signature 직접 사용      │
  └─────────────────────────────────────────────────────────────┘
                    ↓
Step 3: Cell-type Map 생성
  ┌─────────────────────────────────────────────────────────────┐
  │ 각 pixel의 벡터와 centroid 간 Pearson correlation 계산      │
  │ → 최고 상관 centroid의 cell type으로 할당                    │
  │ → min_r threshold로 low-confidence pixel 필터링             │
  │ → per-centroid cell-type map → composite map 병합           │
  └─────────────────────────────────────────────────────────────┘
                    ↓
Step 4: Domain Map 구축 (선택적)
  ┌─────────────────────────────────────────────────────────────┐
  │ Cell-type map을 공간 window로 분할                           │
  │ → window별 cell-type composition 벡터 계산                  │
  │ → Agglomerative clustering → tissue domain 식별             │
  └─────────────────────────────────────────────────────────────┘
```

### 4.2 SSAM Pipeline 구현 (`run_ssam_analysis()`, lines 1360-1503)

**Phase A: KDE Vector Field + Local Maxima Detection**

```python
# 1. 초기화
ds = ssam.SSAMDataset()                          # v1.1+ API (빈 dataset)
analysis = ssam.SSAMAnalysis(ds, ncores=4, verbose=True)

# 2. Transcript 좌표 준비
spots = adata.uns['spots']
gene_mask = ~spots['feature_name'].str.contains('BLANK|NegControl')  # 컨트롤 제거
coords = spots[gene_mask][['x_location', 'y_location', 'feature_name']]
coords[['x', 'y']] -= coords[['x', 'y']].min()  # 원점 이동
coords['x'] /= um_per_px                         # µm → pixel 변환
coords['y'] /= um_per_px

# 3. KDE 실행 (bandwidth를 pixel 단위로 변환)
kde_bw_px = kde_bandwidth_um / um_per_px          # 2.5 µm / 2.0 = 1.25 px
analysis.run_kde(
    locations=locations_df,                        # DataFrame with gene, x, y
    width=width_px, height=height_px, depth=1,
    bandwidth=kde_bw_px,
    sampling_distance=1.0,
    re_run=True,
)

# 4. Thresholds 설정 + local maxima 탐지
analysis.set_thresholds(
    expression_threshold=0.2,                      # 최소 발현 수준
    norm_threshold=5,                              # 최소 L1 norm
)
search_size = max(3, int(round(min_dist / um_per_px)))  # 7µm / 2.0 = 3.5 → 3 or 5 (odd)
analysis.find_localmax(search_size=search_size)

# 5. Normalize + Scale vectors (SSAM 튜토리얼 순서 준수)
analysis.normalize_vectors(normalize_vector=True)
analysis.scale_vectors()                           # ← 2026-03-09 추가 (누락 수정)
```

**Phase B: Clustering + Cell Type Mapping**

```python
# 6. Local maxima에서 AnnData 구축
normalized = ds.normalized_vectors                 # shape: (n_localmax, n_genes)
local_maxs = ds.local_maxs                         # (x_coords, y_coords)
ssam_adata = sc.AnnData(
    np.array(normalized),
    var=pd.DataFrame(index=ssam_genes),
    obs=pd.DataFrame({'x': local_maxs[0], 'y': local_maxs[1]})
)

# 7. Scanpy 클러스터링 파이프라인
sc.pp.normalize_total(ssam_adata, target_sum=1e4)
sc.tl.pca(ssam_adata, svd_solver='arpack')
sc.pp.neighbors(ssam_adata, n_neighbors=15, n_pcs=min(40, n_vars-1, n_obs-1))
sc.tl.umap(ssam_adata, min_dist=0.02, random_state=42)
sc.tl.leiden(ssam_adata, resolution=2.0, random_state=42)

# 8. scRNA-seq 참조 기반 cell type 매핑
mapping_result = _ssam_map_celltypes(analysis, ssam_adata, config, ssam_genes)

# 9. 저장 + 시각화
ssam_adata.write_h5ad(f"{sample_tag}_step2_ssam.h5ad")
_ssam_plots(ssam_adata, ds, mapping_result, output_dir, sample_tag)
```

### 4.3 SSAM Cell Type Mapping (`_ssam_map_celltypes()`, lines 1506-1625)

**Step 1: scRNA-seq Reference 로드**
- 소스: SEA-AD MTG scRNA-seq (auto-download from S3)
- 경로: `config.sc_reference.dest_dir` → `data/scRNAseq/*.h5ad`
- Cell type key: `subclass_label` (configurable)

**Step 2: Reference Signature 구축**
```python
# 각 cell type별 평균 발현 벡터 (common genes만)
for ct in celltypes:
    ct_data = ref_sub[ref_sub.obs[ct_key] == ct]
    ref_signatures[ct] = ct_data.X.mean(axis=0)
ref_signatures = np.log1p(ref_signatures)          # log1p 정규화
```

**Step 3: Leiden Cluster Signature 구축**
```python
# 각 Leiden cluster별 평균 발현 벡터
for cl in leiden_cats:
    mask = ssam_adata.obs['leiden'] == cl
    leiden_signatures[cl] = ssam_sub[mask].X.mean(axis=0)
leiden_signatures = np.log1p(leiden_signatures)     # log1p 정규화
```

**Step 4: Pearson Correlation Matrix**
```python
# (n_celltypes × n_clusters) correlation matrix
for ct in celltypes:
    for cl in leiden_cats:
        r, _ = pearsonr(ref_sig[ct].values, lei_sig[cl].values)
        corr_matrix.loc[ct, cl] = r
```

**Step 5: Cell Type Assignment**
```python
# 각 cluster → 최고 상관 reference cell type
for cl in leiden_cats:
    cluster_assignments[cl] = corr_matrix[cl].astype(float).idxmax()
ssam_adata.obs['leiden_assignment'] = ssam_adata.obs['leiden'].map(cluster_assignments)
```

**Step 6: Pixel-Level Mapping**
```python
# SSAM 내장 함수로 pixel-level cell type map 생성
analysis.map_celltypes(centroids=leiden_signatures.values.T)  # shape: (n_clusters, n_genes)
analysis.filter_celltypemaps(min_norm=3, min_r=0.2)
```

### 4.4 SSAM 시각화 (`_ssam_plots()`, lines 1639-1800)

| # | 출력 파일 | 설명 | 논문 Figure | 조건 |
|---|----------|------|-----------|------|
| 1 | `{tag}_step2_ssam_umap.png` | Leiden 클러스터별 Scanpy UMAP | Ext. Fig. 2c | 항상 |
| 2 | `{tag}_step2_ssam_spatial.png` | Spatial scatter (leiden_assignment 색상) | Ext. Fig. 2a | 항상 |
| 3 | `{tag}_step2_ssam_highest_genes.png` | 상위 20개 발현 유전자 bar plot | - | 항상 |
| 4 | `{tag}_step2_ssam_l1norm_localmax.png` | L1 Norm KDE + Local Maxima 오버레이 (QC 진단) | 튜토리얼 Cell 12 | 항상 |
| 5 | `{tag}_step2_ssam_native_umap.png` | SSAM 내장 UMAP (`ds.plot_umap`) | 튜토리얼 Cell 16 | 항상 |
| 6 | `{tag}_step2_ssam_pca.png` | PCA variance plot (embedding QC) | - | 항상 |
| 7 | `{tag}_step2_ssam_diagnostic_{i}.png` | Per-cluster 3-panel diagnostic (spatial + genes + UMAP) | 튜토리얼 Cell 20 | 항상 (최대 15개) |
| 8 | `{tag}_step2_ssam_celltype_map.png` | Pixel-level SSAM cell type map (`ds.plot_celltypes_map`) | Ext. Fig. 2a | Reference 있을 때 |
| 9 | `{tag}_step2_ssam_umap_celltypes.png` | leiden_assignment UMAP | - | Reference 있을 때 |
| 10 | `{tag}_step2_ssam_correlation.png` | Pearson correlation heatmap (clusters × celltypes) | - | Reference 있을 때 |
| 11 | `{tag}_step2_ssam_kde_3d.png` | 3D wireframe surface (KDE density field) | PDF slide 61 | 항상 |
| 12 | `{tag}_step2_ssam_downsampling.png` | 2-panel: L1 norm heatmap vs L1 maxima scatter | PDF slide 63 | 항상 |
| 13 | `{tag}_step2_ssam_gene_kde.png` | Top 4 유전자별 개별 KDE heatmap | PDF slide 62 | 항상 |
| 14 | `{tag}_step2_ssam_celltype_map_zoomed.png` | Full + Zoomed ROI cell type map | 논문 Fig zoom | Reference 있을 때 |

---

## 5. SSAM 구현 비교: 논문/튜토리얼 vs 노트북 vs 파이프라인

### 5.1 3자 비교표

| 기능/단계 | SSAM 논문 + 박정빈 튜토리얼 | 노트북 (`2_4_brain_ssam.ipynb`) | 파이프라인 (`run_ssam_analysis()`) | 상태 |
|-----------|---------------------------|-------------------------------|-----------------------------------|------|
| **API 버전** | v1.0 (`SSAMDataset(genes, loci, h, w)`) | v1.0 (`SSAMDataset(genes, loci, h, w)`) | v1.1+ (`SSAMDataset()` + `run_kde()`) | 파이프라인이 최신 |
| **KDE 실행** | `run_fast_kde()` 또는 `run_kde()` | `analysis.run_fast_kde(bandwidth=2.5)` | `analysis.run_kde(locations=df, bandwidth=kde_bw_px)` | 동등 (API만 다름) |
| **Local maxima 탐지** | `analysis.find_localmax()` | `analysis.find_localmax()` | `analysis.find_localmax(search_size=N)` | 동등 |
| **벡터 정규화** | `normalize_vectors()` | `normalize_vectors(normalize_vector=True)` | `normalize_vectors(normalize_vector=True)` | 동등 |
| **벡터 스케일링** | `scale_vectors()` (별도 호출) | `scale_vectors()` 호출 없음 | `scale_vectors()` (2026-03-09 추가) | **수정 완료** |
| **De novo 클러스터링** | `analysis.cluster_vectors()` (HDBSCAN) | Scanpy Leiden (resolution=2.0) | `analysis.cluster_vectors()` + Scanpy Leiden (병행) | **수정 완료** (양쪽 모두 실행) |
| **Reference 타입** | smFISH cell-by-gene (동일 플랫폼) | DAPI 사전 annotation + Allen scRNA-seq | SEA-AD scRNA-seq | 모두 유효한 Guided Mode |
| **Signature 구축** | `np.log(signatures + 1)` | `np.log(signatures + 1)` | `np.log1p()` | 동등 (`log1p = log(x+1)`) |
| **Cell type 매핑** | `analysis.map_celltypes(centroids)` | `analysis.map_celltypes(sigs_values)` | `analysis.map_celltypes(centroids=lei_sig.T)` | 동등 |
| **Map 필터링** | `filter_celltypemaps(min_r=0.3)` | `filter_celltypemaps(min_norm=3, min_r=0.2)` | `filter_celltypemaps(min_norm=3, min_r=0.2)` | 동등 (threshold 차이) |
| **Pixel-level map 시각화** | `ds.plot_celltypes_map(colors=...)` | `ds.plot_celltypes_map(colors=..., rotate=3)` | `ds.plot_celltypes_map(colors=..., ax=ax)` | 동등 |
| **UMAP** | `ds.run_umap()` + `ds.plot_umap()` | Scanpy UMAP | `ds.run_umap()` + Scanpy UMAP (병행) | **수정 완료** (양쪽 모두 실행) |
| **Correlation heatmap** | 없음 | 없음 | `sns.heatmap(corr_matrix)` | 파이프라인 추가 |
| **L1 norm 시각화** | `ds.plot_l1norm()` | 없음 | `ds.plot_l1norm()` + `ds.plot_localmax()` 오버레이 | **구현 완료** |
| **Local maxima 시각화** | `ds.plot_localmax(s=0.1, c='red')` | `ds.plot_localmax()` | `ds.plot_localmax(s=0.3, c='red')` (L1 norm 위에 오버레이) | **구현 완료** |
| **SSAM Native UMAP** | `ds.plot_umap(s=1)` | 없음 | `ds.plot_umap(s=1)` | **구현 완료** |
| **Diagnostic plots** | `ds.plot_diagnostic_plot(i, use_embedding='umap')` × N | 없음 | `ds.plot_diagnostic_plot(i, use_embedding='umap')` × 15 | **구현 완료** |
| **Watershed segmentation** | `analysis.run_watershed(vfnorm_thresh_im)` | 없음 | 없음 | 미구현 (선택적) |
| **Domain map** | Agglomerative clustering on windows | 없음 | 없음 | 미구현 (선택적) |
| **zarr v3 호환** | N/A | N/A | `_patch_zarr_v3_compat()` | 파이프라인 전용 |
| **BLANK/NegControl 필터** | N/A | 없음 | `gene_mask` 필터 | 파이프라인 개선 |
| **PCA plot** | 없음 | `sc.pl.pca(adata)` | `sc.pl.pca(ssam_adata)` | **구현 완료** |
| **log1p 전처리** | `sc.pp.log1p(adata)` (reference에 적용) | 없음 | `sc.pp.log1p(ssam_adata)` (scanpy pipeline에 추가) | **수정 완료** |
| **3D KDE surface** | PDF slide 61 (wireframe) | 없음 | `plot_wireframe(ds.vf_norm)` subsampled | **구현 완료** |
| **다운샘플링 비교** | PDF slide 63 (L1 norm vs maxima) | 없음 | 2-panel: `imshow(vf_norm)` + `scatter(local_maxs)` | **구현 완료** |
| **Per-gene KDE heatmap** | PDF slide 62 (개념도) | 없음 | Top 4 genes `imshow(ds.vf[gene_i])` | **구현 완료** |
| **Zoomed cell-type map** | 논문 zoom panels | 없음 | Full + densest ROI crop (500px) | **구현 완료** |
| **h5ad 캐싱** | N/A | N/A | `ssam.h5ad` 캐시 로드 → KDE 재실행 + scanpy/mapping 스킵 | **구현 완료** |

### 5.2 De Novo vs Scanpy Clustering: 이제 병행 실행

**2026-03-09 수정 전**: 파이프라인은 Scanpy Leiden만 사용하고 SSAM 내장 클러스터링을 건너뛰었음.

**2026-03-09 수정 후**: 파이프라인이 **양쪽 모두 실행**:
1. `analysis.cluster_vectors()` + `ds.run_umap()` → SSAM 내장 HDBSCAN 클러스터링 + SSAM UMAP
2. Scanpy `normalize_total → log1p → PCA → neighbors → UMAP → Leiden` → 표준 scRNA-seq 파이프라인

**병행 실행의 이점**:
- `ds.plot_diagnostic_plot(i, use_embedding='umap')`은 SSAM 내장 클러스터링/UMAP이 **반드시 필요**
- `ds.plot_umap(s=1)`도 `ds.run_umap()` 결과가 있어야 작동
- Scanpy Leiden은 `resolution` 파라미터로 클러스터 수를 직접 제어 가능하므로 downstream 분석의 주 클러스터링으로 유지
- SSAM 내장 클러스터링은 진단 시각화 전용, Scanpy Leiden은 실제 분석 결과에 사용

**실행 순서**:
```
normalize_vectors → scale_vectors → cluster_vectors + ds.run_umap (SSAM 내장)
     → AnnData 생성 → normalize_total → log1p → PCA → neighbors → UMAP → Leiden (Scanpy)
```

### 5.3 scale_vectors() 누락 수정 (2026-03-09)

**문제**: SSAM 튜토리얼에서는 `normalize_vectors()` 후 `scale_vectors()`를 별도로 호출한다. 파이프라인에서는 `scale_vectors()`가 누락되어 있었다.

**수정 내용** (line 1441-1442):
```python
# Before (누락)
analysis.normalize_vectors(normalize_vector=True)

# After (수정)
analysis.normalize_vectors(normalize_vector=True)
analysis.scale_vectors()
```

**영향**: `scale_vectors()`는 정규화된 벡터에 원본 L1 norm을 다시 곱하여 발현 크기 정보를 복원한다. 이 단계 없이도 Scanpy 클러스터링은 동작하지만, `ds.normalized_vectors`의 scale이 원본과 달라져 downstream `map_celltypes()` 상관 계산에 영향을 줄 수 있다.

### 5.4 이전 미구현 → 구현 완료 기능 (2026-03-09)

#### 5.4.1 L1 Norm + Local Maxima 시각화 → **구현 완료**

**SSAM 튜토리얼 코드** (Cell 12):
```python
ds.plot_l1norm()
ds.plot_localmax(s=0.1, c='red')
```

**파이프라인 구현** (`_ssam_plots()`, lines 1687-1708):
```python
ds.plot_l1norm(ax=ax)          # L1 norm heatmap
ds.plot_localmax(s=0.3, c='red', ax=ax)  # local maxima overlay
# TypeError fallback: ax= 미지원 SSAM 버전 대응
```
- 출력: `{tag}_step2_ssam_l1norm_localmax.png`
- QC 용도: KDE bandwidth 적절성, local maxima 수 vs 예상 세포 수, 조직 밀도 패턴 검증

#### 5.4.2 SSAM Native UMAP → **구현 완료**

**SSAM 튜토리얼 코드** (Cell 16):
```python
ds.plot_umap(s=1)
```

**파이프라인 구현** (`_ssam_plots()`, lines 1710-1729):
```python
ds.plot_umap(s=1, ax=ax)
```
- 출력: `{tag}_step2_ssam_native_umap.png`
- SSAM 내장 UMAP은 `ds.run_umap()`으로 계산. Scanpy UMAP과 별도 임베딩.
- 전제 조건: `analysis.cluster_vectors()` + `ds.run_umap()`이 먼저 실행되어야 함

#### 5.4.3 Diagnostic Plots → **구현 완료**

**SSAM 튜토리얼 코드** (Cell 20):
```python
for i in range(15):
    plt.figure(figsize=[30, 5])
    ds.plot_diagnostic_plot(i, use_embedding='umap')
```

**파이프라인 구현** (`_ssam_plots()`, lines 1743-1759):
```python
n_diag = min(15, len(ssam_adata.obs['leiden'].cat.categories))
for i in range(n_diag):
    ds.plot_diagnostic_plot(i, use_embedding='umap')
```
- 출력: `{tag}_step2_ssam_diagnostic_{i}.png` (최대 15개)
- 3패널: 공간 분포 | 유전자 발현 bar | UMAP highlight
- 전제 조건: `analysis.cluster_vectors()` + `ds.run_umap()`이 먼저 실행되어야 함

#### 5.4.4 PCA Plot → **구현 완료**

**파이프라인 구현** (`_ssam_plots()`, lines 1731-1741):
```python
sc.pl.pca(ssam_adata, show=False)
```
- 출력: `{tag}_step2_ssam_pca.png`
- Embedding QC: variance ratio 확인

#### 5.4.5 sc.pp.log1p 전처리 → **수정 완료**

**SSAM 튜토리얼** (Cell 25): `sc.pp.normalize_total → sc.pp.log1p → sc.tl.pca`

**수정 전 파이프라인**: `normalize_total → PCA` (log1p 누락)
**수정 후 파이프라인** (line 1474-1475):
```python
sc.pp.normalize_total(ssam_adata, target_sum=1e4)
sc.pp.log1p(ssam_adata)  # 추가
sc.tl.pca(ssam_adata, svd_solver='arpack')
```
- **영향**: log1p 없이 PCA하면 고발현 유전자가 지배적 → 클러스터링 품질 저하
- **연쇄 수정**: `_ssam_map_celltypes()`에서 leiden signature에 이중 log1p 방지 (ssam_adata.X가 이미 log-space)

#### 5.4.6 3D KDE Surface Plot → **구현 완료**

**출처**: 박정빈 SSAM 강의 PDF slide 61 — Gaussian KDE의 3D wireframe surface

**파이프라인 구현**:
```python
from mpl_toolkits.mplot3d import Axes3D
norm_2d = ds.vf_norm.compute().squeeze()
# Subsample to 200x200 (원본 3925x4338은 3D 렌더링에 너무 큼)
step_x = max(1, norm_2d.shape[0] // 200)
sub = norm_2d[::step_x, ::step_y]
ax3d.plot_wireframe(X, Y, sub, rstride=1, cstride=1, linewidth=0.3)
```
- 출력: `{tag}_step2_ssam_kde_3d.png`
- `ds.vf_norm`은 SSAM zarr 저장소에서 L2 norm of vector field를 반환
- 서브샘플링으로 렌더링 속도 보장 (~200x200 = 40K points)

#### 5.4.7 다운샘플링 비교 Plot → **구현 완료**

**출처**: 박정빈 SSAM 강의 PDF slide 63 — "다운샘플링: full vector field vs L1 maxima"

**파이프라인 구현**:
```python
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
ax1.imshow(norm_2d.T, origin='lower', cmap='magma')     # Full field
ax2.scatter(local_maxs[0], local_maxs[1], s=0.3)         # Maxima only
```
- 출력: `{tag}_step2_ssam_downsampling.png`
- 좌: 전체 KDE field (수천만 벡터), 우: local maxima만 (수만 pseudo-cells)
- SSAM의 핵심 효율성을 시각적으로 보여줌

#### 5.4.8 Per-Gene KDE Heatmap → **구현 완료**

**출처**: 박정빈 SSAM 강의 PDF slide 62 — "33 종류의 유전자를 쌓아올림 → 33차원 벡터필드"

**파이프라인 구현**:
```python
for gene_i in top_4_genes:
    gene_field = ds.vf[gene_i].compute().squeeze()
    axes[idx].imshow(gene_field.T, origin='lower', cmap='hot')
```
- 출력: `{tag}_step2_ssam_gene_kde.png`
- 발현량 상위 4개 유전자의 개별 KDE 공간 분포
- 벡터필드가 다차원임을 시각적으로 보여줌 (유전자마다 다른 공간 패턴)

#### 5.4.9 Zoomed Cell-Type Map → **구현 완료**

**파이프라인 구현**:
```python
# 가장 밀도 높은 500x500 ROI 자동 탐지
from scipy.ndimage import uniform_filter
smoothed = uniform_filter(norm_2d, size=500)
peak = np.unravel_index(np.argmax(smoothed), smoothed.shape)
# Full map + ROI rectangle (좌) / Zoomed (우)
```
- 출력: `{tag}_step2_ssam_celltype_map_zoomed.png`
- 전체 map에서 ROI 위치를 흰색 점선으로 표시
- Zoom에서 개별 cell type 영역의 경계 품질을 확인 가능

#### 5.4.10 h5ad 기반 캐싱 → **구현 완료**

**문제**: SSAM KDE 계산이 354개 유전자에 대해 ~10분 소요. 시각화만 수정할 때도 매번 재계산.

**해결 (최종)**: h5ad 파일 기반 캐싱 — scanpy pipeline + reference mapping 결과를 `{sample_tag}_step2_ssam.h5ad`로 저장하고, 이후 실행 시 캐시를 로드하여 scanpy/mapping을 스킵.

> **참고**: 초기에는 zarr store 기반 캐싱 (`SSAMDataset(store=path)`)을 시도했으나, zarr v3에서 `DirectoryStore` 제거 및 `numpy.int64` shape 파싱 에러로 실패. h5ad 기반으로 전환.

```python
ssam_h5ad_path = os.path.join(output_dir, f"{sample_tag}_step2_ssam.h5ad")
if os.path.exists(ssam_h5ad_path):
    # 캐시 로드: scanpy pipeline + reference mapping 결과
    ssam_adata = sc.read_h5ad(ssam_h5ad_path)

    # KDE는 재실행 필요 (ds.plot_* 호출에 필요)
    analysis.run_kde(locations=df, ..., re_run=True)
    analysis.find_localmax(...)
    analysis.normalize_vectors()
    analysis.scale_vectors()
    analysis.cluster_vectors()  # ds.plot_diagnostic_plot 용
    ds.run_umap()               # ds.plot_umap 용

    # 시각화만 실행 (scanpy pipeline + mapping 스킵)
    _ssam_plots(ssam_adata, ds, mapping_result, output_dir, sample_tag)
    return adata
```

**캐시 동작**:
- **첫 실행**: ~22분 (KDE ~10분 + scanpy + mapping + plots) → `ssam.h5ad` 저장
- **이후 실행**: ~7분 (KDE ~10분 재실행, scanpy/mapping 스킵, plots만) — KDE는 `ds.plot_*` API가 `ds.vf` 데이터를 요구하므로 재실행 불가피
- **캐시 위치**: `{output_dir}/{sample_tag}_step2_ssam.h5ad`
- **캐시 내용**: `ssam_adata.X` (log1p-normalized), `.obs['leiden']`, `.obs['leiden_assignment']`, `.obsm['X_umap']`, `.uns['ssam_leiden_celltype_correlations']`

**zarr 캐싱이 실패한 이유**:
1. zarr v3에서 `zarr.DirectoryStore` 클래스가 제거됨 (SSAM 내부에서 사용)
2. `zarr.storage.LocalStore`로 대체 시 `numpy.int64` 값이 shape tuple에 포함되어 `parse_shapelike` 에러 발생
3. SSAM v1.1.3의 zarr v2 종속성이 zarr v3과 근본적으로 비호환

### 5.5 미구현 기능 (의도적 생략)

#### 5.5.1 Watershed Segmentation (선택적 — 후처리 단계)

**SSAM 튜토리얼 코드**:
```python
from skimage import filters
vfnorm_threshold = filters.threshold_local(ds.vf_norm, 35)
vfnorm_thresh_im = (ds.vf_norm > vfnorm_threshold).squeeze().compute()
analysis.run_watershed(vfnorm_thresh_im)
ds.plot_watershed_celltypes_map()
```

**의미**: Vector field norm에 adaptive threshold를 적용하여 세포 영역을 정의하고, watershed 알고리즘으로 개별 세포를 segmentation한다. SSAM의 "segmentation-free" 분석 결과를 활용하여 역으로 세포 경계를 추정하는 후처리.

**파이프라인에서의 필요성**: 낮음. 논문의 주요 결과는 segmentation-free 분석 자체이며, watershed는 SSAM 결과의 검증 목적. Step 3 (Cellpose resegmentation)이 이미 segmentation을 담당.

#### 5.5.2 Domain Map (SSAM Step 4) (선택적)

**SSAM 논문 방법**:
- Cell-type map을 spatial window로 분할
- Window별 cell-type composition 벡터 계산
- Agglomerative clustering → tissue domain 식별

**파이프라인에서의 필요성**: 낮음. Step 5 (Optimal Expansion)에서 Points2Regions 기반 domain 분석이 이미 수행됨.

### 5.6 API 버전 차이 상세

| 항목 | Notebook (SSAM v1.0) | Pipeline (SSAM v1.1+) |
|------|---------------------|----------------------|
| Dataset 생성 | `SSAMDataset(genes, mrna_loci, height, width)` | `SSAMDataset()` (빈 생성) |
| KDE 실행 | `analysis.run_fast_kde(bandwidth=2.5, use_mmap=False, re_run=True)` | `analysis.run_kde(locations=df, width=W, height=H, depth=1, bandwidth=bw_px, sampling_distance=1.0, re_run=True)` |
| 데이터 전달 | 생성자에서 직접 전달 | `run_kde()`의 `locations` DataFrame으로 전달 |
| zarr 백엔드 | zarr v2 | zarr v3 (호환 패치 필요) |
| 유전자 이름 접근 | `ds.genes` | `ds.zarr_group['genes'][:]` |

**zarr v3 호환 패치** (`_patch_zarr_v3_compat()`, lines 1336-1357):
SSAM v1.1.3은 내부적으로 zarr v2 API (`zarr_group.array(name=..., data=...)`)를 사용한다. zarr v3에서는 이 API가 변경되어 `shape`/`dtype`를 필수로 요구한다. 파이프라인은 `zarr.Group.array`를 `create_array`로 리다이렉트하는 monkey-patch를 적용하여 호환성을 확보한다.

---

## 6. Notebook vs Pipeline 구현 비교 (전체)

### 6.1 원본 Notebook들
| 노트북 | Pipeline 함수 | 상태 |
|--------|--------------|------|
| `2_1_batch_processing_distance_to_nuclei_across_samples.ipynb` | `calculate_distance_to_centroid()` + `plot_centroid_distance_analysis()` | 구현됨 |
| `2_3_brain_cell_overlaps.ipynb` | `run_ovrlpy_analysis()` | 구현됨 (선택적) |
| `2_4_brain_ssam.ipynb` | `run_ssam_analysis()` | 구현됨 (선택적) |
| `DEPRECATED_2_2_spots2region_subcellular_exploration_of_output_100cls.ipynb` | 폐기 | 폐기됨 |
| `points2regions/run_p2r.ipynb` | `classify_p2r_subcellular()` | 구현됨 |
| `points2regions/compute_distance.ipynb` | `calculate_distance_to_boundary()` | 구현됨 |
| `points2regions/compute_colors.ipynb` | `assign_p2r_colors()` | 구현됨 |
| `points2regions/figures.ipynb` | `plot_p2r_heatmap()`, `plot_p2r_top_genes()` | 구현됨 |

### 6.2 구현 상태 요약

| 기능 | Notebook | Pipeline | 비고 |
|------|----------|----------|------|
| P2R clustering | 100 cls | 다중 [50,100,200,500] | 개선 |
| P2R subcellular | nuclear/cyto 수동 | `overlaps_nucleus > 0.5` 자동 | 구현됨 |
| Cell type 매핑 | `cell_id.map(cell_id_to_class)` 직접 | `_normalize_cell_ids()` + 3단계 fallback | 개선 (dtype 안전) |
| Centroid distance | `xb.calculating.dispersion()` | 직접 계산 | 동일 결과 |
| Boundary distance (EDT) | 별도 노트북 | signed EDT 통합 | 개선 (signed 부호) |
| Boundary spatial overlay | matplotlib 수동 | `plot_boundary_spatial_overlay()` | 통합 |
| ovrlpy | 별도 노트북 | 선택적 통합 + per-cell incoherence | 구현됨 (확장) |
| SSAM | 별도 노트북 | 선택적 통합 + scRNAseq ref 매핑 | 구현됨 |

---

## 7. SSAM 파라미터 상세 가이드

### 7.1 현재 설정값

```yaml
segmentation_free:
  ssam:
    ssam_vf_um_per_px: 2.0     # Vector field 해상도 (µm/pixel)
    kde_bandwidth_um: 2.5       # KDE bandwidth (µm)
    norm_thres: 5               # 최소 expression L1 norm
    exp_thres: 0.2              # 최소 expression threshold
    min_dist: 7                 # 최소 local maxima 간격 (µm)

sc_reference:
  auto_download: true
  url: "https://sea-ad-single-cell-profiling.s3.us-west-2.amazonaws.com/..."
  dest_dir: "data/scRNAseq"
  celltype_key: "subclass_label"
```

### 7.2 파라미터별 의미와 튜닝

| 파라미터 | 현재값 | 범위 | 영향도 | 상세 설명 |
|----------|--------|------|--------|----------|
| `ssam_vf_um_per_px` | 2.0 | 1.0-5.0 | MEDIUM | Vector field의 공간 해상도. 1.0 µm/px = 고해상도(느림, 메모리 多), 5.0 = 저해상도(빠름). 2.0이 계산 비용과 품질의 균형점. SSAM 튜토리얼에서도 비슷한 해상도 사용. |
| `kde_bandwidth_um` | 2.5 | 1.0-10.0 | **HIGH** | KDE Gaussian kernel의 bandwidth. **가장 중요한 파라미터**. SSAM 논문: "FWHM ≈ 2×bandwidth 가 세포 직경과 비슷해야 함". 뇌 조직 세포 ~10 µm → bw 2.5-5.0 적절. 작은 값 = 개별 세포 해상도, 큰 값 = 조직 구조 수준. |
| `norm_thres` | 5 | 1-20 | LOW | Local maxima의 최소 L1 norm. 낮은 발현 영역의 pseudo-cell을 필터링. 높이면 배경 노이즈 감소, 낮추면 민감도 증가. |
| `exp_thres` | 0.2 | 0.05-0.5 | LOW | 개별 유전자의 최소 발현 수준. SSAM 내부에서 저발현 유전자 필터링에 사용. |
| `min_dist` | 7 | 3-15 µm | LOW | Local maxima 간 최소 거리. 높이면 적은 pseudo-cell 탐지 (빠름, 저해상도). 세포 직경의 50-100% 권장. |
| Leiden `resolution` | 2.0 | 0.5-5.0 | **HIGH** | (코드 내 하드코딩) Scanpy Leiden resolution. 높을수록 더 많은 클러스터. 논문 노트북도 2.0 사용. 조직의 세포 유형 다양성에 따라 조정. |
| `n_neighbors` | 15 | 5-30 | MEDIUM | (코드 내 하드코딩) Scanpy neighbor graph의 k. 높이면 더 넓은 이웃 고려 → 부드러운 클러스터링. |
| `filter min_norm` | 3 | 1-10 | LOW | Pixel-level map 필터링의 최소 norm threshold. |
| `filter min_r` | 0.2 | 0.1-0.5 | MEDIUM | Pixel-level map 필터링의 최소 Pearson correlation. 높이면 더 엄격한 cell type 할당 (빈 영역 증가). SSAM 튜토리얼은 0.3 사용. |

### 7.3 SSAM 논문 vs 파이프라인 파라미터 비교

| 파라미터 | SSAM 논문 추천 | 박정빈 튜토리얼 | 노트북 | 파이프라인 |
|----------|---------------|----------------|--------|----------|
| KDE bandwidth | h=2 (pixels) → FWHM≈10µm | 기본값 사용 | 2.5 µm | 2.5 µm → 1.25 px |
| Grid resolution | 조직 크기에 맞춤 | 1330×1220 px | 원본 좌표 | um_per_px=2.0 으로 계산 |
| Clustering | HDBSCAN (자동) | HDBSCAN | Leiden res=2.0 | Leiden res=2.0 |
| Reference mapping | centroid 기반 Pearson | min_r=0.3 | min_r=0.2 | min_r=0.2 |
| normalize + scale | 둘 다 필수 | 둘 다 호출 | normalize만 | 둘 다 호출 (수정됨) |

---

## 8. 전체 시각화 목록 (Step 2 전체)

| # | 파일명 | 유형 | 논문 Figure | 분석법 |
|---|-------|------|-----------|--------|
| 1 | `p2r_classification.csv` | 데이터 | Ext Fig 3b | P2R 클러스터별 nuclear/cyto 분류 |
| 2 | `p2r_celltype_heatmap.png` | Heatmap | Ext Fig 3b | 행=P2R cluster, 열=cell type, 값=비율 |
| 3 | `p2r_topgenes_{ct}.png` | Stacked bar | Ext Fig 3d | 각 subcellular cluster의 상위 유전자 |
| 4 | `p2r_gene_heatmap_{ct}.png` | Heatmap | - | 유전자×클러스터 발현 매트릭스 |
| 5 | `coherence_map.png` | Spatial map | Fig 1e | z-coherence: 낮은=overlapping |
| 6 | `ovrlpy_pseudocells.png` | Scatter | - | Pseudo-cell coherence |
| 7 | `ovrlpy_tissue.png` | Spatial | - | 조직 레벨 coherence |
| 8 | `ovrlpy_doublets.csv` | 데이터 | - | Overlapping cell 목록 |
| 9 | `cell_incoherence.csv` | 데이터 | - | 세포별 incoherence 점수 |
| 10 | `incoherence_hist.png` | Histogram | - | Incoherence 분포 |
| 11 | `centroid_stripplot.png` | Stripplot | Fig 1f | 극단 유전자의 거리 분포 |
| 12 | `centroid_boxplot.png` | Boxplot | - | 거리 카테고리별 분포 |
| 13 | `centroid_ecdf.png` | ECDF | Fig 2f | 누적 거리 분포 |
| 14 | `centroid_welchs_heatmap.png` | Heatmap | - | Domain 간 거리 차이 |
| 15 | `gene_distance_stats.csv` | 데이터 | - | 유전자별 거리 통계 |
| 16 | `p2r_boundary_{ct}.png` | Boxplot | Ext Fig 3c, Fig 1j | 핵 경계 거리 per P2R cluster |
| 17 | `p2r_mean_distance_scatter.png` | Scatter | - | P2R 클러스터별 평균 signed distance |
| 18 | `distance_boundary_dist.png` | Histogram | - | 2패널: signed EDT + absolute, log scale |
| 19 | `distance_centroid_dist.png` | Histogram | - | Centroid 거리 분포, log scale |
| 20 | `boundary_overview.png` | Spatial+overlay | - | 2패널 전체 뷰 (100K transcripts, 1K boundaries) |
| 21 | `boundary_zoomed_roi.png` | Zoomed overlay | - | 3패널 500µm ROI |
| 22 | `centroid_welchs_pvals.csv` | 매트릭스 | - | Welch's t-test p-value 매트릭스 |
| 23 | `centroid_plotting_subset.csv` | 데이터 | - | 극단 유전자 transcript 서브셋 |
| 24 | `points2regions_k{nc}_bins.h5ad` | AnnData | - | P2R bin-level 데이터 |
| **SSAM 시각화** | | | | |
| 25 | `ssam.h5ad` | AnnData | - | SSAM pseudo-cell 데이터 |
| 26 | `ssam_umap.png` | UMAP | Ext Fig 2c | SSAM Leiden cluster UMAP |
| 27 | `ssam_spatial.png` | Spatial map | Ext Fig 2a | SSAM cell type spatial map |
| 28 | `ssam_highest_genes.png` | Bar | - | Cluster별 상위 유전자 |
| 29 | `ssam_celltype_map.png` | Pixel map | Ext Fig 2a | Pixel-level cell type map (ds.plot) |
| 30 | `ssam_umap_celltypes.png` | UMAP | - | leiden_assignment UMAP |
| 31 | `ssam_correlation.png` | Heatmap | - | Leiden × Reference Pearson correlation |
| 32 | `ssam_l1norm_localmax.png` | Spatial overlay | - | L1 Norm KDE + pseudo-cell local maxima (red) |
| 33 | `ssam_native_umap.png` | UMAP | - | SSAM 내장 UMAP (`ds.plot_umap`, scanpy UMAP과 별도) |
| 34 | `ssam_pca.png` | PCA | - | PCA variance / embedding QC |
| 35 | `ssam_diagnostic_{i}.png` | 3-panel diagnostic | - | Per-cluster: spatial + top genes + UMAP highlight (×15) |
| 36 | `ssam_kde_3d.png` | 3D wireframe | PDF slide 61 | KDE density field 3D surface (subsampled) |
| 37 | `ssam_downsampling.png` | 2-panel comparison | PDF slide 63 | Full vector field vs L1 maxima scatter |
| 38 | `ssam_gene_kde.png` | 4-panel heatmap | PDF slide 62 | Top 4 유전자별 KDE 공간 분포 |
| 39 | `ssam_celltype_map_zoomed.png` | Spatial + zoom | - | Full cell type map + ROI 확대 (Reference 필요) |

---

## 9. Input / Output 상세

### 9.1 Input
| 파일 | 형식 | 출처 |
|------|------|------|
| `{sample_tag}_step1_exploration.h5ad` | AnnData | Step 1 |
| `{sample_tag}_transcripts.parquet` | Parquet | Step 0 |
| Nucleus boundaries (선택) | CSV/Parquet | Xenium 원본 |
| scRNA-seq Reference (선택) | h5ad | SEA-AD (auto-download) |

### 9.2 Output
| 파일 | 형식 | 설명 |
|------|------|------|
| `{sample_tag}_step2_points2regions.h5ad` | AnnData | P2R 도메인 정보 포함 |
| `{sample_tag}_step2_ssam.h5ad` | AnnData | SSAM cell type 포함 (선택) |
| `{sample_tag}_step2_points2regions_k{n}_bins.h5ad` | AnnData | Bin-level 데이터 |
| 31개 시각화/데이터 파일 | PNG/CSV | 분석 결과 (상세: 섹션 8) |

---

## 10. 관련 설정값 정리

```yaml
segmentation_free:
  run_points2regions: true       # P2R 실행 (핵심 분석)
  run_ssam: false                # SSAM 실행 (pip install ssam 필요)
  run_overlaps: false            # ovrlpy 실행 (z_location 필요)
  distance_metric: "both"        # centroid/boundary/both

  points2regions:
    sigma: 3.0                   # Gaussian blur σ
    n_clusters: [50, 100, 200, 500]
    bin_width: 1.0
    min_genes_per_bin: 15

  ssam:
    ssam_vf_um_per_px: 2.0
    kde_bandwidth_um: 2.5
    norm_thres: 5
    exp_thres: 0.2
    min_dist: 7

  overlaps:
    um_per_pixel: 2.0
    radius: 50.0
    bw: 1

sc_reference:
  auto_download: true
  url: "https://sea-ad-single-cell-profiling.s3.us-west-2.amazonaws.com/..."
  dest_dir: "data/scRNAseq"
  celltype_key: "subclass_label"
```

### 10.1 파라미터 조정 가이드

| 파라미터 | 현재값 | 범위/옵션 | 영향도 | 조정 가이드 |
|----------|--------|----------|--------|------------|
| **실행 플래그** | | | | |
| `run_points2regions` | `true` | true/false | - | P2R 핵심 분석. 논문의 주요 결과를 재현하려면 true |
| `run_ssam` | `true` | true/false | - | SSAM 실행. `pip install ssam` 필요. 실행 시간 김 (30분+) |
| `run_overlaps` | `true` | true/false | - | ovrlpy 3D coherence. z_location 컬럼 필요. `pip install ovrlpy` |
| `distance_metric` | `"both"` | centroid/boundary/both | LOW | 계산할 거리 유형. "both" 권장 |
| **P2R 파라미터** | | | | |
| `sigma` | `3.0` | 1.0-10.0 | **HIGH** | Gaussian blur σ (µm). 1.0=세포 수준, 3.0=도메인(기본), 5.0+=대규모 |
| `n_clusters` | `[50,100,200,500]` | 10-1000 리스트 | **HIGH** | KMeans 클러스터 수. 50=대분류, 500=세분류 |
| `bin_width` | `1.0` | 0.5-5.0 µm | MEDIUM | 공간 bin 너비. 작을수록 고해상도(메모리↑) |
| `min_genes_per_bin` | `15` | 5-50 | MEDIUM | 최소 유전자/bin 필터 |
| **SSAM 파라미터** | | | | |
| `ssam_vf_um_per_px` | `2.0` | 1.0-5.0 | MEDIUM | Vector field 해상도. 1.0=고해상도(느림), 5.0=저해상도(빠름) |
| `kde_bandwidth_um` | `2.5` | 1.0-10.0 µm | **HIGH** | KDE bandwidth. 세포 크기와 비례. FWHM≈2×bw≈세포직경 |
| `norm_thres` | `5` | 1-20 | LOW | 최소 발현 norm |
| `exp_thres` | `0.2` | 0.05-0.5 | LOW | Signal threshold |
| `min_dist` | `7` | 3-15 µm | LOW | Sampling point 최소 간격 |
| **ovrlpy 파라미터** | | | | |
| `um_per_pixel` | `2.0` | 1.0-5.0 | MEDIUM | KDE 해상도 |
| `radius` | `50.0` | 20-200 µm | MEDIUM | 분석 반경. 50 µm ≈ 3-5 세포 직경 |
| `bw` | `1` | 1-5 | LOW | KDE bandwidth (pixels) |

---

## 11. 시각화 상세 분석 가이드

### 11.1 P2R Celltype Heatmap (`p2r_celltype_heatmap.png`) → 논문 Extended Data Fig. 3b

**논문 위치**: Extended Data Fig. 3b (p.828) - "P2R cluster × cell type confusion matrix"

**무엇을 봐야 하는가**:
```
              Cell Type:  Oligo  Astro  Neuron  Micro  OPC
    P2R Cluster 0 (nuc)   0.85   0.02   0.05   0.03  0.05
    P2R Cluster 1 (cyto)  0.72   0.10   0.08   0.05  0.05
    P2R Cluster 2 (nuc)   0.05   0.80   0.05   0.05  0.05
    P2R Cluster 3 (cyto)  0.03   0.65   0.15   0.10  0.07

    행 = P2R 클러스터 (nuclear/cytoplasmic 라벨 포함)
    열 = 세포 유형 (Leiden clustering에서 할당)
    값 = 해당 P2R 클러스터의 transcript 중 해당 세포 유형에 속하는 비율

    ① 높은 순도 블록: 하나의 세포 유형에 집중된 P2R 클러스터
    ② 혼재 행: 여러 세포 유형에 분산된 P2R 클러스터 → non-specific domain
    ③ nuclear vs cytoplasmic: 같은 세포 유형의 nuc/cyto 쌍을 비교
```

**해석법**:
- **높은 대각 성분 (>0.7)**: P2R 클러스터가 특정 세포 유형을 잘 대표함
- **낮은 대각 성분 (<0.3)**: P2R 클러스터가 여러 세포 유형의 혼합 영역
- Nuclear 클러스터가 Cytoplasmic 클러스터보다 일반적으로 더 높은 순도를 보임
- 논문: "44 cell-type-specific clusters" 중 nuclear cluster들이 더 세포 유형 특이적

---

### 11.2 P2R Top Genes (`p2r_topgenes_{ct}.png`) → 논문 Extended Data Fig. 3d

**논문 위치**: Extended Data Fig. 3d (p.828)

**무엇을 봐야 하는가**:
```
    각 패널 = 하나의 subcellular cluster (nuclear 또는 cytoplasmic)
    stacked bar = 상위 유전자의 상대적 발현 비율

    ① Nuclear cluster: nuclear-enriched genes (핵 내 mRNA)가 상위
    ② Cytoplasmic cluster: membrane/cytoplasmic genes가 상위
    ③ 같은 세포 유형의 nuclear vs cyto 비교
```

**해석법 (Astrocytes 예시)**:
- **Nuclear cluster**: GFAP, AQP4 등 핵 근처에서 전사되는 유전자가 상위
- **Cytoplasmic cluster**: SLC1A2, GJA1 등 세포 표면/세포질에 localize되는 mRNA가 상위
- **Oligodendrocyte 예시**: Opalin (nuclear), MAG (cytoplasmic/membrane)

---

### 11.3 Coherence Map (`coherence_map.png`) → 논문 Fig. 1e

**논문 위치**: Fig. 1e (p.815)

**해석법**:
- **Coherence > 0.8**: Top/bottom 일관됨 → 정상 세포
- **Coherence < 0.3**: Top/bottom 불일치 → overlapping cells
- 밀집 영역 (cortical layers)에서 더 많은 overlapping
- ovrlpy 결과에서 incoherent cell 비율 5% 이상이면 doublet/overlap 문제 유의

---

### 11.4 Centroid Stripplot (`centroid_stripplot.png`) → 논문 Fig. 1f

**논문 위치**: Fig. 1f (p.815)

**해석법**:
- 상위 5-10 유전자 (centroid에서 먼): cytoplasmic/membrane genes
- 하위 5-10 유전자 (centroid에 가까운): nuclear genes
- **Opalin**: 핵 내 mRNA → 짧은 거리 (논문에서 명시적 언급)
- **MAG, PLP1**: 세포막/myelin sheath → 긴 거리

---

### 11.5 Centroid ECDF (`centroid_ecdf.png`) → 논문 Fig. 2f

**해석법**:
- Nuclear gene ECDF: 5 µm에서 80%+ 도달
- Cytoplasmic gene ECDF: 5 µm에서 40-50%만 도달
- 두 그룹 간의 ECDF 곡선 간격이 클수록 subcellular resolution이 좋음

---

### 11.6 P2R Boundary Distance (`p2r_boundary_{ct}.png`) → 논문 Fig. 1j, Ext. Fig. 3c

**해석법**:
- **Signed EDT convention**: 음수=핵 내부, 양수=핵 외부
- **Nuclear P2R cluster**: 대부분 transcript가 음수 영역
- **Cytoplasmic P2R cluster**: 대부분 transcript가 양수 영역
- 경계 근처 (±2 µm): 핵막 관련 유전자 가능

---

### 11.7 SSAM Spatial Map (`ssam_spatial.png`) → 논문 Extended Data Fig. 2a

**무엇을 봐야 하는가**:
```
    각 점 = SSAM local maximum (pseudo-cell)
    색상 = Reference-mapped cell type (leiden_assignment)
    또는 Leiden cluster (reference 없을 때)

    ① 조직 구조 반영: 뇌의 layer 구조가 보이는지
    ② Leiden 기반 cell type과의 일치: SSAM 결과가 Step 1과 비슷한지
    ③ 해상도: pseudo-cell 밀도가 실제 세포 밀도와 비슷한지
```

**해석법**:
- SSAM은 segmentation bias 없이 조직 구조를 보여줌
- Leiden 클러스터링 결과와 비교하여 두 방법의 일관성 확인
- KDE bandwidth 2.5 µm → 개별 세포 수준의 해상도

---

### 11.8 SSAM Pixel-Level Cell Type Map (`ssam_celltype_map.png`) → 논문 Ext. Fig. 2a

**무엇을 봐야 하는가**:
```
    연속적인 pixel-level cell type map (ds.plot_celltypes_map)
    각 pixel = Pearson correlation 기반 cell type 할당
    filter_celltypemaps(min_norm=3, min_r=0.2)로 low-confidence 제거

    ① 빈 영역: min_r 이하의 correlation → 할당 안 됨 (배경)
    ② 경계 선명도: cell type 간 경계가 명확한지
    ③ 조직 구조 일치: gray/white matter 구분 등
```

**Spatial scatter vs Pixel map 차이**:
- `ssam_spatial.png`: Local maxima (점) 수준, pseudo-cell 위치만 표시
- `ssam_celltype_map.png`: 전체 pixel grid, SSAM의 KDE vector field 기반 연속 맵

---

### 11.9 SSAM UMAP (`ssam_umap.png`, `ssam_umap_celltypes.png`) → 논문 Ext. Fig. 2c

**두 가지 UMAP**:
1. `ssam_umap.png`: Leiden 클러스터 색상 → 클러스터 분리도 확인
2. `ssam_umap_celltypes.png`: Reference-mapped cell type 색상 → 생물학적 의미 확인

**해석법**:
- Leiden UMAP에서 well-separated 클러스터 → 좋은 클러스터링
- Cell type UMAP에서 같은 cell type이 연속적 → 일관된 매핑
- 여러 Leiden 클러스터가 같은 cell type으로 매핑되면 UMAP에서 인접해야 함

---

### 11.10 SSAM Correlation Heatmap (`ssam_correlation.png`) → 파이프라인 추가

**무엇을 봐야 하는가**:
```
    행 = Reference cell types (scRNA-seq subclass_label)
    열 = SSAM Leiden clusters
    값 = Pearson correlation (-0.5 ~ 1.0)
    색상 = RdBu_r (빨강=음의 상관, 파랑=양의 상관)

    ① 각 열(cluster)에서 가장 높은 값 = 해당 cluster의 cell type 할당
    ② 높은 상관 (>0.5): 강한 매칭
    ③ 여러 cell type과 비슷한 상관: 혼합 클러스터 또는 관련 cell type
    ④ 전체적으로 낮은 상관 (<0.3): 유전자 패널 불일치 가능
```

**정상 결과 기준**:
- 대부분의 클러스터에서 최고 상관 > 0.4
- 하나의 reference cell type에 여러 클러스터 매핑 가능 (sub-type 수준 차이)
- 논문 노트북: 50 reference types 중 44가 하나 이상의 클러스터에 매핑됨

---

### 11.11 SSAM Highest Genes (`ssam_highest_genes.png`)

**무엇을 봐야 하는가**:
- 상위 20개 발현 유전자의 비율
- 특정 유전자가 과도하게 지배적이면 KDE가 해당 유전자에 편향될 수 있음
- 미토콘드리아 유전자 (MT-*)가 상위에 있으면 주의 (세포 스트레스/사멸 신호)

---

### 11.12 SSAM L1 Norm + Local Maxima (`ssam_l1norm_localmax.png`) → SSAM 튜토리얼

**출처**: 박정빈 SSAM 튜토리얼의 `ds.plot_l1norm()` + `ds.plot_localmax()` 패턴

**무엇을 봐야 하는가**:
- **L1 Norm 배경**: KDE 밀도 필드의 전체 강도를 보여줌. 밝은 영역 = transcript 밀도 높음
- **빨간 점 (Local Maxima)**: SSAM이 식별한 pseudo-cell 위치. 각 점 = 하나의 pseudo-cell centroid
- **QC 체크**: local maxima가 조직 구조를 따르는지 확인 (gray/white matter 경계에서 밀도 변화)
- **경고 신호**: 빈 영역에 local maxima가 있으면 noise, 조직 내부에 maxima가 없으면 bandwidth 너무 큼

---

### 11.13 SSAM Native UMAP (`ssam_native_umap.png`) → SSAM 튜토리얼 Cell 16

**출처**: 박정빈 SSAM 튜토리얼의 `ds.plot_umap(s=1)` — SSAM 라이브러리 자체 UMAP

**Scanpy UMAP과의 차이**:
- `ssam_umap.png`: Scanpy pipeline (normalize → log1p → PCA → neighbors → UMAP)으로 계산
- `ssam_native_umap.png`: SSAM 내장 `ds.run_umap()`으로 계산, `analysis.cluster_vectors()`의 클러스터 색상 사용
- 두 UMAP은 **다른 알고리즘 경로**이므로 embedding 형태가 다를 수 있음

**무엇을 봐야 하는가**:
- SSAM 자체 클러스터링(`cluster_vectors()`)이 Scanpy Leiden과 유사한 패턴을 보이는지 교차 검증
- SSAM 내장 UMAP에서 클러스터 분리가 더 잘 되거나 못 되면, 전처리 차이(log1p 유무 등)의 영향

---

### 11.14 SSAM PCA (`ssam_pca.png`)

**무엇을 봐야 하는가**:
- PCA variance ratio — 처음 몇 PC가 전체 variance의 상당 부분을 설명해야 함
- Elbow 지점이 명확하면 데이터의 주요 변동 축이 잘 포착된 것
- PC1이 variance의 >30% 설명 시 batch effect 또는 dominant cell type 확인 필요

---

### 11.15 SSAM Per-Cluster Diagnostic (`ssam_diagnostic_{i}.png`) → SSAM 튜토리얼

**출처**: 박정빈 SSAM 튜토리얼의 `ds.plot_diagnostic_plot(i, use_embedding='umap')` 패턴

**3-panel 구성**:
1. **좌측 (Spatial)**: 해당 클러스터의 공간 분포 — 특정 조직 영역에 집중되어야 함
2. **중앙 (Gene expression)**: 해당 클러스터의 상위 발현 유전자 bar chart — marker gene 확인
3. **우측 (UMAP)**: 해당 클러스터가 UMAP 상에서 하이라이트 — 다른 클러스터와 분리도 확인

**해석법**:
- 각 클러스터가 고유한 공간 패턴 + 고유한 유전자 조합을 보여야 정상
- 공간적으로 산재되고 특이 유전자가 없으면 noise cluster일 가능성
- 최대 15개 클러스터만 출력 (파일 bloat 방지)

---

### 11.16 SSAM 3D KDE Surface (`ssam_kde_3d.png`) → PDF Slide 61

**출처**: 박정빈 SSAM 강의 PDF slide 61 — Kernel Density Estimation 3D wireframe

**구현 방법**:
```python
from mpl_toolkits.mplot3d import Axes3D
norm_2d = ds.vf_norm.compute().squeeze()  # L2 norm of vector field
# Subsample to 200x200 for rendering performance
ax3d.plot_wireframe(X, Y, sub, rstride=1, cstride=1, linewidth=0.3)
```

**무엇을 봐야 하는가**:
- **높은 peak**: 특정 위치에 transcript가 집중 → 세포가 있는 위치
- **평탄한 영역**: transcript가 희박한 곳 → 세포 간 공간 또는 조직 외부
- **peak의 너비**: bandwidth 설정에 따라 달라짐 (h=2.5µm → 좁은 peak = 개별 세포 수준)
- **전체 형태**: 조직의 3D 밀도 profile을 직관적으로 보여줌

---

### 11.17 SSAM 다운샘플링 비교 (`ssam_downsampling.png`) → PDF Slide 63

**출처**: 박정빈 SSAM 강의 PDF slide 63 — "다운샘플링"

**2-panel 구성**:
1. **좌 (Full Vector Field)**: L1/L2 norm heatmap, 전체 KDE 그리드 (예: 3925×4338 = 17M vectors)
2. **우 (L1 Maxima)**: Local maxima 위치만 scatter plot으로 표시 (예: 30,000 pseudo-cells)

**무엇을 봐야 하는가**:
- **압축 비율**: 17M → 30K = ~570배 다운샘플링. 이것이 SSAM의 핵심 효율성
- **maxima 분포**: 조직 구조(gray/white matter)가 보존되는지 확인
- **밀도 차이**: heatmap에서 밝은 영역에 maxima가 집중되어야 정상
- **경고 신호**: maxima가 한쪽에 편향되면 threshold 또는 bandwidth 조정 필요

---

### 11.18 SSAM Per-Gene KDE (`ssam_gene_kde.png`) → PDF Slide 62 개념

**출처**: 박정빈 SSAM 강의 PDF slide 62 — "33종류의 유전자를 쌓아올림 → 33차원 벡터필드"

**4-panel 구성**: 발현량 상위 4개 유전자의 개별 KDE heatmap

**무엇을 봐야 하는가**:
- **유전자별 공간 패턴**: 각 유전자가 조직의 특정 영역에 집중되는지 확인
- **Marker gene 검증**: 예를 들어 MBP (oligodendrocyte marker)가 white matter에 집중, GFAP (astrocyte)가 gray matter에 분포
- **Cross-contamination 체크**: 모든 유전자의 KDE가 동일 패턴이면 batch effect 의심

---

### 11.19 SSAM Zoomed Cell-Type Map (`ssam_celltype_map_zoomed.png`)

**2-panel 구성**:
1. **좌 (Full map)**: 전체 pixel-level cell type map + ROI 위치 (흰색 점선 사각형)
2. **우 (Zoomed ROI)**: 가장 밀도 높은 영역 500×500 px 확대

**무엇을 봐야 하는가**:
- **세포 경계**: Zoom에서 개별 cell type 영역이 명확한 경계를 가지는지
- **혼합 영역**: 서로 다른 cell type이 자연스럽게 인접하는지 (biologically plausible)
- **노이즈**: 단일 pixel이 주변과 다른 cell type인 경우 → `filter_celltypemaps` threshold 조정 필요
- **빈 영역(검정)**: min_r threshold에 의해 필터링된 불확실 영역 (적절한 비율: <30%)

---

### 11.20 Boundary Overview/Zoomed (`boundary_overview.png`, `boundary_zoomed_roi.png`)

**해석법**:
- Nuclear mask (cyan contour)가 DAPI signal과 정확히 겹치는지 확인
- Transcript (점)이 핵 내부/외부에 적절히 분포하는지 확인
- 핵 외부의 transcript가 적절한 거리 내에 있는지 (5-10 µm 이내)

---

## 12. 논문 Figure 직접 대응 종합

### 12.1 파이프라인 출력 → 논문 Figure 매핑

| # | 파이프라인 출력 | 논문 Figure | 페이지 | 원문 설명 |
|---|----------------|-----------|--------|----------|
| 1 | `p2r_heatmap_k{n}.png` | Ext. Fig. 3b | 828 | P2R spatial clustering heatmap |
| 2 | `p2r_spatial_k{n}.png` | Fig. 1d | 815 | P2R cluster spatial map |
| 3 | `coherence_3d.png` | Fig. 1e | 815 | 3D coherence of z coordinates |
| 4 | `centroid_strip_plot.png` | Fig. 1f (좌) | 815 | Transcript-centroid distance per gene |
| 5 | `centroid_ecdf.png` | Fig. 2f | 816 | Cumulative distance to centroid |
| 6 | `p2r_boundary_{ct}.png` | Fig. 1j, Ext. Fig. 3c | 815, 828 | Distance to nuclear edge per P2R cluster |
| 7 | `p2r_topgenes_{ct}.png` | Ext. Fig. 3d | 828 | Top genes per subcellular cluster |
| 8 | `ssam_spatial.png` | Ext. Fig. 2a | 827 | SSAM de novo cell type map |
| 9 | `ssam_celltype_map.png` | Ext. Fig. 2a | 827 | SSAM pixel-level cell type map |
| 10 | `ssam_umap.png` | Ext. Fig. 2c | 827 | SSAM cluster UMAP |
| 11 | `ssam_umap_celltypes.png` | - | - | Reference-mapped UMAP (파이프라인 추가) |
| 12 | `ssam_correlation.png` | - | - | Leiden-Reference 상관 (파이프라인 추가) |
| 13 | `ssam_l1norm_localmax.png` | 튜토리얼 Cell 12 | - | L1 Norm + Local Maxima QC (2026-03-09 추가) |
| 14 | `ssam_native_umap.png` | 튜토리얼 Cell 16 | - | SSAM 내장 UMAP (2026-03-09 추가) |
| 15 | `ssam_pca.png` | - | - | PCA embedding QC (2026-03-09 추가) |
| 16 | `ssam_diagnostic_{i}.png` | 튜토리얼 Cell 20 | - | Per-cluster 3-panel diagnostic ×15 (2026-03-09 추가) |

### 12.2 정상 결과 판별 체크리스트

#### SSAM 분석
- [ ] **Local maxima 수**: 예상 세포 수의 50-200% 범위 (brain tissue: 수천-수만)
- [ ] **Leiden 클러스터 수**: 20-80 범위 (resolution=2.0 기준)
- [ ] **Reference 매핑률**: common genes > 50개, 대부분 클러스터의 최고 상관 > 0.3
- [ ] **Spatial map**: 조직 구조 반영 (gray/white matter 구분)
- [ ] **UMAP**: Well-separated 클러스터 (overlapping 최소)
- [ ] **Pixel map**: 빈 영역 < 30% (min_r=0.2 기준)
- [ ] **L1 Norm**: 조직 구조에 맞는 밀도 패턴, local maxima(빨간 점)가 조직 내부에 고르게 분포
- [ ] **PCA**: 상위 PC들이 전체 variance의 상당 부분 설명 (elbow 명확)
- [ ] **Diagnostic plots**: 각 클러스터가 고유 공간 패턴 + 고유 유전자 조합 보유

#### P2R 분석
- [ ] **클러스터-유전자 특이성**: 각 클러스터가 고유 유전자 패턴
- [ ] **Nuclear vs Cyto 분리**: overlaps_nucleus > 0.5 기준으로 명확한 분류
- [ ] **공간 연속성**: 같은 클러스터가 공간적으로 연결

#### Distance 분석
- [ ] **Centroid distance**: Nuclear genes < 5 µm, Cytoplasmic genes > 8 µm (중앙값)
- [ ] **Boundary distance**: Nuclear clusters 중앙값 < 0, Cytoplasmic clusters 중앙값 > 0
- [ ] **ECDF 분리**: Nuclear vs Cytoplasmic gene ECDF 곡선 간 간격 존재

### 12.3 논문 원문 인용 (시각화 관련)

| 시각화 | 논문 원문 인용 | 페이지 |
|--------|--------------|--------|
| P2R / SSAM clustering | "Using one of these approaches, SSAM de novo mode, we identified 44 cell-type-specific clusters" | 814 |
| Z-coherence (Fig. 1e) | "3D coherence of z coordinates across the tissue section" | 815 caption |
| Centroid distance (Fig. 1f) | "we identified some mRNAs enriched in the nucleus, and others in the cytoplasm" | 814 |
| ECDF (Fig. 2f) | "Cumulative proportion of reads by distance from the cell centroid" | 816 caption |
| Subcellular distribution | "Xenium's signal density facilitates the in situ identification of subcellular structures" | 814 |
| Nuclear vs cytoplasmic | "Most nuclear genes" (Opalin, Sox10) vs "Least nuclear genes" (MAG, CAPN2) | 815 Fig.1f |

---

## 13. Summary

Step 2는 **segmentation 없이 공간 패턴을 분석**하는 핵심 단계이다. P2R clustering으로 도메인을 정의하고, distance metrics로 subcellular 분포를 정량화하며, SSAM으로 reference-free cell type mapping을 수행한다.

| 컴포넌트 | 구현 상태 | 비고 |
|---------|---------|------|
| Points2Regions clustering | 완전 | 다중 해상도 지원 |
| P2R subcellular 분류 | 완전 | overlaps_nucleus 기반 |
| Centroid distance | 완전 | 논문 Fig. 1f 대응 |
| Boundary distance (EDT) | 완전 | Signed convention |
| ovrlpy coherence | 완전 (선택적) | z_location 필요 |
| SSAM KDE + local maxima | 완전 | v1.1+ API, zarr v3 패치 |
| SSAM Scanpy clustering | 완전 | Leiden res=2.0 |
| SSAM reference mapping | 완전 | SEA-AD scRNA-seq |
| SSAM pixel-level map | 완전 | map_celltypes + filter |
| SSAM scale_vectors() | **수정 완료** | 2026-03-09 추가 |
| SSAM L1 norm / localmax plot | **구현 완료** | 2026-03-09 추가, QC 진단용 |
| SSAM Native UMAP (`ds.plot_umap`) | **구현 완료** | 2026-03-09 추가, SSAM 내장 UMAP |
| SSAM PCA plot | **구현 완료** | 2026-03-09 추가, embedding QC |
| SSAM per-cluster diagnostic | **구현 완료** | 2026-03-09 추가, 최대 15 클러스터 |
| `sc.pp.log1p` 전처리 | **수정 완료** | 2026-03-09 추가, scanpy 표준 workflow |
| `analysis.cluster_vectors()` | **추가 완료** | 2026-03-09 추가, diagnostic plot 전제조건 |
| `ds.run_umap()` | **추가 완료** | 2026-03-09 추가, SSAM native UMAP 전제조건 |
| SSAM 3D KDE surface | **구현 완료** | 2026-03-09 추가, PDF slide 61 |
| SSAM 다운샘플링 비교 | **구현 완료** | 2026-03-09 추가, PDF slide 63 |
| SSAM Per-Gene KDE heatmap | **구현 완료** | 2026-03-09 추가, PDF slide 62 |
| SSAM Zoomed cell-type map | **구현 완료** | 2026-03-09 추가, ROI 자동 탐지 |
| h5ad 기반 캐싱 | **구현 완료** | 2026-03-09 추가, h5ad 캐시로 scanpy/mapping 스킵 (KDE는 재실행) |
| SSAM diagnostic plots | 미구현 (선택적) | per-celltype 상세 |
| SSAM watershed | 미구현 (선택적) | Step 3에서 대체 |
| SSAM domain map | 미구현 (선택적) | Step 5에서 대체 |

**구현 완성도**: **HIGH** - 논문/튜토리얼의 핵심 기능 모두 구현. 선택적 기능 (diagnostic plots, watershed, domain map)은 다른 Step에서 대체되거나 QC 목적으로만 필요.

**코드 수정 이력**:
- 2026-03-09: `scale_vectors()` 누락 수정 (SSAM 튜토리얼 준수)
