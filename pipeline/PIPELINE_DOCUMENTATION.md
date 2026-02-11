# MASLD Xenium Pipeline Documentation

## 📋 Overview

이 문서는 `/notebooks` 디렉토리의 원본 Jupyter notebooks를 `/pipeline` 디렉토리의 Python 스크립트로 재구현한 **MASLD Xenium 파이프라인**의 전체 프로세스를 설명합니다.

### 핵심 목표
- **Xenium 원시 데이터** → **구조화된 공간 전사체 분석** 파이프라인 자동화
- 7단계 프로세스를 통해 세포 분류, 유전자 발현 분석, 공간 최적화 수행
- 시뮬레이션과 벤치마킹을 통한 방법론 검증

---

## 🔄 Pipeline Architecture Overview

```
Input: Xenium Raw Output (machine output)
    ↓
Step 0: Formatting
    ↓
Step 1: Dataset Exploration
    ↓
Step 2: Segmentation-Free Analysis (Points2Regions)
    ↓
Step 3: Resegmentation (Cellpose)
    ↓
Step 4: Techniques Comparison & Validation
    ↓
Step 5: Optimal Expansion
    ↓
Step 6: Segmentation Benchmark (+ Baysor)
    ↓
Step 7: Simulation & Preprocessing Benchmarking
    ↓
Output: Comprehensive Analysis Results
```

### 데이터 흐름 다이어그램

```
Raw Xenium → [Step 0] → AnnData (h5ad)
                           ↓
                    [Step 1] → Exploration Stats
                           ↓
                    [Step 2] → Points2Regions Domains
                           ↓
                    [Step 3] → Cellpose Segmentation
                           ↓
              ┌─────────────┴─────────────┐
              ↓                           ↓
          [Step 4]              [Step 5] Expansion
          Comparison              KDTree Assignment
              ↓                           ↓
              └─────────────┬─────────────┘
                           ↓
                    [Step 6] Benchmark
                    (Nuclei/Cellpose/Expansion/Baysor)
                           ↓
                    [Step 7] Simulation
                    (scRNAseq → Simulated Xenium)
```

---

## 📊 Steps Detail

### **Step 0: Formatting**
**Notebook Reference:** `0_formatting/` 라인의 `format_xenium_adata_mid_2023()`

#### 목적
Xenium 머신 출력 형식을 표준 AnnData 객체로 변환

#### 입출력
- **입력:**
  - Raw Xenium directory (`cell_feature_matrix.tar.gz`, `cells.csv`, `transcripts.csv` 등)
  - 10X analysis 파일 (UMAP, PCA, clusters)
  - Gene panel JSON

- **출력:**
  - `{sample_tag}.h5ad` (formatted AnnData)
  - `{sample_tag}_transcripts.parquet` (transcript sidecar)
  - Processed DAPI image (TIFF)

#### 주요 기능
1. **데이터 압축 해제**
   - `cell_feature_matrix.tar.gz` → folder 전개
   - `.gz` 파일들 압축 해제

2. **AnnData 구축**
   - Matrix Market 형식 (MTX) 읽기
   - Cell metadata 병합
   - Gene features 파싱 (Gene ID/Name/Reason of Inclusion)

3. **유전자 정규화**
   - Gene panel JSON에서 Ensembl ID 매핑
   - 컨트롤 프로브 제거 (negative/positive controls)

4. **Metadata 통합**
   - 10X UMAP/PCA/Clustering 임포트
   - Transcript 좌표계 보존

5. **QC Filtering**
   - `min_counts` threshold 적용
   - `min_genes` threshold 적용
   - 필터링된 세포 제거

#### 시각화
- 데이터셋 크기 요약 (cells × genes)
- QC 필터링 전후 비교

#### 핵심 코드 패턴
```python
# Decompress cell_feature_matrix
matrix_path = os.path.join(cfm_path, 'matrix.mtx')
a = mmread(matrix_path)
ad = a.todense()

# Parse features dynamically (3-col or 2-col format)
if features_df.shape[1] == 3:
    features_df.columns = ['gene_id', 'gene_name', 'reason_of_inclusion']
else:
    features_df.columns = ['gene_id', 'reason_of_inclusion']
```

---

### **Step 1: Dataset Exploration**
**Notebook Reference:** `1_datasets_exploration/` 라인의 `run_step1()`

#### 목적
데이터셋의 기본 통계, 유전자 분산, 클러스터링, 공간 구조 분석

#### 입출력
- **입력:** Step 0 AnnData + Transcript sidecar
- **출력:**
  - `{sample_tag}_step1_exploration.h5ad`
  - `step1_exploration/` 디렉토리 내 분석 결과

#### 주요 기능

##### 1-1. 일반 통계 (General Statistics)
```
- 총 세포 수, 유전자 수
- 평균 counts per cell
- Sparsity (표현된 유전자 비율)
```

##### 1-2. Transcript Dispersion Analysis
- **Distance Histogram + KDE**
  - 모든 transcript-to-centroid 거리 분포

- **ECDF (Empirical CDF)**
  - 누적 분포 함수 (전체 + 상위 10 유전자 별)

- **Violin Plot**
  - 상위 20 유전자의 거리 분포

- **KS Test**
  - 유전자 쌍 간 거리 분포 차이 검정
  - 유전자 쌍 간 공간 분리도 측정

##### 1-3. Clustering & Marker Annotation
- **PCA + Neighbors**
  - 고변동 유전자 (HVG) 선택
  - k-nearest neighbors 그래프 구축

- **Leiden Clustering**
  - 다중 해상도 (resolutions: 0.5, 1.0, 1.5)

- **Marker Gene Detection**
  - Wilcoxon rank-sum test
  - 클러스터별 상위 marker genes

- **시각화**
  - Dotplot (cluster × top markers)
  - Heatmap (마커 유전자 발현)
  - UMAP 임베딩
  - Spatial scatter plot

##### 1-4. Neighborhood Analysis
- **Spatial Neighbors Graph**
  - Delaunay triangulation 기반 공간 이웃

- **Neighborhood Diversity**
  - Shannon entropy of cell types in neighborhoods

- **Enrichment Analysis**
  - 공간적으로 인접한 세포 유형 쌍 분석

- **Centrality Scores**
  - Betweenness, closeness centrality

#### 시각화
- Distance histograms & KDE plots
- ECDF curves (전체 + 유전자별)
- Violin plots
- Marker gene heatmaps & dotplots
- UMAP embeddings
- Spatial scatter plots
- Neighborhood diversity maps

---

### **Step 2: Segmentation-Free Analysis**
**Notebook Reference:** `2_segmentation_free_analysis/`

#### 목적
세포 경계를 정의하지 않고 공간 transcript 패턴 분석 (Points2Regions)

#### 입출력
- **입력:** Step 1 AnnData + Transcripts
- **출력:**
  - `{sample_tag}_step2_points2regions.h5ad`
  - `points2regions/` 서브디렉토리의 분석 결과

#### 주요 기능

##### 2-1. Points2Regions Clustering
- **도메인 정의**
  - Transcript 좌표 기반 Voronoi/spatial clustering
  - 자동 도메인 수 결정
  - Domain boundaries → GeoJSON 저장

##### 2-2. Distance Metrics
- **Centroid Distance**
  - 각 transcript → domain centroid 거리
  - 계산식: `sqrt((x - centroid_x)^2 + (y - centroid_y)^2)`
  - **주의:** `xb.calculating.dispersion()` = centroid distance (NOT boundary distance)

- **Boundary Distance (EDT)**
  - Euclidean Distance Transform 기반
  - Rasterization → Distance map
  - 각 pixel → 가장 가까운 boundary까지 거리

##### 2-3. 시각화

**극단값 유전자 분석**
```python
# Genes with most extreme dispersion
stripplot + boxplot (최대/최소 centroid distance genes)
```

**거리 기반 색상 코딩 Boxplot**
```python
# Boxplot colored by distance category
- Core (0-25%)
- Intermediate (25-50%)
- Peripheral (50-75%)
- Outer (75-100%)
```

**ECDF Plot**
- 누적 분포 (distance 기반)

**Welch's t-test Heatmap**
- Domain 간 유전자 거리 분포 차이

**Summary Histograms**
- 모든 transcript 거리 분포

#### 핵심 코드
```python
# Points2Regions 도메인 정의
from points2regions import Points2Regions
p2r = Points2Regions(transcript_df, n_regions=auto)
domains = p2r.fit()

# Centroid distance 계산
dist_to_centroid = np.sqrt((transcripts.x - centroid.x)**2 +
                            (transcripts.y - centroid.y)**2)

# EDT (Boundary distance)
edt = distance_transform_edt(binary_mask)
```

---

### **Step 3: Resegmentation (Cellpose)**
**Notebook Reference:** `3_techniques_comparison/3_1, 3_2`

#### 목적
Cellpose를 사용한 새로운 세포 경계 정의 및 transcript 재할당

#### 입출력
- **입력:**
  - DAPI 이미지 (morphology_focus.ome.tif 또는 DAPI.tif)
  - 원본 Transcript CSV
  - Step 2에서 생성된 domain_polygons.json (선택사항)

- **출력:**
  - `{sample_tag}_step3_resegmented.h5ad`
  - `{sample_tag}_step3_transcripts_resegmented.csv`
  - `{sample_tag}_step3_resegmented_masks.tif`
  - `step3_done.txt` (완료 마커)

#### 주요 기능

##### 3-1. Device Detection
```python
# 자동 디바이스 선택: CUDA > MPS (Apple Silicon) > CPU
device = _detect_device()
```

##### 3-2. DAPI 이미지 로딩
- TIF 형식 읽기
- 다중 채널 처리 (필요시 단일 채널 추출)

##### 3-3. Cellpose Nuclei Segmentation
```python
from cellpose import models
model = models.Cellpose(gpu=use_cuda, model_type='nuclei')
masks, flows, styles = model.eval(dapi, diameter=30)
```

주요 파라미터:
- `diameter` : 핵의 예상 직경 (pixels)
- `flow_threshold` : 흐름 기반 분할 임계값
- `cellprob_threshold` : 세포 확률 임계값

##### 3-4. Mask 저장 (TIF)
```python
# 정수 마스크 (각 픽셀 = cell_id)
tifffile.imwrite('masks.tif', masks.astype(np.uint32))
```

##### 3-5. Transcript-to-Cell Mapping
```python
# 각 transcript의 좌표로 마스크 조회
for _, tr in transcripts.iterrows():
    cell_id = masks[int(tr.y_location), int(tr.x_location)]
    transcript_df.loc[_, 'cell_id'] = cell_id
```

##### 3-6. Cell × Gene Matrix 구축
```python
# 재할당된 transcript로부터 cell×gene count matrix 생성
cell_feature_matrix = count_matrix(transcripts_resegmented)
adata = ad.AnnData(X=cell_feature_matrix)
```

##### 3-7. Cell Centroids 추출 (regionprops)
```python
from skimage.measure import regionprops
for region in regionprops(masks):
    adata.obs.loc[region.label, ['centroid_x', 'centroid_y']] = region.centroid
```

##### 3-8. Domain Assignment (선택사항)
- Step 2에서 생성된 `domain_polygons.json` 로드
- 각 cell centroid → domain 할당 (shapely Point-in-Polygon)

#### 시각화
- Cellpose segmentation overlay
- Resegmented cell count distribution
- Domain assignment map

---

### **Step 4: Techniques Comparison & Validation**
**Notebook Reference:** `3_techniques_comparison/3_3-3_7`

#### 목적
원본 vs 재분할 데이터를 효율성(efficiency), 특이성(specificity), 양성률(positivity), 확산(diffusion) 메트릭으로 비교

#### 입출력
- **입력:**
  - Step 3 재분할 데이터 (또는 Step 1 폴백)
  - Step 1 원본 데이터 (비교용)
  - scRNAseq 레퍼런스 (선택사항, specificity 계산용)

- **출력:**
  - `step4_techniques_comparison/` 디렉토리의 비교 결과
  - `step4_done.txt` 마커

#### 주요 기능

##### 4-1. Efficiency Analysis
- **Transcript/Gene per Cell**
  - Histogram: 재분할 vs 원본 세포당 transcript 수
  - Histogram: 세포당 발현 유전자 수

- **Expression Ratio (ST vs scRNAseq)**
  - CPM (Counts Per Million) 정규화
  - Xenium 표현 수준 vs scRNAseq 비교
  - 영역별 효율성 분석 (domain-specific)

##### 4-2. Specificity Analysis
- **Negative Marker Purity**
  - Gene-gene co-expression matrix
  - 특정 유전자 쌍이 동시에 발현되는 세포 비율

- **Gene-Gene Correlation Heatmap**
  - Pearson correlation (원본 vs 재분할)

##### 4-3. Positivity Analysis
- **Positivity Histogram**
  - 각 유전자의 양성 세포 비율 분포

- **Clustering on Positivity**
  - Preprocessing (normalize, log1p, scale)
  - Leiden clustering

- **Violin Plot per Cluster**
  - 클러스터별 양성률 분포

##### 4-4. Diffusion Analysis
- **Transcript-to-Centroid Distance**
  - Pixel → µm 변환 (Xenium: 4.70588 pixel/µm)
  - 원본 vs 재분할 거리 비교

- **Complementary CDF**
  - P(distance > x) 플롯

- **Per-Gene ECDF**
  - 상위 유전자별 누적 분포

- **Gene × Method Distance Heatmap**
  - 유전자 × 방법(원본/재분할) 평균 거리

#### 핵심 코드
```python
# Co-expression 계산
def _coexpression_calculation(exp, min_exp=0):
    for col in exp.columns:
        sel = exp.loc[:, col] > min_exp
        positive_cells = exp.loc[sel, :]
        coexpression.loc[:, col] = np.sum(positive_cells > min_exp) / positive_cells.shape[0]

# Distance 변환 (pixel → µm)
pixel_to_um = 1 / 4.70588  # Xenium pixel size
distances_um = distances_px * pixel_to_um
```

#### 시각화
- Efficiency: Transcript/gene histograms
- Specificity: Co-expression heatmap
- Positivity: Violin plot per cluster
- Diffusion: ECDF, complementary CDF, per-gene distance heatmaps

---

### **Step 5: Optimal Expansion**
**Notebook Reference:** `4_optimal_expansion/4_1`

#### 목적
미할당 transcript를 가장 가까운 주석이 달린 도메인에 할당하고, 상관관계 기반 "turnover" 분석으로 최적 확장 반경 계산

#### 입출력
- **입력:**
  - Step 0 원본 transcript (전체)
  - Step 1/3/4 주석이 달린 세포 데이터

- **출력:**
  - `{sample_tag}_step5_expanded_transcripts.csv`
  - `{sample_tag}_step5_optimal_expansion_summary.csv`
  - `step5_optimal_expansion/` 분석 결과 (plots, CSVs)
  - `{sample_tag}_step5_done.txt` 마커

#### 주요 기능

##### 5-1. Domain Assignment via KDTree
```python
# 모든 원본 transcript를 cell domain에 매핑
from scipy.spatial import cKDTree

tree = cKDTree(cell_centroids)
distances, indices = tree.query(unassigned_transcripts)

# 거리 threshold로 필터링 (선택사항)
assigned = unassigned_transcripts[distances < threshold]
```

##### 5-2. Spatial Mapping
- 할당된 vs 미할당 transcript 시각화
- Distance distribution histogram

##### 5-3. Turnover Analysis
**핵심 개념:** Turnover = background expression이 비핵심 expression과 동일해지는 거리

- **Nuclear vs Background Expression**
  - Core region (domain 내부): "nuclear" profile
  - Outer region (domain 외부): "background" profile

- **Distance-Bin Correlation**
  ```python
  # 거리 범위별로 binning
  # 각 bin에서 이웃 영역과의 유전자 발현 상관성 계산
  correlation_by_distance = [
      pearson_corr(nuclear_profile, background_binX)
      for bin_X in distance_bins
  ]
  ```

- **Turnover Distance Detection**
  - Correlation이 0.5 이상 (또는 설정값)에 도달하는 거리

- **Nuclei Size Measurement**
  ```python
  def dist_nuc(reads_ctdsub):
      """Median distance to nucleus edges via ConvexHull"""
      for cell_id, transcripts in reads.groupby('cell_id'):
          hull = ConvexHull(transcripts[['x_location', 'y_location']])
          nuclei_dist = mean(hull.vertices distances)
  ```

- **Optimal Expansion Calculation**
  ```python
  optimal_expansion = turnover_distance - nuclei_size
  ```

##### 5-4. Summary & Output
- Barplot: 방법별 최적 확장 반경
- CSVs: 유전자별, 도메인별 상세 결과

#### 시각화
- Spatial map with assigned/unassigned transcripts
- Distance to boundary histogram
- Turnover distance vs distance plot
- Optimal expansion barplot (method comparison)

---

### **Step 6: Segmentation Benchmark**
**Notebook Reference:** `5_segmentation_benchmark/`

#### 목적
여러 분할 방법(nuclei, Cellpose, expansion, Baysor)을 하나의 통일된 벤치마크 프레임워크에서 비교

#### 입출력
- **입력:**
  - Step 0: Nuclei segmentation (기준)
  - Step 3: Cellpose resegmentation
  - Step 5: Expansion data
  - Raw Xenium transcripts (Baysor용)
  - (선택) Step 3 TIF masks (Baysor 사전 정보)

- **출력:**
  - `step6_benchmark/` 디렉토리:
    - Concatenated AnnData (모든 방법)
    - 전처리된 벤치마크 데이터
    - 비교 메트릭 (CSV)
    - 시각화 (UMAP, spatial scatter, violin plots)
  - `step6_done.txt` 마커

#### 주요 기능

##### 6-1. Baysor 실행 (선택사항)
```bash
# Xenium 데이터를 Baysor 포맷으로 변환
# Baysor 배치 처리 실행
baysor run -x x -y y -g gene transcripts.csv \
    --prior-segmentation masks.tif \
    -o baysor_output/
```

준비 단계:
- Transcripts CSV → Baysor 포맷 (gene, x, y)
- 필요시 ROI 자르기 (crop)

##### 6-2. 분할 결과 로딩
```python
# 각 방법에서 cell×gene matrix 로드
input_nuclei = load_adata(step0_nuclei.h5ad)      # 기준
input_advanced = load_adata(step3_cellpose.h5ad)  # Cellpose
input_expansion = load_adata(step5_expansion)     # Expansion (CSV → AnnData)
input_baysor = run_baysor(...)                    # Baysor
```

##### 6-3. Concatenation & Preprocessing
```python
# 모든 방법을 batch key로 연결
combined_adata = ad.concat([
    input_nuclei,
    input_advanced,
    input_expansion,
    input_baysor
], axis=0, keys=['nuclei', 'cellpose', 'expansion', 'baysor'])

# 통일된 전처리
sc.pp.normalize_total(combined_adata)
sc.pp.log1p(combined_adata)
sc.pp.highly_variable_genes(combined_adata)
sc.tl.pca(combined_adata)
sc.pp.neighbors(combined_adata)
sc.tl.umap(combined_adata)
sc.tl.leiden(combined_adata)
```

##### 6-4. Annotation Transfer
```python
# 레퍼런스(nuclei)에서 다른 방법으로 cell type 전파
# Majority voting: 각 셀의 이웃 주석으로 할당
for method in ['cellpose', 'expansion', 'baysor']:
    transfer_labels(from=nuclei, to=method, weight='correlation')
```

##### 6-5. Benchmark Metrics
```python
metrics_per_method = {
    'n_cells': len(adata),
    'median_reads_per_cell': adata.obs['n_reads'].median(),
    'assigned_transcripts': sum(adata.obs['n_reads'] > 0),
    'clustering_quality': silhouette_score(...),
    'cell_size_uniformity': cv(adata.obs['n_reads'])
}
```

##### 6-6. Visualizations
- **UMAP per Method**
  - 각 방법별 embedding, cell type 색상화

- **Spatial Scatter**
  - 원래 공간 좌표에서 cell type 분포

- **Cell Type Barplot**
  - 방법별 cell type 구성 비교

- **Transcript Count Violin**
  - Cell type별 reads per cell 분포

#### 핵심 코드
```python
# Baysor 데이터 준비
def prep_xenium_data_for_baysor(xenium_dir, out_dir):
    spots = pd.read_csv(f"{xenium_dir}/transcripts.csv")
    spots_baysor = spots[['feature_name', 'x_location', 'y_location']]
    spots_baysor.columns = ['gene', 'x', 'y']
    spots_baysor.to_csv(f"{out_dir}/transcripts_baysor.csv", index=False)

# Concatenation with batch key
combined = ad.concat(method_list, keys=method_names, axis=0)
combined.obs['method'] = combined.obs.index.get_level_values(0)
```

#### 시각화
- 방법별 UMAP 비교
- 공간 산점도 (cell type 색상)
- Cell type 구성 바플롯
- Transcript count violin plot (클러스터별)

---

### **Step 7: Simulation & Preprocessing Benchmarking**
**Notebook Reference:** `6_simulating_preprocessing/`

#### 목적
scRNAseq 데이터에서 Xenium을 시뮬레이션하고, 다양한 전처리 파라미터 조합을 벤치마킹하여 최적 파이프라인 도출

#### 입출력
- **입력:**
  - (선택) scRNAseq 레퍼런스 데이터 또는 CellxGene Census에서 다운로드
  - 시뮬레이션 설정 (noise level, misseg rate 등)

- **출력:**
  - `step7_simulation/` 디렉토리:
    - `simulated_standard.h5ad` (표준 시뮬레이션)
    - `simulated_high_noise.h5ad` (고 잡음 시뮬레이션)
    - `benchmark_results.csv` (그리드 검색 결과)
    - 벤치마크 시각화 (bar plot, heatmap, 중요도, box plot)
  - `{sample_tag}_step7_done.txt` 마커

#### 주요 기능

##### 7-1. Reference Acquisition
```python
# Option 1: CellxGene Census에서 다운로드
import cellxgene_census
data = cellxgene_census.get_anndata(...)

# Option 2: 로컬 scRNAseq 파일 로드
adata_ref = sc.read_h5ad('reference.h5ad')
```

##### 7-2. Simulation

**Step 7-2a: Subsampling & HVG**
```python
# 크기 조정 (메모리/계산 효율)
adata_sub = adata_ref[adata_ref.obs['cell_type'].value_counts() > min_cells]

# Highly variable genes 선택
sc.pp.highly_variable_genes(adata_sub, n_top_genes=200)
```

**Step 7-2b: Rank Marker Genes**
```python
# Cell type별 marker 검출
sc.tl.rank_genes_groups(adata_sub, groupby='cell_type', method='wilcoxon')
```

**Step 7-2c: Standard Simulation**
```python
# 기본 noise/missegmentation 레벨
simulated_standard = simulate_xenium(
    adata=adata_sub,
    n_reads_per_cell=500,
    dropout_rate=0.1,      # 10% zero inflation
    misseg_rate=0.05,      # 5% misassignment
    noise_level='standard'
)
```

**Step 7-2d: High-Noise Simulation**
```python
# 높은 noise/missegmentation (robust 메서드 테스트)
simulated_high_noise = simulate_xenium(
    ...,
    dropout_rate=0.2,      # 20% zero inflation
    misseg_rate=0.1,       # 10% misassignment
    noise_level='high'
)
```

##### 7-3. Preprocessing Grid Search

**Step 7-3a: 그리드 정의**
```python
preprocessing_grid = {
    'normalize': [True],
    'log1p': [True],
    'hvg_filter': [50, 100, 200],
    'scale': [True],
    'pca_n_comps': [10, 20, 30],
    'leiden_resolutions': [0.3, 0.5, 0.8, 1.0]
}
```

**Step 7-3b: Grid Search 실행**
```python
# 모든 조합 반복 (3 × 3 = 9 combinations)
for hvg, pca_n, res in itertools.product(...):
    adata_proc = preprocess(simulated_standard,
                             hvg_n=hvg,
                             pca_n=pca_n,
                             leiden_res=res)
    metrics = evaluate_clustering(adata_proc)
```

**Step 7-3c: Clustering Quality Metrics**

각 전처리 조합에 대해 다음 메트릭 계산:

1. **Normalized Mutual Information (NMI)**
   - 시뮬레이션 true label과 예측 클러스터의 mutual information
   - 범위: [0, 1], 높을수록 좋음

2. **Adjusted Rand Index (ARI)**
   - True vs 예측 클러스터의 일치도 (adjustment for chance)
   - 범위: [-1, 1], 높을수록 좋음

3. **Fowlkes-Mallows Index (FMI)**
   - Precision and recall의 기하평균
   - 범위: [0, 1], 높을수록 좋음

4. **Variation of Information (VI)**
   - 엔트로피 기반 distance metric
   - 낮을수록 좋음

```python
from sklearn.metrics import (
    normalized_mutual_info_score as NMI,
    adjusted_rand_score as ARI,
    fowlkes_mallows_score as FMI
)

nmi = NMI(true_labels, predicted_labels)
ari = ARI(true_labels, predicted_labels)
fmi = FMI(true_labels, predicted_labels)
vi = VI(true_labels, predicted_labels)
```

##### 7-4. Benchmark Visualizations

**7-4a: Bar Plot (Overall Performance)**
```
x축: 전처리 조합
y축: 평균 메트릭 점수 (NMI, ARI, FMI 정규화)
```

**7-4b: Heatmap (Parameter Sensitivity)**
```
행: HVG 수
열: PCA components
값: 최고 클러스터링 점수
→ 최적 파라미터 조합 식별
```

**7-4c: Feature Importance**
```
각 파라미터(HVG, PCA, Leiden resolution)의
메트릭 변화에 대한 영향도 계산
```

**7-4d: Box Plot (Noise Robustness)**
```
x축: 전처리 조합
y축: 메트릭 점수
그룹: Standard vs High-Noise 시뮬레이션
→ Robust한 파라미터 조합 식별
```

##### 7-5. Output Summary
```csv
preprocessing_combo, hvg, pca_n, leiden_res, NMI, ARI, FMI, VI
combo_1, 100, 20, 0.5, 0.85, 0.80, 0.82, 0.35
combo_2, 100, 20, 0.8, 0.88, 0.83, 0.85, 0.30
...
```

#### 핵심 코드
```python
# 시뮬레이션 함수 (pseudo)
def simulate_xenium(adata, n_reads_per_cell, dropout_rate, misseg_rate):
    """Simulate Xenium from scRNAseq"""
    # Random sampling with Poisson noise
    X_simulated = poisson(X * n_reads_per_cell / X.sum(axis=1))

    # Dropout (zero inflation)
    dropout_mask = np.random.random(X_simulated.shape) < dropout_rate
    X_simulated[dropout_mask] = 0

    # Misassignment (random label swap)
    misseg_idx = np.random.choice(n_cells, size=int(n_cells*misseg_rate), replace=False)
    new_labels = np.random.permutation(adata.obs['cell_type'].values)
    adata.obs.loc[misseg_idx, 'cell_type'] = new_labels[misseg_idx]

    return adata

# 그리드 검색
results = []
for hvg_n in [50, 100, 200]:
    for pca_n in [10, 20, 30]:
        adata_proc = preprocess(simulated, hvg=hvg_n, pca=pca_n)
        nmi = normalized_mutual_info_score(
            simulated.obs['cell_type'],
            adata_proc.obs['leiden']
        )
        results.append({'hvg': hvg_n, 'pca': pca_n, 'NMI': nmi})
```

#### 시각화
- Preprocessing grid 성능 비교 바플롯
- 파라미터별 메트릭 민감도 히트맵
- 파라미터 중요도 바플롯
- Standard vs High-Noise 비교 박스플롯

---

## 🔄 Notebook vs Pipeline 비교

### 구조적 차이

| 측면 | Notebook | Pipeline |
|------|----------|----------|
| **포맷** | Jupyter (.ipynb) | Python 모듈 (.py) |
| **실행 방식** | 셀 단위 interactive | `python pipeline_main.py` 전체 자동 실행 |
| **설정** | Notebook 내 하드코딩 | `config.yaml` 외부 설정 |
| **재현성** | 낮음 (실행 순서, 셀 누락 가능) | 높음 (결정적 실행) |
| **에러 처리** | Manual (try/except 부재) | Robust (로깅, 건너뛰기) |
| **모듈화** | 단일 파일 (Step 0: 0_formatting/, ...) | 8개 독립 모듈 (step0.py, ..., step7.py) |

### 코드 레벨 차이

#### 1. Config 관리
**Notebook:**

```python
# 변수 하드코딩
INPUT_PATH = "/path/to/data/Sample_outs"
OUTPUT_DIR = "../output"
MIN_GENES = 3
```

**Pipeline:**
```yaml
# config.yaml
input_path: "/path/to/data/Sample_outs"
output_dir: "./xenium-output"
filtering:
  min_genes: 3
  min_counts: 10
```

```python
# Python에서 로드
config = yaml.safe_load('config.yaml')
min_genes = config['filtering']['min_genes']
```

#### 2. 에러 처리 & 스킵 로직
**Notebook:**
```python
# 에러 처리 없음, 항상 재실행
adata = sc.read_h5ad(...)
```

**Pipeline:**
```python
step0_output = os.path.join(step0_dir, f"{sample_tag}.h5ad")
if os.path.exists(step0_output):
    print(f"Skipping Step 0: Output found")
else:
    try:
        adata = step0.run_step0(step0_config)
    except Exception as e:
        logger.error(f"Step 0 failed: {e}")
        return
```

#### 3. 임포트 방식
**Notebook:**
```python
# 모든 라이브러리 notebook 상단에서 임포트
import scanpy as sc
import pandas as pd
# ... 수십 줄
```

**Pipeline:**
```python
# 각 step 모듈에서 필요한 것만 임포트
# step0.py: matplotlib, scanpy, shutil, ...
# step1.py: scanpy, squidpy, scipy.stats, ...
```

#### 4. 함수 조직
**Notebook:**
```python
# 상위 레벨 함수 (예: format_xenium_adata_mid_2023)
def format_xenium_adata_mid_2023(path, tag, output_path):
    # 수백 줄 코드
    ...
    return adata
```

**Pipeline:**
```python
# run_step0(config) 진입점
def run_step0(config):
    """Step 0: Formatting"""
    path = config['input_path']
    output_path = config['output_dir']
    adata = format_xenium_adata_mid_2023(path, ..., output_path)
    return adata
```

#### 5. 데이터 흐름
**Notebook (Step 0-1):**
```
Step 0 Notebook
    ↓ (수동으로 출력 저장)
Step 1 Notebook
    ↓ (Step 0 아웃풋 경로 수동 입력)
...
```

**Pipeline:**
```
pipeline_main.py
    → step0_dir = "xenium-output/Sample_outs/step0_formatting"
    → step0_output = step0_dir + "/{sample_tag}.h5ad"
    ↓
    → step1_config['previous_step_adata_path'] = step0_output
    → step1_dir = "xenium-output/Sample_outs/step1_exploration"
    ↓ (자동 체이닝)
```

### 기능적 추가사항

Pipeline에서만 있는 기능:

1. **자동 샘플 이름 추출**
   ```python
   sample_name = os.path.basename(input_path.rstrip(os.sep))
   sample_output_dir = os.path.join(base_output_root, sample_name)
   ```

2. **유연한 데이터 경로 해석**
   - DAPI 이미지 다중 포맷 지원
   - Transcript 사이드카 자동 검출
   - Domain polygons 자동 링크

3. **조건부 Step 4 (Comparison)**
   - Step 3 실행 여부에 따라 comparison 범위 결정
   - Fallback: Step 3 실패 시 Step 1 데이터로 비교

4. **통합 로깅**
   ```python
   logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
   ```

5. **진행 상황 추적**
   ```
   Step 0: Formatting        [완료]
   Step 1: Dataset Exploration [실행 중]
   ...
   ```

---

## 🛠️ 기술 스택 & 핵심 라이브러리

| 라이브러리 | 용도 | 사용 단계 |
|-----------|------|---------|
| **scanpy** | AnnData I/O, QC, clustering | 모든 단계 |
| **squidpy** | Spatial analysis, neighbors | Step 1 |
| **Points2Regions** | Spatial domain clustering | Step 2 |
| **Cellpose** | Nuclei segmentation | Step 3 |
| **scikit-image** | Image processing (regionprops, EDT) | Step 2, 3 |
| **scipy.spatial** | cKDTree, ConvexHull | Step 5 |
| **pandas** | Data manipulation | 모든 단계 |
| **matplotlib/seaborn** | Visualization | 모든 단계 |
| **cellxgene_census** | scRNAseq reference download | Step 7 |
| **sklearn.metrics** | Clustering quality (NMI, ARI, FMI) | Step 7 |

---

## 📊 데이터 형식

### 입력 형식
- **Xenium machine output:**
  - `cell_feature_matrix.tar.gz` (matrix.mtx, barcodes.tsv, features.tsv)
  - `cells.csv` (cell metadata)
  - `transcripts.csv` (x_location, y_location, feature_name)
  - `morphology_focus.ome.tif` (DAPI)
  - `gene_panel.json` (Ensembl mapping)

### 중간 포맷
- **AnnData (.h5ad)**
  ```
  adata.X → cell×gene count matrix (sparse)
  adata.obs → cell metadata (cell_id, n_reads, n_genes, ...)
  adata.var → gene metadata (gene_name, gene_id, ...)
  adata.obsm → embeddings (UMAP, PCA)
  adata.obsl → spatial coordinates (x_coord, y_coord)
  adata.uns → miscellaneous (hvg list, clustering resolution, ...)
  ```

- **Transcript CSV**
  ```
  cell_id, x_location, y_location, feature_name, ...
  ```

- **Domain GeoJSON** (Step 2)
  ```json
  {
    "features": [
      {
        "properties": {"domain_id": 0},
        "geometry": {"type": "Polygon", "coordinates": [...]}
      }
    ]
  }
  ```

### 출력 형식
- **Benchmark Results (CSV)**
  ```
  method, n_cells, median_reads, assigned_prop, clustering_metric
  ```

- **Metrics Summary (CSV)**
  ```
  gene, method, mean_distance_um, std_distance_um, n_transcripts
  ```

---

## ⚙️ 설정 및 실행

### config.yaml 구조
```yaml
# 입출력
input_path: "/path/to/xenium/data"
output_dir: "./xenium-output"
sample_tag: "hbreast"

# QC 필터
filtering:
  min_genes: 3
  min_counts: 10

# Step 별 파라미터
step1:
  leiden_resolutions: [0.5, 1.0, 1.5]

step2:
  n_regions_auto: true

step3:
  cellpose_diameter: 30

step5:
  distance_threshold: 50

step6:
  run_baysor: true

step7:
  n_hvg: [50, 100, 200]
  n_pca: [10, 20, 30]
```

### 실행 방법
```bash
# 기본 설정으로 실행
python pipeline_main.py

# 커스텀 config 지정
python pipeline_main.py --config my_config.yaml
```

---

## 📈 결과 해석 가이드

### Step 0: Formatting
✅ **확인:**
- AnnData shape: (n_cells, n_genes)
- DAPI 이미지 로드됨
- Control probes 제거됨

❌ **문제 해결:**
- `cell_feature_matrix` 압축 해제 실패 → 파일 복사 확인
- Gene panel 매핑 실패 → JSON 형식 확인

### Step 1: Dataset Exploration
✅ **확인:**
- PCA/UMAP 임베딩 생성
- Leiden clustering 성공 (resolutions별 하나 이상 클러스터)
- Marker genes 검출됨

📊 **해석:**
- Distance ECDF: 가파를수록 집중된 유전자 분포
- KS test p-values: 낮을수록 유전자 간 공간 분리도 높음

### Step 2: Points2Regions
✅ **확인:**
- Domain polygons GeoJSON 생성
- Centroid distance vs Boundary distance 메트릭 계산

📊 **해석:**
- Domain 수: 조직 구조 복잡도 반영
- Boundary distance가 높을수록: 세포 경계 외 신호 많음

### Step 3: Resegmentation
✅ **확인:**
- Cellpose 마스크 생성 (TIF)
- 재할당된 cell×gene matrix 생성

📊 **해석:**
- Cell count 급증/감소: 분할 민감도 조정 필요
- Domain assignment: Step 2 도메인과 Cellpose 일관성

### Step 4: Techniques Comparison
📊 **핵심 지표:**
- **Efficiency:** Step 3 > Step 0 (재분할로 더 많은 transcript 할당)
- **Specificity:** Negative marker purity 높을수록 좋음
- **Diffusion:** Centroid distance 짧을수록 좋음

### Step 5: Optimal Expansion
📊 **결과:**
- Turnover distance > Nuclei size → Positive expansion 가능
- Optimal expansion = Turnover - Nuclei size
- 예) Turnover 25µm, Nuclei 15µm → Optimal ≈ 10µm

### Step 6: Benchmark
📊 **비교:**
| 방법 | n_cells | median_reads | 장점 | 단점 |
|------|---------|--------------|------|------|
| Nuclei | 고 | 낮음 | 해석 용이 | 신호 손실 |
| Cellpose | 중 | 중간 | 자동화 | 분할 오류 |
| Expansion | 중 | 중간 | 합리적 | 매개변수 의존 |
| Baysor | 중 | 높음 | 신호 포획 | 계산 비용 |

### Step 7: Simulation
📊 **최적 파라미터:**
```
HVG: 100 (충분한 feature 수)
PCA: 20 (차원 축소, 계산 효율)
Leiden resolution: 0.5-0.8 (안정적 클러스터링)
```

---

## 🐛 트러블슈팅

### Step 0 에러
```
[ERROR] 'cell_feature_matrix' folder and tar.gz not found
```
→ Xenium 원본 폴더 구조 확인

### Step 3 GPU 문제
```
CUDA out of memory
```
→ Cellpose diameter 감소 또는 CPU 모드 사용

### Step 6 Baysor 실패
```
Baysor not found in PATH
```
→ `pip install baysor` 또는 Docker 이용

### Step 7 메모리 부족
```
Memory error during grid search
```
→ HVG 수 감소 또는 시뮬레이션 샘플 크기 축소

---

## 📚 참고 문헌

- **Scanpy:** Wolf, F. A., Angerer, P., & Theis, F. J. (2018). *bioRxiv*
- **Squidpy:** Palla, G., et al. (2022). *Nature Methods*
- **Cellpose:** Stringer, C., et al. (2021). *Nature Methods*
- **Points2Regions:** [https://github.com/YichaoOU/points2regions](GitHub)

---

## 📝 체크리스트

### 파이프라인 실행 전
- [ ] Xenium raw data 준비
- [ ] config.yaml 수정 (input_path, output_dir)
- [ ] 필수 라이브러리 설치 (`pip install -r requirements.txt`)
- [ ] GPU 사용 가능 확인 (Step 3 Cellpose)

### 각 단계별 확인
- [ ] Step 0: AnnData 생성 확인
- [ ] Step 1: UMAP/clustering 시각화 확인
- [ ] Step 2: Domain polygons GeoJSON 생성
- [ ] Step 3: Cellpose 마스크 생성
- [ ] Step 4: 비교 메트릭 계산
- [ ] Step 5: Expansion transcript 저장
- [ ] Step 6: 벤치마크 결과 생성
- [ ] Step 7: 시뮬레이션 메트릭 저장

### 최종 확인
- [ ] 모든 output 디렉토리 생성됨
- [ ] 시각화 PNG/PDF 생성됨
- [ ] 메트릭 CSV 생성됨
- [ ] 로그 파일 검토

---

## 🔗 관련 파일

- **파이프라인 메인:** `pipeline/pipeline_main.py`
- **설정 파일:** `pipeline/config.yaml`
- **스텝 모듈들:** `pipeline/xenium_step*.py` (0-7)
- **유틸리티:** `xb/` 라이브러리
- **벤치마크 유틸:** `pipeline/benchmark_utils/`

---

**문서 버전:** v1.0
**마지막 업데이트:** 2026-02-10
**유지보수:** MASLD Xenium Pipeline Team
