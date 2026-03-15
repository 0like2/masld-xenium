# Step 4: Techniques Comparison & Validation — Metrics Detailed Reference

> **Paper**: "Optimizing Xenium In Situ data utility by quality assessment and best-practice analysis workflows"
> (Salas et al., *Nature Methods*, Volume 22, April 2025, pp. 813–823)
> **Pipeline file**: `pipeline/xenium_step4_techniques_comparison.py`
> **Original notebooks**: `notebooks/3_techniques_comparison/3_3` ~ `3_7`

---

## Overview

Step 4는 **원본(Original) 세그멘테이션 데이터**와 **Step 3에서 재세그멘테이션(Resegmented)된 데이터**를 비교하여 재세그멘테이션의 품질을 정량적으로 평가한다. 논문에서는 Xenium을 포함한 8개 SRT 플랫폼 간 크로스-플랫폼 비교에 사용된 4가지 핵심 메트릭을 정의하며, 파이프라인은 이를 Resegmented vs Original 비교에 적용한다.

논문 Fig. 2a에서 제시한 비교 워크플로우:
```
Image-based ST datasets → Nuclei Segmentation (Cellpose) → Common segmentation
                        → Regional annotation → Region-specific comparable datasets
                        → Comparison: Gene efficiency / Detection specificity /
                                      Read diffusion / *Genes profiled / *Cells profiled /
                                      *Total reads per cell
```

**4가지 핵심 메트릭 카테고리:**

| # | Metric Category | Paper Section | Notebook | Pipeline Function |
|---|----------------|--------------|----------|-------------------|
| 1 | **Efficiency** | Fig. 2b,c,g; Methods | `3_3` | `analyze_efficiency()` |
| 2 | **Specificity (NMP/NCP)** | Fig. 2d, 3e; Methods | `3_4` | `analyze_specificity()` |
| 3 | **Positivity** | Fig. 2b,e; Methods | `3_5` | `analyze_positivity()` |
| 4 | **Diffusion** | Fig. 2f; Methods | `3_6` | `analyze_diffusion()` |

---

## 1. Efficiency (Detection Efficiency)

### 1.1 논문에서의 정의

**"Detection efficiency"** 는 각 SRT 플랫폼이 개별 유전자의 전사체를 얼마나 효율적으로 검출하는지를 정량화한다 (논문 p.815–816, Fig. 2c).

> "We calculated the detection efficiency for individual genes for each technology by comparing read counts obtained for each gene with a reference region-matched scRNA-seq dataset."

논문은 두 가지 접근법을 사용:

1. **Per-gene expression ratio (SRT/scRNA-seq)**: scRNA-seq를 레퍼런스로 삼아 각 유전자의 검출 효율을 비교
2. **Cross-platform direct comparison**: 모든 플랫폼의 데이터를 같은 전처리로 처리 후 클러스터링하여 유전자 발현 수준을 직접 비교

#### 수학적 정의

**Per-gene Expression Ratio:**
```
For each gene g:
  1. positive cells = cells where count(g) > minreads  (minreads = 1)
  2. ST_mean(g) = mean expression of g among positive ST cells
  3. SC_median(g) = median expression of g among positive scRNA-seq cells
  4. Ratio(g) = ST_mean(g) / SC_median(g)
```

- **Ratio ≈ 1**: SRT 플랫폼 검출 효율이 scRNA-seq와 유사
- **Ratio > 1**: SRT에서 과검출 (over-detection)
- **Ratio < 1**: SRT에서 과소검출 (under-detection)

논문 결과: Xenium의 효율은 scRNA-seq (Chromium v2)의 **1.2–1.5배** (Fig. 2c), 뇌 영역에 따라 다름.

### 1.2 노트북 구현 (3_3)

노트북 `3_3_efficiency_between_methods.ipynb`에서의 핵심 로직:

```python
# 각 유전자에 대해 ST mean과 scRNA median 계산
minreads = 1
for each gene g in common_genes:
    st_expr = ST_cells[ST_cells[g] > minreads][g]
    sc_expr = SC_cells[SC_cells[g] > minreads][g]
    ratio(g) = mean(st_expr) / median(sc_expr)
```

노트북 특이사항:
- **3개 뇌 영역별 분석**: Cortex, Hippocampus, Thalamus (Fig. 2g의 Xenium vs Visium scatter)
- 유전자를 lowercase로 통일하여 매칭
- 최소 4개 이상의 데이터셋에서 공통인 유전자만 사용 (`df[0] > 2`)
- `xb.comparing.median_calculator()` 함수 사용

### 1.3 파이프라인 구현

`analyze_efficiency()` 함수 (`xenium_step4_techniques_comparison.py:237–551`)

**Sub-analysis 구성:**

| Sub-step | 설명 | 출력 파일 |
|----------|------|---------|
| **4-2a** | Transcripts/Genes per cell 히스토그램 | `efficiency_transcripts_per_cell_comparison.png`, `efficiency_genes_per_cell_comparison.png` |
| **4-2b** | Expression ratio (ST mean / scRNA median) | `efficiency_expression_ratio.csv`, `efficiency_expression_ratio.png` |
| **4-2c** | Region별 효율 분석 | `efficiency_region_breakdown.csv`, `efficiency_expression_ratio_by_region.csv` |
| 추가 | Reseg vs Original boxplot | `efficiency_reseg_vs_original_boxplot.png` |
| 추가 | ST vs scRNAseq scatter (log-log) | `efficiency_st_vs_sc_scatter.png` |

**주요 계산 흐름:**

```
1. QC metrics 계산 (scanpy.pp.calculate_qc_metrics)
   → total_counts, n_genes_by_counts
2. Resegmented vs Original 히스토그램 비교
3. scRNA-seq reference 사용 시:
   a. Common genes 식별
   b. Per-gene: ST mean (expressing cells, count > minreads)
   c. Per-gene: scRNA median (expressing cells, count > minreads)
   d. Ratio = ST_mean / SC_median
4. Region 정보 존재 시: region별 breakdown
```

**출력 통계 (CSV):**
```
Dataset | median_transcripts_per_cell | mean_transcripts_per_cell |
        | median_genes_per_cell       | mean_genes_per_cell
```

### 1.4 논문과 파이프라인 간 차이점

| 항목 | 논문/노트북 | 파이프라인 |
|------|-----------|----------|
| 비교 대상 | 6개 SRT 플랫폼 간 | Resegmented vs Original (동일 플랫폼 내) |
| minreads | 1 | config `efficiency_minreads` (기본값 1) |
| 정규화 | Raw counts 사용 | Raw counts 사용 (layers['raw'] 우선) |
| 영역 분리 | Cortex/Hippocampus/Thalamus 수동 | obs column 자동 감지 (`region_annotation` 등) |
| scRNA-seq | 필수 (Zhang et al. 2023) | 선택적 (`sc_reference_path`) |

---

## 2. Specificity (Negative Marker Purity, NMP/NCP)

### 2.1 논문에서의 정의

**Negative Co-expression Purity (NCP)** (논문에서는 NMP와 NCP를 혼용): scRNA-seq에서 공동발현하지 않는 유전자 쌍이 SRT에서도 공동발현하지 않는 비율을 측정한다 (논문 p.816–817, Fig. 2d).

> "We implemented a metric called negative co-expression purity (NCP), which quantifies the percentage of non-co-expressed genes in our reference single-cell dataset that do not appear to be coexpressed in each SRT dataset."

#### 수학적 정의 (논문 Methods)

**NMP (Negative Marker Purity) — reads 기반:**

Step 1: 정규화된 평균 발현 계산
```
x̄ₘ_g,c = mean raw expression of gene g in cell type c in modality m (sp or sc)

X̂ₘ_g,c = x̄ₘ_g,c / Σ_c' x̄ₘ_g,c'
```
즉, 각 유전자 g에 대해 cell type별 평균 발현을 모든 cell type 합으로 나눠 정규화.

Step 2: Negative marker 쌍 식별
```
negative_marker_mask = (mean_ct_sc_relative < minimum_exp)
   where minimum_exp = 0.005 (reads) or 0.05 (coexpression)
```
scRNA-seq에서 특정 cell type에서의 상대 발현이 threshold 미만인 (gene, celltype) 쌍을 "negative marker"로 정의.

Step 3: Negative marker 영역의 발현 집계
```
X̂ₛₚ_neg = Σ_{g∈G} Σ_{c∈C^neg_g} X̂ₛₚ_g,c / |P^neg|
X̂ₛc_neg = Σ_{g∈G} Σ_{c∈C^neg_g} X̂ₛc_g,c / |P^neg|
```
여기서 P^neg는 negative marker (gene, celltype) 쌍의 집합.

Step 4: NMP score 계산
```
NMP = { 1 - (X̂ₛₚ_neg - X̂ₛc_neg)   if X̂ₛₚ_neg > X̂ₛc_neg
      { 1                            otherwise
```

- **NMP = 1.0**: SRT에서 negative marker 발현이 scRNA-seq 이하 → 완벽한 specificity
- **NMP < 1.0**: SRT에서 negative marker가 더 많이 검출됨 → read leakage/diffusion 존재
- **NMP ≈ 0**: 심각한 specificity 문제

논문 결과: 모든 플랫폼 NMP > 0.8, HS-ISS와 Molecular Cartography가 가장 높음 (Fig. 2d).

### 2.2 논문의 두 가지 NMP 변형

`metrics.py`에 구현된 두 가지 변형:

#### (A) `negative_marker_purity_reads()` — reads 기반 (논문 주요 메트릭)
```
- Cell type별 mean expression 사용
- 상대 발현 = mean_celltype / mean_over_all_celltypes
- Negative marker threshold: relative expression < 0.005
- 정규화: cell type 합으로 나눔 (sum normalization)
- 결과: 0~1 범위, reads의 cell type 간 leakage 측정
```

#### (B) `negative_marker_purity_cells()` — cells 기반
```
- Cell type별 positive cell fraction 사용 (binary: expressed or not)
- Negative marker threshold: positive cell ratio < 0.005
- 정규화 없음 (이미 비율)
- 결과: 0~1 범위, 세포 수준의 specificity 측정
```

### 2.3 파이프라인 구현 — co-expression 기반 NMP

Step 4 파이프라인은 `metrics.py`의 함수가 아닌 **자체 내장된 co-expression 기반 NMP** (`_negative_marker_purity_coexpression()`)를 사용한다. 이는 노트북 3_4에서 사용하는 `xb.calculating.negative_marker_purity_coexpression()` 함수를 파이프라인 내부에 재구현한 것이다.

`_negative_marker_purity_coexpression()` (`xenium_step4_techniques_comparison.py:51–117`)

**계산 과정:**

```python
# 1. Gene name 소문자 통일 + 공통 유전자 필터링
# 2. Co-expression matrix 계산 (gene × gene)
for each gene col:
    positive_cells = cells where expression(col) > min_exp (default 0)
    coexpression[other_gene, col] = fraction of positive_cells expressing other_gene

# 3. Negative marker 식별
neg_marker_mask = (sc_coexpression < 0.05)

# 4. NMP 계산
lowvals_diff = sp_coexpression[neg_mask] - sc_coexpression[neg_mask]
lowvals_diff[lowvals_diff < 0] = 0   # 음수는 무시 (SP가 더 낮은 건 문제 없음)
NMP = 1 - mean(lowvals_diff)
```

**노트북 3_4와의 대응:**
- 노트북은 `ratio < 10` 필터 (효율이 극단적으로 높은 유전자 제외) → 파이프라인도 동일하게 `efficiency_expression_ratio.csv`에서 `ratio < 10` 필터 적용
- 노트북은 cortex 영역만 분석 → 파이프라인은 전체 데이터 사용

### 2.4 Specificity 분석 출력

`analyze_specificity()` (`xenium_step4_techniques_comparison.py:553–689`)

| Sub-step | 설명 | 출력 파일 |
|----------|------|---------|
| **4-3a** | NMP score (scalar) | `specificity_nmp_score.txt` |
| **4-3a-1** | Per-gene purity CSV | `specificity_nmp_per_gene_{method}.csv` |
| **4-3a-2** | NMP per-gene boxplot | `specificity_nmp_per_gene_boxplot.png` |
| **4-3a-3** | Efficiency vs Specificity scatter | `specificity_vs_efficiency_scatter_{method}.png` |
| **4-3b** | Gene-gene correlation heatmap | `specificity_gene_correlation_{method}.png` |

**Efficiency vs Specificity scatter** (논문 노트북 3_4 마지막 셀):
- X축: Expression ratio (efficiency)
- Y축: Purity score (specificity)
- 목적: 효율과 특이성 간 trade-off 시각화

### 2.5 Gene-Gene Correlation (scRNA-seq 없을 때의 대체 메트릭)

scRNA-seq reference가 없는 경우, 파이프라인은 **gene-gene correlation heatmap**을 대체 specificity 지표로 사용:

```python
top_genes = top 50 expressed genes
corr_matrix = np.corrcoef(X[:, top_genes], rowvar=False)
```

- 높은 off-diagonal correlation: 잠재적 read diffusion (인접 세포 간 전사체 누출)
- Resegmented vs Original 비교 시 correlation 감소 → specificity 개선 의미

---

## 3. Positivity (Gene Detection Rate)

### 3.1 논문에서의 정의

**Positivity**는 각 유전자가 검출되는 세포의 비율(fraction)을 의미한다. 논문에서는 이를 플랫폼 간 세포당 유전자 수(genes per cell)와 전사체 수(transcripts per cell) 비교의 맥락에서 사용한다 (Fig. 2b,e).

> "To further validate these observations, we independently clustered the cells from each dataset using standardized analysis pipelines, identified shared populations and compared gene expression levels across technologies."

### 3.2 노트북 구현 (3_5)

`3_5_Computing_positivity_after_preprocessing_for_all_ST_techs.ipynb`

**핵심 전처리 파이프라인 (모든 플랫폼 동일):**
```python
def main_preprocessing(adata):
    adata.layers['raw'] = adata.X.copy()
    sc.pp.filter_cells(adata, min_counts=10)
    sc.pp.filter_cells(adata, min_genes=3)
    sc.pp.normalize_total(adata, target_sum=None)
    sc.pp.log1p(adata)
    sc.pp.neighbors(adata, n_neighbors=8, n_pcs=0)
    sc.tl.leiden(adata, resolution=2.2, key_added='leiden_2_2')
    sc.tl.leiden(adata, resolution=1.4, key_added='leiden_1_4')
    sc.tl.leiden(adata, resolution=0.6, key_added='leiden_0_6')
    sc.tl.umap(adata, min_dist=0.1)
    return adata
```

**분석 내용:**
1. **Per-gene positivity**: `n_cells_by_counts / n_obs` = 해당 유전자가 발현된 세포 비율
2. **Cluster별 발현 비교**: Leiden 클러스터링 후, 특정 마커 유전자의 최고 발현 클러스터를 자동 식별
3. **Violin plot**: 각 유전자에 대해 최적 클러스터에서의 발현 분포를 플랫폼 간 비교

### 3.3 파이프라인 구현

`analyze_positivity()` (`xenium_step4_techniques_comparison.py:692–821`)

**전처리 파라미터 (config.yaml 대응):**

| 파라미터 | 노트북 기본값 | Config 키 | 설명 |
|---------|-----------|----------|------|
| `n_neighbors` | 8 | `comparison.positivity.n_neighbors` | KNN 그래프 이웃 수 |
| `n_pcs` | 0 | `comparison.positivity.n_pcs` | PCA 차원 (0 = 전체 발현 사용) |
| `leiden_resolution` | 2.2 (primary) | `comparison.positivity.leiden_resolution` | 주요 Leiden 해상도 |
| `leiden_resolutions` | [2.2, 1.4, 0.6] | `comparison.positivity.leiden_resolutions` | 다중 해상도 |
| `umap_min_dist` | 0.1 | `comparison.positivity.umap_min_dist` | UMAP 최소 거리 |
| `min_counts` | 10 | `comparison.positivity.min_counts` | 최소 전사체/세포 |
| `min_genes` | 3 | `comparison.positivity.min_genes` | 최소 유전자/세포 |

**n_pcs = 0의 의미**: PCA를 수행하지 않고 전체 유전자 발현 matrix를 직접 KNN 그래프 구축에 사용. Xenium 패널이 ~248–541 유전자로 상대적으로 작기 때문에 차원 축소 없이도 효과적.

**Sub-analysis 구성:**

| Sub-step | 설명 | 출력 파일 |
|----------|------|---------|
| **4-4a** | Positivity 분포 히스토그램 | `{sample}_positivity_dist_comparison.png` |
| **4-4b** | 전처리 + Leiden 클러스터링 (3개 해상도) | — |
| **4-4c** | 클러스터별 Violin plot (top 10 양성 유전자) | `{sample}_positivity_violin_clusters.png` |
| **4-4d** | UMAP per top gene (top 4) | `{sample}_positivity_umap_top_genes.png` |
| 추가 | Optimal cluster per gene table | `{sample}_positivity_optimal_cluster.csv` |

**Positivity 계산:**
```python
# Positivity = fraction of cells expressing each gene
adata.var['positivity'] = adata.var['n_cells_by_counts'] / adata.n_obs
```

- **positivity ≈ 1.0**: 거의 모든 세포에서 발현되는 유전자 (housekeeping)
- **positivity ≈ 0.0**: 극소수 세포에서만 발현 (cell type-specific marker)

**Optimal cluster 식별:**
```python
for each gene:
    optimal_cluster = argmax(mean_expression_per_cluster)
```
각 유전자가 가장 높게 발현되는 Leiden 클러스터를 자동 식별하여, 세포 타입 특이적 마커의 클러스터 배치를 검증.

---

## 4. Diffusion (Read Diffusion / Distance to Centroid)

### 4.1 논문에서의 정의

**Read diffusion**은 각 전사체(read/transcript)가 할당된 세포의 중심(centroid)으로부터 얼마나 떨어져 있는지를 측정한다 (논문 p.817, Fig. 2f).

> "With the aim of exploring the diffusion of the different technologies, reads identified outside the nuclei were assigned to their closest cell."

Diffusion은 세그멘테이션 정확도와 전사체 확산(RNA diffusion) 수준을 반영하는 지표:
- **짧은 거리** → 전사체가 세포 내부에 잘 국한됨 → 좋은 세그멘테이션
- **긴 거리** → 전사체가 세포 경계 밖으로 확산 or 잘못된 세포에 할당

### 4.2 수학적 정의

**Transcript-to-Centroid Distance:**
```
For each transcript t assigned to cell c:
    d(t) = sqrt((x_t - cx_c)² + (y_t - cy_c)²)

    where:
      (x_t, y_t) = transcript position
      (cx_c, cy_c) = centroid of assigned cell c
```

**단위 변환 (pixel → µm):**
```
distance_um = distance_px / conversion_factor

Technology    | Factor (px/µm)
Xenium        | 4.70588
CosMx         | 8.3333
Vizgen        | 9.20586
MERFISH       | 9.28
HybrISS       | 3.11
Resolved Bio  | 7.24
```

논문에서의 결과 (Fig. 2f):
- MERFISH, ISS 계열: 전사체가 세포 중심 근처에 집중
- CosMx, MC, MERSCOPE: 전사체가 더 넓게 분포 (더 큰 세포 or 더 많은 diffusion)
- Xenium: MERFISH와 유사한 짧은 distance

### 4.3 노트북 구현 (3_6)

`3_6_Diffussion_on_resegmented_data.ipynb`

핵심 계산:
```python
# 1. 거리 계산 (pixel 단위)
rd2['distance_to_centroid'] = sqrt(
    (rd2['x'] - rd2['closest_cell_y'])^2 +
    (rd2['y'] - rd2['closest_cell_x'])^2
)

# 2. 플랫폼별 pixel→µm 변환
rd2.loc[method=='Xenium', 'distance_to_centroid'] /= 4.70588

# 3. Complementary CDF plot
sns.ecdfplot(data=rd2, x='distance_to_centroid', hue='method', complementary=True)
```

노트북의 추가 분석:
- 10% 랜덤 샘플링 (`sample(list(transcripts.index), int(transcripts.shape[0]*0.1))`)
- **Per-gene ECDF subplots**: 상위 9개 유전자별 거리 분포
- **Gene × Method heatmap**: 유전자별 평균 거리를 플랫폼 간 비교
- **Assigned reads barplot**: 세포에 할당된 전사체 비율

### 4.4 파이프라인 구현

`analyze_diffusion()` (`xenium_step4_techniques_comparison.py:824–1103`)

**Sub-analysis 구성:**

| Sub-step | 설명 | 출력 파일 |
|----------|------|---------|
| **4-5a** | Transcript-to-centroid 거리 계산 (px→µm) | `reseg_diffusion_stats.csv`, `original_diffusion_stats.csv` |
| **4-5a+** | Per-gene 거리 요약 | `diffusion_per_gene_summary.csv` |
| **4-5b** | Complementary CDF plot | `diffusion_complementary_cdf_comparison.png` |
| **4-5c** | Per-gene ECDF subplots (top 9 유전자) | `diffusion_per_gene_ecdf.png` |
| **4-5d** | Gene × Method mean distance heatmap | `diffusion_gene_method_heatmap.png`, `diffusion_gene_method_mean_distances.csv` |
| **4-5e** | Assigned reads stacked barplot | `diffusion_assigned_reads_barplot.png` |

**좌표 단위 자동 감지 로직:**

파이프라인은 centroid과 transcript 좌표의 단위가 다를 수 있는 상황을 자동 처리:
```python
cx_range = centroid_x.max() - centroid_x.min()
tx_range = transcript_x.max() - transcript_x.min()
ratio = cx_range / tx_range

if ratio > 2.0:
    # Centroids in pixels, transcripts in µm → convert centroids
    cx_vals /= conversion_factor
elif ratio < 0.5:
    # Centroids in µm, transcripts in pixels → convert transcripts
    x_vals /= conversion_factor
```

**Cell ID 매칭 전략:**

transcript의 cell_id와 adata.obs의 인덱스가 다른 형식일 수 있으므로 (예: numeric index vs Xenium barcode string), 파이프라인은 자동으로 적절한 key를 선택:
```python
if transcript cell_ids are non-numeric strings and adata.obs has 'cell_id' column:
    use adata.obs['cell_id'] as centroid key
else:
    use adata.obs.index as centroid key
```

### 4.5 Complementary CDF (1-CDF) 해석

Complementary CDF plot은 "해당 거리 이상에 위치한 전사체의 비율"을 보여준다:

```
Y축 값 = P(Distance > x) = 1 - CDF(x)
```

- **커브가 빠르게 떨어짐** → 대부분의 전사체가 centroid 가까이 → 좋은 세그멘테이션
- **커브가 천천히 떨어짐** → 전사체가 넓게 분포 → 더 많은 diffusion
- **Resegmented 커브 < Original 커브** → 재세그멘테이션이 전사체를 더 정확한 세포에 할당

### 4.6 Assigned Reads Proportion 해석

```
Proportion = (transcripts assigned to cells) / (total transcripts)
```

논문 (p.814): Xenium 기본 세그멘테이션에서 평균 76.8%의 reads가 세포에 할당됨. 재세그멘테이션 후 이 비율의 변화는 세그멘테이션 전략의 보수적/적극적 수준을 반영.

---

## 5. 메트릭 간 관계 및 종합 해석

### 5.1 Efficiency vs Specificity Trade-off

논문에서 강조하는 핵심 관계:

```
                    High Specificity (NMP → 1.0)
                           ↑
                           |
     Conservative seg.     |     ★ Ideal
     (fewer reads/cell,    |     (high capture,
      high purity)         |      high purity)
                           |
    ──────────────────────────────────────→ High Efficiency (Ratio → 1.0)
                           |
     Poor performance      |     Aggressive seg.
     (low capture,         |     (more reads/cell,
      low purity)          |      possible contamination)
                           |
```

- **재세그멘테이션의 목표**: Efficiency를 높이면서 Specificity를 유지
- 노트북 3_4에서 ratio < 10 필터를 적용하는 이유: 극단적으로 높은 efficiency ratio를 가진 유전자는 NMP 계산을 왜곡

### 5.2 Diffusion과 Expansion의 관계

논문 Fig. 3a,b에서 보여주는 핵심 통찰:

```
Distance to centroid ↔ Cell type specificity
   ≤ 5.06 µm (nuclei radius): Nuclear signature 우세
   5.06–10.71 µm: 혼합 영역 (nuclear + cytoplasmic)
   > 10.71 µm: Background/domain signature 우세
```

Step 5 (Optimal Expansion)에서 이 관계를 이용하여 최적 확장 거리를 결정.

### 5.3 전체 품질 지표 요약

| Metric | 좋은 결과 | 나쁜 결과 | 논문 기준값 (Xenium) |
|--------|---------|---------|-----------------|
| **Efficiency ratio** | ≈ 1.0 | ≪ 1 or ≫ 1 | 1.2–1.5 (vs scRNA-seq) |
| **NMP score** | → 1.0 | → 0 | > 0.8 |
| **Positivity** | 마커 유전자가 해당 클러스터에서 높은 발현 | 고르게 분포 | 패널 의존적 |
| **Mean distance to centroid** | ≤ 5 µm | > 15 µm | ~10.71 µm (with expansion) |
| **Assigned reads %** | 70–90% | < 50% | 76.8% (default seg.) |

---

## 6. Config 설정 참조

```yaml
# pipeline/config.yaml — Step 4 관련 설정

comparison:
  run_efficiency: true          # Efficiency 분석 활성화
  run_specificity: true         # Specificity (NMP) 분석 활성화
  run_positivity: true          # Positivity 분석 활성화
  run_diffusion: true           # Diffusion 분석 활성화
  technology: "xenium"          # pixel→µm 변환 플랫폼

  sc_reference_path: null       # scRNAseq .h5ad 경로 (NMP + Efficiency에 필요)
  efficiency_minreads: 1        # 최소 read 수 (expressing cell 판정)

  positivity:
    n_neighbors: 8              # KNN 이웃 수
    n_pcs: 0                    # PCA 차원 (0 = no PCA)
    leiden_resolution: 2.2      # Primary Leiden resolution
    leiden_resolutions: [2.2, 1.4, 0.6]  # Multi-resolution
    umap_min_dist: 0.1          # UMAP 최소 거리
    min_counts: 10              # 세포 필터링: 최소 transcript 수
    min_genes: 3                # 세포 필터링: 최소 유전자 수
```

---

## 7. 출력 디렉토리 구조

```
xenium-output/
└── figures/
    └── 4_techniques_comparison/
        ├── spatial_roi_map.png
        │
        ├── # Efficiency (4-2)
        ├── efficiency_transcripts_per_cell_comparison.png
        ├── efficiency_genes_per_cell_comparison.png
        ├── efficiency_metrics_comparison.csv
        ├── efficiency_expression_ratio.csv
        ├── efficiency_expression_ratio.png
        ├── efficiency_st_vs_sc_scatter.png
        ├── efficiency_region_breakdown.csv
        ├── efficiency_region_breakdown.png
        ├── efficiency_expression_ratio_by_region.csv
        ├── efficiency_expression_ratio_by_region.png
        ├── efficiency_ratio_boxplot_log2.png
        ├── efficiency_reseg_vs_original_boxplot.png
        │
        ├── # Specificity (4-3)
        ├── specificity_nmp_score.txt
        ├── specificity_nmp_per_gene_resegmented.csv
        ├── specificity_nmp_per_gene_original.csv
        ├── specificity_nmp_per_gene_all.csv
        ├── specificity_nmp_per_gene_boxplot.png
        ├── specificity_vs_efficiency_scatter_resegmented.png
        ├── specificity_gene_correlation_resegmented.png
        ├── specificity_gene_correlation_original.png
        │
        ├── # Positivity (4-4)
        ├── {sample}_resegmented_gene_positivity.csv
        ├── {sample}_original_gene_positivity.csv
        ├── {sample}_positivity_dist_comparison.png
        ├── {sample}_positivity_violin_clusters.png
        ├── {sample}_positivity_umap_top_genes.png
        ├── {sample}_positivity_optimal_cluster.csv
        │
        ├── # Diffusion (4-5)
        ├── reseg_diffusion_stats.csv
        ├── original_diffusion_stats.csv
        ├── diffusion_per_gene_summary.csv
        ├── diffusion_complementary_cdf_comparison.png
        ├── diffusion_per_gene_ecdf.png
        ├── diffusion_gene_method_heatmap.png
        ├── diffusion_gene_method_mean_distances.csv
        └── diffusion_assigned_reads_barplot.png
```

---

## 8. 핵심 수식 요약

### Efficiency Ratio
$$\text{Ratio}(g) = \frac{\text{mean}_{c \in \text{ST}}(x_{g,c} \mid x_{g,c} > 1)}{\text{median}_{c \in \text{SC}}(x_{g,c} \mid x_{g,c} > 1)}$$

### Negative Marker Purity (NMP)
$$\hat{X}^{(m)}_{g,c} = \frac{\bar{x}^{(m)}_{g,c}}{\sum_{c' \in C} \bar{x}^{(m)}_{g,c'}}$$

$$\text{NMP} = \begin{cases} 1 - (\hat{X}^{(sp)}_{neg} - \hat{X}^{(sc)}_{neg}) & \text{if } \hat{X}^{(sp)}_{neg} > \hat{X}^{(sc)}_{neg} \\ 1 & \text{otherwise} \end{cases}$$

### Positivity
$$\text{Positivity}(g) = \frac{|\{c : x_{g,c} > 0\}|}{N_{\text{cells}}}$$

### Diffusion Distance
$$d(t, c) = \sqrt{(x_t - \text{cx}_c)^2 + (y_t - \text{cy}_c)^2} \div \text{px\_per\_µm}$$

---

*Document generated: 2026-02-20*
*Based on: Paper (Salas et al., Nat. Methods 2025), Notebooks 3_3–3_6, Pipeline Step 4 implementation*
