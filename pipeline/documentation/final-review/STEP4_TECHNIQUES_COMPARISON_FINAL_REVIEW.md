# Step 4: Techniques Comparison & Validation - Final Review

## 1. Overview

### 1.1 What This Step Does
Step 4는 Step 3의 **재분할(resegmented)** 데이터와 Step 0/1의 **원본(nuclei)** 데이터를 4가지 메트릭으로 정량적으로 비교한다:
1. **Efficiency** (효율성): 세포당 transcript/gene 수 비교, scRNAseq 대비 발현 효율
2. **Specificity** (특이성): Negative Marker Purity (NMP) 기반 비특이적 co-expression 측정
3. **Positivity** (양성률): 유전자별 양성 세포 비율 분석 및 클러스터링
4. **Diffusion** (확산): Transcript-to-centroid 거리 분포 비교

### 1.2 Paper Context
논문의 "Xenium detection efficiency matches ISH" 섹션 (p.816-817)과 Fig. 2에서 이 비교 분석을 수행한다:
- "We calculated the detection efficiency for individual genes for each technology" → Efficiency
- "We implemented a metric called negative co-expression purity (NCP)" → Specificity (NMP/NCP)
- "We aimed to explore the diffusion of the different technologies" → Diffusion
- "To systematically identify subcellular mRNA clusters" → Positivity (preprocessing)

### 1.3 논문 Figure 매칭
| 논문 Figure | 설명 | Step 4 구현 |
|------------|------|------------|
| **Fig. 2a** | Comparison workflow diagram | `run_step4()` 전체 workflow |
| **Fig. 2b** | Transcripts/cell, Genes/cell boxplot (플랫폼 간) | `analyze_efficiency()` |
| **Fig. 2c** | SRT/SC gene efficiency ratio boxplot | `analyze_efficiency()` → expression ratio |
| **Fig. 2d** | Gene specificity (NCP) boxplot | `analyze_specificity()` → NMP score |
| **Fig. 2e** | Violin: transcript counts per gene per platform | `analyze_efficiency()` |
| **Fig. 2f** | Cumulative distance to centroid (플랫폼 간) | `analyze_diffusion()` → ECDF |
| **Fig. 2g** | Xenium vs Visium scatter | 미구현 (visium_path=null) |
| **Extended Data Fig. 4b** | Assigned reads proportion | `analyze_diffusion()` → assigned_reads |
| **Extended Data Fig. 4c** | SRT/SC efficiency 비교 (hippocampus, thalamus) | `analyze_efficiency()` |
| **Extended Data Fig. 4d** | NCP low-efficiency genes | `analyze_specificity()` |
| **Extended Data Fig. 4f** | Density: distance to centroid per gene | `analyze_diffusion()` |
| **Extended Data Fig. 4g** | Resegmented datasets ROI | Step 3 masks overlay |

---

## 2. Pipeline Process Flow

```
Step 3 Output (재분할)                    Step 0/1 Output (원본)
├── step3_resegmented.h5ad               ├── {tag}.h5ad
├── step3_transcripts_resegmented.csv    ├── transcripts.parquet
                    ↓                              ↓
              ┌─────────────────────────────────────────┐
              │        [Step 4: Techniques Comparison]     │
              └─────────────────────────────────────────┘
                                ↓
    ┌───── 4-1. Load Data ──────────────────────────────────┐
    │  원본 adata + 재분할 adata + scRNAseq reference (선택)  │
    └────────────────────────────────────────────────────────┘
         ↓
    ┌───── 4-2. Efficiency Analysis ────────────────────────┐
    │  Transcripts/cell histogram                            │
    │  Genes/cell histogram                                  │
    │  ST/SC expression ratio                                │
    │  Region-specific efficiency                            │
    └────────────────────────────────────────────────────────┘
         ↓
    ┌───── 4-3. Specificity (NMP) ──────────────────────────┐
    │  Negative marker pairs 식별                            │
    │  Co-expression 계산                                    │
    │  NMP score (0-1)                                       │
    └────────────────────────────────────────────────────────┘
         ↓
    ┌───── 4-4. Positivity Analysis ────────────────────────┐
    │  Gene positivity (양성 세포 비율)                       │
    │  Preprocessing + Leiden clustering                     │
    │  Violin plot per cluster                               │
    └────────────────────────────────────────────────────────┘
         ↓
    ┌───── 4-5. Diffusion Analysis ─────────────────────────┐
    │  Transcript-centroid 거리 (pixel→µm)                   │
    │  Complementary CDF                                     │
    │  Per-gene ECDF                                         │
    │  Gene×Method distance heatmap                          │
    └────────────────────────────────────────────────────────┘
         ↓
├── figures/4_techniques_comparison/ 출력 디렉토리
│   ├── efficiency_*.png/csv
│   ├── specificity_*.png/csv
│   ├── positivity_*.png/csv
│   ├── positivity_key_markers*.png/csv  (NEW: brain marker comparison)
│   └── diffusion_*.png/csv
```

---

## 3. Sub-step 상세 분석

### 3.1 Efficiency Analysis (`analyze_efficiency()`, lines 246-618)

**목적**: 재분할 vs 원본의 transcript 포착 효율 비교

**3.1.1 Transcripts/cell 비교**
```python
# 두 방법의 세포당 transcript 수 히스토그램 비교
reseg_counts = adata_reseg.obs['n_counts']
orig_counts = adata_orig.obs['n_counts']
# → 재분할이 더 많은 transcript를 포착하는지 확인
```

**3.1.2 Genes/cell 비교**
```python
# 세포당 발현 유전자 수 비교
reseg_genes = adata_reseg.obs['n_genes']
orig_genes = adata_orig.obs['n_genes']
```

**3.1.3 ST/scRNAseq Expression Ratio** (scRNAseq reference 필요)
```python
# Raw counts 기반 유전자별 발현량 비교 (CPM 아님 — notebook 방식)
# Per-gene: ST mean over expressing cells / scRNA median over expressing cells
# Paper: "cells with more than one read (count > 1)" — strict > minreads 사용
st_expr = st_col[st_col > minreads]  # minreads = config.efficiency_minreads (default 1)
sc_expr = sc_col[sc_col > minreads]
st_means = np.mean(st_expr)   # ST: 발현 세포의 평균
sc_medians = np.median(sc_expr)  # SC: 발현 세포의 중앙값
ratio = st_means / sc_medians
# ratio > 1: Xenium에서 더 잘 탐지, ratio < 1: scRNAseq에서 더 잘 탐지
```

**3.1.4 Region-specific Efficiency**
도메인/영역별로 효율을 분리하여, 특정 조직 영역에서의 segmentation 품질을 평가.

**출력 파일**:
- `efficiency_transcripts_per_cell_comparison.png` - 히스토그램
- `efficiency_genes_per_cell_comparison.png` - 히스토그램
- `efficiency_expression_ratio.png` - 2-panel: Boxplot(log2) + Histogram(log2, cumulative %)
- `efficiency_expression_ratio_by_region.png` - 영역별 ratio
- `efficiency_metrics_comparison.csv` - 정량적 요약

### 3.2 Specificity Analysis - NMP (`analyze_specificity()`, lines 621-757)

**핵심 개념: Negative Marker Purity (NMP)**

NMP는 논문에서 독자적으로 정의한 메트릭으로, 비특이적 co-expression의 정도를 측정한다.

**정의**:
```
NMP = 1 - (X_neg^(sp) - X_neg^(sc)) / X_neg^(sc)   if X_neg^(sp) > X_neg^(sc)
    = 1                                               otherwise
```
여기서:
- `X_neg^(sp)`: spatial 데이터에서 negative marker의 평균 발현
- `X_neg^(sc)`: scRNA-seq reference에서 negative marker의 평균 발현
- Negative marker: 특정 세포 유형에서 발현되지 않아야 하는 유전자

**해석**:
- NMP = 1.0: 완벽한 순도 (비특이적 발현 없음)
- NMP = 0.0: 심각한 누출 (leakage) - segmentation 경계 문제
- NMP > 0.8: 양호 (논문 기준)

**처리 과정**:
```python
def _negative_marker_purity_coexpression(adata, adata_sc):
    # 1. scRNAseq에서 비-co-expressed 유전자 쌍 식별
    # 2. Spatial 데이터에서 해당 쌍의 co-expression 비율 계산
    # 3. NMP = co-expressed 비율이 낮을수록 높음
```

**출력**:
- `specificity_nmp_score.txt` - 전체 NMP 점수
- `specificity_nmp_per_gene_{method}.csv` - 유전자별 NMP
- `specificity_nmp_per_gene_boxplot.png` - NMP 분포
- `specificity_gene_correlation_{method}.png` - Gene-gene correlation heatmap

### 3.3 Positivity Analysis (`analyze_positivity()`, lines 760-1026)

**논문 원래 목적 (노트북 3_5)**:
논문 Methods (p.12): "datasets from all platforms were preprocessed, clustered and annotated. Populations consistently identified across technologies were further used for comparison"
→ 여러 SRT 플랫폼(Xenium, CosMx, Vizgen, MERFISH 등)에서 **동일한 전처리를 적용한 후**, 각 기술에서 독립적으로 클러스터링하여 **같은 세포 유형을 식별**하고, 해당 세포 유형에서의 **유전자별 positivity(양성 세포 비율)를 기술 간 비교**하는 것이 원래 목적. Fig. 2e에 해당.

**파이프라인 현재 구현**:
본 파이프라인은 Xenium 단일 데이터셋만 사용하므로, 크로스 플랫폼 비교는 불가능.
대신 **재분할 vs 원본 비교**를 통해 클러스터링 품질과 유전자별 positivity 패턴을 평가:
- Positivity distribution: Reseg vs Original 비교 (히스토그램)
- Clustering + Violin + UMAP: **양쪽 데이터셋** 전처리 후 나란히 비교 (side-by-side)
- Brain cell type marker 분석: 알츠하이머 패널 주요 마커 발현 비교 (NEW)

**처리 과정**:
```
1. Gene positivity 계산 (Reseg + Original 각각)
   positivity_rate(gene_g) = count(cells with gene_g > 0) / total_cells
   → positivity_dist_comparison.png (두 데이터셋 히스토그램 비교)

2. 양쪽 데이터셋 동일 파이프라인 전처리 (_preprocess_for_positivity)
   for label, ad in datasets:  # Reseg + Original 모두
     filter_cells(min_counts=10, min_genes=3) → normalize_total → log1p
     → neighbors(n_neighbors=8, n_pcs=0) → Leiden (다중 해상도) → UMAP

3. 클러스터별 발현 시각화 (양쪽 비교)
   → violin_clusters.png (positivity 상위 6 유전자 × 데이터셋 side-by-side)
   → umap_top_genes.png (상위 4 유전자 × 데이터셋 side-by-side)
   → optimal_cluster.csv (Reseg만: 유전자별 최고 발현 클러스터)

4. Brain cell type marker 비교 (NEW)
   brain_markers = {Astrocytes, Oligodendrocytes, Neurons(exc/inh), Microglia, OPC}
   → 각 마커의 최적 클러스터에서 데이터셋 간 발현 비교
   → key_markers.csv + key_markers_violin.png
```

**설정값**:
```yaml
comparison:
  positivity:
    n_neighbors: 8           # k-NN 그래프
    n_pcs: 0                 # PCA 사용 안 함 (직접 expression 사용)
    leiden_resolution: 2.2   # 주 해상도
    leiden_resolutions: [2.2, 1.4, 0.6]  # 다중 해상도
    umap_min_dist: 0.1
    min_counts: 10
    min_genes: 3
```

**Brain Cell Type Markers** (하드코딩):
```python
brain_markers = {
    'Astrocytes': ['GFAP', 'AQP4', 'ALDOC', 'S100B', 'SLC1A2', 'SLC1A3'],
    'Oligodendrocytes': ['OLIG2', 'MBP', 'MOG', 'CLDND1', 'ERMN', 'PLP1'],
    'Neurons (excitatory)': ['SYT1', 'SNAP25', 'SLC17A7', 'CAMK2A'],
    'Neurons (inhibitory)': ['GAD1', 'GAD2', 'PVALB', 'SST', 'VIP'],
    'Microglia': ['CSF1R', 'TREM2', 'CX3CR1', 'P2RY12'],
    'OPC': ['VCAN', 'PDGFRA'],
}
```

**시각화**:
- `positivity_dist_comparison.png` - 원본 vs 재분할 positivity 분포 (히스토그램)
- `positivity_violin_clusters.png` - **양쪽 비교**: 상위 6 유전자 × 데이터셋 side-by-side violin (단일 데이터셋이면 상위 10 유전자)
- `positivity_key_markers.csv` - **양쪽 비교**: Brain marker별 최적 클러스터 발현 데이터 (NEW)
- `positivity_key_markers_violin.png` - **양쪽 비교**: Brain marker violin (cell type별 subplot) (NEW)
- `positivity_umap_top_genes.png` - **양쪽 비교**: 상위 4 유전자 × 데이터셋 UMAP
- `positivity_optimal_cluster.csv` - **Reseg만**: 유전자별 최적 클러스터 할당

### 3.4 Diffusion Analysis (`analyze_diffusion()`, lines 1029-1359)

**목적**: Transcript가 세포 centroid로부터 얼마나 멀리 분산되어 있는지 비교

**핵심 - pixel → µm 변환**:
```python
# Xenium: 4.70588 pixels/µm
pixel_to_um = 1 / 4.70588  # ≈ 0.2125 µm/pixel
distances_um = distances_px * pixel_to_um
```

**기술별 변환 팩터**:
| 기술 | pixels/µm | 설명 |
|------|-----------|------|
| Xenium | 4.70588 | 기본값 |
| CosMx | 1.0 | 이미 µm 단위 |
| Vizgen | 1.0 | 이미 µm 단위 |
| MERFISH | 9.2203 | 레퍼런스 기반 |
| HybRISS | 2.5641 | 레퍼런스 기반 |

**시각화**:

| Plot | 설명 | 분석법 |
|------|------|--------|
| `diffusion_complementary_cdf_comparison.png` | P(distance > x) 플롯 | 재분할이 원본보다 더 concentrated한지 확인 (5µm/10µm 참조선 포함) |
| `diffusion_per_gene_ecdf.png` | 상위 유전자별 ECDF | Nuclear vs cytoplasmic 유전자의 diffusion 패턴 |
| `diffusion_gene_method_heatmap.png` | Gene×Method 평균 거리 | 방법 간 유전자별 거리 차이 |
| `diffusion_assigned_reads_barplot.png` | Nuclear capture rate stacked barplot | 핵 내부/외부 transcript 비율 비교 |

**핵심 분석법 - Complementary CDF**:
- X축: distance to centroid (µm)
- Y축: P(distance > x), 즉 해당 거리 이상에 있는 reads 비율
- 재분할 곡선이 원본보다 아래 = 더 concentrated (좋은 segmentation)
- 곡선이 위 = 더 diffused (reads가 멀리 분산됨)

---

## 4. Notebook vs Pipeline 구현 비교

### 4.1 원본 Notebook들
| 노트북 | Pipeline 함수 | 상태 |
|--------|--------------|------|
| `3_3_efficiency_between_methods.ipynb` | `analyze_efficiency()` | 구현됨 |
| `3_4_negative_marker_purity_for_specificity.ipynb` | `analyze_specificity()` | 구현됨 (간소화) |
| `3_5_Computing_positivity_after_preprocessing_for_all_ST_techs.ipynb` | `analyze_positivity()` | 구현됨 |
| `3_6_Diffussion_on_resegmented_data.ipynb` | `analyze_diffusion()` | 구현됨 |
| `3_7_Xenium_vs_Visium_comparison.ipynb` | 미구현 | 미구현 |

### 4.2 구현 상태 요약

| 분석 | Notebook 기능 | Pipeline 구현 | 심각도 |
|------|-------------|-------------|--------|
| Efficiency - histogram | 완전 | 완전 (transcripts + genes per cell) | OK |
| Efficiency - ST/SC ratio | NxN scatter 포함 | 2-panel (boxplot+histogram) + scatter + log2 ratio + region breakdown (raw counts, ST mean/SC median) | OK |
| Efficiency - reseg vs orig | 비교 | boxplot 비교 포함 | OK |
| Specificity - NMP | 유전자별 상세 | inlined NMP + per-gene CSV + boxplot | MEDIUM (약간 간소화) |
| Specificity - efficiency scatter | - | NMP vs Efficiency scatter 추가 | 개선 |
| Positivity - clustering | 다중 해상도 | 다중 해상도 (2.2, 1.4, 0.6), 양쪽 비교 | OK |
| Positivity - brain markers | - | Brain cell type marker 발현 비교 (6 types) | 개선 (NEW) |
| Diffusion - ECDF | 완전 | 완전 (complementary CDF + per-gene ECDF, paper-style) | OK |
| Diffusion - heatmap | - | Gene×Method mean distance heatmap 추가 | 개선 |
| Diffusion - assigned reads | 비율 확인 | nuclear capture rate stacked barplot | OK |
| Xenium vs Visium (3_7) | 완전 | 미구현 (visium_path=null) | LOW (별도 분석) |

### 4.3 알려진 이슈

1. ~~**출력 디렉토리 불일치**~~: **수정됨** — `analyze_efficiency()`도 이제 `figs_dir`에 저장. 모든 Step 4 출력이 동일 디렉토리에 위치하여 specificity의 ratio CSV 탐색도 정상 작동.
2. **sc_reference_path 접근**: ~~혼용~~ **정상** — `analyze_specificity`에 전달되는 `config`는 실제로 `comp_config` (`config['comparison']`)이므로 `config.get('sc_reference_path')` = `comparison.sc_reference_path`로 정상 작동.
3. **Gene name 컬럼 탐지**: `feature_name` → `gene` → `gene_name` → `target` 순으로 시도. Xenium 데이터는 항상 `feature_name`이 존재하여 실질적 문제 없음. 비-Xenium 데이터에서는 'unknown' fallback 가능.
4. **좌표 불일치 탐지**: pixel/µm 좌표 감지 heuristic이 `tx_range < 20000` 기준 사용. Xenium 데이터(µm ~8000-12000, px ~40000-60000)에서는 적절하나, 비표준 데이터에서 부정확할 수 있음.

---

## 5. 전체 출력 목록

### 5.1 Efficiency 분석 출력 (13개)
| # | 파일명 | 유형 | 논문 Figure | 분석법 |
|---|-------|------|-----------|--------|
| 1 | `spatial_roi_map.png` | Spatial | - | 분석 영역 (region별 색상) |
| 2 | `efficiency_transcripts_per_cell_comparison.png` | Histogram | Fig 2b | 재분할의 reads 포착 개선 여부 |
| 3 | `efficiency_genes_per_cell_comparison.png` | Histogram | Fig 2b | 유전자 diversity 비교 |
| 4 | `efficiency_metrics_comparison.csv` | CSV | - | Median/mean 정량적 효율 요약 |
| 5 | `efficiency_expression_ratio.csv` | CSV | - | ST/SC ratio per gene |
| 6 | `efficiency_expression_ratio.png` | 2-Panel (Box+Hist) | Fig 2c | ST/SC ratio: Boxplot(log2) + Histogram(cumulative %) |
| 7 | `efficiency_st_vs_sc_scatter.png` | Scatter | Fig 2c | Log-log ST vs SC with identity line |
| 8 | `efficiency_expression_ratio_by_region.csv` | CSV | Ext Fig 4c | Per-gene per-region ratio |
| 9 | `efficiency_expression_ratio_by_region.png` | Boxplot | Ext Fig 4c | 영역별 효율 차이 |
| 10 | `efficiency_region_breakdown.csv` | CSV | - | Region별 median/mean transcripts/genes |
| 11 | `efficiency_region_breakdown.png` | Boxplot | - | Side-by-side region 비교 |
| 12 | `efficiency_ratio_boxplot_log2.png` | Boxplot | - | log2 ratio 분포 |
| 13 | `efficiency_reseg_vs_original_boxplot.png` | Boxplot | - | 재분할 vs 원본 직접 비교 |

### 5.2 Specificity 분석 출력 (6개)
| # | 파일명 | 유형 | 논문 Figure | 분석법 |
|---|-------|------|-----------|--------|
| 14 | `specificity_nmp_score.txt` | Text | Fig 2d | NMP 전체 점수 (per dataset) |
| 15 | `specificity_nmp_per_gene_{label}.csv` | CSV | - | 유전자별 NMP (per dataset) |
| 16 | `specificity_nmp_per_gene_all.csv` | CSV | - | 유전자별 NMP (모든 dataset 합쳐서) |
| 17 | `specificity_nmp_per_gene_boxplot.png` | Boxplot | Fig 2d | 유전자별 NMP 분포 |
| 18 | `specificity_vs_efficiency_scatter_{label}.png` | Scatter | - | NMP vs Efficiency 관계 |
| 19 | `specificity_gene_correlation_{label}.png` | Heatmap | - | Top 50 gene-gene correlation |

### 5.3 Positivity 분석 출력 (7개)
| # | 파일명 | 유형 | 논문 Figure | 분석법 |
|---|-------|------|-----------|--------|
| 20 | `{tag}_{label}_gene_positivity.csv` | CSV | - | Top 50 gene positivity fractions (Reseg + Original 각각) |
| 21 | `{tag}_positivity_dist_comparison.png` | Histogram | - | Positivity rate 분포 (Reseg vs Original) |
| 22 | `{tag}_positivity_violin_clusters.png` | Violin | - | 클러스터별 positivity (양쪽 side-by-side, 상위 6 유전자) |
| 23 | `{tag}_positivity_key_markers.csv` | CSV | - | Brain cell type marker 발현 데이터 (NEW) |
| 24 | `{tag}_positivity_key_markers_violin.png` | Violin | - | Brain marker 최적 클러스터 발현 비교 (NEW) |
| 25 | `{tag}_positivity_umap_top_genes.png` | UMAP | - | Top 4 유전자 발현 공간 패턴 (양쪽 side-by-side) |
| 26 | `{tag}_positivity_optimal_cluster.csv` | CSV | - | Per-gene 최적 cluster 할당 (Reseg만) |

### 5.4 Diffusion 분석 출력 (8개)
| # | 파일명 | 유형 | 논문 Figure | 분석법 |
|---|-------|------|-----------|--------|
| 27 | `reseg_diffusion_stats.csv` | CSV | - | 재분할 거리 통계 |
| 28 | `original_diffusion_stats.csv` | CSV | - | 원본 거리 통계 |
| 29 | `diffusion_per_gene_summary.csv` | CSV | - | Per-gene mean/median distance |
| 30 | `diffusion_complementary_cdf_comparison.png` | CDF | Fig 2f | 아래=더 concentrated (5µm/10µm 참조선 포함) |
| 31 | `diffusion_per_gene_ecdf.png` | ECDF (9 subplots) | Ext Fig 4f | 유전자별 거리 분포 |
| 32 | `diffusion_gene_method_heatmap.png` | Heatmap | - | Gene×Method mean distance |
| 33 | `diffusion_gene_method_mean_distances.csv` | CSV | - | Gene×Method pivot table |
| 34 | `diffusion_assigned_reads_barplot.png` | Stacked Barplot | Ext Fig 4b | Nuclear capture rate (In Nucleus / Outside Nucleus) |

**총 출력**: 34+ 파일 (21 PNG, 14 CSV, 1 TXT)

---

## 6. Input / Output 상세

### 6.1 Input
| 파일 | 형식 | 설명 |
|------|------|------|
| Step 3 `resegmented.h5ad` | AnnData | 재분할 세포 데이터 |
| Step 3 `transcripts_resegmented.csv` | CSV | 재할당된 transcript |
| Step 0/1 `*.h5ad` | AnnData | 원본 세포 데이터 |
| Step 0 `transcripts.parquet` | Parquet | 원본 transcript |
| scRNAseq reference (선택) | h5ad | NMP/효율 비교용 |

### 6.2 Output
| 파일 | 형식 | 설명 |
|------|------|------|
| `figures/4_techniques_comparison/` | Directory | 모든 분석 결과 |
| 34개+ 시각화/데이터 파일 | PNG/CSV/TXT | 비교 분석 결과 (위 시각화 목록 참조) |
| `step4_done.txt` | Marker | 완료 표시 (pipeline_main에서 생성) |

---

## 7. 관련 설정값 정리

```yaml
comparison:
  run_efficiency: true         # Efficiency 분석 실행
  run_specificity: true        # NMP 분석 실행
  run_positivity: true         # Positivity 분석 실행
  run_diffusion: true          # Diffusion 분석 실행
  technology: "xenium"         # pixel→µm 변환 팩터 결정
  sc_reference_path: null      # scRNAseq reference 경로
  visium_path: null            # Visium 비교 (미구현)
  efficiency_minreads: 1       # 양성 세포 최소 reads

  positivity:
    n_neighbors: 8
    n_pcs: 0
    leiden_resolution: 2.2
    leiden_resolutions: [2.2, 1.4, 0.6]
    umap_min_dist: 0.1
    min_counts: 10
    min_genes: 3
```

### 7.1 파라미터 조정 가이드

| 파라미터 | 현재값 | 범위/옵션 | 영향도 | 조정 가이드 |
|----------|--------|----------|--------|------------|
| **분석 실행 플래그** | | | | |
| `run_efficiency` | `true` | true/false | - | Efficiency 분석. 재분할 vs 원본 transcript 포착률 비교 |
| `run_specificity` | `true` | true/false | - | NMP 분석. scRNAseq reference 필요. 없으면 자동 건너뜀 |
| `run_positivity` | `true` | true/false | - | Positivity 분석. 클러스터링 기반 유전자 양성률 |
| `run_diffusion` | `true` | true/false | - | Diffusion 분석. Transcript-centroid 거리 비교 |
| **핵심 파라미터** | | | | |
| `technology` | `"xenium"` | xenium/cosmx/vizgen/merfish/hybriss | **HIGH** | Pixel→µm 변환 팩터 결정. 잘못 설정하면 모든 거리 계산 오류. Xenium=4.70588 px/µm |
| `sc_reference_path` | `null` | 파일 경로/null | MEDIUM | scRNAseq reference h5ad. NMP (Specificity) 계산에 필수. null이면 NMP 건너뜀 |
| `efficiency_minreads` | `1` | 1-10 | LOW | 양성 세포 판정 최소 reads. 1=어떤 read든 발현으로 간주. 높이면 더 엄격한 양성 기준 |
| **Positivity 서브 파라미터** | | | | |
| `positivity.n_neighbors` | `8` | 5-30 | MEDIUM | Positivity k-NN 그래프. 8=노트북 기본. 높이면 더 부드러운 클러스터, 낮으면 세밀한 구분 |
| `positivity.n_pcs` | `0` | 0/10-50 | LOW | PCA 사용 여부. 0=PCA 없이 직접 expression 사용 (노트북 기본). Xenium 패널은 유전자가 적어 PCA 불필요 |
| `positivity.leiden_resolution` | `2.2` | 0.5-5.0 | MEDIUM | 주 Leiden 해상도. 2.2=노트북 기본. 높으면 더 많은 클러스터, 낮으면 대분류 |
| `positivity.leiden_resolutions` | `[2.2, 1.4, 0.6]` | 리스트 | LOW | 다중 해상도 비교용. 첫 번째 값이 주 분석에 사용 |
| `positivity.umap_min_dist` | `0.1` | 0.01-0.5 | LOW | UMAP 최소 거리. 0.01=촘촘한 클러스터, 0.5=넓은 분포. 시각화만 영향 |
| `positivity.min_counts` | `10` | 1-50 | MEDIUM | Positivity 분석용 셀 필터링. Step 0과 독립 |
| `positivity.min_genes` | `3` | 1-20 | MEDIUM | Positivity 분석용 셀 필터링 |

**튜닝 팁**:
- `technology`는 **반드시 올바르게 설정**. Xenium 데이터면 `"xenium"`, CosMx면 `"cosmx"`. 잘못 설정하면 diffusion 거리 단위가 틀어짐
- `sc_reference_path`를 설정하면 NMP 분석으로 segmentation 품질을 정량적으로 평가 가능. SEA-AD reference (`data/scRNAseq/`)를 사용하면 뇌 데이터에 적합
- Positivity 파라미터들은 노트북 기본값이 대부분 적절. `leiden_resolution`만 세포 유형 수에 맞춰 조정

---

## 8. 핵심 로직의 의미

### 8.1 NMP (Negative Marker Purity)
- **높은 NMP (>0.8)**: Segmentation 품질 양호 - 세포 간 transcript 누출 최소
- **낮은 NMP (<0.6)**: Segmentation 경계 불량 - 이웃 세포의 transcript 혼입
- 논문에서 "all the different SRT technologies presented a mean high specificity (NCP > 0.8)" (p.816)

### 8.2 ST/SC Efficiency Ratio
- **Ratio > 1**: Xenium이 scRNAseq보다 해당 유전자를 더 잘 탐지
- **Ratio = 1**: 동일한 탐지 효율
- **Ratio < 1**: scRNAseq에서 더 잘 탐지 (Xenium 패널에 없거나 탐지 한계)
- 논문: "detection efficiency was found to be between 1.2 and 1.5 times higher than that of scRNA-seq" (p.816)

### 8.3 Diffusion Distance
- 논문에서 "transcripts located more than 10.71 µm, on average, from the cell centroid exhibited a higher gene expression correlation with domain-specific background signatures" (p.817)
- 이 거리 이상의 transcript는 잘못 할당된 것일 가능성이 높음

---

## 10. 시각화 상세 분석 가이드

> **총 19개 PNG 시각화** 전체에 대한 해석 가이드.
> 코드 내 생성 순서에 따라 번호를 매김.

### 10.1 Spatial ROI Map (`spatial_roi_map.png`)

**코드 위치**: `run_step4()`, lines 174-200

**무엇을 봐야 하는가**:
```
    Y (µm)
    ↑
    │ ● ● ●     Region A (파랑)
    │   ● ● ●
    │ ▲ ▲ ▲ ▲   Region B (주황)
    │   ▲ ▲ ▲
    │ ■ ■ ■     Region C (초록)
    │   ■ ■
    └──────────────→ X (µm)

    각 점 = 하나의 세포 centroid
    색상 = 영역 annotation (region_annotation/spatial_annotation/domain)

    ① 영역 분포: 각 region이 공간적으로 잘 분리되어 있는지
    ② 경계선: 영역 간 경계가 명확한지 vs 혼재되어 있는지
    ③ 밀도: 세포 밀도가 영역별로 다른지
```

**해석법**:
- **영역 분리가 명확**: 정상. Region annotation이 공간 패턴을 잘 반영
- **영역이 혼재**: Annotation이 부정확하거나, 실제 조직 구조가 복잡
- **빈 영역**: 해당 region에서 세포가 적거나 segmentation이 해당 영역을 놓침
- **이 plot은 이후 region-based efficiency 분석의 기반** → 영역이 잘 정의되지 않으면 region 분석 결과 해석에 주의

---

### 10.2 Efficiency - Transcripts per Cell (`efficiency_transcripts_per_cell_comparison.png`) → 논문 Fig. 2b

**논문 위치**: Fig. 2b (p.816) - "Transcripts per cell by technology/segmentation"

**무엇을 봐야 하는가**:
```
    Count
    ↑
    │    Nuclei (파랑)    Resegmented (주황)
    │      ██              ████
    │    ██████          ████████
    │  ████████████    ████████████████
    │████████████████████████████████████▓▓
    └──────────────────────────────────────→ Transcripts/cell
      0    50   100  150  200  250  300

    ① 분포 이동 (shift): 재분할이 오른쪽으로 이동 = 더 많은 reads 포착
    ② 중앙값 비교: Reseg > Nuclei = 확장으로 인한 추가 reads
    ③ 꼬리 비교: 긴 오른쪽 꼬리 = 일부 세포가 매우 많은 reads 포착
    ④ 왼쪽 꼬리 (< 10): 저품질 세포 비율 비교
```

**해석법**:
- **재분할 > 원본**: 정상적. 확장으로 인해 cytoplasmic reads가 추가됨
- **재분할 >> 원본 (2배 이상)**: 과도한 확장 → misassignment 의심
- **재분할 ≈ 원본**: 확장이 효과가 없거나, 확장 거리가 너무 짧음
- **논문 수치**: Xenium 평균 ~150-250 reads/cell (패널에 따라 다름)
- **NMP와 교차 분석**: reads 증가 + NMP 유지 = 좋은 확장. reads 증가 + NMP 감소 = misassignment

---

### 10.3 Efficiency - Genes per Cell (`efficiency_genes_per_cell_comparison.png`) → 논문 Fig. 2b

**논문 위치**: Fig. 2b (p.816) - Transcripts/cell과 동일 패널 계열

**무엇을 봐야 하는가**:
```
    Count
    ↑
    │    Nuclei (파랑)    Resegmented (주황)
    │      ██              ████
    │    ██████          ████████
    │  ████████████    ████████████████
    │████████████████████████████████████▓▓
    └──────────────────────────────────────→ Genes/cell
      0    20    40    60    80   100

    ① 분포 이동: 재분할이 오른쪽으로 = 더 다양한 유전자 탐지
    ② 중앙값: Reseg > Nuclei = 확장이 유전자 diversity에 기여
    ③ 최대값: Xenium 패널 크기 (~300-500 genes)에 근접 여부
    ④ 저검출 세포 (< 5 genes): 저품질 세포/debris 비율
```

**해석법**:
- **재분할에서 증가**: cytoplasmic reads 추가로 새로운 유전자가 탐지됨
- **재분할에서 큰 변화 없음**: 확장 영역의 reads가 이미 탐지된 유전자에 집중
- **Transcripts/cell은 증가했는데 Genes/cell은 불변**: 특정 유전자의 발현만 높아짐 (housekeeping 등)
- **Genes/cell 감소**: 있을 수 없으나, 세포 filtering 기준이 다른 경우 가능

---

### 10.4 Efficiency - ST/SC Expression Ratio (`efficiency_expression_ratio.png`) → 논문 Fig. 2c

**논문 위치**: Fig. 2c (p.816) - "SRT/scRNAseq gene efficiency ratio"

**구현**: 2-패널 Figure (Paper Fig. 2c 스타일)
- **좌측 패널**: Boxplot + Stripplot (log2 스케일) — 각 점 = 유전자, summary stats 주석
- **우측 패널**: Histogram + Cumulative % overlay (twin axis)
- 색상: Xenium purple (#56018f)

**무엇을 봐야 하는가**:
```
    좌측 패널 (Boxplot)              우측 패널 (Histogram + Cumulative %)
    log2(SRT/SC ratio)               Number of Genes          Cum %
    ↑                                ↑                         ↑
    │     ·                          │     ████                │ 100%
    2.0 │  ╭─╮                       │   ████████         ───  │  80%
    │     │ │  ← summary stats       │ ████████████     ╱     │  60%
    1.0 │  │━│     n=xxx genes        │████████████████╱       │  40%
    │     │ │     Median: x.xx        │█████████████╱          │  20%
    0.0 │──│─│── Ratio=1 (빨강)       │─ˡ──────────╱──────→   │   0%
    │     ╰─╯     Over(>1): xx%      │  ↑        ↑
   -1.0 │                              빨강     주황 (median)
    │                                 log2(SRT/SC ratio)
    └──Xenium──→

    ① 좌측 boxplot: 전체 분포 요약 + 개별 유전자 점
    ② 우측 histogram: log2 ratio 분포 + 누적 % 오버레이 (녹색 곡선)
    ③ 빨간 실선 (y=0 / x=0): log2(1) = 0 = 동일 효율 기준
    ④ 주황 점선: 중앙값 위치
    ⑤ Summary stats: n, median, mean, Over/Under-detected %, Within 0.5-2x %
```

**해석법**:
- **중앙값 > 0 (log2)**: Xenium ISH가 대부분 유전자에서 scRNAseq보다 높은 탐지 효율
- **논문 결과**: "detection efficiency was between 1.2 and 1.5 times higher than scRNA-seq" (p.816) → log2 중앙값 ≈ 0.26~0.58 기대
- **Over-detected (>1) 비율이 높음**: 대부분 유전자에서 Xenium이 우세
- **우측 누적 곡선 (녹색)**: 50% 지점이 0 이상이면 과반수 유전자가 Xenium 우세
- **분포가 넓음**: 유전자별 탐지 효율 편차가 크다 → 프로브 품질 차이
- **참고**: 단독 log2 boxplot은 10.8 (`efficiency_ratio_boxplot_log2.png`)에서 별도 제공

---

### 10.5 Efficiency - ST vs SC Scatter (`efficiency_st_vs_sc_scatter.png`) → 논문 Fig. 2c

**논문 위치**: Fig. 2c (p.816) - ST/scRNAseq expression comparison

**무엇을 봐야 하는가**:
```
    log10(ST mean, raw)
    ↑
    4 │                    ·  ·
      │                 · · ·
    3 │              · · · ·
      │           · · · ·      ← ST 과탐지 영역
    2 │        · · · · ·
      │     · · · · ·  /
    1 │  · · · · · · / ← identity line (빨간 점선)
      │  · · · · · /
    0 │  · · · · /
      │  · · · /  ← SC 과탐지 영역
      └────────────────────→ log10(scRNA median, raw)
         0    1    2    3    4

    각 점 = 하나의 유전자
    ① Identity line 위: Xenium이 더 효율적으로 탐지
    ② Identity line 아래: scRNAseq가 더 효율적
    ③ 점 분산: 넓으면 유전자별 효율 차이 큼
    ④ 기울기: 1보다 크면 전반적으로 Xenium이 우세
```

**해석법**:
- **대부분 identity 위**: Xenium ISH가 scRNAseq보다 유전자 탐지 효율이 높음 (논문 결론)
- **일부 유전자가 크게 아래**: 해당 유전자의 Xenium 프로브 효율이 낮거나, scRNAseq에서 특이적으로 높은 발현
- **전체 패턴**: 강한 양의 상관 = 두 기술 간 일관성, 분산 크면 기술 간 bias 존재
- **log-log 스케일**: 저발현/고발현 유전자 모두 균등하게 시각화

---

### 10.6 Efficiency - Region Expression Ratio (`efficiency_expression_ratio_by_region.png`) → 논문 Ext Fig. 4c

**논문 위치**: Extended Data Fig. 4c (p.829) - "SRT/SC efficiency by brain region"

**무엇을 봐야 하는가**:
```
    log2(ST/SC ratio)
    2.0 ┤
        │  ╭─╮     ╭─╮     ╭─╮
    1.0 ┤  │━│     │━│     │ │
        │  │ │     │ │     │━│
    0.0 ┤──│─│─────│─│─────│─│──── ← 동일 효율 기준선
        │  │ │     │ │     │ │
   -1.0 ┤  ╰─╯     ╰─╯     ╰─╯
        └──Cortex──Hippo──Thalamus──→

    ① 각 boxplot = 해당 영역의 유전자별 ratio 분포
    ② 빨간 점선 (y=0): ratio = 1 (동일 효율) 기준선
    ③ 중앙값 위치: 0 위면 해당 영역에서 Xenium이 우세
    ④ 영역 간 차이: 특정 영역에서 효율이 더 높은지
```

**해석법**:
- **모든 영역이 0 이상**: 전체적으로 Xenium이 효율적 (논문 일치)
- **특정 영역만 0 이하**: 해당 영역에서 segmentation 품질 이슈 가능
- **영역 간 ratio 차이 큼**: 조직 유형에 따라 탐지 효율이 다름 (세포 크기, 밀도 영향)
- **논문**: "efficiency was consistent across hippocampus, cortex, and thalamus" (Ext Fig. 4c)

---

### 10.7 Efficiency - Region Breakdown (`efficiency_region_breakdown.png`)

**무엇을 봐야 하는가**:
```
    Transcripts/cell              Genes/cell
    ↑                             ↑
    │  ╭─╮  ╭─╮  ╭─╮            │  ╭─╮  ╭─╮  ╭─╮
    │  │ │  │ │  │━│            │  │━│  │ │  │ │
    │  │━│  │━│  │ │            │  │ │  │━│  │━│
    │  │ │  │ │  │ │            │  │ │  │ │  │ │
    │  ╰─╯  ╰─╯  ╰─╯            │  ╰─╯  ╰─╯  ╰─╯
    └──Ctx──Hippo──Thal──→       └──Ctx──Hippo──Thal──→

    좌측: 영역별 세포당 transcript 수 boxplot
    우측: 영역별 세포당 유전자 수 boxplot
    ※ Resegmented 데이터만 대상 (영역 annotation 필요)
```

**해석법**:
- **영역별 차이**: 조직 구조에 따른 자연스러운 차이 (신경 vs 글리아 비율 등)
- **높은 중앙값 영역**: 세포가 크거나 전사 활성이 높은 영역
- **이상치**: 극단적 세포 = 합쳐진 세포(doublet) 또는 debris
- **Region ROI map과 교차 확인**: 공간적 위치와 발현 패턴의 일치 여부

---

### 10.8 Efficiency - Log2 Ratio Boxplot (`efficiency_ratio_boxplot_log2.png`)

**무엇을 봐야 하는가**:
```
    log2(ST mean / scRNA median)
    ↑
    │        ·
    │     ╭──╮
    2.0 │  │  │     ← 상위 25%: Xenium이 2배 이상 효율
    │     │  │
    1.0 │  │━━│     ← 중앙값: 대부분 유전자에서 Xenium 우세
    │     │  │
    0.0 │──│──│──── ← 동일 효율 기준선
    │     ╰──╯
   -1.0 │     ·     ← 이상치: scRNAseq에서 더 잘 탐지되는 유전자
    │
    └──────────→

    Box: IQR (25-75 percentile)
    Whiskers: 1.5×IQR
    Stripplot: 개별 유전자 점
    빨간 점선: log2(1) = 0 (동일 효율)
```

**해석법**:
- **중앙값 > 0**: 전반적으로 Xenium이 scRNAseq보다 효율적 (논문: log2(1.2)~log2(1.5) ≈ 0.26~0.58)
- **분포 폭**: 좁으면 유전자 간 균일, 넓으면 프로브별 편차 큼
- **음수 이상치**: 해당 유전자의 Xenium 프로브 성능이 낮음
- **10.4의 histogram과 동일 데이터를 boxplot으로** → 요약 통계에 집중

---

### 10.9 Efficiency - Reseg vs Original Boxplot (`efficiency_reseg_vs_original_boxplot.png`)

**무엇을 봐야 하는가**:
```
    Counts/cell                 Genes/cell
    ↑                           ↑
    │  ╭─╮     ╭─╮             │  ╭─╮     ╭─╮
    │  │ │     │ │             │  │ │     │ │
    │  │ │     │━│             │  │━│     │━│
    │  │━│     │ │             │  │ │     │ │
    │  │ │     │ │             │  │ │     │ │
    │  ╰─╯     ╰─╯             │  ╰─╯     ╰─╯
    └──Reseg───Original──→     └──Reseg───Original──→

    좌측: 세포당 총 counts 비교
    우측: 세포당 탐지 유전자 수 비교
    ① Reseg > Original: 재분할이 더 많은 reads/genes를 포착
    ② 중앙값 차이: 확장 효과의 크기
    ③ 분포 폭: 재분할에서 넓어지면 세포 간 변동성 증가
```

**해석법**:
- **Reseg 중앙값 > Original**: 정상. cytoplasmic 확장으로 추가 reads 포착
- **큰 차이 (>2배)**: 과도한 확장 가능 → NMP 확인 필요
- **Genes/cell이 거의 같은데 Counts/cell만 증가**: 이미 알려진 유전자의 reads만 추가 (depth 증가)
- **둘 다 증가**: 새로운 유전자도 추가 탐지 (diversity + depth 모두 개선)
- **10.2/10.3 histogram의 요약 버전** → 빠른 비교에 적합

---

### 10.10 Specificity - NMP Boxplot (`specificity_nmp_per_gene_boxplot.png`) → 논문 Fig. 2d

**논문 위치**: Fig. 2d (p.816) - "Gene specificity (NCP/NMP) by technology"

**무엇을 봐야 하는가**:
```
    NMP Score
    1.0 ┤────────────── 완벽한 순도 ──
        │     ╭─╮
    0.9 ┤     │━│      ← 대부분 유전자가 0.8-1.0
        │     │ │
    0.8 ┤─────│─│────── 논문 기준선 (NCP > 0.8 = 양호) ──
        │     ╰─╯
    0.7 ┤        ·     ← 이상치 (outlier): 특정 유전자가 낮은 NMP
        │
    0.6 ┤
        │
    0.5 ┤              ← 0.5 이하 = 심각한 누출
        └──Nuclei──Reseg──→

    ① 중앙값: 0.8 이상이면 양호
    ② 하위 이상치: 어떤 유전자가 낮은 NMP인지 식별
    ③ Nuclei vs Reseg: 재분할이 NMP를 악화시키지 않는지 확인
```

**해석법**:
- **NMP > 0.8**: 논문 기준 양호 → segmentation 경계가 적절
- **NMP < 0.6**: 심각한 transcript 누출 → 세포 간 경계 불량
- **Nuclei NMP > Reseg NMP**: 확장이 과도하여 이웃 세포의 reads 혼입
- **Nuclei NMP < Reseg NMP**: 드문 경우. 재분할이 더 정확한 경계를 설정
- **낮은 NMP 유전자 분석**: 해당 유전자가 diffuse하게 발현되는지 확인 (생물학적 원인 vs 기술적 원인)

---

### 10.11 Specificity - Efficiency vs Specificity Scatter (`specificity_vs_efficiency_scatter_{label}.png`)

**무엇을 봐야 하는가**:
```
    Specificity (Purity)
    1.0 ┤─────────────────────── 완벽한 순도
        │  · · · · · · · · ·
    0.9 │  · · · · · · · · ·   ← 고효율 + 고순도 = 이상적
        │  · · · · · · · ·
    0.8 │──·─·─·─·─·─·────── NMP 기준선
        │    · · · ·
    0.7 │      · ·        ← 고효율 + 저순도 = 확장 과도
        │        ·
    0.6 │  ·               ← 저효율 + 저순도 = 프로브 문제
        │
    0.5 │
        └────────────────────→ Efficiency Ratio (ST/scRNAseq)
           0    1    2    3

    각 점 = 하나의 유전자
    빨간 점선 (x=1): 동일 효율 기준
    회색 점선 (y=1): 완벽한 순도 기준
```

**해석법**:
- **우상단 (높은 효율 + 높은 순도)**: 이상적 유전자 — Xenium에서 잘 탐지되면서 특이적
- **우하단 (높은 효율 + 낮은 순도)**: 해당 유전자가 많이 탐지되지만 비특이적 co-expression 존재 → 확장에 의한 leakage 의심
- **좌상단 (낮은 효율 + 높은 순도)**: 탐지량은 적지만 특이적 → 프로브 효율은 낮으나 정확
- **좌하단 (낮은 효율 + 낮은 순도)**: 프로브 문제 — 탐지도 낮고 특이성도 낮음
- **Reseg vs Original 별도 생성**: 재분할 후 유전자들이 우하단으로 이동하면 확장이 과도

---

### 10.12 Specificity - Gene Correlation Heatmap (`specificity_gene_correlation_{label}.png`)

**무엇을 봐야 하는가**:
```
              Gene1  Gene2  Gene3  Gene4  Gene5
    Gene1     1.0    0.8    0.1   -0.2    0.3
    Gene2     0.8    1.0    0.2   -0.1    0.4
    Gene3     0.1    0.2    1.0    0.7    0.0
    Gene4    -0.2   -0.1    0.7    1.0   -0.1
    Gene5     0.3    0.4    0.0   -0.1    1.0

    색상: 빨강=양의 상관, 파랑=음의 상관, 흰색=무상관
    대상: 발현량 상위 50개 유전자
    Pearson correlation (cell × gene matrix)
```

**해석법**:
- **강한 양의 상관 블록**: 같은 세포 유형에서 co-expressed되는 유전자 그룹
- **강한 음의 상관**: 서로 다른 세포 유형의 marker (교차 발현 없음 = 좋은 segmentation)
- **Reseg vs Original 비교**: 재분할 후 양의 상관이 전반적으로 증가하면 → leakage로 인한 비특이적 co-expression 증가
- **모든 유전자 간 양의 상관**: ambient RNA 오염 또는 과도한 확장 신호
- **NMP가 reference가 없을 때 대안적 proxy**: Gene-gene correlation 패턴으로 segmentation 품질 간접 평가

---

### 10.13 Positivity - Distribution Comparison (`{tag}_positivity_dist_comparison.png`)

**무엇을 봐야 하는가**:
```
    Number of Genes
    ↑
    │     Nuclei (파랑)     Reseg (주황)
    │  ████
    │  ████████           ████
    │  ████████████       ████████
    │  ████████████████   ████████████████
    └──────────────────────────────────────→ Fraction of Positive Cells
      0    0.2   0.4   0.6   0.8   1.0

    X축: 양성 세포 비율 (positivity rate = 해당 유전자 발현 세포 / 전체 세포)
    Y축: 해당 비율을 가진 유전자의 수

    ① 왼쪽 피크 (0-0.2): 세포 유형 특이적 유전자 (소수 세포에서만 발현)
    ② 오른쪽 피크 (0.8-1.0): housekeeping 유전자 (대부분 세포에서 발현)
    ③ 중간대 (0.3-0.7): 여러 세포 유형에 걸쳐 발현
    ④ Nuclei vs Reseg 비교: 재분할로 positivity가 이동하는지
```

**해석법**:
- **재분할에서 전체적으로 오른쪽 이동**: 확장으로 더 많은 세포에서 유전자 탐지 = 정상
- **큰 이동 (>0.2)**: 과도한 확장 → 비특이적 유전자 할당 가능
- **bimodal 분포**: 세포 유형 특이적 + housekeeping 유전자가 잘 분리됨 = 좋은 데이터
- **unimodal (한쪽 치우침)**: 대부분 유전자가 비슷한 positivity = 세포 유형 구분 약함

---

### 10.14 Positivity - Cluster Violin (`{tag}_positivity_violin_clusters.png`)

**논문 맥락**: 노트북 3_5에서는 여러 SRT 플랫폼을 각각 클러스터링한 후, 같은 세포 유형(예: astrocyte)을 식별하여 기술 간 positivity를 비교 (Fig. 2e). 본 파이프라인은 Xenium 단일 데이터이므로, **재분할 vs 원본 데이터의 클러스터링 비교** 용도로 활용.

**구현**: 양쪽 데이터셋 동일 파이프라인 전처리 후 side-by-side 비교.
- 2 데이터셋: `n_genes_violin(=6) × len(processed)(=2)` 그리드 (행=유전자, 열=데이터셋)
- 1 데이터셋 fallback: `10 × 1` 세로 배열

**무엇을 봐야 하는가**:
```
    [Resegmented]                          [Original]
    Gene1 (pos=0.85)                      Gene1 (pos=0.85)
    Expression                             Expression
    ↑  ╭──╮  ╭╮    ╭──╮                  ↑  ╭──╮  ╭╮    ╭──╮
    │  │  │  ││    │  │                  │  │  │  ││    │  │
    │  ╰──╯  ╰╯    ╰──╯                  │  ╰──╯  ╰╯    ╰──╯
    └──C0──C1──C2──C3──→                  └──C0──C1──C2──C3──→

    Gene2 (pos=0.12)                      Gene2 (pos=0.12)
    ...                                    ...

    6개 유전자 × 2 데이터셋 side-by-side (양쪽 비교 시)
    X축: Leiden 클러스터 (primary resolution, 기본 2.2)
    Y축: log1p 정규화 발현량
    Violin 폭 = 해당 클러스터 내 발현 분포
```

**해석법**:
- **특정 클러스터에서만 고발현**: 세포 유형 특이적 marker → 좋은 clustering
- **모든 클러스터에서 균일**: housekeeping gene 또는 clustering 해상도 부적절
- **Reseg vs Original 비교**: 같은 유전자가 두 데이터셋에서 유사한 클러스터 패턴을 보이는지
- **Reseg에서만 추가 발현**: 확장으로 인해 추가 transcript가 포착된 경우
- **bimodal violin**: 해당 클러스터 내에 두 하위 집단이 존재 → 해상도를 높여야 할 수 있음
- **다중 해상도 (2.2/1.4/0.6)**: obs에 `leiden_2_2`, `leiden_1_4`, `leiden_0_6` 모두 저장되나, violin은 primary (첫 번째 = 2.2)만 사용

### 10.14b Positivity - Key Brain Markers Violin (`{tag}_positivity_key_markers_violin.png`) (NEW)

**무엇을 봐야 하는가**:
```
    GFAP              AQP4              OLIG2             MBP
    (Astrocytes)      (Astrocytes)      (Oligo)           (Oligo)
    ┌──────┐          ┌──────┐          ┌──────┐          ┌──────┐
    │ ╭╮ ╭╮│          │ ╭╮ ╭╮│          │ ╭╮ ╭╮│          │ ╭╮ ╭╮│
    │ ││ │││          │ ││ │││          │ ││ │││          │ ││ │││
    │ ╰╯ ╰╯│          │ ╰╯ ╰╯│          │ ╰╯ ╰╯│          │ ╰╯ ╰╯│
    └Res─Orig┘        └Res─Orig┘        └Res─Orig┘        └Res─Orig┘

    각 subplot = 하나의 brain marker 유전자
    X축: Resegmented vs Original
    Y축: 최적 클러스터에서의 발현량
    색상: Reseg=#56018f (보라), Original=#4DA1A9 (청록)
```

**해석법**:
- **Reseg에서 발현 증가**: 확장으로 인해 해당 세포 유형의 transcript가 추가 포착
- **Reseg에서 발현 감소**: 드문 경우, 클러스터 구성이 변경됨
- **Cell type별 일관성**: 같은 cell type의 marker들이 유사한 패턴을 보여야 정상
- **최적 클러스터**: 각 유전자별로 평균 발현이 가장 높은 클러스터에서의 세포만 사용

---

### 10.15 Positivity - UMAP Top Genes (`{tag}_positivity_umap_top_genes.png`)

**무엇을 봐야 하는가**:
```
    [Resegmented]
    Gene1 [Reseg]    Gene2 [Reseg]    Gene3 [Reseg]    Gene4 [Reseg]
    ┌──────┐          ┌──────┐          ┌──────┐          ┌──────┐
    │··▓▓··│          │ ▓▓···│          │···▓▓·│          │·▓▓▓··│
    │·▓▓▓··│          │▓▓▓···│          │··▓▓▓·│          │··▓▓··│
    └──────┘          └──────┘          └──────┘          └──────┘

    [Original]
    Gene1 [Orig]     Gene2 [Orig]     Gene3 [Orig]     Gene4 [Orig]
    ┌──────┐          ┌──────┐          ┌──────┐          ┌──────┐
    │··▓▓··│          │ ▓▓···│          │···▓▓·│          │·▓▓▓··│
    │·▓▓▓··│          │▓▓▓···│          │··▓▓▓·│          │··▓▓··│
    └──────┘          └──────┘          └──────┘          └──────┘

    n_datasets(=2) 행 × 4열 그리드 (양쪽 비교 시 2×4=8 subplots)
    각 점 = 하나의 세포 (UMAP 공간)
    색상 강도 = 해당 유전자의 발현량
    각 subplot 제목: "{gene} [{dataset}]"

    ① 발현 패턴: 특정 UMAP 영역에 집중 vs 전체 분산
    ② 데이터셋 간 비교: 같은 유전자가 Reseg vs Original에서 유사한 공간 패턴인지
    ③ 발현 수준: 강도가 클러스터 내에서 균일한지 vs 이질적
```

**해석법**:
- **특정 클러스터에 집중**: 세포 유형 특이적 marker → 해당 클러스터의 identity 파악 가능
- **넓게 분산**: housekeeping gene 또는 많은 세포 유형에서 공통 발현
- **Reseg vs Original 비교**: 같은 유전자의 UMAP 패턴이 일치하면 재분할이 세포 유형 구조를 유지한 것
- **Reseg에서 더 넓은 발현**: 확장으로 인해 추가 세포에서 유전자 탐지 → positivity 증가
- **UMAP 구조와 violin plot 교차 확인**: violin에서 고발현인 클러스터가 UMAP의 같은 영역인지 확인

---

### 10.16 Diffusion - Complementary CDF (`diffusion_complementary_cdf_comparison.png`) → 논문 Fig. 2f

**논문 위치**: Fig. 2f (p.816) - "Cumulative proportion of reads by distance"

**무엇을 봐야 하는가**:
```
    P(distance > x)
    1.0 ├──╲
        │    ╲  Nuclei (파랑)
    0.8 │     ╲╲
        │       ╲╲  Reseg (주황)
    0.6 │        ╲╲
        │          ╲╲
    0.4 │            ╲╲
        │              ╲──╲
    0.2 │                ╲───╲
        │                  ╲────╲
    0.0 ├────────────────────╲────╲───
        └──────────────────────────→ Distance to centroid (µm)
           0     5    10    15    20

    핵심: 곡선이 아래에 있을수록 transcript가 centroid에 더 가까움 (좋은 segmentation)

    ① Reseg 곡선이 Nuclei보다 위 = 재분할이 더 diffused (확장 효과)
    ② 10.71 µm 지점 확인: 이 거리에서의 P(d>x) 값
    ③ 두 곡선의 교차점: 어디서 패턴이 반전되는지
```

**해석법**:
- **Complementary CDF = 1 - ECDF**: 해당 거리 이상에 있는 reads의 비율
- **곡선이 빠르게 0에 접근**: reads가 centroid 근처에 집중 (compact segmentation)
- **곡선이 천천히 감소**: reads가 넓게 분산 (diffuse segmentation 또는 과도한 확장)
- **Reseg vs Nuclei 비교**: Reseg가 약간 위에 있는 것은 정상 (cytoplasmic reads 포함). 하지만 크게 위에 있으면 misassignment
- **논문의 turnover distance (10.71 µm)**: 이 거리에서 P(d>x)가 아직 높으면 주의

---

### 10.17 Diffusion - Per Gene ECDF (`diffusion_per_gene_ecdf.png`) → 논문 Extended Data Fig. 4f

**논문 위치**: Extended Data Fig. 4f (p.829) - "Density: distance to centroid per gene"

**무엇을 봐야 하는가**:
```
    ECDF
    1.0 ├─────────────────────────────
        │   Nuclear ╱────── fast rise
    0.8 │          ╱
        │         ╱
    0.6 │        ╱  Cytoplasmic
        │       ╱         ╱───── slow rise
    0.4 │      ╱         ╱
        │     ╱         ╱
    0.2 │    ╱         ╱
        │   ╱         ╱
    0.0 ├──╱─────────╱
        └────────────────────→ Distance (µm)
           0    5   10   15

    각 곡선 = 하나의 유전자
    ① Nuclear gene: 왼쪽으로 치우침 (centroid 근처)
    ② Cytoplasmic gene: 오른쪽으로 치우침 (먼 거리)
    ③ Nuclei vs Reseg 비교: 같은 유전자의 두 방법 비교
```

**해석법**:
- 같은 유전자의 Nuclei vs Reseg ECDF 비교: Reseg에서 곡선이 오른쪽으로 이동하면 확장으로 인한 추가 reads
- Nuclear gene의 이동이 작고, Cytoplasmic gene의 이동이 크면 → 적절한 확장
- 모든 유전자의 ECDF가 동일하게 이동하면 → 비특이적 확장 (misassignment 가능)

---

### 10.18 Diffusion - Gene×Method Heatmap (`diffusion_gene_method_heatmap.png`)

**무엇을 봐야 하는가**:
```
              Nuclei    Reseg    (Baysor)
    Gene1      5.2       6.8      5.5      ← Nuclear gene: 변화 적음
    Gene2      7.8       12.1     8.2      ← Cyto gene: 재분할에서 증가
    Gene3      4.1       5.0      4.5
    Gene4      9.5       14.3     10.1     ← 가장 큰 변화 = 확장 민감
    Gene5      6.3       8.2      6.8

    값 = 평균 distance to centroid (µm)
    색상: 어두움 = 짧은 거리, 밝음 = 긴 거리

    ① 방법 간 패턴: Nuclear gene은 방법에 무관하게 일관
    ② Reseg에서 크게 증가하는 유전자: cytoplasmic/diffuse 유전자
    ③ 전체 패턴: 재분할이 모든 유전자에서 균일하게 증가하면 비특이적 확장
```

**해석법**:
- **Nuclear genes (Opalin 등)**: 모든 방법에서 5-7 µm → 방법에 무관하게 일관
- **Cytoplasmic genes (MAG 등)**: Nuclei에서 8-10 µm → Reseg에서 12-15 µm로 증가
- **차이 패턴**: cytoplasmic gene만 증가하고 nuclear gene은 유지 → 좋은 확장
- **모든 유전자 균일 증가**: 비특이적 확장 또는 큰 expansion distance → NMP 확인 필요

---

### 10.19 Assigned Reads Barplot (`diffusion_assigned_reads_barplot.png`) → 논문 Extended Data Fig. 4b

**논문 위치**: Extended Data Fig. 4b (p.829) - "Proportion of assigned reads"

**구현**: Nuclear capture rate (Stacked Barplot)
- **Resegmented**: `in_cell > 0` (Cellpose nuclear mask)
- **Original**: `overlaps_nucleus > 0` (Xenium native nuclear overlap flag)
- 색상: In Nucleus=#4DA1A9 (청록), Outside Nucleus=#e8e8e8 (회색)

**무엇을 봐야 하는가**:
```
    Fraction
    1.0 ├────────────────────────
        │  ░░░░░░░░   ░░░░░░░░░░░░  ← Outside Nucleus (회색)
    0.8 │  ░░░░░░░░   ░░░░░░░░░░░░
        │  ████████   ████████████
    0.6 │  ████████   ████████████
        │  ████████   ████████████  ← In Nucleus (청록)
    0.4 │  ████████   ████████████
        │  ████████   ████████████
    0.2 │  ████████   ████████████
        └──Reseg─────Original──────→

    제목: "Proportion of Reads Inside Nuclei (Ext. Data Fig. 4b style)"
    ① In Nucleus 비율: 핵 내부에 위치한 transcript 비율
    ② Stacked bar: In Nucleus + Outside Nucleus = 1.0
    ③ Reseg vs Original 비교: 핵 포착률 차이
```

**해석법**:
- **In Nucleus 10-30%**: 정상 — 핵 내부에만 위치한 transcript는 전체의 소수
- **Reseg < Original**: 정상. 확장으로 인해 핵 외부 reads가 추가되어 핵 비율은 상대적으로 감소
- **Reseg ≈ Original**: 확장이 적거나 핵 밖 transcript가 적은 경우
- **In Nucleus > 50%**: 패널 유전자가 대부분 nuclear-enriched이거나 확장이 매우 작은 경우
- **참고**: 이 metric은 "cell assignment rate"가 아닌 "nuclear capture rate"임. 논문의 76.8% 수치는 cell assignment(세포 할당) 비율이며, 본 barplot은 핵 overlap 비율을 보여줌

---

## 9. Summary

Step 4는 **정량적 품질 비교**의 핵심 단계로, 4가지 메트릭 (Efficiency, Specificity, Positivity, Diffusion)을 통해 재분할의 개선 효과를 평가한다.

| 분석 | 구현 상태 | 심각도 |
|------|---------|--------|
| Efficiency - histograms | 완전 (transcripts + genes) | OK |
| Efficiency - ST/SC ratio | 완전 (scatter + boxplot + log2 + region) | OK |
| Efficiency - reseg vs orig | 완전 (boxplot 비교) | OK |
| Specificity - NMP | 구현됨 (per-gene CSV + boxplot) | MEDIUM (약간 간소화) |
| Specificity - NMP vs Efficiency | 완전 (scatter plot) | OK |
| Positivity - clustering | 완전 (다중 해상도, 양쪽 비교) | OK |
| Positivity - brain markers | 완전 (6 cell type, 양쪽 비교) | OK (NEW) |
| Diffusion - CDF/ECDF | 완전 (paper-style, 5/10µm 참조선) | OK |
| Diffusion - heatmap | 완전 (Gene×Method) | OK |
| Diffusion - assigned reads | 완전 (nuclear capture rate stacked barplot) | OK |
| Xenium vs Visium (3_7) | 미구현 (visium_path=null) | LOW (별도 분석) |

**구현 완성도**: **HIGH** - 핵심 기능 모두 구현. 34+ 출력 파일 생성. Xenium vs Visium 비교만 미구현 (별도 데이터 필요).

**알려진 이슈**:
1. ~~출력 디렉토리 불일치~~ — **수정됨**: 모든 출력이 `figures/4_techniques_comparison/`에 통합
2. ~~sc_reference_path config 접근 경로 혼용~~ — **수정됨**: `comp_config.get('sc_reference_path')`로 정상 작동
3. 좌표 단위 (µm/pixel) 자동 감지 heuristic (centroid/tx range ratio 기반) — LOW (비표준 데이터에서 부정확할 수 있음)
4. Brain marker 유전자 목록이 하드코딩 — LOW (알츠하이머 패널 전용, 다른 데이터에서는 누락 가능)

---

## 11. 논문 Figure 직접 대응 및 시각화 정상 판별 종합 가이드

### 11.1 파이프라인 출력 → 논문 Figure 매핑 종합표

| # | 파이프라인 출력 파일 | 논문 Figure | 논문 페이지 | 논문 원문 설명 |
|---|---|---|---|---|
| 1 | `efficiency_transcripts_per_cell_comparison.png` | **Fig. 2b** (p.816) 좌측 | 816 | "Box plot showing the numbers of transcripts per cell and genes per cell for each dataset" |
| 2 | `efficiency_genes_per_cell_comparison.png` | **Fig. 2b** (p.816) 우측 | 816 | 세포당 유전자 수 비교 |
| 3 | `efficiency_expression_ratio.png` | **Fig. 2c** (p.815) | 815 | 2-panel: Boxplot(log2) + Histogram(cumulative %). Raw counts: ST mean / SC median |
| 4 | `efficiency_st_vs_sc_scatter.png` | **Fig. 2c** 관련 | 815 | ST vs SC 발현량 log-log scatter |
| 5 | `efficiency_expression_ratio_by_region.png` | **Extended Data Fig. 4c** | 829 | 뇌 영역별 효율 비교 |
| 6 | `specificity_nmp_per_gene_boxplot.png` | **Fig. 2d** (p.816) | 816 | "Box plot showing NCP scores, ranging from 0 to 1" |
| 7 | `specificity_vs_efficiency_scatter_{label}.png` | 직접 대응 없음 (파이프라인 추가) | - | NMP vs Efficiency 관계 scatter |
| 8 | `specificity_gene_correlation_{label}.png` | 직접 대응 없음 | - | Gene-gene correlation heatmap |
| 9 | `positivity_dist_comparison.png` | **Fig. 2e** 관련 | 816 | Positivity rate 분포 비교 |
| 10 | `positivity_violin_clusters.png` | **Fig. 2e** (p.816) 관련 | 816 | "Violin plot of transcripts detected per gene across datasets" (양쪽 side-by-side) |
| 10b | `positivity_key_markers_violin.png` | 직접 대응 없음 (파이프라인 추가) | - | Brain cell type marker 발현 비교 (NEW) |
| 11 | `positivity_umap_top_genes.png` | 직접 대응 없음 | - | 상위 유전자 발현 공간 패턴 (양쪽 side-by-side) |
| 12 | `diffusion_complementary_cdf_comparison.png` | **Fig. 2f** (p.816) | 816 | "Cumulative proportion of reads by distance from the cell centroid" |
| 13 | `diffusion_per_gene_ecdf.png` | **Extended Data Fig. 4f** | 829 | "Density: distance to centroid per gene" |
| 14 | `diffusion_gene_method_heatmap.png` | 직접 대응 없음 (파이프라인 추가) | - | Gene × Method mean distance |
| 15 | `diffusion_assigned_reads_barplot.png` | **Extended Data Fig. 4b** | 829 | "Assigned reads proportion" |
| 16 | `efficiency_reseg_vs_original_boxplot.png` | 직접 대응 없음 (파이프라인 추가) | - | 재분할 vs 원본 직접 비교 |
| 17 | `efficiency_ratio_boxplot_log2.png` | 직접 대응 없음 (파이프라인 추가) | - | log2 ratio 요약 boxplot |

### 11.2 정상 결과 판별 체크리스트

#### 시각화 1: Efficiency - Transcripts/cell (`efficiency_transcripts_per_cell_comparison.png`) → Fig. 2b
- [ ] **중앙값**: 150-250 reads/cell (논문: "186.6 reads per cell")
- [ ] **재분할 > 원본**: 재분할 분포가 오른쪽으로 이동 (더 많은 reads 포착)
- [ ] **중앙값 차이**: 재분할이 20-50% 더 높음이 정상
- **왜 정상인가**: 논문 Fig. 2b에서 확장(expansion)이 reads/cell을 증가시킴. "resulting in very different numbers of reads per cell" (p.816). Cellpose resegmentation + expansion이 원본보다 더 많은 transcript를 포착하는 것이 기대됨.
- **읽는법**: X축=방법(Original/Reseg), Y축=reads/cell. Boxplot의 중앙선=중앙값, 박스=IQR. 재분할의 박스가 원본보다 위에 있으면 효율 개선.
- **비정상 신호**: 재분할 >> 원본 (2배 이상) → 과도한 확장, misassignment 의심

#### 시각화 2: Efficiency - ST/SC Ratio (`efficiency_expression_ratio.png`) → Fig. 2c
- [ ] **피크 위치**: ratio 1.0 이상 (Xenium이 scRNAseq보다 효율적)
- [ ] **대부분 유전자**: ratio 1.2-1.5 범위
- [ ] **빨간 점선 (ratio=1)**: 오른쪽에 더 많은 유전자
- **왜 정상인가**: 논문 "detection efficiency was found to be between 1.2 and 1.5 times higher than that of scRNA-seq (Chromium v2)" (p.815). Xenium의 ISH 기반 탐지가 scRNAseq의 PCR 기반보다 특정 유전자에서 효율적.
- **읽는법**: X축=SRT/SC ratio (log2 또는 linear), Y축=빈도. ratio>1이면 Xenium이 더 효율적. 특정 유전자(highly expressed)에서 ratio가 높은 것이 정상.
- **비정상 신호**: 피크 < 0.5 → 프로브 문제 또는 reference 불일치

#### 시각화 3: Specificity - NMP (`specificity_nmp_per_gene_boxplot.png`) → Fig. 2d
- [ ] **중앙값 > 0.8**: 양호한 segmentation 품질
- [ ] **이상치 < 0.6**: 소수의 유전자만
- [ ] **Reseg NMP ≈ Original NMP**: 재분할이 NMP를 악화시키지 않음
- **왜 정상인가**: 논문 "all the different SRT technologies presented a mean high specificity (NCP > 0.8)" (p.816-817). NMP(Negative Marker Purity)는 세포 유형 특이적이지 않은 유전자가 해당 세포에서 낮게 발현되는 정도를 측정. 높을수록 segmentation이 정확.
- **읽는법**: Y축=NMP score (0-1). 높은 값=해당 세포 유형에서 발현되지 않아야 할 유전자가 실제로 낮게 발현(좋은 segmentation). 여러 방법 비교 시 높은 NMP=더 정확한 segmentation.
- **비정상 신호**: NMP < 0.6 → 심각한 transcript 누출, 확장 과도

#### 시각화 4: Diffusion - Complementary CDF (`diffusion_complementary_cdf_comparison.png`) → Fig. 2f
- [ ] **Reseg 곡선이 Nuclei보다 약간 위**: 정상 (cytoplasmic reads 포함으로 약간 퍼짐)
- [ ] **10 µm에서 P(d>x)**: 0.2 미만이 정상
- [ ] **두 곡선의 간격**: 적당 (너무 크지 않음)
- **왜 정상인가**: 논문 Fig. 2f의 complementary CDF. Nuclei-only는 핵 근처 reads만 포함하므로 CDF가 빠르게 0에 수렴. Resegmented는 확장된 영역 reads 포함으로 약간 높은 CDF. "MERFISH and ISS methods displayed a higher concentration of reads near the cellular centroid" (p.817).
- **읽는법**: X축=centroid 거리, Y축=P(d≥x) (거리 이상인 reads 비율). 곡선이 빨리 0에 가까워질수록 reads가 centroid 근처에 집중. 여러 방법의 곡선 비교로 segmentation 품질 비교 가능.
- **비정상 신호**: Reseg 곡선이 크게 위 → 과도한 확장으로 먼 reads 포함

#### 시각화 5: Diffusion - Assigned Reads (`diffusion_assigned_reads_barplot.png`) → ExtData Fig. 4b
- [ ] **Nuclear capture rate**: In Nucleus 비율이 합리적 범위 (10-30%)
- [ ] **Reseg ≤ Original**: 확장으로 핵 외부 reads 추가 → 핵 비율 상대적 감소가 정상
- [ ] **Stacked bar 합계 = 1.0**: In Nucleus + Outside Nucleus
- **왜 정상인가**: 핵 내부 transcript는 전체의 일부 (ISH 기반 기술은 세포질 reads도 다수). 확장(expansion) 시 핵 외부 reads가 추가되므로 핵 비율은 상대적으로 감소.
- **읽는법**: Stacked barplot. 청록(In Nucleus) + 회색(Outside Nucleus) = 1.0. Resegmented는 `in_cell`(Cellpose nuclear mask), Original은 `overlaps_nucleus`(Xenium native) 사용.
- **주의**: 이 metric은 "cell assignment rate"가 아닌 "nuclear capture rate". 논문의 76.8%는 세포 할당률이며, 본 barplot은 핵 overlap 비율을 측정하므로 수치가 다름.
- **비정상 신호**: In Nucleus > 60% → 비정상적으로 높은 핵 포착 (확장이 거의 없거나 데이터 이상)

### 11.3 논문 원문 인용 (시각화 관련)

| 시각화 | 논문 원문 인용 | 페이지 |
|---|---|---|
| Efficiency (Fig. 2b) | "resulting in very different numbers of reads per cell, with CosMx yielding the highest number" | 816 |
| ST/SC ratio (Fig. 2c) | "detection efficiency was found to be between 1.2 and 1.5 times higher than that of scRNA-seq (Chromium v2), depending on the metric and region analyzed" | 815 |
| NMP/NCP (Fig. 2d) | "all the different SRT technologies presented a mean high specificity (NCP > 0.8), with HS-ISS and Molecular Cartography being the most specific technologies" | 816-817 |
| Xenium NMP 특이성 | "The specificity of Xenium was slightly lower than that of other commercial platforms, but was consistently higher than that of CosMx, which presented the lowest values" | 817 |
| Diffusion (Fig. 2f) | "MERFISH and ISS methods displayed a higher concentration of reads near the cellular centroid, whereas ISH-based commercial platforms (CosMx, Molecular Cartography and MERSCOPE) had reads positioned farther away" | 817 |
| Turnover distance | "transcripts located more than 10.71 µm, on average, from the cell centroid exhibited a higher gene expression correlation with domain-specific background signatures" | 817 |
| Assigned reads | "Using Xenium's default segmentation, an average of 186.6 reads per cell was observed throughout the datasets, with 76.8% of reads being assigned to cells" | 813 |
| 확장 영향 | "Cell mask expansion after segmentation can negatively impact the characterization of cell populations by misassigning reads to neighboring cells" | 821 |
