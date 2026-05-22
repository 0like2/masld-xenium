# Step 1: Dataset Exploration - Final Review

## 1. Overview

### 1.1 What This Step Does
Step 1은 Xenium 데이터셋의 **기본 특성을 탐색하고 세포 유형을 식별**하는 단계이다. 총 세포 수, 유전자 발현 분포, transcript 분산(dispersion) 분석, Leiden 클러스터링, marker gene 탐지, 공간적 이웃 분석(neighborhood analysis)을 수행한다.

### 1.2 Paper Context
논문의 "Xenium datasets offer high-quality tissue population data" 섹션 (p.813-814)에서 데이터셋 특성을 탐색한 내용이 이 단계에 해당한다. Cell-by-gene matrix의 log-transformation, normalization, 40 principal components 기반 neighborhood graph, Leiden clustering으로 세포 유형을 식별하였다.

### 1.3 논문 Figure 매칭
| 논문 Figure | 설명 | Step 1 구현 |
|------------|------|------------|
| **Fig. 1b** | 데이터셋별 세포 수, 유전자 수, reads/cell 요약 | `calculate_general_stats()` |
| **Fig. 1c** | UMAP - cell type 색상화 (mouse brain) | `perform_clustering_and_annotation()` → UMAP plot |
| **Fig. 1d** | Spatial map - cell type 분포 | `perform_clustering_and_annotation()` → spatial scatter |
| **Fig. 1f** | Distance to centroid boxplot (nuclear/cyto genes) | `calculate_transcript_dispersion()` |
| **Fig. 2b** | Transcripts/cell, Genes/cell violin plot | `calculate_general_stats()` |
| **Fig. 2f** | Cumulative proportion of reads by distance to centroid | `calculate_transcript_dispersion()` → ECDF |
| **Extended Data Fig. 1c** | Violin: transcripts/cell, genes/cell 비교 | `calculate_general_stats()` |
| **Extended Data Fig. 1f** | Density plot: distance to centroid | `calculate_transcript_dispersion()` |

---

## 2. Pipeline Process Flow

```
Step 0 Output
├── {sample_tag}.h5ad                    ← AnnData (cells × genes)
├── {sample_tag}_transcripts.parquet     ← Transcript sidecar
         ↓
    [Step 1: Dataset Exploration]
         ↓
    ┌──── 1-1. General Statistics ────┐
    │  총 세포/유전자 수, sparsity,     │
    │  reads/cell 분포, 히트맵          │
    └─────────────────────────────────┘
         ↓
    ┌──── 1-2. Transcript Dispersion ─┐
    │  Distance histogram + KDE        │
    │  ECDF (전체 + 유전자별)           │
    │  Violin (상위 유전자)             │
    │  KS Test (유전자 쌍)             │
    └─────────────────────────────────┘
         ↓
    ┌──── 1-3. Clustering ────────────┐
    │  HVG 선택 → PCA → Neighbors     │
    │  Leiden (다중 해상도)             │
    │  Marker gene detection           │
    │  UMAP + Spatial scatter          │
    └─────────────────────────────────┘
         ↓
    ┌──── 1-4. Neighborhood Analysis ─┐
    │  Spatial neighbors graph         │
    │  Neighborhood enrichment         │
    │  Centrality scores               │
    └─────────────────────────────────┘
         ↓
├── {sample_tag}_step1_exploration.h5ad  ← 분석 결과 통합 AnnData
├── step1_exploration/                   ← 통계/플롯 출력 디렉토리
```

---

## 3. Sub-step 상세 분석

### 3.1 General Statistics (`calculate_general_stats()`, lines 90-189)

**목적**: 데이터셋의 기본 품질 지표를 계산하고 시각화

**계산 항목**:
| 지표 | 계산법 | 의미 |
|------|--------|------|
| Total cells | `adata.n_obs` | 전체 세포 수 |
| Total genes | `adata.n_vars` | 프로파일된 유전자 수 |
| Mean counts/cell | `adata.obs['n_counts'].mean()` | 평균 reads per cell |
| Median counts/cell | `adata.obs['n_counts'].median()` | 중앙값 reads per cell |
| Mean genes/cell | `adata.obs['n_genes'].mean()` | 평균 발현 유전자 수 |
| Sparsity | `1 - (adata.X > 0).sum() / total_entries` | 0인 엔트리 비율 |

**출력 파일**:
- `{sample_tag}_step1_stats.txt` - 텍스트 요약
- `{sample_tag}_step1_stats.csv` - CSV 형식 통계
- `{sample_tag}_step1_stats_heatmap.png` - 유전자별 통계 히트맵

### 3.2 Transcript Dispersion Analysis (`calculate_transcript_dispersion()`, lines 258-367)

**목적**: 각 transcript의 세포 centroid로부터의 거리를 계산하여 subcellular 분포 패턴 분석

**핵심 개념 - Dispersion (분산)**:
```
dispersion = distance(transcript_xy, cell_centroid_xy)
           = sqrt((x_transcript - x_centroid)² + (y_transcript - y_centroid)²)
```

**주의**: `xb.calculating.dispersion()` 함수의 "dispersion"은 **centroid distance**를 의미한다. Boundary distance와 다름.

**처리 과정**:
1. Transcript sidecar에서 cell_id > 0인 할당된 transcript 추출
2. cells.csv에서 세포별 centroid 좌표 매칭
3. 유클리드 거리 계산 (pixel 단위 → µm 변환 가능)
4. 유전자별 거리 분포 통계 (mean, median, std)

**시각화**:
| Plot | 설명 | 분석법 |
|------|------|--------|
| `dispersion_dist.png` | 전체 거리 히스토그램 + KDE | 분포 형태 확인: 정규/이봉/우치우 편향 |
| `dispersion_ecdf.png` | ECDF (경험적 누적분포) | 거리 임계값에서의 transcript 비율 확인 |
| `dispersion_ecdf_genes.png` | 상위 10 유전자별 ECDF | nuclear vs cytoplasmic gene 구분 |
| `dispersion_violin.png` | 상위 20 유전자 거리 violin | 유전자 간 공간 분포 차이 |
| `dispersion_metrics.txt` | 통계 요약 | 정량적 비교용 |

**분석법**:
- ECDF에서 급격히 상승하는 유전자 = nuclear-enriched (핵 가까이 위치)
- ECDF가 완만한 유전자 = diffuse/cytoplasmic (넓게 분포)
- 논문 Fig. 1f에서 이 패턴을 확인: Opalin(nuclear) vs MAG(cytoplasmic)

### 3.3 KS Test (`run_ks_tests()`, lines 369-453)

**목적**: 유전자 쌍 간 거리 분포의 통계적 차이를 Kolmogorov-Smirnov 검정으로 측정

```python
# 유전자 A와 B의 거리 분포가 같은 분포에서 왔는지 검정
ks_stat, p_value = scipy.stats.ks_2samp(distances_geneA, distances_geneB)
```

**출력**:
- `ks_pvalues.csv` - 유전자 쌍별 p-value 매트릭스
- `ks_statistics.csv` - KS 통계량 매트릭스
- `ks_heatmap.png` - KS 통계량 히트맵

**분석법**: KS 통계량이 높은 유전자 쌍 = 공간적으로 서로 다른 분포 (nuclear vs cytoplasmic)

### 3.4 Clustering & Marker Annotation (`perform_clustering_and_annotation()`, lines 458-588)

**처리 과정**:
```
normalize_total(target_sum=1e4)  ← library-size normalization
        ↓
    log1p()                      ← log transformation
        ↓
highly_variable_genes()          ← 고변동 유전자 선택
        ↓
    PCA (n_comps=40)             ← 차원 축소
        ↓
neighbors(n_neighbors=15)        ← k-NN 그래프 구축
        ↓
Leiden(resolutions=[0.5,0.8,1.0,1.5])  ← 클러스터링 (다중 해상도)
        ↓
rank_genes_groups(method='wilcoxon')   ← marker gene 탐지
        ↓
    UMAP                         ← 2D 임베딩
```

**설정값 의미**:
| 파라미터 | 기본값 | 의미 |
|----------|--------|------|
| `resolutions` | [0.5, 0.8, 1.0, 1.5] | Leiden 해상도: 높을수록 더 많은 클러스터 |
| `primary_resolution` | 1.0 | marker/시각화에 사용되는 주 해상도 |
| `n_neighbors` | 15 | k-NN 그래프의 이웃 수 |

**시각화**:
| Plot | 설명 | 분석법 |
|------|------|--------|
| `hvg.png` | Highly variable genes 선택 결과 | 분산-평균 관계에서 선택된 유전자 확인 |
| `markers_res{r}.csv` | 클러스터별 marker gene 목록 | 상위 DEG로 세포 유형 추정 |
| `markers_dotplot.png` | 클러스터×marker dotplot | 점 크기=발현 비율, 색=발현 강도 |
| `markers_heatmap.png` | Marker gene heatmap | 클러스터별 유전자 발현 패턴 |
| `umap_res{r}.png` | UMAP 임베딩 (클러스터 색상) | 세포 유형 분리 확인 |
| `spatial_res{r}.png` | 공간 scatter (클러스터 색상) | 조직 내 세포 유형 배치 |

### 3.5 Neighborhood Analysis (`analyze_neighborhoods()`, lines 593-667)

**목적**: Squidpy를 사용한 공간적 이웃 구조 분석

**처리 과정**:
1. **Spatial Neighbors Graph**: Delaunay triangulation 기반 공간 이웃 그래프
2. **Neighborhood Enrichment**: 특정 세포 유형 쌍이 예상보다 더 자주/덜 자주 이웃하는지 검정
3. **Centrality Scores**: 그래프 이론 기반 중심성 (betweenness, closeness)

```python
# Squidpy 공간 이웃 분석
sq.gr.spatial_neighbors(adata, coord_type='generic', radius=100.0)
sq.gr.nhood_enrichment(adata, cluster_key=f'leiden_{primary_res}')
sq.gr.centrality_scores(adata, cluster_key=f'leiden_{primary_res}')
```

**시각화**:
| Plot | 설명 | 분석법 |
|------|------|--------|
| `grad_enrichment.png` | Neighborhood enrichment z-score heatmap | 양수=공동위치, 음수=배타적 |
| `centrality.png` | Centrality 점수 (betweenness, closeness) | 높은 centrality = hub 역할 세포 유형 |

---

## 4. Notebook vs Pipeline 구현 비교

### 4.1 원본 Notebook들
| 노트북 | Pipeline 함수 | 상태 |
|--------|--------------|------|
| `1_1_Statistics_all_samples_using_txsim.ipynb` | `calculate_general_stats()` | 구현됨 |
| `1_2_Celltype_identification_PART1_msbrain_nuclei.ipynb` | `perform_clustering_and_annotation()` | 구현됨 |
| `1_3_Read_specific_dispersion_metrics_Xenium.ipynb` | `calculate_transcript_dispersion()` | 구현됨 |
| `1_5_Celltype_identification_mscoronal_multisection-EXPANDED-cells.ipynb` | 부분 구현 (expanded cells 미지원) | 부분 |
| `1_6_Batch_preprocessing_real_Xenium_datasets.ipynb` | `perform_clustering_and_annotation()` | 통합됨 |
| `ADDITIONAL_1_7_Celltype_architecture_structure_scores_msbrain_PART2.ipynb` | `analyze_neighborhoods()` | 구현됨 |

### 4.2 상세 비교

| 기능 | Notebook | Pipeline | 차이점 |
|------|----------|----------|--------|
| 통계 계산 | txsim 라이브러리 사용 | 직접 구현 | txsim 의존성 제거 |
| Dispersion | `xb.calculating.dispersion()` | 직접 유클리드 거리 계산 | 동일 결과 |
| KS Test | 유전자 쌍별 수동 | `scipy.stats.ks_2samp` 루프 | 동일 |
| Clustering | 단일 해상도 | 다중 해상도 지원 | 개선 |
| Marker detection | Wilcoxon | Wilcoxon | 동일 |
| Neighborhood | Squidpy | Squidpy | 동일 |
| 해상도 설정 | 하드코딩 | config.yaml | 개선 |

---

## 5. 전체 시각화 목록 및 분석법

| # | 시각화 파일 | 유형 | 의미 | 분석법 |
|---|-----------|------|------|--------|
| 1 | `stats.txt/csv` | 통계 | 데이터셋 기본 품질 | 세포/유전자 수 확인, sparsity 평가 |
| 2 | `stats_heatmap.png` | Heatmap | 유전자별 통계 | 발현 패턴 개요 확인 |
| 3 | `dispersion_dist.png` | Histogram+KDE | Transcript-centroid 거리 분포 | 분포 형태로 segmentation 품질 평가 |
| 4 | `dispersion_ecdf.png` | ECDF | 누적 거리 분포 | 50% 지점 = 중앙 transcript 거리 |
| 5 | `dispersion_ecdf_genes.png` | ECDF/gene | 유전자별 거리 분포 | Nuclear vs cytoplasmic gene 구분 |
| 6 | `dispersion_violin.png` | Violin | 상위 유전자 거리 | 유전자 간 분포 폭 비교 |
| 7 | `dispersion_metrics.txt` | 텍스트 | Dispersion 요약 통계 | 정량 지표 확인 |
| 8 | `ks_pvalues.csv` | 매트릭스 | 유전자 쌍 KS p-value | p<0.05: 유의미한 분포 차이 |
| 9 | `ks_statistics.csv` | 매트릭스 | KS 통계량 | 높을수록 큰 분포 차이 |
| 10 | `ks_heatmap.png` | Heatmap | KS 통계량 시각화 | 밝은 영역 = 큰 차이 유전자 쌍 |
| 11 | `hvg.png` | Scatter | HVG 선택 결과 | 선택된 유전자가 합리적인지 확인 |
| 12 | `markers_res{r}.csv` | 테이블 | 클러스터별 marker gene | 상위 DEG로 세포 유형 추정 |
| 13 | `markers_dotplot.png` | Dotplot | Marker×Cluster | 특이적 marker 패턴 확인 |
| 14 | `markers_heatmap.png` | Heatmap | Marker 발현 | 클러스터 간 발현 차이 |
| 15 | `umap_res{r}.png` | UMAP | 세포 클러스터링 | 클러스터 분리도 확인 |
| 16 | `spatial_res{r}.png` | Spatial scatter | 공간 클러스터 배치 | 조직 구조와의 일치 확인 |
| 17 | `umap_res{r}_noMTRNR.png` | UMAP | MTRNR 제외 클러스터링 | MTRNR (미토콘드리아 rRNA) 유전자 제거 후 재클러스터링. MTRNR이 dominant하면 세포 유형 분리를 방해할 수 있으므로 제외 후 비교 |
| 18 | `spatial_res{r}_noMTRNR.png` | Spatial scatter | MTRNR 제외 공간 배치 | MTRNR 제외 후 공간 클러스터 분포 변화 확인 |
| 19 | `markers_dotplot_noMTRNR.png` | Dotplot | MTRNR 제외 Marker | MTRNR 제외 후 marker gene 변화. 원본과 비교하여 MTRNR 영향도 평가 |
| 20 | `grad_enrichment.png` | Heatmap | Neighborhood enrichment | 세포 유형 공동위치 패턴 |
| 21 | `centrality.png` | Bar/Scatter | Centrality 점수 | Hub 세포 유형 식별 |

---

## 6. Input / Output 상세

### 6.1 Input
| 파일 | 형식 | 출처 |
|------|------|------|
| `{sample_tag}.h5ad` | AnnData | Step 0 |
| `{sample_tag}_transcripts.parquet` | Parquet | Step 0 |
| 또는 `transcripts.csv` | CSV | 원본 Xenium |

### 6.2 Output
| 파일 | 형식 | 설명 |
|------|------|------|
| `{sample_tag}_step1_exploration.h5ad` | AnnData | 클러스터링 + 이웃 분석 결과 포함 |
| `step1_exploration/` 내 통계/플롯 | CSV/PNG/TXT | 21개 이상의 분석 결과 파일 |

### 6.3 출력 AnnData 구조 변경사항
```
adata.obs 추가 컬럼:
  - leiden_0.5, leiden_0.8, leiden_1.0, leiden_1.5  ← 다중 해상도 클러스터
  - n_counts, n_genes                               ← QC 메트릭

adata.obsm 추가:
  - X_pca   ← PCA 좌표
  - X_umap  ← UMAP 좌표

adata.uns 추가:
  - rank_genes_groups  ← marker gene 분석 결과
  - nhood_enrichment   ← neighborhood enrichment 결과
  - centrality_scores  ← centrality 결과
```

---

## 7. 관련 설정값 정리

```yaml
exploration:
  neighbor_radius: 100.0           # Squidpy 공간 이웃 반경 (µm)
  run_transcript_dispersion: true  # Dispersion 분석 실행 여부
  resolutions: [0.5, 0.8, 1.0, 1.5]  # Leiden 해상도 목록
  primary_resolution: 1.0         # 주 분석용 해상도
```

### 7.1 파라미터 조정 가이드

| 파라미터 | 현재값 | 범위/옵션 | 영향도 | 조정 가이드 |
|----------|--------|----------|--------|------------|
| `neighbor_radius` | `100.0` | 50-500 µm | MEDIUM | Squidpy 공간 이웃 탐색 반경. 작을수록 미세 구조, 클수록 대규모 패턴. 뇌 조직: 100-200 µm 권장 |
| `run_transcript_dispersion` | `true` | true/false | LOW | Dispersion 분석 실행 여부. 시간 절약 시 false 가능. 논문 Fig. 1f 재현에 필수 |
| `resolutions` | `[0.5, 0.8, 1.0, 1.5]` | 0.1-3.0 리스트 | **HIGH** | Leiden clustering 해상도 목록. 값이 클수록 더 많은 클러스터 생성. 0.5(대분류)~1.5(세분류) |
| `primary_resolution` | `1.0` | resolutions 목록 중 하나 | **HIGH** | Marker gene 탐지, UMAP 시각화에 사용되는 주 해상도. 세포 유형 수에 맞춰 조정 |

**튜닝 팁**:
- `primary_resolution`이 가장 중요. 예상 세포 유형 수에 맞춰 설정: 5-10개 유형=0.5, 10-20개=1.0, 20-30개=1.5
- `resolutions` 리스트에 여러 값을 포함하면 각 해상도에서의 결과를 비교하여 최적 값을 찾을 수 있음
- `neighbor_radius`: 조직 밀도에 따라 조정. 세포 밀도가 높은 cortex=100µm, 분산된 white matter=200µm

---

## 8. 핵심 로직의 의미

### 8.1 Leiden Resolution
- **0.5**: 대분류 (major cell types: neurons, glia, vascular)
- **0.8**: 중분류 (sub-types 일부 분리)
- **1.0**: 세분류 (논문 기본 사용 해상도)
- **1.5**: 과분류 가능성 (rare population 탐지용)

### 8.2 Neighborhood Enrichment Z-score
- **Z > 2**: 해당 세포 유형 쌍이 공간적으로 유의미하게 공동위치
- **Z < -2**: 해당 세포 유형 쌍이 공간적으로 배타적
- **|Z| < 2**: 유의미하지 않은 공간적 관계

### 8.3 KS Statistic
- **0에 가까움**: 두 유전자의 거리 분포가 유사 (같은 subcellular 구획)
- **1에 가까움**: 두 유전자의 거리 분포가 완전히 다름 (다른 구획)

---

## 10. 시각화 상세 분석 가이드

### 10.1 Dispersion Distribution (`dispersion_dist.png`) → 논문 Fig. 1f 관련

**논문 위치**: Fig. 1f (p.815) - "Distance of transcripts to their assigned cell centroid"

**무엇을 봐야 하는가**:
```
    Count
    ↑
    │  ██
    │  ████
    │  ██████
    │  ████████▓▓
    │  ████████████▓▓▓▓
    │  ████████████████████▓▓▓▓▓▓
    └──────────────────────────────→ Distance to centroid (µm)
       0    5   10   15   20   25

    확인 사항:
    ① 피크 위치 (mode): 5-8 µm이 정상 (핵 반경 근처)
    ② 오른쪽 꼬리 (tail): 긴 꼬리 = cytoplasmic reads 많음
    ③ 이봉 분포: 두 번째 피크 = segmentation 오류 또는 이웃 세포 혼입
```

**해석법**:
- **피크 ≈ 5-7 µm**: 대부분의 transcript가 핵 근처에 위치. 정상적인 Xenium 데이터
- **피크 > 10 µm**: 핵 segmentation이 세포 centroid와 크게 어긋남. 또는 확장이 과도함
- **이봉 분포 (bimodal)**: 첫 번째 피크 = nuclear reads, 두 번째 피크 = cytoplasmic reads 또는 misassigned reads
- **논문 결과**: "transcripts located more than 10.71 µm from centroid exhibited higher correlation with background" (p.817)
- **KDE 곡선**: 히스토그램 위의 부드러운 밀도 곡선으로 분포 형태를 더 명확히 파악

**좋은 결과 vs 나쁜 결과**:
- 좋음: 단봉(unimodal), 피크 5-7 µm, 20 µm 이상에서 급격히 감소
- 나쁨: 이봉, 피크 > 15 µm, 긴 오른쪽 꼬리 (50 µm+)

---

### 10.2 Dispersion ECDF (`dispersion_ecdf.png`) → 논문 Fig. 2f

**논문 위치**: Fig. 2f (p.816) - "Cumulative proportion of reads by distance to centroid"

**무엇을 봐야 하는가**:
```
    Cumulative
    Proportion
    1.0 ├─────────────────────────────────── ←── 100% reads
        │                          ╱──────
        │                    ╱────╱
    0.5 │               ╱───╱           ←── 50% reads at ~7 µm (좋음)
        │          ╱───╱
        │     ╱───╱
    0.0 ├───╱
        └──────────────────────────────→ Distance to centroid (µm)
           0     5    10    15    20

    핵심 판독점:
    ① 50% 도달 거리: 7 µm 이하 = 좋은 segmentation
    ② 90% 도달 거리: 15 µm 이하 = 적절한 확장
    ③ 곡선 기울기: 급격 = concentrated, 완만 = diffused
```

**해석법**:
- **ECDF 급상승 (steep)**: transcript가 centroid 근처에 집중 → 좋은 segmentation
- **ECDF 완만 (gradual)**: transcript가 넓게 퍼짐 → diffuse한 발현 또는 segmentation 문제
- **50% 지점이 5 µm 이하**: 매우 좋은 segmentation (핵 반경 내에 반이 위치)
- **50% 지점이 10 µm 이상**: segmentation 경계가 넓거나 misassignment 존재
- 논문에서는 플랫폼 간 비교에 사용: Xenium vs CosMx vs MERFISH

---

### 10.3 유전자별 Dispersion ECDF (`dispersion_ecdf_genes.png`) → 논문 Fig. 1f

**논문 위치**: Fig. 1f (p.815) 좌측 - boxplot과 함께 유전자별 centroid distance 비교

**무엇을 봐야 하는가**:
```
    Cumulative
    Proportion
    1.0 ├──────────────────────────────
        │    Opalin ╱─────── ← Nuclear gene (핵에 가까움)
        │          ╱  Mbp ╱────
    0.5 │    ╱────╱      ╱
        │   ╱      ╱────╱  MAG ╱── ← Cytoplasmic gene (핵에서 멀음)
        │  ╱  ╱───╱      ╱───╱
    0.0 ├─╱──╱───────────╱
        └────────────────────────→ Distance (µm)
           0    5   10   15   20

    ① 왼쪽으로 치우친 유전자 = Nuclear-enriched (핵 내 mRNA)
    ② 오른쪽으로 치우친 유전자 = Cytoplasmic/membrane (세포질 mRNA)
    ③ 유전자 간 간격 = subcellular compartment 차이의 크기
```

**해석법**:
- **Nuclear genes** (Opalin, Mog 등): ECDF가 빠르게 상승 → 핵 가까이 분포
- **Cytoplasmic genes** (MAG, PLP1 등): ECDF가 느리게 상승 → 넓게 분포
- 두 그룹 간의 ECDF 차이가 클수록 segmentation이 subcellular 정보를 잘 보존
- **KS test 결과와 연계**: ECDF 차이가 큰 유전자 쌍 = KS 통계량이 높음

---

### 10.4 Dispersion Violin (`dispersion_violin.png`)

**무엇을 봐야 하는가**:
```
    Distance (µm)
    30 ┤
       │    ╭─╮          ╭───╮
    20 ┤    │ │    ╭─╮   │   │
       │  ╭─│ │╮   │ │   │   │
    10 ┤  │ │━│ │  │━│   │ ━ │   ← Median 위치가 핵심
       │  │ │ │ │  │ │   │   │
     0 ┤  ╰─╯ ╰─╯  ╰─╯   ╰───╯
       └──Gene1──Gene2──Gene3──→

    ① Median (가로선) 위치: Nuclear gene은 낮음, Cyto gene은 높음
    ② Violin 폭: 넓음 = 많은 transcript, 좁음 = 적은 transcript
    ③ 분포 대칭성: 비대칭 = 특정 방향으로 편향된 분포
```

**해석법**:
- 상위 20 유전자의 거리 분포를 한눈에 비교
- Median이 5 µm 이하인 유전자 = nuclear-enriched
- Median이 10 µm 이상인 유전자 = cytoplasmic/membrane
- Violin이 매우 넓은 유전자 = 다양한 세포 유형에서 다른 subcellular 분포

---

### 10.5 KS Heatmap (`ks_heatmap.png`)

**무엇을 봐야 하는가**:
```
         Gene1  Gene2  Gene3  Gene4
    Gene1  0.0   0.15   0.45   0.52
    Gene2  0.15  0.0    0.38   0.48
    Gene3  0.45  0.38   0.0    0.12   ← 같은 구획의 유전자 = 낮은 KS
    Gene4  0.52  0.48   0.12   0.0

    색상: 밝음(노랑) = 높은 KS 통계량 = 큰 분포 차이
          어두움(보라) = 낮은 KS 통계량 = 유사한 분포

    ① 밝은 블록: nuclear vs cytoplasmic 유전자 쌍
    ② 어두운 블록: 같은 subcellular 구획의 유전자 쌍
    ③ 대각선은 항상 0 (자기 자신과의 비교)
```

**해석법**:
- KS 통계량 > 0.3: 유의미한 분포 차이 (다른 subcellular 구획)
- KS 통계량 < 0.1: 유사한 분포 (같은 구획)
- 블록 패턴이 보이면 유전자 그룹이 distinct subcellular compartment에 존재

---

### 10.6 HVG 선택 (`hvg.png`)

**무엇을 봐야 하는가**:
```
    Dispersion
    (normalized)
    ↑         ● ●  ← Highly variable genes (선택됨)
    │       ● ●  ●
    │     ● ●●● ● ●
    │   ○○○○●●●○○○
    │  ○○○○○○○○○○     ← Non-HVG (선택 안 됨)
    │ ○○○○○○○○○
    └─────────────────→ Mean expression

    ● = HVG (선택됨), ○ = Non-HVG
    ① 선택된 유전자 수가 합리적인지 (보통 전체의 10-30%)
    ② Cell type marker gene들이 HVG에 포함되는지
    ③ 너무 높은 mean expression의 유전자(housekeeping)가 제외되는지
```

---

### 10.7 Marker Dotplot (`markers_dotplot.png`) → 논문 Fig. 1c 관련

**무엇을 봐야 하는가**:
```
    Cluster  0    1    2    3    4    5
    Gene1    ●                        ← Gene1은 Cluster 0에 특이적
    Gene2         ●●                  ← Gene2는 Cluster 1에 특이적
    Gene3              ●   ●         ← Gene3은 Cluster 2,3에서 발현
    Gene4                   ●●●      ← Gene4는 Cluster 4에서 강하게 발현
    Gene5                        ○   ← Gene5는 약하게 발현

    점 크기 = 발현 비율 (해당 클러스터에서 양성인 세포 %)
    점 색상 = 평균 발현 강도 (진할수록 높은 발현)

    좋은 marker:
    ① 하나의 클러스터에서만 큰+진한 점 (특이적 발현)
    ② 다른 클러스터에서는 작거나 없는 점 (비특이적 발현 없음)
```

**해석법**:
- **이상적인 marker**: 한 클러스터에서 80%+ 발현, 나머지에서 10% 미만
- **비특이적 marker**: 여러 클러스터에서 비슷한 크기의 점 → 해당 유전자로는 클러스터 구분 불가
- Xenium 패널에서는 curated marker 유전자들이므로, 대부분 특이적 패턴을 보여야 함
- 알츠하이머 뇌 데이터에서 기대되는 marker: GFAP(astrocyte), MBP(oligo), SNAP25(neuron)

---

### 10.8 UMAP (`umap_res{r}.png`) → 논문 Fig. 1c

**논문 위치**: Fig. 1c (p.815) - "UMAP showing cell types"

**무엇을 봐야 하는가**:
```
    ┌───────────────────────────────┐
    │     ○○○           ●●●        │
    │    ○○○○○         ●●●●●       │  ← 잘 분리된 클러스터
    │     ○○○           ●●●        │
    │                               │
    │  ▲▲▲     □□□                  │
    │ ▲▲▲▲▲   □□□□□     ◇◇◇       │  ← 각 기호 = 다른 클러스터
    │  ▲▲▲     □□□      ◇◇◇       │
    └───────────────────────────────┘

    확인 사항:
    ① 클러스터 간 분리도: 잘 떨어져 있으면 좋음
    ② 클러스터 내 밀도: 촘촘하면 동질적 (homogeneous)
    ③ 브릿지 세포: 클러스터 사이의 세포 = transitional state
    ④ 해상도별 비교: 0.5→대분류, 1.0→세분류, 1.5→과분류
```

**해석법**:
- **잘 분리된 클러스터**: 세포 유형이 명확히 구분됨 → 데이터 품질 양호
- **혼재된 클러스터**: 세포 유형 구분이 어려움 → segmentation 품질 문제 또는 해상도 부적절
- **Resolution 0.5**: 5-10개 클러스터 (major cell types)
- **Resolution 1.0**: 10-20개 클러스터 (sub-types 포함)
- **Resolution 1.5**: 20-30개 클러스터 (rare population 탐지, 과분류 위험)

---

### 10.9 Spatial Scatter (`spatial_res{r}.png`) → 논문 Fig. 1d

**논문 위치**: Fig. 1d (p.815) - "Spatial distribution of cell types"

**무엇을 봐야 하는가**:
```
    ┌───────────────────────────────┐
    │ ●●●●●●    ○○○○○○○    ▲▲▲▲   │  ← Layer structure
    │ ●●●●●     ○○○○○○     ▲▲▲    │     (cortex layers)
    │                               │
    │ □□□□□□□□□□□□□□□□□□□          │  ← White matter
    │ □□□□□□□□□□□□□□□□□□□          │
    │                               │
    │ ◇◇◇◇◇    ●●●●●●             │  ← Mixed regions
    └───────────────────────────────┘

    ① 조직 구조와 클러스터 배치의 일치 확인
    ② 특정 영역에 편중된 세포 유형 식별
    ③ 무작위 분포 = 품질 문제 또는 해상도 부적절
```

**해석법**:
- **뇌 조직**: Cortical layers (I-VI), white matter, meninges 등 구조가 보여야 함
- **Neuron (SNAP25+)**: Gray matter에 집중
- **Oligodendrocyte (MBP+)**: White matter에 집중
- **Astrocyte (GFAP+)**: 넓게 분포하되 gray/white matter 경계에 많음
- **Microglia (AIF1+)**: 전체적으로 산재
- **공간적 무작위 분포**: 클러스터링이 세포 유형을 제대로 포착하지 못함

---

### 10.10 Neighborhood Enrichment (`grad_enrichment.png`)

**무엇을 봐야 하는가**:
```
              Cluster 0  1    2    3    4
    Cluster 0   --     +3.2 -1.5  +5.1  0.2
    Cluster 1  +3.2    --   +2.8  -0.3  +1.5
    Cluster 2  -1.5   +2.8  --    -2.1  +0.8
    Cluster 3  +5.1   -0.3 -2.1   --    +4.2
    Cluster 4   0.2   +1.5  +0.8  +4.2  --

    색상: 빨강/양수 = 공동위치 (co-localization)
          파랑/음수 = 배타적 (spatial exclusion)
          흰색/0 = 무작위

    ① 강한 양수 (Z > 3): 항상 함께 있는 세포 유형
    ② 강한 음수 (Z < -3): 절대 함께 있지 않는 세포 유형
    ③ 대각선: 같은 유형끼리의 self-enrichment (보통 양수)
```

**해석법 (뇌 조직)**:
- **Neuron-Neuron**: 양수 (같은 layer에 모여 있음)
- **Neuron-Oligodendrocyte**: 음수 (gray vs white matter 분리)
- **Astrocyte-모든 유형**: 약한 양수 (넓게 분포)
- **Microglia-모든 유형**: ~0 (무작위 분포)
- 알츠하이머 뇌에서는 Microglia-Amyloid plaque 근처에 enrichment 가능

---

### 10.11 Centrality Scores (`centrality.png`)

**무엇을 봐야 하는가**:
```
    Betweenness     Closeness     Degree
    Centrality      Centrality    Centrality
    ↑               ↑             ↑
    │  ██            │  ████       │  ████
    │  ████          │  ████       │  ████████
    │  ██████        │  ██████     │  ████████
    └──────→         └──────→      └──────→
    Cell Types       Cell Types    Cell Types

    ① Betweenness 높음: 다른 세포 유형 간의 "다리" 역할
    ② Closeness 높음: 조직의 중심부에 위치
    ③ Degree 높음: 많은 이웃과 접촉
```

**해석법**:
- **높은 Betweenness**: 해당 세포 유형이 서로 다른 tissue domain을 연결 (예: astrocytes가 gray/white matter 경계에서)
- **높은 Closeness**: 조직 중심에 위치하는 세포 유형
- **높은 Degree**: 밀집 영역에 있는 세포 유형

---

## 9. Summary

Step 1은 데이터셋의 **탐색적 분석(EDA)** 단계로, QC 지표 계산, transcript 공간 분석, 세포 클러스터링, 공간 이웃 분석을 수행한다.

| 평가 항목 | 상태 |
|----------|------|
| General statistics | 완전 구현 |
| Transcript dispersion | 완전 구현 |
| KS test | 완전 구현 |
| Clustering (multi-resolution) | 완전 구현 |
| Marker detection | 완전 구현 |
| Neighborhood analysis | 완전 구현 |
| Expanded cell analysis (1_5) | 미구현 (Step 5에서 처리) |

**구현 완성도**: **HIGH** - 핵심 분석 기능이 모두 구현됨

---

## 11. 논문 Figure 직접 대응 및 시각화 정상 판별 종합 가이드

### 11.1 파이프라인 출력 → 논문 Figure 매핑 종합표

| # | 파이프라인 출력 파일 | 논문 Figure | 논문 페이지 | 논문 원문 설명 |
|---|---|---|---|---|
| 1 | `dispersion_dist.png` | **Fig. 1f** (p.815) 좌측 boxplot | 815 | Distance to centroid 분포 - nuclear vs cytoplasmic mRNAs |
| 2 | `dispersion_ecdf.png` | **Fig. 2f** (p.816) | 816 | "Cumulative proportion of reads by distance from the cell centroid" |
| 3 | `dispersion_ecdf_genes.png` | **Fig. 1f** (p.815) | 815 | 유전자별 centroid 거리 분포: nuclear gene(Opalin)은 가깝고, cytoplasmic gene(MAG)은 멀다 |
| 4 | `dispersion_violin.png` | **Fig. 1f** 관련 | 815 | 상위 유전자의 거리 분포 비교 |
| 5 | `ks_heatmap.png` | 직접 대응 없음 | - | 유전자 쌍 간 거리 분포 차이의 정량화 (KS 통계량) |
| 6 | `hvg.png` | 직접 대응 없음 | - | HVG 선택은 표준 Scanpy 전처리 단계 |
| 7 | `markers_dotplot.png` | **Fig. 1c** 관련 | 815 | 클러스터별 marker gene 발현 패턴 |
| 8 | `umap_res{r}.png` | **Fig. 1c** (p.815) | 815 | "UMAP of cells colored by cell type" |
| 9 | `spatial_res{r}.png` | **Fig. 1d** (p.815) | 815 | "Spatial map of cell types; replicate 1 is shown" |
| 10 | `grad_enrichment.png` | 직접 대응 없음 (노트북 1_7) | - | Squidpy neighborhood enrichment - 세포 유형 간 공간 관계 |
| 11 | `centrality.png` | 직접 대응 없음 (노트북 1_7) | - | Graph centrality - hub 세포 유형 식별 |

### 11.2 정상 결과 판별 체크리스트

#### 시각화 1: `dispersion_dist.png` (Transcript-Centroid 거리 분포)
- [ ] **피크 위치**: 5-8 µm에 위치 (핵 반경 근처)
- [ ] **분포 형태**: 단봉(unimodal) right-skewed
- [ ] **20 µm 이후**: 급격히 감소
- [ ] **이봉(bimodal) 아님**: 두 번째 피크가 없어야 정상
- **왜 정상인가**: 논문 Fig. 1f에서 대부분의 transcript가 핵 근처(5-7 µm)에 위치. "transcripts located more than 10.71 µm from centroid exhibited higher correlation with background" (p.817). 생물학적으로 대부분의 mRNA는 핵에서 전사 후 핵 주변에 존재.
- **읽는법**: X축=centroid까지 거리(µm), Y축=transcript 수. 피크가 왼쪽(작은 거리)에 있으면 transcript가 핵 근처에 잘 위치한 것.
- **비정상 신호**: 피크 > 15 µm (segmentation 불량), 이봉 분포 (misassignment)

#### 시각화 2: `dispersion_ecdf.png` (누적 거리 분포)
- [ ] **50% 도달 거리**: 7 µm 이하
- [ ] **90% 도달 거리**: 15 µm 이하
- [ ] **곡선 형태**: S자형 (sigmoid)
- **왜 정상인가**: 논문 Fig. 2f에서 Xenium의 ECDF가 다른 플랫폼(CosMx, MERSCOPE)보다 centroid 근처에 집중됨. Xenium의 ISH 기반 탐지가 읽기의 위치 정밀도가 높기 때문.
- **읽는법**: X축=거리, Y축=누적 비율. 곡선이 빨리 1.0에 가까워질수록 reads가 centroid에 가까운 것. 가파른 상승=좋은 결과.
- **비정상 신호**: 50% 지점이 10 µm 이상 → segmentation 경계 과도

#### 시각화 3: `dispersion_ecdf_genes.png` (유전자별 ECDF)
- [ ] **Nuclear gene** (Opalin, Mog 등): ECDF가 빠르게 상승 (5 µm에서 80%+)
- [ ] **Cytoplasmic gene** (MAG, PLP1 등): ECDF가 느리게 상승
- [ ] **두 그룹 간 ECDF 간격**: 명확히 구분됨
- **왜 정상인가**: 논문 Fig. 1f에서 "Most nuclear genes" (Opalin, Sox10 등)와 "Least nuclear genes" (MAG, CAPN2 등)의 명확한 거리 차이 확인. 이는 subcellular localization의 생물학적 차이를 반영.
- **읽는법**: 각 곡선은 하나의 유전자. 왼쪽(작은 거리)에서 빠르게 상승=핵 내 유전자, 오른쪽에서 천천히 상승=세포질 유전자. 두 그룹이 분리되면 Xenium이 subcellular resolution을 제공한다는 증거.
- **비정상 신호**: 모든 유전자의 ECDF가 겹침 → subcellular resolution 부족

#### 시각화 4: `umap_res{r}.png` (UMAP 클러스터링)
- [ ] **클러스터 분리**: 주요 세포 유형(Neuron, Oligo, Astro, Micro)이 명확히 분리됨
- [ ] **클러스터 수**: resolution 1.0에서 10-20개가 정상 (뇌 조직)
- [ ] **브릿지 세포**: 소수 존재 가능 (transitional state)
- **왜 정상인가**: 논문 Fig. 1c에서 50개 세포 유형이 UMAP에서 명확히 분리. "we identified 50 cell types" (p.814). Xenium의 높은 capture efficiency로 세포 유형 구분이 가능.
- **읽는법**: 각 점=1세포, 색상=클러스터. 잘 분리된 클러스터=뚜렷한 세포 유형. UMAP은 지역 구조 보존에 초점이므로, 클러스터 간 거리는 해석에 주의.
- **비정상 신호**: 모든 세포가 하나의 클러스터 → 전처리 오류, 무수한 작은 클러스터 → 해상도 과도

#### 시각화 5: `spatial_res{r}.png` (공간 클러스터 배치)
- [ ] **조직 구조 보존**: 뇌의 cortical layers, white matter, meninges 등이 보여야 함
- [ ] **Neuron**: gray matter에 집중
- [ ] **Oligodendrocyte**: white matter에 집중
- [ ] **Astrocyte**: gray/white matter 경계에 많음
- **왜 정상인가**: 논문 Fig. 1d에서 세포 유형이 해부학적 조직 구조와 일치. "When assigning these cells to anatomical tissue domains" (p.814). 뇌의 laminar 구조는 well-established anatomy.
- **읽는법**: 각 점=1세포의 물리적 위치 (x, y 좌표). 색상=클러스터. 같은 색의 세포가 공간적으로 모여 있으면 클러스터링이 조직 구조를 반영하는 것.
- **비정상 신호**: 세포 유형이 공간적으로 무작위 분포 → 클러스터링이 세포 유형을 못 포착

#### 시각화 6: `grad_enrichment.png` (Neighborhood Enrichment)
- [ ] **자기 자신 enrichment**: 대각선이 양수 (같은 유형끼리 모임)
- [ ] **Neuron-Oligo**: 음수 (gray vs white matter 분리)
- [ ] **Astrocyte-대부분**: 약한 양수 (넓게 분포)
- **왜 정상인가**: 뇌 조직에서 세포 유형 간 공간적 segregation은 well-known biological pattern. 대각선 양수=같은 유형끼리 이웃하는 경향.
- **읽는법**: 히트맵의 행/열=세포 유형. 빨강=enrichment(이웃할 확률 높음), 파랑=depletion(이웃할 확률 낮음). 대각선은 자기 유형 간 enrichment.
- **비정상 신호**: 모든 값이 0 근처 → 공간적 관계 없음

### 11.3 논문 원문 인용 (시각화 관련)

| 시각화 | 논문 원문 인용 | 페이지 |
|---|---|---|
| Dispersion (Fig. 1f) | "we identified some mRNAs enriched in the nucleus, and others in the cytoplasm" | 814 |
| Dispersion turnover | "transcripts located more than 10.71 µm, on average, from the cell centroid exhibited a higher gene expression correlation with domain-specific background signatures" | 817 |
| ECDF 플랫폼 비교 (Fig. 2f) | "MERFISH and ISS methods displayed a higher concentration of reads near the cellular centroid, whereas ISH-based commercial platforms (CosMx, Molecular Cartography and MERSCOPE) had reads positioned farther away" | 817 |
| UMAP cell types (Fig. 1c) | "we identified 50 cell types that could be mapped onto the tissue to create a cell-type map" | 814 |
| Spatial map (Fig. 1d) | "When assigning these cells to anatomical tissue domains (Methods), we observed a consistent distribution of domain-specific cell types" | 814 |
| Reads/cell 기본값 | "Using Xenium's default segmentation, an average of 186.6 reads per cell was observed throughout the datasets, with 76.8% of reads being assigned to cells" | 813 |
| Subcellular resolution | "Xenium's signal density facilitates the in situ identification of subcellular structures" | 814 |
| P2R clusters | "Using one of these approaches, SSAM de novo mode, we identified 44 cell-type-specific clusters" | 814 |
