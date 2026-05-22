# 파이프라인 전체 단계별 검토 보고서

**날짜:** 2026-02-18
**브랜치:** step3-refactor
**논문:** Salas et al., "Optimizing Xenium In Situ data utility", Nature Methods, 2025년 4월

---

## 개요

전체 8단계(Step 0~7)를 다음 기준으로 검토:
1. **분석 정확성** - 데이터 무결성, 통계 결과
2. **시각화 품질** - 가독성, 완성도, 미적 요소
3. **논문 Figure 매칭** - 파이프라인 출력 ↔ 논문 Figure 대응 관계

### 출력 현황 요약

| 단계 | 이름 | 상태 | 출력물 |
|------|------|------|--------|
| **Step 0** | Formatting | 정상 | h5ad 1개 |
| **Step 1** | Dataset Exploration | 경미한 이슈 | h5ad 1개 + 플롯 13개 + 통계 |
| **Step 2** | Segmentation Free (P2R) | **이슈 있음** | h5ad 5개 + 플롯 7개 |
| **Step 3** | Resegmentation (Cellpose) | **이슈 있음** | h5ad 1개 + 플롯 2개 |
| **Step 4** | Techniques Comparison | 정상 (Step 3 이슈 반영) | 플롯/통계 16개 |
| **Step 5** | Optimal Expansion | **이슈 있음** | crossover 플롯 100+개 + 요약 |
| **Step 6** | Segmentation Benchmark | **이슈 있음** | h5ad 1개 + 플롯 11개 |
| **Step 7** | Simulation | **무출력 (실패)** | done marker만 존재 |

---

## Step 0: Formatting

### 데이터 요약
- **44,943 cells x 354 genes** (QC 필터링 후)
- 필터링 기준: `min_counts=10`, `min_genes=3`
- Counts/cell: median 216, mean 234.6 (범위: 10~1,016)
- Genes/cell: median 28, mean 28.4 (범위: 3~76)

### 코드 검토
- Control probe 식별 및 제거 (NegControlProbe, NegControlCodeword, antisense, BLANK)
- 10X 분석 결과 import (UMAP, PCA, graph clusters, kmeans 2~10)
- Transcript parquet sidecar 저장 (h5ad 크기 최적화)
- Spatial coordinates → `obsm['spatial']`
- `gene_panel.json` → Ensembl ID 매핑

### 논문 Figure 매칭
- 직접 대응 Figure 없음 (데이터 전처리 단계)

### 발견된 이슈
없음. Step 0은 정상 동작.

---

## Step 1: Dataset Exploration

### 통계 요약 (stats.txt)
- 44,943 cells, 354 genes
- Reads QV>20 비율: 85.6% (논문: 72~91% 범위 내)
- Panel 내 reads 비율: 99.8%
- Cell assigned reads: 75.7% (논문 Fig 1b: 평균 76.8%)
- Control reads: median 0%, mean 0.02% (매우 낮음 - 정상)

### 시각화 검토

| 파이프라인 출력 | 상태 | 논문 Figure 매칭 |
|---|---|---|
| UMAP (Leiden 1.0) | 클러스터 잘 분리되나 **색상이 거의 보이지 않음** - 점 크기/투명도 문제 | Fig 1c (cell types UMAP) |
| Spatial Map (Leiden 1.0) | 조직 구조 보이고 17개 클러스터 표시 | Fig 1d (spatial map) |
| Dispersion distribution | Histogram + KDE 정상, median 6.22 um | Fig 2f 관련 (distance to centroid) |
| Dispersion ECDF | 부드러운 CDF 곡선, 정상 | Fig 2f (cumulative distance) |
| ECDF per gene (Top 10) | 유전자별 분포 차이 보임 - **범례 잘림** | Fig 2f 관련 |
| Violin per gene (Top 20) | 유전자별 distance 분포 적절 | Extended Data Fig 관련 |
| KS heatmap | Pairwise KS test, 유전자 간 유의미한 차이 | 논문 직접 대응 없음 (추가 분석) |
| Markers dotplot | 클러스터별 마커 발현 패턴 | Extended Data Fig 1h,i 관련 |
| Markers heatmap | Dendrogram + 발현 히트맵 정상 | Extended Data Fig 1h,i 관련 |
| HVG plot | Mean-dispersion 관계 정상 | Fig 4 preprocessing 관련 |
| Stats heatmap | 단일 샘플 summary - **정보량 제한적** | Fig 1b (dataset table) |
| Centrality | 클러스터별 centrality scores | 논문 직접 대응 없음 |
| Neighborhood enrichment | 클러스터 간 공간적 enrichment | Extended Data Fig 관련 |

### 발견된 이슈
1. **UMAP 시각화** - 점이 너무 희미해서 클러스터 색상이 거의 안 보임. `alpha`, `size` 조정 필요
2. **ECDF per gene** - 범례(gene name)가 잘려서 어떤 유전자인지 식별 어려움
3. **Stats heatmap** - 단일 샘플이라 히트맵이 1행만 있어 시각적 가치가 낮음

---

## Step 2: Segmentation Free Analysis (Points2Regions)

### 시각화 검토

| 파이프라인 출력 | 상태 | 논문 Figure 매칭 |
|---|---|---|
| Distance to Centroid dist (log) | 정상. Bell-shaped, peak ~5-10 um | Fig 2f 관련 |
| Distance to Boundary (signed) | 정상. Nuclear edge 기준 분포 | Fig 1j 관련 |
| Centroid Boxplot | **좋음.** Nuclear(빨강) vs Cyto(파랑) 유전자 분리 명확 | Fig 1f (subcellular distribution) |
| Centroid ECDF | **좋음.** Nuclear/cytoplasmic 유전자 ECDF 분리 | Fig 2f (cumulative distance) |
| Centroid Stripplot | 정상. Extreme genes 분포 | 보조 시각화 |
| Welch's t-test Heatmap | **좋음.** Nuclear vs Cyto 블록 구조 명확 | 추가 분석 |
| P2R Mean Distance Scatter | **좋음.** Nuclei→Cyto 그래디언트 정상 | Fig 1j (distance to nuclei edge) |
| P2R Celltype Heatmap | **문제: 모든 클러스터가 "Unknown"** | Fig 1i (P2R cluster → celltype) |
| P2R Boundary Boxplot | Nuclear/cyto 분류 자체는 정상 | Fig 1j 관련 |
| P2R Top Genes | **문제: 너무 넓어서 판독 불가** (13755x595px) | Fig 1k (DEGs per cluster) |
| P2R Gene Heatmap | **문제: 너무 넓어서 판독 불가** (11409x889px) | 보조 시각화 |

### 발견된 이슈
1. **P2R 클러스터 cell type이 모두 "Unknown"** - scRNA-seq 참조 데이터와의 annotation transfer가 실패한 것으로 보임. 논문 Fig 1i에서는 Oligo, Astrocyte 등 구체적 cell type이 매핑되어야 함
2. **Top genes / Gene heatmap이 너무 넓음** - 50개 클러스터 x nuclear/cyto = ~100개 패널이 가로로 나열되어 해상도가 부족
3. **Boundary distance가 1400 um까지** - 핵에서 먼 unassigned transcript들이 포함된 것으로 보임

---

## Step 3: Resegmentation (Cellpose)

### 데이터 검토
- **39,100 cells x 539 genes** (Step 0: 44,943 x 354)
- **539 genes가 354보다 많음** - 원본 Xenium output의 모든 feature 사용 (control probe 포함 가능성)
- **Counts/cell 크게 감소**: median 35 (Step 0: 216) - 약 84% 감소
- **min counts=1, min genes=1** - QC 필터링이 적용되지 않음
- `obsm['spatial']` 없음 - 좌표가 obs에만 저장됨

### 시각화 검토

| 파이프라인 출력 | 상태 | 논문 Figure 매칭 |
|---|---|---|
| Mask QC (expanded + nuclei) | 정상. 39,583 labels, 조직 형태 유지 | Fig 3c (segmentation comparison) |
| Transcript overlay | **이슈: 가장자리에 큰 artifact 원형** (좌상단 흰색 영역) | Fig 3c 관련 |

### 발견된 이슈
1. **Gene 수 불일치** (539 vs 354) - Control probe가 포함된 것으로 보임
2. **Counts/cell 급감** - median 216→35, 약 84% 감소. Expansion distance나 transcript 할당 로직 확인 필요
3. **QC 필터링 미적용** - min_counts=1인 셀이 존재
4. **Transcript overlay 가장자리 artifact** - 조직 밖 영역에 큰 원형 mask들

---

## Step 4: Techniques Comparison

### 핵심 지표
- **NMP Score**: Resegmented 0.979 > Original 0.960 (nuclei-only이므로 더 높은 순도)
- **Read 할당률**: Resegmented ~10% vs Original ~75% (**심각**)

### 시각화 검토

| 파이프라인 출력 | 상태 | 논문 Figure 매칭 |
|---|---|---|
| Expression ratio histogram | 정상. ratio~1 중심 분포 | Fig 2c 관련 (detection efficiency) |
| ST vs scRNA scatter | 대부분 identity 아래 - ST 검출 낮음 | Fig 2c, 2g (ST/SC ratio) |
| Transcripts/cell comparison | **핵심 비교**: Reseg peak~20 vs Original peak~200 | Fig 2b (transcripts per cell) |
| Genes/cell comparison | Reseg median~14 vs Original median~28 | Fig 2b (genes per cell) |
| Reseg vs Original boxplot | **명확한 감소** 확인 | 보조 시각화 |
| Gene correlation (Original) | 블록 구조 정상 | Fig 2d 관련 |
| Gene correlation (Resegmented) | 약한 correlation - counts 부족 영향 | Fig 2d 관련 |
| NMP per gene boxplot | **좋음.** Reseg NMP=0.979 > Original 0.960 | Fig 3e (NMP scatter) |
| Diffusion assigned reads | **심각: Reseg 10% vs Original 75%** | 논문 직접 대응 없음 |
| Diffusion complementary CDF | Reseg 꼬리 ~80um vs Original ~20um | Fig 2f (cumulative distance) |
| Diffusion gene-method heatmap | **좋음.** 유전자별 방법 간 거리 비교 | 보조 시각화 |
| Diffusion per-gene ECDF | 9개 유전자 개별 비교, 명확한 차이 | 보조 시각화 |
| Positivity distribution | Reseg에서 0% positive 유전자 급증 | 보조 시각화 |
| Positivity UMAP top genes | 발현 패턴 정상 | 보조 시각화 |
| Positivity violin clusters | 유전자별 클러스터 분포 정상 | 보조 시각화 |

### 핵심 발견
1. **Resegmentation의 transcript 할당률이 ~10%로 극히 낮음** (Original 75%). Cellpose nuclei-only 모델이 expansion 없이 핵만 감지해서 대부분의 cytoplasmic transcript를 놓치는 것으로 보임
2. NMP는 오히려 Reseg가 높음 (0.979 vs 0.960) - nuclei-only이므로 더 순수한 signal
3. 이 결과 자체는 논문의 핵심 메시지와 일치: **"nuclei segmentation captures purer signal but loses reads"**

---

## Step 5: Optimal Expansion

### 핵심 결과
- **Optimal expansion = -4.4 um (음수!)** - expansion이 필요 없다는 의미
- Mean nuclei size: 4.83 um
- Mean turnover: 0.43

### 시각화 검토

| 파이프라인 출력 | 상태 | 논문 Figure 매칭 |
|---|---|---|
| Optimal expansion txt | **문제: -4.4 um (음수)** | Fig 3b (predicted optimal expansion) |
| Expansion distances histogram | 분포 정상 (대부분 0 근처) | 보조 시각화 |
| Turnover barplot | Cell type별 turnover - 생물학적으로 합리적 | Fig 3a 관련 (PCC vs distance) |
| Turnover filtered (>5 domains) | 7개 cell type 통과 | Fig 3b 관련 |
| Expansion map | 조직 위 색상 분포 정상 | 보조 시각화 |
| Reads vs centroids | 조직 커버리지 정상 | 보조 시각화 |
| Crossover plots (~100개) | **좋음.** corr_nuc vs corr_back crossover 패턴 | Fig 3a (PCC vs distance) |

### 발견된 이슈
1. **Optimal expansion = -4.4 um (음수)** - 논문에서는 cell type별로 다른 양수 optimal distance를 보여줌 (Fig 3b). 전체 평균이 음수인 것은 비정상적
2. **Turnover distance 값이 매우 작음** (대부분 0~1 um) - crossover point가 너무 일찍 발생하거나 correlation 계산에 문제 가능성
3. **Crossover plot** - Oligodendrocyte 예시에서 corr_nuc는 감소, corr_back은 증가하는 올바른 패턴이지만, correlation 값이 0.97~1.0 범위로 매우 좁아서 crossover 감지가 어려움

---

## Step 6: Segmentation Benchmark

### 시각화 검토

| 파이프라인 출력 | 상태 | 논문 Figure 매칭 |
|---|---|---|
| UMAP benchmark (3 panels) | **문제: 너무 작고 흐림**, segmentation method 패널은 OK | Fig 3f (UMAP Baysor vs Nuclei) |
| Spatial segmentation map | **Baysor 크기가 비정상적으로 작음** (crop region만) | Fig 1d, 3c 관련 |
| Spatial celltype map | 같은 이슈 - Baysor 패널만 작음 | 보조 시각화 |
| ARI heatmap | **expansion vs baysor ARI=0.015 - 매우 낮음** | Fig 3d (ARI heatmap) |
| Counts violin by method | **좋음.** Baysor median~6500 >> 나머지. Nuclei~200, Cellpose~50 | Fig 3g (violin cell counts) |
| Cell type frequency barplot | **좋은 시각화.** Baysor에서 Astrocyte 과대, VLMC 과대 | Fig 3h (cell type barplot) |
| DEG dotplot by segmentation | **좋음.** 방법별 마커 유전자 차이 | 보조 시각화 |
| Marker genes by celltype | 세부 panel 많지만 readable | 보조 시각화 |
| Marker genes by segmentation | 방법별 top DEG 비교 | 보조 시각화 |
| HVG selection plot | 정상 | Fig 4 관련 |
| PCA scree plot | 정상. PC9에서 90% 도달 | 보조 시각화 |
| QC histogram | **문제: Counts 분포가 0 근처 밀집** (대부분 methods의 counts가 낮음) | 보조 시각화 |

### 발견된 이슈
1. **Baysor가 crop region에서만 실행됨** - config에서 `crop.enabled: true, size_um: 2000`으로 설정. 정상 동작이지만 다른 방법들과 공정한 비교가 안 됨
2. **ARI가 expansion vs baysor = 0.015** - 사실상 random. crop region 차이가 원인일 수 있음
3. **UMAP이 너무 작고 cluster가 97개** - 해상도가 너무 높거나 combined dataset이 너무 heterogeneous
4. **Cellpose counts가 nuclei보다 낮음** (violin plot) - Step 3의 문제가 여기서도 반영됨

---

## Step 7: Simulation

### 상태: **무출력 (실패)**
- `step7_done.txt` marker만 존재하고 **실제 출력(plot, h5ad)이 전혀 없음**
- `cellxgene_census` 미설치로 reference 다운로드 실패 → 전체 Step 7이 skip된 것으로 추정
- `pipeline_main.py`에서 try/except로 예외를 잡고 done marker를 쓰고 있음 → **실패했지만 done으로 마킹된 상태**

### 논문 Figure 매칭
- Fig 4a-e (simulation workflow, preprocessing heatmap, ARI comparison)에 해당하는 출력이 모두 누락

### 발견된 이슈
1. **Silent failure** - done marker가 실제 성공 여부와 무관하게 기록됨
2. **Missing dependency** - `cellxgene_census` 패키지 필요
3. **pipeline_main.py 395~397행** - try/except 내에서 done marker를 쓰므로 실패가 은폐됨

---

## 우선순위별 이슈 요약 (수정 현황: 2026-02-18)

### 심각 (논문 재현성 차단)
1. **Step 7 무출력** - ✅ 수정완료: done marker가 성공 후에만 기록, import guard 추가
2. **Step 3 낮은 transcript 할당률** (~10%) - ✅ 수정완료: `expand_labels()` 적용, dual assignment 구현
3. **Step 2 P2R 클러스터 모두 "Unknown"** - ✅ 수정완료: cell_id 타입 캐스팅(str) + spatial nearest-neighbor fallback 추가

### 높음 (시각화 품질)
4. **Step 2 P2R top genes/heatmap 판독 불가** - ✅ 수정완료: grid layout (max 10열/행) + heatmap 최대 너비 30inch 제한
5. **Step 6 Baysor crop vs 전체 조직 비교** - ✅ 수정완료: Dual benchmark mode (full tissue + crop 비교)
6. **Step 5 음수 optimal expansion** - ✅ 수정완료: true crossover detection (smoothing + above→below 전환 요구) + nanmedian aggregation
7. **Step 1 UMAP 너무 희미** - ✅ 수정완료: size 자동 스케일 (120000/n_cells) + alpha=0.7 + frameon=False

### 보통 (외형적/경미)
8. **Step 1 ECDF 범례 잘림** - ✅ 수정완료: figsize 확대 (12,6) + fontsize=8 + framealpha=0.9
9. **Step 1 Stats heatmap** - ✅ 수정완료: 단일 샘플 시 bar chart로 자동 전환
10. **Step 3 gene 수 539 vs 354** - ✅ 수정완료: NegControl/BLANK/antisense 필터링 추가
11. **Step 3 QC 필터링 미적용** - ✅ 수정완료: sc.pp.filter_cells(min_counts, min_genes) + edge cell 감지
12. **Step 6 UMAP 97개 클러스터** - ✅ 수정완료: 40개 초과 시 resolution 자동 0.5x 감소 + 점 크기 자동 스케일

### 추가 수정사항
13. **Step 2 Boundary distance 1400µm까지 확장** - ✅ 수정완료: 1st/99th percentile 기준 x축 clipping
14. **Step 3 Background label(0) 포함 버그** - ✅ 수정완료: cell_id_reseg > 0 필터 추가
15. **Step 6 UMAP 점 크기** - ✅ 수정완료: auto-scale (80000/n_cells) + alpha=0.7

---

## 논문 Figure ↔ 파이프라인 매칭 종합표

| 논문 Figure | 설명 | 파이프라인 단계 | 상태 |
|---|---|---|---|
| Fig 1a | 워크플로우 개요 | 해당 없음 | 해당 없음 |
| Fig 1b | 데이터셋 통계 테이블 | Step 1 (stats_heatmap) | 부분적 |
| Fig 1c | UMAP cell types | Step 1 (umap_res1.0) | 희미함 |
| Fig 1d | Spatial map | Step 1 (spatial_res1.0) | 정상 |
| Fig 1e | 3D coherence (ovrlpy) | Step 2 (ovrlpy) | 미실행 (config) |
| Fig 1f | Subcellular distribution | Step 2 (centroid_boxplot) | **좋음** |
| Fig 1g,h | Spatial transcript maps | 미구현 | 누락 |
| Fig 1i | P2R clusters | Step 2 (p2r_celltype_heatmap) | **모두 Unknown** |
| Fig 1j | Distance to nuclei edge | Step 2 (p2r_mean_distance_scatter) | **좋음** |
| Fig 1k | P2R cluster별 DEG | Step 2 (p2r_topgenes) | 판독 불가 |
| Fig 2a | 비교 워크플로우 | 해당 없음 | 해당 없음 |
| Fig 2b | Transcripts/genes per cell | Step 4 (comparison 플롯) | 정상 |
| Fig 2c | Detection efficiency | Step 4 (expression_ratio) | 정상 |
| Fig 2d | Gene specificity NCP | Step 4 (nmp_per_gene_boxplot) | **좋음** |
| Fig 2e | 유전자별 violin | Step 1 (dispersion_violin) | 정상 |
| Fig 2f | Cumulative distance to centroid | Step 1/2 (ECDF 플롯) | **좋음** |
| Fig 2g | Xenium vs Visium | 해당 없음 (단일 데이터셋) | 해당 없음 |
| Fig 3a | Distance + PCC | Step 5 (crossover_plots) | **좋음** |
| Fig 3b | Cell type별 optimal expansion | Step 5 (turnover_barplot) | **음수 값** |
| Fig 3c | Segmentation 비교 이미지 | Step 3 (mask_qc) + Step 6 | 부분적 |
| Fig 3d | ARI heatmap | Step 6 (rand_index_heatmap) | 낮은 ARI |
| Fig 3e | NMP scatter | Step 4 (nmp_per_gene_boxplot) | **좋음** |
| Fig 3f | UMAP Baysor vs Nuclei | Step 6 (umap_benchmark) | 너무 작음 |
| Fig 3g | Cell counts violin | Step 6 (counts_violin_by_method) | **좋음** |
| Fig 3h | Cell type barplot | Step 6 (celltype_frequency_barplot) | **좋음** |
| Fig 4a-e | Simulation/Preprocessing | Step 7 | **실패** |
| Fig 5a-h | Imputation/Domain ID | 미구현 | 해당 없음 |
