# MASLD Xenium Pipeline: 종합 검토 보고서

> **최종 업데이트:** 2026-02-12
> **브랜치:** `step3-refactor`
> **논문 참조:** Marco Salas et al., "Optimizing Xenium In Situ data utility by quality assessment and best-practice analysis workflows", *Nature Methods* 22, 813–823 (April 2025). [DOI: 10.1038/s41592-025-02617-2](https://doi.org/10.1038/s41592-025-02617-2)

---

## 1. 프로젝트 개요

### 1.1 논문이 제시하는 목표

이 논문은 10x Genomics Xenium In Situ 플랫폼에서 생성되는 공간 전사체 데이터에 대한 **독립적 품질 평가와 최적 분석 워크플로우**를 제시한다. 25개 Xenium 데이터셋(14개 실험, 총 12억 reads, 600만 cells)을 분석하여 다음을 수행한다:

1. **데이터 품질 특성화** — 검출 효율, 특이성, 재현성 평가 (Fig. 1)
2. **세포 분할 비교를 위한 공통 재분할** — Cellpose 기반 nuclei segmentation 후 expansion (Fig. 2a, 3c)
3. **기술 간 비교** — Xenium vs CosMx, MERSCOPE, MERFISH, HS-ISS, MC, Visium (Fig. 2)
4. **최적 확장 반경 산출** — Nuclear vs background 유전자 발현 상관관계 교차점(turnover) 분석 (Fig. 3a,b)
5. **분할 전략 벤치마크** — Baysor+Cellpose가 최적 (NMP + assigned reads 기준) (Fig. 3d,e)
6. **전처리 최적화** — 시뮬레이션 기반 전처리 그리드 서치로 최적 워크플로우 도출 (Fig. 4)
7. **SVF 식별 및 도메인 탐색** — SpatialDE, Squidpy, Hotspot 등 비교 (Fig. 4f-j)
8. **유전자 임퓨테이션** — SpaGE, Tangram 등 7개 방법 벤치마크 (Fig. 5)

### 1.2 파이프라인의 역할

이 파이프라인(`pipeline/`)은 위 논문의 핵심 분석을 **재현 가능한 end-to-end Python 워크플로우**로 구현한다. 원본 Jupyter 노트북(`notebooks/`)의 분석 로직을 자동화하고, 새로운 Xenium 데이터셋에 바로 적용할 수 있도록 한다.

### 1.3 Step-Notebook-논문 매핑

| Pipeline Step | Notebook Source | 논문 Figure/Section | 핵심 목적 |
|---|---|---|---|
| **Step 0** — Formatting | `0_formatting/` | Methods | Xenium 원시 출력 → AnnData 변환 |
| **Step 1** — Exploration | `1_datasets_exploration/` | Fig. 1b-f | 데이터 품질 통계, 분산(dispersion) 분석, 클러스터링 |
| **Step 2** — Segmentation-Free | `2_segmentation_free_analysis/` | Fig. 1e-k | Points2Regions 기반 세포 미할당 분자 패턴 분석 |
| **Step 3** — Resegmentation | `3_techniques_comparison/3_1` | Fig. 2a, 3c | Cellpose nuclei segmentation + label expansion |
| **Step 4** — Comparison | `3_techniques_comparison/3_3-3_7` | Fig. 2c-g | 효율성/특이성/양성률/확산 비교 분석 |
| **Step 5** — Optimal Expansion | `4_optimal_expansion/` | Fig. 3a,b | Correlation 기반 turnover → 최적 확장 반경 |
| **Step 6** — Benchmark | `5_segmentation_benchmark/` | Fig. 3d-h | 분할 전략별 클러스터링 품질 비교 |
| **Step 7** — Simulation | `6_simulating_preprocessing/` | Fig. 4a-e | 시뮬레이션 기반 전처리 파라미터 최적화 |

> **참고:** 노트북 `7_domain_exploration/`, `8_SVF_identification/`, `9_spatial_domain_annotation/`, `10_gene_imputation/`은 현재 파이프라인에 포함되지 않음 (향후 확장 후보).

---

## 2. Step별 상세 검토 (최신 코드 기준)

### Step 0: Formatting (`xenium_step0_formatting.py`)

#### 구현 상태: **양호 (85%)**

| 기능 | 구현 | 논문/노트북 대응 |
|------|:---:|---|
| Xenium raw output → AnnData (sparse CSR) | **O** | Methods |
| 압축 파일 자동 해제 (.gz, .tar.gz) | **O** | — |
| 음성 대조군 분리 (NegControlProbe, BLANK) | **O** | Fig. 1b QC metrics |
| Transcripts → sidecar .parquet 저장 | **O** | 노트북보다 개선 |
| Analysis 결과 로딩 (UMAP/PCA/tSNE/KMeans) | **O** | — |
| gene_panel.json → Ensembl ID 매핑 | **O** | — |
| OME-TIFF → plain TIFF 변환 | **O** | — |
| QC 필터링 (min_counts, min_genes) | **O** | Fig. 1b (0.21% cells <10 reads 제외) |
| `overlaps_nucleus=1` 핵 필터링 | **config** | 노트북 `0_3` (config `filter_nuclei_only` 옵션) |
| R용 CSV 내보내기 (Seurat 호환) | X | 노트북 `0_1`, `0_2` — Python 파이프라인에서 불필요 |

#### 미구현 (낮은 우선순위)
- R(Seurat) 포맷 내보내기 — Python 전용 파이프라인이므로 선택사항

---

### Step 1: Dataset Exploration (`xenium_step1_dataset_exploration.py`)

#### 구현 상태: **양호 (80%)**

**논문 대응:** 이 단계는 논문 Fig. 1b의 데이터셋 특성 요약 및 Fig. 1f의 subcellular patterning 분석의 기초를 담당한다.

| 기능 | 구현 | 논문 대응 |
|------|:---:|---|
| QC 통계 (n_cells, n_genes, median, QV>20 등) | **O** | Fig. 1b 요약 테이블 |
| Transcript dispersion (distance-to-centroid) | **O** | Fig. 1f (subcellular distribution) |
| KS 검정 (유전자 간 거리 분포 비교) | **O** | 노트북 `1_3` |
| ECDF 플롯 (전체 + 유전자별 top 10) | **O** | Fig. 2f (cumulative proportion) |
| Violin plot (top 20 유전자) | **O** | 노트북 `1_3` |
| HVG 선별 + Leiden 다중 해상도 클러스터링 | **O** | Fig. 1c cell type identification |
| Marker gene 추출 (Wilcoxon, dotplot + heatmap) | **O** | Fig. 1c,d |
| Neighborhood enrichment (Squidpy) | **O** | Fig. 1d 조직 구조 |
| Neighborhood centrality scores | **O** | 노트북 `1_7` |
| 통계 요약 히트맵 | **O** | Fig. 1b normalized metrics |
| KS p-value 히트맵 | **O** | — |

#### 미구현 (낮은 우선순위)
- 다수 데이터셋 간 배치 비교 (논문은 25개 데이터셋 사용, 파이프라인은 단일 데이터셋 처리)
- Cell architecture structure scores (노트북 `1_7` ADDITIONAL)

---

### Step 2: Segmentation-Free Analysis (`xenium_step2_segmentation_free_analysis.py`)

#### 구현 상태: **보통 (65%)**

**논문 대응:** 논문의 "Xenium retains key 3D and subcellular cell information" 섹션 (Fig. 1e-k). Points2Regions(Ref. 10)로 세포 분할 없이 분자 클러스터를 식별하고, subcellular 패턴을 분석.

| 기능 | 구현 | 논문 대응 |
|------|:---:|---|
| Points2Regions 비지도 공간 클러스터링 | **O** | Fig. 1i,j (Points2Regions clusters) |
| Distance-to-centroid (유전자별) | **O** | Fig. 1f (distance-to-centroid boxplot) |
| Distance-to-boundary (EDT 기반) | **O** | Fig. 1j (nuclear edge distance) |
| 유전자별 거리 통계 CSV | **O** | — |
| Extreme genes stripplot/boxplot (inward/outward) | **O** | Fig. 1f (most nuclear/least nuclear) |
| ECDF by gene | **O** | — |
| Welch's t-test p-value 히트맵 | **O** | 노트북 `2_1` |
| SSAM 분석 | X | Extended Data Fig. 2a (SSAM de novo mode) |
| Overlaps (ovrlpy) 분석 | **placeholder** | Fig. 1e (z-dimension overlaps) |

#### 잔여 이슈
- **SSAM**: config에 파라미터 존재하지만 실행 코드 미구현. 논문에서 SSAM은 44개 세포타입 클러스터를 식별하는 데 사용 (Extended Data Fig. 2a)
- **Overlaps**: `ovrlpy` placeholder만 존재 — config 기본값이 `run_overlaps: false`로 설정됨
- Points2Regions 결과를 nuclear/cytoplasmic/extracellular로 분류하는 로직 없음 (Extended Data Fig. 3b)

---

### Step 3: Resegmentation (`xenium_step3_resegmentation.py`)

#### 구현 상태: **양호 (90%)**

**논문 대응:** 논문의 핵심 전략 — "To facilitate a fair comparison, cells were resegmented using a common segmentation algorithm (Cellpose), and reads were reassigned to individual cells" (Fig. 2a). 논문은 nuclei-only segmentation + 제한적 expansion (<10-30% reads assigned)을 권장.

| 기능 | 구현 | 논문 대응 |
|------|:---:|---|
| Cellpose nuclei segmentation (GPU/MPS/CPU) | **O** | Fig. 2a, Methods |
| Cellpose v4.0.8 API 호환 (pretrained_model) | **O** | — |
| `expand_labels()` configurable expansion | **O** | Fig. 3b (expansion distances 1-15µm) |
| Dual transcript assignment (`in_cell`/`closest_cell`) | **O** | 노트북 `3_1` 핵심 로직 |
| Nuclei mask + expanded mask 별도 저장 (.tif) | **O** | — |
| Regionprops: centroid, area, perimeter | **O** | — |
| Cell area px→µm² 변환 | **O** | Config `um_per_pixel_inv` |
| Distance-to-centroid per transcript | **O** | Fig. 2f, 3a |
| Mask + transcript overlay 시각화 | **O** | — |
| Domain assignment (JSON polygons) | **O** | Fig. 1d (anatomical regions) |
| Mask QC 시각화 (nuclei vs expanded) | **O** | — |

#### 최근 수정사항 (2026-02-12)
- Cellpose v4 API 호환: `model_type` → `pretrained_model`, `channels` 제거, 반환값 3개
- `pipeline_main.py`: Step 3 marker 파일이 실제 output h5ad 존재 시에만 생성되도록 수정

---

### Step 4: Techniques Comparison (`xenium_step4_techniques_comparison.py`)

#### 구현 상태: **양호 (75%)**

**논문 대응:** 논문 Fig. 2c-g의 핵심 분석. "Xenium detection efficiency matches ISH" 섹션에서 6개 SRT 기술의 검출 효율/특이성/확산을 비교. 여기서는 원본(nuclei) vs 재분할(Cellpose) 비교에 초점.

| 기능 | 구현 | 논문 대응 |
|------|:---:|---|
| **Efficiency** — transcripts/genes per cell histogram | **O** | Fig. 2b |
| **Efficiency** — ST/scRNAseq expression ratio (CPM) | **O** | Fig. 2c (SRT/SC ratio) |
| **Efficiency** — region-based breakdown + CSV | **O** | Fig. 2g (per-region scatter) |
| **Specificity** — NMP (Negative Marker Purity) | **O** | Fig. 3e (NMP scatter) |
| **Specificity** — per-gene purity CSV + boxplot | **O** | Fig. 2d (NCP score) |
| **Specificity** — NMP vs Efficiency scatter | **O** | Fig. 3e |
| **Specificity** — gene-gene correlation heatmap | **O** | — |
| **Positivity** — gene detection rate histogram | **O** | 노트북 `3_5` |
| **Positivity** — Leiden clustering + violin plots | **O** | 노트북 `3_5` (per-cluster expression) |
| **Positivity** — UMAP per top gene + optimal cluster CSV | **O** | — |
| **Diffusion** — px→µm 단위 변환 | **O** | Fig. 2f (distance in µm) |
| **Diffusion** — Complementary CDF (1-CDF) plot | **O** | Fig. 2f |
| **Diffusion** — per-gene ECDF subplots | **O** | Extended Data Fig. 4f |
| **Diffusion** — Gene×Method distance heatmap | **O** | — |
| **Diffusion** — per-gene summary CSV | **O** | — |
| Spatial ROI map (domain-colored cells) | **O** | — |
| `spatial_utils.py` 공유 변환 함수 | **O** | — |

#### 잔여 이슈
- Efficiency의 expression ratio 계산은 scRNAseq reference가 필요 (`sc_reference_path` in config)
- 6개 기술 간 비교는 단일 데이터셋 파이프라인 구조상 미지원 (논문은 25개 데이터셋 배치 비교)
- Xenium vs Visium 비교 (노트북 `3_7`) 미구현

---

### Step 5: Optimal Expansion (`xenium_step5_optimal_expansion.py`)

#### 구현 상태: **양호 (85%)**

**논문 대응:** 논문의 "Nuclear expansion influences cell-type expression profiles" 섹션 (Fig. 3a,b). 핵심 발견: "transcripts located more than 10.71 µm from the cell centroid exhibited a higher correlation with domain-specific background signatures". 최적 확장 반경은 세포타입별로 다름 (Fig. 3b).

| 기능 | 구현 | 논문 대응 |
|------|:---:|---|
| cKDTree 기반 domain assignment (청크 처리) | **O** | Methods |
| Nuclear vs background expression 프로필 구축 | **O** | Fig. 3a (PCC curves) |
| Distance-bin correlation curves | **O** | Fig. 3a (corr vs distance plot) |
| Turnover distance 검출 (corr 교차점) | **O** | Fig. 3a (crossover point) |
| Nuclei size 계산 (ConvexHull, `dist_nuc()`) | **O** | Fig. 3b (nuclei edge) |
| Cell size 계산 | **O** | Fig. 3b (cell edge) |
| `optimal_expansion = turnover - nuclei_size` | **O** | Fig. 3b 핵심 공식 |
| Domain별 반복 분석 (min_reads_per_domain 필터) | **O** | Methods (>5000 reads) |
| Celltype별 turnover summary CSV | **O** | — |
| Turnover barplot (cell_size/nuclei_size overlay) | **O** | Fig. 3b |
| Per-domain crossover line plots | **O** | Fig. 3a (per cell type curves) |
| `domain_key` vs `celltype_key` 분리 | **O** | — |

#### 최근 수정사항 (2026-02-12)
- `domain_key`와 `celltype_key`가 동일 소스 컬럼을 사용하던 HIGH 버그 수정
- Fallback path `step4_reseg` → `step3_reseg` 수정
- `dist_nuc()` missing distance column 경고 추가
- `pipeline_main.py`: Step 5 marker 파일이 실제 output CSV 존재 시에만 생성

---

### Step 6: Segmentation Benchmark (`xenium_step6_segmentation_benchmark.py`)

#### 구현 상태: **양호 (80%)**

**논문 대응:** 논문의 "Baysor and Cellpose outperform standard Xenium segmentation" 섹션 (Fig. 3c-h). Baysor(BA2 P0.8) + Cellpose nuclei segmentation이 최적 전략으로 도출됨. "Although cells defined by Baysor had a higher count per cell, the identified cellular populations were the same across both segmentation strategies" (Fig. 3f-h).

| 기능 | 구현 | 논문 대응 |
|------|:---:|---|
| **Baysor 실행** — real execution + dry_run mode | **O** | Fig. 3c,e (Baysor segmentation) |
| **Baysor** — early binary check (shutil.which) | **O** | — |
| **Baysor** — prior segmentation TIF 지원 | **O** | Fig. 3e (prior segm. confidence) |
| Transcript CSV → AnnData aggregation | **O** | — |
| Multi-method loading (nuclei/cellpose/expansion/baysor) | **O** | Fig. 3d (52 strategies) |
| Annotation transfer (kNN majority voting, k=15) | **O** | Fig. 3f (UMAP cell types) |
| Per-cell annotation confidence scores | **O** | — |
| Preprocessing: normalize/log1p/HVG/PCA/Leiden | **O** | Methods (Louvain→Leiden) |
| Preprocessing params from config | **O** | Fig. 4c (target_sum=100, etc.) |
| QC visualizations (counts/genes hist, HVG, PCA scree) | **O** | — |
| UMAP (3 panels: method/leiden/celltype) | **O** | Fig. 3f |
| Spatial scatter map (per method, per celltype) | **O** | Fig. 3c (ROI comparison) |
| Cell type frequency barplot | **O** | Fig. 3h |
| Counts violin by method | **O** | Fig. 3g |
| Clustering quality: Silhouette/CH/DB scores | **O** | — |
| NMP metrics (cells + reads, if scRNA ref) | **O** | Fig. 3e (NMP axis) |
| Rand Index (pairwise method comparison) | **O** | Fig. 3d (ARI heatmap) |
| Marker gene ranking (Wilcoxon) + DEG dotplot | **O** | — |
| 5th percentile metrics | **O** | — |

#### 최근 수정사항 (2026-02-12)
- Baysor `shutil.which()` 체크를 데이터 준비 전으로 이동 (불필요한 transcript 로딩 방지)
- Silhouette/CH/DB 클러스터링 품질 메트릭 추가
- Annotation confidence scores 추가

#### 잔여 이슈
- 논문은 52개 segmentation strategy를 비교 (Fig. 3d) — 파이프라인은 4개(nuclei/cellpose/expansion/baysor)
- Baysor 바이너리가 서버에 설치되어 있어야 함

---

### Step 7: Simulation & Preprocessing Optimization (`xenium_step7_simulation.py`)

#### 구현 상태: **보통 (65%)**

**논문 대응:** 논문 Fig. 4 전체. "Preparing Xenium data: best practices in preprocessing" 섹션. CellxGene Census에서 scRNAseq를 다운로드하여 Xenium과 유사한 데이터를 시뮬레이션하고, 전처리 조합을 그리드 서치하여 최적 워크플로우 도출. **논문의 최적 워크플로우:** (1) library-size normalization (target_sum=100), (2) log1p, (3) scaling, (4) all PCs + 16 neighbors, (5) Louvain clustering.

| 기능 | 구현 | 논문 대응 |
|------|:---:|---|
| CellxGene Census scRNAseq 다운로드 | **O** | Fig. 4a |
| 데이터셋 필터링 (max_cells, min/max_celltypes) | **O** | Methods (8k-20k, 2-30 types) |
| Marker gene 선별 (t-test/Wilcoxon) | **O** | Fig. 4a (subsetting profiled genes) |
| Noise simulation (random ±1) | **O** | Fig. 4a (introducing unspecific reads) |
| Missegmentation simulation (cell mixing) | **O** | Fig. 4a (simulating mis-segmented cells) |
| Dual simulation (standard + high noise) | **O** | — |
| Preprocessing grid search | **O** | Fig. 4b,c (ranked workflows) |
| 4 metrics: ARI, NMI, FMI, VI | **O** | Fig. 4d |
| Perturbation analysis (one-at-a-time) | **O** | Fig. 4d (parameter sensitivity) |
| ARI/NMI/FMI bar charts | **O** | Fig. 4b |
| ARI heatmap (n_neighbors × n_pcs) | **O** | — |
| Parameter importance bar chart | **O** | — |
| Perturbation grouped bar chart | **O** | — |

#### 잔여 이슈
- 논문은 다수 데이터셋(tissues)에 걸친 cross-dataset ARI 히트맵 생성 (Fig. 4b) — 파이프라인은 단일 tissue
- Preprocessing grid 조합 수가 config에 의존 (기본 ~50개, 논문은 유사)
- Real data에 대한 parameter-tuning 검증 (Fig. 4e) 미구현

---

### pipeline_main.py Orchestrator

#### 구현 상태: **양호 (90%)**

| 기능 | 구현 |
|------|:---:|
| Step 0→7 순차 실행 | **O** |
| 이전 스텝 결과 재사용 (h5ad/marker 존재 시 skip) | **O** |
| Step 간 데이터 경로 자동 전달 | **O** |
| DAPI 이미지 경로 자동 탐색 (3개 후보) | **O** |
| `transcripts.parquet` 폴백 지원 | **O** |
| Step 3 marker 조건부 생성 (output 확인) | **O** |
| Step 5 marker 조건부 생성 (output 확인) | **O** |
| Domain map 자동 연결 (Step 2 → Step 3) | **O** |
| Baysor prior TIF 자동 연결 (Step 3 → Step 6) | **O** |
| Config 기반 server/local 경로 전환 | **O** |

---

## 3. 전체 요약: 수정 이력 및 잔여 이슈

### 3.1 해결된 CRITICAL/HIGH 이슈 (2026-02 리팩토링)

| # | 이전 상태 | 수정 내용 |
|---|---|---|
| ~~C1~~ | Step 2: distance-to-boundary만 구현 | **해결** — distance-to-centroid 추가, config `distance_metric` 옵션 |
| ~~C2~~ | Step 5: turnover/crossover 미구현 | **해결** — correlation-based turnover 알고리즘 완전 구현 |
| ~~C3~~ | Step 6: Baysor 시뮬레이션만 | **해결** — 실제 subprocess 실행 + graceful skip + dry_run 모드 |
| ~~C4~~ | Step 4: Diffusion px→µm 변환 누락 | **해결** — `spatial_utils.py` PIXEL_TO_UM_FACTORS 적용 |
| ~~H1~~ | Step 4: Efficiency = QC metric만 | **해결** — expression ratio (ST/scRNAseq) 계산 추가 |
| ~~H2~~ | Step 6: 전처리 파라미터 하드코딩 | **해결** — config에서 모든 파라미터 로드 |
| ~~H3~~ | Step 1: ECDF, violin, KS 누락 | **해결** — ECDF/violin/KS test/heatmap 모두 구현 |
| ~~H4~~ | Step 3: expand_labels 미구현 | **해결** — configurable expansion_distance 추가 |
| ~~H5~~ | Step 5: domain_key = celltype_key 버그 | **해결** — 별도 소스 컬럼으로 분리 |
| ~~H6~~ | Step 3: Cellpose v4 API 비호환 | **해결** — pretrained_model, 반환값 3개, channels/tile 제거 |
| ~~H7~~ | pipeline_main: 스킵된 step에 marker 생성 | **해결** — output 파일 존재 확인 후 marker 생성 |

### 3.2 잔여 이슈 (우선순위별)

#### MEDIUM — 기능 보완

| # | Step | 이슈 | 논문 참조 |
|---|---|---|---|
| M1 | Step 2 | SSAM 분석 미구현 (placeholder) | Extended Data Fig. 2a |
| M2 | Step 2 | Points2Regions 결과를 nuclear/cytoplasmic/extracellular로 분류 | Extended Data Fig. 3b |
| M3 | Step 7 | Real data parameter-tuning 검증 | Fig. 4e |
| M4 | Step 4 | Xenium vs Visium 비교 (노트북 `3_7`) | Fig. 2g |

#### LOW — 선택적 확장

| # | Step | 이슈 | 논문 참조 |
|---|---|---|---|
| L1 | — | SVF 식별 (SpatialDE, Squidpy, Hotspot 등) | Fig. 4f-j |
| L2 | — | Domain exploration (SpaGCN, BANKSY 등) | Fig. 5g,h |
| L3 | — | Gene imputation (SpaGE, Tangram 등) | Fig. 5a-f |
| L4 | Step 0 | R(Seurat) 포맷 내보내기 | — |
| L5 | Step 2 | Overlaps (ovrlpy) 실제 구현 | Fig. 1e |

---

## 4. 논문의 핵심 결론과 파이프라인 구현 상태

| 논문 핵심 결론 | 파이프라인 | 상태 |
|---|---|---|
| Xenium 검출 효율은 scRNA-seq의 1.2-1.5x (Chromium v2) | Step 4 efficiency ratio | **O** |
| NCP > 0.8로 높은 특이성 (HS-ISS, MC와 유사) | Step 4/6 NMP score | **O** |
| Baysor + Cellpose nuclei = 최적 분할 전략 (BA2 P0.8) | Step 6 benchmark | **O** |
| Cellpose expansion ≠ rigid expansion; 세포타입별 최적 반경이 다름 | Step 5 turnover analysis | **O** |
| 최적 전처리: normalize(100) + log1p + scale + all PCs + 16 neighbors + Louvain | Step 7 grid search | **O** |
| 3D subcellular 정보 보존: nuclear vs cytoplasmic gene 패턴 | Step 2 distance analysis | **O** |
| Points2Regions로 segmentation-free 세포타입 식별 가능 | Step 2 Points2Regions | **O** |
| SVF 알고리즘 간 일관성이 낮음 (Kendall tau 상관관계 낮음) | 미구현 (향후 L1) | X |
| Gene imputation: SpaGE가 가장 높은 PCC 달성 | 미구현 (향후 L3) | X |
| Domain identification: binning 기반 방법이 수동 annotation과 가장 유사 | 미구현 (향후 L2) | X |

---

## 5. 아키텍처 및 기술 참고

### 5.1 Config 구조 (`config.yaml`)
- 서버 경로: `/data/project/lyrak/masld-xenium/data/...`
- 로컬 Mac 개발 경로는 주석으로 보존
- Baysor 바이너리 경로: `baysor.executable_path` (서버에서 `which baysor`로 확인 필요)
- 각 Step별 enable/disable 플래그 제공 (`run_resegmentation`, `run_expansion` 등)

### 5.2 공유 유틸리티
- `pipeline/utils/spatial_utils.py` — `PIXEL_TO_UM_FACTORS` dict, 좌표 변환 함수
- `pipeline/benchmark_utils/metrics.py` — ARI, NMI, FMI, VI 계산
- `xb/calculating.py` — `dispersion()`, `dist_nuc()` 등 원본 라이브러리

### 5.3 의존성
- `pyproject.toml`로 PEP 621 패키징 (139개 의존성)
- 진입점: `xenium-pipeline = "pipeline.pipeline_main:main"`
- 설치: `pip install -e .` → `xenium-pipeline --config pipeline/config.yaml`

---

## 6. 검증 방법

```bash
# 서버에서 실행
git pull origin step3-refactor
pip install -e .
rm -f xenium-output/*/step*_done.txt   # 기존 buggy marker 제거
xenium-pipeline --config pipeline/config.yaml

# Step별 개별 검증
python -c "from pipeline import xenium_step3_resegmentation as s3; print('OK')"

# Import 전체 검증
python -c "
import sys; sys.path.insert(0, 'pipeline')
for m in ['xenium_step0_formatting','xenium_step1_dataset_exploration',
          'xenium_step2_segmentation_free_analysis','xenium_step3_resegmentation',
          'xenium_step4_techniques_comparison','xenium_step5_optimal_expansion',
          'xenium_step6_segmentation_benchmark','xenium_step7_simulation']:
    __import__(m); print(f'  {m}: OK')
"
```

---

*이 보고서는 논문 Marco Salas et al. (Nature Methods 2025)의 분석 워크플로우를 기준으로 파이프라인 구현 상태를 검토한 것입니다.*
