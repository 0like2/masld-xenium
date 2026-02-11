# Pipeline vs Notebook 전체 검토 보고서

## Context
`notebooks/` 디렉토리의 Jupyter 노트북들을 `pipeline/` 디렉토리에서 Python 스크립트로 재구현하고 있음. 각 스텝별로 빠진 부분, 로직 오류, 시각화 누락 등을 철저히 검토함.

**매핑 관계:**
| Pipeline Script | Notebook Source |
|---|---|
| step0_formatting | notebooks/0_formatting/ |
| step1_dataset_exploration | notebooks/1_datasets_exploration/ |
| step2_segmentation_free_analysis | notebooks/2_segmentation_free_analysis/ |
| step3_resegmentation | notebooks/3_techniques_comparison/3_1_resegmentation_notebooks/ |
| step4_techniques_comparison | notebooks/3_techniques_comparison/3_3~3_7 |
| step5_optimal_expansion | notebooks/4_optimal_expansion/ |
| step6_segmentation_benchmark | notebooks/5_segmentation_benchmark/ |
| step7_simulation | notebooks/6_simulating_preprocessing/ |

---

## Step 0: Formatting (`xenium_step0_formatting.py`)

### 구현 상태: **양호 (80%)**

**잘 구현된 부분:**
- Xenium 원시 데이터(matrix.mtx, features.tsv, cells.csv) → AnnData 변환
- 압축 파일(.gz, .tar.gz) 자동 해제
- 음성 대조군(NegControlProbe 등) 처리 (노트북에 없는 추가 기능)
- Transcripts를 sidecar .parquet로 저장 (노트북보다 개선된 방식)
- Analysis 폴더에서 UMAP/PCA/tSNE/KMeans 결과 로딩

**누락된 기능:**
1. **핵(nuclei) 필터링 미구현** — 노트북 `0_3`에서 `overlaps_nucleus=1` 필터링을 하지만 파이프라인에 없음
2. **R용 CSV 내보내기 미구현** — 노트북 `0_1`, `0_2` (Seurat용 포맷). 이 부분은 R 전용이므로 Python 파이프라인에서는 선택사항
3. **시뮬레이션 데이터 포매팅 미구현** — 노트북 `0_2`

**수정 권장:**
- nuclei 필터링 옵션 추가 (config에서 `filter_nuclei_only: true/false`로 제어)

---

## Step 1: Dataset Exploration (`xenium_step1_dataset_exploration.py`)

### 구현 상태: **보통 (60%)**

**잘 구현된 부분:**
- 기본 QC 통계 계산 (n_cells, n_genes, median_genes/counts 등)
- Transcript dispersion 분석 (거리 계산 + 히스토그램)
- HVG 계산 및 Leiden 클러스터링 (다중 해상도)
- Marker gene 추출 (Wilcoxon test, top 5/cluster)
- Neighborhood enrichment 분석 (Squidpy 사용)

**누락된 시각화 (중요!):**
1. **ECDF(누적 분포) 플롯** — 노트북 `1_3`에서 유전자/샘플별 ECDF 커브가 핵심 결과물
2. **Violin plot** — 거리 분포의 바이올린 플롯 (소스별 비교)
3. **유전자×세포타입 거리 히트맵** — 세포 내 위치 분석의 핵심
4. **KS 검정(통계 테스트)** — 샘플 소스 간 Kolmogorov-Smirnov 검정 완전 누락
5. **Marker gene 랭킹 히트맵** — dotplot만 있고, 상세 히트맵 누락
6. **통계 요약 히트맵** — 노트북은 데이터셋 그룹별 히트맵 생성, 파이프라인은 텍스트 파일만 저장

**수정 권장:**
- `step1_1`에 ECDF plot 추가 (`sns.ecdfplot` 또는 수동 구현)
- Violin plot 추가
- KS 검정 결과 CSV 저장 추가
- 출력 포맷을 TXT → CSV로 변경 (분석 재사용성)

---

## Step 2: Segmentation-Free Analysis (`xenium_step2_segmentation_free_analysis.py`)

### 구현 상태: **미흡 (40%)**

**잘 구현된 부분:**
- Points2Regions 분석 구현
- Distance Transform(EDT) 기반 boundary 거리 계산
- 거리 히스토그램 시각화

### **치명적 로직 차이:**
- **파이프라인: distance-to-boundary (핵막까지 거리)**
- **노트북: distance-to-centroid (핵 중심까지 거리)**
- 이 두 메트릭은 완전히 다른 의미를 가짐. 노트북은 `xb.calculating.dispersion()`으로 centroid 거리를 계산

**누락된 기능:**
1. **유전자별 거리 분석** — 노트북 `2_1`의 핵심: inward/outward 유전자 분류 (stripplot, boxplot)
2. **Welch's t-test** — 유전자 간 거리 비교 p-value 매트릭스
3. **ECDF 플롯** — 유전자별 누적 분포 곡선
4. **Overlaps 분석 미구현** — `ovrlpy` import만 하고 실제 호출 없음 (placeholder만 존재, line 249: `pass`)
5. **SSAM 분석 미구현** — config에 파라미터 있지만 실행 코드 없음
6. **유전자별 통계 CSV 내보내기** — 노트북은 유전자별 거리 통계를 CSV로 저장

**수정 권장 (우선순위 순):**
1. distance-to-centroid 옵션 추가 (가장 중요)
2. 유전자별 시각화 (stripplot, ECDF) 추가
3. ovrlpy 호출 구현 또는 config에서 `run_overlaps: false`로 기본값 변경
4. SSAM은 별도 스텝으로 분리하거나 명시적 미구현 표시

---

## Step 3: Resegmentation (`xenium_step3_resegmentation.py`)

### 구현 상태: **양호 (70%)**

**핵심 기능 구현됨** (Cellpose 기반 재분할). 별도의 상세 노트북 비교 필요하나, 기본 로직은 구현되어 있음.

---

## Step 4: Techniques Comparison (`xenium_step4_techniques_comparison.py`)

### 구현 상태: **미흡 (35%)**

**구현된 부분:**
- QC 메트릭 기반 효율성 비교 (히스토그램)
- NMP 점수 계산 (참조 데이터 있을 때)
- Gene-gene correlation 히트맵 (proxy)
- Positivity 비율 계산
- Diffusion 거리 계산 + ECDF 플롯

### **치명적 문제:**

#### 4-1. Efficiency (노트북 3.3)
- **노트북**: ST/scRNAseq 발현 비율(expression ratio) 계산 → 6개 방법 비교
- **파이프라인**: 단순 QC 메트릭(total_counts, n_genes) 비교 → **완전히 다른 메트릭**
- 영역별(Cortex/Hippocampus/Thalamus) 분석 누락
- 프로브 수 vs 효율성 상관 분석 누락

#### 4-2. Specificity (노트북 3.4)
- 효율성 비율 기반 유전자 사전필터링 누락 (ratio < 10 조건)
- 세포타입별 purity 분석 누락 (파이프라인은 단일 NMP 점수만)
- Efficiency vs Specificity 2D scatter plot 누락

#### 4-3. Positivity (노트북 3.5)
- 전처리 파이프라인 완전 누락 (normalize, log1p, clustering, UMAP)
- Leiden 클러스터링 없음 → 유전자별 최적 클러스터 식별 불가
- Violin plot 누락 (유전자 발현 분포)

#### 4-4. Diffusion (노트북 3.6) — **단위 변환 버그**
- **노트북**: 픽셀→마이크로미터 변환 적용 (방법별 변환 계수 사용)
  - Xenium: `/4.70588`, CosMx: `/8.3333`, etc.
- **파이프라인**: 단위 변환 없음 → 거리값이 임의 픽셀 단위
- 유전자별 ECDF 서브플롯 누락
- Gene×Method 거리 히트맵 누락
- Complementary CDF (1-CDF) 대신 표준 CDF 사용

**수정 권장:**
- Diffusion에 단위 변환 로직 추가 (최소한 Xenium 변환 계수)
- 실제 expression ratio 기반 efficiency 계산으로 교체
- 전처리 파이프라인 추가 (positivity 분석용)

---

## Step 5: Optimal Expansion (`xenium_step5_optimal_expansion.py`)

### 구현 상태: **미흡 (45%)**

**구현된 부분:**
- cKDTree 기반 도메인 할당 ✓ (잘 구현됨)
- 청크 처리로 대용량 데이터 지원 ✓
- 거리 임계값 필터링 ✓
- 확장 결과 시각화 (scatter plot, 거리 히스토그램) ✓

### **핵심 누락: Turnover/Crossover 분석 미완성**

노트북의 핵심 알고리즘이 구현되지 않음:

1. **Correlation 기반 turnover 계산 누락**
   - 노트북: 각 거리 bin에서 `corr(expression[dist], nucleus_expression)` vs `corr(expression[dist], background_expression)` 계산
   - 두 correlation이 교차하는 지점 = turnover distance
   - 파이프라인: correlation 계산 없음, `dist_nuc()` (ConvexHull)만 계산

2. **최적 확장 반경 계산 공식 누락**
   - `optimal_expansion = mean_turnover - mean_nuclei_size`
   - 파이프라인은 nuclei_size만 계산하고 turnover를 계산하지 않으므로 최종 결과를 도출할 수 없음

3. **유전자별 Crosstab 분석 누락**
   - 노트북: `pd.crosstab(distance, feature_name)` → 거리별 유전자 발현 패턴
   - 파이프라인: 이 분석 전체 누락

4. **도메인별 반복 분석 누락**
   - 노트북: 각 세포타입 × 각 도메인(>5000 reads)에 대해 반복
   - 파이프라인: 세포타입별만 반복 (도메인 내부 루프 없음)

5. **시각화 누락**
   - 노트북: correlation difference vs distance 라인플롯 (turnover 지점 표시)
   - 노트북: 세포타입별 nuclei_size(빨간) + cell_size(검정) + turnover(바) 오버레이 플롯

**수정 권장:**
- turnover 계산 핵심 알고리즘 구현 (correlation-based crossover)
- 도메인별 루프 추가
- 최종 optimal expansion 계산식 구현

---

## Step 6: Segmentation Benchmark (`xenium_step6_segmentation_benchmark.py`)

### 구현 상태: **미흡 (30%)**

### **치명적 문제:**

1. **Baysor 실행이 시뮬레이션됨** (line 131-140)
   - `subprocess.run(cmd, check=True)` 가 주석 처리됨
   - 빈 CSV 파일을 생성하여 파이프라인 계속 진행
   - 실제 Baysor 세분화 결과 없이 벤치마크 무의미

2. **전처리 파라미터 불일치**
   - 노트북: `target_sum=100, min_counts=40, min_genes=15, n_neighbors=15, umap_min_dist=0.1`
   - 파이프라인: 모든 파라미터 기본값 사용 (결과 품질 저하)

3. **클러스터링 알고리즘 불일치**
   - 노트북: **Louvain** (resolution=1.0)
   - 파이프라인: **Leiden** (resolution 미지정) → 다른 결과

4. **Annotation Transfer 미구현**
   - 노트북: reference adata에서 majority voting으로 세포타입 전이
   - 파이프라인: leiden 클러스터만 생성, 세포타입 annotation 없음

5. **누락된 시각화:**
   - Spatial map (세포타입별 색상)
   - Cell type frequency barplot (방법 간 비교)
   - Violin plot (방법별 counts 분포)

6. **UMAP 저장 로직 불안정** (line 298-301)
   - `figures/umap.png` 경로 의존 (scanpy 기본 저장) → 불안정

**수정 권장:**
- Baysor 실제 실행 활성화 (또는 명시적 skip 옵션)
- 전처리 파라미터를 config에서 로드하도록 수정 + Louvain 옵션 추가
- Annotation transfer 구현
- 시각화 추가

---

## Step 7: Simulation (`xenium_step7_simulation.py`)

### 구현 상태: **미흡 (35%)**

**구현된 부분:**
- CellxGene Census에서 scRNAseq 다운로드
- 기본 시뮬레이션 (noise, missegmentation)
- 4개 메트릭 계산 (ARI, NMI, FMI, VI) ✓

**누락된 기능:**
1. **데이터셋 필터링 로직** — 노트북: 8k-20k 세포, 2-30 세포타입 필터링 → 파이프라인: 전체 다운로드
2. **전처리 그리드 축소** — 노트북: 30+ 조합 → 파이프라인: 9개(3×3)만
3. **파라미터 범위 축소** — normalization 방법, target_sum 변형, scale 변형 등 누락
4. **Perturbation 분석 누락** — 최적 워크플로우에서 파라미터 1개씩 변경 효과 테스트
5. **Cross-dataset 히트맵 누락** — 다수 데이터셋에 걸친 ARI 히트맵
6. **시각화 대폭 축소** — 노트북 8종류 → 파이프라인 2종류(ARI bar, NMI bar)만

---

## pipeline_main.py 오케스트레이션

### 구현 상태: **양호**
- Step 0→7 순차 실행
- 이전 스텝 결과 재사용 (파일 존재 시 skip)
- Step 간 데이터 경로 전달
- DAPI 이미지 경로 자동 탐색

**문제점:**
- Step 6 import가 `run_step6` 내부에서 동적 import (line 309) — 일관성 부족
- `from pipeline.benchmark_utils import metrics` (step6 line 10) — 상대 import 문제 가능

---

## 전체 요약: 심각도별 분류

### CRITICAL (즉시 수정 필요)
1. **Step 2**: distance-to-boundary vs distance-to-centroid 로직 차이
2. **Step 5**: turnover/crossover 핵심 알고리즘 미구현 → 최적 확장 반경 계산 불가
3. **Step 6**: Baysor 실행이 시뮬레이션 (빈 CSV 생성)
4. **Step 4**: Diffusion 단위 변환 누락 (픽셀 vs 마이크로미터)

### HIGH (기능 완성도에 큰 영향)
5. **Step 4**: Efficiency가 QC 메트릭 비교로 대체됨 (expression ratio 아님)
6. **Step 6**: 전처리 파라미터/알고리즘 불일치 (Louvain→Leiden, 파라미터 기본값)
7. **Step 7**: 전처리 그리드 30+ → 9개로 대폭 축소
8. **Step 1**: KS 검정, ECDF, violin plot 등 핵심 통계/시각화 누락
9. **Step 2**: ovrlpy 분석 placeholder만 존재 (pass)

### MEDIUM (시각화/분석 보완 필요)
10. **Step 1**: 통계 요약 히트맵 → 텍스트 파일로 단순화
11. **Step 2**: 유전자별 거리 분석 (stripplot, boxplot) 누락
12. **Step 4**: Positivity 분석에서 전처리/클러스터링 미수행
13. **Step 5**: 도메인별 반복 분석 + correlation 시각화 누락
14. **Step 6**: Annotation transfer 미구현
15. **Step 7**: Perturbation 분석 + cross-dataset 히트맵 누락

### LOW (선택적 개선)
16. **Step 0**: nuclei 필터링 옵션 추가
17. **Step 0**: R용 CSV 내보내기 (Python 파이프라인에서는 불필요할 수 있음)
18. **Step 6**: UMAP 저장 경로 안정화

---

## 수정 계획

### Phase 1: Critical 수정 (4건)
1. **Step 2** — distance-to-centroid 계산 함수 추가 + 기존 boundary 방식과 병행 옵션
2. **Step 5** — correlation-based turnover 알고리즘 구현 + optimal expansion 공식 완성
3. **Step 6** — Baysor 실제 실행 활성화 (주석 해제 + 에러 핸들링)
4. **Step 4** — Diffusion에 Xenium 단위 변환 계수 추가

### Phase 2: High 수정 (5건)
5. **Step 4** — expression ratio 기반 efficiency 계산으로 교체
6. **Step 6** — 전처리 파라미터를 config에서 로드하도록 수정 + Louvain 옵션 추가
7. **Step 7** — 전처리 그리드 확장 (최소 20+ 조합)
8. **Step 1** — ECDF plot, violin plot, KS test 추가
9. **Step 2** — ovrlpy placeholder를 실제 구현으로 교체하거나 기본값 false로 변경

### Phase 3: Medium 시각화 보강 (6건)
10~15번 항목의 시각화 및 분석 기능 추가

### 검증 방법
- 각 스텝을 개별 실행하여 출력물(h5ad, CSV, PNG) 생성 확인
- 시각화 결과물이 노트북 출력과 유사한 형태인지 비교
- `python pipeline/pipeline_main.py` 전체 파이프라인 실행 테스트
