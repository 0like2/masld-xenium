# Step 4 Renewal: Resegmentation QC (Cellpose 재분할 품질 검증)

> **작성일**: 2026-02-23
> **목적**: Step 4를 cross-platform 비교에서 **Cellpose 재분할 품질 검증(QC)**으로 리뉴얼
> **기반**: 논문 vs 노트북 vs 파이프라인 구조 분석 결과

---

## 1. 왜 리뉴얼이 필요한가?

### 1.1 현재 Step 4의 문제점

현재 Step 4는 **논문 Section D (Fig. 2)의 cross-platform 비교 메트릭**을 가져와서 **Reseg vs Original 비교**에 적용하고 있음. 이는 **맥락 불일치(context mismatch)**:

```
논문의 원래 맥락 (Section D, Fig. 2):
┌──────────────────────────────────────────────────────┐
│  Xenium vs CosMx vs Vizgen vs MERFISH vs HybRISS    │
│  → 6개 SRT 플랫폼 간 효율/특이성/확산 비교             │
│  → scRNAseq를 공통 reference로 사용                    │
│  → notebooks 3_3 ~ 3_7                               │
└──────────────────────────────────────────────────────┘

파이프라인이 적용한 맥락:
┌──────────────────────────────────────────────────────┐
│  Cellpose Resegmented vs Original (Nuclei-only)      │
│  → 동일 기술의 segmentation 방법 간 비교               │
│  → scRNAseq reference 없을 수도 있음                   │
│  → 논문 Section F (Fig. 3)에 해당                      │
└──────────────────────────────────────────────────────┘
```

### 1.2 구체적 문제

| 현재 분석 | 문제 | 판정 |
|-----------|------|------|
| **Efficiency** (ST/SC ratio) | scRNAseq reference 필요, cross-platform용 메트릭 | 제거 (기본 stats만 유지) |
| **Co-expression NMP** | scRNAseq 필요, O(n²) 느림, cross-platform용 | 제거 |
| **Positivity** (clustering) | cross-platform 비교용 (Fig. 2e), Step 6에서 더 잘 수행 | 제거 |
| **Diffusion** (거리 분석) | Step 6에 없는 유일한 분석, expansion 품질 직접 검증 | **핵심 유지** |
| **Gene correlation** | scRNAseq 불필요, misassignment 간접 감지 | 유지 (선택적) |
| **Basic stats** (counts/genes per cell) | 간단하고 유용, 의존성 없음 | 유지 |

### 1.3 Step 4 vs Step 6 중복 제거

```
Step 6이 이미 하는 것 (cell type annotation 후):
├── Counts/cell violin (Fig. 3g)           → Step 4 Efficiency histogram과 중복
├── Cell type frequency bar (Fig. 3h)       → Step 4 Positivity와 중복
├── UMAP by segmentation (Fig. 3f)          → Step 4 Positivity UMAP과 중복
├── ARI (Fig. 3d)                           → Step 4에 없음
├── NMP scatter (Fig. 3e)                   → Step 4 NMP와 중복 (더 정확)
└── Clustering quality (silhouette 등)       → Step 4에 없음

Step 6에 없는 것 (Step 4만의 고유 가치):
├── ★ Diffusion analysis (transcript-to-centroid distance)
├── Gene-gene correlation heatmap
└── Quick stats summary (세포 수, assigned reads 비율)
```

---

## 2. 리뉴얼 구조

### 2.1 새로운 Step 4 역할 정의

```
기존: "Techniques Comparison" (기술 간 비교)
리뉴얼: "Resegmentation QC" (재분할 품질 검증)
```

**핵심 질문**: "Cellpose 재분할이 원본보다 나은가? expansion이 적절한가?"

### 2.2 리뉴얼 후 Pipeline Flow

```
Step 3: Cellpose Resegmentation
├── step3_resegmented.h5ad (cell × gene matrix)
├── step3_transcripts_resegmented.csv (transcript 좌표 + cell assignment)
└── cell centroids (x_centroid, y_centroid)
         ↓
┌─────────────────────────────────────────────────────────┐
│        [Step 4 (Renewal): Resegmentation QC]            │
│                                                         │
│  Input 요구사항:                                         │
│  - Step 3 출력 (adata + transcripts) ← 필수             │
│  - Step 0/1 출력 (원본 adata + transcripts) ← 비교용     │
│  - scRNAseq reference ← 불필요 (제거)                    │
└─────────────────────────────────────────────────────────┘
         ↓
    ┌───── 4-1. Quick Stats Comparison ────────────────────┐
    │  세포 수, counts/cell, genes/cell, assigned reads    │
    │  → 재분할 효과를 숫자 하나로 요약                       │
    │  → 의존성: adata만 필요                               │
    └──────────────────────────────────────────────────────┘
         ↓
    ┌───── 4-2. Diffusion Analysis (핵심) ─────────────────┐
    │  Transcript-to-centroid distance 분포                 │
    │  → expansion 품질 직접 검증                            │
    │  → 의존성: transcripts + centroids                    │
    └──────────────────────────────────────────────────────┘
         ↓
    ┌───── 4-3. Gene Correlation (선택) ───────────────────┐
    │  Gene-gene correlation heatmap 비교                   │
    │  → misassignment에 의한 비특이적 co-expression 감지    │
    │  → 의존성: adata만 필요 (scRNAseq 불필요)              │
    └──────────────────────────────────────────────────────┘
         ↓
├── figures/4_resegmentation_qc/ 출력 디렉토리
│   ├── qc_stats_comparison.csv
│   ├── qc_stats_comparison.png
│   ├── diffusion_complementary_cdf_comparison.png
│   ├── diffusion_per_gene_ecdf.png
│   ├── diffusion_gene_method_heatmap.png
│   ├── diffusion_per_gene_summary.csv
│   └── gene_correlation_{label}.png (선택)
         ↓
Step 5: Optimal Expansion (diffusion 결과가 expansion 파라미터 선택에 연결)
         ↓
Step 6: Full Segmentation Benchmark (cell type annotation 후 NMP, ARI 등)
```

---

## 3. Sub-step 상세 설계

### 3.1 Quick Stats Comparison (유지/간소화)

**목적**: 재분할 전후 기본 통계 비교 — 한눈에 파악

**데이터 요구사항**: `adata.X` (or `layers['raw']`) + `adata.obs`

**메트릭**:

| 메트릭 | 계산 | 의미 |
|--------|------|------|
| Cell count | `adata.shape[0]` | 재분할로 세포 수가 어떻게 변했는가 |
| Median counts/cell | `np.median(np.sum(X, axis=1))` | 세포당 transcript 포착량 변화 |
| Median genes/cell | `np.median(np.sum(X > 0, axis=1))` | 세포당 유전자 diversity 변화 |
| 5th percentile counts | `np.percentile(counts, 5)` | 저품질 세포의 최소 포착량 |
| Proportion assigned reads | `sum(X) / total_spots` | 전체 transcript 중 세포에 할당된 비율 |

**출력**:
- `qc_stats_comparison.csv` — 정량적 요약 테이블
- `qc_stats_comparison.png` — side-by-side boxplot (counts/cell + genes/cell)

**해석 가이드**:
```
                 Reseg    Original   해석
Cell count:      12,500   8,200      → Cellpose가 더 많은 세포 탐지 (+52%)
Median counts:   145      68         → 확장으로 2.1배 reads 포착 증가
Median genes:    52       38         → 유전자 diversity 1.4배 증가
5th pctile:      15       8          → 저품질 세포도 개선
Assigned reads:  0.72     0.45       → 72% reads가 세포에 할당 (vs 45%)

✓ 모든 지표에서 재분할이 우수 → expansion이 효과적
⚠ counts/cell이 3배 이상 증가하면 → 과도한 확장 의심 (diffusion 확인 필요)
```

### 3.2 Diffusion Analysis (핵심 — 유지/강화)

**목적**: Transcript가 세포 centroid에서 얼마나 멀리 분산되어 있는지 비교

**이것이 핵심인 이유**:
1. **Step 6에 없는 유일한 분석** — Step 4만의 고유 가치
2. **Expansion 품질 직접 검증** — 물리적 거리로 측정
3. **Fig. 3a와 직결** — 논문의 expansion distance 개념
4. **Step 5 optimal expansion 선택에 연결** — 최적 확장 거리 결정의 근거

**논문 근거**:
> "transcripts located more than 10.71 µm, on average, from the cell centroid exhibited
> a higher gene expression correlation with domain-specific background signatures" (p.817)
> → 10.71 µm = misassignment 임계값

**데이터 요구사항**: `transcripts.csv` (좌표 + cell_id) + `adata.obs` (centroids)

**처리 과정**:
```python
# 1. 각 transcript의 할당된 세포 centroid까지 유클리드 거리 계산
distance = sqrt((tx_x - centroid_x)² + (tx_y - centroid_y)²)

# 2. Pixel → µm 변환
distance_um = distance_px / conversion_factor  # Xenium: 4.70588 px/µm

# 3. Reseg vs Original 비교
```

**시각화 (4개)**:

#### 3.2.1 Complementary CDF (`diffusion_complementary_cdf_comparison.png`)

**논문 위치**: Fig. 2f 스타일 (맥락만 reseg vs original로 변경)

```
    P(distance > x)
    1.0 ┤ ╲╲
        │  ╲ ╲
    0.8 ┤   ╲  ╲
        │    ╲   ╲     Reseg (주황)
    0.6 ┤     ╲    ╲── Original (파랑)
        │      ╲     ╲
    0.4 ┤       ╲      ╲
        │        ╲       ╲
    0.2 ┤         ╲        ╲
        │          ╲         ╲
    0.0 ┤           ╲          ╲
        └──┤──┤──┤──┤──┤──┤──┤──→ Distance (µm)
           0  5  10 15 20 25 30

    하단 annotation:
    Reseg: 85% ≤5µm, 95% ≤10µm, med=3.2µm
    Original: 92% ≤5µm, 98% ≤10µm, med=2.1µm
```

**해석법**:
- **곡선이 아래** = 더 concentrated (transcripts가 centroid 가까이)
- **Original이 아래**: 정상 — 핵 내부 reads만이므로 centroid에 가까움
- **Reseg 곡선이 약간 위**: 정상 — expansion으로 cytoplasmic reads 추가
- **Reseg 곡선이 매우 위 (5µm에서 < 60%)**: 과도한 확장 경고
- **10.71 µm 기준선**: 이 이상의 reads 비율이 높으면 misassignment 위험

#### 3.2.2 Per-Gene ECDF (`diffusion_per_gene_ecdf.png`)

```
    Gene1 (n=15000)     Gene2 (n=12000)     Gene3 (n=8000)
    ┌──────────┐         ┌──────────┐         ┌──────────┐
    │╲╲        │         │╲  ╲      │         │╲╲        │
    │  ╲ ╲     │         │ ╲   ╲    │         │  ╲╲      │
    │    ╲  ╲  │         │  ╲    ╲  │         │    ╲╲    │
    │     ╲   ╲│         │   ╲     ╲│         │      ╲╲  │
    └──────────┘         └──────────┘         └──────────┘
    상위 9개 유전자의 개별 ECDF
```

**해석법**:
- **핵 유전자 (예: MALAT1)**: 두 방법 모두 짧은 거리 → 핵에 집중
- **세포질 유전자**: Reseg에서 더 긴 꼬리 → 확장으로 포착
- **유전자 간 차이**: 유전자별 subcellular localization 패턴 반영

#### 3.2.3 Gene × Method Heatmap (`diffusion_gene_method_heatmap.png`)

```
              Reseg   Original
    Gene1     3.2     2.1      ← 차이 작음 = 주로 핵 발현
    Gene2     8.5     2.3      ← 차이 큼 = cytoplasmic transcript
    Gene3     4.1     3.8      ← 비슷 = 핵/세포질 균일 분포
    ...
    (상위 30개 유전자, 색상: 평균 거리 µm)
```

**해석법**:
- **Reseg-Original 차이가 큰 유전자**: cytoplasmic mRNA (expansion 효과 큼)
- **차이가 거의 없는 유전자**: 핵에 집중 발현 (expansion 영향 적음)
- **Reseg에서 >10 µm**: 해당 유전자의 transcript가 매우 멀리 할당됨 → 확인 필요

#### 3.2.4 Per-Gene Summary CSV (`diffusion_per_gene_summary.csv`)

| feature_name | method | mean_distance_um | median_distance_um | std_distance_um | n_transcripts |
|-------------|--------|------------------|--------------------|-----------------|---------------|
| GFAP | Resegmented | 4.2 | 3.1 | 3.5 | 15234 |
| GFAP | Original | 2.1 | 1.8 | 1.9 | 8912 |
| MBP | Resegmented | 6.8 | 5.2 | 4.1 | 12456 |
| ... | | | | | |

### 3.3 Gene Correlation Heatmap (선택적 유지)

**목적**: scRNAseq reference 없이도 misassignment를 간접적으로 감지

**원리**: 과도한 expansion → 이웃 세포의 reads 혼입 → 원래 상관 없는 유전자 간 양의 상관 증가

**데이터 요구사항**: `adata.X` (cell × gene matrix) — scRNAseq **불필요**

```
              Original                          Resegmented
              Gene1 Gene2 Gene3 Gene4            Gene1 Gene2 Gene3 Gene4
    Gene1     1.0   0.1  -0.3   0.8            Gene1  1.0   0.3  -0.1   0.7
    Gene2     0.1   1.0   0.7  -0.2     →      Gene2  0.3   1.0   0.6   0.1
    Gene3    -0.3   0.7   1.0   0.0            Gene3 -0.1   0.6   1.0   0.2
    Gene4     0.8  -0.2   0.0   1.0            Gene4  0.7   0.1   0.2   1.0

    ⚠ 전반적인 양의 상관 증가 → leakage 가능성
```

**해석법**:
- **Reseg에서 음의 상관이 0에 가까워짐**: 서로 다른 세포 유형의 marker가 혼재 → leakage
- **강한 상관 블록 유지**: 세포 유형 특이적 co-expression 패턴 보존 → 양호
- **전체 평균 상관 증가**: ambient RNA 또는 과도한 확장 신호

**출력**:
- `gene_correlation_Resegmented.png` — top 50 gene correlation heatmap
- `gene_correlation_Original.png` — 비교용

---

## 4. 제거 항목 상세 근거

### 4.1 Efficiency - ST/SC Expression Ratio → 제거

| 근거 | 설명 |
|------|------|
| **scRNAseq 의존** | reference 없으면 실행 불가, 있어도 gene name 매칭 복잡 |
| **맥락 불일치** | 원래 cross-platform efficiency 비교용 (Fig. 2c) |
| **Step 6 중복** | counts/cell 비교는 Step 6 violin에서 더 정확하게 수행 |
| **대체** | Quick stats의 counts/genes per cell로 충분 |

### 4.2 Co-expression NMP → 제거

| 근거 | 설명 |
|------|------|
| **scRNAseq 의존** | `_coexpression_calculation(exp_sc)` 필수 |
| **성능 문제** | O(n_genes²) for loop → 매우 느림 |
| **맥락 불일치** | cross-platform co-expression 비교용 (Fig. 2d) |
| **해석 어려움** | co-expression 패턴 차이가 reseg 맥락에서 뭘 의미하는지 불명확 |
| **대체** | Gene correlation heatmap이 scRNAseq 없이 유사한 정보 제공 |
| **참고** | Cell-type NMP는 Step 6에서 cell type annotation 후 정확하게 계산 |

**코드 분석 노트**:
`_negative_marker_purity_coexpression()` (line 51-117)은 `key='celltype'` 파라미터를 받지만 **함수 내부에서 key를 사용하지 않음**. 기술적으로는 cell type 없이 작동하지만, scRNAseq reference (`adata_sc`)는 반드시 필요.

### 4.3 Positivity Analysis → 제거

| 근거 | 설명 |
|------|------|
| **맥락 불일치** | 원래 multi-platform positivity 비교용 (Fig. 2e, notebook 3_5) |
| **Step 6 완전 중복** | UMAP, clustering, cell type bar → Step 6에서 전부 수행 |
| **불필요한 전처리** | Leiden clustering + UMAP → Step 4에서 할 필요 없음 |
| **대체** | Step 6의 joint clustering이 더 정확한 비교 제공 |

### 4.4 Region-specific Efficiency → 제거

| 근거 | 설명 |
|------|------|
| **의존성** | `region_annotation` 필요 (Step 2 domain output) |
| **맥락 불일치** | cross-platform region 비교용 (Ext Fig. 4c) |
| **대체** | Quick stats로 전체적 비교 충분 |

---

## 5. 제거 vs 유지 요약

### 5.1 Before → After

```
현재 Step 4 (34+ 출력 파일)              리뉴얼 Step 4 (~10 출력 파일)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━           ━━━━━━━━━━━━━━━━━━━━━━━━━━━━
├─ Efficiency (13 파일) ──── ✂ 제거     ├─ Quick Stats (2 파일) ── 유지/간소화
│  ├─ ST/SC ratio                      │  ├─ qc_stats_comparison.csv
│  ├─ Region breakdown                 │  └─ qc_stats_comparison.png
│  ├─ Expression ratio                 │
│  └─ Boxplots × 5                     ├─ Diffusion (6 파일) ── ★ 핵심 유지
│                                      │  ├─ complementary_cdf.png
├─ Specificity (6 파일) ── ✂ 제거      │  ├─ per_gene_ecdf.png
│  ├─ Co-expression NMP                │  ├─ gene_method_heatmap.png
│  └─ Correlation heatmap ── ♻ 이동    │  ├─ per_gene_summary.csv
│                                      │  ├─ reseg_diffusion_stats.csv
├─ Positivity (7 파일) ── ✂ 제거       │  └─ original_diffusion_stats.csv
│  ├─ Clustering + UMAP                │
│  ├─ Brain markers                    └─ Gene Correlation (2 파일) ── 선택적
│  └─ Violin plots                        ├─ gene_correlation_Reseg.png
│                                         └─ gene_correlation_Orig.png
├─ Diffusion (8 파일) ── ★ 유지
│  ├─ CDF, ECDF, heatmap
│  └─ Summary CSV
│
└─ 총: 34+ 파일                         총: ~10 파일
   scRNAseq 필요: O                      scRNAseq 필요: X
   실행 시간: 긴                          실행 시간: 짧음
```

### 5.2 Config 변경

```yaml
# 현재
comparison:
  run_efficiency: true       # → 제거
  run_specificity: true      # → 제거
  run_positivity: true       # → 제거
  run_diffusion: true        # → 유지
  sc_reference_path: null    # → 제거 (불필요)
  technology: "xenium"       # → 유지 (pixel→µm 변환)
  positivity:                # → 전체 제거
    n_neighbors: 8
    leiden_resolution: 2.2
    ...

# 리뉴얼
resegmentation_qc:
  run_quick_stats: true          # 기본 통계 비교
  run_diffusion: true            # 핵심: transcript-to-centroid 거리
  run_gene_correlation: false    # 선택: gene-gene correlation heatmap
  technology: "xenium"           # pixel→µm 변환 팩터
  correlation_top_n_genes: 50    # gene correlation heatmap 대상 유전자 수
```

---

## 6. Step 간 연결

### 6.1 Step 4의 새로운 위치

```
Step 3 (Cellpose Resegmentation)
  │ "세포를 다시 나눴다"
  ↓
Step 4 (Resegmentation QC)     ← 리뉴얼
  │ "잘 나눴는지 빠르게 확인"
  │ ├── Quick stats: 세포 수, counts/cell 숫자 확인
  │ ├── Diffusion: expansion이 적절한지 거리로 확인
  │ └── Correlation: misassignment가 없는지 간접 확인
  ↓
Step 5 (Optimal Expansion)
  │ "expansion 거리를 최적화"
  │ ← Step 4 diffusion 결과가 여기 판단의 근거
  ↓
Step 6 (Segmentation Benchmark)
  │ "cell type annotation 후 본격 비교"
  │ ├── NMP (cell-type 기반) — 정확한 specificity 측정
  │ ├── ARI — segmentation 간 일치도
  │ ├── Clustering quality — silhouette, CH, DB
  │ └── Violin, bar, UMAP, spatial maps
  ↓
  "어떤 segmentation이 가장 좋은지 결론"
```

### 6.2 Step 4 → Step 5 연결점

Step 4 Diffusion 결과에서:
- **Median distance ≤ 5 µm**: expansion이 적당 → Step 5에서 현재 distance 유지
- **Median distance > 10 µm**: expansion이 과도 → Step 5에서 distance 줄여야 함
- **Per-gene 차이 큼**: 특정 유전자만 멀리 할당 → gene-specific 분석 필요

### 6.3 Step 4 → Step 6 연결점

Step 4에서 확인한 것:
- Quick stats로 "재분할이 원본보다 좋다" 1차 확인
- Diffusion으로 "expansion이 적절하다" 2차 확인

Step 6에서 최종 확인:
- NMP로 "transcript leakage가 없다" 정밀 확인
- ARI, clustering으로 "세포 분류가 일관적이다" 확인

---

## 7. 논문 Figure 매칭 (리뉴얼 후)

| 논문 Figure | 설명 | 리뉴얼 Step 4 | 비고 |
|------------|------|-------------|------|
| **Fig. 2b** | Transcripts/cell boxplot (cross-platform) | ✂ 제거 | Quick stats가 간접 대체 |
| **Fig. 2c** | ST/SC efficiency ratio | ✂ 제거 | scRNAseq 의존 |
| **Fig. 2d** | Gene specificity (NCP) | ✂ 제거 | Step 6 NMP가 대체 |
| **Fig. 2e** | Positivity by platform | ✂ 제거 | Step 6 clustering이 대체 |
| **Fig. 2f** | Distance to centroid (cross-platform) | ★ **유지** (맥락 변경) | Reseg vs Original로 적용 |
| **Fig. 3a** | Expansion distance 개념 | ★ **연결** | Diffusion이 이 개념 직접 검증 |
| **Fig. 3g** | Counts/cell violin (segmentation 간) | Quick stats로 간소 | Step 6에서 정확한 violin |

---

## 8. 출력 파일 목록 (리뉴얼 후)

### 8.1 Quick Stats (2개)

| # | 파일명 | 유형 | 설명 |
|---|-------|------|------|
| 1 | `qc_stats_comparison.csv` | CSV | 정량적 요약 (cell count, median counts/genes, assigned reads) |
| 2 | `qc_stats_comparison.png` | PNG | Side-by-side boxplot (counts/cell + genes/cell) |

### 8.2 Diffusion Analysis (6개)

| # | 파일명 | 유형 | 설명 |
|---|-------|------|------|
| 3 | `diffusion_complementary_cdf_comparison.png` | PNG | Complementary CDF (5µm/10µm 참조선) |
| 4 | `diffusion_per_gene_ecdf.png` | PNG | 상위 9 유전자별 ECDF |
| 5 | `diffusion_gene_method_heatmap.png` | PNG | Gene × Method mean distance |
| 6 | `diffusion_per_gene_summary.csv` | CSV | Per-gene distance 통계 |
| 7 | `reseg_diffusion_stats.csv` | CSV | Reseg 거리 통계 |
| 8 | `original_diffusion_stats.csv` | CSV | Original 거리 통계 |

### 8.3 Gene Correlation (선택, 2개)

| # | 파일명 | 유형 | 설명 |
|---|-------|------|------|
| 9 | `gene_correlation_Resegmented.png` | PNG | Top 50 gene-gene correlation heatmap |
| 10 | `gene_correlation_Original.png` | PNG | 비교용 |

**총 출력**: 10 파일 (기존 34+ → 70% 감소)

---

## 9. 구현 우선순위

| 순서 | 작업 | 난이도 | 설명 |
|------|------|--------|------|
| 1 | Efficiency/Positivity/NMP 코드 제거 | 낮음 | config flag 비활성화 또는 함수 삭제 |
| 2 | Quick Stats 함수 구현 | 낮음 | 기존 `analyze_efficiency()`에서 histogram 부분만 추출 |
| 3 | Diffusion 함수 유지/정리 | 낮음 | 기존 `analyze_diffusion()` 거의 그대로 유지 |
| 4 | Gene Correlation 분리 | 낮음 | 기존 `analyze_specificity()`에서 correlation 부분만 추출 |
| 5 | Config 구조 변경 | 낮음 | `comparison` → `resegmentation_qc` |
| 6 | `run_step4()` 리팩토링 | 중간 | 새 구조에 맞게 entry point 재작성 |

---

## 10. 핵심 요약

```
┌─────────────────────────────────────────────────────────────────┐
│                    Step 4 리뉴얼 핵심 3줄 요약                    │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  1. Cross-platform 메트릭 (NMP, efficiency ratio, positivity)   │
│     → 전부 제거. Step 6에서 cell type 기반으로 더 정확하게 수행    │
│                                                                 │
│  2. Diffusion (transcript-to-centroid distance)                 │
│     → 핵심 유지. Step 6에 없는 유일한 분석. expansion 품질 직접 검증│
│                                                                 │
│  3. Quick Stats + Gene Correlation                              │
│     → 보조. scRNAseq 불필요. 빠른 1차 확인용                      │
│                                                                 │
│  결과: 34+ 파일 → 10 파일, scRNAseq 의존성 제거, 실행 시간 단축    │
└─────────────────────────────────────────────────────────────────┘
```
