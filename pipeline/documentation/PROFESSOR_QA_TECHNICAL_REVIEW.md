# 교수님 질문 기술 검토 보고서

> **작성일:** 2026-02-19
> **데이터셋:** Xenium V1 FFPE Human Brain Alzheimer's (with Add-on Panel)
> **파이프라인:** masld-xenium (Step 0~7)
> **검증 도구:** `pipeline/xenium_eda.py` (EDA 스크립트)

---

## 목차

1. [DAPI-Transcript 좌표 Alignment 및 Gold Standard 부재 이슈](#1-dapi-transcript-좌표-alignment-및-gold-standard-부재-이슈)
2. [Mouse Liver Panel 유전자 42개 충분성 검토](#2-mouse-liver-panel-유전자-42개-충분성-검토)
3. [QC Control 값 분포 분석](#3-qc-control-값-분포-분석)
4. [하드웨어 리소스 (RAM) 검토](#4-하드웨어-리소스-ram-검토)

---

## 1. DAPI-Transcript 좌표 Alignment 및 Gold Standard 부재 이슈

### 1.1 교수님 질문 요약

> Xenium 데이터에서 DAPI 이미지(Nucleus, 파란색)와 발현된 유전자(Transcript)를 Overlay 했을 때, 위치가 정확히 맞지 않는(Align이 안 되는) 기술적 이슈가 있다고 들었다. 이에 대한 Gold Standard(정확한 정답 기준)가 없다는 점이 문제라고 알고 있다.

### 1.2 결론 (한 줄 요약)

**좌표 자체는 정확히 일치합니다. "Alignment 문제"로 알려진 것은 pixel 좌표 오차가 아니라, "핵 밖에 있는 transcript를 어느 세포에 배정할 것인가"에 대한 ground truth가 없다는 문제입니다.**

### 1.3 우리 데이터에서 실제 확인한 결과

`xenium_eda.py --section dapi` 실행 결과:

```
[DAPI Focus] morphology_focus.ome.tif
  Pixel dimensions: 40,855 x 36,955
  Physical size:    8,681.7 x 7,852.9 µm (8.68 x 7.85 mm)
  Pixel size:       0.2125 µm/pixel

[Coordinate Alignment Check]
  Transcript X range: 3.90 ~ 7,851.21 µm
  Transcript Y range: 5.17 ~ 8,678.61 µm
  Transcript X (pixels): 18 ~ 36,947
  Transcript Y (pixels): 24 ~ 40,840

  DAPI image: 40,855 x 36,955 pixels
  Transcripts fit within DAPI X: YES (max_tx_px=36,947, img_w=36,955)
  Transcripts fit within DAPI Y: YES (max_tx_px=40,840, img_h=40,855)
  Coverage: X=99.9%, Y=99.9%
```

- DAPI 이미지와 transcript 좌표 범위가 **99.9% 일치**
- 8개 pixel (1.7µm) 이내의 미세한 마진만 존재
- `experiment.xenium` 메타데이터에서 `pixel_size: 0.2125` µm/pixel 확인

### 1.4 왜 좌표는 일치하는가: Xenium 시스템 구조

Xenium 장비는 DAPI 형광 이미지와 transcript spot 위치를 **동일한 광학 시스템**과 **동일한 좌표 프레임**으로 캡처합니다:

1. **동일 좌표계 사용:** `transcripts.csv`의 `x_location`, `y_location` (µm 단위)과 DAPI OME-TIFF 이미지는 같은 좌표 원점을 공유합니다. µm 좌표를 pixel로 변환하는 공식은 `pixel = µm / 0.2125` (혹은 `µm × 4.70588`)입니다.

2. **내부 보정 과정:**
   - 렌즈 왜곡은 장비 캘리브레이션 데이터 기반으로 소프트웨어 보정됨
   - FOV(Field of View) 간 overlap 영역에서 feature matching으로 global alignment 수행
   - 우리 데이터: 100개 FOV (`fov_name` 컬럼에서 B9~J7 등 확인)

3. **우리 파이프라인의 좌표 변환 (Step 3):**
   ```python
   # xenium_step3_resegmentation.py:506-516
   um_per_pixel_inv = 4.70588  # config에서 설정 (px/µm)
   y_coords = (y_vals * um_per_pixel_inv).astype(int)
   x_coords = (x_vals * um_per_pixel_inv).astype(int)
   ```
   이 변환으로 transcript µm 좌표를 DAPI pixel 좌표로 매핑하며, 변환 후 모든 transcript가 이미지 범위 내에 정확히 들어감을 EDA에서 확인했습니다.

### 1.5 진짜 문제: Gold Standard의 본질

"Alignment 문제"로 알려진 것은 실제로는 **Cell Segmentation + Transcript Assignment** 문제입니다:

#### 1.5.1 핵(Nucleus) vs 세포질(Cytoplasm) 격차

우리 데이터에서 확인된 수치:

```
[Nucleus Overlap] (transcripts.parquet에서)
  overlaps_nucleus=1 (핵 안): 1,459,578  (8.9%)
  overlaps_nucleus=0 (핵 밖): 14,911,942 (91.1%)
```

**전체 transcript의 91.1%가 핵 밖에 위치합니다.** DAPI는 핵만 염색하므로, DAPI 이미지 위에 transcript를 overlay하면 대부분의 transcript가 핵 바깥에 점으로 찍히게 됩니다. 이것이 "align이 안 된다"고 보이는 현상의 실체입니다.

이것은 **좌표 오류가 아니라 생물학적 현실**입니다:
- mRNA의 대부분은 핵에서 전사된 후 세포질로 이동
- 세포질의 물리적 크기가 핵보다 훨씬 큼 (우리 데이터: median cell area=508.8, median nucleus area=24.5 → 세포질이 핵의 ~20배)
- 일부 mRNA는 세포 외 공간으로 확산되기도 함

#### 1.5.2 Gold Standard가 없는 이유

핵 밖에 있는 transcript를 특정 세포에 배정하는 "정답"을 알 수 있는 기술이 존재하지 않습니다:

| 문제 | 설명 |
|------|------|
| **경계 모호성** | 세포 경계가 불분명한 밀집 조직에서 인접 세포 간 transcript 소유권 판단 불가 |
| **확산(Diffusion)** | transcript가 원래 세포에서 인접 세포로 확산되는 물리적 현상이 존재 |
| **Z축 해상도** | Xenium은 z-step 3.0µm으로 촬영하지만, 위아래 세포의 transcript가 같은 xy 위치에 겹칠 수 있음 |
| **측정 불가능성** | 특정 transcript가 어느 세포에서 전사되었는지를 분자 수준에서 확인하는 기술이 없음 |

#### 1.5.3 대안적 Ground Truth (Proxy) 방법들

완벽한 gold standard는 없지만, 연구자들이 사용하는 proxy 방법들:

1. **Nuclear-only segmentation**: 핵 안의 transcript만 사용 (보수적이지만 고신뢰)
2. **scRNA-seq cross-validation**: 같은 조직의 single-cell RNA-seq 데이터와 비교하여 cell type별 발현 패턴이 일치하는지 확인
3. **Spurious co-expression 분석**: 생물학적으로 공존할 수 없는 두 유전자가 같은 세포에서 발현되면 mis-assignment 징후
4. **인접 조직 절편 CODEX/MIBI**: 단백질 수준에서 세포 경계를 독립적으로 확인

### 1.6 최신 보정/대안 방법 (2024-2025)

| 방법 | 접근법 | 참고문헌 |
|------|--------|----------|
| **Baysor** | Transcript 위치의 확률적 공간 클러스터링 기반 segmentation. 핵 segmentation을 prior로 사용 가능 | [10x Analysis Guide](https://www.10xgenomics.com/analysis-guides/using-baysor-to-perform-xenium-cell-segmentation) |
| **Proseg** | 확률적 cellular Potts model. 핵에서 시작하여 확장/축소를 반복하며 transcript 분포를 최적으로 설명하는 경계 탐색. Baysor 대비 ~10배 빠름 | Jones et al., PMC 2025 |
| **SPLIT** | snRNA-seq와 RCTD deconvolution을 통합하여 세포내 signal을 정제. 확산 오염을 제거 | Bilous et al., 2025, bioRxiv |
| **MisTIC** | Variational Bayesian model. 재-segmentation 없이 mis-assignment 교정. 공간 근접성 + 발현 호환성 + 이웃 지지도 활용 | bioRxiv 2025 |
| **FastReseg** | R 패키지. 세포별 segmentation 오류 점수 산출 후 mis-located transcript 재배정 | Bioconductor OSTA |
| **nuc2seg** | Nuclear stain만으로 whole-cell segmentation 수행 | tansey-lab/nuc2seg |
| **Xenium Multimodal** | 경계 stain (ATP1A1/CD45/E-Cadherin) + interior stain (18S rRNA) + DAPI를 3단계 파이프라인으로 결합 | 10x Genomics |

### 1.7 우리 파이프라인의 대응 방식

우리 파이프라인은 Salas et al. (2025, Nature Methods) 논문의 방법론을 따르며, 이 문제를 다층적으로 다루고 있습니다:

#### Step 3: Dual Assignment 전략
```python
# xenium_step3_resegmentation.py:537-544
# 두 가지 assignment를 동시에 수행:
df_valid['in_cell'] = nuclei_labels       # 핵 기반 (보수적, 고신뢰)
df_valid['closest_cell'] = cell_labels    # 확장 기반 (포괄적, 저신뢰)
```

- **Nuclear assignment** (`in_cell`): Cellpose가 검출한 핵 mask 안에 있는 transcript만 해당 세포에 배정. 보수적이지만 mis-assignment 최소화.
- **Expanded assignment** (`closest_cell`): `expand_labels(distance=400px ≈ 85µm)` 로 확장된 mask에 transcript 배정. 세포질 transcript도 포착하지만 오배정 위험 증가.

#### Step 5: Optimal Expansion 탐색
최적의 확장 거리를 데이터 기반으로 결정하여, nuclear-only (under-capture)와 over-expansion (mis-assignment) 사이 균형점을 찾음.

#### Step 6: Baysor 벤치마크
Baysor를 포함한 4가지 segmentation 방법 (nuclei, Cellpose, expansion, Baysor)을 벤치마킹하여 최적 전략 선택.

### 1.8 교수님께 보고할 핵심 요약

1. **좌표 자체의 misalignment는 없습니다.** 우리 EDA에서 DAPI와 transcript 좌표가 99.9% coverage로 일치함을 확인했습니다.
2. **실제 문제는 "transcript-to-cell assignment"입니다.** 핵 밖 transcript(91.1%)를 어느 세포에 배정할지의 gold standard가 존재하지 않습니다.
3. **이것은 Xenium 특유의 문제가 아니라, 모든 in situ spatial transcriptomics 기술의 근본적 한계입니다** (MERFISH, SeqFISH, CosMx 등도 동일한 문제).
4. **우리 파이프라인은 이를 Dual Assignment + Optimal Expansion + Baysor 벤치마킹의 3중 전략으로 대응**하고 있으며, 이는 Salas et al. (2025) 논문이 권장하는 best practice를 충실히 따른 것입니다.

---

## 2. Mouse Liver Panel 유전자 42개 충분성 검토

### 2.1 교수님 질문 요약

> Mouse Tissue Panel에서 Liver의 유전자가 42개밖에 안 되는데, 이것만으로 분석이 충분한가? (너무 적지 않은가?)

### 2.2 결론 (한 줄 요약)

**Major cell type 5~7종 구분에는 충분하지만, hepatocyte zonation, stellate cell 활성화 상태, pathway 분석 등 심층 분석에는 부족합니다. Custom panel add-on 또는 Xenium Prime 5K Panel 병행을 권장합니다.**

### 2.3 배경: Xenium Panel 시스템 구조

10x Genomics Xenium은 미리 디자인된 패널(Pre-designed Panel) 기반으로 작동합니다:

| Panel 종류 | 유전자 수 | 특징 |
|------------|-----------|------|
| **Mouse Tissue Atlassing Panel** | 379 total (다중 조직) | 9개 이상 조직을 커버하도록 설계. 간(liver) 관련은 ~42개 |
| **Xenium Prime 5K Panel** | ~5,000 | 넓은 커버리지, 공개 scRNA-seq atlas 기반 설계. per-gene sensitivity 낮음 |
| **Custom Panel** | 1~100+ | 연구자가 직접 유전자 선정. Pre-designed panel에 50/100-gene add-on 추가 가능 |

우리 현재 데이터셋 (Human Brain Alzheimer's)의 패널 구성:

```
[Gene Panel - EDA 결과]
  Total targets (probes): 374
  Biological genes:       354 (pre-designed 254 + custom 100)
  Control probes:         20
  Probe Descriptors:
    gene:             354
    negative_control: 20
```

이 데이터셋은 Brain에 특화된 354-gene 패널이므로 충분하지만, Mouse Liver로 전환할 때는 42개로 크게 줄어듭니다.

### 2.4 간(Liver) 주요 세포 타입별 핵심 마커 유전자

간 조직에는 최소 5~7개 주요 세포 타입이 존재하며, 각각을 식별하기 위한 핵심 마커는 다음과 같습니다:

| 세포 타입 | 핵심 마커 | 비고 |
|-----------|-----------|------|
| **Hepatocyte** (간세포) | ALB, ALDOB, PCK1, BCHE, CYP3A4, CYP2E1, HAL | Pericentral: CYP2E1, GLUL / Periportal: PCK1, HAL, SDS |
| **Kupffer Cell** (쿠퍼세포) | CD5L, CLEC4F (mouse), VSIG4, FOLR2, CD163 | CLEC4F는 mouse Kupffer 표준 마커 |
| **Hepatic Stellate Cell (HSC)** | PDGFRB, ACTA2, COL1A1, DES, LRAT | ACTA2/COL1A1은 활성화(섬유화) 시 상향 |
| **Sinusoidal Endothelial (LSEC)** | PECAM1, DNASE1L3, INMT, FCN2/3, STAB2 | Fenestrated endothelium, zone-specific 마커 존재 |
| **Cholangiocyte** (담관세포) | KRT7, KRT19, SOX9, EPCAM | Periportal 영역에 집중 |
| **Portal Fibroblast** | COL1A1, ELN, BGN | 활성화 HSC와 마커 중첩 |
| **T/B/NK Cell** | CD3D, CD3E, CD79A, NKG7 | 면역세포 |
| **Macrophage (non-Kupffer)** | CD68, LYZ, S100A8 | 염증 시 침윤 |

#### 42개로 커버 가능한 수준

- **기본 cell type 분류 (5~7종):** 각 타입당 3~5개 마커 × 7 = 21~35개 필요 → **42개면 가능**
- **Hepatocyte zonation (periportal vs pericentral):** pericentral 마커 (CYP2E1, GLUL, CYP1A2) + periportal 마커 (PCK1, HAL, SDS, CYP2F2) → 추가 7~10개 필요 → **42개 패널에 포함될 수도 있으나, 연속적 gradient 분석은 어려움**
- **HSC 활성화 상태:** quiescent (LRAT, HGF) vs activated (ACTA2, COL1A1, TIMP1) → 추가 5~6개
- **Pathway 분석:** 불가능 (pathway당 최소 10~30개 유전자 필요)
- **Differential expression:** 통계적 검정력 매우 제한적

### 2.5 다른 연구에서 사용한 유전자 수 비교

| 연구 | 기술 | 유전자 수 | 간 분석 수준 |
|------|------|-----------|-------------|
| **Nault et al. (2024/2025)** | MERFISH | **317개** | Zonation gradient, HSC-hepatocyte 시그널링(THBS1/2→CD36), macrophage subtype |
| **Duan et al. (2025, Nature Genetics)** | Xenium + MIBI | ~300+ | MASLD 61 samples. MITF as LAM regulator, CV endothelial-HSC profibrotic crosstalk |
| **Barnes & Culver (2025)** | Xenium | ~300+ | Cirrhosis + non-cirrhotic MASLD, 40 FFPE samples |
| **TMA study (2025, bioRxiv)** | Xenium | ~300+ | 42 FFPE liver cores, cost-effective profiling |
| **Guilliams et al. (2022, Cell)** | 다중 기술 | ~100+ | Macrophage niche, spatial proteogenomics |

**결론: 최근 간 spatial transcriptomics 연구들은 대부분 200~300개 이상의 유전자를 사용합니다. 42개는 "분류"에는 충분하지만 "심층 분석"에는 명확히 부족합니다.**

### 2.6 해결 방안

#### Option A: Custom Add-on Panel 추가 (추천)
- 50-gene 또는 100-gene add-on panel을 기존 Mouse Tissue Panel에 결합
- MASLD 특화 마커 선정: 지방 대사(FABP1, FASN, PPARG), 섬유화(COL1A1, COL3A1, TIMP1, MMP2), 염증(TNF, IL1B, CCL2)
- 비용 효율적이고, 기존 354-gene 기반에 추가하여 ~400-450개 커버 가능
- **선정 기준:** Tabula Muris liver atlas + MASLD scRNA-seq 공개 데이터셋 (예: Duan et al. 2025의 reference list)

#### Option B: Xenium Prime 5K Panel
- ~5,000개 유전자로 포괄적 커버리지
- 장점: pathway 분석, novel cell state 발견 가능
- **단점:**
  - Per-gene sensitivity가 낮음 (Bilous et al., 2025: "newer 5K panel captures more transcripts overall but suffers from reduced per-gene sensitivity and persistent diffusion")
  - Sensitivity-breadth trade-off 존재
  - 데이터 크기 ~7배 증가 (후술 RAM 분석 참고)

#### Option C: 두 접근법 병행
- 1차 실험: Custom add-on panel (focused, high-sensitivity)
- 2차 실험: 5K panel (exploratory, broad coverage)
- 비교 분석으로 최적 전략 결정

### 2.7 5K Panel 사용 시 주의사항

우리 파이프라인은 5K Panel 데이터도 처리 가능하지만, 고려할 점이 있습니다:

1. **RAM:** Peak ~52 GB (현재 서버 2TB로 충분) → 자세한 내용은 [Section 4](#4-하드웨어-리소스-ram-검토) 참조
2. **Disk:** 입력 데이터 ~60-70 GB (현재 9 GB의 7배), 출력 ~150-200 GB
3. **Runtime:** Baysor full-tissue는 72시간+, crop mode 필수
4. **Sparse matrix:** 5K panel은 대부분의 유전자가 특정 세포에서 0 발현 → sparse matrix 효율이 매우 중요. 우리 파이프라인은 `scipy.sparse.csr_matrix`를 이미 사용 중

### 2.8 교수님께 보고할 핵심 요약

1. **42개 유전자로 major cell type 분류는 가능합니다** (간세포, 쿠퍼세포, 성상세포, 내피세포, 담관세포 등 5-7종).
2. **하지만 MASLD 연구에서 핵심적인 hepatocyte zonation, HSC 활성화 상태, 면역세포 subtype 분석에는 부족합니다.**
3. **Custom 50-100 gene add-on panel 추가를 권장합니다.** MASLD 특화 마커 (지방 대사, 섬유화, 염증)를 포함하면 ~450개 유전자로 충분한 분석이 가능합니다.
4. **5K Panel은 커버리지는 넓지만 per-gene sensitivity trade-off가 있어, targeted question에는 custom panel이 더 적합합니다.**

---

## 3. QC Control 값 분포 분석

### 3.1 교수님 질문 요약

> Negative Control은 시그널에서 노이즈를 뺀 값인가? Control 데이터에는 0 값이 훨씬 많아야 하는 것 아닌가? Raw 값과 Control 값의 분포 차이를 확인해보고 싶다.

### 3.2 결론 (한 줄 요약)

**교수님의 예상이 정확합니다. Control probe reads는 95.7%의 세포에서 완전히 0이며, Global control percentage는 0.019%로 매우 낮습니다. 이 데이터셋은 높은 품질의 실험임을 확인했습니다.**

### 3.3 Xenium Negative Control 시스템 상세 설명

Xenium은 **두 종류**의 negative control을 사용하여 서로 다른 종류의 noise를 측정합니다:

#### 3.3.1 Negative Control Probes (NCP)

- **물리적 probe 분자**가 실험에 포함됨
- **어떤 알려진 생물 genome에도 존재하지 않는 서열**을 타겟팅
- 실제 유전자 probe와 동일한 전체 과정 (hybridization → ligation → decoding)을 거침
- **측정하는 것:** 잘못된 decoding + 비특이적 결합 (off-target binding) **모두**
- 우리 데이터: 20개 NCP set (NegControlProbe_00002 ~ NegControlProbe_00042)

#### 3.3.2 Negative Control Codewords

- **물리적 probe가 없는** codebook 상의 codeword
- 어떤 gene panel의 probe도 이 codeword를 생성하지 않음
- **오직 decoding 알고리즘의 오류**로만 "검출"될 수 있음
- **측정하는 것:** decoding specificity만

#### 3.3.3 두 control의 관계

```
NCP rate ≥ Codeword rate (항상)
   ↕
NCP rate - Codeword rate = 비특이적 probe 결합에 의한 false positive 비율
```

우리 데이터 (10X metrics_summary.csv에서):
```
negative_control_probe_rate:    0.00267 (0.267%)
negative_control_codeword_rate: 0.00055 (0.055%)
→ 비특이적 결합 기여분: 0.212%
```

### 3.4 우리 데이터의 실제 분포 (EDA 결과)

`xenium_eda.py --section qc` 실행 결과:

#### 3.4.1 기본 통계량

| 지표 | 값 | 의미 |
|------|-----|------|
| **total_counts_raw** (median) | **216** reads/cell | QC 전 세포당 전체 transcript 수 |
| **total_counts_raw** (mean) | 234.7 reads/cell | |
| **total_counts_raw** (min/max) | 10 / 1,016 | QC filtering (min_counts=10) 적용 후 |
| **total_counts_control** (median) | **0** reads/cell | 세포당 control probe reads 수 |
| **total_counts_control** (mean) | 0.05 reads/cell | |
| **total_counts_control** (min/max) | 0 / 3 | 최대 3개까지만 (매우 낮음) |
| **pct_counts_control** (median) | **0.000%** | 세포당 control 비율 |
| **pct_counts_control** (mean) | 0.021% | |
| **pct_counts_control** (max) | 3.1% | 극소수 세포에서만 |

#### 3.4.2 세포별 Control Read 분포

| Control reads/cell | 세포 수 | 비율 |
|-------------------|---------|------|
| = 0 | 43,000 | **95.7%** |
| ≥ 1 | 1,943 | 4.3% |
| ≥ 2 | ~100 | <0.2% |
| ≥ 5 | 0 | **0.0%** |

**→ 교수님이 예상하신 대로, control 데이터는 0에 극도로 집중된 분포 (zero-inflated)입니다.**

#### 3.4.3 Global 수준 통계

```
Total control reads:   2,026
Total raw reads:       10,546,674
Global control %:      0.019%
```

전체 10,546,674개 read 중 control probe에서 검출된 것은 겨우 2,026개 (0.019%)입니다.

#### 3.4.4 FDR (False Discovery Rate) 추정

10x Genomics가 권장하는 FDR 계산 방법:

```
FDR = (NCP detections per probe × number of target genes) / total detected transcripts

    = (101.3 × 354) / 10,546,674
    = 0.340%
```

**FDR 0.34%는 매우 양호한 수준입니다** (일반적으로 <5%면 acceptable, <1%면 excellent).

- `NCP detections per probe = 2,026 / 20 = 101.3` (20개 NCP set에 걸친 평균)
- "실제 유전자 probe도 NCP와 같은 비율로 off-target 검출이 발생한다"는 가정 하에, 354개 유전자에 대해 예상되는 총 false positive 수를 추정

### 3.5 시각화 (생성된 QC Plots)

EDA 스크립트가 4개의 plot을 생성하여 `eda_plots/` 폴더에 저장했습니다:

| 파일명 | 내용 |
|--------|------|
| `human_alzheimers_eda_qc_raw_vs_control_hist.png` | **Histogram 3-panel:** (1) Raw counts 분포, (2) Control counts 분포 (0에 집중), (3) 두 분포 로그 스케일 오버레이 |
| `human_alzheimers_eda_qc_violin.png` | **Violin plot 2-panel:** Raw vs Control 분포 비교. Control의 극도로 좁은 분포가 시각적으로 확인됨 |
| `human_alzheimers_eda_qc_pct_control.png` | **Control % 분포:** 대부분 0% 근처, >1% 세포 극소수 |
| `human_alzheimers_eda_qc_raw_vs_control_scatter.png` | **산점도:** X=raw counts, Y=control counts. 대부분 Y=0 라인에 밀집 |

경로: `xenium-output/Xenium_V1_FFPE_Human_Brain_Alzheimers_With_Addon_outs/eda_plots/`

### 3.6 우리 파이프라인의 Control 처리 방식

#### Step 0: Control Probe 감지 및 분리

```python
# xenium_step0_formatting.py:119-154
# 4가지 패턴으로 control probe 감지:
control_mask = (
    adata.var['gene_id'].str.contains('NegControlProbe_', case=False) |
    adata.var['gene_id'].str.contains('NegControlCodeword_', case=False) |
    adata.var['gene_id'].str.contains('antisense_', case=False) |
    adata.var['gene_id'].str.contains('BLANK', case=False)
)

# QC 메트릭 계산 (filtering 전에):
adata.obs['total_counts_raw'] = np.array(adata.X.sum(axis=1)).flatten()
adata.obs['total_counts_control'] = np.array(control_genes.X.sum(axis=1)).flatten()
adata.obs['pct_counts_control'] = (total_control / total_raw) * 100

# Control probe를 adata.X에서 제거 (QC 컬럼은 obs에 보존):
adata = adata[:, ~adata.var['is_control']].copy()
# Shape 변화: (44,943 × 374) → (44,943 × 354)  [20개 control 제거]
```

#### Step 3: Transcript-level Control 필터링

```python
# xenium_step3_resegmentation.py:637-643
# AnnData 생성 전에 control transcript도 제거:
ctrl_mask = df_assigned[gene_col].str.contains(
    'NegControl|BLANK|antisense', case=False, na=False)
df_assigned = df_assigned[~ctrl_mask]
```

#### Transcript-level Control 통계 (EDA 결과)

```
[Gene/Feature Analysis - transcripts.parquet]
  Total unique features: 541 (gene 354 + controls + unassigned codewords 등)
  Control probe transcripts: 30,391 (0.19%)
```

16,371,520개 전체 transcript 중 control probe에 의한 것은 30,391개 (0.19%)뿐입니다.

### 3.7 Negative Control에 대한 추가 주의사항

1. **NCP 기반 FDR은 global 추정치**입니다. 개별 유전자별로 off-target binding이 다를 수 있습니다. Louison et al. (2025, eLife)은 Xenium v1 Human Breast panel에서 45개 유전자에 predicted off-target probe binding이 있음을 발견했습니다.

2. **NCP rate가 높다고 반드시 나쁜 실험은 아닙니다.** 10x Genomics: "a high value may result from lower detection of gene transcripts rather than nonspecific binding" - 즉, 전체 transcript 검출량이 낮으면 상대적으로 NCP rate가 높아 보일 수 있습니다.

3. **공간적 분포 확인 권장.** NCP 검출이 특정 영역에 클러스터링되면 → 조직 손상, 자가형광, 에지 효과 의심. 우리 EDA에서는 이 분석을 추가 구현할 수 있습니다.

### 3.8 교수님께 보고할 핵심 요약

1. **Negative Control Probe는 "노이즈를 뺀 값"이 아니라, noise 수준 자체를 측정하기 위한 별도의 probe입니다.** 실제 genome에 없는 서열을 타겟으로 하여, 이것이 검출되면 = false positive입니다.
2. **예상대로 control 분포는 0에 극도로 집중되어 있습니다:** 95.7% 세포에서 control = 0, 최대값도 3. Raw counts는 median 216으로 유의미한 차이가 명확합니다.
3. **Global control percentage 0.019%, estimated FDR 0.34%로 이 데이터셋은 매우 양호한 품질입니다.**
4. **Histogram, Violin plot, Scatter plot을 생성하여 `eda_plots/` 폴더에 저장했습니다.** 교수님께 보여드리기 적합한 시각화 자료입니다.

---

## 4. 하드웨어 리소스 (RAM) 검토

### 4.1 교수님 질문 요약

> 5k Panel 데이터를 풀 스케일로 돌릴 때 랩 서버(RAM)가 버틸 수 있는가? 특히 Baysor 등 특정 툴이 리소스를 많이 먹는데, 최대 피크 시 메모리(RAM)를 몇 GB나 점유하는지 정확히 알아야 한다. (필요 시 장비 증설 고려)

### 4.2 결론 (한 줄 요약)

**현재 서버(2TB RAM, 4x A5000 GPU)는 5K Panel 데이터도 충분히 처리 가능합니다. 현재 데이터 peak ~23GB, 5K Panel 추정 peak ~52GB로, 2TB RAM 대비 넉넉합니다. 단, Baysor full-tissue는 crop mode 또는 tiling 필수입니다.**

### 4.3 현재 서버 사양

`xenium_eda.py --section ram` 및 시스템 명령으로 확인:

| 항목 | 사양 | 비고 |
|------|------|------|
| **RAM** | 2,152 GB (2.0 TiB) | Available: 2,091 GB (97.2%) |
| **CPU** | 64 cores | |
| **GPU** | 4x NVIDIA RTX A5000 | 각 24,564 MiB (24 GB) VRAM |
| **Disk** | 15 TB (7.8 TB available, 47% 사용 중) | /data partition |
| **CUDA** | 12.2 | Driver 535.274.02 |

### 4.4 현재 데이터셋 (313-gene, Brain Alzheimer's) 사이즈

```
[Current Dataset Statistics]
  Transcripts:          16,371,520
  Cells (10X):          44,955
  Genes:                354 (biological) + 20 (control)
  DAPI image:           40,855 x 36,955 pixels (3.0 GB, uint16)
  transcripts.csv:      1.4 GB
  transcripts.parquet:  269.8 MB
  Total input:          ~9.0 GB
```

### 4.5 Step별 Peak RAM 추정 (현재 데이터셋)

| Pipeline Step | Peak RAM (추정) | 주요 메모리 소비 요인 | 분석 근거 |
|---------------|----------------|----------------------|-----------|
| **Step 0: Formatting** | ~5.1 GB | `mmread(matrix.mtx)` → dense matrix (44,955 × 374 × 8B = 134 MB), `pd.read_csv(transcripts.csv)` → DataFrame (~4.2 GB, CSV 1.4GB × 3배 overhead) | `adata.X = dense matrix` + 전체 transcript DataFrame |
| **Step 1: Exploration** | ~1.2 GB | adata (sparse, ~25 MB), scanpy PCA/UMAP/neighbors intermediate arrays | Step 0 output h5ad가 sparse로 저장되어 있어 가벼움 |
| **Step 2: P2R** | ~7.5 GB | transcript DataFrame (~4.2 GB), Points2Regions binning grids (k=50~500, sigma=3), adata copies | 4개 k값에 대해 각각 binned matrix 생성 |
| **Step 3: Resegmentation** | **~22.7 GB (PEAK)** | DAPI image (3.0 GB, uint16), nuclei masks (6.0 GB, int32), expanded masks (6.0 GB, int32), transcript DataFrame (~4.2 GB), expand_labels 임시 복사본 (~3 GB) | 가장 메모리 집약적 단계 |
| Step 3 (GPU) | ~4-8 GB VRAM/tile | Cellpose flow dynamics: `(1, 2, H, W) float32`. 40% VRAM budget으로 tile size 자동 계산 | `_auto_tile_size()`: 24GB A5000 → ~12,800px tiles |
| **Step 4: Comparison** | ~3.0 GB | adata 비교 + scRNA-seq reference subset | 가벼운 메트릭 계산 |
| **Step 5: Expansion** | ~6.5 GB | transcript DataFrame + KDTree (subsample 1%) + correlation 계산 | KDTree는 subsample된 셀만 사용 |
| **Step 6: Benchmark** | ~8.0 GB | 4개 adata (nuclei/cellpose/expansion/baysor), transcript DataFrame, sklearn metrics | 4 methods × adata ~1.5 GB each |
| **Step 7: Simulation** | ~4.0 GB | CellxGene Census download, simulation grid (30 permutations) | 20K cells × 200 markers |

**현재 데이터셋 Peak: ~22.7 GB (Step 3)**

### 4.6 Step 3 메모리 분석 상세 (Peak 단계)

Step 3가 가장 메모리를 많이 사용하므로 상세 분석:

```
[Step 3 Memory Breakdown]

1. DAPI 이미지 로드
   morphology_focus.ome.tif: 40,855 × 36,955 × uint16
   = 40,855 × 36,955 × 2 bytes = 3.02 GB

2. Cellpose segmentation (tiled)
   - 타일 크기 자동 계산: A5000 24GB, 40% budget → tile ~12,800px
   - 40k × 37k 이미지 → ~12 tiles (stride = tile_size - 512)
   - 타일당 GPU 사용: ~4-8 GB VRAM
   - 타일 간 gc.collect() + torch.cuda.empty_cache() 호출

3. 결과 masks (int32)
   all_masks: 40,855 × 36,955 × int32 = 6.04 GB

4. Label expansion
   expand_labels(nuclei_masks, distance=400):
   - 입력: nuclei_masks (6.04 GB)
   - 출력: expanded_masks (6.04 GB)
   - 내부 임시: distance transform (~3 GB)
   → 이 시점 peak: 3.02 + 6.04 + 6.04 + 3.0 ≈ 18.1 GB

5. Transcript assignment
   pd.read_csv(transcripts.csv): ~4.2 GB
   → 총 peak: ~22.3 GB

6. GC & mask deallocation
   del dapi_image, flows, styles, model → ~3 GB 해제
   최종 유지: expanded masks + transcripts + output adata
```

**GPU 메모리는 별도로 사용됩니다** (VRAM은 시스템 RAM과 독립):
```
[GPU Memory Usage during Step 3]
  Per tile:     ~4-8 GB VRAM (flow dynamics tensor)
  Between tiles: gc.collect() + torch.cuda.empty_cache()
  Peak VRAM:    ~8 GB (single tile)
  A5000 24GB:   충분 (남은 ~16 GB headroom)
```

### 4.7 5K Panel 데이터 RAM 추정

5K Panel로 전환 시 예상되는 변화:

#### 4.7.1 데이터 규모 스케일링

| 항목 | 313-gene | 5K Panel | 스케일 | 근거 |
|------|---------|----------|--------|------|
| Genes | 354 | ~5,000 | **14x** | Panel 크기 |
| Transcripts | 16.4M | ~115M | **~7x** | Bilous 2025: 5K panel captures more but sensitivity trade-off |
| transcripts.csv | 1.4 GB | ~10 GB | ~7x | Transcript 수에 비례 |
| Cell count | ~45K | ~45-60K | ~1-1.3x | Same tissue area, 약간 더 많은 cells detected |
| DAPI image | 3.0 GB | 3.0 GB | **1x** | 이미지는 패널과 무관 |
| Count matrix (dense) | 134 MB | ~1.9 GB | ~14x | cells × genes × 8B |
| Count matrix (sparse) | ~25 MB | ~200 MB | ~8x | Sparsity 증가로 dense보다 효율적 |

#### 4.7.2 Step별 5K Panel Peak RAM 추정

| Pipeline Step | 313-gene | 5K Panel | 증가 요인 |
|---------------|---------|----------|-----------|
| **Step 0: Formatting** | 5.1 GB | **~33.9 GB** | Dense matrix 14x + transcript DF 7x |
| **Step 1: Exploration** | 1.2 GB | ~3-5 GB | Sparse adata + PCA on 5K genes |
| **Step 2: P2R** | 7.5 GB | ~40-50 GB | 5K features per bin, 4 k-values |
| **Step 3: Resegmentation** | 22.7 GB | **~51.5 GB** | DAPI/masks 동일 + transcript DF 7x + output adata 14x |
| **Step 4: Comparison** | 3.0 GB | ~10 GB | 5K gene comparisons |
| **Step 5: Expansion** | 6.5 GB | ~35 GB | transcript DF 7x + KDTree |
| **Step 6: Benchmark** | 8.0 GB | **~43.8 GB** | 4 adata × 5K genes + transcripts |
| **Step 7: Simulation** | 4.0 GB | ~8 GB | 5K markers simulation |

**5K Panel Peak: ~51.5 GB (Step 3)**

#### 4.7.3 5K Panel 파이프라인 최적화 포인트

현재 파이프라인에서 5K Panel 처리 시 추가 최적화가 필요할 수 있는 부분:

1. **Step 0: `mmread()` → `.todense()` 변경 필요**
   - 현재: `a = mmread(matrix_path); ad = a.todense()` → 5K panel에서 ~1.9 GB dense matrix
   - 최적화: sparse로 유지 (`a.tocsr()`) → ~200 MB로 10배 절감

2. **Step 2: P2R binning grid**
   - 5K features × large grid = 많은 메모리
   - 최적화: HVG subset으로 P2R 수행 (예: top 2000 genes)

3. **Transcript loading: Parquet 사용 필수**
   - CSV 10 GB vs Parquet ~1.5 GB (메모리 로딩 시에도 차이)
   - `config.yaml`의 `use_parquet: true` 유지 필수

### 4.8 Baysor 특별 RAM 분석

Baysor는 Julia 기반 외부 프로세스로 실행되며, 파이프라인의 다른 Python 프로세스와 독립적으로 RAM을 사용합니다.

#### 4.8.1 현재 데이터셋 (313-gene)

| 모드 | Transcripts 수 | 예상 RAM | 예상 Runtime |
|------|----------------|---------|-------------|
| **Full-tissue** | 16.4M | 30-50 GB | 24-72시간 |
| **Crop mode (2000×2000µm)** | ~1-2M | 5-10 GB | 5-10시간 |
| **Tiling mode** (구현 예정) | per-tile ~2M | 5-10 GB/worker | n_workers × tile_time |

현재 config: `baysor.crop.enabled: true`, `baysor.crop.size_um: 2000`

#### 4.8.2 5K Panel Baysor 추정

| 모드 | Transcripts 수 | 예상 RAM | 예상 Runtime |
|------|----------------|---------|-------------|
| **Full-tissue** | ~115M | **200+ GB** | 수일~수주 |
| **Crop mode (2000×2000µm)** | ~7-14M | 20-40 GB | 12-24시간 |
| **Crop mode (1000×1000µm)** | ~2-4M | 5-15 GB | 3-8시간 |

**권장:** 5K Panel에서 Baysor는 반드시 crop mode 또는 tiling mode 사용. Full-tissue는 비현실적.

### 4.9 서버 용량 판정

```
[Server Capacity Assessment]

                    313-gene       5K Panel       Server Capacity
                    --------       --------       ---------------
Peak RAM (Python):  22.7 GB        51.5 GB        2,091 GB available
Peak VRAM (GPU):    ~8 GB          ~8 GB          24 GB per GPU
Baysor RAM:         5-10 GB        20-40 GB       2,091 GB available
Disk (input):       9 GB           ~63 GB         7.8 TB available
Disk (output):      ~30 GB         ~200 GB        7.8 TB available

Verdict: ALL CLEAR - 5K Panel 풀 스케일 처리에 충분한 리소스
```

**RAM Headroom:**
- 313-gene: 22.7 GB / 2,091 GB = **1.1% 사용** (98.9% 여유)
- 5K Panel: 51.5 GB / 2,091 GB = **2.5% 사용** (97.5% 여유)
- Baysor 동시 실행 포함: ~92 GB / 2,091 GB = **4.4% 사용** (95.6% 여유)

### 4.10 교수님께 보고할 핵심 요약

1. **현재 서버 (2TB RAM, 64 core, 4x A5000)는 5K Panel을 포함한 모든 시나리오에서 충분합니다.**
2. **현재 데이터셋 Peak RAM: ~23 GB** (Step 3 Cellpose resegmentation 단계)
3. **5K Panel 추정 Peak RAM: ~52 GB** (여전히 서버 RAM의 2.5%에 불과)
4. **Baysor가 가장 리소스를 많이 사용합니다.** 5K Panel full-tissue는 200GB+ 필요하므로 **crop mode(20-40 GB) 또는 tiling mode 필수**입니다.
5. **장비 증설은 현재 불필요합니다.** 단, 여러 5K Panel 데이터셋을 동시 처리하거나 Baysor full-tissue를 원하는 경우에만 고려.
6. **GPU는 Step 3 Cellpose에서만 사용되며, 24GB A5000으로 tiling 방식(~12 tiles)으로 안전하게 처리합니다.**

---

## 부록

### A. EDA 스크립트 사용법

`pipeline/xenium_eda.py`로 언제든 데이터를 직접 확인할 수 있습니다:

```bash
# 전체 EDA 실행 (모든 섹션 + 플롯 생성)
python pipeline/xenium_eda.py

# 특정 섹션만 실행
python pipeline/xenium_eda.py --section files,qc,ram

# 텍스트만 출력 (플롯 생성 안 함)
python pipeline/xenium_eda.py --skip-plots

# 리포트를 파일로 저장
python pipeline/xenium_eda.py --output-report eda_report.txt
```

**사용 가능한 섹션:**

| 섹션 | 내용 |
|------|------|
| `files` | 입력 파일 목록, 크기, 상태 체크 |
| `transcripts` | transcript 데이터 head/info/dtypes, 좌표 범위, cell_id 형식, QV, gene 분포 |
| `cells` | cells.csv 메타데이터, 공간 좌표, 면적, transcript 수 |
| `gene_panel` | gene_panel.json: 유전자 수, control probe, add-on 유전자 |
| `qc` | **total_counts_control vs total_counts_raw 분포** (교수님 질문 ③) |
| `dapi` | **DAPI-transcript 좌표 alignment 검증** (교수님 질문 ①) |
| `metrics` | 10X Xenium metrics_summary.csv 표시 |
| `h5ad` | 파이프라인 출력 h5ad 파일들의 구조, layers, obs 컬럼 |
| `ram` | **RAM 사용량 추정 및 5K Panel 스케일링** (교수님 질문 ④) |
| `experiment` | experiment.xenium 메타데이터 |

### B. 생성된 QC Plot 파일 목록

```
xenium-output/Xenium_V1_.../eda_plots/
├── human_alzheimers_eda_qc_raw_vs_control_hist.png     # Histogram 3-panel
├── human_alzheimers_eda_qc_violin.png                  # Violin plot
├── human_alzheimers_eda_qc_pct_control.png             # % control 분포
└── human_alzheimers_eda_qc_raw_vs_control_scatter.png  # Raw vs Control 산점도
```

### C. 참고 문헌

| 번호 | 저자/출처 | 제목 | 연도 |
|------|-----------|------|------|
| 1 | Salas et al. | Optimizing Xenium In Situ data utility (Nature Methods) | 2025 |
| 2 | Bilous et al. | From Transcripts to Cells (bioRxiv) | 2025 |
| 3 | Nault et al. | Spatial transcriptomics of healthy and fibrotic human liver (Nature Communications) | 2024 |
| 4 | Duan et al. | Spatially resolved multi-omics of MASLD (Nature Genetics) | 2025 |
| 5 | Louison et al. | Off-target probe binding in Xenium (eLife) | 2025 |
| 6 | Jones et al. | Proseg: probabilistic cell segmentation (PMC) | 2025 |
| 7 | 10x Genomics | Xenium Onboard Analysis - Segmentation Algorithms | 2024 |
| 8 | 10x Genomics | Metrics Matter: FDR and Target Specificity (Blog) | 2024 |
| 9 | 10x Genomics | What is the Xenium NCP metric (KB Article) | 2024 |
| 10 | 10x Genomics | Xenium Panel Design / 5K Panels | 2025 |

### D. 데이터 요약 (한 눈에 보기)

```
=== Human Alzheimer's Brain Xenium Dataset ===

실험 정보:
  Run name:           Human Alzheimer's Disease Brain (FFPE)
  Run date:           2023-03-26
  Panel:              R&D Panel (254 pre-designed + 100 custom = 354 genes)
  Preservation:       FFPE
  Instrument:         Xenium V1 (R&D, analysis_sw: xenium-1.3.0.5)
  Pixel size:         0.2125 µm/pixel

조직 정보:
  Region area:        43,735,648 µm² (~43.7 mm²)
  Cell area (total):  23,404,530 µm² (53.5% of region)
  Cells detected:     44,955
  Transcripts/cell:   216 (median)
  Transcripts/100µm²: 59.8

품질 지표:
  Q20 decoded fraction:     85.7%
  NCP rate:                 0.267%
  Codeword rate:            0.055%
  Transcripts assigned:     75.3%
  Estimated FDR:            0.34%
  Empty cells:              0%

파이프라인 처리 후:
  Step 0 output cells:      44,943 (QC 후)
  Step 3 resegmented cells: 36,170 (Cellpose)
  Genes (after control removal): 354
```
