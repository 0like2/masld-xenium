# Xenium 분석 파이프라인 - 노트북 파일 매핑

## 📊 Step별 참고 노트북 정리

### Step 0: 데이터 포맷팅

| 단계 | 노트북 파일 | 경로 | 주요 함수 | 선택사항 |
|------|-----------|------|----------|----------|
| **0.0** | `0_0_Formatting xenium to anndata` | `notebooks/0_formatting/` | `format_xenium_adata()` | 포맷 버전 선택 |
| | | | `format_xenium_adata_2023()` | (2023년 초) |
| | | | `format_xenium_adata_mid_2023()` | (2024 권장 ✓) |
| | | | `format_background()` | 배경 이미지 추출 |

**입력**: Xenium 원본 파일 (CSVs, MTX)
**출력**: `{tag}.h5ad` (AnnData 객체)

---

### Step 1: 데이터셋 탐색 & 전처리

| 단계 | 노트북 파일 | 경로 | 주요 함수 | 선택사항 |
|------|-----------|------|----------|----------|
| **1.1** | `1_1_Statistics_all_samples_using_txsim` | `notebooks/1_datasets_exploration/` | `all_quality_metrics()` | 샘플 필터링 |
| **1.6** | `1_6_Batch_preprocessing_real_Xenium_datasets` | | `keep_nuclei_and_quality()` | QV > 20 기준 |
| | | | `main_preprocessing()` | target_sum=100 ✓ |
| | | | `allcombs()` | 618가지 조합 |

**입력**: `{tag}.h5ad` (Step 0)
**출력**: 전처리된 AnnData, 품질 메트릭 CSV

**핵심 파라미터**:
```python
main_preprocessing(
    adata,
    target_sum=100,      # ✓ 권장
    norm=True,           # ✓ 정규화
    lg=True,             # ✓ 로그 변환
    scale=False,         # 권장 (False)
    hvg=False            # 권장 (False)
)
```

---

### Step 2: 세포 분할 없이 분석 (Segmentation-Free)

| 단계 | 노트북 파일 | 경로 | 주요 함수 | 선택사항 |
|------|-----------|------|----------|----------|
| **2.1** | `2_1_batch_processing_distance_to_nuclei_across_samples` | `notebooks/2_segmentation_free_analysis/` | `dispersion()` | 리드 거리 계산 |
| | | | `dist_nuc()` | 핵까지 거리 |
| **2.4** | `2_4_brain_ssam` | | SSAM 분석 | 세포 자동 식별 |

**입력**: Step 0 (미처리) + Step 1 (처리된 데이터)
**출력**:
- 거리-발현 데이터 CSV
- 핵 vs 세포질 분류
- 통계 그래프 (boxplot, ECDF)

**핵심 메트릭**:
- 리드-세포중심 거리
- 핵 겹침 여부 (overlaps_nucleus)
- 유전자별 평균 거리

---

### Step 3: 플랫폼 간 비교 (선택사항)

| 단계 | 노트북 파일 | 경로 | 주요 함수 | 선택사항 |
|------|-----------|------|----------|----------|
| **3.1** | `3_1_resegmentation_notebooks/` | `notebooks/3_techniques_comparison/` | `cellpose_segment()` | 재분할 알고리즘 |
| **3.2** | `3_2_cell_to_domain_assignment_notebooks/` | | 도메인 할당 | ROI별 분석 |
| **3.3** | `3_3_efficiency_between_methods` | | `combine_med()` | 효율 계산 |
| | | | `median_calculator()` | |
| **3.5** | `3_5_Computing_positivity_after_preprocessing_for_all_ST_techs` | | 양성률 계산 | |

**입력**: 여러 플랫폼 데이터 (CosMx, MERFISH, Visium 등)
**출력**: 플랫폼 간 검출 효율 비교 그래프

---

### Step 4: 최적 세포 확장 거리

| 단계 | 노트북 파일 | 경로 | 주요 함수 | 선택사항 |
|------|-----------|------|----------|----------|
| **4.1** | `4_1_Optimal_expansion_multisection` | `notebooks/4_optimal_expansion/` | `dist_nuc()` | 거리 범위 설정 |
| | | | cKDTree 근접 이웃 | |
| | | | 상관계수 계산 | 임계값 |

**입력**: Step 2 (거리 정보)
**출력**:
- 최적 확장 거리 (µm)
- 세포 유형별 최적값 Table
- 거리-상관계수 곡선

**결론**:
```
Xenium 기본값: 15 µm (너무 큼)
최적값: 5.65 ~ 10.71 µm
권장: 중심에서 10.71 µm 또는 핵경계에서 5.65 µm
```

---

### Step 5: 세포 분할 벤치마킹 (선택사항)

| 단계 | 노트북 파일 | 경로 | 주요 함수 | 선택사항 |
|------|-----------|------|----------|----------|
| **5.1** | `5_1_Compare_Clustering on_different_segmentations` | `notebooks/5_segmentation_benchmark/` | ARI, NMI, NMP 메트릭 | 알고리즘 선택 |
| | | | | Prior Confidence |

**입력**: 다양한 세포 분할 결과
**출력**: 분할 알고리즘 성능 비교

**권장 설정**:
```python
# Baysor + Prior Confidence = 0.8 (최고 성능)
- NMP 점수 최고
- 핵 정보 신뢰도 높음 (가중치 0.8)
```

---

### Step 6: 전처리 시뮬레이션

| 단계 | 노트북 파일 | 경로 | 주요 함수 | 선택사항 |
|------|-----------|------|----------|----------|
| **6.1** | `6_1_extract_scRNAseq_from_Census_cellxgene` | `notebooks/6_simulating_preprocessing/` | scRNA-seq 추출 | 데이터셋 선택 |
| **6.2** | `6_2_extracting_characteristics_simulated_datasets` | | 시뮬레이션 파라미터 | |
| **6.3** | `6_3_Simulated_Xenium_different_preprocessing_python` | | `allcombs_simulated()` | 618가지 조합 |
| | | | `main_preprocessing()` | |
| **6.4** | `6_4_Assessing_simulated_clusters` | | `compute_vi()` | 평가 메트릭 |
| | | | `compute_fmi()` | |

**입력**: scRNA-seq 참고 데이터
**출력**: 최적 전처리 경로 (Golden Path)

**최종 권장**:
```
✓ 정규화: target_sum = 100
✓ 로그 변환: Yes
✓ 표준화: No
✓ 고변이 유전자: No
✓ 클러스터링: Leiden/Louvain
```

---

## 📁 xb 모듈 함수 위치

```python
# formatting.py
from xb.formatting import (
    format_xenium_adata,
    format_xenium_adata_2023,
    format_xenium_adata_mid_2023,
    format_background,
    keep_nuclei_and_quality
)

# preprocessing.py
from xb.preprocessing import (
    main_preprocessing,
    preprocess_adata
)

# simulating.py
from xb.simulating import (
    allcombs,
    allcombs_simulated,
    missegmentation_simulation,
    noise_adder,
    subset_of_single_cell
)

# calculating.py
from xb.calculating import (
    dispersion,
    dist_nuc,
    entropy,
    compute_vi,
    compute_fmi
)

# comparing.py
from xb.comparing import (
    combine_med,
    median_calculator
)
```

---

## 🔄 필수 vs 선택 단계

| Step | 필수 | 이유 |
|------|------|------|
| 0 | ✓ | 데이터 포맷 변환 필수 |
| 1 | ✓ | 기본 품질 검증 필수 |
| 2 | ✗ | 핵-세포질 분석 필요시만 |
| 3 | ✗ | 벤치마킹 필요시만 |
| 4 | △ | 세포 경계 최적화 권장 |
| 5 | △ | 분할 알고리즘 비교 선택 |
| 6 | ✓ | 전처리 최적화 중요 |

---

## 📝 데이터 흐름 요약

```
원본 파일
├─ transcripts.csv
├─ cell_feature_matrix/
│  ├─ matrix.mtx
│  ├─ barcodes.tsv
│  └─ features.tsv
└─ cells.csv
    ↓ [Step 0: format_xenium_adata_mid_2023()]
{tag}.h5ad
├─ .X: (n_cells × n_genes)
├─ .obs: 세포 메타데이터
├─ .var: 유전자 정보
├─ .obsm['spatial']: 공간 좌표
└─ .uns['spots']: 리드 정보
    ↓ [Step 1: main_preprocessing()]
{tag}_preprocessed.h5ad
├─ 클러스터링 정보 추가
├─ 정규화된 발현값
└─ 품질 메트릭
    ↓ [Step 2,4: dispersion(), dist_nuc()]
거리-발현 분석
├─ 유전자별 거리 정보
└─ 핵 vs 세포질 분류
    ↓ [Step 6: allcombs()]
최적 전처리 경로 결정
    ↓
최종 분석 완료
```

---

## ⚙️ 핵심 파라미터 체크리스트

### 필수 설정
- [ ] Step 0: 포맷 버전 = `mid_2023` (현재 Xenium)
- [ ] Step 1: `target_sum = 100`
- [ ] Step 1: `norm = True`, `lg = True`
- [ ] Step 1: `scale = False`, `hvg = False`

### 권장 설정
- [ ] Step 4: 거리 범위 = 5-20 µm
- [ ] Step 5: 알고리즘 = Baysor, Prior Confidence = 0.8
- [ ] Step 6: 618가지 조합으로 검증

### 선택 설정
- [ ] Step 2: 핵-세포질 분석 필요?
- [ ] Step 3: 플랫폼 간 비교 필요?
- [ ] Step 5: 세포 분할 벤치마킹 필요?

