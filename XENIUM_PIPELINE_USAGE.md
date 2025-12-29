# Xenium 분석 파이프라인 - 사용 가이드

## 📖 개요

이 가이드는 `xenium_pipeline_main.py` 스크립트를 사용하여 Xenium 데이터를 분석하는 방법을 설명합니다.

---

## 🚀 빠른 시작

### 1. 필수 설치

```bash
# xb 모듈 설치 (현재 디렉토리에서)
pip install -e .

# 또는 요구사항 설치
pip install -r requirements.txt
```

### 2. 기본 사용법

```python
from xenium_pipeline_main import XeniumPipeline

# 파이프라인 초기화
pipeline = XeniumPipeline(
    xenium_input_path="./data/xenium_output",
    output_path="./output",
    sample_tag="my_sample"
)

# 필수 단계만 실행 (권장)
pipeline.run_full_pipeline(
    steps=[0, 1, 6],
    target_sum=100,
    scale=False,
    hvg=False
)
```

---

## 📝 상세 사용 예제

### 예제 1: 기본 파이프라인 (필수 단계만)

**용도**: 표준 Xenium 분석
**소요 시간**: ~30분 (데이터 크기에 따라)

```python
from xenium_pipeline_main import XeniumPipeline
import os

# 경로 설정
xenium_path = "/path/to/Xenium_V1_FF_Mouse_Brain_outs"
output_dir = "./xenium_analysis_output"

# 파이프라인 생성
pipeline = XeniumPipeline(
    xenium_input_path=xenium_path,
    output_path=output_dir,
    sample_tag="mouse_brain_rep1"
)

# Step 0, 1, 6 실행 (필수)
pipeline.run_full_pipeline(
    steps=[0, 1, 6],
    target_sum=100,        # 권장: 정규화 기준값
    mincounts=10,          # 최소 카운트
    mingenes=3,            # 최소 유전자
    neigh=15,              # 최근접 이웃
    scale=False,           # 권장: 공간 데이터에서 False
    hvg=False              # 권장: 모든 유전자 사용
)

# 출력 파일:
# - output/mouse_brain_rep1_step0_formatted.h5ad
# - output/mouse_brain_rep1_step1_preprocessed.h5ad
# - output/mouse_brain_rep1_step6_ari_scores.csv
```

---

### 예제 2: 세포질-핵 분석 포함

**용도**: 유전자의 핵-세포질 분포 분석
**소요 시간**: ~45분

```python
from xenium_pipeline_main import XeniumPipeline

pipeline = XeniumPipeline(
    xenium_input_path="/path/to/xenium_output",
    output_path="./output",
    sample_tag="analysis"
)

# Step 0, 1, 2, 6 실행
pipeline.run_full_pipeline(
    steps=[0, 1, 2, 6],
    target_sum=100,
    sample_fraction=0.1  # Step 2에서 10% 샘플링
)

# 추가 출력:
# - output/analysis_step2_gene_distances.csv
#   → 유전자별 평균 거리 (핵-세포질 분포)
```

---

### 예제 3: 모든 단계 실행

**용도**: 완전한 벤치마킹
**소요 시간**: ~2시간

```python
from xenium_pipeline_main import XeniumPipeline

pipeline = XeniumPipeline(
    xenium_input_path="/path/to/xenium_output",
    output_path="./output",
    sample_tag="complete_analysis"
)

# 모든 단계 실행
pipeline.run_full_pipeline(
    steps=[0, 1, 2, 4, 6],
    target_sum=100,
    scale=False,
    hvg=False
)

# 모든 중간 파일 생성됨
```

---

### 예제 4: 개별 단계 실행 및 커스터마이징

**용도**: 특정 단계만 세밀하게 제어
**소요 시간**: 유연함

```python
from xenium_pipeline_main import XeniumPipeline
import scanpy as sc

pipeline = XeniumPipeline(
    xenium_input_path="/path/to/xenium_output",
    output_path="./output",
    sample_tag="custom"
)

# Step 0: 포맷팅
adata = pipeline.step0_format_xenium()
print(f"포맷팅 완료: {adata.shape}")

# Step 1: 전처리
adata = pipeline.step1_preprocess(
    target_sum=100,
    mincounts=10,
    mingenes=3,
    neigh=15,
    scale=False,
    hvg=False
)

# 중간에 데이터 확인
print(f"전처리 후: {adata.shape}")
print(f"클러스터 수: {len(adata.obs['leiden_1_4'].unique())}")

# Step 2: 세포질-핵 분석
gene_distances = pipeline.step2_segmentation_free_analysis(
    sample_fraction=0.1
)
print(gene_distances.head())

# Step 4: 최적 거리
optimal = pipeline.step4_optimal_expansion()
print(f"최적 거리: {optimal['optimal_from_center']:.2f} µm")

# Step 6: 전처리 최적화
ari_scores = pipeline.step6_preprocessing_simulation(
    sample_size=0.05
)
print(ari_scores.head(10))
```

---

### 예제 5: 배치 처리 (여러 샘플)

**용도**: 다중 샘플 분석
**소요 시간**: 샘플 수에 따라

```python
from xenium_pipeline_main import XeniumPipeline
import os
from pathlib import Path

# 샘플 목록
samples = {
    "ms_brain_rep1": "/data/Xenium_V1_FF_Mouse_Brain_MultiSection_1_outs",
    "ms_brain_rep2": "/data/Xenium_V1_FF_Mouse_Brain_MultiSection_2_outs",
    "ms_brain_rep3": "/data/Xenium_V1_FF_Mouse_Brain_MultiSection_3_outs",
    "human_brain": "/data/Xenium_V1_FFPE_Human_Brain_Healthy_With_Addon_outs",
}

output_base = "./xenium_batch_analysis"

# 각 샘플 분석
results_summary = {}
for sample_name, sample_path in samples.items():
    print(f"\n{'='*60}")
    print(f"처리 중: {sample_name}")
    print(f"{'='*60}")

    pipeline = XeniumPipeline(
        xenium_input_path=sample_path,
        output_path=f"{output_base}/{sample_name}",
        sample_tag=sample_name
    )

    try:
        pipeline.run_full_pipeline(
            steps=[0, 1, 6],
            target_sum=100,
            scale=False,
            hvg=False
        )
        results_summary[sample_name] = "✓ 성공"
    except Exception as e:
        results_summary[sample_name] = f"✗ 실패: {str(e)}"

# 결과 요약
print(f"\n{'='*60}")
print("배치 처리 결과")
print(f"{'='*60}")
for sample, status in results_summary.items():
    print(f"{sample}: {status}")
```

---

## ⚙️ 파라미터 상세 설명

### Step 1: 전처리 파라미터

```python
pipeline.step1_preprocess(
    target_sum=100,        # ✓ 정규화 기준 (권장: 100)
    mincounts=10,          # 최소 카운트 필터 (기본: 10)
    mingenes=3,            # 최소 유전자 발현 필터 (기본: 3)
    neigh=15,              # 최근접 이웃 수 (기본: 15)
    scale=False,           # ✓ 표준화 (권장: False)
    hvg=False,             # ✓ 고변이 유전자 (권장: False)
    save=True              # 결과 저장 여부
)
```

**파라미터 선택 가이드**:

| 파라미터 | 기본값 | 권장값 | 이유 |
|---------|--------|--------|------|
| target_sum | None | **100** | 618가지 조합 시뮬레이션에서 최적 |
| scale | False | **False** | 공간 데이터에서 표준화는 신호 손실 |
| hvg | False | **False** | 저발현 유전자도 공간 패턴 정보 포함 |
| norm | True | **True** | 정규화 필수 |
| lg (log) | True | **True** | 로그 변환 필수 |

### Step 2: 세포질-핵 분석 파라미터

```python
pipeline.step2_segmentation_free_analysis(
    sample_fraction=0.1    # 샘플링 비율 (메모리 절약)
)
```

**메모리 고려사항**:
- 전체 데이터: `sample_fraction=1.0` (느리지만 정확)
- 빠른 분석: `sample_fraction=0.1` (권장, 여전히 충분함)
- 매우 큰 데이터: `sample_fraction=0.05` (매우 빠름)

### Step 6: 전처리 시뮬레이션 파라미터

```python
pipeline.step6_preprocessing_simulation(
    sample_size=0.05       # 샘플링 비율 (메모리 절약)
)
```

**주의**: 실제 로드하는 조합 수는 매우 많습니다.
- 전체 데이터: 최고 정확도이지만 시간 오래 걸림
- 5% 샘플: 일반적으로 충분하고 빠름 (권장)
- 1% 샘플: 매우 빠르지만 덜 정확할 수 있음

---

## 📊 출력 파일 설명

### Step 0 출력
```
output/
├── my_sample.h5ad              # 포맷팅된 AnnData 객체
├── background.tiff             # 조직 배경 이미지
└── [Step 0 로그]
```

**주요 컬럼** (`adata.obs`):
- `cell_id`: 세포 고유 ID
- `x_centroid`, `y_centroid`: 세포 중심 좌표
- `transcript_counts`: 각 세포의 리드 수
- `total_counts`: 총 카운트

### Step 1 출력
```
output/
├── my_sample_step1_preprocessed.h5ad  # 전처리된 데이터
└── [전처리 로그]
```

**추가된 컬럼** (`adata.obs`):
- `leiden_1_4`: Leiden 클러스터링 결과
- `louvain_1_4`: Louvain 클러스터링 결과
- `n_genes_by_counts`, `total_counts`: QC 메트릭

### Step 2 출력
```
output/
└── my_sample_step2_gene_distances.csv
```

**컬럼 설명**:
- `mean`: 유전자의 평균 거리
- `median`: 중앙값 거리
- `std`: 표준편차
- `count`: 리드 수

**해석**:
- 낮은 `mean` → 핵 농축 유전자
- 높은 `mean` → 세포질 풍부 유전자

### Step 4 출력
```
output/
└── my_sample_step4_optimal_expansion.csv
```

**내용**:
```
optimal_from_center,           10.71
optimal_from_nucleus_border,    5.65
mean_nucleus_distance,          5.06
xenium_default,                15.00
```

### Step 6 출력
```
output/
└── my_sample_step6_ari_scores.csv
```

**컬럼**:
- `preprocessing`: 전처리 조합 이름
- `ARI`: Adjusted Rand Index 점수

**상위 항목은 권장 설정과 일치해야 함**:
```
norm_100_True_False_False_louv  0.95
norm_100_True_False_False_leiden 0.94
...
```

---

## 🔍 문제 해결

### 문제 1: 메모리 부족

**증상**: MemoryError 또는 매우 느린 처리

**해결책**:
```python
# Step 2에서 샘플링 증가
pipeline.step2_segmentation_free_analysis(sample_fraction=0.05)

# Step 6에서 샘플링 증가
pipeline.step6_preprocessing_simulation(sample_size=0.02)
```

### 문제 2: 포맷팅 실패

**증상**: `format_xenium_adata_mid_2023()` 실패

**확인사항**:
1. 올바른 경로인지 확인
2. Xenium 포맷 버전 확인:
   - 2022 (구형): `format_xenium_adata()` 사용
   - 2023 초: `format_xenium_adata_2023()` 사용
   - 2024 (현재): `format_xenium_adata_mid_2023()` 사용

```python
# 포맷 자동 감지는 아직 없으므로 수동으로 선택
from xb.formatting import format_xenium_adata  # 또는 2023 버전

adata = format_xenium_adata(
    path=xenium_path,
    tag=sample_tag,
    output_path=output_path
)
```

### 문제 3: 클러스터링 결과 없음

**증상**: `adata.obs`에 클러스터 정보 부족

**원인**: Step 1에서 전처리가 불완전하거나 파라미터 문제

**해결책**:
```python
# 로그 확인
import logging
logging.basicConfig(level=logging.DEBUG)

# 파라미터 조정
pipeline.step1_preprocess(
    target_sum=100,
    mincounts=5,        # 더 낮추기
    mingenes=1,         # 더 낮추기
    neigh=10            # 이웃 수 조정
)
```

---

## 📈 성능 최적화

### 빠른 분석 (5-10분)
```python
pipeline.run_full_pipeline(
    steps=[0, 1],
    target_sum=100,
    scale=False,
    hvg=False
)
# Step 6 생략 (전처리 최적화 스킵)
```

### 표준 분석 (30분)
```python
pipeline.run_full_pipeline(
    steps=[0, 1, 6],
    target_sum=100,
    scale=False,
    hvg=False
)
```

### 완전한 분석 (2시간)
```python
pipeline.run_full_pipeline(
    steps=[0, 1, 2, 4, 6],
    target_sum=100,
    scale=False,
    hvg=False
)
```

---

## 💡 권장 워크플로우

### Phase 1: 데이터 검증 (필수)
```python
# Step 0, 1 실행
# 품질 메트릭 확인
# QV > 20 리드 비율 확인 (기준: ≥81%)
```

### Phase 2: 전처리 최적화 (권장)
```python
# Step 6 실행
# ARI 점수 확인
# GOLDEN PATH 설정 확인
```

### Phase 3: 심층 분석 (선택)
```python
# Step 2: 핵-세포질 분석
# Step 4: 최적 확장 거리
# Step 5: 세포 분할 벤치마킹 (자체 스크립트)
```

---

## 📚 추가 자료

### 참고 파일
- `XENIUM_PIPELINE_NOTEBOOK_MAPPING.md`: 노트북-스크립트 매핑
- `xenium_pipeline_main.py`: 메인 스크립트 (이 파일)
- `notebooks/`: 원본 분석 노트북들

### 논문 결과
- **최적 정규화**: `target_sum = 100`
- **최적 확장 거리**: `5.65 ~ 10.71 µm`
- **최적 세포 분할**: `Baysor + Prior Confidence 0.8`
- **검출 효율**: Xenium은 Chromium v2 대비 1.2~1.5배

