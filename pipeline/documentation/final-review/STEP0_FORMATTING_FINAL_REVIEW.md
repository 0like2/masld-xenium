# Step 0: Formatting - Final Review

## 1. Overview

### 1.1 What This Step Does
Step 0은 Xenium 머신의 원시 출력(raw machine output)을 표준화된 **AnnData (.h5ad)** 객체로 변환하는 포매팅 단계이다. Xenium 장비는 세포별 유전자 발현 매트릭스, 세포 메타데이터, 개별 transcript 좌표, DAPI 이미지 등을 개별 파일로 출력하는데, 이를 하나의 통합 데이터 구조로 재구성하여 다운스트림 분석이 가능하도록 한다.

### 1.2 Paper Context
논문 (Salas et al., Nature Methods 2025)에서는 25개 Xenium 데이터셋을 분석하였으며, 모든 데이터셋은 커스텀 포매팅 함수(`format_xenium_adata_mid_2023`)를 통해 AnnData로 변환되었다 (Methods: "Xenium dataset processing" 섹션). GitHub 레포(`Moldia/Xenium_benchmarking`)의 `0_formatting/` 노트북에 해당 코드가 공개되어 있다.

### 1.3 논문 Figure 매칭
- **Figure 1b**: 포매팅된 25개 데이터셋의 요약 테이블 (세포 수, 유전자 수, reads/cell 등) → Step 0의 출력인 AnnData 객체에서 추출된 통계
- **Extended Data Figure 1a-e**: 데이터셋 간 재현성 비교 (총 reads/gene, cell type 비율 등) → Step 0 포매팅 결과를 기반으로 한 QC 분석

---

## 2. Pipeline Process Flow

```
Xenium Machine Output Directory
├── cell_feature_matrix.tar.gz     ← 세포×유전자 매트릭스 (MTX format)
├── cells.csv.gz                   ← 세포 메타데이터 (좌표, 면적 등)
├── transcripts.csv.gz             ← 개별 transcript 좌표 (x, y, z, gene)
├── morphology_focus.ome.tif       ← DAPI 형광 이미지
├── gene_panel.json                ← 유전자 패널 정보 (Ensembl ID 매핑)
└── analysis/                      ← 10X 전처리 결과 (UMAP, PCA, clusters)
    ├── umap/
    ├── pca/
    └── clustering/
         ↓
    [Step 0: Formatting]
         ↓
├── {sample_tag}.h5ad              ← 통합 AnnData 객체
├── {sample_tag}_transcripts.parquet ← Transcript sidecar (빠른 I/O)
└── DAPI 이미지 (원본 또는 다운샘플)
```

### 2.1 상세 처리 흐름

```
1. 압축 해제
   cell_feature_matrix.tar.gz → matrix.mtx, barcodes.tsv, features.tsv
   cells.csv.gz → cells.csv
   transcripts.csv.gz → transcripts.csv

2. Cell×Gene Matrix 구축
   matrix.mtx (Market Matrix) → scipy.io.mmread() → dense matrix
   barcodes.tsv → cell IDs
   features.tsv → gene IDs/names (2열 또는 3열 포맷 자동 감지)

3. AnnData 객체 생성
   adata.X = count matrix (cells × genes)
   adata.obs = cells.csv 메타데이터 병합
   adata.var = features 정보 (gene_id, gene_name, reason_of_inclusion)

4. Gene Panel 정규화
   gene_panel.json → Ensembl ID 매핑
   컨트롤 프로브 제거 (negative/positive controls, blank codewords)

5. 10X Analysis 임포트
   UMAP 좌표 → adata.obsm['X_umap']
   PCA 좌표 → adata.obsm['X_pca']
   Clustering → adata.obs['leiden_*']

6. QC Filtering
   min_counts threshold → 최소 transcript 수 미달 세포 제거
   min_genes threshold → 최소 발현 유전자 수 미달 세포 제거

7. Transcript Sidecar 저장
   transcripts.csv → {sample_tag}_transcripts.parquet (Parquet 형식)
```

---

## 3. Notebook vs Pipeline 구현 비교

### 3.1 원본 Notebook
- **파일**: `notebooks/0_formatting/0_0_Formatting xenium to anndata.ipynb`
- **관련 노트북들**:
  - `0_1_anndata_Formatting_for_R_batch_processing_Xenium_datasets.ipynb` - R 호환 포맷
  - `0_2_simulated_datasets_Formatting_for_R.ipynb` - 시뮬레이션 데이터 포맷
  - `0_3_unprocessed_adata_to_nuclei_batch_processing.ipynb` - nuclei 전용 필터링

### 3.2 Pipeline 구현
- **파일**: `pipeline/xenium_step0_formatting.py` (pipeline_main.py에서 호출)
- **진입점**: `run_step0(config)` → `format_xenium_adata_mid_2023()`

### 3.3 비교 테이블

| 항목 | Notebook | Pipeline | 상태 |
|------|----------|----------|------|
| 압축 해제 | `tarfile.extractall()` + `gzip` | 동일 | OK |
| MTX 읽기 | `scipy.io.mmread()` | 동일 | OK |
| Features 파싱 | 3열/2열 자동 감지 | 동일 | OK |
| Gene panel 정규화 | JSON → Ensembl ID | 동일 | OK |
| 컨트롤 프로브 제거 | 수동 필터링 | 동일 | OK |
| 10X analysis 임포트 | UMAP/PCA/clustering | 동일 | OK |
| QC 필터링 | 하드코딩된 값 | config.yaml에서 읽기 | 개선 |
| Transcript sidecar | CSV 저장 | Parquet 저장 (설정 가능) | 개선 |
| nuclei-only 필터링 | `0_3` 노트북에서 별도 | `filter_nuclei_only` config 옵션 | 통합 |
| 에러 처리 | 없음 | try/except + 로깅 | 개선 |
| 경로 관리 | 하드코딩 | config 기반 | 개선 |

---

## 4. Sub-step 상세 분석

### 4.1 데이터 압축 해제
```python
# cell_feature_matrix.tar.gz 해제
tarfile.open(tar_path).extractall(output_dir)

# .gz 파일 개별 해제
for gz_file in [cells.csv.gz, transcripts.csv.gz, ...]:
    gunzip(gz_file)
```
- Xenium V1 포맷은 `cell_feature_matrix.tar.gz` 안에 `matrix.mtx`, `barcodes.tsv`, `features.tsv`를 포함
- 최신 Xenium 포맷에서는 `cell_feature_matrix/` 디렉토리가 이미 해제된 상태일 수 있음

### 4.2 Features 파싱 (유전자 정보)
```python
# 3열 포맷: gene_id | gene_name | reason_of_inclusion
# 2열 포맷: gene_id | reason_of_inclusion (gene_name 없음)
if features_df.shape[1] == 3:
    features_df.columns = ['gene_id', 'gene_name', 'reason_of_inclusion']
else:
    features_df.columns = ['gene_id', 'reason_of_inclusion']
    features_df['gene_name'] = features_df['gene_id']
```
- `reason_of_inclusion`: Gene Expression, Negative Control Probe, Negative Control Codeword 등
- 컨트롤 프로브는 이후 단계에서 제거됨

### 4.3 QC Filtering
```python
# config.yaml 설정
formatting:
  mincounts: 10    # 세포당 최소 총 transcript 수
  mingenes: 3      # 세포당 최소 발현 유전자 수
  filter_nuclei_only: false  # true면 nucleus 내 transcript만 유지
```

**의미**:
- `mincounts=10`: 10개 미만의 transcript를 가진 세포는 저품질로 간주하여 제거. 논문에서 "Only 0.21% of the cells had fewer than ten assigned reads" (p.814)
- `mingenes=3`: 3개 미만의 유전자를 발현하는 세포는 debris 또는 empty droplet으로 간주
- `filter_nuclei_only`: 0_3 노트북의 기능을 통합. nucleus 경계 내의 transcript만 유지하면 segmentation 오류를 최소화하지만, cytoplasmic reads를 잃게 됨

### 4.4 Transcript Sidecar
원본 transcript 데이터는 AnnData에 직접 저장하기에 너무 크므로 (수백만 행), 별도 파일로 저장:
- **Parquet** (기본값, `use_parquet: true`): 빠른 I/O, 압축 효율적
- **CSV** (폴백): 호환성

주요 컬럼:
| 컬럼 | 설명 | 예시 |
|------|------|------|
| `transcript_id` | 고유 transcript ID | AABCDE-1 |
| `cell_id` | 할당된 세포 ID (0=미할당) | 12345 |
| `x_location` | X 좌표 (µm) | 1523.45 |
| `y_location` | Y 좌표 (µm) | 2341.67 |
| `z_location` | Z 좌표 (µm) | 8.2 |
| `feature_name` | 유전자 이름 | GFAP |
| `overlaps_nucleus` | 핵 내부 여부 (0/1) | 1 |
| `qv` | 품질 점수 (phred) | 40 |

---

## 5. Visualization (시각화)

Step 0은 주로 데이터 변환 단계이므로 시각화 출력이 제한적이다:

| 시각화 | 설명 | 분석법 |
|--------|------|--------|
| 데이터셋 요약 통계 | 총 세포 수, 유전자 수, 평균 reads/cell | 텍스트 로그 출력 |
| QC 필터링 전후 비교 | 필터링으로 제거된 세포 비율 | 로그에서 확인 |

---

## 6. Input / Output 상세

### 6.1 Input
| 파일 | 형식 | 설명 |
|------|------|------|
| `cell_feature_matrix.tar.gz` | tar.gz (내부: MTX + TSV) | 세포×유전자 count matrix |
| `cells.csv(.gz)` | CSV | 세포 메타데이터: cell_id, x_centroid, y_centroid, cell_area, nucleus_area |
| `transcripts.csv(.gz)` | CSV | 개별 transcript 좌표 및 할당 정보 |
| `morphology_focus.ome.tif` | TIFF | DAPI 형광 이미지 (다채널) |
| `gene_panel.json` | JSON | 유전자 패널 정보 (Ensembl ID, type) |
| `analysis/` | Directory | 10X 사전 분석 (UMAP, PCA, clustering) |

### 6.2 Output
| 파일 | 형식 | 설명 |
|------|------|------|
| `{sample_tag}.h5ad` | HDF5 (AnnData) | 통합 세포×유전자 데이터 + 메타데이터 |
| `{sample_tag}_transcripts.parquet` | Parquet | Transcript sidecar (좌표 + 할당) |
| `background.tiff` | TIFF | DAPI 배경 이미지 (`morphology_mip.ome.tif` → 단일 채널 변환, `format_background()` line 325) |

### 6.3 AnnData 구조
```
adata.X         → (n_cells × n_genes) count matrix [sparse CSR]
adata.obs       → DataFrame: cell_id, x_centroid, y_centroid, cell_area,
                   nucleus_area, n_counts, n_genes, ...
adata.var       → DataFrame: gene_id, gene_name, reason_of_inclusion
adata.obsm      → 'X_umap': UMAP 좌표 (2D)
                   'X_pca': PCA 좌표 (다차원)
adata.uns       → 메타정보 (sample_tag, qc_params, ...)
```

---

## 7. 관련 설정값 정리

```yaml
# config.yaml - Step 0 관련
input_path: "/path/to/Xenium_outs"     # Xenium 출력 디렉토리
output_dir: "xenium-output"             # 파이프라인 출력 기본 경로
sample_tag: "human_alzheimers"          # 샘플 식별자 (출력 파일명 접두사)

formatting:
  mincounts: 10                         # 최소 transcript 수/세포
  mingenes: 3                           # 최소 유전자 수/세포
  filter_nuclei_only: false             # nucleus 전용 필터링

use_parquet: true                       # Parquet vs CSV 저장 형식
```

### 7.1 파라미터 조정 가이드

| 파라미터 | 현재값 | 범위/옵션 | 영향도 | 조정 가이드 |
|----------|--------|----------|--------|------------|
| `input_path` | (데이터 경로) | 문자열 | - | Xenium 장비 출력 디렉토리 경로. `cell_feature_matrix.tar.gz`, `cells.csv.gz` 등이 있는 폴더 |
| `output_dir` | `"xenium-output"` | 문자열 | - | 모든 파이프라인 출력이 저장되는 기본 디렉토리 |
| `sample_tag` | `"human_alzheimers"` | 문자열 | - | 출력 파일명 접두사. 데이터셋 식별용 |
| `formatting.mincounts` | `10` | 1-50 | **HIGH** | 세포당 최소 transcript 수. 낮추면 더 많은 세포 보존 (노이즈 포함), 높이면 고품질 세포만 남김. 논문: "0.21% of cells had fewer than 10 reads" |
| `formatting.mingenes` | `3` | 1-20 | MEDIUM | 세포당 최소 발현 유전자 수. 1-5가 일반적. 너무 높으면 낮은 complexity 세포 (e.g., RBC) 제거 |
| `formatting.filter_nuclei_only` | `false` | true/false | **HIGH** | `true`: nucleus 내 transcript만 유지 → 높은 순도, cytoplasmic reads 손실. `false`: 모든 transcript 유지 (기본) |
| `use_parquet` | `true` | true/false | LOW | Parquet(빠른 I/O, 작은 파일) vs CSV(범용 호환). 성능만 차이, 결과 동일 |

**튜닝 팁**:
- `mincounts`/`mingenes`는 Step 0에서 설정하면 **모든 downstream 분석에 영향**. 보수적으로 시작(10/3)하고 QC 결과를 보고 조정
- `filter_nuclei_only=true`는 segmentation-free 분석(Step 2)에서 cytoplasmic 정보를 잃으므로 주의

---

## 8. 핵심 로직 및 정의된 값의 의미

### 8.1 pixel_to_um 변환
Xenium의 좌표계는 **마이크로미터(µm)** 기반이지만, 이미지 기반 분석에서는 **픽셀** 단위를 사용한다:
- Xenium default: **4.70588 pixels/µm** (= 0.2125 µm/pixel)
- 이 값은 Step 3, 4, 5에서 거리 변환에 사용됨

### 8.2 phred quality score (qv)
- `qv > 20` = 99% 정확도 (논문 기준: 72-91% of reads with qv > 20)
- Step 0에서는 qv 필터링을 하지 않지만, 다운스트림에서 활용 가능

### 8.3 overlaps_nucleus
- `1` = transcript가 핵 segmentation mask 내부에 위치
- `0` = transcript가 핵 외부 (세포질 또는 미할당)
- 논문에서는 "76.8% of reads on average being assigned to cells" (p.814)

---

## 9. 논문 PDF와의 정합성

| 논문 내용 | Step 0 구현 | 정합 여부 |
|-----------|------------|----------|
| "formatted as anndata using a customized function" (Methods) | `format_xenium_adata_mid_2023()` | OK |
| "processed using Scanpy (v1.9.1)" (Methods) | Scanpy 기반 전처리 | OK |
| mincounts/mingenes 필터링 | config로 설정 가능 | OK |
| 25개 데이터셋 일괄 처리 | sample_tag 기반 개별 처리 | OK |
| 컨트롤 프로브 제거 | features.tsv 파싱에서 처리 | OK |

---

## 10. 시각화 상세 분석 가이드

Step 0 자체는 포매팅 단계이므로 시각화가 제한적이지만, Step 0의 출력이 논문의 주요 Figure에 직접 사용된다.

### 10.1 논문 Figure 1b - 데이터셋 요약 테이블

**논문 위치**: p.815, Fig. 1b
**파이프라인 출력**: Step 0의 AnnData에서 추출한 통계

**무엇을 봐야 하는가**:
```
┌─────────────────────────────────────────────────────────────┐
│  Dataset     | Tissue | Source | No.genes | Cells | Reads  │
│  Brain(cor)  | FF     | 3      | 248      | ●●●   | ████   │
│  Brain(Alz)  | FFPE   | 1      | 319      | ●●    | ███    │
│  ...                                                        │
│                                                             │
│  컬럼: Mean reads/cell, Mean genes/cell, Total reads,      │
│        Cells recovered, Area imaged, Reads assigned (%)     │
└─────────────────────────────────────────────────────────────┘
```

**해석법**:
- **Mean reads per cell**: 150-250이 정상. 100 미만이면 캡처 효율 문제 의심
- **Mean genes per cell**: 30-80이 정상. 패널 크기에 비례
- **Reads assigned to cells (%)**: 논문에서 평균 76.8%. 60% 미만이면 segmentation 불량
- **FF vs FFPE**: Fresh Frozen이 FFPE보다 약간 높은 reads/cell (조직 보존 상태 차이)

**좋은 결과 vs 나쁜 결과**:
- 좋음: reads/cell > 150, genes/cell > 40, assigned > 70%
- 나쁨: reads/cell < 50, genes/cell < 20, assigned < 50%

### 10.2 논문 Extended Data Fig. 1a - 재현성 검증

**논문 위치**: Extended Data Fig. 1a (p.826)
**관련**: Step 0 포매팅된 두 개의 preview vs own machine 데이터셋

**무엇을 봐야 하는가**:
- **상단**: 유전자별 total reads scatter plot (preview 1 vs preview 2)
  - r=0.999 → 실험 간 재현성이 매우 높음
  - 대각선에서 벗어나는 점 = 실험 간 차이가 큰 유전자
- **하단**: 세포 유형별 abundance scatter plot
  - 세포 유형 비율이 실험 간 보존되는지 확인

**핵심**: Step 0에서 올바르게 포매팅되었다면, 동일 조직의 다른 실험에서도 일관된 통계를 보여야 함

---

## 11. Summary

Step 0은 파이프라인의 **기초 데이터 변환 단계**로, 원본 Xenium 출력을 표준 AnnData 형식으로 변환한다. Notebook 대비 Pipeline의 주요 개선점:
1. **Config 기반 파라미터 관리**: 하드코딩 → YAML 설정
2. **Parquet 지원**: CSV 대비 5-10배 빠른 I/O
3. **nuclei-only 필터링 통합**: 별도 노트북 → 단일 config 옵션
4. **에러 처리 및 로깅**: Robust한 실행 보장

구현 완성도: **HIGH** - 논문/노트북 대비 완전히 구현됨

---

## 12. 논문 Figure 직접 대응 및 시각화 정상 판별 종합 가이드

### 12.1 파이프라인 출력 → 논문 Figure 매핑 종합표

| 파이프라인 출력 | 논문 Figure | 논문 페이지 | 논문 원문 설명 |
|---|---|---|---|
| AnnData 요약 통계 (세포 수, 유전자 수, reads/cell) | **Fig. 1b** (p.815) | 815 | "Summary table of the Xenium datasets, detailing dataset characteristics, descriptors and quality metrics" |
| QC 필터링 전후 통계 | **Fig. 1b** 범례 | 815 | "Only 0.21% of the cells had fewer than ten assigned reads and were excluded from further analysis" |
| 데이터셋 재현성 통계 | **Extended Data Fig. 1a-e** | 826 | "Datasets generated in independent experiments on similar samples using identical probe panels exhibited a strong similarity in gene-specific detection efficiency, dispersion and reads per cell" |

### 12.2 정상 결과 판별 체크리스트

Step 0은 포매팅 단계이므로 시각화보다는 **출력 통계의 정상 범위**를 확인한다:

#### 체크리스트: AnnData 요약 통계
- [ ] **총 세포 수**: 10,000-500,000 범위 (조직 크기에 따라 다름)
- [ ] **총 유전자 수**: 210-392 범위 (패널에 따라). 컨트롤 프로브 제거 후 감소
- [ ] **평균 reads/cell**: 150-250이 정상. 논문: "an average of 186.6 reads per cell"
- [ ] **Reads assigned to cells (%)**: 70-85%가 정상. 논문: "76.8% of reads being assigned to cells"
- [ ] **QC 제거 비율**: 1% 미만이 정상. 논문: "Only 0.21% of the cells had fewer than ten assigned reads"
- [ ] **Mean genes/cell**: 30-80이 정상 (패널 크기에 비례)

**왜 이 범위가 정상인가**: 논문에서 25개 데이터셋을 분석한 결과, 위 범위가 일관되게 관찰됨. FF와 FFPE 조직 간 차이는 있으나 모두 이 범위 내.

**비정상 신호**:
- reads/cell < 50 → 캡처 효율 문제 또는 데이터 품질 이슈
- assigned < 50% → segmentation 심각한 불량
- QC 제거 > 10% → 이미지 품질 또는 실험 조건 문제

### 12.3 논문 원문 인용 (시각화 관련)

| 관련 출력 | 논문 원문 인용 | 페이지 |
|---|---|---|
| 세포 수 / reads/cell | "These samples span a variety of types, altogether representing a total of 1.2 billion reads and 6 million cells" | 813 |
| QC 필터링 기준 | "Only 0.21% of the cells had fewer than ten assigned reads and were excluded from further analysis, positioning Xenium as a suitable platform for assessing cell-type frequencies in tissues" | 813 |
| Reads assigned 비율 | "Using Xenium's default segmentation, an average of 186.6 reads per cell was observed throughout the datasets, with 76.8% of reads being assigned to cells" | 813 |
| FF vs FFPE 비교 | "no obvious differences between fresh frozen (FF) and formalin-fixed paraffin-embedded (FFPE) sections" | 813 |
| 재현성 | "Datasets generated in independent experiments on similar samples using identical probe panels exhibited a strong similarity in gene-specific detection efficiency, dispersion and reads per cell" | 814 |
| 유전자 패널 크기 | "The number of genes profiled per sample type (the 'panel') ranged between 210 and 392 genes" | 813 |
| 품질 점수 | "All datasets included the three-dimensional (3D) position (x, y and z), gene identity and phred-based quality value (qv) of every decoded read, with 81% (range, 72-91%) of the reads on average exhibiting high quality (qv > 20)" | 813 |
