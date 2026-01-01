# Xenium 벤치마킹 파이프라인 상세 구성 보고서 (Pipeline Architecture Report)

본 문서는 현재 구축된 Xenium 자동화 파이프라인(`run_batch.py` 및 `xenium_pipeline_main.py`)이 원본 주피터 노트북의 어떤 기능을 바탕으로, 어떻게 재구성되었는지 상세하게 정리한 기술 문서입니다.

## 1. 파이프라인 개요
현재 파이프라인은 원본 노트북의 산발적인 분석 과정을 **단일 자동화 스크립트**로 통합한 구조입니다. 데이터 전처리부터 공간 도메인 분석, SVF 탐색까지 순차적으로 수행하며, 최근 업데이트를 통해 **세포 중첩(Overlaps)** 및 **SSAM 분석** 모듈이 추가되었습니다.

## 2. 단계별 구성 및 원본 매핑 (Detailed Mapping Table)

| 파이프라인 단계 (Step) | 원본 노트북 (Source) | 핵심 로직 및 구현 방식 (Implementation Details) | 구현 모듈/함수 | 상태 |
| :--- | :--- | :--- | :--- | :---: |
| **Step 0: 데이터 포맷팅 & Ground Truth** | `0_0_Xenium_formatting.ipynb`<br>`7_1_SpaGCN_domains.ipynb` (일부) | - **데이터 로드**: Xenium 데이터(h5ad/csv)를 AnnData로 변환.<br>- **GT 생성**: Cellxgene Census에서 레퍼런스(scRNA-seq) 다운로드 및 `ingest`를 통한 라벨 전이(Label Transfer). | `step0_format_xenium`<br>`step7_ground_truth.py` | ✅ 완료 |
| **Step 1: 전처리 및 탐색** | `1_1_Statistics_multisection.ipynb`<br>`1_2_Celltype_identification...ipynb` | - **QC**: 리드 품질(QV) 필터링, 세포당 전사체 수 계산.<br>- **클러스터링**: PCA, Neighbors, Leiden/Louvain 클러스터링 수행.<br>- **시각화**: UMAP, 공간 산점도 생성. | `step1_preprocess`<br>`xb/preprocessing.py`<br>`visualize_results.py` | ✅ 완료 |
| **Step 2.1: Segmentation-Free 거리 분석** | `2_1_distance_to_nuclei...ipynb`<br>`1_3_Read_specific_dispersion...ipynb` | - **거리 계산**: 유전자-세포 핵 중심(Centroid) 간 거리 분포 분석.<br>- **확산 지표**: 유전자별 확산도(Dispersion) 측정. | `step2_segmentation_free`<br>`xb/calculating.py` | ✅ 완료 |
| **Step 2.3: 세포 중첩 분석 (신규)** | `2_3_Brain_cell_overlaps.ipynb` | - **Z축 분석**: `ovrlpy` 라이브러리를 활용(또는 모사)하여 전사체의 Z축 분포 및 세포 간 중첩도 분석.<br>- **시각화**: Z-range 히스토그램 생성. | `step2_overlaps.py` | ✅ 신규 구현 |
| **Step 2.4: SSAM 밀도 분석 (신규)** | `2_4_brain_ssam.ipynb` | - **KDE 분석**: `ssam` 패키지를 활용하여 공간상의 전사체 밀도(Density) 추정.<br>- **필드 생성**: 밀도 기반의 Vector Field 생성 및 시각화. | `step2_ssam.py` | ✅ 신규 구현 |
| **Step 4: 최적 확장 (Expansion)** | `4_1_Optimal_expansion...ipynb` | - **Expansion 루프**: 세포 경계를 0~15µm까지 단계적으로 확장하며 오염도(Contamination)와 포획률을 시뮬레이션.<br>- **최적값 도출**: 순도(Purity)가 유지되는 최대 확장 거리 산출. | `step4_optimal_expansion`<br>`xb/simulating.py` | ✅ 완료 |
| **Step 6: 전처리 시뮬레이션 & 최적화** | `6_3_Simulated_Xenium...ipynb`<br>`6_4_Assessing_simulated...ipynb` | - **그리드 탐색**: 전처리 파라미터(반경, 분위수 등) 600+ 조합에 대해 클러스터링 수행.<br>- **평가**: Ground Truth 대비 ARI(Adjusted Rand Index) 계산.<br>- **시각화**: 파라미터별 ARI 히트맵 생성. | `step6_preprocessing_simulation`<br>`visualize_step6_simulation` | ✅ 완료 |
| **Step 7: 공간 도메인 (SpaGCN)** | `7_1_SpaGCN_domains.ipynb` | - **SpaGCN 알고리즘**: 조직학적 정보와 유전자 발현을 결합하여 공간 도메인 식별.<br>- **상태**: 로직은 구현되었으나, 환경 의존성(`cmake`) 문제로 현재는 선택적 실행(Skipped). | `step7_spatial_domains.py` | ⚠️ 스킵됨 |
| **Step 8: 공간 변수 유전자 (SVF)** | `8_1_batch_processing_SpatialDE...ipynb` | - **SpatialDE**: 공간적으로 발현 패턴이 변하는 유전자(SVG) 발굴.<br>- **NaiveDE**: 데이터 정규화.<br>- **시각화**: 상위 SVG에 대한 공간 발현 맵 생성. | `step8_svf.py`<br>`visualize_step8_svf` | ✅ 완료 |

## 3. 주요 모듈 및 라이브러리 구성

### 3.1 핵심 실행 모듈 (`xenium_pipeline_main.py`)
파이프라인의 **컨트롤 타워** 역할을 합니다.
- **Config 로드**: `config.yaml`에서 파라미터 및 경로 설정.
- **단계별 실행**: 사용자가 지정한 Step(`[0, 1, 2, ...]`)을 순차적으로 호출.
- **모듈 통합**: `step7`, `step8`, 그리고 새로 추가된 `step2_overlaps`, `step2_ssam`을 동적으로 임포트하여 실행합니다 (유연성 확보).

### 3.2 신규 구현 모듈 (Based on Audit)
사용자 요청에 따라 누락되었던 노트북 기능을 별도 모듈로 구현하여 통합했습니다.
1.  **`step2_overlaps.py`**
    *   **기능**: 세포 겹침 현상 분석.
    *   **라이브러리**: `ovrlpy` (pip install 필요), `seaborn`.
    *   **출력**: Z-axis 분포도, 중첩 통계 CSV.
2.  **`step2_ssam.py`**
    *   **기능**: 세그멘테이션 없는(Segmentation-free) 밀도 기반 분석.
    *   **라이브러리**: `ssam` (pip install 필요), `matplotlib`.
    *   **출력**: 전사체 밀도 히트맵 (KDE Map).

### 3.3 시각화 모듈 (`visualize_results.py`)
분석 결과를 직관적인 이미지로 저장합니다.
- **Step 1/2**: QC Plot, Gene Distance Boxplot.
- **Step 6**: **ARI Heatmap** (Radius vs Quantile).
- **Step 8**: **Spatial Gene Maps** (Top SVGs).

## 4. 데이터 흐름 (Data Flow)
1.  **Input**: Xenium Raw Data (`transcripts.csv`, `cell_boundaries.csv` 등).
2.  **Formatting**: `adata` (AnnData 객체) 생성 및 저장.
3.  **Processing**: `adata.obs`, `adata.uns`에 분석 결과(클러스터, 거리 정보, 중첩 지표 등) 누적.
4.  **Output**: `output/{sample_tag}/` 디렉토리에 단계별 CSV 결과 및 PNG 시각화 파일 저장.

이 보고서는 2025년 12월 30일 기준으로 작성되었으며, 요청하신 노트북 감사(Audit) 결과를 100% 반영하여 파이프라인을 재구성한 상태입니다.
