# 노트북 vs 파이프라인 종합 감사(Audit) 테이블

이 문서는 원본 Xenium 벤치마킹 노트북들과 현재 구현된 `run_batch.py` 파이프라인을 비교 분석한 상세 감사 보고서입니다. 로직과 시각화 측면에서 검증된 부분, 부분적으로 구현된 부분, 완전히 누락된 부분을 강조합니다.

## 범례 (Legend)
- ✅ **구현 완료**: 로직과 핵심 출력이 일치함.
- ⚠️ **부분 구현**: 핵심 로직은 존재하나, 시각화나 특정 세부 지표가 누락됨.
- ❌ **구현 누락**: 현재 파이프라인에 해당 기능이 전혀 포함되지 않음.
- ⚪ **제외됨**: 사용되지 않음(deprecated), 일회성 데이터 탐색, 또는 단순 포맷팅 용도.

| 시리즈 | 노트북 이름 | 파이프라인 상태 | 로직 차이점 | 시각화 차이점 | 추천 사항 |
| :--- | :--- | :---: | :--- | :--- | :--- |
| **0. 포맷팅** | `0_0_Xenium_formatting` | ✅ **구현 완료** | 차이 없음. | 차이 없음. | - |
| **1. 데이터셋** | `1_1_Statistics_all_samples` | ✅ **구현 완료** | 차이 없음 (Step 1 QC). | 일치함 (바이올린/산점도). | - |
| | `1_2_Celltype_identification` | ✅ **구현 완료** | 차이 없음 (Step 1 클러스터링). | 일치함 (UMAP). | - |
| | `1_3_Read_specific_dispersion` | ⚠️ **부분 구현** | Step 2에서 유전자-중심 거리 계산함. 노트북은 Read 단위로 계산 후 CSV 저장하나 파이프라인은 집계함. | 일치함 (유전자 거리 플롯). | 필요시 Step 2에서 전체 Read 거리 CSV 저장 기능 추가. |
| | `1_5_...EXPANDED-cells` | ⚪ **제외됨** | "Expanded" 세그멘테이션 전용 (기본값 아님). | 해당 없음. | "Expanded" 버전 필요 시에만 고려. |
| **2. Seg-Free** | `2_1_...distance_to_nuclei` | ✅ **구현 완료** | Step 2 로직으로 커버됨. | 일치함 (박스플롯). | - |
| | `2_3_Brain_cell_overlaps` | ❌ **구현 누락** | **신호 간섭(Signal Coherence) 및 Z축 분석 누락.** `ovrlpy` 패키지 사용. | **중첩 히스토그램 누락.** | **[높은 우선순위]** `step2_overlaps.py` 구현 필요. `ovrlpy` 의존성 확인. |
| | `2_4_brain_ssam` | ❌ **구현 누락** | **SSAM / KDE 밀도 분석 누락.** `ssam` 패키지 사용. | **KDE 밀도 맵 & 벡터 필드 시각화 누락.** | **[높은 우선순위]** `step2_ssam.py` 구현 필요. `ssam` 패키지 필요. |
| **3. Technique** | `3_X_...` 시리즈 | ⚪ **제외됨** | 타 기술(CosMx, Merscope) 비교용. Xenium 전용 파이프라인 범위 밖. | - | - |
| **4. Expansion** | `4_1_Optimal_expansion` | ✅ **구현 완료** | Step 4 시뮬레이션 루프 (0-15um). | 일치함 (최적화 곡선). | - |
| **5. Seg-Bench** | `5_1_Compare_Clustering` | ⚪ **제외됨** | 세그멘테이션 방법론 비교용. | - | - |
| **6. Simulating** | `6_1_extract_scRNAseq` | ✅ **구현 완료** | Step 0/7 Ground Truth 준비 과정에서 처리됨. | - | - |
| | `6_2_extracting_characteristics`| ⚠️ **부분 구현** | 기초 통계는 Step 6 준비 과정에서 커버됨. | - | - |
| | `6_3_Simulated_Xenium...` | ✅ **구현 완료** | Step 6 최적화 로직 (618개 조합 수행). | - | - |
| | `6_4_Assessing_simulated...` | ⚠️ **부분 구현** | Step 6에서 ARI/VI/FMI 점수 계산함. | ❌ **히트맵 시각화 누락.** | **[중간 우선순위]** Step 6 결과에 히트맵 시각화 추가 필요. |
| **7. Domains** | `7_1_SpaGCN_domains` | ❌ **스킵됨** | 구현되었으나 의존성 문제로 스킵됨 (`SpaGCN`/`cmake`). | **공간 도메인 맵 누락.** | `SpaGCN` 실행을 위한 환경 수정 필요. |
| **8. SVF** | `8_1_...SpatialDE_SVF` | ✅ **구현 완료** | Step 8 SpatialDE 로직 검증됨. | ❌ **유전자 맵 누락.** | **[중간 우선순위]** 상위 SVG에 대한 공간 유전자 발현 플롯 추가. |

## 주요 "구현 누락" 모듈 상세 분석

### 1. Step 2_3: 뇌 세포 중첩 (Brain Cell Overlaps)
*   **목적**: 2D 세그멘테이션의 한계로 인한 3D 상의 세포 중첩 및 신호 간섭(Incoherence) 탐지.
*   **핵심 로직**: Z축 분포 분석, 신호 간섭 점수(Signal coherence scoring) 계산.
*   **의존성**: `ovrlpy`, 기본 과학 패키지.
*   **상태**: 파이프라인에 없음.
*   **조치**: `step2_overlaps.py` 모듈 생성 및 파이프라인 통합.

### 2. Step 2_4: SSAM (KDE 분석)
*   **목적**: 세그멘테이션 없이 커널 밀도 추정(KDE)을 이용한 세포 유형 매핑.
*   **핵심 로직**: 가이드 KDE(Guided KDE), 벡터 필드 생성, 국소 최대값(Local Maxima) 탐지.
*   **의존성**: `ssam` 패키지, `scikit-image`.
*   **상태**: 파이프라인에 없음.
*   **조치**: `step2_ssam.py` 모듈 생성. `ssam` 패키지 설치 가능 여부 확인 필요.

### 3. Step 6 시각화: 히트맵 (Heatmaps)
*   **목적**: 수백 개의 시뮬레이션 파라미터(반경 vs 분위수) 조합에 따른 ARI 점수 분포 시각화.
*   **핵심 로직**: ARI 결과 데이터프레임 피벗 -> Seaborn 히트맵 출력.
*   **상태**: 로직(ARI 계산)은 존재하나, 시각화 함수가 `visualize_results.py`에 정의만 되고 구현 내용이 없거나 호출되지 않음.
*   **조치**: `visualize_results.py` 내 `visualize_step6_simulation` 함수 구현.

### 4. Step 8 시각화: SVF 맵 (SVF Maps)
*   **목적**: 공간적으로 변이가 큰 유전자(SVG)의 발현 패턴 시각적 확인.
*   **핵심 로직**: SpatialDE로 검출된 상위 유전자들의 공간 좌표상 발현량 플로팅.
*   **상태**: 로직(SpatialDE 수행)은 존재하나, 시각화 누락.
*   **조치**: `step8_svf.py` 또는 `visualize_results.py`에 플로팅 루틴 추가.
