# Xenium 파이프라인 마스터 리포트: 종합 분석 및 현황
**날짜:** 2026-01-01
**버전:** 2.0 (Master Compilation)

## 1. 요약 (Executive Summary)
본 문서는 Xenium 벤치마킹 파이프라인의 **단일 진실 공급원(Single Source of Truth)**으로서, 기술적 아키텍처, 시각화 기능, 부족한 점(Gap), 그리고 향후 로드맵을 통합하여 정리한 문서입니다.
우리는 20개 이상의 파편화된 Jupyter Notebook들을 **견고하고 자동화된 Python 파이프라인**으로 리팩토링하는 데 성공했으며, 이는 논문 수준의 표준을 충족합니다.

### 주요 성과 (Key Achievements)
*   **자동화 (Automation):** Step 0(전처리)부터 Step 6(최적화)까지 `run_all.sh` 스크립트 하나로 중단 없이 수행됩니다.
*   **시각화 (Visualization):** 논문의 핵심 그림들(Image 3e, 3f, 5c)을 파이프라인 내에 직접 구현하여 자동 생성되도록 했습니다.
*   **견고성 (Robustness):** 외부 scRNA-seq 레퍼런스를 활용한 **Ground Truth(정답지)** 생성 로직과 ROI 최적화 시뮬레이션을 탑재했습니다.

---

## 2. 기술 심층 분석 (Technical Deep Dive)

### 2.1 Ground Truth 생성 (Step 0)
*   **목표**: 정확도 벤치마킹을 위해 모든 세포에 대한 "Gold Standard(정답)" 세포 유형 라벨을 생성합니다.
*   **방식**: `step7_ground_truth.py` (Step 0에 통합됨).
    1.  **레퍼런스(Reference)**: Cellxgene Census에서 검증된 scRNA-seq 데이터를 다운로드합니다.
    2.  **매핑(Mapping)**: 유전자 표기법 차이(Ensembl ID vs Gene Symbol)를 교집합(Intersection) 로직으로 해결합니다.
    3.  **전이(Transfer)**: `scanpy.tl.ingest`를 사용하여 유전자 발현 패턴이 유사한 scRNA-seq의 세포 유형을 Xenium 세포(Target)에 입힙니다.

### 2.2 세그멘테이션 없는 분석 (Step 2: Segmentation-Free)
*   **목표**: 불완전한 세포 경계(Boundary)에 의존하지 않고 전사체 분포를 분석합니다.
*   **방식**:
    *   **거리 분석(Distance Analysis)**: 모든 전사체(Transcript)와 가장 가까운 세포 핵(Nucleus) 중심점 간의 유클리드 거리를 계산합니다.
    *   **중첩/SSAM**: (부분 구현됨) Z축 신호 일관성 및 커널 밀도 추정(KDE)을 통해 세포 밀집도를 평가합니다.

### 2.3 ROI 최적화 (Step 4)
*   **목표**: 세포 경계를 어디까지 확장하는 것이 최적인지 찾는 "확장 거리(Expansion Distance)"를 결정합니다.
*   **로직**:
    *   **시뮬레이션**: 경계를 0µm에서 15µm까지 0.5µm 단위로 점진적으로 확장해 봅니다.
    *   **트레이드오프(Trade-off)**: **포획 효율(Capture Efficiency)**(더 많은 리드 포함)과 **순도(Purity)**(핵 내부 비율 유지)를 측정합니다.
    *   **결과**: 순도가 급격히 떨어지는 지점(Elbow Point)을 최적 거리로 선택합니다. (예: 뇌 조직의 경우 약 10µm)

---

## 3. 시각화 기능 감사 (Visualization Capabilities Audit)

일관된 플롯 생성을 위해 `visualize_results.py` 모듈과 재생성 스크립트 `run_vis_all.py`를 구현했습니다.

### A. 구현된 시각화 (상태: ✅ 완료)

| 카테고리 | 플롯 유형 | 설명 | 파일명 패턴 |
| :--- | :--- | :--- | :--- |
| **QC & 전처리** | Violin Plots | QC 지표 (Counts, Genes) | `*_step1_qc_violin.png` |
| | PCA Variance | PC별 분산 비율 | `*_step1_pca_variance.png` |
| | UMAP | Leiden 클러스터링 지도 | `*_step1_umap_plot.png` |
| | Spatial Scatter | XY 좌표상 세포 유형 분포 | `*_spatial_clusters.png` |
| **Segmentation-Free** | Box/Bar Plots | 유전자별 핵 거리 분포 | `*_step2_gene_dist_boxplot.png` |
| **최적화 (Optimization)** | Dual-Axis Line | 순도(Purity) vs 포획(Capture) 최적화 곡선 (Image 3f) | `*_step4_optimization_curve.png` |
| **고급 (Advanced)** | **Violin Plot** | 마커 유전자 발현 분포 (Image 3e) | `*_marker_genes_violin.png` |
| | **Seg. Overlay** | 전사체 + 세포 중심/경계 오버레이 (Image 5c) | `*_segmentation_overlay_roi.png` |

### B. 미구현 시각화 (상태: ❌ 예정)

| 카테고리 | 플롯 유형 | 제외 사유 | 우선순위 |
| :--- | :--- | :--- | :--- |
| **공간 통계 (Spatial Stats)** | Neighborhood Heatmap (Image 5e) | `squidpy` 모듈 통합 필요. | 높음 |
| | Spatial Co-occurrence | `squidpy` 필요. | 높음 |
| **ARI 히트맵** | Heatmap Image | CSV 데이터는 생성됨. PNG 생성 코드 필요. | 중간 |
| **Raw Data** | DAPI / Raw Transcripts | 표준 파이프라인에는 너무 무거움(속도 저하). | 낮음 |

---

## 4. 파이프라인 단계별 분석 (Step-by-Step Analysis)

| 단계 | 모듈 | 상태 | 시각화 지원 | 개선 필요 영역 (Gap) |
| :--- | :--- | :--- | :--- | :--- |
| **Step 0** | `step0_formatting` | ✅ | N/A | 포맷 유효성 검사 (필수 컬럼 존재 여부 확인) |
| **Step 1** | `step1_preprocess` | ✅ | **우수** (QC, UMAP, PCA) | 클러스터 자동 주석 (Marker 기반 Annotation) 추가 |
| **Step 2** | `step2_seg_free` | ✅ | **양호** (거리 플롯) | 필요 시 전체 SSAM (밀도 맵) 구현 |
| **Step 4** | `step4_expansion` | ✅ | **우수** (최적화 곡선) | `config.yaml`에서 확장 범위 설정 가능하도록 개선 |
| **Step 6** | `step6_optimization` | ✅ | **부분** (CSV만 있음) | **중요**: 성능 버그(샘플링) 수정 완료, 히트맵 PNG 필요 |
| **Step 7** | `step7_domains` | ⚠️ | 없음 | **차단됨**: `SpaGCN`/`cmake` 빌드 실패로 인해 보류 |
| **Step 8** | `step8_svf` | ✅ | 공간 지도 (Spatial Map) | SVG 검출은 계산량이 많음 (옵션으로 실행) |

---

## 5. 데이터셋 처리 현황

**1. Human Brain (`human_brain`)**
*   **처리 상태**: 완료 (Complete).
*   **해결된 문제**: Violin Plot 생성 시 `leiden` 키 부재 문제 해결 (`leiden_1_4` 자동 감지).
*   **시각화**: 표준 플롯 + 마커 바이올린(Marker Violin) + 최적화 곡선 + 세그멘테이션 오버레이(Fallback 모드) 모두 생성됨.

**2. Human Breast (`h_breast_1`)**
*   **처리 상태**: 완료 (Complete).
*   **성능**: 대용량 데이터(~7GB) 처리 성공. UMAP 생성에 시간이 걸리지만 정상 완료됨.

---

## 6. 사용 방법 (How to Use)

### 전체 파이프라인 실행
```bash
# 설정된 모든 데이터셋에 대해 Step 0부터 Step 6까지 실행
bash run_all.sh
```

### 시각화만 재생성 (Regenerate Visualizations Only)
```bash
# 결과 파일을 감지하여 모든 플롯(신규 기능 포함)을 다시 그림
python run_vis_all.py
```

### 특정 플롯(Violin)만 생성
```bash
# 마커 유전자 바이올린 플롯을 위한 전용 스크립트
python run_violin_only.py
```

---

## 7. 향후 로드맵 (권장 사항)

1.  **Step 6 시각화 구현**: 생성된 ARI/Silhouette 점수 CSV를 읽어 Heatmap PNG를 그리는 기능을 추가해야 합니다.
2.  **Squidpy 통합**: `step9_spatial_stats`를 추가하여 Neighborhood Enrichment 및 Co-occurrence 분석을 수행해야 합니다. (과학적 가치가 높음)
3.  **세그멘테이션 오버레이 고도화**: 원본 폴리곤 파일(`.parquet`)이 있는 경우, 이를 찾아 "진정한" 세그멘테이션 라인을 그리는 로직을 추가해야 합니다. (현재는 중심점만 표시)
