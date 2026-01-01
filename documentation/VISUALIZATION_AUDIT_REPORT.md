# 노트북 시각화 검토 보고서 (Visualization Audit Report)

사용자가 요청한 특정 노트북 목록에 대해 시각화 로직의 파이프라인 포함 여부를 검토한 결과입니다.

## 1. 개요
대부분의 핵심 시각화(QC, UMAP, 거리 분포, 최적화 곡선)는 구현되어 있습니다. 그러나 **비교 분석(Comparison)** 및 **심화 구조 분석(Structure Scores)** 관련 시각화는 현재 파이프라인의 범위(단일 샘플 처리 중심)를 벗어나거나 아직 구현되지 않았습니다.

## 2. 노트북별 상세 검토

### 1. 데이터셋 탐색 (Dataset Exploration)
| 노트북 | 파이프라인 상태 | 주요 시각화 & 누락 사항 |
| :--- | :---: | :--- |
| **`1_1_Statistics_all_samples...`** | ✅ **구현됨** | - **QC Plots**: Read counts, Genes per cell 바이올린 플롯 (Step 1).<br>- **구현**: `visualize_step1_qc`. |
| **`1_2_Celltype_identification...`** | ✅ **구현됨** | - **Cluster Maps**: UMAP 및 공간 상의 클러스터 분포 (Step 1).<br>- **구현**: `visualize_pca_and_clustering`. |
| **`1_3_Read_specific_dispersion...`** | ⚠️ **부분 구현** | - **Dispersion**: 유전자별 확산도 바/박스 플롯 (Step 2).<br>- **차이**: 리드 단위의 세밀한 "Dispersion Metrics" 비교 시각화는 단순화됨 (거리 분포로 대체). |
| **`1_5_...EXPANDED-cells...`** | ⚪ **제외됨** | - **내용**: "Expanded" 세그멘테이션과 원본 비교.<br>- **현황**: 파이프라인은 단일 세그멘테이션 입력을 가정하므로, **비교 시각화**는 없음. |
| **`1_6_Batch_preprocessing...`** | ⚪ **제외됨** | - **내용**: 여러 배치의 통합 및 보정 전후 비교.<br>- **현황**: 파이프라인은 **단일 샘플(Single Sample)** 처리에 최적화됨. 배치 효과 시각화는 다중 샘플 통합 분석 단계에서 필요. |
| **`1_7_...structure_scores...`** | ❌ **누락됨** | - **내용**: **세포 구조 점수(Structure Scores)** 및 조직학적 구조 시각화.<br>- **분석**: 현재 파이프라인 Step 1~8 어디에도 해당 로직(Neighborhood enrichment 등) 없음.<br>- **추천**: 추후 "조직 구조 분석" 단계로 추가 고려 가능. |

### 2. Segmentation-Free 분석
| 노트북 | 파이프라인 상태 | 주요 시각화 & 누락 사항 |
| :--- | :---: | :--- |
| **`2_1_...distance_to_nuclei...`** | ✅ **구현됨** | - **Distance Plots**: 유전자-핵 거리 박스플롯 (Step 2).<br>- **구현**: `visualize_step2_gene_distances`. |
| **`2_3_brain_cell_overlaps`** | ✅ **신규 구현** | - **Overlap Hist**: 세포 중첩도 및 Z축 분포 히스토그램.<br>- **현황**: `step2_overlaps.py` 모듈 추가로 해결됨. |
| **`2_4_brain_ssam`** | ✅ **신규 구현** | - **Density Map**: SSAM (KDE) 밀도 맵.<br>- **현황**: `step2_ssam.py` 모듈 추가로 해결됨. |
| `DEPRECATED_2_2...` | ⚪ 제외됨 | Deprecated 노트북. |

### 3. 기술 비교 (Techniques Comparison)
| 노트북 | 파이프라인 상태 | 주요 시각화 & 누락 사항 |
| :--- | :---: | :--- |
| **`3_X` 전체 시리즈** | ⚪ **범위 외** | - **내용**: Xenium vs Visium/CosMx/Merscope 기술 간 성능/일치도 비교.<br>- **현황**: 벤치마킹 파이프라인이 "단일 Xenium 샘플 분석"에 초점을 맞추고 있어, **타 기술 데이터와의 1:1 비교 시각화**는 포함되어 있지 않음. |
| `3_7_Xenium_vs_Visium...` | ⚪ **범위 외** | 상동. |

### 4. 최적 확장 (Optimal Expansion)
| 노트북 | 파이프라인 상태 | 주요 시각화 & 누락 사항 |
| :--- | :---: | :--- |
| **`4_1_Optimal_expansion...`** | ✅ **구현됨** | - **Optimization Curve**: 확장 거리에 따른 순도/포획률 변화 곡선.<br>- **구현**: `visualize_step4_optimal_expansion`. |

### 5. 세그멘테이션 벤치마크
| 노트북 | 파이프라인 상태 | 주요 시각화 & 누락 사항 |
| :--- | :---: | :--- |
| **`5_1_Compare_Clustering...`** | ⚪ **제외됨** | - **내용**: 서로 다른 세그멘테이션 방법(Cellpose, Baysor 등) 간의 클러스터링 일치도 비교.<br>- **현황**: 단일 세그멘테이션 입력 기반이므로 비교 시각화 없음. |

## 3. 결론 및 제언
1.  **핵심 기능 일치**: Step 1, 2, 4, 6의 주요 분석 시각화는 모두 구현되었습니다.
2.  **구조 분석 누락**: **`1_7` (Structure Scores)** 관련 분석 및 시각화는 현재 없습니다. 조직의 공간적 패턴(Neighborhood Analysis)이 중요하다면 추가 개발이 필요합니다.
3.  **비교 분석 제외**: 시리즈 3, 5와 같은 "방법론 간 비교"는 현재 파이프라인의 목적(단일 샘플 파이프라인 실행)과 다르므로 제외되었습니다.

이 보고서를 바탕으로 추가적으로 구현하고 싶은 시각화(예: `1_7`의 구조 점수 등)가 있다면 말씀해 주세요.
