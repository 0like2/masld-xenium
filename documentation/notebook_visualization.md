# Notebook Visualization Analysis Report (Detailed)

사용자가 제공한 모든 코드 스니펫에 대해 파이프라인(`visualize_results.py`, `step2_*.py`) 구현 여부를 하나하나 대조한 상세 리포트입니다.

## 1. Notebook 1_5: Celltype Identification & Space (Series 1)

| # | Code Snippet (요약) | 설명 (Purpose) | 구현 상태 | 비고 / 대체 구현 |
|---|---|---|---|---|
| **1** | `sns.scatterplot(spotscorr...)` | **Replicate Correlation**: 2개 실험 간 유전자 발현 상관관계. | ❌ 미구현 | 단일 샘플 분석 구조상 기술적 재현성 분석은 제외됨. |
| **2** | `sc.pl.umap(color=['leiden_1_4'])` | **Leiden UMAP**: 클러스터링 결과 시각화. | ✅ **구현됨** | `visualize_pca_and_clustering` 함수. |
| **3** | `sc.pl.umap(color='dataset_source')` | **Batch UMAP**: 데이터 출처에 따른 분포 확인. | ❌ 미구현 | 단일 샘플에서는 불필요 (Source가 1개). |
| **4** | `sns.heatmap(gene_sorted)` | **Marker Gene Heatmap**: 클러스터별 마커 유전자 발현량. | ❌ **미구현** | **[중요]** 향후 추가 권장 항목. (현재 클러스터 번호만 출력됨) |
| **5** | `sc.pl.highly_variable_genes` | **HVG Plot**: 유전자 변동성 및 평균 발현량 분포. | ❌ 미구현 | QC용 플롯으로, 필수적이지 않아 제외됨. (`n_genes`, `counts` 바이올린 플롯은 있음) |
| **6** | `sc.pl.umap(color=['leiden_1_4'])` | **Leiden UMAP (반복)**: 위와 동일. | ✅ **구현됨** | 위와 동일. |
| **7** | `sc.pl.umap(color=['replicate'])` | **Replicate UMAP**: 실험 반복별 분포 확인. | ❌ 미구현 | 단일 샘플 분석에는 해당 없음. |
| **8** | `sc.pl.umap(color='Gfap')` | **Gene UMAP**: 특정 유전자(`Gfap`) 발현 맵. | ⚠️ **부분 구현** | Step 8에서 **자동으로 상위 SVG(공간 변수 유전자) 9개**를 찍어주지만, 사용자가 지정한 유전자를 찍는 기능은 없음. |
| **9** | `sc.pl.spatial(color='region_annotation')` | **Spatial Region Map**: 조직 영역(Annotations) 시각화. | ❌ 미구현 | **Step 7 (Spatial Domains)** 기능이 라이브러리 문제로 비활성화되어, 영역(Region) 정보가 생성되지 않음. |
| **10** | `sc.pl.spatial(color='region_level4')` | **Spatial Level 4**: 세부 조직 영역 시각화. | ❌ 미구현 | 위와 동일 (Step 7 의존). |
| **11** | `sc.pl.umap(...outline=True)` | **Region Outline UMAP**: UMAP 위에 조직 영역 윤곽선 표시. | ❌ 미구현 | 위와 동일 (Step 7 의존). |
| **12** | `sc.pl.umap(...adata_oligo...)` | **Subset UMAP**: 특정 세포군(Oligo/Astro)만 하이라이트. | ❌ 미구현 | 특정 가설 검증용 커스텀 플롯이라 일반 파이프라인에는 포함되지 않음. |

---

## 2. Notebook 2_3: Cell Overlaps (Step 2)

| # | Code Snippet (요약) | 설명 (Purpose) | 구현 상태 | 비고 / 대체 구현 |
|---|---|---|---|---|
| **13** | `plt.scatter(df.x,df.y)` | **Raw Transcript Scatter**: 전체 전사체 좌표 점찍기. | ⚠️ **대체됨** | `step2_ssam`에서 **Density Map (밀도 맵)**으로 더 보기 좋게 대체 구현됨. 점(Scatter)은 데이터가 너무 커서 렌더링이 느림. |
| **14** | `correlations_**15` (Softmax) | **Signal Coherence**: 신호 일관성 계산 로직. | ❌ **생략됨** | 매우 무거운 매트릭스 연산으로, 벤치마킹 파이프라인의 속도를 위해 제외됨. |
| **15** | `plt.imshow(distance_...winter_r)` | **Cosine Sim Heatmap**: Z축 상/하단 유사도 맵. | ❌ **단순화됨** | `step2_overlaps`에서 **"Z-Range Histogram"**으로 단순화하여 Z축 데이터 퀄리티를 체크함. |
| **16** | `vis.plot_instance` (Rectangles) | **ROI Visualization**: 특정 영역(ROI) 확대 시각화. | ❌ **미구현** | 자동화된 파이프라인에서는 특정 ROI를 지정하기 어려워 전체 통계(Histogram)로 대체함. |
| **17** | `ax1.bar(...)` | **Incoherence Histogram**: 비일관성 수치 분포. | ✅ **유사 구현** | `step2_overlaps`의 `z_range_distribution.png` (Z축 범위 분포)가 유사한 목적을 수행함. |

---

## 3. Notebook 2_4: SSAM (Step 2-4)

| # | Code Snippet (요약) | 설명 (Purpose) | 구현 상태 | 비고 / 대체 구현 |
|---|---|---|---|---|
| **18** | `plt.scatter(...nipy_spectral)` | **Class Scatter**: 세포 타입별 분포. | ❌ **미구현** | SSAM의 전체 워크플로우(De novo clustering)가 수행되지 않아 클래스 정보가 없음. |
| **19** | `ds.plot_localmax()` | **Local Maxima**: 밀도 최고점(세포 중심) 시각화. | ❌ **생략됨** | `step2_ssam`은 계산 부하를 줄이기 위해 **Local Max 탐색 과정을 생략**하고 밀도 맵(Density Map)만 생성함. |
| **20** | `ds.plot_celltypes_map(...)` | **SSAM Celltype Map**: 벡터 기반 세포 지도. | ❌ **미구현** | `ssam` 패키지의 무거운 연산(벡터 필드, 로컬 맥스)이 필요하여 제외됨. |
| **21** | `sc.pl.highest_expr_genes` | **Top Genes**: SSAM 결과의 상위 발현 유전자. | ❌ **미구현** | SSAM AnnData를 생성하지 않으므로 불가능. |
| **22** | `sc.pl.pca`, `sc.pl.umap` | **SSAM PCA/UMAP**: SSAM 결과의 차원 축소. | ❌ **미구현** | 위와 동일. |
| **23** | `ds.plot_celltypes_map(rotate=3)` | **De Novo Map**: SSAM 기반 De novo 클러스터링 지도. | ❌ **미구현** | 위와 동일. |

---

## 4. Notebook 3-2: Custom Segmentation & Assignment (New)

이 부분은 **Raw Image Processing (DAPI)** 및 **Custom Segmentation (Cellpose)** 영역입니다.

| # | Code Snippet (요약) | 설명 (Purpose) | 구현 상태 | 비고 / 대체 구현 |
|---|---|---|---|---|
| **24** | `plt.imshow(full_mask)` | **Mask Visualization**: 세포 영역(Mask) 시각화. | ❌ **제외됨 (Out of Scope)** | 현재 파이프라인은 Xenium 장비가 제공한 기본 Segmentation(`cell_boundaries.csv`)을 신뢰하고 사용합니다. 직접 DAPI 이미지를 다시 세그멘테이션하는 과정은 벤치마킹 범위를 벗어납니다. |
| **25** | `plt.scatter(adata.obs.x...)` | **Centroid Plot**: 추출된 세포 중심점 플롯. | ✅ **유사 구현** | `visualize_step8_svf` 등에서 세포 플롯(`sc.pl.spatial`)을 통해 중심점을 기반으로 한 시각화를 수행합니다. |
| **26** | `plt.imshow(dapi_image)` | **DAPI Image Overlay**: DAPI 핵 염색 이미지 확인. | ❌ **미구현** | 대용량 TIFF/OME-TIFF 이미지 처리 모듈이 없습니다. |
| **27** | `model = models.Cellpose(...)` | **Nuclei Segmentation**: Cellpose 딥러닝 모델로 핵 분할. | ❌ **미구현** | **커스텀 세그멘테이션**은 파이프라인에 포함되어 있지 않습니다. (장비 제공 결과 사용) |
| **28** | `closest_cell.append(...)` | **Read Assignment**: 전사체를 세그멘테이션 마스크에 할당. | ✅ **기본 기능** | Xenium 포맷팅 단계(Step 0)에서 이미 할당된 데이터(`cell_feature_matrix`)를 로드합니다. 파이프라인 내에서 재할당 로직은 수행하지 않습니다. |
| **29** | `plt.scatter(read_positions['x']...)` | **Reads vs Nuclei**: 전사체와 핵 위치 함께 시각화. | ⚠️ **대체됨** | `step2_ssam/transcript_density_map.png`에서 전사체의 밀도 분포를 확인하는 것으로 대체합니다. |

---

## 5. 종합 결론 (Overall Conclusion)

1.  **시각화 커버리지**:
    *   **분석용 시각화 (UMAP, QC)**: 대부분 구현되어 있습니다.
    *   **검증용 시각화 (Heatmap)**: **구현 필요 (Marker Gene Heatmap)**.
    *   **탐색용 시각화 (Raw Reads, DAPI)**: 대용량 데이터 렌더링 문제 및 파이프라인 목적(Downstream Analysis) 상 **제외**되었습니다.

2.  **Segmentation (Notebook 3-2)**:
    *   보내주신 Notebook 3-2는 **"Raw Data부터 다시 Segmentation을 수행"**하는 과정입니다.
    *   현재 파이프라인은 **"이미 Segmentation된 데이터(AnnData)를 받아 분석"**하는 것에 초점이 맞춰져 있으므로, DAPI 이미지 처리나 Cellpose 모델 구동 코드는 포함되어 있지 않습니다.
