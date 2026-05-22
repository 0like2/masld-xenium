# Step4 Techniques Comparison Detailed Review (Requested Format)

## 1. 4-1 Load resegmented + original data
- Notebook 의도: 리세그멘테이션 결과와 원본 데이터를 함께 불러와 이후 비교 분석의 공통 입력으로 사용.
- 논문: 공정 비교를 위해 공통 기준(세그멘테이션/영역 정렬/참조 데이터)을 맞춘 뒤 비교 수행.
- Pipeline 대응: `run_step4()`에서 resegmented/original AnnData와 transcript CSV를 로드.
  - 근거: `pipeline/xenium_step4_techniques_comparison.py:110`
- 차이: 큰 방향은 일치하나, 원본 transcript 포맷에 따라 일부 하위 분석(특히 diffusion)이 자동 스킵될 수 있음.

## 2. 4-2 Efficiency analysis
- Notebook 의도: 기술 간/영역 간 유전자 검출 효율 비교(특히 scRNA 대비 ratio), region별 분해 비교.
  - 근거: `notebooks/3_techniques_comparison/3_3_efficiency_between_methods.ipynb`
- 논문: positive-cell 기반 median 중심 효율 비교(SRT/scRNA), 공통 영역 기준 비교.
- Pipeline 대응:
  - QC 분포 비교(histogram, 요약통계): `pipeline/xenium_step4_techniques_comparison.py:225`
  - ST/scRNA ratio 계산(옵션): `pipeline/xenium_step4_techniques_comparison.py:256`
  - region breakdown(옵션): `pipeline/xenium_step4_techniques_comparison.py:303`
- 차이:
  - Pipeline은 mean CPM ratio(`target_sum=1e6`) 중심이라 논문/노트북의 median 기반 정의와 직접 동일하지 않음.
  - region 키 탐색에 `region_annotation`이 없어 Step3 도메인 결과를 놓칠 수 있음.

## 3. 4-3 Specificity analysis
- Notebook 의도: NMP/NCP 기반 specificity 비교 + efficiency와의 관계 분석.
  - 근거: `notebooks/3_techniques_comparison/3_4_negative_marker_purity_for_specificity.ipynb`
- 논문: NCP/NMP를 핵심 specificity 지표로 사용하여 플랫폼/전략 비교.
- Pipeline 대응:
  - NMP-like 계산(참조 scRNA 필요): `pipeline/xenium_step4_techniques_comparison.py:349`
  - 보조적으로 gene-gene correlation heatmap: `pipeline/xenium_step4_techniques_comparison.py:369`
- 차이:
  - Pipeline은 논문 full protocol 대비 단순화된 proxy 성격.
  - top-gene 선택/교차 subset 시 gene set 불일치에 취약할 수 있음.

## 4. 4-4 Positivity analysis
- Notebook 의도: 전처리 후 positivity, 클러스터별 marker/UMAP 시각화, 다중 해상도 비교.
  - 근거: `notebooks/3_techniques_comparison/3_5_Computing_positivity_after_preprocessing_for_all_ST_techs.ipynb`
- 논문: 전처리/표현형 분해에 따른 세포군 해석력 평가 맥락.
- Pipeline 대응:
  - positivity histogram + 상위 유전자 테이블: `pipeline/xenium_step4_techniques_comparison.py:423`
  - Leiden/UMAP/violin + gene별 optimal cluster: `pipeline/xenium_step4_techniques_comparison.py:432`
- 차이:
  - Notebook의 파라미터/다중 resolution 탐색 대비 Pipeline은 단일 단순 경로(Leiden 1.0).

## 5. 4-5 Diffusion analysis
- Notebook 의도: 기술별 transcript-to-centroid 거리 분포(ECDF), gene별 확산 차이 비교.
  - 근거: `notebooks/3_techniques_comparison/3_6_Diffussion_on_resegmented_data.ipynb`
- 논문: 플랫폼별 subcellular diffusion 특성 비교.
- Pipeline 대응:
  - px->um 변환 적용, complementary ECDF, per-gene ECDF, gene×method heatmap.
  - 근거: `pipeline/xenium_step4_techniques_comparison.py:498`
- 차이:
  - Pipeline은 주로 resegmented vs original 내부 비교 축.
  - notebook의 다중 플랫폼 통합 비교와는 범위가 다름.

## 6. Cross-cutting issue: Region key compatibility
- Notebook 의도: region 기반 비교를 강하게 전제(CTX/HPF/TH 등).
- 논문: 공통 해부학 영역 정렬이 핵심 전제.
- Pipeline 대응: region 후보를 `spatial_annotation/region/tissue_region/domain`만 검색.
  - 근거: `pipeline/xenium_step4_techniques_comparison.py:305`
- 차이: Step3이 쓰는 `region_annotation`(`pipeline/xenium_step3_resegmentation.py:191`)이 step4에서 자동 인식되지 않음.

## 7. 3_7_Xenium_vs_Visium_comparison.ipynb
- Notebook 의도: area-normalized gene counts, X/V ratio, region별 side-by-side.
  - 근거: `notebooks/3_techniques_comparison/3_7_Xenium_vs_Visium_comparison.ipynb`
- 논문: Xenium가 Visium 대비 tissue-level sensitivity 우수(핵심 결과).
- Pipeline 대응: 없음.
- 차이: 핵심 missing block.

---

## 최종 판정
- Step4는 `3_3~3_6`의 핵심 프레임을 부분적으로 재현.
- 논문 충실 재현 관점에서는 아직 미완성:
  1. `3_7 Xenium vs Visium` 미구현
  2. efficiency 정의 불일치(mean CPM vs median-based)
  3. region 키 정합성(`region_annotation`) 보강 필요

