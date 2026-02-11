# Step3 Techniques Comparison Review (Requested Format)

## 1. 3_1_resegmentation_notebooks
- Notebook 의도: Cellpose 기반 핵 분할 후, 확장 마스크(`expand_labels`)를 사용해 transcript를 `in_cell`/`closest_cell`로 할당.
  - 근거: `notebooks/3_techniques_comparison/3_1_resegmentation_notebooks/batch_segmentation_cellpose-Xenium.ipynb`
- 논문: 공통 segmentation 전략(Cellpose 등)으로 기술 간 비교 기반을 정리하고, diffusion 성격 파악을 위해 centroid 거리 관점 분석 수행.
- Pipeline 대응:
  - Cellpose nuclei segmentation 수행
  - transcript -> cell 라벨 매핑
  - cell×gene matrix 생성 및 AnnData 저장
  - 근거: `pipeline/xenium_step3_resegmentation.py:75`
- 차이:
  - Notebook의 `closest_cell`(확장 라벨 기반) 흐름과 달리 Pipeline은 mask 라벨 기반 할당 중심.
  - 확장 기반 nearest assignment를 동일 형태로 export하지 않음.

## 2. 3_2_cell_to_domain_assignment_notebooks
- Notebook 의도: 영역 polygon(json 등)을 이용해 세포를 domain/region에 매핑.
  - 근거: `notebooks/3_techniques_comparison/3_2_cell_to_domain_assignment_notebooks/*`
- 논문: 공통 해부학 영역 정렬을 통해 기술 간 공정 비교 수행.
- Pipeline 대응:
  - Step3에서 optional domain assignment 구현(`domain_map_path` 기반 polygon 포함 판정)
  - 결과를 `region_annotation`에 저장
  - 근거: `pipeline/xenium_step3_resegmentation.py:177`
- 차이:
  - Notebook은 데이터셋별 domain assignment 워크플로우가 더 풍부함.
  - Pipeline은 단일 optional 경로 중심.

## 3. Segmentation model and runtime behavior
- Notebook 의도: Cellpose를 다양한 노트북 변형으로 실행하며 결과 비교/실험.
- 논문: Cellpose/Baysor 등 여러 전략을 benchmark해 최적 조합 탐색.
- Pipeline 대응:
  - 디바이스 자동 감지(CUDA>MPS>CPU)
  - Cellpose `nuclei` 모델 실행, 대형 이미지 타일 옵션
  - 근거: `pipeline/xenium_step3_resegmentation.py:35`
- 차이:
  - Step3 자체는 Cellpose 중심 재분할 구현이고, Baysor/다전략 벤치마크는 Step6 책임.

## 4. Transcript assignment output semantics
- Notebook 의도: 후속 효율/특이성/diffusion 비교를 위한 assignment 테이블(예: `closest_cell`) 축적.
- 논문: 비교 목적에 맞게 assignment 규칙을 명확히 구분(효율/특이성 vs diffusion).
- Pipeline 대응:
  - `cell_id_reseg > 0` transcript만 `step3_transcripts_resegmented.csv`로 저장
  - 근거: `pipeline/xenium_step3_resegmentation.py:157`
- 차이:
  - diffusion 분석 맥락에서 notebook/논문의 외부 read 처리 전략과 결과가 달라질 수 있음.

## 5. Data products (what Step3 actually produces)
- Notebook 의도: segmentation/assignment 중간산출물을 후속 비교 분석에 활용.
- 논문: segmentation 품질과 downstream 영향(효율, 특이성, diffusion) 평가.
- Pipeline 대응:
  - resegmented mask TIF
  - resegmented transcript CSV
  - resegmented AnnData(H5AD)
  - optional domain-annotated AnnData
  - 근거: `pipeline/xenium_step3_resegmentation.py:97`
- 차이:
  - 산출물 구성은 실무적으로 충분하지만, notebook의 실험형 중간 결과(다수 변형)는 축약됨.

## 6. Step3 -> Step4 연결 관점 이슈
- Notebook 의도: domain/region/assignment 정보를 step-by-step 비교 지표로 직접 연결.
- 논문: region-aligned 비교가 핵심.
- Pipeline 대응:
  - Step3가 `region_annotation`을 생성 가능
  - Step4의 region 탐색 후보에 `region_annotation`이 없음
  - 근거: `pipeline/xenium_step3_resegmentation.py:191`
- 차이:
  - Step3에서 생성한 region 정보가 Step4 지역 분석에 자동 반영되지 않을 수 있음.

## 7. Overall parity judgment for Step3
- Notebook 의도: 재분할 + 도메인 할당 + 후속 비교를 위한 데이터 구조 준비.
- 논문: segmentation 전략 차이가 downstream 결과에 미치는 영향 검증.
- Pipeline 대응: Step3는 “Cellpose 기반 재분할 파이프라인”으로 기능적으로 잘 구현됨.
- 차이:
  - 논문/노트북의 전체 segmentation benchmark 관점(다전략/세부 파라미터 비교)은 Step3 단독으로는 미포함.
  - closest-cell 기반 diffusion 친화적 assignment semantics는 동일 재현 아님.

---

## 최종 판정
- Step3는 재분할 파이프라인으로서는 구현 완성도가 높음.
- 다만, 논문/노트북과의 엄밀 동등성 관점에서는 다음 보강이 필요:
  1. closest-cell assignment 옵션 제공(특히 diffusion 정합성 목적)
  2. Step4와 region 키(`region_annotation`) 연동 강화
  3. 논문식 다전략 benchmark 결과와의 연결은 Step6 포함 통합 검증 필요

