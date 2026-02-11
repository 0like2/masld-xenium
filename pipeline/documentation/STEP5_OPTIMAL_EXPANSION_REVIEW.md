
# Step5 Optimal Expansion Review (Requested Format)

## 1. 4_1_Optimal_expansion_multisection.ipynb (overall objective)
- Notebook 의도: domain annotation이 있는 read/cell을 기반으로, cell type × domain별 turnover 지점을 구해 최적 expansion을 추정.
  - 근거: `notebooks/4_optimal_expansion/4_1_Optimal_expansion_multisection.ipynb`
- 논문: nucleus signature와 domain background signature의 correlation crossover를 이용해 optimal expansion 정의.
- Pipeline 대응: Step5 전체가 해당 notebook의 KDTree assignment + turnover 계산 로직을 반영.
  - 근거: `pipeline/xenium_step5_optimal_expansion.py:1`
- 차이: notebook의 실험형/수동 분석 흐름을 production 파이프라인 형태로 단순화.

## 2. 5-1 Load original reads from Step0
- Notebook 의도: 원본 reads를 기준으로 domain 재할당 및 turnover 분석 입력 구성.
- 논문: Xenium segmentation/expansion 결과를 바탕으로 nuclear/background signature 비교.
- Pipeline 대응:
  - Step0 `.h5ad`의 `uns['spots']`에서 reads 로딩
  - Step0 파일 탐색 fallback 포함
  - 근거: `pipeline/xenium_step5_optimal_expansion.py:331`
- 차이:
  - Step0 산출물(`spots`) 의존성이 강함.
  - `spots` 누락 시 Step5 전체가 중단됨.

## 3. 5-2 Load annotated cells (domain source)
- Notebook 의도: cell metadata(`Class`, `spatial_annotation`)를 읽어 reads에 initial annotation/domain 부여.
- 논문: cell type/domain 단위로 최적 expansion 추정.
- Pipeline 대응:
  - `previous_step_adata_path` 우선 사용
  - 실패 시 내부 경로 fallback 탐색
  - 근거: `pipeline/xenium_step5_optimal_expansion.py:354`
- 차이:
  - fallback 경로에 `step4_resegmentation`이 하드코딩되어 있는데, 현재 파이프라인 실제 폴더명은 `step3_resegmentation`.
  - 즉 fallback 경로 일부가 불일치.

## 4. 5-3 Map domain assignments to reads
- Notebook 의도: `cell_id -> domain`, `cell_id -> initial_annotation` 매핑 후 unassigned reads 분리.
  - 근거: notebook code에서 `domain`, `initial_annotation`, `nancells/annotatedcells` 분리.
- 논문: domain-specific background와 cell-type signature 비교를 위한 전처리.
- Pipeline 대응:
  - domain key 후보 탐색 후 mapping
  - `reads_original['domain']`, `reads_original['initial_annotation']` 생성
  - 근거: `pipeline/xenium_step5_optimal_expansion.py:381`
- 차이:
  - domain key 우선순위가 `['spatial_annotation','Class','leiden','cluster','graph_clusters']`로 제한됨.
  - `region_annotation` 등 다른 유효 키는 자동 사용되지 않을 수 있음.

## 5. 5-4 KDTree nearest-domain expansion
- Notebook 의도: domain이 있는 reads 일부를 anchor로 샘플링 후 KDTree로 unassigned reads의 nearest domain 추정.
  - 근거: notebook에서 `cKDTree`, 1% subsampling.
- 논문: 최적 expansion 계산 전 background/domain signature 구성의 기반 단계.
- Pipeline 대응:
  - `subsample_fraction` 기반 anchor 샘플링
  - `cKDTree.query(k=1)`로 nearest domain 할당
  - 근거: `pipeline/xenium_step5_optimal_expansion.py:428`
- 차이:
  - notebook/논문의 분석 맥락은 동일하나, 샘플링 비율/threshold에 따라 결과 변동성이 큼.

## 6. 5-5 Optional distance threshold
- Notebook 의도: 기본 nearest 할당 후 결과 검토.
- 논문: method 설명에서는 turnover 기반 기준이 핵심이며 hard threshold는 본질 아님.
- Pipeline 대응:
  - `distance_threshold` 설정 시 너무 먼 read를 NaN 처리 가능
  - 근거: `pipeline/xenium_step5_optimal_expansion.py:466`
- 차이:
  - 논문 핵심 알고리즘에는 없는 운영용 옵션(실무용 guardrail).

## 7. 5-6 / 5-7 Visualization and output export
- Notebook 의도: 도메인 재할당 결과 시각화 및 파일 저장.
- 논문: 결과 해석을 위한 시각화(Fig. 3 계열) 제공.
- Pipeline 대응:
  - spatial scatter + distance histogram 저장
  - expanded transcripts CSV 저장
  - 근거: `pipeline/xenium_step5_optimal_expansion.py:492`
- 차이:
  - 시각화 스타일은 논문 figure를 완전히 복제하지 않고 운영형 요약 플롯으로 제공.

## 8. 5-8 Turnover/Crossover analysis
- Notebook 의도:
  - distance-bin × gene matrix 구성
  - nuclear vs background correlation 곡선 계산
  - `diff < 0.1` 최초 거리 = turnover
  - celltype별 turnover 집계
- 논문:
  - cell type–domain pair(>=5,000 reads)에서 distance별 signature를 nuclear/background와 Pearson 상관 비교
  - background correlation이 nuclear correlation을 넘어서는 최소 거리로 optimal point 정의
- Pipeline 대응:
  - `calculate_turnover()`에서 동일한 핵심 로직 구현
  - `min_reads_per_domain`, `diff_threshold` 파라미터화
  - 근거: `pipeline/xenium_step5_optimal_expansion.py:73`
- 차이:
  - `overlaps_nucleus`가 없으면 median-distance 기반 proxy로 대체.
  - 논문의 1-um interval 명시는 코드상 `round(0)` 기반으로 근사 적용.

## 9. Optimal expansion formula and outputs
- Notebook 의도: 평균 turnover와 nuclei size를 비교해 expansion from nucleus border 산출.
- 논문: optimal expansion = turnover distance - nuclei edge distance.
- Pipeline 대응:
  - `optimal_expansion = mean(turnover) - mean(nuclei_size)` 계산
  - summary matrix/per-celltype/txt 결과 저장
  - 근거: `pipeline/xenium_step5_optimal_expansion.py:255`
- 차이:
  - 핵심 수식은 일치.
  - 입력 데이터 품질(centroid, overlaps_nucleus, domain key)에 따라 안정성이 달라짐.

## 10. Step5 -> Step6 handoff semantics
- Notebook 의도: expansion 결과를 downstream 비교/검증에 활용.
- 논문: segmentation/expansion이 downstream 품질 지표에 미치는 영향 평가.
- Pipeline 대응:
  - Step5 CSV를 Step6 입력으로 연결
  - 근거: `pipeline/pipeline_main.py:303`
- 차이:
  - 파이프라인 주석에도 명시된 대로 Step5는 현재 domain 할당 중심이며, cell-level reassignment semantics와 완전 동일하진 않음.

---

## 최종 판정
- Step5는 `4_1` notebook과 논문의 optimal expansion 핵심 아이디어를 **상당히 충실하게 재구현**했습니다.
- 특히 turnover/crossover 계산의 구조는 논문 Methods와 정합성이 높습니다.

다만 운영/정합성 관점의 보강 포인트:
1. fallback 경로명 불일치(`step4_resegmentation` -> `step3_resegmentation`) 수정 필요
2. domain key 탐색 후보에 `region_annotation` 등 추가 필요
3. `overlaps_nucleus` 부재 시 proxy 사용 여부를 결과 메타데이터에 명시 권장

