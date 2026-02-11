# Step6 Segmentation Benchmark Review (Requested Format)

## 1. 5_segmentation_benchmark (overall objective)
- Notebook 의도: 다양한 segmentation 방법(Cellpose, Watershed, MESMER, Baysor, Clustermap, binning 등)을 같은 데이터 축에서 비교하고, 클러스터/지표를 통해 성능을 평가.
  - 근거: `notebooks/5_segmentation_benchmark/README.md:1`
- 논문: segmentation 전략 전반(Baysor, MESMER, Watershed, Cellpose, Clustermap, binning + expansion)을 benchmark하고 NMP/assigned reads/reads-per-cell/genes-per-cell 등으로 비교.
- Pipeline 대응: Step6에서 nuclei/cellpose/expansion/baysor를 합쳐 benchmark 수행.
  - 근거: `pipeline/xenium_step6_segmentation_benchmark.py:1`
- 차이: notebook/논문 대비 방법 커버리지가 축소되어 있음(특히 Watershed/MESMER/Clustermap/binning 미포함).

## 2. 6-1 Baysor execution and prior handling
- Notebook 의도: Baysor를 prior segmentation 유무/하이퍼파라미터(예: prior confidence, scale)와 함께 실행해 결과 비교.
  - 근거: `notebooks/5_segmentation_benchmark/run_baysor.py:92`
- 논문: Baysor prior segmentation confidence를 포함해 여러 설정을 비교.
- Pipeline 대응:
  - Xenium transcript를 Baysor 입력 포맷으로 변환 후 실행
  - prior segmentation TIF를 config로 주입 가능
  - dry-run 지원
  - 근거: `pipeline/xenium_step6_segmentation_benchmark.py:83`
- 차이:
  - Pipeline은 단일 Baysor 실행 경로 중심이며, notebook/논문처럼 광범위한 파라미터 sweep/다중 조건 반복 비교는 기본 제공되지 않음.

## 3. 6-2 Load segmentation result sets
- Notebook 의도: 여러 방법의 assignment 결과를 count matrix로 만들어 한 축에서 비교.
  - 근거: `notebooks/5_segmentation_benchmark/gen_counts.py:61`
- 논문: 여러 segmentation 방법을 공통 지표로 벤치마킹.
- Pipeline 대응:
  - 입력 소스: nuclei(h5ad), cellpose(h5ad), expansion(csv), baysor(csv)
  - expansion/baysor CSV는 crosstab 집계로 AnnData 변환
  - 근거: `pipeline/xenium_step6_segmentation_benchmark.py:505`
- 차이:
  - notebook 생태계의 methods(Clustermap, Watershed, MESMER, binning) 입력 경로가 step6 기본 흐름에는 없음.

## 4. 6-3 Concatenate and preprocess
- Notebook 의도: 결합된 세포들을 동일 전처리/클러스터링으로 비교(스크립트/노트북에서 louvain 기반 파라미터 사용).
  - 근거: `notebooks/5_segmentation_benchmark/5_1_Compare_Clustering on_different_segmentations.ipynb`
- 논문: 공정 비교를 위한 동일 분석 파이프라인 적용.
- Pipeline 대응:
  - QC filter(min_counts/min_genes), normalize/log1p, PCA, neighbors, UMAP, Leiden 수행
  - config 기반 파라미터화
  - 근거: `pipeline/xenium_step6_segmentation_benchmark.py:297`
- 차이:
  - notebook은 `louvain` 키/명칭으로 분석한 흔적이 강함.
  - pipeline은 Leiden 고정 구현(주석으로 차이 인지).

## 5. 6-4 Annotation transfer
- Notebook 의도: reference와의 cluster/class 대응으로 segmentation 간 세포군 비교.
  - 근거: notebook에서 crosstab 기반 cluster-to-class 대응
- 논문: scRNA 참조와 결합해 세그멘테이션 품질/NMP 해석.
- Pipeline 대응:
  - reference AnnData가 있으면 PCA-kNN majority voting으로 `celltype_majority`, `celltype_cluster` 생성
  - 근거: `pipeline/xenium_step6_segmentation_benchmark.py:234`
- 차이:
  - notebook의 수동/규칙 기반 cluster 대응과는 방식이 다름.
  - 참조 데이터 미지정 시 해당 단계 스킵.

## 6. 6-5 Benchmark metrics
- Notebook 의도: assigned reads, n_cells, reads/cell, genes/cell, NMP 등 핵심 지표를 비교.
  - 근거: `notebooks/5_segmentation_benchmark/metrics.py:9`
- 논문: 특히 assigned proportion + reads/genes per cell의 median 및 5th percentile, 그리고 NMP를 핵심으로 사용.
- Pipeline 대응:
  - 현재 저장 지표: `n_cells`, `median_reads`, `median_genes` (+ 조건부 assigned_prop)
  - 근거: `pipeline/xenium_step6_segmentation_benchmark.py:584`
- 차이:
  - 5th percentile 지표 미계산.
  - NMP 미계산(관련 함수는 `pipeline/benchmark_utils/metrics.py`에 있으나 step6에서 호출 안 함).
  - `assigned_prop`는 `subset.uns['spots']`가 있어야 계산되는데, concat/변환 경로에서 대부분 누락되어 실제로는 자주 비어 있을 가능성 큼.

## 7. 6-6 Visualizations
- Notebook 의도: UMAP, segmentation 간 세포군 분포, counts 분포 등 시각적 비교.
- 논문: UMAP/비율/분포 기반 figure로 비교 결과 제시.
- Pipeline 대응:
  - UMAP(method/leiden), spatial scatter, celltype barplot, counts violin 생성
  - 근거: `pipeline/xenium_step6_segmentation_benchmark.py:343`
- 차이:
  - 논문 Figure 수준의 다중 조건/다중 하이퍼파라미터 시각화는 축약됨.

## 8. Segmentation method coverage gap
- Notebook 의도: 방법군 자체를 폭넓게 탐색(run_segmentation, run_clustermap, run_baysor 스크립트 체계).
  - 근거: `notebooks/5_segmentation_benchmark/run_segmentation.py:47`
- 논문: Baysor/MESMER/Watershed/Cellpose/Clustermap/binning(+expansion)까지 평가.
- Pipeline 대응: nuclei/cellpose/expansion/baysor 중심.
- 차이: benchmark breadth가 논문 대비 좁아, “최적 방법 탐색”보다는 “선택된 방법 비교”에 가까움.

## 9. Step5 expansion input semantics (step6 주석과 정합)
- Notebook 의도: segmentation 결과와 assignment semantics를 맞춰 공정 비교.
- 논문: 세포 수준 assignment 품질이 NMP/assigned 비율에 직접 영향.
- Pipeline 대응:
  - step6는 step5 CSV를 expansion 입력으로 사용 가능
  - 주석에서 step5가 domain assignment 중심임을 명시
  - 근거: `pipeline/pipeline_main.py:303`
- 차이:
  - expansion 입력의 cell-level 의미가 다른 방법들과 완전히 동등하지 않을 수 있음.

---

## 최종 판정
- Step6는 “실행 가능한 segmentation benchmark 파이프라인”으로 잘 구성되어 있으며, Baysor 포함 비교 흐름을 자동화했다는 점은 강점.
- 다만 논문/원 notebook 대비로는 **부분 재현**:
  1. 방법 커버리지 축소(Watershed/MESMER/Clustermap/binning 부재)
  2. 핵심 지표 일부 미반영(NMP, 5th percentile)
  3. assigned proportion 계산 경로가 데이터 구조에 따라 실질적으로 비활성화될 수 있음

보강 우선순위:
1. step6 metrics에 `percentile_5th_reads_cells`, `percentile_5th_genes_cells`, NMP 추가
2. segmentation method 입력 확장(Clustermap/Watershed/MESMER/binning)
3. `assigned_prop` 계산을 위해 total transcript 수를 별도 메타데이터로 유지

