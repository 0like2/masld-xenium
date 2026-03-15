# Pipeline Logic Structure Review: 논문 vs 노트북 vs 파이프라인

> **작성일**: 2026-02-22
> **목적**: 논문, 노트북, 파이프라인 간 논리 구조 비교 및 불일치 진단

---

## 1. 논문의 큰 이야기 (왜 이 연구를 했는가?)

> "Xenium이라는 새 현미경이 얼마나 좋은지 테스트하고,
> 데이터 분석하는 최고의 방법을 찾아주자!"

논문은 크게 **두 가지**를 합니다:
1. **Xenium 평가**: 다른 기술들과 비교해서 얼마나 좋은가?
2. **분석 도구 평가**: 어떤 분석 방법이 가장 좋은가?

---

## 2. 논문의 분석 순서 (A~K)

```
┌─────────────────────────────────────────────────────────────────┐
│ A. 데이터 준비: 25개 Xenium 데이터셋 정리                         │
│    "실험 결과물을 컴퓨터가 읽을 수 있게 정리하자"                    │
└─────────┬───────────────────────────────────────────────────────┘
          │
    ┌─────┴─────┬──────────────┬──────────────┐
    ▼           ▼              ▼              ▼
┌────────┐ ┌──────────┐ ┌──────────┐ ┌──────────────────┐
│ B.세포  │ │ C.세포없이│ │ D.6개    │ │ G.전처리 방법    │
│ 분류   │ │ 분석     │ │ 기술비교 │ │ 비교 (시뮬레이션)│
│        │ │(SSAM등)  │ │(cross-   │ │                  │
│"이 세포│ │"세포경계 │ │platform) │ │"데이터 정리 순서 │
│는 뭐?" │ │없이도   │ │"Xenium이 │ │  어떻게 하면    │
│        │ │분석 가능"│ │ 제일?" │ │  가장 정확?"    │
└───┬────┘ └──────────┘ └──────────┘ └──────────────────┘
    │ (세포 분류 결과 필요)
    │
    ├──────────────┬──────────────┐
    ▼              ▼              ▼
┌────────┐  ┌──────────┐  ┌──────────────┐
│ E.확장  │  │ H.공간   │  │ I.유전자      │
│ 거리   │  │ 유전자   │  │ 예측         │
│ 최적화 │  │ 찾기     │  │ (imputation) │
│        │  │ (SVF)    │  │              │
│"핵에서 │  │"어디서   │  │"빠진 유전자를 │
│몇 µm   │  │특별히   │  │ 추측하자"    │
│확장?"  │  │발현?"   │  │              │
└───┬────┘  └──────────┘  └──────────────┘
    │
    ▼
┌──────────────┐     ┌──────────────┐
│ F.세그멘테이션│     │ J.영역 구분   │
│ 벤치마크     │     │ (도메인)      │
│              │     │              │
│"Baysor vs    │     │"뇌의 어디가  │
│ Cellpose vs  │     │ 어느 영역?"  │
│ 기본 핵?"    │     │              │
└──────────────┘     └──────────────┘
          │                  │
          └───────┬──────────┘
                  ▼
         ┌──────────────┐
         │ K. 최종 권장  │
         │ 파이프라인    │
         │              │
         │"이렇게 분석  │
         │ 하세요!"     │
         └──────────────┘
```

**핵심 포인트**: D(cross-platform 비교)와 F(segmentation 비교)는 **별개의 분석**입니다.

---

## 3. 논문 섹션별 상세 설명

### A. 데이터 준비 (Fig. 1a-b)
- **입력**: 25개 Xenium 실험 원본 데이터 (mouse brain, human cancer 등)
- **분석**: QC 통계 (reads/cell, genes/cell, QV>20 비율)
- **결론**: 81% reads가 QV>20, 평균 186.6 reads/cell, 0.21%만 저품질

### B. 세포 분류 (Fig. 1c-d)
- **입력**: Mouse brain 7개 연속 절편
- **분석**: Xenium 기본 segmentation → normalize → PCA → Leiden → 50개 세포 유형 식별
- **결론**: 실험 간 세포 유형 비율이 일관됨, 재현성 확인

### C. 세포 없이 분석 (Fig. 1e-k)
- **입력**: Xenium 3D read 좌표 (x, y, z)
- **분석**: SSAM (segmentation-free), ovrlpy (z축 coherence), Points2Regions
- **결론**: 세포 경계 없이도 분석 가능, 3D 해상도가 핵심 장점

### D. Cross-Platform 비교 (Fig. 2) ★
- **입력**: 6개 SRT 플랫폼 데이터 (CosMx, HS-ISS, MERFISH, MERSCOPE, Xenium, MC) + scRNA-seq
- **분석**:
  - 모든 플랫폼을 **Cellpose로 동일하게 재분할** (공정한 비교를 위해)
  - **Efficiency**: transcripts/cell, genes/cell → Fig. 2b (box plot, 6개 기술)
  - **Specificity**: NCP score → Fig. 2d
  - **Diffusion**: distance to centroid → Fig. 2f
  - Xenium vs Visium → Fig. 2g
- **결론**: Xenium이 ISH 기반 기술과 동등, scRNA-seq보다 1.2-1.5배 높은 효율

### E. 확장 거리 최적화 (Fig. 3a-b)
- **입력**: Mouse brain dataset 1 + 세포 유형 주석
- **분석**: 핵에서 거리별 correlation (nuclear vs background signature)
- **결론**: 최적 확장 ~5.64 µm, 세포 유형마다 다름 → 고정 15µm 확장은 과도

### F. Segmentation 벤치마크 (Fig. 3c-h) ★
- **입력**: Mouse brain Xenium + DAPI + 다양한 segmentation 알고리즘
- **분석**:
  - 52개 segmentation 조합 테스트
  - **핵심 메트릭 2개**: Proportion of assigned reads (X축) vs NMP (Y축) → Fig. 3e
  - ARI heatmap → Fig. 3d
  - Counts/cell violin → Fig. 3g
  - Cell type frequency bar → Fig. 3h
- **결론**: **Baysor + Xenium nuclei prior (BA2 P0.8)이 최고**

### G. 전처리 방법 비교 (Fig. 4a-e)
- **입력**: scRNA-seq → Xenium 시뮬레이션 데이터 생성
- **분석**: 전처리 조합 grid search (normalization, log, HVG, PCs, neighbors, clustering)
- **결론**: Library-size norm(100) + log + scale + all PCs + 16 neighbors + Louvain

### H-J. SVF, Imputation, Domain (Fig. 4f-j, 5)
- 각각 독립적인 벤치마크 트랙
- SVF: Hotspot 추천 (false positive 최소)
- Imputation: SpaGE 추천
- Domain: binning-based가 의외로 최고

---

## 4. 노트북 구조 (논문 저자의 실제 코드)

```
노트북 0: 데이터 포맷팅 ─────────────────────────────── 논문 A
    │
노트북 1: 데이터 탐색 + 세포 분류 ──────────────────── 논문 B
    │
    ├── 노트북 2: 세포 없이 분석 (P2R, ovrlpy, SSAM) ─ 논문 C
    │
    ├── 노트북 3: ★ 6개 기술 비교 (cross-platform) ─── 논문 D
    │     ├─ 3_1: 6개 기술 전부 Cellpose로 재분할
    │     ├─ 3_2: 뇌 영역별 세포 배정
    │     ├─ 3_3: Efficiency (검출 효율) ──────────── Fig 2b (box plot, 6개 기술)
    │     ├─ 3_4: NMP (특이성) ────────────────────── Fig 2d
    │     ├─ 3_5: Positivity ─────────────────────── Extended Data
    │     ├─ 3_6: Diffusion (확산) ────────────────── Fig 2f
    │     └─ 3_7: Xenium vs Visium ────────────────── Fig 2g
    │
    ├── 노트북 4: 확장 거리 최적화 ──────────────────── 논문 E (Fig 3a-b)
    │
    ├── 노트북 5: ★ Segmentation 벤치마크 ────────────── 논문 F
    │     └─ 5_1: Baysor vs Nuclei 비교
    │           ├─ UMAP ──────────────────────────── Fig 3f
    │           ├─ Counts violin ─────────────────── Fig 3g
    │           ├─ Cell type bar ─────────────────── Fig 3h
    │           └─ (+ ARI, NMP scatter) ──────────── Fig 3d-e
    │
    ├── 노트북 6: 전처리 시뮬레이션 ─────────────────── 논문 G (Fig 4a-e)
    ├── 노트북 7: 도메인 탐색 ───────────────────────── 논문 J (Fig 5g-h)
    └── 노트북 8: SVF 비교 ──────────────────────────── 논문 H (Fig 4f-j)
```

### 노트북 간 데이터 의존성

```
0_0 (Format) ──────────────────────────────────────────────────┐
  ├── 0_3 (Nuclei filter)                                      │
  │     └── 1_2 (Cell typing + domains) ──────────────────┐    │
  │           ├── 2_3 (ovrlpy overlaps)                   │    │
  │           ├── 2_4 (SSAM)                              │    │
  │           ├── 4_1 (Optimal expansion)                 │    │
  │           ├── 5_1 (Segmentation benchmark)            │    │
  │           ├── 7_1 (Domain detection) → 7_2            │    │
  │           └── 3_2 (Domain assign) ┐                   │    │
  │                                    │                   │    │
  ├── 1_1 (Dataset statistics) ─── Fig 1B                 │    │
  │                                    │                   │    │
  ├── 3_1 (Resegmentation 6 techs) ───┤                   │    │
  │     └── 3_6 (Diffusion)           │                   │    │
  │                                    │                   │    │
  │              3_2 (Domain assign) ──┤                   │    │
  │                ├── 3_3 (Efficiency)── 3_4 (NMP)       │    │
  │                ├── 3_5 (Positivity)                   │    │
  │                └── 3_7 (Xenium vs Visium)             │    │
  │                                                        │    │
  ├── 8_1 (SVF batch) → 8_2 (SVF comparison)              │    │
  │                                                        │    │
  └── 6_1 (scRNAseq download)                             │    │
        └── 6_3 (Simulate) → 6_4 (Assess)                 │    │
```

---

## 5. 우리 파이프라인 구조

```
Step 0: 포맷팅 ─────────────────────────────────────── 노트북 0
    │
Step 1: 데이터 탐색 + 세포 분류 ────────────────────── 노트북 1
    │
Step 2: 세포 없이 분석 (P2R, ovrlpy) ──────────────── 노트북 2
    │
Step 3: Cellpose 재분할 ───────────────────────────── 노트북 3_1 (일부만)
    │
Step 4: ⚠️ "기술 비교" ────────────────────────────── 노트북 3_3~3_7
    │     ├─ Efficiency (histogram, reseg vs original)
    │     ├─ NMP (reseg vs original)
    │     ├─ Positivity (reseg vs original)
    │     └─ Diffusion (reseg vs original)
    │
Step 5: 확장 거리 최적화 ──────────────────────────── 노트북 4
    │
Step 6: Segmentation 벤치마크 ─────────────────────── 노트북 5
    │     ├─ UMAP, Violin, Bar plot
    │     ├─ ARI, NMP, Assigned reads
    │     └─ Clustering quality (silhouette 등)
    │
Step 7: 전처리 시뮬레이션 ─────────────────────────── 노트북 6
```

### 파이프라인 데이터 흐름

```
Raw Xenium Output
       │
       ▼
   [Step 0: Formatting]
       │
       ├─→ step0_adata.h5ad ───────────────────────────────────┐
       ├─→ transcripts.csv ────────────────────────────────────┤
       │                                                       │
       ▼                                                       │
   [Step 1: Exploration]                                       │
       │                                                       │
       ├─→ step1_adata.h5ad ──┬─→ Step 4 (original baseline)  │
       │                      │                                │
       ▼                      │                                │
   [Step 2: Seg-Free]         │                                │
       │                      │                                │
       ├─→ domain_polygons ───┤                                │
       │                      │                                │
       ▼                      │                                │
   [Step 3: Resegmentation]   │                                │
       │                      │                                │
       ├─→ step3_adata.h5ad ──┼─→ Step 4 (resegmented) ───┐   │
       ├─→ step3_masks.tif ───┼───────────────────────────────┐│
       │                      │                            │  ││
       ▼                      │                            │  ││
   [Step 4: Comparison] ⚠️    │  (순수 분석, 데이터 출력 없음) │  ││
       │                      │                            │  ││
       ▼                      │                            │  ││
   [Step 5: Expansion]  ◄─────┘                            │  ││
       │                                                   │  ││
       ├─→ expanded_transcripts.csv                        │  ││
       │                                                   │  ││
       ▼                                                   ▼  ▼▼
   [Step 6: Benchmark] ◄── nuclei(S0) + cellpose(S3) + expansion(S5) + baysor(inline)
       │
       ▼
   [Step 7: Simulation] (독립 - Census 데이터 사용)
```

---

## 6. 핵심 문제: 뭐가 어떻게 꼬였는가

### 비유로 설명

```
논문 저자의 원래 의도:

  노트북 3 = "삼성 vs 애플 vs LG 카메라 화질 비교" (6개 브랜드 비교)
  노트북 5 = "같은 사진을 포토샵 A vs B로 보정"   (같은 사진, 다른 도구)

파이프라인이 한 것:

  Step 4 = "카메라 화질 비교 기준"으로 "포토샵 A vs B"를 비교 ⚠️
  Step 6 = 원래 의도대로 "포토샵 A vs B" 비교 ✅
```

### 구체적 문제점

```
┌──────────────────────────────────────────────────────────────────────┐
│                    논문의 원래 구조                                    │
│                                                                      │
│  노트북 3_3~3_7 ──→ 6개 SRT 기술 간 비교 (Fig 2) ← "cross-platform" │
│                     CosMx vs Xenium vs MERFISH vs ...                │
│                                                                      │
│  노트북 5_1    ──→ Baysor vs Nuclei 비교 (Fig 3)  ← "segmentation"  │
│                     같은 데이터 안에서 방법만 다르게                     │
└──────────────────────────────────────────────────────────────────────┘

                              ↓ 파이프라인 변환 과정에서...

┌──────────────────────────────────────────────────────────────────────┐
│                  파이프라인의 현재 구조                                 │
│                                                                      │
│  Step 4  ──→ 노트북 3_3~3_7 메트릭을 Reseg vs Original에 적용 ⚠️    │
│              (cross-platform용 메트릭을 segmentation 비교에 강제 적용) │
│                                                                      │
│  Step 6  ──→ 노트북 5_1 방식대로 segmentation 벤치마크 ✅             │
│              (이미 violin, bar, UMAP, ARI, NMP 전부 구현)             │
└──────────────────────────────────────────────────────────────────────┘
```

---

## 7. Step 4 vs Step 6 중복 분석

### Step 6이 이미 구현하고 있는 논문 Fig. 3 요소들

| 논문 Figure | Step 6 구현 | 코드 위치 |
|---|---|---|
| **Fig. 3g** (counts/cell violin) | `_save_counts_violin()` | L1352 |
| **Fig. 3h** (cell type frequency bar) | `_save_celltype_barplot()` | L1309 |
| **Fig. 3f** (UMAP by segmentation) | `_save_umap()` | L975 |
| **Fig. 3d** (ARI heatmap) | `_compute_rand_index()` | L1589 |
| **Fig. 3e** (NMP vs assigned reads scatter) | `_save_reads_vs_nmp_scatter()` | L1391 |
| + clustering quality | silhouette, CH, DB scores | L1896 |
| + spatial maps | `_save_spatial_map()` | L1256 |
| + DAPI overlays | `_save_dapi_*()` | L1036+ |

### Step 4 메트릭 판정

| Step 4 메트릭 | Step 6에 있음? | 필요? | 이유 |
|---|---|---|---|
| **4-2 Efficiency** (histogram) | ✅ violin + median_reads | **중복** | Step 6이 더 정확 |
| **4-3 NMP** | ✅ NMP per method + scatter | **중복** | Step 6이 더 정확 |
| **4-4 Positivity** | ❌ 없음 | **불필요** | cross-platform 전용 |
| **4-5 Diffusion** | ❌ 없음 | **불필요** | 논문도 segmentation 간 비교에 안 씀 |
| **4-2b Expression ratio** (ST/SC) | ❌ 없음 | **불필요** | cross-platform 전용 |

---

## 8. Figure 매핑 오류 상세

| 파이프라인 출력 파일 | 현재 문서 매핑 | 실제 해당 논문 Figure | 원본 노트북 |
|---|---|---|---|
| `efficiency_transcripts_per_cell_comparison.png` | ~~Fig. 2b~~ | **Fig. 3g** (violin) | **Notebook 5_1** |
| `efficiency_genes_per_cell_comparison.png` | ~~Fig. 2b~~ | **Fig. 3g** 변형 | **Notebook 5_1** |
| Fig. 2b (cross-platform box plot) | — | 파이프라인에 **미구현** | **Notebook 3_3** |

### 시각화 차이

| | 논문 Fig. 2b | 파이프라인 Step 4 | 논문 Fig. 3g |
|---|---|---|---|
| **Plot type** | Box plot | Histogram | Violin plot |
| **비교 대상** | 6개 SRT 플랫폼 | Reseg vs Original | Baysor vs Nuclei |
| **Metric** | Raw transcripts/cell | Raw transcripts/cell | Counts per cell |
| **데이터** | 6개 플랫폼 h5ad | 단일 Alzheimer's dataset | 단일 mouse brain |

---

## 9. 최종 정리: 무엇이 필요하고 무엇이 불필요한가

| 파이프라인 Step | 논문 섹션 | 매핑 상태 | 비고 |
|---|---|---|---|
| **Step 0** (포맷팅) | A | ✅ 정상 | |
| **Step 1** (탐색+분류) | B | ✅ 정상 | |
| **Step 2** (세포 없이 분석) | C | ✅ 정상 | |
| **Step 3** (Cellpose 재분할) | F 준비단계 | ✅ 정상 | |
| **Step 4** (4개 메트릭 비교) | ~~D~~ | ⚠️ **문제** | cross-platform 메트릭을 잘못 적용, Step 6과 중복 |
| **Step 5** (확장 최적화) | E | ✅ 정상 | |
| **Step 6** (벤치마크) | F | ✅ **정상** | 이미 Fig 3 전체 구현 |
| **Step 7** (시뮬레이션) | G | ✅ 정상 | |

### 미구현 논문 섹션

| 논문 섹션 | 구현 여부 | 이유 |
|---|---|---|
| D. Cross-platform 비교 | ❌ 미구현 | 6개 플랫폼 데이터 필요 (단일 데이터셋 파이프라인으로 불가) |
| H. SVF 비교 | ❌ 미구현 | 다수 SVF 알고리즘 + 다수 데이터셋 필요 |
| I. Gene imputation | ❌ 미구현 | 다수 imputation 알고리즘 벤치마크 필요 |
| J. Domain identification | ❌ 미구현 | 다수 domain 알고리즘 + manual annotation 필요 |

---

## 10. 올바른 파이프라인 흐름 (논문 기준)

```
Step 3 (resegmentation) → Step 5 (optimal expansion) → Step 6 (benchmark)
                                                          ↑ 여기서 비교 완료
```

Step 4는 Step 3 결과와 원본을 비교하는 "중간 검증" 역할을 하려 했지만,
Step 6이 이미 nuclei/cellpose/expansion/baysor 전체를 한꺼번에 비교하므로
Step 4를 거치지 않고 바로 Step 6에서 확인하는 것이 논문에 맞는 방식입니다.

---

## 부록: 논문 Figure → 노트북 → 파이프라인 Step 전체 매핑

| Paper Figure | 논문 섹션 | 노트북 | 파이프라인 Step | 상태 |
|---|---|---|---|---|
| Fig 1a-b | A (데이터 준비) | 1_1 | Step 0-1 | ✅ |
| Fig 1c-d | B (세포 분류) | 1_2 | Step 1 | ✅ |
| Fig 1e-k | C (seg-free) | 2_1, 2_3, 2_4 | Step 2 | ✅ |
| Fig 2a | D (cross-platform) | 3_1 | Step 3 (부분) | ⚠️ |
| **Fig 2b** | **D (efficiency)** | **3_3** | **Step 4 (잘못 적용)** | **⚠️ cross-platform 전용** |
| **Fig 2d** | **D (specificity)** | **3_4** | **Step 4 (잘못 적용)** | **⚠️ cross-platform 전용** |
| **Fig 2f** | **D (diffusion)** | **3_6** | **Step 4 (잘못 적용)** | **⚠️ cross-platform 전용** |
| Fig 2g | D (Xe vs Visium) | 3_7 | 미구현 | — |
| Fig 3a-b | E (expansion) | 4_1 | Step 5 | ✅ |
| Fig 3c | F (seg visual) | 5_1 | Step 6 (DAPI) | ✅ |
| **Fig 3d** | **F (ARI)** | **5_1** | **Step 6** | **✅** |
| **Fig 3e** | **F (NMP scatter)** | **5_1** | **Step 6** | **✅** |
| **Fig 3f** | **F (UMAP)** | **5_1** | **Step 6** | **✅** |
| **Fig 3g** | **F (violin)** | **5_1** | **Step 6** | **✅** |
| **Fig 3h** | **F (bar)** | **5_1** | **Step 6** | **✅** |
| Fig 4a-e | G (preprocessing) | 6_3, 6_4 | Step 7 | ✅ |
| Fig 4f-j | H (SVF) | 8_1, 8_2 | 미구현 | — |
| Fig 5a-f | I (imputation) | 10_1 | 미구현 | — |
| Fig 5g-h | J (domain) | 7_1, 7_2 | 미구현 | — |
