# MTRNR2L12 / MTRNR2L8 유전자 분석 보고서

> 작성일: 2026-02-27
> 데이터셋: Xenium V1 FFPE Human Brain Alzheimer's (With Addon)

---

## 1. Gene Panel 정보

| 항목 | MTRNR2L12 | MTRNR2L8 |
|------|-----------|----------|
| **Ensembl ID** | ENSG00000269028 | ENSG00000255823 |
| **Codeword** | 432 | 421 |
| **Probe 수 (coverage)** | 8 | 8 |
| **Panel 종류** | **Addon** (hBrain_100g) | **Addon** (hBrain_100g) |
| **Descriptor** | "gene" (일반 유전자, 컨트롤 아님) | "gene" |

### Panel 구성
- **Panel 이름:** PD_349_LX16_hBrain_100g_PD326-addon-alzheimers
- **Base panel:** hBrain_254g (254개 유전자)
- **Addon panel:** hBrain_100g (100개 추가 유전자, 알츠하이머 특화)
- **총 타겟:** 354개 유전자
- **Chemistry:** DE_24 (10x Xenium)
- **Gene panel 파일:** `data/Xenium_V1_FFPE_Human_Brain_Alzheimers_With_Addon_outs/gene_panel.json`

---

## 2. 생물학적 정체

- **Mitochondrial rRNA-like pseudogenes** (미토콘드리아 ribosomal RNA 유사 위유전자)
- 실제 미토콘드리아 유전자(MT-RNR2, 16S rRNA)에서 유래한 핵 내 복제본 (NUMTs: Nuclear Mitochondrial DNA segments)
- 모든 세포 타입에서 높게 발현 (housekeeping-like, 세포 타입 특이성 없음)
- Xenium에서 특히 지배적인 이유: 미토콘드리아가 세포당 수백~수천 개 존재하여 transcript copy number가 극히 높음

---

## 3. 이 데이터에서의 발현 비율

| 지표 | 값 |
|------|-----|
| **MTRNR2L12** | 전체 reads의 **45.22%** (약 4.77M / 10.54M counts) |
| **MTRNR2L8** | 전체 reads의 **31.00%** (약 3.27M / 10.54M counts) |
| **Combined** | **76.22%** |
| **Per-cell 중앙값 (MTRNR2L12)** | ~47% of cell total counts |
| **Per-cell 중앙값 (MTRNR2L8)** | ~30% of cell total counts |
| **Dispersion (핵 거리)** | MTRNR2L12: 9.67 µm, MTRNR2L8: 9.31 µm (cytoplasmic) |

> 354개 유전자 중 단 2개가 전체 reads의 3/4을 차지

### 시각화 출력
- `human_alzheimers_step1_mtrnr_proportion.png` (좌: Top 20 barplot, 우: per-cell violin)
- `human_alzheimers_step1_dispersion_localization.png` (nuclear vs cytoplasmic 분류에서 cytoplasmic 극우측)

---

## 4. Clustering/UMAP에 대한 영향

### 문제: Raw counts PCA 시 MTRNR이 지배

정규화 없이 PCA를 수행하면:
- PC1이 사실상 "MTRNR 발현량 축"이 됨
- 모든 downstream (neighbors, UMAP, Leiden)이 왜곡
- Leiden res=1.0에서 17개 클러스터 생성, 대부분의 top marker가 MTRNR2L12/L8

### 해결: normalize_total + log1p 전처리

정규화 적용 후:
- **셀별 library size 정규화** (target_sum=1e4): MTRNR 고발현 셀의 영향 완화
- **log1p**: 고발현 유전자 효과 압축 (log 스케일에서 MTRNR의 지배력 감소)
- Leiden res=1.0에서 **8개 클러스터**, 생물학적 마커가 정상적으로 검출

| 비교 항목 | 수정 전 (raw counts PCA) | 수정 후 (normalize + log1p) |
|-----------|-------------------------|---------------------------|
| **클러스터 수** (res=1.0) | 17개 (과도 분절) | **8개** (적절) |
| **UMAP 구조** | blob 형태, 분리 불량 | 클러스터간 분리 명확 |
| **Top markers** | MTRNR2L12/L8이 거의 모든 클러스터의 top marker | CRYAB, SNCA, GFAP, GJA1 등 생물학적 마커 |

---

## 5. 파이프라인에서의 처리 전략

### Step 1: Dataset Exploration
- `normalize_total(target_sum=1e4)` + `log1p`로 발현량 압축 → 클러스터링 왜곡 완화
- MTRNR proportion 시각화로 문제 크기 정량화
- MTRNR 제거 없이도 정규화만으로 클러스터링 품질 개선

### Step 5: Optimal Expansion Analysis
- Background reads의 **15% 초과 유전자 자동 필터** (`bck_frac > 0.15`)
- MTRNR 제거 시 correlation dynamic range: ~0.03 → ~0.20으로 확대
- 세포 타입별 turnover 패턴이 드러남

### HVG Selection
- `highly_variable_genes(min_mean=0.3, max_mean=7)` 파라미터에서 `max_mean=7` 조건에 의해 MTRNR 유전자(mean ~8.4)가 **HVG에서 자동 제외**
- PCA는 HVG 기반으로 수행되므로 MTRNR의 직접적 영향 차단

---

## 6. 정규화 후 클러스터 마커 해석

| Cluster | Top 3 Markers | 추정 Cell Type | 셀 수 (%) |
|---------|--------------|---------------|-----------|
| 0 | CRYAB, SPP1, CA2 | Oligodendrocytes | 10,204 (22.7%) |
| 1 | SNCA, CLU, NRGN | Neurons (excitatory) | 8,479 (18.9%) |
| 2 | GFAP, MTRNR2L12, SPP1 | Astrocytes (reactive) | 6,450 (14.4%) |
| 3 | GJA1, AQP4, CLU | Astrocytes (protoplasmic) | 5,673 (12.6%) |
| 4 | HSPA1B, FLT1, VIM | Endothelial / Stress | 5,529 (12.3%) |
| 5 | MTRNR2L12, MTRNR2L8, GPNMB | Microglia/Immune | 4,512 (10.0%) |
| 6 | CRYAB, QDPR, SPP1 | Oligodendrocytes (subtype) | 3,911 (8.7%) |
| 7 | VCAN, BCAN, SEMA5A | OPC (Oligo Progenitor) | 185 (0.4%) |

> Cluster 5는 여전히 MTRNR이 top marker로 남아 있음. 이 클러스터는 상대적으로 MTRNR 고발현 + 낮은 다른 유전자 발현 셀로 구성되며, GPNMB/PTPRC/LYVE1 등 면역 세포 마커를 동시에 발현. MTRNR 제거 후 re-clustering 시 면역 세포 클러스터로 더 깨끗하게 분리 가능.

---

## 7. 참고 문헌 및 관련 파일

- **Gene panel:** `data/Xenium_V1_FFPE_Human_Brain_Alzheimers_With_Addon_outs/gene_panel.json`
- **MTRNR proportion plot:** `xenium-output/.../step1_exploration/human_alzheimers_step1_mtrnr_proportion.png`
- **Localization violin:** `xenium-output/.../step1_exploration/human_alzheimers_step1_dispersion_localization.png`
- **Step 5 background filter:** `pipeline/xenium_step5_optimal_expansion.py` (lines 298-311)
- **Paper:** Salas et al., "Optimizing Xenium In Situ data utility", Nature Methods, April 2025
