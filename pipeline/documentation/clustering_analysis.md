# Clustering Analysis: 뇌 세포 유형과 Leiden 클러스터 매핑

> **Dataset**: Xenium V1 FFPE Human Brain Alzheimer's (10x Genomics)
> **Resolution**: Leiden 1.0 | **총 세포 수**: 44,943 | **유전자 수**: 354
> **마커 선정**: Wilcoxon rank-sum test (각 클러스터 vs 나머지, adjusted p-value)

---

## 1. 인간 뇌 조직의 주요 세포 유형

인간 뇌는 크게 다음 세포 유형들로 구성된다:

| 세포 유형 | 영문 | 역할 | 뇌 내 비율 |
|-----------|------|------|-----------|
| 뉴런 (신경세포) | Neurons | 전기 신호 전달, 정보 처리 | ~50% (수적으로는 glia보다 적음) |
| 성상세포 | Astrocytes | BBB 유지, 대사 지원, 시냅스 조절 | glia 중 가장 많음 |
| 올리고덴드로사이트 | Oligodendrocytes | 수초(myelin) 형성, 축삭 절연 | white matter에 풍부 |
| 올리고덴드로사이트 전구세포 | OPCs | 올리고덴드로사이트로 분화 | 전체의 5-8% |
| 미세아교세포 | Microglia | 뇌 면역 감시, 시냅스 가지치기 | 전체의 5-10% |
| 혈관내피세포 | Endothelial cells | 혈관벽 구성, BBB 형성 | 소수 |
| 혈관주위세포 | Pericytes / Perivascular | 혈관 안정화 | 극소수 |

알츠하이머 뇌에서는 reactive astrogliosis(반응성 성상세포 증식), disease-associated microglia(DAM) 활성화, 올리고덴드로사이트 스트레스 등이 특징적으로 나타난다.

---

## 2. 클러스터별 세포 유형 어노테이션

### Cluster 0 — Mature Oligodendrocytes (성숙 올리고덴드로사이트)

| 항목 | 값 |
|------|-----|
| **세포 수** | 10,204 (22.7%) |
| **Top 3 마커** | CRYAB, SPP1, CA2 |
| **Top 5 (rank 4-5)** | QDPR, EDIL3 |
| **Median counts/genes** | 214 / 31 |

**마커 유전자 해석:**
- **CRYAB** (Crystallin Alpha B): 소분자 heat shock protein. White matter에 풍부하며 올리고덴드로사이트에서 높게 발현. 수초 유지 및 세포 스트레스 보호 역할. AD에서 아밀로이드 플라크 주변에 축적되는 것으로 알려짐.
- **SPP1** (Osteopontin): 세포외기질 당단백질. Disease-associated glia에서 상향조절. AD 뇌에서 플라크 주변 glia 활성과 관련.
- **CA2** (Carbonic Anhydrase 2): 탄산탈수효소. **올리고덴드로사이트의 canonical marker**로, 수초 형성 시 pH 조절에 핵심적 역할.
- **QDPR** (Dihydropteridine Reductase): BH4 재활용 효소. 올리고덴드로사이트에서 발현 보고.
- **EDIL3** (EGF-like Repeats and Discoidin I-like Domains 3): 혈관신생 및 세포 부착 관련.

**판정**: CA2가 교과서적 올리고덴드로사이트 마커이며, CRYAB도 white matter glia 특이적. 가장 큰 클러스터로, 뇌 조직의 올리고덴드로사이트 풍부 영역(white matter)을 반영한다.

---

### Cluster 1 — Excitatory Neurons (흥분성 뉴런)

| 항목 | 값 |
|------|-----|
| **세포 수** | 8,479 (18.9%) |
| **Top 3 마커** | SNCA, CLU, NRGN |
| **Top 5 (rank 4-5)** | BEX1, HSP90AA1 |
| **Median counts/genes** | 305 / 33 |

**마커 유전자 해석:**
- **SNCA** (Alpha-Synuclein): 시냅스 전 말단에 풍부한 단백질. 뉴런 특이적 마커. 시냅스 소포 순환과 도파민 조절에 관여. Parkinson's/Lewy body 병리와 직결되지만, 정상 뉴런에서도 높게 발현.
- **CLU** (Clusterin/ApoJ): **AD GWAS 상위 risk gene**. 분비형 샤페론으로 Aβ clearance에 관여. 뉴런과 astrocyte 모두에서 발현되지만, 이 클러스터에서 뉴런 마커들과 함께 나타남.
- **NRGN** (Neurogranin): **흥분성 뉴런의 대표 마커**. 수상돌기 스파인에 위치하며 CaMKII/calmodulin 신호전달 조절. CSF에서 AD 바이오마커로 활용됨.
- **BEX1** (Brain Expressed X-linked 1): 뇌 특이적 발현 유전자. 신경 발달과 분화에 관여.
- **HSP90AA1**: Heat shock protein. 스트레스 반응에 광범위하게 발현되지만 뉴런에서 높음.

**판정**: SNCA + NRGN 조합은 흥분성 뉴런을 강력하게 지지한다. CLU가 AD risk gene으로 함께 나타나는 것은 알츠하이머 뇌 샘플의 특성을 반영. 전체에서 두 번째로 큰 클러스터로, 뇌의 neuronal population을 대표한다. **Median counts가 305로 가장 높은 것**은 뉴런의 높은 전사 활성을 반영.

---

### Cluster 2 — Reactive Astrocytes (반응성 성상세포)

| 항목 | 값 |
|------|-----|
| **세포 수** | 6,450 (14.4%) |
| **Top 3 마커** | GFAP, MTRNR2L12, SPP1 |
| **Top 5 (rank 4-5)** | MOBP, CRYAB |
| **Median counts/genes** | 197 / 25 |

**마커 유전자 해석:**
- **GFAP** (Glial Fibrillary Acidic Protein): **성상세포의 가장 대표적인 canonical marker**. 중간섬유 단백질로 성상세포 세포골격의 핵심 구성. AD에서 반응성 성상세포(reactive astrogliosis) 시 극적으로 상향조절됨.
- **MTRNR2L12**: 미토콘드리아 rRNA 유사유전자. 생물학적 의미보다는 미토콘드리아 transcript 오염 가능성. (아래 주의사항 참조)
- **SPP1** (Osteopontin): Cluster 0에서도 나타남. 반응성 glia 전반에서 상향조절.
- **MOBP** (Myelin-associated Oligodendrocyte Basic Protein): 수초 관련 — 일부 올리고덴드로사이트 혼재 가능성.
- **CRYAB**: Cluster 0과 공유. 스트레스 반응 glia에서 공통 발현.

**판정**: GFAP가 1순위 마커로 나온 것은 성상세포 클러스터를 명확히 지시한다. 다만 MOBP, CRYAB 등 올리고 마커가 하위에 있어, 일부 glia 혼재 가능성이 있다. AD 뇌에서 reactive astrocyte는 GFAP 발현이 크게 증가하므로 별도 클러스터로 분리된 것이 합리적이다.

---

### Cluster 3 — Homeostatic Astrocytes (항상성 성상세포)

| 항목 | 값 |
|------|-----|
| **세포 수** | 5,673 (12.6%) |
| **Top 3 마커** | GJA1, AQP4, CLU |
| **Top 5 (rank 4-5)** | GLUL, ALDOC |
| **Median counts/genes** | 274 / 32 |

**마커 유전자 해석:**
- **GJA1** (Connexin 43): 성상세포 gap junction의 핵심 구성 단백질. 성상세포 네트워크 형성과 이온/대사물 교환에 필수적. **성상세포 canonical marker**.
- **AQP4** (Aquaporin 4): 성상세포 endfeet에 특이적으로 발현되는 수분 채널. 뇌척수액 순환(glymphatic system)의 핵심. AD에서 AQP4 재분포가 Aβ clearance 저하와 연관.
- **CLU** (Clusterin): Cluster 1과 공유. 성상세포에서도 분비.
- **GLUL** (Glutamine Synthetase): 글루타메이트를 글루타민으로 전환. **성상세포 특이적 효소**. 시냅스 글루타메이트 재활용의 핵심.
- **ALDOC** (Aldolase C/Fructose-bisphosphate Aldolase C): 해당과정 효소. **성상세포 특이적 동형체**로, astrocyte-neuron lactate shuttle에 관여.

**판정**: **가장 깔끔한 클러스터.** GJA1, AQP4, GLUL, ALDOC 전부가 교과서적 성상세포 마커다. Cluster 2(reactive)와 달리 항상성 기능 마커들이 주를 이루어, **homeostatic astrocyte**로 분류된다. 두 성상세포 클러스터의 분리는 AD 뇌에서 reactive vs homeostatic 상태 공존을 반영한다.

---

### Cluster 4 — Endothelial / Vascular Cells (혈관내피/혈관계 세포)

| 항목 | 값 |
|------|-----|
| **세포 수** | 5,529 (12.3%) |
| **Top 3 마커** | HSPA1B, FLT1, VIM |
| **Top 5 (rank 4-5)** | MTRNR2L8, IFITM3 |
| **Median counts/genes** | 194 / 24 |

**마커 유전자 해석:**
- **HSPA1B** (HSP70 family): Heat shock protein. 스트레스 반응에서 광범위하게 발현되지만, 혈관내피 세포에서도 높은 발현 보고.
- **FLT1** (VEGFR1, Vascular Endothelial Growth Factor Receptor 1): **혈관내피세포의 canonical marker**. VEGF 신호를 수용하여 혈관신생 조절. BBB(혈뇌장벽) 구성 세포에서 핵심적.
- **VIM** (Vimentin): 중간섬유 단백질. 내피세포, 중간엽 세포, 반응성 glia에서 발현. 혈관계 세포의 구조 단백질.
- **IFITM3** (Interferon-induced Transmembrane Protein 3): 선천면역 반응 단백질. 내피세포와 면역세포에서 발현. AD에서 γ-secretase 조절자로 주목받음.

**판정**: FLT1이 핵심 마커로, 혈관내피세포 클러스터로 확실하다. VIM + IFITM3 조합은 BBB 구성 세포 및 혈관 주변 면역 반응을 반영. 뇌 조직의 혈관계(vasculature)를 대표한다.

---

### Cluster 5 — Microglia / Immune Cells (미세아교세포/면역세포)

| 항목 | 값 |
|------|-----|
| **세포 수** | 4,512 (10.0%) |
| **Top 3 마커** | MTRNR2L12, MTRNR2L8, GPNMB |
| **Top 5 (rank 4-5)** | PTPRC, LYVE1 |
| **Median counts/genes** | 139 / 18 |

**마커 유전자 해석:**
- **MTRNR2L12, MTRNR2L8**: 미토콘드리아 rRNA 유사유전자. 생물학적 마커라기보다 technical artifact 가능성이 높음 (아래 주의사항 참조). 다만 이 유전자들의 log fold change가 매우 낮아(0.11, 0.08) 실질적 차별 발현은 아님.
- **GPNMB** (Glycoprotein NMB): **Disease-Associated Microglia (DAM)의 핵심 마커**. AD에서 활성화된 미세아교세포에서 특이적으로 상향조절. 리소좀 기능, 식세포작용, 지질 대사 관련.
- **PTPRC** (CD45): **Pan-immune marker**. 모든 조혈계 면역세포에서 발현. 뇌에서는 주로 미세아교세포와 침윤 면역세포를 표지.
- **LYVE1** (Lymphatic Vessel Endothelial Hyaluronan Receptor 1): 뇌막 림프관 내피세포와 **perivascular macrophage (혈관주위 대식세포)** 마커.

**판정**: GPNMB + PTPRC + LYVE1 조합은 미세아교세포/면역세포 클러스터를 강력히 지지한다. 특히 GPNMB는 AD-specific DAM 마커로, 알츠하이머 뇌의 면역 활성 상태를 반영. Median counts(139)와 genes(18)가 가장 낮은 것은 미세아교세포의 상대적으로 작은 세포 크기와 낮은 transcript 수를 반영한다. **MTRNR 유전자 제외 시 GPNMB, PTPRC, LYVE1이 top 마커로 올라옴** (noMTRNR dotplot 참조).

---

### Cluster 6 — Stressed Oligodendrocytes (스트레스 올리고덴드로사이트)

| 항목 | 값 |
|------|-----|
| **세포 수** | 3,911 (8.7%) |
| **Top 3 마커** | CRYAB, QDPR, SPP1 |
| **Top 5 (rank 4-5)** | MTRNR2L12, FTH1 |
| **Median counts/genes** | 148 / 21 |

**마커 유전자 해석:**
- **CRYAB**: Cluster 0과 공유. 스트레스 보호 역할.
- **QDPR**: Cluster 0과 공유. 올리고덴드로사이트 계열.
- **SPP1**: Disease-associated glia 마커.
- **FTH1** (Ferritin Heavy Chain 1): 철 저장 단백질. 올리고덴드로사이트는 뇌에서 철 함량이 가장 높은 세포 유형. AD에서 철 항상성 이상(iron dysregulation)과 연관.

**판정**: Cluster 0과 마커가 겹치지만 별도로 분리된 하위집단. p-value와 log fold change가 Cluster 0보다 약하고, median counts/genes도 낮다. **스트레스 상태이거나 탈수초화(demyelination) 과정에 있는 올리고덴드로사이트 하위집단**으로 해석된다. FTH1의 등장은 AD에서의 철 대사 이상을 시사한다.

---

### Cluster 7 — OPCs (올리고덴드로사이트 전구세포)

| 항목 | 값 |
|------|-----|
| **세포 수** | 185 (0.4%) |
| **Top 3 마커** | VCAN, BCAN, SEMA5A |
| **Top 5 (rank 4-5)** | OLIG2, DNER |
| **Median counts/genes** | 259 / 40 |

**마커 유전자 해석:**
- **VCAN** (Versican): 대형 chondroitin sulfate proteoglycan. **OPC의 대표적 세포외기질 마커**. OPC가 분비하여 이동과 분화를 조절하는 환경을 형성.
- **BCAN** (Brevican): 뇌 특이적 proteoglycan. **OPC 및 perineuronal net의 핵심 구성요소**. 시냅스 가소성과 축삭 성장 조절.
- **SEMA5A** (Semaphorin 5A): 축삭 유도 분자. OPC의 이동(migration)과 위치 결정에 관여. 수초화 과정의 조절자.
- **OLIG2** (Oligodendrocyte Transcription Factor 2): **올리고덴드로사이트 계열의 master 전사인자**. OPC에서 성숙 올리고까지 전 과정에서 발현되지만, OPC에서 특히 높음.
- **DNER** (Delta/Notch-like EGF Repeat Containing): Notch 신호 경로 관련. 신경-glia 상호작용과 분화 조절.

**판정**: **매우 깔끔한 OPC 시그니처**. VCAN + BCAN + OLIG2 + SEMA5A 조합은 OPC를 의심의 여지 없이 가리킨다. 전체의 0.4%(185개)로 극소수이지만, 마커의 log fold change가 매우 높아(VCAN: 7.87, BCAN: 6.70, OLIG2: 6.51) 통계적으로 매우 유의미하다. Median genes(40)가 모든 클러스터 중 가장 높은 것은 OPC의 활발한 전사 활성(분화 준비 상태)을 반영한다.

---

## 3. 전체 요약: 클러스터 → 세포 유형 매핑

| Cluster | 세포 유형 | 세포 수 (%) | 핵심 마커 | 신뢰도 |
|---------|----------|------------|----------|--------|
| 0 | Mature Oligodendrocytes | 10,204 (22.7%) | CA2, CRYAB | 높음 |
| 1 | Excitatory Neurons | 8,479 (18.9%) | SNCA, NRGN | 높음 |
| 2 | Reactive Astrocytes | 6,450 (14.4%) | GFAP | 높음 |
| 3 | Homeostatic Astrocytes | 5,673 (12.6%) | GJA1, AQP4, GLUL, ALDOC | 매우 높음 |
| 4 | Endothelial / Vascular | 5,529 (12.3%) | FLT1, VIM | 높음 |
| 5 | Microglia / Immune | 4,512 (10.0%) | GPNMB, PTPRC, LYVE1 | 높음 |
| 6 | Stressed Oligodendrocytes | 3,911 (8.7%) | CRYAB, QDPR, FTH1 | 중간 |
| 7 | OPCs | 185 (0.4%) | VCAN, BCAN, OLIG2 | 매우 높음 |

---

## 4. 공간적 분포 패턴 (Neighborhood Enrichment)

Neighborhood Enrichment Z-score 히트맵에서 관찰되는 주요 패턴:

- **Cluster 0 자가응집 (z >> 100)**: 올리고덴드로사이트가 white matter 영역에 밀집. 공간적으로 가장 강한 자가응집을 보임.
- **Cluster 0 ↔ 1 배제 (z < 0)**: 올리고(white matter)와 뉴런(gray matter)이 공간적으로 분리 — 해부학적으로 정확한 패턴.
- **Cluster 2 ↔ 0 양의 상관**: Reactive astrocyte가 올리고덴드로사이트 근처에 분포 — white matter gliosis를 반영.
- **Cluster 3 ↔ 4 양의 상관**: 항상성 성상세포가 혈관 주변에 위치 — astrocyte endfeet-BBB 구조를 반영 (AQP4가 BBB endfeet 마커).
- **Cluster 5 ↔ 4 양의 상관**: 미세아교세포가 혈관 근처에 분포 — perivascular immune surveillance.

Spatial Map에서도 Cluster 0(파란색, 올리고)이 조직 하단-좌측에 집중, Cluster 1(뉴런)이 상단-우측에 분포하는 등 gray matter/white matter 경계가 관찰된다.

---

## 5. 알츠하이머 특이적 소견

이 클러스터링 결과에서 AD 관련 주요 소견:

1. **성상세포가 두 클러스터로 분리 (2 vs 3)**: Reactive(GFAP↑) vs Homeostatic(GJA1, AQP4) 상태 공존은 AD의 reactive astrogliosis를 직접 반영한다.

2. **올리고덴드로사이트도 두 클러스터로 분리 (0 vs 6)**: Cluster 6의 FTH1(철 저장)과 낮은 transcript 수는 AD에서의 탈수초화(demyelination) 및 iron dysregulation을 시사한다.

3. **DAM 마커 GPNMB의 등장 (Cluster 5)**: Disease-Associated Microglia가 별도 클러스터로 분리될 만큼 활성화되어 있다.

4. **CLU가 뉴런 클러스터에서 2순위 (Cluster 1)**: AD GWAS top risk gene이 뉴런에서 높게 발현되는 것은 Aβ clearance 부담을 반영.

5. **IFITM3 (Cluster 4)**: 최근 AD 연구에서 γ-secretase 조절자로 보고된 유전자가 혈관 클러스터에서 나타남 — neurovascular unit에서의 AD 병리 관여 가능성.

---

## 6. 주의사항: MTRNR 유전자

MTRNR2L12, MTRNR2L8 등 미토콘드리아 rRNA 유사유전자(pseudogene)가 모든 클러스터에서 top expressed gene으로 나타난다. 이는:

- Xenium 패널에 포함된 프로브가 미토콘드리아 transcript를 잡는 것으로, **기술적 특성(technical artifact)**일 가능성이 높다.
- Cluster 5에서 top 마커로 나타났지만 **log fold change가 0.11/0.08로 극히 낮아**, 실질적 차별 발현이 아니다.
- 파이프라인에서 **MTRNR 제외 버전 dotplot**(`_noMTRNR`)을 별도 생성하여 이 문제를 보완하고 있다.
- MTRNR 제외 시 Cluster 5의 진정한 마커인 GPNMB, PTPRC, LYVE1이 상위로 올라온다.

---

## 7. 클러스터링 품질 평가

**긍정적 평가:**
- 뇌의 6대 주요 세포 유형(Neuron, Astrocyte, Oligodendrocyte, OPC, Microglia, Endothelial)이 **모두 포착됨**
- Astrocyte와 Oligodendrocyte가 각각 2개씩 하위 상태(reactive/homeostatic, mature/stressed)로 분리 — 생물학적으로 의미 있는 해상도
- Cluster 3(Astrocyte)과 7(OPC)의 마커가 **교과서적으로 깔끔**
- Spatial 분포가 해부학적 구조(gray/white matter 분리)와 일치
- 전체 세포 유형 비율이 알려진 뇌 조직 구성과 대체로 일치

**개선 가능 사항:**
- Inhibitory neuron(억제성 뉴런, GABAergic)이 별도 클러스터로 분리되지 않음 — 더 높은 resolution이 필요할 수 있음
- Cluster 6(Stressed Oligo)의 정체성이 Cluster 0과 겹침 — 추가 검증 필요
- MTRNR 유전자가 마커 선정에 노이즈를 줌 — 사전 필터링 권장
