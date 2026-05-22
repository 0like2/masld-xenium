# Step 5: 최적 확장 — 초등학생도 알기 쉬운 상세 가이드

> **작성일**: 2026-02-23
> **목적**: Step 5의 전체 과정을 누구나 이해할 수 있도록 상세히 풀어 설명

---

## 먼저, 왜 이걸 하는 걸까?

```
세포 = 풍선
풍선 안에 = 핵(nucleus) ← Cellpose가 찾은 것
풍선 밖 바닥에 = 떨어진 쪽지(transcript) ← 아직 주인 없음

문제: 쪽지가 핵 밖에 떨어져 있다!
     → 누구 세포 건지 모른다!

해결: 핵을 조금 부풀려서(expand) 주변 쪽지를 주워오자!
     → 그런데 얼마나 부풀려야 할까?
     → 너무 적게 부풀리면 → 쪽지를 못 주움
     → 너무 많이 부풀리면 → 옆집 세포의 쪽지를 훔침!
```

**Step 5의 목표: "딱 알맞게" 부풀리는 거리를 찾는 것!**

---

## 전체 과정 한눈에 보기

```
[5-1] 데이터 불러오기 (쪽지 전부 + 세포 정보)
  ↓
[5-2] 각 세포가 어느 "동네"에 있는지 정하기 (Domain Mapping)
  ↓
[5-3] 주인 없는 쪽지를 가장 가까운 세포에 배달하기 (KDTree)
  ↓
[5-4] 거리 필터 (너무 먼 쪽지는 배달 취소)
  ↓
[5-5] 그림 그리기 (어디에 뭐가 있는지 지도)
  ↓
[5-6] 확장된 쪽지 저장
  ↓
[5-7] Turnover 분석 (핵심!!) — "어디까지가 내 쪽지이고, 어디부터 남의 쪽지인지 찾기"
  ↓
[결과] 최적 확장 거리 = turnover 거리 - 핵 반경
```

---

## [5-1] 데이터 불러오기

### 뭘 불러오나?

**두 가지**를 불러옵니다:

#### (1) 전체 transcript (쪽지) — Step 0에서 만든 것

```python
# Step 0에서 만든 h5ad 파일을 열어서
adata_step0 = sc.read_h5ad(step0_file)

# 그 안에 있는 전체 transcript 테이블을 꺼냄
reads_original = adata_step0.uns['spots']
```

이 테이블은 이렇게 생겼습니다:

```
transcript_id | cell_id | x_location | y_location | feature_name
─────────────┼─────────┼────────────┼────────────┼─────────────
  tx_00001   |   352   |  1234.5    |  5678.2    |  GFAP
  tx_00002   |    -1   |  1240.3    |  5680.1    |  MBP     ← 주인 없음! (-1)
  tx_00003   |   352   |  1235.1    |  5679.0    |  AQP4
  tx_00004   |    -1   |  2100.0    |  3200.5    |  SYT1    ← 주인 없음!
  ...
```

- `cell_id = 352` → 352번 세포에 속한 쪽지
- `cell_id = -1` → **주인 없는 쪽지** (이것들을 배달해야 함!)

#### (2) 세포 정보 — Step 1 또는 Step 3에서 만든 것

```python
# Step 3의 resegmented h5ad, 또는 Step 1의 exploration h5ad
adata_annotated = sc.read_h5ad(input_adata_file)
```

이 안에는:
- 각 세포의 **위치** (x_centroid, y_centroid)
- 각 세포의 **종류** (celltype: 뉴런, 성상세포, 미세아교세포...)
- 각 세포의 **동네** (domain: 이 세포가 뇌의 어느 영역에 있는지)

### 만약 세포 종류(celltype)가 없으면?

scRNA-seq 레퍼런스 데이터를 가져와서 **kNN 라벨 전달**을 합니다 (line 580-668):

```
쉽게 말하면:
"이 세포가 어떤 유전자를 많이 발현하는지 보고,
 scRNA-seq에서 가장 비슷한 세포를 15개 찾아서,
 그 중 가장 많은 유형으로 이름을 붙여주자!"

예: 15개 중 12개가 "Astrocyte"이면 → 이 세포도 "Astrocyte"
    confidence = 12/15 = 0.80
```

---

## [5-2] 각 세포가 어느 "동네"에 있는지 정하기 (Domain Mapping)

### 왜 "동네"가 필요한가?

```
뇌 조직을 위에서 내려다보면:

    ┌───────────────────────────────┐
    │  피질(Cortex)                 │  ← 동네 A
    │   뉴런이 많은 곳              │
    ├───────────────────────────────┤
    │  해마(Hippocampus)            │  ← 동네 B
    │   기억과 관련된 곳             │
    ├───────────────────────────────┤
    │  시상(Thalamus)               │  ← 동네 C
    │   신호 중계소                  │
    └───────────────────────────────┘

같은 "동네"의 배경 쪽지들은 비슷한 패턴을 가짐
→ 동네별로 "배경 신호"를 따로 계산해야 정확함!
```

### 4단계 우선순위로 동네를 정함

코드는 아래 순서로 동네 정보를 찾습니다 (line 804-846):

```
우선순위 1: spatial_annotation (이미 누가 표시해둔 영역 정보)
  ↓ 없으면
우선순위 2: P2R Domain Mapping (Step 2에서 만든 Points2Regions 사용)
  ↓ 없으면
우선순위 3: leiden (유전자 발현 기반 클러스터링 결과)
  ↓ 없으면
우선순위 4: 새로 leiden 클러스터링 돌림 (최후의 수단)
```

### 우선순위 2: P2R Domain Mapping 상세 (`_load_p2r_spatial_domains()`, line 60-177)

이것이 가장 핵심적인 동네 배정 방법입니다.

#### P2R이 뭔가?

Step 2에서 **Points2Regions**이라는 도구로 transcript 밀도를 분석해서, 조직을 작은 정사각형 격자(bin)로 나누고, 비슷한 bin끼리 묶어서 클러스터를 만들었습니다.

```
조직을 격자로 나눈 모습:

    ┌──┬──┬──┬──┬──┬──┐
    │A │A │A │B │B │B │
    ├──┼──┼──┼──┼──┼──┤
    │A │A │B │B │B │C │     A, B, C = P2R 클러스터
    ├──┼──┼──┼──┼──┼──┤     (비슷한 transcript 패턴끼리 묶음)
    │A │B │B │C │C │C │
    ├──┼──┼──┼──┼──┼──┤
    │B │B │C │C │C │C │
    └──┴──┴──┴──┴──┴──┘
```

#### 과정:

**Step a: P2R 파일 찾기**

```python
# "k50" 파일을 먼저 찾고, 없으면 가장 작은 k 파일을 선택
p2r_path = f"{sample_tag}_step2_points2regions_k{n_clusters}_bins.h5ad"

# 못 찾으면 → k=100, k=200 등 가장 작은 걸 자동으로 탐색
```

**Step b: Meta-clustering (필요시)**

요청한 k(예: 50)보다 파일의 k(예: 200)가 크면, **너무 잘게 쪼개져 있으니 합치는 작업**을 합니다:

```python
# bin의 공간 좌표(x, y)를 기준으로 KMeans 클러스터링
# → 물리적으로 가까운 bin끼리 묶어서 50개 "큰 동네"로 합침
km = KMeans(n_clusters=50, random_state=42)
spatial_labels = km.fit_predict(bin_coords)
```

```
합치기 전 (k=200):           합치기 후 (k=50):
┌──┬──┬──┬──┬──┬──┐         ┌──────────┬──────────┐
│23│24│25│67│68│69│         │          │          │
├──┼──┼──┼──┼──┼──┤         │ region_0 │ region_1 │
│23│24│67│67│69│89│   →     │          │          │
├──┼──┼──┼──┼──┼──┤         ├──────────┼──────────┤
│24│67│67│89│89│89│         │ region_2 │ region_3 │
└──┴──┴──┴──┴──┴──┘         └──────────┴──────────┘
200개의 작은 조각            50개의 큰 동네
```

**Step c: 세포를 가장 가까운 P2R bin에 매핑**

```python
# P2R bin들의 좌표로 KDTree 구축
tree = cKDTree(bin_coords)   # bin 위치로 나무(tree) 만들기

# 각 세포의 중심점(centroid)에서 가장 가까운 bin 찾기
dists, idxs = tree.query(cell_coords, k=1)

# 해당 bin의 클러스터 라벨 = 세포의 동네!
domains = bin_clusters[idxs]

# 너무 먼 세포(>100um)는 동네 없음(NaN)으로 처리
domains[dists > max_dist_um] = np.nan
```

```
비유로 설명:

    P2R bin = 동네 표지판
    세포 = 아이

    각 아이가 가장 가까운 표지판을 찾아가면 → 그 아이의 동네가 정해짐!
    단, 표지판에서 100m 이상 떨어진 아이는 → "동네 미정"
```

**Step d: 좌표 단위 자동 감지**

세포 좌표가 pixel이고 P2R bin이 um이면 규모가 다르므로 자동 변환:

```python
# 세포 좌표 범위가 bin 좌표 범위의 2배 이상이면 → pixel→um 변환
if cell_range > bin_range * 2:
    cell_coords = cell_coords / scale  # 예: pixel / 4.7 = um
```

---

## [5-3] 주인 없는 쪽지를 가장 가까운 세포에 배달하기 (KDTree)

### "읽기(transcript) 분류" — 주인 있음 vs 주인 없음

```python
# domain이 있는 transcript = 세포에 소속됨 (주인 있음)
annotatedcells = reads_original[reads_original['domain'].notna()]

# domain이 없는 transcript = 미할당 (주인 없음)
nancells = reads_original[reads_original['domain'].isna()]
```

```
예시:
  주인 있는 쪽지: 50만 개 (세포에 이미 소속)
  주인 없는 쪽지: 20만 개 (바닥에 떨어져 있음)
  → 이 20만 개를 배달해야 함!
```

### KDTree란?

```
KDTree = "주변 검색 나무"

보통 방법: 모든 세포와 거리를 일일이 계산 → 매우 느림
KDTree 방법: 공간을 나무처럼 쪼개서 빠르게 검색 → 매우 빠름

비유:
  택배 기사가 주소를 찾을 때
  X 온 동네를 다 돌아다니며 찾기 (느림)
  O "이 지역 → 이 블록 → 이 집" 순으로 좁혀가기 (빠름)
```

### 서브샘플링 — 왜 전부 안 쓰나?

```python
subsample_frac = 0.01  # 전체의 1%만 사용
n_sample = int(n_assigned * subsample_frac)  # 50만 x 0.01 = 5000개

# 5000개의 "대표 쪽지"만 뽑아서 KDTree를 만듦
annotated_sub = annotatedcells.sample(n=n_sample, random_state=42)
tree = cKDTree(annotated_sub[['x_location', 'y_location']].values)
```

```
비유:
  50만 개의 집 주소를 전부 나무에 넣으면 → 메모리 폭발
  5000개의 대표 주소만 넣어도 → 주변 동네를 충분히 찾을 수 있음
  (1%만 써도 결과는 거의 동일)
```

### 배달 실행!

```python
# 주인 없는 쪽지 20만 개의 좌표
coords2 = nancells[['x_location', 'y_location']].values

# KDTree에 물어보기: "이 쪽지에서 가장 가까운 대표 쪽지는 누구?"
dists, idxs = tree.query(coords2, k=1)

# 가장 가까운 대표 쪽지의 동네 = 이 쪽지의 동네!
assigned_domains = anchor_domains[idxs]
```

```
비유:
  쪽지: "나는 (1240, 5680) 좌표에 있어!"
  KDTree: "가장 가까운 대표 쪽지는... 352번 세포의 '피질' 동네 것이네!"
  쪽지: "그러면 나도 '피질' 동네로 배달!"
```

100만 개가 넘으면 **100만 개씩 나눠서(chunk)** 처리:

```python
chunk_size = 1000000
for i in range(n_chunks):
    chunk = coords2[start:end]
    dists, idxs = tree.query(chunk, k=1)
```

---

## [5-4] 거리 필터 (너무 먼 쪽지는 배달 취소)

```python
if dist_threshold:  # 예: 15um
    mask_too_far = concat_distances > dist_threshold
    assigned_domains[mask_too_far] = np.nan  # "너무 멀어! 배달 취소!"
```

```
비유:
  "가장 가까운 집이 15m 이상 떨어져 있으면
   그건 아무 집의 것도 아닌 거야.
   길바닥에 떨어진 쓰레기(배경 노이즈)일 가능성이 높아!"
```

(기본값은 `null` = 거리 제한 없이 무조건 가장 가까운 곳에 배달)

---

## [5-5] 그림 그리기

### (1) Expansion Map (`expansion_map.png`) → 논문 Fig. 3a 좌측

모든 transcript를 지도 위에 **동네별 색상**으로 표시:

```python
sns.scatterplot(data=plot_df, x='x_location', y='y_location',
                hue='domain', s=1)  # 각 점 = 1개 transcript, 색 = 동네
```

```
    ┌───────────────────────────────┐
    │ ** **   ++ ++ ++              │  * = 동네A (피질)
    │ * * *   + + ++ ++             │  + = 동네B (해마)
    │ * *    ++ ++  ## ##           │  # = 동네C (시상)
    │ *     ++ +   ## ## ##         │
    │       ++    ## ## ## ##       │  전에는 주인없던 쪽지도 이제
    └───────────────────────────────┘  색깔이 있다 = 동네가 정해짐!
```

### (2) Distance Histogram (`expansion_distances.png`)

미할당 쪽지 → 가장 가까운 세포까지의 거리 분포:

```python
sns.histplot(concat_distances, bins=50, kde=True)
```

```
    빈도
    |  ####
    |  ######
    |  ########
    |  ##############
    |  ####################....
    └──────────────────────────────→ 거리 (um)
      0    5   10   15   20   25

    "대부분의 주인없는 쪽지가 세포에서 5um 이내에 있구나!"
```

### (3) Reads vs Centroids QC (`reads_vs_centroids.png`)

transcript 위치(파랑)와 세포 중심점(빨강)을 겹쳐 그림:

```python
ax.scatter(reads['x_location'], reads['y_location'], s=1)        # 파랑: 쪽지
ax.scatter(centroids['x_cell'], centroids['y_cell'], color='red') # 빨강: 세포 중심
```

→ "쪽지들이 세포 중심 주변에 잘 몰려있나?" QC 확인용

---

## [5-6] 확장된 쪽지 저장

```python
reads_original.to_csv(f"{sample_tag}_step5_expanded_transcripts.csv")
```

이제 모든 transcript에 `domain`과 `initial_annotation`(세포 유형)이 붙어있음!

---

## [5-7] Turnover 분석 (핵심!!!)

**이것이 Step 5의 진짜 목적입니다.**

### 핵심 질문:

> "세포 중심에서 얼마나 멀어지면, 그 쪽지가 더 이상 '이 세포의 것'이 아니라 '배경 소음'이 되는가?"

### 비유로 이해하기:

```
 내 집 (= 핵, nucleus)
 내 마당 (= 핵 바깥, cytoplasm)
 이웃집 담장 넘어 (= 배경, background)

집 안에서 발견된 물건 → 100% 내 거!
마당에서 발견된 물건 → 아마 내 거... 거리에 따라 다름
담장 넘어 발견된 물건 → 이웃 거거나 쓰레기!

질문: "마당의 어디까지가 '내 영역'인가?"
답: turnover distance! (이 거리에서 "내 것 같은 정도"와 "쓰레기 같은 정도"가 같아짐)
```

### 상세 과정:

#### Step 7-1: reads 분류 — "세포에 할당됨" vs "미할당(배경)"

```python
# cell_id가 -1이면 → 미할당 (배경)
mask_not_assigned = reads_original['cell_id'] == -1

reads_not_assigned = reads_original[mask_not_assigned]    # 배경 쪽지
reads_assigned_all = reads_original[~mask_not_assigned]   # 세포 소속 쪽지
```

#### Step 7-2: Negative control 제거

```python
# BLANK, NegControlProbe, antisense 같은 가짜 유전자는 제거
_ctrl = reads['feature_name'].str.contains('BLANK|NegControl|antisense')
reads = reads[~_ctrl]
```

이들은 실험에서 일부러 넣은 "가짜 유전자"로, 노이즈를 추가하므로 제거.

#### Step 7-3: 각 transcript → 자기 세포 중심까지 거리 계산

```python
# 각 쪽지의 위치와 소속 세포의 중심(centroid) 사이 거리
reads['distance'] = sqrt(
    (쪽지_x - 세포중심_x)^2 + (쪽지_y - 세포중심_y)^2
)
reads['distance'] = reads['distance'].round(0)  # 1um 단위로 반올림
```

```
예시:
  쪽지 위치: (1240, 5680)
  세포 중심: (1235, 5678)
  거리 = sqrt((5)^2 + (2)^2) = sqrt(29) = 5um

  → 이 쪽지는 세포 중심에서 5um 떨어져 있다
```

#### Step 7-4: overlaps_nucleus (핵 내부 여부) — 없으면 proxy 생성

```python
# "이 쪽지가 핵 안에 있나?" 정보가 필요한데, 없을 수 있음
# → 없으면 "거리의 중앙값보다 가까우면 핵 안"으로 간주

med_dist = reads['distance'].median()  # 예: 4.5um
reads['overlaps_nucleus'] = (reads['distance'] < med_dist).astype(int)
# 거리 < 4.5um → 1 (핵 안), 거리 >= 4.5um → 0 (핵 밖)
```

#### Step 7-5: Background expression profile (배경 유전자 패턴) 만들기

```python
# 미할당 쪽지들을 동네(domain)별로 모아서, 유전자별 개수를 세기
background_express = pd.crosstab(
    reads_not_assigned['domain'],       # 동네
    reads_not_assigned['feature_name']  # 유전자 이름
)
```

```
결과 (배경 유전자 패턴):
              GFAP   MBP   SYT1   AQP4  ...
동네A(피질)     120   80    200    95   ...
동네B(해마)     150   60    300    110  ...
동네C(시상)      90   40    150    70   ...

→ 각 동네의 "바닥에 떨어진 쪽지의 유전자 패턴"
→ 이것이 "배경 소음"의 지문!
```

#### Step 7-6: 메인 루프 — 세포 유형 x 동네별 분석

```python
for celltype in ['Neuron', 'Astrocyte', 'Microglia', ...]:     # 각 세포 유형마다
    for domain in ['동네A', '동네B', '동네C', ...]:             # 각 동네마다
        # 이 세포유형 + 이 동네에 해당하는 쪽지들만 모음
        reads_ctd = reads_assigned[(celltype) & (domain)]
```

단, reads가 `min_reads_per_domain`(기본 5000) 이하인 (세포유형, 동네) 쌍은 건너뜀 — 데이터가 너무 적으면 통계가 불안정하니까.

#### Step 7-7: Nuclear expression profile (핵 유전자 패턴) 만들기

```python
# 핵 안에 있는 쪽지만 모아서 유전자별 개수 세기
reads_nuclear = reads_ctd[reads_ctd['overlaps_nucleus'] == 1]
nuclear_profile = crosstab(reads_nuclear, gene별 개수)
```

```
핵 내부 유전자 패턴 (Neuron, 피질):
  GFAP: 5    ← 뉴런에서는 GFAP 적음
  MBP: 3     ← 뉴런에서는 MBP 적음
  SYT1: 250  ← 뉴런 marker! 핵에서 많이 나옴
  SNAP25: 180
```

#### Step 7-8: Dominant housekeeping gene 제거

```python
# 배경에서 전체 reads의 >15%를 차지하는 유전자 제거
bck_frac = bck_total / bck_sum
dominant_mask = bck_frac > 0.15
# 예: MTRNR2L12가 배경의 52%를 차지 → 제거!
```

```
왜?
  MTRNR2L12 같은 유전자가 어디서나 엄청 많으면
  → 모든 거리에서 상관관계가 0.97로 나옴
  → 차이가 0.03밖에 안 되어서 turnover를 못 찾음!

  이런 유전자를 빼면:
  → 상관관계 범위가 0.03 → 0.20으로 넓어짐
  → turnover를 정확히 찾을 수 있음!
```

#### Step 7-9: 거리별 상관관계 계산 (핵심 중의 핵심!)

```python
# 거리 0um, 1um, 2um, ... 에서 각각
for dist in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, ...]:

    # 이 거리에 있는 쪽지들의 유전자 패턴
    dist_profile = expression_at_distance[dist]

    # 핵 패턴과 비교 (Pearson 상관계수)
    corr_nuclear = corrcoef(dist_profile, nuclear_profile)

    # 배경 패턴과 비교 (Pearson 상관계수)
    corr_background = corrcoef(dist_profile, background_profile)
```

```
결과 테이블:
거리  | 핵 상관  | 배경 상관 | 차이(diff)
──────┼─────────┼──────────┼──────────
 0um  |  0.92   |   0.15   |  0.77    ← 핵에서 가까울수록 "핵 것" 같음
 2um  |  0.88   |   0.20   |  0.68
 4um  |  0.75   |   0.32   |  0.43
 6um  |  0.60   |   0.48   |  0.12
 8um  |  0.45   |   0.55   | -0.10    ← 여기부터 "배경"이 더 비슷!
10um  |  0.35   |   0.62   | -0.27
12um  |  0.28   |   0.68   | -0.40

→ 약 7-8um 부근에서 "뒤바뀜(turnover)"이 일어남!
```

```
그래프로 보면:

    상관관계(PCC)
    1.0 |--\
        |    \  핵 (파랑) — 가까울수록 높음
    0.8 |     \
        |      \
    0.6 |       \\
        |        \\  ← ** 교차점 = Turnover Distance! **
    0.4 |         //
        |       //  배경 (주황) — 멀수록 높음
    0.2 |     //
        |   //
    0.0 |--/
        └──────────────────→ 거리 (um)
           0  2  4  6  8  10 12

    교차점(~7-8um) = 이 거리에서 "내 것"과 "배경"이 반반!
    → 이 이상은 "배경 소음"이므로 확장하면 안 됨!
```

#### Step 7-10: Turnover distance 감지 — Half-life 방식

```python
# 1. 약간 smoothing (노이즈 제거, 윈도우 3)
diff_smooth = summary['diff'].rolling(window=3, center=True).mean()

# 2. diff의 최대값(peak) 찾기
peak_diff = diff_smooth.max()       # 예: 0.77
peak_dist = diff_smooth.idxmax()    # 예: 0um

# 3. Half-life threshold = peak의 10%
halflife_thresh = 0.1 * 0.77 = 0.077

# 4. peak 이후에서 threshold 아래로 떨어지는 첫 지점 = turnover!
after_peak = diff_smooth[거리 >= peak_dist]
turnover = after_peak[after_peak < 0.077].index.min()  # 예: 8um
```

```
비유:
  "핵 신호와 배경 신호의 차이가 최대(0.77)인 지점에서 시작해서,
   그 차이가 원래의 10%(0.077) 이하로 줄어드는 지점을 찾자!
   → 그 지점이 turnover distance!"
```

**왜 "Half-life" 방식이 좋은가?**

```
이전 방식: diff < 0.05 (고정 threshold)
  → 상관관계 범위가 작은 세포유형에서는 영원히 못 찾음
  → 범위가 큰 세포유형에서는 너무 일찍 찾음

Half-life 방식: diff < peak x 0.1 (적응형 threshold)
  → peak가 0.80이면 threshold = 0.08 (큰 신호)
  → peak가 0.05이면 threshold = 0.005 (작은 신호)
  → 자동으로 맞춰짐!
```

#### Step 7-11: 핵 크기 측정 — ConvexHull

```python
def dist_nuc(reads):
    for cell_id, cell_reads in reads.groupby('cell_id'):
        # ConvexHull = 모든 점을 감싸는 가장 작은 볼록 다각형
        hull = ConvexHull(cell_reads[['x_location', 'y_location']])

        # 꼭짓점들에서 중심까지의 평균 거리 = 핵 반경
        nuclei_radius = mean(꼭짓점_거리들)
```

```
비유:
  핵 안의 쪽지들을 고무밴드로 감싸면:

     ,------,
    | .  . . |    . = 핵 안 쪽지
    |  .  .  |    고무밴드 = ConvexHull
    | .   .  |    중심에서 고무밴드까지 평균 거리 = 핵 반경!
     '------'
```

#### Step 7-12: 최적 확장 계산!!!

```python
optimal_expansion = mean_turnover - mean_nuclei_size

# 논문 결과:
# 10.71um (turnover) - 5.06um (핵 반경) = 5.64um (최적 확장!)
```

```
시각적으로:

    <───── 10.71um (turnover distance) ──────>
    <── 5.06um ──><──── 5.64um ────>

    [### 핵 ###]   [..... 확장 .....]  |  배경
    < 핵 반경 >    < 이만큼만         |  (여기는
                     확장하면 됨! >    |   쓰레기!)

    Xenium 기본값은 15um 확장인데,
    실제 최적은 5.64um → 기본값이 너무 과도했다!
```

**음수 보호**: turnover가 핵 반경보다 작으면(비정상) → 0으로 clamping:

```python
if optimal_expansion < 0:
    optimal_expansion = 0.0  # "확장 안 하는 게 낫다"
```

---

## [5-7e] 결과 저장 — 그래프 + CSV + TXT

### (1) Turnover Barplot (`turnover_barplot.png`) → 논문 Fig. 3b

```python
# 세포 유형별 수평 막대 그래프
sns.barplot(y='cluster', x='score', palette=colors)         # 막대: per-domain turnover
sns.scatterplot(y='cluster', x='cell_size', color='black')  # 검정 점: 전체 세포 크기
sns.scatterplot(y='cluster', x='nuclei_size', color='#D83066')  # 핑크 점: 핵 크기
```

```
                     Distance (um)
                     0    5   10   15   20
    Neuron      ==================== *  ^
    Oligo       ============== *  ^
    Astrocyte   ================ *  ^
    Microglia   ========== *  ^

    막대 = turnover distance (per-domain 평균 +/- SD)
    * = 전체 세포 크기 (ConvexHull)
    ^ = 핵 크기 (ConvexHull)
    막대 끝 - ^ = 최적 확장!
```

### (2) Crossover Plots (`crossover_plots/*.png`) → 논문 Fig. 3a 우측

각 (세포유형 x 동네) 쌍별로 상관관계 곡선:

```
    Neuron_동네A.png    Neuron_동네B.png    Astrocyte_동네A.png ...
```

### (3) CSV/TXT 파일들

```
turnover_summary.csv       — 전체 matrix (세포유형 x 동네)
turnover_per_celltype.csv  — 세포유형별 요약 (turnover, 핵크기, 세포크기)
optimal_expansion.txt      — 최종 답: "최적 확장 = X.XX um"
```

---

## 전체 요약: Step 5가 답하는 질문

```
Q: "Cellpose로 핵을 찾았는데, 핵 밖의 쪽지는 어떻게 하지?"
A: "가장 가까운 세포에 배달하되, 세포마다 '적절한 거리'가 있다!"

Q: "적절한 거리는 어떻게 찾아?"
A: "핵에서 멀어질수록 '핵 패턴'과의 닮은 정도가 줄고 '배경 패턴'과의 닮은 정도가 늘어.
    둘이 같아지는 지점(turnover) = 세포의 실질적 경계!
    최적 확장 = turnover - 핵 크기"

Q: "논문의 답은?"
A: "turnover 10.71um - 핵 5.06um = 최적 확장 5.64um
    (Xenium 기본 15um보다 훨씬 짧다!)"
```

---

## 코드 위치 참조 (빠른 찾기용)

| 기능 | 함수/위치 | 라인 |
|------|----------|------|
| ConvexHull 핵 크기 | `dist_nuc()` | 35-53 |
| P2R Domain Mapping | `_load_p2r_spatial_domains()` | 60-177 |
| Turnover 분석 전체 | `calculate_turnover()` | 182-575 |
| kNN 라벨 전달 | `_transfer_celltype_labels()` | 580-668 |
| Step 5 메인 entry | `run_step5()` | 673-1193 |
| 데이터 로딩 (5-1) | `run_step5()` 내부 | 693-778 |
| Domain 할당 (5-2) | 4-level 우선순위 | 804-846 |
| KDTree 배달 (5-3/4) | cKDTree 구축+쿼리 | 932-993 |
| 시각화 (5-5) | scatter + histplot | 997-1033 |
| Turnover 호출 (5-7) | `calculate_turnover()` | 1040-1183 |
