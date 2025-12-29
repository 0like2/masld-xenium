# Pipeline 시각화 가이드 (숨겨진 기능!)

## 📊 Pipeline의 시각화 기능

**좋은 소식**: Pipeline에는 `xb/plotting.py`라는 시각화 모듈이 있습니다!

```
xb/plotting.py
├── map_of_clusters()          ⭐ 공간 클러스터 시각화
├── generate_hex_colors()      색상 생성
├── plot_cell_counts()         품질 히스토그램
└── plot_domains()             도메인 시각화
```

---

## 1️⃣ 공간 클러스터 시각화 (map_of_clusters)

**가장 유용한 함수!**

### 기본 사용법

```python
from xb.plotting import map_of_clusters

# 방법 1: 모든 클러스터를 한 그림에 (권장)
map_of_clusters(
    adata=adata,
    key='leiden_1_4',           # 클러스터링 키
    clusters='all',              # 모두 표시
    size=8,                       # 점 크기
    background='white',           # 배경색
    figuresize=(10, 7),
    save=None                    # './output' (저장 경로)
)

# 방법 2: 클러스터별 개별 그림
map_of_clusters(
    adata=adata,
    key='leiden_1_4',
    clusters='individual',       # 각 클러스터마다 별도 그림
    size=8,
    save='./output'             # 저장 폴더
)

# 방법 3: 특정 클러스터만 강조
map_of_clusters(
    adata=adata,
    key='leiden_1_4',
    clusters=['0', '1', '5'],    # 클러스터 0, 1, 5만 표시
    size=8,
    save='./output'
)
```

### 파라미터 설명

| 파라미터 | 타입 | 설명 | 기본값 |
|---------|------|------|--------|
| `adata` | AnnData | 전처리된 데이터 | 필수 |
| `key` | str | `adata.obs`의 클러스터링 키 | 'leiden' |
| `clusters` | str or list | 'all' / 'individual' / ['0','1'] | 'all' |
| `size` | int | 점의 크기 (픽셀) | 8 |
| `background` | str | 배경색 ('white', 'black' 등) | 'white' |
| `figuresize` | tuple | 그림 크기 (가로, 세로) | (10, 7) |
| `save` | str | 저장 폴더 경로 (None=표시만) | None |
| `format` | str | 저장 형식 ('pdf', 'png', 'jpg') | 'pdf' |

### 출력 파일명

```
저장 옵션이 활성화되면:
├─ map_all_clusters_{size}_{background}_{key}.pdf       # clusters='all'
├─ map_individual_cluster_0_{size}{background}_{key}.pdf # 클러스터 0
├─ map_individual_cluster_1_{size}{background}_{key}.pdf # 클러스터 1
└─ map_group_of_clusters_{012...}_{size}{background}_{key}.pdf  # 선택 클러스터들
```

---

## 2️⃣ 품질 히스토그램 (plot_cell_counts)

### 기본 사용법

```python
from xb.plotting import plot_cell_counts

clustering_params = {
    'min_counts_x_cell': 10,      # 최소 카운트 기준선
    'min_genes_x_cell': 3         # 최소 유전자 기준선
}

plot_cell_counts(
    adata=adata,
    plot_path='./output/',
    save=True,
    clustering_params=clustering_params
)
```

### 출력

```
2개 패널의 히스토그램:
├─ 왼쪽: 세포당 카운트 분포
│   └─ 빨간 선: min_counts 기준값
└─ 오른쪽: 세포당 유전자 수 분포
    └─ 빨간 선: min_genes 기준값

저장 파일: cell_counts_histogram.png
```

### 시각적 예시

```
세포당 카운트                    세포당 유전자 수
|                                |
| █                              | █
| █ █                            | █ █
| █ █ █                          | █ █ █
|_█_█_█_|_ (min=10)             |_█_█_█_|_ (min=3)
  0  100  200                      0  100  200
```

---

## 3️⃣ 도메인 시각화 (plot_domains)

### 기본 사용법

```python
from xb.plotting import plot_domains

plot_domains(
    adata=adata,
    groupby='leiden_1_4'  # 도메인/클러스터 키
)
```

### 특징

- Scanpy의 spatial 플롯 사용
- 여러 샘플이 있으면 자동으로 각각 표시
- 로컬 공간 구조를 시각화

---

## 4️⃣ 색상 생성 (generate_hex_colors)

### 기본 사용법

```python
from xb.plotting import generate_hex_colors

# 70개의 랜덤 색상 생성
colors = generate_hex_colors(num_colors=70)
# 출력: ['#a3f2d1', '#f9a2e1', '#2c9e3a', ...]

# 사용자 정의 개수
colors = generate_hex_colors(num_colors=30)
```

---

## 🎯 Step별 시각화 권장

### Step 1 후 권장 시각화

```python
from xenium_pipeline_main import XeniumPipeline
from xb.plotting import map_of_clusters, plot_cell_counts

# Step 1 실행
pipeline = XeniumPipeline(
    xenium_input_path="./data/unprocessed_adata/ms_brain_rep1.h5ad",
    output_path="./output",
    sample_tag="ms_brain_rep1",
    load_preprocessed=True
)

adata = pipeline.step1_preprocess(
    target_sum=100,
    scale=False,
    hvg=False,
    save=True
)

# 1️⃣ 품질 확인
plot_cell_counts(
    adata=adata,
    plot_path='./output/',
    clustering_params={
        'min_counts_x_cell': 10,
        'min_genes_x_cell': 3
    }
)

# 2️⃣ 클러스터 공간 분포 (전체)
map_of_clusters(
    adata=adata,
    key='leiden_1_4',
    clusters='all',
    size=6,
    background='white',
    save='./output',
    format='png'
)

# 3️⃣ 클러스터별 개별 시각화 (선택)
map_of_clusters(
    adata=adata,
    key='leiden_1_4',
    clusters='individual',
    size=8,
    save='./output'
)

# 4️⃣ 특정 클러스터 강조 (선택)
map_of_clusters(
    adata=adata,
    key='leiden_1_4',
    clusters=['0', '1', '2', '3'],  # 상위 4개 클러스터
    size=10,
    save='./output'
)
```

---

## 🚀 Pipeline + Squidpy 조합 사용법

**최고의 시각화 조합!**

```python
import scanpy as sc
import squidpy as sq
from xb.plotting import map_of_clusters, plot_cell_counts
from xenium_pipeline_main import XeniumPipeline

# 1단계: Pipeline으로 최적 전처리
pipeline = XeniumPipeline(
    xenium_input_path="./data/ms_brain_rep1.h5ad",
    output_path="./output",
    sample_tag="ms_brain",
    load_preprocessed=True
)

adata = pipeline.step1_preprocess(
    target_sum=100,
    scale=False,
    hvg=False
)

# 2단계: Pipeline 시각화 (기본)
plot_cell_counts(adata=adata, plot_path='./output/')
map_of_clusters(adata=adata, key='leiden_1_4', save='./output')

# 3단계: Squidpy 추가 분석 (고급)
# UMAP 생성
sc.tl.umap(adata)

# 마커 검색
sc.tl.rank_genes_groups(adata, groupby='leiden_1_4', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=10)

# 공간 분석
sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True)
sq.gr.nhood_enrichment(adata, cluster_key='leiden_1_4')
sq.pl.nhood_enrichment(adata, cluster_key='leiden_1_4')

# 4단계: 마커 기반 세포 유형 주석
marker_dict = {
    "Brain_Neuron": ["SYN1", "SNAP25", "NRGN"],
    "Brain_Glia": ["GFAP", "AQP4", "S100B"],
    "Brain_Astro": ["GLAST", "GLUTAMATE", "GFAP"],
}

for label, genes in marker_dict.items():
    genes_present = [g for g in genes if g in adata.var_names]
    if genes_present:
        sc.tl.score_genes(adata, gene_list=genes_present, score_name=f"score_{label}")

# 5단계: 결과 저장
adata.write('./output/ms_brain_annotated.h5ad')
```

---

## 📈 시각화 결과 예시

### map_of_clusters (clusters='all')

```
공간 분포도:
  0  1  2  3  4
5 🔵🟢🔵🟢🔵
4 🟢🔵🟢🔵🟢
3 🔵🟢🔵🟢🔵
2 🟢🔵🟢🔵🟢
1 🔵🟢🔵🟢🔵
  0  1  2  3  4

범례: 🔵=클러스터0, 🟢=클러스터1, 등...
```

### map_of_clusters (clusters='individual')

```
클러스터 0:          클러스터 1:
  강조 영역          강조 영역
  🔵🔵🔵             🟢🟢🟢
  🔵🔵🔵             🟢🟢🟢
  ⚪⚪⚪ (배경)      ⚪⚪⚪ (배경)
```

### plot_cell_counts

```
세포당 카운트 분포        세포당 유전자 분포
높이 |                  높이 |
     | ___              | ___
     ||   |___          |     |___
     ||   |   |___      ||   |   |___
     ||___|___|___|     ||___|___|___|
     0  100 200 300 400 0  50 100 150 200
```

---

## 🔧 커스터마이징 팁

### 1. 색상 커스터마이징

```python
# 자동 색상 생성 대신 사용자 정의
import matplotlib.colors as mcolors

# 특정 색상 사용
custom_colors = {
    0: '#FF0000',  # 빨강
    1: '#00FF00',  # 초록
    2: '#0000FF',  # 파랑
}

# adata.uns에 저장
adata.uns['leiden_1_4_colors'] = [custom_colors.get(i, '#CCCCCC')
                                  for i in range(max(adata.obs['leiden_1_4'])+1)]
```

### 2. 고해상도 저장

```python
# format 파라미터로 조절
map_of_clusters(
    adata=adata,
    key='leiden_1_4',
    save='./output',
    format='png',  # 더 높은 품질
    figuresize=(20, 15)  # 더 큰 크기
)

# 또는 수동으로 DPI 조절
import matplotlib.pyplot as plt
plt.savefig('./output/custom.png', dpi=300)  # 300 DPI
```

### 3. 특정 영역만 확대

```python
# map_of_clusters의 문제: 전체 좌표 사용
# 해결: 부분 데이터셋으로 별도 분석

# 상단 영역만
adata_top = adata[adata.obs['y_centroid'] > adata.obs['y_centroid'].median()]
map_of_clusters(adata=adata_top, key='leiden_1_4', save='./output',
                format='png', figuresize=(15, 10))
```

---

## ❌ Pipeline 시각화의 한계

| 기능 | Pipeline | Squidpy |
|-----|---------|---------|
| 공간 클러스터 | ✓ | ✓ (더 예쁨) |
| UMAP/t-SNE | ✗ | ✓ |
| 마커 검색 | ✗ | ✓ |
| 공간 이웃 분석 | ✗ | ✓ |
| Ligand-Receptor | ✗ | ✓ |
| 공간 자기상관 | ✗ | ✓ |
| 상호작용 대시보드 | ✗ | ✓ (scanviz) |

---

## 💡 결론

### Pipeline의 시각화 강점
✓ **공간 클러스터 시각화** (기본적이지만 빠름)
✓ **품질 확인 히스토그램** (QC용)
✓ **간단하고 직관적** (추가 설정 불필요)

### 완벽한 조합
```
Pipeline: 최적 전처리 + 기본 시각화
         ↓
Squidpy: 심화 분석 + 고급 시각화
         ↓
최종 결과: 정확성 + 완성도 모두 확보! ✓
```

---

## 📁 전체 코드 예시

**complete_analysis.py**:
```python
#!/usr/bin/env python
"""Pipeline + Squidpy 완전 분석 스크립트"""

from xenium_pipeline_main import XeniumPipeline
from xb.plotting import map_of_clusters, plot_cell_counts
import scanpy as sc
import squidpy as sq

# 1단계: Pipeline 전처리
print("Step 1: 전처리 실행...")
pipeline = XeniumPipeline(
    xenium_input_path="./data/ms_brain_rep1.h5ad",
    output_path="./output",
    sample_tag="ms_brain_rep1",
    load_preprocessed=True
)

adata = pipeline.step1_preprocess(target_sum=100, scale=False, hvg=False)

# 2단계: Pipeline 시각화
print("Step 2: Pipeline 시각화...")
plot_cell_counts(adata=adata, plot_path='./output/')
map_of_clusters(adata=adata, key='leiden_1_4', clusters='all', save='./output')

# 3단계: Squidpy 추가
print("Step 3: Squidpy 분석...")
sc.tl.umap(adata)
sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True)
sq.gr.nhood_enrichment(adata, cluster_key='leiden_1_4')

# 4단계: 저장
adata.write('./output/final_analysis.h5ad')
print("✓ 완료! ./output 폴더에 결과가 저장되었습니다.")
```

실행:
```bash
python complete_analysis.py
```

---

**참고 파일**:
- `/data1/project/20rak/masld_xenium/Xenium_benchmarking/xb/plotting.py`
- `/data1/project/20rak/masld_xenium/Xenium_benchmarking/xenium_pipeline_main.py`

**작성일**: 2025-12-29
