# Step 3: Resegmentation (Cellpose) - Final Review

## 1. Overview

### 1.1 What This Step Does
Step 3은 **Cellpose** 딥러닝 모델을 사용하여 DAPI 이미지에서 세포 핵을 재분할(resegmentation)하고, 핵 마스크를 확장(expansion)하여 transcript를 새로운 세포에 재할당한다. Xenium의 기본 segmentation 대신 독립적인 segmentation을 수행하여, 이후 Step 4에서 품질 비교의 기준이 된다.

### 1.2 Paper Context
논문의 "Baysor and Cellpose outperform standard Xenium segmentation" 섹션 (p.817):
- "The Cellpose (v2.2.3) deep-learning models 'nuclei' (CPn) and 'cyto' (CPc) were applied on the DAPI channels with diameter parameters of none, 20, 30 and 40"
- Cellpose cytoplasmic model이 segmentation benchmark에서 주요 비교 대상
- Fig. 2a에서 workflow 설명: Nuclei Segmentation → Cellpose (cyto mod.) → Comparison

### 1.3 논문 Figure 매칭
| 논문 Figure | 설명 | Step 3 구현 |
|------------|------|------------|
| **Fig. 2a** | 비교 workflow (Nuclei Seg → Cellpose → Comparison) | `run_step3()` 전체 |
| **Fig. 3c** | 6가지 segmentation ROI 비교 (Cellpose 포함) | Step 3 masks → Step 6에서 비교 |
| **Fig. 3d** | ARI heatmap (segmentation 전략 간) | Step 3 결과 → Step 6에서 계산 |
| **Extended Data Fig. 4g** | Resegmented datasets ROI | Step 3 masks overlay |
| **Extended Data Fig. 5a-b** | Segmentation ROI 비교 (다양한 방법) | Step 3 masks → Step 6 |
| **Extended Data Fig. 5c** | Segmentation metrics heatmap | Step 3 → Step 6 metrics |

---

## 2. Pipeline Process Flow

```
Step 0/1/2 Outputs
├── DAPI Image (morphology_focus.ome.tif)
├── transcripts.csv (원본 transcript 좌표)
├── Step 2 domain_polygons.json (선택)
         ↓
    [Step 3: Resegmentation]
         ↓
    ┌───── 3-1. Device Detection ────────────────────────────┐
    │  CUDA GPU → MPS (Apple) → CPU 자동 선택               │
    └────────────────────────────────────────────────────────┘
         ↓
    ┌───── 3-2. DAPI Loading ────────────────────────────────┐
    │  Multi-channel TIF → single channel DAPI 추출          │
    │  다운샘플링 (overview용)                                │
    └────────────────────────────────────────────────────────┘
         ↓
    ┌───── 3-3. Cellpose Segmentation ──────────────────────┐
    │  model_type='nuclei', diameter=auto                    │
    │  Tiled processing (대용량 이미지 지원)                   │
    │  → nuclei masks (integer label image)                  │
    └────────────────────────────────────────────────────────┘
         ↓
    ┌───── 3-4. Mask Expansion ─────────────────────────────┐
    │  expand_labels(distance=400px)                         │
    │  nuclei → expanded cell masks                          │
    └────────────────────────────────────────────────────────┘
         ↓
    ┌───── 3-5. Transcript Assignment ──────────────────────┐
    │  각 transcript (x,y) → expanded mask pixel 조회        │
    │  → cell_id 할당                                        │
    └────────────────────────────────────────────────────────┘
         ↓
    ┌───── 3-6. AnnData Construction ──────────────────────┐
    │  Cell×Gene count matrix 구축 (nuclear labels 사용)     │
    │  Cell centroids + area/perimeter (regionprops)         │
    └────────────────────────────────────────────────────────┘
         ↓
    ┌───── 3-6b. Edge Artifact Flagging ───────────────────┐
    │  이미지 경계 접촉 세포 감지 → is_edge_cell 플래그       │
    │  (QC 필터링 전에 수행하여 인덱스 일치 보장)             │
    └────────────────────────────────────────────────────────┘
         ↓
    ┌───── 3-6c. QC Filtering ─────────────────────────────┐
    │  min_counts, min_genes 기반 저품질 세포 제거           │
    └────────────────────────────────────────────────────────┘
         ↓
    ┌───── 3-7. Domain Assignment (선택) ───────────────────┐
    │  Cell centroid → domain polygon (Point-in-Polygon)     │
    └────────────────────────────────────────────────────────┘
         ↓
├── {sample_tag}_step3_resegmented_masks.tif    ← 확장된 cell masks
├── {sample_tag}_step3_nuclei_masks.tif         ← 원본 nuclei masks
├── {sample_tag}_step3_resegmented.h5ad         ← 재분할 AnnData
├── {sample_tag}_step3_transcripts_resegmented.csv ← transcript 재할당
└── 시각화 (overlay, QC, ROI)
```

---

## 3. Sub-step 상세 분석

### 3.1 Device Detection (`_detect_device()`, lines 120-141)

```python
# 우선순위: CUDA GPU > MPS (Apple Silicon) > CPU
if torch.cuda.is_available():
    device = torch.device('cuda')
    use_gpu = True
elif hasattr(torch.backends, 'mps') and torch.backends.mps.is_available():
    device = torch.device('mps')
    use_gpu = True
else:
    device = torch.device('cpu')
    use_gpu = False
```

**GPU 성능 영향**:
- CUDA GPU: ~10분 (전체 Xenium DAPI)
- CPU: ~2-6시간 (이미지 크기에 따라)

### 3.2 DAPI Image Loading (`_load_dapi_downsampled()`, lines 145-162)

- 다중 채널 TIF에서 DAPI 채널 추출 (보통 채널 0)
- 메모리 효율적 로딩 (tifffile 사용)
- Overview용 다운샘플링

### 3.3 Cellpose Tiled Segmentation (`_run_cellpose_tiled()`, lines 45-117)

**왜 Tiled Processing이 필요한가**:
Xenium DAPI 이미지는 보통 20,000 × 20,000+ 픽셀로, 한 번에 GPU 메모리에 로드할 수 없다. Tiled 처리로 이미지를 분할하여 순차 처리한다.

```python
# 타일 크기 자동 결정 (_auto_tile_size(), lines 31-42)
# GPU VRAM의 40%를 flow dynamics tensor용으로 예산, 4096-16384px 범위
tile_size = _auto_tile_size(image_shape, gpu_memory)

# 각 타일에 대해 Cellpose 실행 (512px overlap으로 경계 세포 보존)
for tile in tiles:
    masks_tile, flows, styles = model.eval(tile, diameter=diameter)
    # core region만 사용하여 stitching, label offset으로 고유 ID 보장
```

**핵심 파라미터**:
| 파라미터 | 기본값 | 의미 |
|----------|--------|------|
| `pretrained_model` | 'nuclei' | Cellpose 모델: nuclei(핵만) / cyto(세포체). (v4 API: `model_type` → `pretrained_model`) |
| `diameter` | None (auto) | 예상 핵 직경 (pixels). None=자동 감지 |
| `tile_overlap` | 512 | 타일 간 겹침 영역 (경계 세포 보존) |
| `batch_size` | auto | GPU 메모리에 따라 자동 결정: 32 (≥30GB), 16 (≥16GB), 8 (default) |

**참고**: Cellpose v4에서 API가 변경됨 — `model_type` 대신 `pretrained_model` 사용, `channels` 파라미터 불필요, `eval()` 반환값이 4개→3개로 변경.

### 3.4 Mask Expansion

**개념**: Cellpose는 핵만 탐지하므로, cytoplasmic transcript를 포함하려면 마스크를 확장해야 한다.

```python
from skimage.segmentation import expand_labels
expanded_masks = expand_labels(nuclei_masks, distance=expansion_distance)
# expansion_distance=400 pixels ≈ 85 µm (at 4.70588 px/µm)
```

**논문 맥락**: "Xenium's nuclear segmentation is followed by a default radius expansion of 15 µm" (p.817). 파이프라인의 400 pixel ≈ 85 µm은 논문의 15 µm보다 훨씬 크지만, 이는 설정 가능하다.

**설정값**:
```yaml
resegmentation:
  expansion_distance: 400    # pixels (0이면 확장 없이 nuclei만)
```

### 3.5 Transcript-to-Cell Mapping (lines 478-596)

```python
# 좌표 컬럼 자동 감지: x_global_px/y_global_px → global_x/global_y → x_location/y_location
# x_location/y_location인 경우 µm → pixel 변환 필요
pixel_x = int(transcript_x * um_per_pixel_inv)
pixel_y = int(transcript_y * um_per_pixel_inv)

# CosMx 감지: technology=='cosmx'이면 x/y 스왑 (lines 495-498)

# 확장된 mask에서 cell_id 조회 (이미지 좌표계: masks[y, x])
cell_labels = masks[y_coords, x_coords]
# cell_labels = 0이면 미할당 (어떤 세포에도 속하지 않음)

# Dual Assignment (expansion 활성화 시):
# in_cell = nuclei_masks[y, x]       → nuclear label
# closest_cell = expanded_masks[y, x] → expanded label
```

### 3.6 AnnData Construction

재할당된 transcript로부터 새로운 cell×gene matrix 구축:

```python
# 0. Control probe 제거 (NegControl, BLANK, antisense 등)
df_assigned = df_assigned[~df_assigned[gene_col].str.contains(
    r'NegControl|BLANK|antisense', case=False, na=False)]

# 1. Cell×Gene count matrix (nuclear labels 기반 - expansion 활성화 시)
if nuclei_masks is not None:
    cell_gene_matrix = pd.crosstab(df_nuclear['in_cell'], df_nuclear[gene_col])
else:
    cell_gene_matrix = pd.crosstab(df_for_matrix['cell_id_reseg'], df_for_matrix[gene_col])

# 2. Cell centroids (from mask regionprops)
from skimage.measure import regionprops
for region in regionprops(expanded_masks):
    centroids[region.label] = (region.centroid[1], region.centroid[0])  # (y, x)

# 3. AnnData 생성 + morphology
adata = sc.AnnData(cell_gene_matrix)
adata.obs['y_centroid'], adata.obs['x_centroid'] = centroids
adata.obs['cell_area_um2'] = cell_area_px * (pixel_to_um ** 2)
adata.obs['cell_perimeter_um'] = cell_perimeter_px * pixel_to_um

# 4. Edge artifact flagging (QC 전에 수행 → 인덱스 일치 보장)
if nuclei_masks is not None and expansion_distance > 0:
    edge_labels = set(masks on all 4 image borders)
    adata.obs['is_edge_cell'] = [idx in edge_labels for idx in cell_indices]

# 5. QC filtering (edge flagging 후)
sc.pp.filter_cells(adata, min_counts=mincounts)
sc.pp.filter_cells(adata, min_genes=mingenes)
```

**핵심 순서**: Edge flagging이 QC filtering보다 **먼저** 수행된다. `cell_indices`가 AnnData 생성 시점의 전체 세포 목록이므로, QC filtering 전에 `is_edge_cell`을 할당해야 인덱스가 일치한다. QC filtering이 `filter_cells()`로 행을 제거하면 `is_edge_cell`도 자동으로 서브셋된다.

**Dual Assignment**: expansion 활성화 시, `in_cell` (nuclear label)과 `closest_cell` (expanded label)이 모두 기록된다. AnnData의 count matrix는 nuclear label 기반으로 구축되어, 핵 내부 transcript만 세포별 발현량에 반영된다.

### 3.7 Domain Assignment (선택)

Step 2에서 생성된 `domain_polygons.json`이 있으면, 각 cell centroid를 domain에 할당:

```python
from shapely.geometry import Point
for cell_id, (cx, cy) in centroids.items():
    for domain_id, polygon in domain_polygons.items():
        if polygon.contains(Point(cx, cy)):
            adata.obs.loc[cell_id, 'domain'] = domain_id
            break
```

---

## 4. Notebook vs Pipeline 구현 비교

### 4.1 원본 Notebook들
| 노트북 | Pipeline 함수 | 상태 |
|--------|--------------|------|
| `3_1_resegmentation_notebooks/batch_segmentation_cellpose-Xenium.ipynb` | `_run_cellpose_tiled()` + `run_step3()` | 구현됨 |
| `3_2_cell_to_domain_assignment_notebooks/Xenium Dataset_domain_assignment-newsegmentation.ipynb` | Domain assignment in `run_step3()` | 구현됨 |

### 4.2 상세 비교

| 기능 | Notebook | Pipeline | 비고 |
|------|----------|----------|------|
| Cellpose 실행 | `model.eval()` 직접 호출 | Tiled processing 지원 | 개선: 대용량 이미지 |
| Mask 확장 | `expand_labels()` | `expand_labels()` | 동일 |
| Transcript 할당 | pixel lookup | pixel lookup | 동일 |
| Count matrix | pandas groupby | pandas groupby | 동일 |
| Centroids | `regionprops` | `regionprops` | 동일 |
| Domain 할당 | shapely Point-in-Polygon | shapely Point-in-Polygon | 동일 |
| 디바이스 감지 | 수동 설정 | 자동 (CUDA>MPS>CPU) | 개선 |
| QC | 하드코딩 | config 기반 | 개선 |
| 시각화 | 별도 | 통합 QC/overlay 플롯 | 개선 |

### 4.3 주요 차이점 및 개선사항

1. **expand_labels 사용**: 파이프라인은 `expand_labels()`을 사용하여 핵 마스크를 균일하게 확장. 노트북에서도 동일 방식 사용.

2. **CosMx 좌표 스왑**: 노트북의 CosMx 버전에서 x/y 좌표를 스왑하는 로직이 파이프라인에도 포함 (`technology=='cosmx'` 검출). Xenium에서는 작동하지 않음.

3. **Tiled processing**: 노트북에 없는 기능으로, 파이프라인에서 추가하여 대용량 이미지 처리를 가능하게 함. `_auto_tile_size()`가 GPU VRAM 기반으로 타일 크기 자동 결정.

4. **Cellpose v4 호환**: `model_type` → `pretrained_model` API 변경, `eval()` 반환값 4→3개 변경 대응.

5. **Dual Assignment**: 노트북에 없는 기능. `in_cell` (nuclear label)과 `closest_cell` (expanded label) 동시 기록으로, AnnData는 nuclear label 기반 matrix + expanded label 기반 할당 정보를 모두 보존.

6. **Edge Artifact Flagging**: 노트북에 없는 기능. 이미지 경계 접촉 세포를 `is_edge_cell` 플래그로 표시. QC 전에 수행하여 인덱스 일치 보장 (이전 버그 수정).

7. **Control Probe Filtering**: AnnData matrix 구축 전에 `NegControl|BLANK|antisense` 패턴의 control probe를 제거.

---

## 5. 시각화 목록 및 분석법

| # | 파일명 | 유형 | 패널 수 | 설명 | 분석법 |
|---|-------|------|---------|------|--------|
| 1 | `{tag}_step3_dapi_overview.png` | Image | 1 | DAPI 전체 이미지 (vmax=99.5 percentile) | 이미지 품질 확인 |
| 2 | `{tag}_step3_dapi_mask_overlay.png` | Overlay | 2-3 | DAPI + expanded mask (cyan) + nuclei mask (orange) | Segmentation 정확도 시각적 확인 |
| 3 | `{tag}_step3_dapi_zoomed_roi.png` | Zoomed | 3 | ROI 확대: DAPI / 경계선 / colored mask blend | 세포 경계 정밀도 확인 |
| 4 | `{tag}_step3_mask_qc.png` | Label2RGB | 1-2 | Mask 색상화 (expanded full + nuclei zoomed) | 비정상적 크기 세포 탐지 |
| 5 | `{tag}_step3_mask_transcript_overlay.png` | Overlay | 1 | Binary mask (α=0.5) + transcript scatter (red, 50K sample) | Transcript 할당 정확도 시각 검증 |
| 6 | `{tag}_step3_dapi_transcript_overlay.png` | Overlay | 1 | DAPI + transcript scatter (red, 50K sample) | Transcript 공간 분포 확인 |
| 7 | `{tag}_step3_domain_polygon_overlay.png` | Overlay | 1 | Domain polygon + cell centroids (선택, domain_map_path 필요) | Domain 할당 확인 |

**분석법 상세**:
- `dapi_mask_overlay.png`: Cyan contour (expanded)와 Orange contour (nuclei)가 DAPI 밝은 영역과 정확히 일치하는지 확인. 불일치 = segmentation 오류
- `mask_qc.png`: label2rgb로 색상화된 mask. 극단적으로 큰 마스크 = over-segmentation 실패, 극단적으로 작은 마스크 = debris
- `mask_transcript_overlay.png` / `dapi_transcript_overlay.png`: 최근 추가된 시각화. Mask 경계 내에 transcript가 적절히 위치하는지 확인

---

## 6. Input / Output 상세

### 6.1 Input
| 파일 | 형식 | 설명 |
|------|------|------|
| `morphology_focus.ome.tif` | TIFF (multi-channel) | DAPI 형광 이미지 |
| `transcripts.csv` | CSV | 원본 transcript 좌표 (x, y, gene, cell_id) |
| `domain_polygons.json` (선택) | GeoJSON | Step 2 도메인 경계 |
| Step 1/2 `*.h5ad` | AnnData | 이전 단계 결과 (셀 메타데이터) |

### 6.2 Output
| 파일 | 형식 | 설명 |
|------|------|------|
| `{tag}_step3_resegmented_masks.tif` | TIFF (uint32) | 확장된 cell mask (pixel=cell_id) |
| `{tag}_step3_nuclei_masks.tif` | TIFF (uint32) | 원본 nuclei mask (확장 전) |
| `{tag}_step3_resegmented.h5ad` | AnnData | 재분할 cell×gene 데이터 |
| `{tag}_step3_transcripts_resegmented.csv` | CSV | Transcript 재할당 결과 (cell_id_reseg, in_cell, closest_cell, distance_to_centroid 포함) |
| 7개 시각화 파일 | PNG | QC/overlay 플롯 (위 시각화 목록 참조) |

### 6.3 AnnData 구조
```
adata.X       → (n_cells × n_genes) count matrix [nuclear labels 기반]
adata.obs     → cell_id, x_centroid, y_centroid, n_counts, n_genes,
                cell_area_px, cell_area_um2, cell_perimeter_px, cell_perimeter_um,
                is_edge_cell (bool), region_annotation (선택, domain 할당 시)
adata.var     → gene_name (features)
adata.uns     → segmentation_method='cellpose', expansion_distance=400
```

### 6.4 Transcripts CSV 구조
```
cell_id_reseg  → 확장된 mask의 cell ID (0=미할당)
in_cell        → nuclei mask의 cell ID (nuclear label)
closest_cell   → expanded mask의 cell ID (= cell_id_reseg)
distance_to_centroid → transcript에서 가장 가까운 세포 centroid까지 거리
region_annotation    → domain 할당 (선택)
```

---

## 7. 관련 설정값 정리

```yaml
resegmentation:
  run_resegmentation: true           # Step 3 실행 여부
  expansion_distance: 400            # Mask 확장 거리 (pixels). 0=확장 없이 nuclei만
  um_per_pixel_inv: 4.70588          # Pixel/µm 변환 (Xenium default)

  cellpose:
    # Cellpose v4 API: model_type → pretrained_model (코드에서 자동 처리)
    model_type: 'nuclei'             # 'nuclei' 또는 'cyto' (config에서는 model_type으로 설정)
    diameter: null                   # null=자동 감지
    channels: [0, 0]                 # Grayscale DAPI (v4에서는 불필요하지만 호환성 유지)
    # batch_size: GPU 메모리에 따라 자동 (32/16/8)
    # tile_size: VRAM 40% 기반 자동 계산 (4096-16384px)
    # tile_overlap: 512px (하드코딩)

  domain_assignment:
    n_domains: 5                     # 예상 도메인 수 (사용되지 않음)
    domain_map_path: null            # Domain polygon JSON 경로

formatting:
  mincounts: 10                      # QC: 최소 transcript/cell
  mingenes: 3                        # QC: 최소 genes/cell

comparison:
  technology: "xenium"               # CosMx x/y 스왑 감지에 사용
```

### 7.1 파라미터 조정 가이드

| 파라미터 | 현재값 | 범위/옵션 | 영향도 | 조정 가이드 |
|----------|--------|----------|--------|------------|
| `run_resegmentation` | `true` | true/false | - | Step 3 실행 여부. false면 Step 4-6에서 원본 segmentation만 사용 |
| `expansion_distance` | `400` | 0-1000 px | **TOP 1** | **가장 중요한 파라미터**. 핵 마스크 확장 거리 (pixels). 0=확장 없음(nuclei only). 400px ≈ 85µm (매우 공격적). 논문 최적: ~27px (5.64µm). **권장: 50-100px (10-21µm)** |
| `um_per_pixel_inv` | `4.70588` | 기술 의존 | LOW | Pixel/µm 변환 팩터. Xenium=4.70588 (0.2125 µm/pixel). 변경 불필요 (장비 고정값) |
| `cellpose.model_type` | `'nuclei'` | nuclei/cyto | **HIGH** | Cellpose 모델. `nuclei`=핵만 탐지(확장 필수), `cyto`=세포체 포함(확장 불필요하나 DAPI only에서 정확도 낮음). 논문: 두 모델 비교 |
| `cellpose.diameter` | `null` | null/10-80 px | **HIGH** | 예상 핵 직경. `null`=자동 감지(권장). 수동 설정 시 생물학적 핵 크기에 맞춤: 뇌 ~30-40px |
| `cellpose.channels` | `[0, 0]` | [0,0] | LOW | Grayscale DAPI. Cellpose v4에서는 불필요 (호환성 유지) |
| `domain_assignment.domain_map_path` | `null` | 파일 경로/null | LOW | Domain polygon JSON. null=도메인 할당 건너뜀. Step 2에서 생성한 polygon이 있을 때만 사용 |
| `formatting.mincounts` | `10` | 1-50 | MEDIUM | QC 필터링에 사용 (Step 0과 공유). 재분할된 세포에도 동일 기준 적용 |
| `formatting.mingenes` | `3` | 1-20 | MEDIUM | QC 필터링에 사용 (Step 0과 공유) |
| `comparison.technology` | `"xenium"` | xenium/cosmx/vizgen 등 | LOW | 좌표 변환 및 CosMx x/y 스왑 감지에 사용. Xenium 데이터는 변경 불필요 |

**튜닝 팁**:
- `expansion_distance`가 **파이프라인 전체에서 가장 중요한 파라미터**. Step 5의 turnover 분석 결과(최적 ~5.64µm = ~27px)를 참고하여 설정
- 현재 기본값 400px(≈85µm)는 논문의 최적값(5.64µm)보다 15배 큼 → Step 5 결과를 먼저 확인 후 조정 권장
- `cellpose.diameter=null` (자동 감지)이 대부분의 경우 잘 작동. 비정상적 segmentation 시 수동 설정 시도
- `model_type='cyto'`는 DAPI만 있는 Xenium에서는 핵/세포체 구분이 어려워 `nuclei`보다 불리할 수 있음

---

## 8. 핵심 로직의 의미

### 8.1 expansion_distance = 400 pixels
- 400 pixels ÷ 4.70588 px/µm ≈ **85 µm**
- Xenium 기본 확장: **15 µm**
- 85 µm은 매우 공격적인 확장 → 세포 간 transcript 오할당(misassignment) 위험
- 논문에서는 다양한 확장 거리 (1, 2, 5, 10, 15 µm)를 비교 (Extended Data Fig. 5c)
- **권장**: 실제 분석에서는 config로 적절한 값 설정 필요

### 8.2 Cellpose 'nuclei' vs 'cyto' model
- **nuclei**: 핵만 탐지 → 확장이 필수
- **cyto**: 세포체(cytoplasm)까지 탐지 → 확장 불필요하지만 DAPI만으로는 정확도 떨어짐
- 논문에서는 두 모델을 모두 평가 (Fig. 3c, Extended Data Fig. 5c)

### 8.3 Tiled Processing의 tile_overlap
- 겹침 영역이 없으면 타일 경계의 세포가 분리됨
- 512 pixel 겹침으로 경계 세포를 양쪽 타일에서 탐지 후 병합

---

## 10. 시각화 상세 분석 가이드

### 10.1 DAPI Overview (`dapi_overview.png`)

**무엇을 봐야 하는가**:
```
    ┌───────────────────────────────────────┐
    │ ████████████████████████████████████  │
    │ ██ 밝은 핵들 (DAPI 형광) ████████████  │
    │ ████████████████████████████████████  │
    │ ██                             █████  │
    │ ████  어두운 영역 (세포질/배경) ██████  │
    │ ████████████████████████████████████  │
    └───────────────────────────────────────┘

    확인 사항:
    ① 이미지 품질: 밝기/대비가 적절한지
    ② 조직 구조: 뇌의 gray/white matter 경계가 보이는지
    ③ 아티팩트: 기포, 접힘, 파손 등이 없는지
    ④ 핵 밀도: 영역별 핵 밀도 차이 (cortex > white matter)
```

**해석법**:
- **좋은 DAPI 이미지**: 개별 핵이 명확히 구분, 균일한 배경, 조직 구조 보존
- **나쁜 DAPI 이미지**: 흐릿한 핵, 불균일 배경, 과도한 자가형광
- 이 이미지가 Cellpose의 입력이므로, 이미지 품질이 segmentation 품질을 결정

---

### 10.2 DAPI + Mask Overlay (`dapi_mask_overlay.png`) → 논문 Extended Data Fig. 4g 관련

**논문 위치**: Extended Data Fig. 4g (p.829) - "Resegmented datasets ROI"

**무엇을 봐야 하는가**:
```
    ┌─────────────────────────────────┐
    │     ╭cyan╮    ╭cyan╮           │
    │     │ ●● │    │ ●● │ ← mask   │  DAPI (회색) + Cellpose mask (cyan)
    │     ╰────╯    ╰────╯           │
    │                                 │
    │  ╭────╮  ← 작은 핵도 탐지?      │
    │  │ ●  │                         │
    │  ╰────╯                         │
    │                                 │
    │     ╭───────────────╮           │
    │     │    DAPI but    │ ← 놓친 핵? │
    │     │   no mask!     │           │
    │     ╰───────────────╯           │
    └─────────────────────────────────┘

    핵심 확인 사항:
    ① cyan contour가 DAPI 밝은 영역과 정확히 일치하는지
    ② 놓친 핵 (DAPI 있지만 mask 없음) = under-segmentation
    ③ 잘못된 탐지 (mask 있지만 DAPI 없음) = over-segmentation
    ④ 병합된 핵 (두 개가 하나로) = insufficient resolution
```

**해석법**:
- **정확한 segmentation**: 거의 모든 밝은 DAPI 영역에 cyan contour 존재
- **Under-segmentation**: DAPI가 밝지만 mask가 없는 영역 → diameter 파라미터 조정 필요
- **Over-segmentation**: 하나의 핵이 여러 mask로 분할됨 → diameter가 너무 작음
- **Merged nuclei**: 인접한 핵이 하나의 mask로 합쳐짐 → diameter가 너무 큼 또는 overlap resolution 부족
- 논문에서 "Cellpose (v2.2.3) with diameter=None (auto-detection)" → auto가 대부분 잘 작동

---

### 10.3 Zoomed ROI (`dapi_zoomed_roi.png`) → 논문 Fig. 3c 관련

**논문 위치**: Fig. 3c (p.817) - "ROI comparison of segmentation methods"

**무엇을 봐야 하는가**:
```
    ┌─────────────── 100 µm ──────────────┐
    │                                       │
    │   ╭──╮  ╭──╮  ╭──╮                   │
    │   │○ │  │○ │  │○ │ ← 개별 핵 경계     │
    │   ╰──╯  ╰──╯  ╰──╯                   │
    │      ·· · ·· ·  · ··                  │  ← Transcript dots
    │   ╭──╮                                │
    │   │○ │  ← 핵 + 확장 경계              │
    │   ╰──╯                                │
    │                                       │
    │  확인: 핵 크기 ~10 µm, 간격 적절      │
    └───────────────────────────────────────┘

    핵심 확인:
    ① 핵 크기가 생물학적으로 합리적인지 (직경 5-15 µm)
    ② 인접 핵 간 간격이 적절한지
    ③ 확장된 mask가 이웃 세포와 겹치지 않는지
    ④ Transcript 점이 올바른 세포에 할당되는지
```

**해석법**:
- **핵 직경 ~10 µm**: 정상적인 뇌 세포 핵 크기
- **직경 < 3 µm**: debris이거나 over-segmentation
- **직경 > 30 µm**: 병합되었거나 혈관 등 비세포 구조
- 이 ROI 뷰는 Step 6의 segmentation 비교에서도 동일한 영역으로 사용됨

---

### 10.4 Mask QC Histogram (`mask_qc.png`)

**무엇을 봐야 하는가**:
```
    Count
    ↑
    │        ██
    │      ██████
    │    ██████████
    │  ████████████████
    │████████████████████▓▓
    └──────────────────────────→ Mask area (pixels²)
      0   100  500  1000  5000

    ① 피크 위치: 정상 핵 크기에 대응하는 면적
    ② 왼쪽 꼬리: 매우 작은 mask = debris 또는 fragmentation
    ③ 오른쪽 꼬리: 매우 큰 mask = merged nuclei 또는 artifact
    ④ 이봉 분포: 두 종류의 세포 크기 (neuron vs glia)
```

**해석법**:
- **정규 분포**: 균일한 세포 크기 → 좋은 segmentation
- **오른쪽 극단값**: log-transform 후 확인. 매우 큰 mask는 품질 문제
- **왼쪽 극단값**: debris → QC 필터링 (min_counts, min_genes)으로 제거
- **뇌 조직 기대값**:
  - Neuron 핵: 직경 ~12-15 µm (면적 ~115-180 µm²)
  - Glia 핵: 직경 ~5-10 µm (면적 ~20-80 µm²)
  - 2배 이상 큰 mask: 세포 병합 의심

---

### 10.5 Mask + Transcript Overlay (`mask_transcript_overlay.png`)

**무엇을 봐야 하는가**:
```
    ┌─────────────────────────────┐
    │  ╭───╮ ·  ·  · ╭───╮      │
    │  │·· │ ·    ·   │·  │      │
    │  │ · │  ·  ·    │ · │      │  Mask contour + transcript dots
    │  ╰───╯   ·      ╰───╯      │
    │      ·         ·  ·         │  ← 미할당 transcript
    │   ·     ·                   │     (mask 밖)
    └─────────────────────────────┘

    ① Mask 내부의 transcript 밀도가 외부보다 높은지
    ② Mask 경계 근처에 많은 transcript → 경계 정확도 중요
    ③ Mask 외부의 고립된 transcript → 미할당 reads (Step 5에서 처리)
```

**해석법**:
- **좋은 할당**: transcript 대부분이 mask 내부에 위치
- **누출 (leakage)**: mask 경계 바로 밖에 많은 transcript → 확장 (expansion)이 필요
- **미할당 비율**: 20-30%가 정상 (논문: "76.8% assigned"). 50% 이상 미할당이면 문제

---

## 9. Summary

Step 3은 **Cellpose 기반 독립적 세포 재분할**을 수행하는 핵심 단계이다. 생성된 segmentation masks와 재할당된 AnnData는 Step 4 (품질 비교), Step 5 (최적 확장), Step 6 (벤치마크)의 입력이 된다.

| 컴포넌트 | 구현 상태 | 심각도 |
|---------|---------|--------|
| Device detection | 완전 (CUDA > MPS > CPU 자동) | OK |
| DAPI loading | 완전 | OK |
| Cellpose segmentation | 완전 (tiled, v4 호환) | OK |
| Mask expansion | 완전 (`expand_labels`) | OK |
| Transcript assignment | 완전 (dual assignment: nuclear + expanded) | OK |
| Control probe filtering | 완전 (NegControl/BLANK/antisense) | OK |
| AnnData construction | 완전 (nuclear label 기반 matrix) | OK |
| Edge artifact flagging | 완전 (QC 전 수행, 인덱스 일치 보장) | OK |
| QC filtering | 완전 (config 기반 min_counts/min_genes) | OK |
| Domain assignment | 완전 (선택적, Shapely Point-in-Polygon) | OK |
| 시각화 (7개) | 완전 (DAPI overlay, ROI, transcript overlay 포함) | OK |

**구현 완성도**: **HIGH** - 논문/노트북 대비 완전히 구현. 추가 기능: tiled processing, dual assignment, edge flagging, control probe filtering, Cellpose v4 호환

---

## 11. 논문 Figure 직접 대응 및 시각화 정상 판별 종합 가이드

### 11.1 파이프라인 출력 → 논문 Figure 매핑 종합표

| # | 파이프라인 출력 파일 | 논문 Figure | 논문 페이지 | 논문 원문 설명 |
|---|---|---|---|---|
| 1 | `dapi_overview.png` | 직접 대응 없음 (입력 확인용) | - | DAPI 전체 이미지 품질 확인 |
| 2 | `dapi_mask_overlay.png` | **Extended Data Fig. 4g** (p.829) | 829 | "Resegmented datasets ROI" - Cellpose mask와 DAPI 오버레이 |
| 3 | `dapi_zoomed_roi.png` | **Fig. 3c** (p.818) 관련 | 818 | "ROI comparison of 6 segmentation methods" - 확대 뷰 |
| 4 | `mask_qc.png` | 직접 대응 없음 (QC용) | - | Mask 크기 분포 히스토그램 |
| 5 | `mask_transcript_overlay.png` | 직접 대응 없음 (검증용) | - | Mask + transcript scatter overlay |
| 6 | `dapi_transcript_overlay.png` | 직접 대응 없음 (검증용) | - | DAPI + transcript scatter overlay |
| 7 | `domain_polygon_overlay.png` | 직접 대응 없음 (선택) | - | Domain polygon + cell centroids |
| 8 | `resegmented_masks.tif` | → Step 6 **Fig. 3c-d** | 818 | Step 6 벤치마크의 비교 대상 |
| 9 | `resegmented.h5ad` | → Step 4, 5, 6 입력 | - | 재분할 AnnData (Step 4-6 공유) |

### 11.2 정상 결과 판별 체크리스트

#### 시각화 1: `dapi_overview.png` (DAPI 전체 이미지)
- [ ] **개별 핵 가시성**: 확대 없이도 밝은 핵이 보이는 영역 존재
- [ ] **조직 구조**: gray/white matter 경계, cortical layers 등 보임
- [ ] **아티팩트 없음**: 기포, 접힘, 과도한 자가형광 없음
- [ ] **밝기 균일성**: 이미지 전체에서 비슷한 배경 밝기
- **왜 정상인가**: Cellpose의 입력이므로, DAPI 이미지 품질이 segmentation 품질을 결정. 좋은 DAPI 이미지는 개별 핵이 명확히 구분되고 배경이 균일.
- **읽는법**: 밝은 점/영역=핵(DNA 염색). 밝기가 높고 균일할수록 Cellpose가 정확하게 segmentation 가능. 어두운 영역=세포가 없거나 접근 불가.

#### 시각화 2: `dapi_mask_overlay.png` (DAPI + Mask Overlay) → ExtData Fig. 4g
- [ ] **Cyan contour가 DAPI 밝은 영역과 일치**: 핵 위치에 mask 존재
- [ ] **놓친 핵 < 5%**: DAPI 있지만 mask 없는 영역 (under-segmentation)
- [ ] **잘못된 탐지 < 3%**: mask 있지만 DAPI 없는 영역 (over-segmentation)
- [ ] **병합 핵 < 5%**: 인접 핵이 하나의 mask로 합쳐진 경우
- **왜 정상인가**: 논문 Extended Data Fig. 4g에서 Cellpose resegmented ROI 확인. Cellpose nuclei model은 DAPI에서 높은 정확도로 핵을 탐지. "Cellpose (v2.2.3) deep-learning models 'nuclei' (CPn) and 'cyto' (CPc) were applied on the DAPI channels" (p.817 Methods).
- **읽는법**: 회색 배경=DAPI, Cyan 윤곽=expanded mask, Orange 윤곽=nuclei mask. Cyan이 DAPI 밝은 영역을 정확히 둘러쌀수록 좋음.
- **비정상 신호**: 대규모 미검출 → diameter 파라미터 조정 필요, 과도한 병합 → diameter가 너무 큼

#### 시각화 3: `dapi_zoomed_roi.png` (확대 ROI) → Fig. 3c 관련
- [ ] **핵 직경**: 5-15 µm (생물학적으로 합리적)
- [ ] **인접 핵 간격**: 겹침 없이 적절한 간격
- [ ] **확장 mask**: 이웃 세포와 과도하게 겹치지 않음
- [ ] **Transcript 점이 올바른 세포에 할당**
- **왜 정상인가**: 논문 Fig. 3c에서 6가지 segmentation 방법의 ROI를 비교. Cellpose nuclei가 개별 핵을 정확히 탐지하는 모습을 보여줌. 뇌 조직에서 neuron 핵 ~12-15 µm, glia 핵 ~5-10 µm이 정상.
- **읽는법**: 3패널 구성 - DAPI / 경계선 / 색상 mask. 경계선이 핵 경계를 정확히 따르는지 확인. 색상 mask에서 각 세포가 고유 색상으로 분리되는지 확인.
- **비정상 신호**: 직경 < 3 µm (debris), 직경 > 30 µm (병합)

#### 시각화 4: `mask_qc.png` (Mask 크기 분포)
- [ ] **단봉 분포**: 주 피크가 정상 핵 크기에 대응
- [ ] **왼쪽 꼬리 적음**: 극단적으로 작은 mask(debris) 비율 낮음
- [ ] **오른쪽 꼬리 적음**: 극단적으로 큰 mask(병합) 비율 낮음
- **왜 정상인가**: Cellpose의 nuclei model은 핵 크기를 학습했으므로, 정상적인 크기 범위의 mask를 생성해야 함. QC 필터링(min_counts, min_genes)이 debris 제거.
- **읽는법**: X축=mask 면적(pixels²), Y축=빈도. 피크가 정상 핵 크기(glia ~20-80 µm², neuron ~115-180 µm²)에 위치해야 함. 이봉이면 neuron과 glia 크기 차이 반영.

#### 시각화 5: `mask_transcript_overlay.png` (Mask + Transcript)
- [ ] **Mask 내부 transcript 밀도 > 외부**: 대부분의 reads가 세포 내
- [ ] **Mask 경계 근처의 transcript**: 경계에서의 전이가 자연스러움
- [ ] **Mask 외부 reads < 30%**: 논문 "76.8% assigned"
- **왜 정상인가**: 정확한 segmentation이면 transcript의 대부분이 mask 내부에 위치. 외부 reads는 미할당(extracellular) 또는 segmentation이 놓친 세포의 것.
- **읽는법**: 반투명 mask 위에 빨간 점(transcript). 점이 mask 안에 많으면 좋은 할당. 점이 mask 밖에 클러스터를 이루면 놓친 세포 존재.

### 11.3 논문 원문 인용 (시각화 관련)

| 시각화 | 논문 원문 인용 | 페이지 |
|---|---|---|
| Cellpose 적용 | "The Cellpose (v2.2.3) deep-learning models 'nuclei' (CPn) and 'cyto' (CPc) were applied on the DAPI channels with diameter parameters of none, 20, 30 and 40" | 822 (Methods) |
| Resegmented ROI (ExtData Fig. 4g) | "Resegmented datasets: ROI showing nuclei detected by Cellpose" | 829 caption |
| Segmentation 비교 (Fig. 3c) | "ROI of the six segmentation methods on mouse brain cortex dataset" | 818 caption |
| 확장 (expansion) | "Xenium's nuclear segmentation is followed by a default radius expansion of 15 µm" | 817 |
| Reads assignment | "76.8% of reads being assigned to cells" | 813 |
| Mask quality | "Baysor and Cellpose outperform standard Xenium segmentation" | 817 |
| Best pipeline | "the optimal algorithm involves two steps: first, identifying nuclei using Cellpose and second, assigning reads to individual cells using Baysor" | 821 |
