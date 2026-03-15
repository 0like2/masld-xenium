# Step 2 & 3 Segmentation 검토 보고서: Cellpose Overlap, 85um Expansion, 결과 비교

> 작성일: 2026-02-27
> 교수님 피드백 3가지에 대한 기술적 검토

---

## 1. Cellpose는 Overlap 없이 동작한다 (Non-overlapping)

### Cellpose 논문 근거 (Stringer et al., Nature Methods 2020)

Cellpose의 핵심 메커니즘은 **gradient flow tracking**이다:

1. **Neural network 출력**: 수평 gradient, 수직 gradient, inside/outside 확률맵 (3개)
2. **Mask recovery 과정**:
   - probability > 0.5인 픽셀만 선택
   - 각 픽셀에서 gradient 방향을 따라 200 iterations 동안 이동 (dynamical system)
   - **같은 fixed point(attractor)로 수렴하는 픽셀들 = 하나의 세포**
   - "all pixels belonging to a given cell can be routed to its center"
3. **결론**: 각 픽셀은 **정확히 하나의 gradient path**를 따르므로, **하나의 세포에만 할당**됨

> **Cellpose 출력은 integer label image로, 각 픽셀 값 = 세포 ID (0=배경). Overlap 불가능.**

### Cellpose 논문 핵심 그림 가이드

> 출처: Stringer, C. et al. "Cellpose: a generalist algorithm for cellular segmentation."
> *Nature Methods* 18, 100-106 (2021). https://doi.org/10.1038/s41592-020-01018-x

#### Fig. 1a — Heat Diffusion Simulation (학습 데이터 생성 원리)

**무엇을 보여주는가**: 수동 annotation된 세포 mask를 neural network이 학습 가능한 **vector flow 표현**으로 변환하는 과정

**그림 구성**:
- 왼쪽: 세포 mask (이진 형태, 불규칙한 모양)
- 중앙: mask 중심에 heat source를 놓고 **열확산(heat diffusion) 시뮬레이션** 실행 → 에너지 함수의 등고선(heatmap)이 생성됨
- 오른쪽: 에너지 함수의 X gradient, Y gradient를 화살표(arrows)로 표시

**핵심 포인트**:
- 열확산의 결과로 모든 gradient가 **하나의 꼭대기(peak = 세포 중심)**를 향함
- 불규칙한 모양의 세포도 gradient가 "모서리를 돌아서(indirectly around corners)" 중심으로 향함
- X/Y gradient를 0-360도 각도로 변환 → **단일 방향 표현** (이것이 network가 예측할 target)
- **결론**: 세포 중심이 유일한 global maximum이므로, gradient를 따라가면 반드시 하나의 중심에 도달

**교수님께 설명 시**: "이 그림은 Cellpose가 학습하는 것이 세포 경계가 아니라 **세포 중심으로 향하는 흐름장(flow field)**임을 보여줍니다. 하나의 세포 안에서는 모든 화살표가 같은 중심을 가리킵니다."

---

#### Fig. 1b — 다양한 세포 형태에 대한 Flow 시각화

**무엇을 보여주는가**: 학습 데이터의 여러 세포 형태에 대한 spatial flow 예시

**그림 구성**:
- 여러 세포 이미지 위에 gradient flow를 **색상 원형(color wheel)**으로 표현
- 각 픽셀의 색상 = gradient 방향 (0-360도를 HSV/sinebow 색상환으로 매핑)
  - 예: 오른쪽을 향하는 flow = 빨강, 위쪽 = 초록, 왼쪽 = 파랑 등
- 밝기 = gradient의 크기(magnitude)

**핵심 포인트**:
- 둥근 세포, 긴 세포, L자 모양 등 **모든 형태**에서 flow가 중심을 향함
- 세포 경계 근처: 색상이 급격히 변함 (인접 세포가 다른 중심을 향하므로)
- 같은 세포 내부: 색상이 부드럽게 변화 (같은 중심을 향하므로)

**교수님께 설명 시**: "같은 색상 영역 = 같은 세포. 색상이 급변하는 경계 = 세포 간 경계. 이 flow를 neural network이 예측하는 것입니다."

---

#### Fig. 1c-d — Neural Network 구조 (U-Net 기반)

**무엇을 보여주는가**: Cellpose network의 architecture

**그림 구성**:
- U-Net 형태의 encoder-decoder 구조 다이어그램
- 입력: 원본 이미지 (DAPI 등)
- 출력: **3개 채널** — (1) X flow, (2) Y flow, (3) cell probability map
- Style vector: 가장 낮은 해상도에서 추출한 이미지 스타일 정보 → 모든 upsampling 층에 전달

**핵심 포인트**:
- 출력이 3개이므로 "어디가 세포인가(확률맵)" + "세포 내부에서 중심이 어디인가(flow)" 동시 예측
- Cell probability > 0.5인 픽셀만 flow tracking 대상 → 배경 픽셀은 애초에 제외

---

#### ★ Fig. 1e — Pixel Dynamics: 핵심 그림 (Non-overlapping 증명) ★

**무엇을 보여주는가**: Test time에 예측된 flow field 위에서 **각 픽셀이 이동하는 과정(dynamical system)**

**그림 구성**:
- 예측된 flow field 위에 여러 색상의 **점(dot)** 또는 **궤적(trajectory)**이 표시됨
- 각 점은 하나의 픽셀을 나타냄
- 화살표/궤적은 해당 픽셀이 flow를 따라 이동하는 경로를 보여줌
- **같은 세포에 속하는 픽셀들 → 같은 색상** (같은 fixed point로 수렴하므로)
- 서로 다른 세포의 픽셀들 → 서로 다른 색상 (다른 fixed point로 수렴)

**Mask recovery 알고리즘 (이 그림이 보여주는 4단계)**:
```
Step 1: cell probability map에서 threshold(0.5) 적용 → 세포 후보 픽셀만 선택
Step 2: 각 후보 픽셀 위치에서 시작, predicted flow field를 따라 이동
        (finite differences, step size = 1, 200 iterations)
Step 3: 200회 반복 후 각 픽셀이 도달한 최종 위치(fixed point) 기록
Step 4: 같은 fixed point에 도달한 픽셀들 → 같은 cell ID 부여
```

**핵심 포인트**:
- 각 픽셀은 **정확히 하나의 경로**를 따름 → **정확히 하나의 fixed point**에 도달
- 하나의 경로가 두 개의 fixed point로 분기하는 것은 **수학적으로 불가능** (deterministic dynamical system)
- 따라서 각 픽셀은 **정확히 하나의 세포**에 할당됨 → **overlap 원천 불가능**
- Fixed point = 세포의 중심(attractor), Basin of attraction = 세포 영역

**교수님께 설명 시**: "이 그림에서 같은 색 점들이 하나의 중심으로 모이는 것을 보시면 됩니다. 한 점이 동시에 두 중심으로 갈 수는 없으므로, Cellpose는 구조적으로 overlap이 불가능합니다. 이것이 Cellpose가 integer label mask를 출력하는 이유입니다."

---

#### Fig. 1f — 최종 Mask 결과: Non-overlapping 확인

**무엇을 보여주는가**: Fig. 1e의 pixel dynamics 결과로 생성된 **최종 segmentation mask**

**그림 구성**:
- 원본 이미지 위에 세포별로 **다른 색상**의 mask가 오버레이됨
- 각 색상 = 하나의 cell ID (integer label)
- 색상 간 **겹침 없음** — 모든 세포 픽셀은 정확히 하나의 색상
- 배경(cell probability < 0.5)은 label = 0 (색상 없음)

**핵심 포인트**:
- Fig. 1e에서 같은 fixed point로 수렴한 픽셀들 → Fig. 1f에서 같은 색상
- 세포 경계가 **sharp하고 겹치지 않음** → Voronoi와 유사하지만 flow 기반
- 출력 형식: `np.ndarray (H x W)`, dtype=int, 값 = {0, 1, 2, ..., N_cells}

**교수님께 설명 시**: "최종 결과물이 이 그림입니다. 각 픽셀에 정수 ID가 부여되어 있고, 한 픽셀이 두 개의 세포에 동시에 속하는 것은 불가능합니다."

---

#### Fig. 5d-e — 3D에서도 동일 원리 적용 (참고용)

**무엇을 보여주는가**: 2D로 학습된 Cellpose를 3D 데이터에 확장하는 방법

**그림 구성**:
- Fig. 5d: XY, XZ, YZ 슬라이스에서 각각 2D flow 예측 → 3개의 flow를 pairwise 평균하여 3D flow field (X, Y, Z) 생성
- Fig. 5e: 3D pixel dynamics — 3D 공간에서 각 픽셀(voxel)이 flow를 따라 fixed point로 수렴

**핵심 포인트**:
- 3D에서도 **같은 non-overlapping 원리** 적용 — 각 voxel은 하나의 fixed point로만 수렴
- 단, 파이프라인의 Xenium 데이터는 **2D DAPI 이미지**만 사용하므로 3D 적용은 해당 없음
- 이 그림은 "Cellpose가 3D overlap을 고려할 수 있다"는 것을 보여주지만, 현재 파이프라인에서는 활용하지 않음

**교수님께 설명 시**: "Cellpose 자체는 3D도 지원하지만, Xenium 데이터에서는 2D DAPI만 사용하므로 조직의 z축 겹침은 고려하지 못합니다. 이것이 read-based 방법(Baysor)이 필요한 이유 중 하나입니다."

---

#### Supplementary Fig. 1 (S1) — Cellpose GUI의 Flow 시각화 (참고용)

**무엇을 보여주는가**: Cellpose GUI에서 실제 output을 시각화한 화면

**그림 구성**:
- 탭 전환으로 (1) 원본 이미지, (2) flow field (RGB 색상환), (3) cell probability map, (4) mask outlines 확인 가능
- Flow field: 색상환(sinebow) 방식으로 각 픽셀의 flow 방향을 색상으로 표현
  - `arctan2(dy, dx)` → 각도 → RGB 색상 (0도/120도/240도 cosine offset)
- Cell probability: 밝은 영역 = 세포 내부(높은 확률), 어두운 영역 = 배경
- Mask outlines: 빨간색 윤곽선으로 최종 segmentation 경계 표시

**핵심 포인트**:
- 실제 Cellpose를 돌렸을 때 나오는 output이 어떻게 생겼는지 보여줌
- Flow field의 색상 패턴: 같은 세포 내부는 부드러운 색상 전환, 세포 경계에서 급격한 색상 변화

---

### 그림 활용 요약

| 설명 목적 | 추천 그림 | 핵심 관찰 포인트 |
|----------|----------|----------------|
| Cellpose가 학습하는 것이 무엇인가 | **Fig. 1a** | 화살표가 세포 중심을 향함 (경계가 아닌 중심!) |
| Flow가 세포 형태에 무관하게 작동함 | **Fig. 1b** | 다양한 형태에서 색상이 중심으로 수렴 |
| **Overlap 불가능의 핵심 증명** | **★ Fig. 1e ★** | 각 점이 하나의 fixed point로만 수렴 (분기 없음) |
| 최종 출력이 integer label mask | **Fig. 1f** | 색상 겹침 없음, 각 픽셀 = 하나의 cell ID |
| 3D 확장 가능성 vs 현재 한계 | **Fig. 5d-e** | 3D 지원되지만 파이프라인은 2D만 사용 |

> **교수님 발표/설명 시 최소 필수 그림**: Fig. 1a (원리) + Fig. 1e (non-overlapping 증명) + Fig. 1f (결과)

### expand_labels: Cellpose와 완전히 다른 원리, 다른 라이브러리

#### Cellpose vs expand_labels 비교

| | **Cellpose** (핵 검출) | **expand_labels** (핵 확장) |
|---|---|---|
| **라이브러리** | `cellpose` (MouseLand/cellpose) | `scikit-image` (`skimage.segmentation`) |
| **내부 핵심 함수** | 자체 neural network + `dynamics.py` | `scipy.ndimage.distance_transform_edt` |
| **원리** | 딥러닝이 gradient flow 예측 → pixel dynamics (200 iterations) | 단순 유클리드 거리 계산 → nearest-label lookup |
| **입력** | DAPI 이미지 (raw pixels) | Cellpose가 출력한 nuclei_masks (integer labels) |
| **하는 일** | 핵이 어디에 있는지 **찾는 것** (segmentation) | 찾아놓은 핵을 주변으로 **넓히는 것** (post-processing) |
| **비유** | "이 사진에서 사람 얼굴을 찾아라" (AI 인식) | "찾은 얼굴 주변에 원을 그려라" (기계적 도형 연산) |
| **Non-overlap 이유** | 각 픽셀이 하나의 flow path만 따르므로 | 각 배경 픽셀의 nearest label이 하나뿐이므로 |

```
파이프라인 Step 3 내 역할 분담:

  DAPI 이미지
      |
      ↓
  [cellpose 라이브러리]        ← 딥러닝: 핵을 "찾는" 단계
  model = CellposeModel(pretrained_model='nuclei')
  masks, flows, styles = model.eval(dapi_image)
      |
      ↓
  nuclei_masks (핵 영역만 있는 integer label image)
      |
      ↓
  [scikit-image 라이브러리]    ← 기계적 연산: 핵을 "넓히는" 단계
  masks = expand_labels(nuclei_masks, distance=400)
      |
      ↓
  expanded_masks (핵 + 주변 Voronoi 영역)
```

**핵심**: 이 두 단계는 **완전히 별개의 라이브러리, 별개의 알고리즘**이다.
- Cellpose는 어디까지나 **핵(nuclei) 검출**만 담당
- `expansion_distance=400`은 Cellpose에 전달되는 것이 아니라, Cellpose **이후에** scikit-image의 거리 계산 함수에 전달됨

#### expand_labels 내부 동작: 왜 기계적 확장인데 겹치지 않는가

#### 왜 기계적 확장인데 겹치지 않는가? — 알고리즘 내부 동작

`expand_labels`는 "각 label을 1px씩 동시에 키우다가 부딪히면 멈추는" 방식이 **아니다**.
실제로는 **distance transform + nearest-label lookup**을 한 번에 계산한다:

```
알고리즘 (scikit-image 내부 구현):

Step 1: 모든 배경 픽셀(=0)에 대해 distance transform 계산
        → 각 배경 픽셀에서 "가장 가까운 label 픽셀까지 몇 px인가" 기록
        → 동시에 "그 가장 가까운 label의 ID가 무엇인가"도 기록
        (scipy.ndimage.distance_transform_edt 사용)

Step 2: distance ≤ expansion_distance인 배경 픽셀만 선택

Step 3: 해당 픽셀에 nearest label ID를 부여
```

**핵심: 각 배경 픽셀은 "가장 가까운 label 하나"만 갖는다.**

```
예시: 세포 A와 B 사이의 빈 공간에 있는 픽셀 P

        A까지 12px     B까지 8px
  [A]-------- P ----------[B]

  → P는 B에 할당 (거리가 더 가까우므로)
  → A와 B가 동시에 P를 차지하는 것은 불가능 (nearest가 하나뿐)
```

정확히 **등거리(equidistant)** 지점이 있어도, 구현상 하나의 label만 선택된다 (tie-breaking):

```
확장 전:                         확장 후 (expansion_distance 내):

  .......A.......                AAAAAAA|BBBBBBB
  .................              AAAAAAA|BBBBBBB
  ...........B...                AAAAA|BBBBBBBBB
                                      ↑
                                 Voronoi 경계 (등거리선)
                                 이 선 위의 픽셀도 둘 중 하나에만 할당됨
```

**결과가 Voronoi tessellation과 동일한 이유**:
- Voronoi diagram = "각 점에 가장 가까운 seed를 기준으로 공간을 분할"
- expand_labels = "각 배경 픽셀에 가장 가까운 label을 부여"
- 수학적으로 **동일한 연산** (seed = label centroid가 아니라 label boundary이므로 정확히는 "일반화된 Voronoi")

따라서 expand_labels의 출력도 Cellpose와 마찬가지로 **integer label image**이며, 한 픽셀이 두 세포에 동시에 속하는 것은 알고리즘 구조상 불가능하다.

```python
# xenium_step3_resegmentation.py:402-404
nuclei_masks = masks  # Keep original nuclei masks for dual assignment
if expansion_distance > 0:
    masks = expand_labels(nuclei_masks, distance=expansion_distance)
# 결과: 여전히 integer label image, 각 픽셀 = 하나의 세포 ID
```

### 전체 파이프라인 흐름

```
DAPI 이미지
    |
[Cellpose 'nuclei' model]  -> nuclei_masks (non-overlapping integer labels)
    |
[expand_labels(distance=400px)]  -> expanded_masks (Voronoi 분할, non-overlapping)
    |
Transcript 할당:
  - in_cell = nuclei_masks[y, x]        <- 핵 내부 transcript만
  - closest_cell = expanded_masks[y, x]  <- 확장 영역 포함
```

### 교수님 질문에 대한 답변

> "오블랩은 있는지 모르겠고 싱글로 가정해가지고 이렇게 그냥 파티션을 맞는다. 맞아 틀려?"

**맞습니다.** Cellpose + expand_labels 조합은:
- 3D overlap을 고려하지 않음 (2D 평면만)
- 각 픽셀은 하나의 세포에만 할당 (single assignment)
- expand_labels는 Voronoi-like partition으로 빈 공간을 채움
- 실제 조직의 3D 겹침은 무시됨 -> 이것이 Baysor 같은 read-based 방법이 필요한 이유

**한계**: Xenium 논문 (Salas et al., 2025) p.822에서도 인정:
> "implementing new segmentation algorithms that account for the 3D structure of the data and incorporating additional staining for cellular membranes would facilitate the correct identification of individual cells"

---

## 2. 85um Expansion: 논문 근거 없음, 수정 필요

### 현재 파이프라인 설정

| 항목 | 값 |
|------|-----|
| `config.yaml` expansion_distance | 400 pixels |
| um_per_pixel_inv | 4.70588 px/um |
| **실제 확장 거리** | **400 / 4.70588 = 85 um** |

### Xenium 논문 (Salas et al., 2025) 수치

| 항목 | 값 | 출처 |
|------|-----|------|
| Xenium default expansion | **15 um** | p.817 "default radius expansion of 15 um" |
| 평균 nuclei 반경 | 5.06 um | p.817 "nuclei...presented a radius of 5.06 um" |
| 배경 신호 전환점 (turnover) | 10.71 um | p.817 Fig. 3a "transcripts located more than 10.71 um" |
| **논문 최적 expansion** | **5.64 um** | p.817 "ideal expansion should be 5.64 um" (= 10.71 - 5.06) |
| 테스트한 범위 | 1, 2, 5, 10, 15 um | p.817 Fig. 3e |

### 비교

| 설정 | 거리 | 픽셀 | 논문 최적 대비 |
|------|------|------|-------------|
| **현재 파이프라인** | **85 um** | **400 px** | **15.1배 큼** |
| 논문 default | 15 um | ~71 px | 2.7배 큼 |
| **논문 최적값** | **5.64 um** | **~27 px** | **기준** |

### expansion_distance 값의 코드 내 흐름 (변수 추적)

이 값이 config에서 읽혀 최종 mask에 적용되기까지의 전체 경로:

```
[1] config.yaml (line 122)
    resegmentation:
      expansion_distance: 400    ← 설정값 (단위: pixels)
      um_per_pixel_inv: 4.70588  ← 픽셀/um 변환 계수

            ↓ config dict로 로드

[2] xenium_step3_resegmentation.py:323
    reseg_config = config.get('resegmentation', {})

            ↓

[3] xenium_step3_resegmentation.py:333
    expansion_distance = reseg_config.get('expansion_distance', 400)
    # 변수명: expansion_distance (int, 단위: pixels)
    # default fallback = 400

            ↓ Cellpose 실행 (expansion과 무관)

[4] xenium_step3_resegmentation.py:380-390
    model = models.CellposeModel(gpu=use_gpu, pretrained_model='nuclei', device=device)
    masks, flows, styles = _run_cellpose_tiled(model, dapi_image, ...)
    # Cellpose 출력: masks (integer label image, 각 픽셀 = 핵 ID)
    # *** Cellpose 자체는 expansion_distance를 사용하지 않음 ***

            ↓ skimage label()로 재라벨링

[5] xenium_step3_resegmentation.py:398
    masks = label(masks)
    # connected component labeling (정수 라벨 재부여)

            ↓ 원본 핵 마스크 보존

[6] xenium_step3_resegmentation.py:401
    nuclei_masks = masks   ← 핵 전용 마스크 (확장 전 원본)

            ↓ ★ expansion_distance가 적용되는 핵심 지점 ★

[7] xenium_step3_resegmentation.py:402-404
    if expansion_distance > 0:
        masks = expand_labels(nuclei_masks, distance=expansion_distance)
    #                                       ^^^^^^^^^^^^^^^^^^^^^^^^
    #   skimage.segmentation.expand_labels의 `distance` 파라미터로 전달
    #   distance=400 → 각 label을 최대 400 pixels까지 Voronoi 확장

            ↓ 확장된 masks로 transcript 할당

[8] xenium_step3_resegmentation.py:511-519  (transcript → cell 매핑)
    y_coords = (y_vals * um_per_pixel_inv).astype(int)  # um → pixel 변환
    x_coords = (x_vals * um_per_pixel_inv).astype(int)

    xenium_step3_resegmentation.py:524-527  (dual assignment)
    in_cell     = nuclei_masks[y, x]   ← 핵 내부만 (확장 전)
    closest_cell = masks[y, x]         ← 확장 영역 포함 (Voronoi)

            ↓ edge cell 판별에도 사용

[9] xenium_step3_resegmentation.py:674
    if nuclei_masks is not None and expansion_distance > 0:
        # 이미지 가장자리에 닿는 expanded label → edge artifact로 flagging
```

**핵심 요약**:
- `expansion_distance`는 **Cellpose에 직접 전달되지 않음** (Cellpose는 핵만 검출)
- Cellpose 출력(`nuclei_masks`) 이후 **별도 단계**에서 `expand_labels(distance=expansion_distance)`로 적용
- 즉, Cellpose의 segmentation 품질과는 무관하고, **확장 후처리(post-processing)**의 파라미터

| 단계 | 코드 위치 | 변수명 | 설명 |
|------|----------|--------|------|
| Config 로드 | `config.yaml:122` | `expansion_distance` (YAML key) | 400 (pixels) |
| Python 변수 | `step3:333` | `expansion_distance` (local var) | `reseg_config.get('expansion_distance', 400)` |
| expand_labels 호출 | `step3:404` | `distance` (함수 인자) | `expand_labels(nuclei_masks, distance=expansion_distance)` |
| 변환 계수 | `config.yaml:123` | `um_per_pixel_inv` | 4.70588 px/um |
| 실제 거리 | 계산값 | — | 400 / 4.70588 = **85 um** |

### 왜 85um가 문제인가

논문 Fig. 3a 기반 분석:
- centroid에서 **10.71um 이상** 떨어진 transcript -> 핵 signature보다 **배경 signature와 더 높은 상관관계**
- 85um까지 확장하면 대부분의 공간이 **배경 noise**로 채워짐
- cell type별 optimal expansion이 다름 (Fig. 3b) -> 고정 85um는 모든 cell type에 부적절

### 85um의 출처

원본 노트북 (`notebooks/3_techniques_comparison/3_1_resegmentation_notebooks/batch_segmentation_cellpose-Xenium.ipynb`):
```python
expanded_nuclei = expand_labels(labeled_nuclei, distance=400)  # 주석/근거 없음
```

**결론**: 노트북에 `distance=400`이 하드코딩되어 있으나 **근거 문서가 없음**. 논문에서 나온 값이 아니며, 빈 공간을 최대한 채우려는 보수적 초기값으로 추정.

### 권장 수정 (코드 변경 사항)

**파일**: `pipeline/config.yaml` (line 122)

```yaml
# 현재 (문제):
expansion_distance: 400  # 85um - 논문 근거 없음

# 수정안 1 (논문 default):
expansion_distance: 71   # 15um (논문 Xenium default)

# 수정안 2 (논문 최적값):
expansion_distance: 27   # 5.64um (논문 Fig. 3b optimal)

# 수정안 3 (Step 5 결과 기반 - 권장):
# Step 5 optimal expansion 분석 실행 후 cell type별 최적값 사용
```

---

## 3. 기계적 분할 vs 딥러닝 분할: 차이가 있어야 정상

### 논문의 분류 체계 (Fig. 3c-e)

| 카테고리 | 방법 | 원리 |
|---------|------|------|
| **Staining-based** | Xenium default, Cellpose, MESMER, Watershed | DAPI 형태 기반 |
| **Read-based** | Baysor, Clustermap | transcript 밀도/조성 기반 |
| **Mixed** | Baysor + Cellpose prior | 형태 + transcript 조합 |

### 왜 차이가 있어야 하는가

논문 Fig. 3d (ARI heatmap):
- **같은 카테고리** 내 방법들은 높은 ARI (유사한 결과)
- **다른 카테고리** 간에는 낮은 ARI (다른 결과)
- Staining-based끼리 cluster, Read-based끼리 cluster -> **정보 소스가 다르면 결과도 달라야 정상**

논문 p.817:
> "Staining-based strategies using DAPI generated similar outputs, with cell expansion being the force driving their differences."

### 차이의 의미

1. **차이가 없다면** -> 두 방법이 같은 정보만 사용 -> 하나가 불필요
2. **차이가 있다면** -> 각 방법이 다른 측면을 포착 -> **상보적 분석 가능**
3. **논문 결론**: Baysor + Cellpose nuclei prior 조합이 최고 성능 (BA2 P0.8)
   - Cellpose가 핵 형태 정보 제공
   - Baysor가 transcript 밀도 기반으로 세포 경계 결정
   - 두 정보가 합쳐져 가장 정확한 segmentation

### dapi_zoomed_roi.png이 보여주는 것

**파일**: `xenium_step3_resegmentation.py:238-311` (`_save_dapi_zoomed_roi()`)

3-panel 시각화:
1. **DAPI (zoomed ROI)**: 원본 DAPI 형광 이미지
2. **DAPI + boundaries**: cyan = expanded masks 경계, orange = nuclei masks 경계
3. **DAPI + colored masks**: 반투명 label 색상 오버레이

이 그림에서 확인해야 할 것:
- **cyan 경계가 핵(orange)보다 훨씬 넓음** -> 85um expansion의 시각적 증거
- **핵 사이 빈 공간이 모두 Voronoi로 분할됨** -> overlap 없이 partition
- **확장된 영역이 다른 세포의 transcript를 포함할 수 있음** -> misassignment 위험

---

## 검증 방법

1. config.yaml에서 expansion_distance 변경 후 Step 3 재실행
2. `dapi_zoomed_roi.png` 재생성하여 expansion 범위 시각적 확인
3. Step 4 (techniques comparison)에서 NMP, ARI 메트릭으로 segmentation 품질 비교
4. Step 5 optimal expansion으로 cell type별 최적값 산출

---

## 요약 (교수님 답변용)

| 피드백 | 결론 |
|-------|------|
| 85um expansion 근거 | **논문 근거 없음.** 논문 최적값은 5.64um, default도 15um. 노트북 하드코딩값(400px)이 검증 없이 사용됨 |
| Cellpose overlap 여부 | **Overlap 없음.** Gradient tracking -> 각 픽셀이 하나의 attractor로 수렴 -> non-overlapping. expand_labels도 Voronoi partition (overlap 불가) |
| 분할 결과 비교 | **차이가 있어야 정상.** Staining-based(Cellpose)와 Read-based(Baysor)는 다른 정보를 사용하므로 결과가 달라야 하며, 조합이 최고 성능 |
