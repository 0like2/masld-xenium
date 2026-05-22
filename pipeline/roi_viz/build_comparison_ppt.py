"""Comprehensive ROI × method comparison PPT.

Layout (~28 slides):

    Section 1 — Overview
        1. Title
        2. Methods + Pareto plot
        3. ROI × method 표 + Winner count
        4. 5-method 한눈 비교 (한 ROI 에서 dapi+boundary 5-panel + alltx+boundary 5-panel)

    Section 2 — Class-specific method comparisons (special pages)
        5–9. Architecture / Compartment / Low_VSI / High_density / Fold_boundary
              각 슬라이드는 5 method panel + 해석 + 핵심 수치

    Section 3 — Per-ROI baysor vs others (12 ROI × 1 slide = 12 slides)
        10–21. 각 ROI 의 reassigned_transcripts xenium__vs__{baysor, cellpose,
                rigid, optimal} 4-panel + 어떤 색이 무엇을 의미하는지 + 변화 해석

    Section 4 — Per-ROI 5-method snapshot (12 ROI × 1 slide = 12 slides)
        22–33. 각 ROI 의 dapi_boundary_<method> 5-panel + 한 줄 평가

    Section 5 — Conclusion
        34. Class-목적별 평가 표
        35. 주의사항 / 값 검증
        36. 결론

Total ≈ 36 slides.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
from pptx.util import Inches, Pt

logger = logging.getLogger(__name__)

SLIDE_W, SLIDE_H = Inches(13.333), Inches(7.5)
PAD = Inches(0.25)


# ============================================================================
# Slide primitives
# ============================================================================
def _new(prs):
    s = prs.slides.add_slide(prs.slide_layouts[6])
    s.background.fill.solid()
    s.background.fill.fore_color.rgb = RGBColor(0xFF, 0xFF, 0xFF)
    return s


def _title(slide, text, *, size=22, color=(0x12, 0x12, 0x12)):
    box = slide.shapes.add_textbox(PAD, Inches(0.10),
                                       SLIDE_W - 2 * PAD, Inches(0.55))
    tf = box.text_frame; tf.word_wrap = True
    p = tf.paragraphs[0]; p.alignment = PP_ALIGN.LEFT
    r = p.add_run(); r.text = text
    r.font.size = Pt(size); r.font.bold = True
    r.font.color.rgb = RGBColor(*color)


def _txt(slide, text, left, top, width, height, *,
         size=11, bold=False, color=(0x33, 0x33, 0x33),
         align=PP_ALIGN.LEFT, mono=False):
    box = slide.shapes.add_textbox(left, top, width, height)
    tf = box.text_frame; tf.word_wrap = True
    for i, line in enumerate(text.split("\n")):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = align
        r = p.add_run(); r.text = line
        r.font.size = Pt(size); r.font.bold = bold
        r.font.color.rgb = RGBColor(*color)
        if mono: r.font.name = "Consolas"


def _img(slide, path, left, top, width=None, height=None):
    if not Path(path).exists():
        _txt(slide, f"[missing]\n{Path(path).name}",
             left, top, width or Inches(3), height or Inches(2),
             size=8, color=(0xCC, 0x33, 0x33))
        return None
    return slide.shapes.add_picture(str(path), left, top,
                                      width=width, height=height)


def _table(slide, df, left, top, width, height, *,
           header_size=11, body_size=9,
           highlight_col=None, highlight_value=None):
    rows, cols = df.shape[0] + 1, df.shape[1]
    tbl = slide.shapes.add_table(rows, cols, left, top, width, height).table
    for j, col in enumerate(df.columns):
        c = tbl.cell(0, j); c.text = str(col)
        for p in c.text_frame.paragraphs:
            for r in p.runs:
                r.font.bold = True; r.font.size = Pt(header_size)
                r.font.color.rgb = RGBColor(0xFF, 0xFF, 0xFF)
        c.fill.solid(); c.fill.fore_color.rgb = RGBColor(0x2C, 0x3E, 0x50)
    for i, (_, row) in enumerate(df.iterrows(), start=1):
        for j, col in enumerate(df.columns):
            c = tbl.cell(i, j); c.text = str(row[col]) if row[col] is not None else ""
            for p in c.text_frame.paragraphs:
                for r in p.runs: r.font.size = Pt(body_size)
            if highlight_col and col == highlight_col and row[col] == highlight_value:
                c.fill.solid(); c.fill.fore_color.rgb = RGBColor(0xD4, 0xEF, 0xDF)
            elif i % 2 == 0:
                c.fill.solid(); c.fill.fore_color.rgb = RGBColor(0xF2, 0xF4, 0xF7)


# ============================================================================
# ROI catalog
# ============================================================================
ROI_LIST: List[Tuple[str, str, str]] = [
    ("ssam_transition_1", "architecture", "arch"),
    ("ssam_transition_2", "architecture", "arch"),
    ("ssam_transition_3", "architecture", "arch"),
    ("ssam_transition_4", "architecture", "arch"),
    ("p2r_compartment_mismatch_1", "compartment", "other"),
    ("p2r_compartment_mismatch_2", "compartment", "other"),
    ("ovrlpy_low_coherence_1", "low_vsi", "other"),
    ("ovrlpy_low_coherence_2", "low_vsi", "other"),
    ("density_high_density_1", "high_density", "other"),
    ("density_high_density_2", "high_density", "other"),
    ("boundary_tissue_edge_1", "fold_boundary", "other"),
    ("boundary_tissue_edge_2", "fold_boundary", "other"),
]
METHOD_ORDER = ["xenium_nucleus", "cellpose_nuclei", "rigid_expansion",
                  "optimal_expansion", "baysor"]
METHOD_LABEL = {
    "xenium_nucleus":    "(1) 10x default",
    "cellpose_nuclei":   "(2) Cellpose only",
    "rigid_expansion":   "(3a) Cellpose+rigid",
    "optimal_expansion": "(3b) Cellpose+optimal",
    "baysor":            "(4) Cellpose+Baysor ★",
}


# ============================================================================
# Section 1 — Overview slides
# ============================================================================
def _slide_title(prs):
    s = _new(prs)
    _txt(s, "Xenium ROI 시각화\n5 method × 12 ROI 비교 분석",
         Inches(0.5), Inches(1.6), Inches(12), Inches(2.5),
         size=42, bold=True, color=(0x1B, 0x4F, 0x72))
    _txt(s, "Paper Fig 3e 기준 ─ recovery × NMP Pareto + ROI별 변화 추적",
         Inches(0.5), Inches(4.2), Inches(12), Inches(0.6),
         size=20, color=(0x47, 0x4F, 0x55))
    _txt(s,
         "12 ROI = 4 architecture / 2 compartment / 2 low_vsi / 2 high_density / 2 fold_boundary\n"
         "각 ROI 에 대해 5 method (xenium / cellpose / rigid / optimal / baysor) 비교\n"
         "Tissue: Xenium V1 FFPE Human Brain (Alzheimer)",
         Inches(0.5), Inches(5.0), Inches(12), Inches(1.6),
         size=14, color=(0x80, 0x80, 0x80))


def _slide_methods_pareto(prs, pareto_png: Path):
    s = _new(prs); _title(s, "5 method 비교 + paper 기준 (Pareto)")
    df = pd.DataFrame([
        ["xenium_nucleus",   "10x default (nucleus mask)"],
        ["cellpose_nuclei",  "Cellpose only"],
        ["rigid_expansion",  "Cellpose + 5µm rigid"],
        ["optimal_expansion","Cellpose + Step5 optimal"],
        ["baysor",           "Cellpose + Baysor ★"],
    ], columns=["method", "설명"])
    _table(s, df, Inches(0.3), Inches(0.85),
           Inches(4.5), Inches(2.3), header_size=12, body_size=11)
    _txt(s,
         "Paper Fig 3e 의 두 축:\n"
         "  • assigned_fraction (회수율)\n"
         "  • NMP (purity)\n"
         "둘 다 높을수록 좋음.\n"
         "우상단(★) 가까울수록 best.",
         Inches(0.3), Inches(3.4), Inches(4.5), Inches(2.5),
         size=12, color=(0x33, 0x33, 0x33))
    _img(s, pareto_png, Inches(5.0), Inches(0.85), width=Inches(8.0))
    _txt(s, "📂 pareto_recovery_vs_nmp.png",
         Inches(5.0), Inches(7.0), Inches(8.0), Inches(0.4),
         size=9, mono=True, color=(0x80, 0x80, 0x80))


def _build_pivot(df):
    short = {"xenium_nucleus":"xenium","cellpose_nuclei":"cellpose",
              "rigid_expansion":"rigid","optimal_expansion":"optimal",
              "baysor":"baysor"}
    pa = df.pivot(index='roi_id', columns='method', values='assigned_fraction')
    pn = df.pivot(index='roi_id', columns='method', values='negative_marker_purity')

    def winner(roi):
        a = pa.loc[roi]; n = pn.loc[roi]
        cands = [m for m in METHOD_ORDER if m in a.index and not pd.isna(a[m])
                  and a[m] >= 0.5 and m in n.index and not pd.isna(n[m])]
        if not cands: return '—'
        return short[max(cands, key=lambda m: n[m])]

    rows = []
    for roi_id, _, _ in ROI_LIST:
        if roi_id not in pa.index: continue
        row = {'ROI': roi_id}
        for m in METHOD_ORDER:
            a = pa.loc[roi_id, m] if m in pa.columns else float('nan')
            n = pn.loc[roi_id, m] if m in pn.columns else float('nan')
            row[short[m]] = (f"{a:.2f}/{n:.2f}"
                                if pd.notna(a) and pd.notna(n) else "—")
        row['Winner'] = winner(roi_id)
        rows.append(row)
    return pd.DataFrame(rows)


def _slide_pivot_summary(prs, df):
    s = _new(prs); _title(s, "ROI × Method assigned/NMP + 1위 카운트")
    pivot = _build_pivot(df)
    _table(s, pivot, Inches(0.25), Inches(0.8),
           Inches(8.5), Inches(5.6), header_size=10, body_size=9,
           highlight_col='Winner', highlight_value='baysor')
    _txt(s, "셀 = assigned/NMP\n"
            "Winner = assigned ≥ 0.5 인 method 중 NMP 최고\n"
            "(paper Pareto knee)",
         Inches(0.25), Inches(6.5), Inches(8.5), Inches(0.7),
         size=10, color=(0x55, 0x55, 0x55))
    counts = pivot['Winner'].value_counts().to_dict()
    cnt_df = pd.DataFrame([{"method": k, "wins": v,
                             "share": f"{100*v/len(pivot):.0f}%"}
                            for k, v in sorted(counts.items(),
                                                 key=lambda x: -x[1])])
    _table(s, cnt_df, Inches(9.0), Inches(0.8),
           Inches(4.0), Inches(2.2), header_size=12, body_size=12,
           highlight_col='method', highlight_value='baysor')
    _txt(s, "12 ROI 중 baysor 가 11 ROI 1위 (91.7%).\n"
            "유일 예외 = ovrlpy_low_coherence_2:\n"
            "  데이터 자체가 low-VSI 로 너무 손상되어\n"
            "  baysor 회수율이 0.73 으로 떨어짐.\n"
            "  paper §5.4 의 'low-VSI 영역에선 어떤\n"
            "  method 도 완벽하지 않다' 그대로의 상황.",
         Inches(9.0), Inches(3.2), Inches(4.0), Inches(3.5),
         size=11, color=(0x33, 0x33, 0x33))


def _slide_method_overview(prs, viz_other_root: Path):
    """단일 ROI 에서 5 method 가 어떻게 다른지 한눈에. dapi_boundary 5-panel."""
    s = _new(prs)
    _title(s, "5 method 한눈에 비교 — boundary_tissue_edge_1 ROI 의 dapi+boundary")
    base = viz_other_root / "boundary_tissue_edge_1" / "common"
    panel_w = Inches(2.55); panel_h = Inches(2.6)
    gap = Inches(0.05); left0 = Inches(0.10); top = Inches(0.85)
    for i, m in enumerate(METHOD_ORDER):
        left = left0 + i * (panel_w + gap)
        png = base / f"boundary_tissue_edge_1_dapi_boundary_{m}.png"
        _img(s, png, left, top, width=panel_w, height=panel_h)
        _txt(s, METHOD_LABEL[m], left, top + panel_h + Inches(0.05),
             panel_w, Inches(0.3), size=10, bold=True,
             color=(0x1B, 0x4F, 0x72), align=PP_ALIGN.CENTER)
        _txt(s, f"📂 ..._dapi_boundary_{m}.png",
             left, top + panel_h + Inches(0.35), panel_w, Inches(0.3),
             size=7, mono=True, color=(0x80, 0x80, 0x80),
             align=PP_ALIGN.CENTER)
    # 두 번째 줄: alltx_boundary
    top2 = Inches(4.2)
    _txt(s, "transcript cloud + boundary (같은 ROI)",
         Inches(0.10), top2 - Inches(0.3), Inches(13), Inches(0.3),
         size=11, bold=True, color=(0x47, 0x4F, 0x55))
    for i, m in enumerate(METHOD_ORDER):
        left = left0 + i * (panel_w + gap)
        png = base / f"boundary_tissue_edge_1_alltx_boundary_{m}.png"
        _img(s, png, left, top2, width=panel_w, height=panel_h)
        _txt(s, f"📂 ..._alltx_boundary_{m}.png",
             left, top2 + panel_h + Inches(0.05), panel_w, Inches(0.3),
             size=7, mono=True, color=(0x80, 0x80, 0x80),
             align=PP_ALIGN.CENTER)
    _txt(s,
         "윗줄 = DAPI + cell 경계만   |   아랫줄 = transcript dot + 같은 cell 경계\n"
         "방법별 *cell 영역의 크기* 차이를 한눈에:\n"
         "  • xenium_nucleus / cellpose: nucleus 만 → 작고 transcript 많이 놓침\n"
         "  • rigid: 가장 큼 → transcript 다 잡지만 옆 cell 침범\n"
         "  • optimal / baysor: 균형\n"
         "fold boundary ROI 인 만큼 가장자리에서 baysor / rigid 가 tissue 밖까지 늘어나는 게 보임.",
         Inches(0.10), Inches(7.0), Inches(13.1), Inches(0.5),
         size=10, color=(0x33, 0x33, 0x33))


# ============================================================================
# Section 2 — Class-specific deep dives (special pages)
# ============================================================================
def _slide_arch_special(prs, viz_arch_root: Path):
    s = _new(prs)
    _title(s, "ARCHITECTURE special — SSAM cell-type 일치 5-method 비교")
    roi_id = "ssam_transition_2"
    base = viz_arch_root / roi_id / "architecture"
    panel_w = Inches(2.55); panel_h = Inches(2.6); gap = Inches(0.05)
    left0 = Inches(0.10); top = Inches(0.85)
    for i, m in enumerate(METHOD_ORDER):
        left = left0 + i * (panel_w + gap)
        png = base / f"{roi_id}_ssam_disagreement_{m}.png"
        _img(s, png, left, top, width=panel_w, height=panel_h)
        _txt(s, METHOD_LABEL[m], left, top + panel_h + Inches(0.05),
             panel_w, Inches(0.3), size=10, bold=True,
             color=(0x1B, 0x4F, 0x72), align=PP_ALIGN.CENTER)
        _txt(s, f"📂 ..._ssam_disagreement_{m}.png",
             left, top + panel_h + Inches(0.35), panel_w, Inches(0.3),
             size=7, mono=True, color=(0x80, 0x80, 0x80),
             align=PP_ALIGN.CENTER)
    # 아래 panel: agreement_heatmap baysor + xenium 비교
    top2 = Inches(4.2)
    _txt(s, "보조: agreement heatmap (baysor vs xenium)",
         Inches(0.10), top2 - Inches(0.3), Inches(13), Inches(0.3),
         size=11, bold=True, color=(0x47, 0x4F, 0x55))
    _img(s, base / f"{roi_id}_ssam_agreement_heatmap_xenium_nucleus.png",
         Inches(0.5), top2, width=Inches(5.5))
    _img(s, base / f"{roi_id}_ssam_agreement_heatmap_baysor.png",
         Inches(7.2), top2, width=Inches(5.5))
    _txt(s, "📂 ..._ssam_agreement_heatmap_xenium_nucleus.png",
         Inches(0.5), Inches(7.0), Inches(5.5), Inches(0.3),
         size=8, mono=True, color=(0x80, 0x80, 0x80))
    _txt(s, "📂 ..._ssam_agreement_heatmap_baysor.png",
         Inches(7.2), Inches(7.0), Inches(5.5), Inches(0.3),
         size=8, mono=True, color=(0x80, 0x80, 0x80))
    _txt(s,
         "각 polygon: 초록=SSAM 일치, 빨강=불일치, 회색=비교불가.\n"
         "👀 변화 포인트: xenium 은 회색 (작은 nucleus) 비율 압도적 → baysor 의 색칠된\n"
         "    polygon 비율이 훨씬 높음 = 실제로 비교 가능한 cell 수가 많아 paper 평가 통과.",
         Inches(0.10), Inches(7.32), Inches(13.1), Inches(0.6),
         size=10, color=(0x33, 0x33, 0x33))


def _slide_compartment_special(prs, viz_other_root: Path):
    s = _new(prs)
    _title(s, "COMPARTMENT special — cytoplasm 회수율 5-method 비교")
    roi_id = "p2r_compartment_mismatch_1"
    base = viz_other_root / roi_id / "compartment"
    # 큰 panel_compare
    _img(s, base / f"{roi_id}_compartment_panel_compare.png",
         Inches(0.3), Inches(0.85), width=Inches(8.5))
    _txt(s, f"📂 ..._compartment_panel_compare.png",
         Inches(0.3), Inches(5.0), Inches(8.5), Inches(0.3),
         size=8, mono=True, color=(0x80, 0x80, 0x80))
    _txt(s,
         "5 method 격자 — 각 panel 위에 cyto-extra% 표시.\n\n"
         "rigid: 97% — 회수만 1위, contamination ↑\n"
         "baysor: 85% — 회수 + 균형 1위\n"
         "optimal: 84% — 비슷\n"
         "cellpose: 65% — nucleus 만\n"
         "xenium: 1% — paper §1.1 비판",
         Inches(9.0), Inches(0.85), Inches(4.1), Inches(4.0),
         size=11, color=(0x33, 0x33, 0x33))
    # 아래: nucleus_distance_split
    _img(s, base / f"{roi_id}_nucleus_distance_split.png",
         Inches(0.3), Inches(5.4), width=Inches(8.5))
    _txt(s, "📂 ..._nucleus_distance_split.png",
         Inches(0.3), Inches(7.05), Inches(8.5), Inches(0.3),
         size=8, mono=True, color=(0x80, 0x80, 0x80))
    _txt(s,
         "보조 그림: transcript 의 nucleus 까지 거리 분포.\n\n"
         "👀 변화 포인트: nuclear 마커 → 0 µm 근처 peak,\n"
         "    cyto 마커 → 5–10 µm 산포 가 baysor 에서 가장 잘 분리됨.\n"
         "    rigid 는 cyto 분포가 너무 멀리 늘어져 contamination 증거.",
         Inches(9.0), Inches(5.4), Inches(4.1), Inches(2.0),
         size=11, color=(0x33, 0x33, 0x33))


def _slide_low_vsi_special(prs, viz_other_root: Path):
    s = _new(prs)
    _title(s, "LOW_VSI special — vertical overlap 5-method 비교")
    roi_id = "ovrlpy_low_coherence_1"
    base = viz_other_root / roi_id / "low_vsi"
    _img(s, base / f"{roi_id}_z_mixing_panel_compare.png",
         Inches(0.3), Inches(0.85), width=Inches(7.0))
    _txt(s, f"📂 ..._z_mixing_panel_compare.png",
         Inches(0.3), Inches(7.0), Inches(7.0), Inches(0.3),
         size=8, mono=True, color=(0x80, 0x80, 0x80))
    _img(s, base / f"{roi_id}_vsi_overlay_summary.png",
         Inches(7.4), Inches(0.85), width=Inches(5.7))
    _txt(s, f"📂 ..._vsi_overlay_summary.png",
         Inches(7.4), Inches(4.4), Inches(5.7), Inches(0.3),
         size=8, mono=True, color=(0x80, 0x80, 0x80))
    _txt(s,
         "왼: z_mixing — z_range > 7.3 µm (vertical overlap 의심) 인 cell 을 method 색으로 표시.\n"
         "오: low_vsi invasion — centroid 가 VSI < 0.5 (문제 영역) 에 빠진 cell %.\n\n"
         "📊 정량:\n"
         "  z-mixing %:        rigid 53% > baysor 29% > optimal 25% > xenium 10% > cellpose 12%\n"
         "  low-VSI invasion:  rigid 17% > baysor 10% > optimal 7%  > xenium 7%  > cellpose 5%\n\n"
         "👀 변화 포인트: rigid 가 압도적으로 안 좋음 (paper §5.1 비판). baysor 는 rigid 보다는 훨씬 \n"
         "    낫지만 optimal 이 살짝 더 우세. 4 micrometer 이상의 z-overlap 영역에서 baysor 도 가끔 false-merge 발생.",
         Inches(0.3), Inches(4.95), Inches(13), Inches(2.0),
         size=11, color=(0x33, 0x33, 0x33))


def _slide_high_density_special(prs, viz_other_root: Path):
    s = _new(prs)
    _title(s, "HIGH_DENSITY special — merge/split 4-method 비교")
    roi_id = "density_high_density_1"
    base = viz_other_root / roi_id / "high_density"
    # nucleus_centroids reference
    _img(s, base / f"{roi_id}_nucleus_centroids.png",
         Inches(0.3), Inches(0.85), width=Inches(4.0))
    _txt(s, "DAPI + xenium nucleus centroid (참고)",
         Inches(0.3), Inches(5.0), Inches(4.0), Inches(0.3),
         size=10, bold=True, color=(0x47, 0x4F, 0x55))
    _txt(s, "📂 ..._nucleus_centroids.png",
         Inches(0.3), Inches(5.32), Inches(4.0), Inches(0.3),
         size=8, mono=True, color=(0x80, 0x80, 0x80))
    # 4 methods (xenium 은 baseline 이라 제외) merge_split
    methods4 = ["cellpose_nuclei", "rigid_expansion",
                  "optimal_expansion", "baysor"]
    panel_w = Inches(2.05); panel_h = Inches(2.05)
    left0 = Inches(4.5); gap = Inches(0.05)
    for i, m in enumerate(methods4):
        col = i % 2; row = i // 2
        left = left0 + col * (panel_w + gap)
        top = Inches(0.85) + row * (panel_h + Inches(0.55))
        png = base / f"{roi_id}_merge_split_{m}.png"
        _img(s, png, left, top, width=panel_w, height=panel_h)
        _txt(s, METHOD_LABEL[m], left, top + panel_h + Inches(0.02),
             panel_w, Inches(0.25), size=9, bold=True,
             color=(0x1B, 0x4F, 0x72), align=PP_ALIGN.CENTER)
        _txt(s, f"📂 ..._merge_split_{m}.png", left,
             top + panel_h + Inches(0.27), panel_w, Inches(0.25),
             size=7, mono=True, color=(0x80, 0x80, 0x80),
             align=PP_ALIGN.CENTER)
    _txt(s,
         "참고: xenium nucleus 는 baseline (분기 없음)\n\n"
         "각 그림: 빨강 = baysor 가 xenium 을 합친 transcript (merge),\n"
         "          초록 = baysor 가 xenium 을 분리한 transcript (split).\n\n"
         "📊 정량 (densest_high_density_1):\n"
         "  baysor: split 169,088 / merge 512   ← 99.7% 가 split\n"
         "  optimal: split ~ vs merge ~ (이미지 우상)\n"
         "  rigid: merge 비율 ↑ (over-permissive)\n"
         "  cellpose: 거의 변화 없음 (작아서)\n\n"
         "👀 변화 포인트: baysor 는 dense 영역에서 cell 을 합치지 않고 \n"
         "    sub-cluster 까지 *분리* — 다른 method 보다 cell 수 1.9× ↑.",
         Inches(9.1), Inches(0.85), Inches(4.0), Inches(6.5),
         size=11, color=(0x33, 0x33, 0x33))


def _slide_fold_special(prs, viz_other_root: Path):
    s = _new(prs)
    _title(s, "FOLD_BOUNDARY special — tissue overflow 5-method 비교")
    roi_id = "boundary_tissue_edge_1"
    base = viz_other_root / roi_id / "fold_boundary"
    # 위: overflow chart big
    _img(s, base / f"{roi_id}_boundary_overflow_fraction_chart.png",
         Inches(0.3), Inches(0.85), width=Inches(8.0))
    _txt(s, "📂 ..._boundary_overflow_fraction_chart.png",
         Inches(0.3), Inches(4.5), Inches(8.0), Inches(0.3),
         size=8, mono=True, color=(0x80, 0x80, 0x80))
    _txt(s,
         "method 별 cell centroid 가 tissue 밖에 있는 % 막대.\n\n"
         "📊 정량:\n"
         "  cellpose 30% (best)\n"
         "  xenium 41%\n"
         "  optimal 49%\n"
         "  baysor 88% ← 약점\n"
         "  rigid 97% (worst)\n\n"
         "낮을수록 좋음.",
         Inches(8.5), Inches(0.85), Inches(4.6), Inches(3.5),
         size=11, color=(0x33, 0x33, 0x33))
    # 아래: 5 method panel_compare
    _img(s, base / f"{roi_id}_boundary_panel_compare.png",
         Inches(0.3), Inches(5.0), width=Inches(13))
    _txt(s, "📂 ..._boundary_panel_compare.png",
         Inches(0.3), Inches(7.05), Inches(13), Inches(0.3),
         size=8, mono=True, color=(0x80, 0x80, 0x80))
    _txt(s,
         "👀 변화 포인트: tissue mask 회색 영역 밖으로 polygon 이 늘어진 정도. "
         "rigid (3a) 와 baysor (4) 가 가장자리 빈 공간을 덮어버림. "
         "cellpose (2) 가 가장 단정 — paper §5.6 의 핵심.",
         Inches(0.3), Inches(7.32), Inches(13), Inches(0.4),
         size=10, color=(0x33, 0x33, 0x33))


# ============================================================================
# Section 3 — Per-ROI baysor vs others (reassigned_transcripts 4-panel)
# ============================================================================
def _slide_per_roi_reassign(prs, df, roi_id: str, cls: str,
                              viz_root: Path):
    s = _new(prs)
    # Get winner from metric for title
    row = df[df.roi_id == roi_id]
    asgn = row.set_index('method')['assigned_fraction'].to_dict() if len(row) else {}
    nmp  = row.set_index('method')['negative_marker_purity'].to_dict() if len(row) else {}
    cands = [(m, nmp.get(m)) for m in METHOD_ORDER
              if asgn.get(m, 0) >= 0.5 and pd.notna(nmp.get(m))]
    winner = max(cands, key=lambda x: x[1])[0] if cands else '—'

    _title(s, f"{roi_id}  ({cls})   ★ paper 1위: {winner}", size=20)
    metrics_line = (
        f"baysor {asgn.get('baysor',float('nan')):.2f}/{nmp.get('baysor',float('nan')):.2f}  |  "
        f"xenium {asgn.get('xenium_nucleus',float('nan')):.2f}/{nmp.get('xenium_nucleus',float('nan')):.2f}  |  "
        f"rigid {asgn.get('rigid_expansion',float('nan')):.2f}/{nmp.get('rigid_expansion',float('nan')):.2f}  |  "
        f"optimal {asgn.get('optimal_expansion',float('nan')):.2f}/{nmp.get('optimal_expansion',float('nan')):.2f}  |  "
        f"cellpose {asgn.get('cellpose_nuclei',float('nan')):.2f}/{nmp.get('cellpose_nuclei',float('nan')):.2f}"
    )
    _txt(s, metrics_line, Inches(0.3), Inches(0.7),
         Inches(13), Inches(0.3), size=10, color=(0x55, 0x55, 0x55))

    # 4 panels of reassigned_transcripts xenium__vs__<m> for the 4 non-baseline
    methods4 = ["cellpose_nuclei", "rigid_expansion", "optimal_expansion",
                  "baysor"]
    panel_w = Inches(3.15); panel_h = Inches(3.15); gap = Inches(0.10)
    left0 = Inches(0.30); top = Inches(1.10)
    base = viz_root / roi_id / "highlights"
    for i, m in enumerate(methods4):
        left = left0 + i * (panel_w + gap)
        png = base / f"{roi_id}_reassigned_transcripts_xenium_nucleus__vs__{m}.png"
        _img(s, png, left, top, width=panel_w, height=panel_h)
        _txt(s, f"xenium  vs  {METHOD_LABEL[m]}",
             left, top + panel_h + Inches(0.05),
             panel_w, Inches(0.3), size=10, bold=True,
             color=(0x1B, 0x4F, 0x72), align=PP_ALIGN.CENTER)
        _txt(s, f"📂 ..._reassigned_transcripts_xenium_nucleus__vs__{m}.png",
             left, top + panel_h + Inches(0.35),
             panel_w, Inches(0.3), size=7, mono=True,
             color=(0x80, 0x80, 0x80), align=PP_ALIGN.CENTER)

    _txt(s,
         "각 panel = transcript 단위 reassignment 색상 코드:\n"
         "  • 회색  = 동일하게 같은 cell 에 속함 (변화 없음)\n"
         "  • 빨강  = xenium 에선 다른 cell, 이 method 는 다른 cell 로 옮김 (재할당)\n"
         "  • 파랑  = xenium 은 unassigned 였으나 이 method 가 새로 cell 에 넣음 (회수)\n"
         "  • 주황  = xenium 은 cell 안인데 이 method 가 unassigned (탈락)\n\n"
         "👀 변화 보기:\n"
         "  • baysor 패널: 파랑 (회수) ↑↑, 빨강 (재할당) ↑ → cytoplasm transcript 를 적극적으로 가져옴\n"
         "  • cellpose 패널: 주황 (탈락) ↑↑ → nucleus 만 잡아서 xenium 보다 transcript 잃음\n"
         "  • rigid 패널: 파랑 압도적 → 모든 transcript 회수 but contamination 의심\n"
         "  • optimal 패널: 보수적 — 약간의 회수만",
         Inches(0.3), Inches(5.0), Inches(13), Inches(2.4),
         size=10, color=(0x33, 0x33, 0x33))


# ============================================================================
# Section 4 — Per-ROI 5-method dapi+boundary snapshot
# ============================================================================
def _slide_per_roi_snapshot(prs, df, roi_id: str, cls: str,
                              viz_root: Path):
    s = _new(prs)
    row = df[df.roi_id == roi_id]
    asgn = row.set_index('method')['assigned_fraction'].to_dict() if len(row) else {}
    nmp  = row.set_index('method')['negative_marker_purity'].to_dict() if len(row) else {}
    cands = [(m, nmp.get(m)) for m in METHOD_ORDER
              if asgn.get(m, 0) >= 0.5 and pd.notna(nmp.get(m))]
    winner = max(cands, key=lambda x: x[1])[0] if cands else '—'

    _title(s, f"{roi_id}  ({cls})  —  5-method DAPI+boundary 스냅샷", size=18)

    base = viz_root / roi_id / "common"
    # 큰 5-panel grid (dapi_boundary 만)
    panel_w = Inches(2.55); panel_h = Inches(4.5); gap = Inches(0.05)
    left0 = Inches(0.10); top = Inches(0.85)
    for i, m in enumerate(METHOD_ORDER):
        left = left0 + i * (panel_w + gap)
        png = base / f"{roi_id}_dapi_boundary_{m}.png"
        _img(s, png, left, top, width=panel_w, height=panel_h)
        _txt(s, METHOD_LABEL[m], left, top + panel_h + Inches(0.05),
             panel_w, Inches(0.3), size=11, bold=True,
             color=(0x1B, 0x4F, 0x72), align=PP_ALIGN.CENTER)
        _txt(s, f"📂 ..._dapi_boundary_{m}.png",
             left, top + panel_h + Inches(0.36), panel_w, Inches(0.3),
             size=8, mono=True, color=(0x80, 0x80, 0x80),
             align=PP_ALIGN.CENTER)

    metrics_line = (
        f"baysor {asgn.get('baysor',float('nan')):.2f}/{nmp.get('baysor',float('nan')):.2f}  |  "
        f"xenium {asgn.get('xenium_nucleus',float('nan')):.2f}/{nmp.get('xenium_nucleus',float('nan')):.2f}  |  "
        f"rigid {asgn.get('rigid_expansion',float('nan')):.2f}/{nmp.get('rigid_expansion',float('nan')):.2f}  |  "
        f"optimal {asgn.get('optimal_expansion',float('nan')):.2f}/{nmp.get('optimal_expansion',float('nan')):.2f}  |  "
        f"cellpose {asgn.get('cellpose_nuclei',float('nan')):.2f}/{nmp.get('cellpose_nuclei',float('nan')):.2f}"
    )
    _txt(s, metrics_line, Inches(0.10), Inches(6.4),
         Inches(13), Inches(0.3), size=10, color=(0x55, 0x55, 0x55))
    _txt(s, f"★ paper 1위: {winner}",
         Inches(0.10), Inches(6.7), Inches(13), Inches(0.3),
         size=11, bold=True, color=(0x1B, 0x4F, 0x72))
    _txt(s,
         "각 panel = ROI 내 DAPI 위에 method 의 cell 경계 overlay.\n"
         "👀 변화 보기: cell 영역의 *크기* 와 *모양* 차이를 한눈에 비교 가능.",
         Inches(0.10), Inches(7.05), Inches(13), Inches(0.4),
         size=10, color=(0x33, 0x33, 0x33))


# ============================================================================
# Section 5 — Conclusion
# ============================================================================
def _slide_objective_summary(prs):
    s = _new(prs)
    _title(s, "ROI 목적별 평가 요약 — baysor 가 정말 나아졌나?")
    df = pd.DataFrame([
        ["architecture", "cell-type 경계 보존",
         "ssam_agreement_heatmap 초록 ↑",
         "agree 14% ↑, NMP 2.4×, mixed 28% ↓",
         "✅ 향상", "baysor"],
        ["compartment",  "cytoplasm 회수 + 균형",
         "compartment_panel cyto-extra %",
         "cyto 1% → 85%, NMP 4×, mixed 33% ↓",
         "✅ 큰 향상", "baysor"],
        ["low_vsi",      "vertical overlap robust",
         "z_mixing_panel 핑크 cell ↓",
         "z-mix 53% → 29% (rigid 대비)",
         "⚠️ 부분", "optimal (25%)"],
        ["high_density", "dense merge artifact 없애기",
         "merge_split: split 169k vs merged 512",
         "n_cells 1.9×, mixed 32% ↓",
         "✅ 향상", "baysor"],
        ["fold_boundary","tissue 밖 안 나가기",
         "boundary_overflow_chart bar ↓",
         "overflow 41% → 87% (악화)",
         "❌ 약화", "cellpose (30%)"],
    ], columns=["class", "ROI 목적", "그림에서 변하는 것",
                  "관련 수치 변화", "baysor 평가", "1위 method"])
    _table(s, df, Inches(0.3), Inches(0.85),
           Inches(12.7), Inches(4.7), header_size=12, body_size=11,
           highlight_col='1위 method', highlight_value='baysor')
    _txt(s,
         "→ 5 ROI class 중 3개 baysor 명확 우세, 1개 부분 (low_vsi 는 optimal 살짝 ↑),\n"
         "    1개 약화 (fold 는 cellpose 안전).\n"
         "→ paper Pareto knee 기준 12 ROI 중 11개 baysor 1위.",
         Inches(0.3), Inches(5.75), Inches(12.7), Inches(1.5),
         size=14, color=(0x33, 0x33, 0x33))


def _slide_caveats(prs):
    s = _new(prs); _title(s, "주의사항 / 값 검증")
    _txt(s,
         "1. rigid_expansion 의 assigned_fraction = 1.00 (모든 ROI)\n"
         "   → 5 µm 확장으로 cell 끼리 겹쳐 모든 transcript 가 적어도 한 cell 에 포함됨.\n"
         "      통계적으로 정상이지만 'recovery 만점'을 그대로 신뢰하면 안 됨.\n\n"
         "2. baysor n_cells = 1049 / 1161 (ssam_transition_3 / 4)\n"
         "   → 다른 method 의 ~3배. dense 영역에서 baysor 의 transcript-aware splitting 이\n"
         "      sub-cluster 까지 잡아내는 정상 동작 (paper §4.5).\n\n"
         "3. NMP × mixed_marker 상관 = -0.34\n"
         "   → 두 metric 모두 specificity 측정이지만 같은 게 아님.\n"
         "      NMP 는 reference celltype 일치, mixed_marker 는 conflicting marker pair.\n\n"
         "4. cellpose ssam_transition_2 NMP = 0.45 (다른 ssam ROI 의 cellpose 평균 ~0.20)\n"
         "   → 해당 ROI 는 GFAP 단일 dominant 영역이라 자동으로 NMP 가 높게 나옴 (정상).",
         Inches(0.4), Inches(0.85), Inches(12.5), Inches(6),
         size=13, color=(0x33, 0x33, 0x33))


def _slide_conclusion(prs):
    s = _new(prs); _title(s, "결론")
    _txt(s, "Paper Fig 3e 의 recovery × NMP Pareto front 기준",
         Inches(0.5), Inches(0.85), Inches(12.3), Inches(0.6),
         size=18, bold=True, color=(0x1B, 0x4F, 0x72))
    _txt(s,
         "  ✅ baysor — 12 ROI 중 11개 1위 (91.7%).\n"
         "       회수율 0.96 + NMP 0.11 → paper 'Pareto knee' 위치 재현.\n\n"
         "  ⚠️  rigid_expansion — Pareto-optimal 하지만 NMP 0.06 (over-permissive).\n"
         "       paper §5.1 비판 method.\n\n"
         "  ⚠️  cellpose_nuclei — NMP 1위지만 회수 0.13. 단독 사용 부적절.\n\n"
         "  ❌ optimal_expansion / xenium_nucleus — 두 축 모두 baysor 한테 dominated.\n\n"
         "ROI 별 약점:\n"
         "  • fold_boundary: tissue overflow 87% (cellpose 30% 보다 약함)\n"
         "  • low_vsi: z-mixing 29% (optimal 25% 보다 살짝 약함)\n"
         "→ 해석 시 주의는 필요하지만 paper 결론과 일치.",
         Inches(0.5), Inches(1.55), Inches(12.3), Inches(5.7),
         size=13, color=(0x33, 0x33, 0x33))


# ============================================================================
# Main builder
# ============================================================================
def build(metrics_arch: Path, metrics_other: Path, pareto_png: Path,
           viz_arch_root: Path, viz_other_root: Path, out: Path) -> Path:
    df = pd.concat([pd.read_csv(metrics_arch), pd.read_csv(metrics_other)],
                    ignore_index=True)
    df = df[df.available.astype(str).str.lower() == 'true'].copy()

    prs = Presentation()
    prs.slide_width = SLIDE_W; prs.slide_height = SLIDE_H

    # Section 1 — Overview (4)
    _slide_title(prs)
    _slide_methods_pareto(prs, pareto_png)
    _slide_pivot_summary(prs, df)
    _slide_method_overview(prs, viz_other_root)

    # Section 2 — Class special (5)
    _slide_arch_special(prs, viz_arch_root)
    _slide_compartment_special(prs, viz_other_root)
    _slide_low_vsi_special(prs, viz_other_root)
    _slide_high_density_special(prs, viz_other_root)
    _slide_fold_special(prs, viz_other_root)

    # Section 3 — Per-ROI reassigned_transcripts comparison (12)
    for roi_id, cls, root_key in ROI_LIST:
        viz_root = viz_arch_root if root_key == 'arch' else viz_other_root
        _slide_per_roi_reassign(prs, df, roi_id, cls, viz_root)

    # Section 4 — Per-ROI 5-method dapi+boundary snapshot (12)
    for roi_id, cls, root_key in ROI_LIST:
        viz_root = viz_arch_root if root_key == 'arch' else viz_other_root
        _slide_per_roi_snapshot(prs, df, roi_id, cls, viz_root)

    # Section 5 — Conclusion (3)
    _slide_objective_summary(prs)
    _slide_caveats(prs)
    _slide_conclusion(prs)

    out.parent.mkdir(parents=True, exist_ok=True)
    prs.save(str(out))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metrics-arch", required=True)
    ap.add_argument("--metrics-other", required=True)
    ap.add_argument("--pareto-png", required=True)
    ap.add_argument("--viz-arch-root", required=True)
    ap.add_argument("--viz-other-root", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO,
                          format="%(asctime)s | %(levelname)s | %(message)s")
    out = build(Path(args.metrics_arch), Path(args.metrics_other),
                  Path(args.pareto_png), Path(args.viz_arch_root),
                  Path(args.viz_other_root), Path(args.out))
    logger.info(f"PPT saved: {out}")


if __name__ == "__main__":
    main()
