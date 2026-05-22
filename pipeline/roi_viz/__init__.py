"""ROI resegmentation comparison visualization package.

Modular pipeline that, for each Step 2 ROI, generates:
- Common images (context, DAPI crop, reference crop, metric card, per-method
  DAPI+boundary / all-tx+boundary / markers+boundary / celltype map).
- Difference highlight images (boundary overlap, difference mask, transcript
  reassignment) between a baseline (xenium_nucleus) and each compared method.
- ROI-class-specific add-on images (easy_control / architecture / compartment /
  low_vsi / high_density / fold_boundary).

Built around an explicit cache layer so that partial edits do not trigger
re-computation of the whole pipeline.
"""

__version__ = "0.2.0"
