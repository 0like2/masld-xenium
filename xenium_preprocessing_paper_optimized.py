#!/usr/bin/env python3
"""
Xenium Data Preprocessing Pipeline - Paper Optimized Version (Option C)
========================================================================

Pure paper-based optimization: Fixed parameters based on Notebook 1_6 benchmarking.

This is the SIMPLEST, MOST REPRODUCIBLE version with parameters fixed to
the values that showed best performance in the optimization study.

Key Features:
1. target_sum=100 (Xenium-optimized, NOT 1e4)
2. min_counts=40, min_genes=15 (corrected filtering)
3. HVG selection with optimized parameters (min/mean/disp)
4. NO scaling (shown to be detrimental)
5. All PCA components (n_pcs=0)
6. Single primary Leiden clustering (resolution=1.0)
7. min_dist=0.1 for UMAP (tighter clusters)
8. No extra complexity - just the essentials

References:
    - Notebook 1_2: Cell type identification basics
    - Notebook 1_6: Batch preprocessing optimization (MAIN SOURCE)
    - Notebook 6_3: Simulated data validation

Usage:
    # Basic usage
    python xenium_preprocessing_paper_optimized.py --input data.h5ad --output preprocessed.h5ad

    # Generate and preprocess simulated data
    python xenium_preprocessing_paper_optimized.py --simulate --output simulated.h5ad

    # From Xenium directory
    python xenium_preprocessing_paper_optimized.py --xenium-dir /path/to/xenium --output out.h5ad
"""

import argparse
import json
import warnings
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore", category=UserWarning)

# Set random seed for reproducibility
np.random.seed(0)

# Configure scanpy
sc.settings.set_figure_params(dpi=120, facecolor="white")
sc.settings.verbosity = 2


class XeniumPreprocessorPaperOptimized:
    """
    Paper-optimized preprocessing class for Xenium data.

    All parameters are fixed based on Notebook 1_6 benchmarking results.
    No customization needed - this is the recommended pipeline.
    """

    # FIXED PARAMETERS (from Notebook 1_6 optimization)
    TARGET_SUM = 100  # Xenium-optimized (NOT 1e4)
    MIN_COUNTS = 40
    MIN_GENES = 15
    MIN_CELLS = 30

    # HVG parameters (from Notebook 1_6)
    HVG_MIN_MEAN = 0.0125
    HVG_MAX_MEAN = 3
    HVG_MIN_DISP = 0.5

    # Dimensionality reduction (from Notebook 1_6)
    N_NEIGHBORS = 15
    N_PCS = 0  # Use all components
    UMAP_MIN_DIST = 0.1
    UMAP_SPREAD = 1.0

    # Clustering (primary resolution)
    LEIDEN_RESOLUTION = 1.0

    def __init__(
        self,
        input_file: Optional[str] = None,
        output_dir: str = "./output",
        library_id: str = "xenium_sample",
    ):
        """Initialize preprocessor with fixed parameters."""
        self.input_file = input_file
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.library_id = library_id
        self.adata = None
        self.log = []

    def log_message(self, msg: str) -> None:
        """Log and print message."""
        print(f"[INFO] {msg}")
        self.log.append(msg)

    def load_data(self) -> "XeniumPreprocessorPaperOptimized":
        """Load Xenium data from H5AD file."""
        self.log_message(f"Loading data from {self.input_file}...")

        if self.input_file.endswith(".h5ad"):
            self.adata = sc.read_h5ad(self.input_file)
        else:
            raise ValueError("Input file must be .h5ad format")

        self.log_message(f"Loaded: {self.adata.shape[0]} cells × {self.adata.shape[1]} genes")
        return self

    def load_xenium_directory(self, base_dir: str) -> "XeniumPreprocessorPaperOptimized":
        """Load from Xenium output directory."""
        self.log_message(f"Loading Xenium data from {base_dir}...")

        base_dir = Path(base_dir)
        counts_path = base_dir / "cell_feature_matrix"

        try:
            # Load count matrix
            self.adata = sc.read_10x_mtx(counts_path, var_names="gene_symbols", make_unique=True)

            # Load cell metadata
            cells = pd.read_parquet(base_dir / "cells.parquet").set_index("cell_id")
            common = self.adata.obs_names.intersection(cells.index)
            self.adata = self.adata[common].copy()
            self.adata.obs = cells.loc[common]

            # Add spatial coordinates
            self.adata.obsm["spatial"] = self.adata.obs[["x_centroid", "y_centroid"]].to_numpy()

            # Add library metadata
            self.adata.obs["library_id"] = self.library_id
            self.adata.layers["counts"] = self.adata.X.copy()

            # Setup spatial metadata
            self.adata.uns.setdefault("spatial", {})[self.library_id] = {
                "metadata": {"coordinate_system": "xenium"},
                "scalefactors": {"spot_diameter_fullres": 1.0},
            }

            self.log_message(f"Loaded: {self.adata.shape[0]} cells × {self.adata.shape[1]} genes")
            return self

        except Exception as e:
            raise RuntimeError(f"Failed to load Xenium directory: {e}")

    def quality_control(self) -> "XeniumPreprocessorPaperOptimized":
        """Compute QC metrics and filter cells/genes."""
        self.log_message("Computing QC metrics...")

        # Calculate QC metrics
        sc.pp.calculate_qc_metrics(self.adata, inplace=True)

        # Print QC statistics
        qc_summary = self.adata.obs[["total_counts", "n_genes_by_counts"]].describe(
            percentiles=[0.5, 0.9, 0.99]
        )
        self.log_message(f"\nQC Summary:\n{qc_summary}")

        # Filter cells using FIXED thresholds
        max_counts = float(np.quantile(self.adata.obs["total_counts"], 0.99) * 4)

        self.log_message(
            f"Filtering: min_counts={self.MIN_COUNTS}, "
            f"min_genes={self.MIN_GENES}, max_counts={max_counts}"
        )

        # Filter cells
        n_before = self.adata.n_obs
        self.adata = self.adata[
            (self.adata.obs["total_counts"] >= self.MIN_COUNTS)
            & (self.adata.obs["total_counts"] <= max_counts)
            & (self.adata.obs["n_genes_by_counts"] >= self.MIN_GENES)
        ].copy()
        n_after = self.adata.n_obs
        self.log_message(f"Filtered cells: {n_before} → {n_after} ({100*n_after/n_before:.1f}% retained)")

        # Filter genes
        n_genes_before = self.adata.n_vars
        sc.pp.filter_genes(self.adata, min_cells=self.MIN_CELLS)
        n_genes_after = self.adata.n_vars
        self.log_message(
            f"Filtered genes: {n_genes_before} → {n_genes_after} "
            f"({100*n_genes_after/n_genes_before:.1f}% retained)"
        )

        # Store raw
        self.adata.raw = self.adata.copy()
        self.adata.layers["raw"] = self.adata.X.copy()

        return self

    def normalize_and_hvg(self) -> "XeniumPreprocessorPaperOptimized":
        """Normalize, log-transform, and select HVGs (FIXED PARAMETERS)."""
        self.log_message(f"Normalizing with target_sum={self.TARGET_SUM}...")

        # Ensure counts are in layers
        if "counts" not in self.adata.layers:
            self.adata.layers["counts"] = self.adata.X.copy()

        # FIXED: target_sum=100 (Xenium optimized)
        sc.pp.normalize_total(self.adata, target_sum=self.TARGET_SUM)

        # Log transform
        self.log_message("Log-transforming...")
        sc.pp.log1p(self.adata)

        # Select HVGs with FIXED parameters
        self.log_message(
            f"Selecting HVGs (min_mean={self.HVG_MIN_MEAN}, "
            f"max_mean={self.HVG_MAX_MEAN}, min_disp={self.HVG_MIN_DISP})..."
        )
        sc.pp.highly_variable_genes(
            self.adata,
            min_mean=self.HVG_MIN_MEAN,
            max_mean=self.HVG_MAX_MEAN,
            min_disp=self.HVG_MIN_DISP,
            subset=True,
        )
        self.log_message(f"Selected genes: {self.adata.shape[1]}")

        # NOTE: Scaling NOT applied (Notebook 1_6 showed it's detrimental)

        return self

    def dimensionality_reduction(self) -> "XeniumPreprocessorPaperOptimized":
        """PCA and UMAP (FIXED PARAMETERS)."""
        self.log_message("Computing PCA (all components)...")
        sc.tl.pca(self.adata)

        self.log_message(
            f"Computing neighbors (n_neighbors={self.N_NEIGHBORS}, n_pcs={self.N_PCS})..."
        )
        sc.pp.neighbors(
            self.adata,
            n_neighbors=self.N_NEIGHBORS,
            n_pcs=self.N_PCS,
            metric="cosine"
        )

        self.log_message(f"Computing UMAP (min_dist={self.UMAP_MIN_DIST})...")
        sc.tl.umap(self.adata, min_dist=self.UMAP_MIN_DIST, spread=self.UMAP_SPREAD)

        return self

    def clustering(self) -> "XeniumPreprocessorPaperOptimized":
        """Leiden clustering with FIXED resolution."""
        self.log_message(f"Running Leiden clustering (resolution={self.LEIDEN_RESOLUTION})...")
        sc.tl.leiden(self.adata, resolution=self.LEIDEN_RESOLUTION, key_added="leiden")

        n_clusters = self.adata.obs["leiden"].nunique()
        self.log_message(f"Found {n_clusters} clusters")

        return self

    def marker_scoring(self) -> "XeniumPreprocessorPaperOptimized":
        """Score basic cell type markers."""
        marker_dict = {
            "Neuronal": ["SYN1", "RBFOX3", "NEFL"],
            "Glial": ["GFAP", "S100B", "AQP4"],
            "Immune": ["PTPRC", "CD3E", "CD79A"],
            "Endothelial": ["PECAM1", "VWF", "CDH5"],
        }

        self.log_message("Computing cell type marker scores...")

        # Filter to available genes
        marker_dict = {
            k: [g for g in v if g in self.adata.var_names]
            for k, v in marker_dict.items()
        }
        marker_dict = {k: v for k, v in marker_dict.items() if len(v) > 0}

        # Score each cell type
        for label, genes in marker_dict.items():
            sc.tl.score_genes(self.adata, gene_list=genes, score_name=f"score_{label}")

        # Assign cell type hint
        score_cols = [c for c in self.adata.obs.columns if c.startswith("score_")]
        if score_cols:
            score_df = self.adata.obs[score_cols]
            self.adata.obs["celltype_hint"] = score_df.idxmax(axis=1).str.replace("score_", "")
        else:
            self.adata.obs["celltype_hint"] = "Unknown"

        self.log_message(f"Cell type hints: {self.adata.obs['celltype_hint'].value_counts().to_dict()}")

        return self

    def spatial_neighbors(self) -> "XeniumPreprocessorPaperOptimized":
        """Compute spatial neighbor graph."""
        self.log_message("Computing spatial neighbors (Delaunay)...")
        try:
            sq.gr.spatial_neighbors(self.adata, coord_type="generic", delaunay=True)
            self.log_message("Spatial neighbor graph computed")
        except Exception as e:
            self.log_message(f"Warning: Could not compute spatial neighbors: {e}")

        return self

    def rank_genes_groups(self) -> "XeniumPreprocessorPaperOptimized":
        """Rank differentially expressed genes."""
        self.log_message("Ranking genes by cluster (Wilcoxon test)...")
        sc.tl.rank_genes_groups(self.adata, groupby="leiden", method="wilcoxon")
        self.log_message("DEG analysis complete")

        return self

    def save_results(self, output_h5ad: Optional[str] = None) -> str:
        """Save preprocessed data and metadata."""
        if output_h5ad is None:
            output_h5ad = self.output_dir / "xenium_preprocessed_paper_optimized.h5ad"

        self.log_message(f"Saving to {output_h5ad}...")
        self.adata.write_h5ad(output_h5ad, compression="gzip")

        # Save metadata
        metadata = {
            "n_obs": self.adata.n_obs,
            "n_vars": self.adata.n_vars,
            "n_clusters": self.adata.obs["leiden"].nunique(),
            "parameters": {
                "target_sum": self.TARGET_SUM,
                "min_counts": self.MIN_COUNTS,
                "min_genes": self.MIN_GENES,
                "hvg_min_mean": self.HVG_MIN_MEAN,
                "hvg_max_mean": self.HVG_MAX_MEAN,
                "hvg_min_disp": self.HVG_MIN_DISP,
                "n_neighbors": self.N_NEIGHBORS,
                "leiden_resolution": self.LEIDEN_RESOLUTION,
                "umap_min_dist": self.UMAP_MIN_DIST,
            },
            "preprocessing_log": self.log,
        }

        metadata_file = self.output_dir / "metadata.json"
        with open(metadata_file, "w") as f:
            json.dump(metadata, f, indent=2)

        self.log_message(f"Metadata saved to {metadata_file}")

        # Save figures
        self._save_figures()

        return str(output_h5ad)

    def _save_figures(self) -> None:
        """Save QC and analysis figures."""
        self.log_message("Saving figures...")

        # QC metrics
        fig, axes = plt.subplots(1, 2, figsize=(12, 4))
        self.adata.obs["total_counts"].hist(bins=50, ax=axes[0])
        axes[0].set_xlabel("Total counts")
        axes[0].set_ylabel("Frequency")

        self.adata.obs["n_genes_by_counts"].hist(bins=50, ax=axes[1])
        axes[1].set_xlabel("Number of genes")
        axes[1].set_ylabel("Frequency")
        plt.tight_layout()
        plt.savefig(self.output_dir / "01_qc_metrics.png", dpi=120, bbox_inches="tight")
        plt.close()

        # UMAP colored by cluster
        fig = sc.pl.umap(
            self.adata,
            color=["leiden", "total_counts"],
            wspace=0.3,
            ncols=2,
            show=False,
            return_fig=True,
        )
        fig.savefig(self.output_dir / "02_umap_clusters.png", dpi=120, bbox_inches="tight")
        plt.close()

        # Spatial plot colored by cluster
        if "leiden" in self.adata.obs.columns:
            try:
                sq.pl.spatial_scatter(
                    self.adata,
                    color="leiden",
                    size=1,
                    library_key="library_id",
                    spatial_key="spatial",
                    save=str(self.output_dir / "03_spatial_clusters.png"),
                )
            except Exception as e:
                self.log_message(f"Warning: Could not save spatial plot: {e}")

        self.log_message("Figures saved")

    def run_pipeline(self, output_h5ad: Optional[str] = None) -> "XeniumPreprocessorPaperOptimized":
        """Run complete preprocessing pipeline with FIXED parameters."""
        self.log_message("=" * 70)
        self.log_message("XENIUM PREPROCESSING - PAPER OPTIMIZED (OPTION C)")
        self.log_message("=" * 70)
        self.log_message("Pure paper-based optimization with fixed parameters")
        self.log_message("Based on Notebook 1_6 benchmarking study")
        self.log_message("=" * 70)

        self.quality_control()
        self.normalize_and_hvg()
        self.dimensionality_reduction()
        self.clustering()
        self.marker_scoring()
        self.rank_genes_groups()
        self.spatial_neighbors()
        self.save_results(output_h5ad=output_h5ad)

        self.log_message("=" * 70)
        self.log_message("PREPROCESSING COMPLETE!")
        self.log_message("=" * 70)

        return self


def generate_simulated_data(
    n_cells: int = 5000,
    n_genes: int = 2000,
    n_types: int = 5,
    output_file: Optional[str] = None,
) -> object:
    """Generate simulated Xenium-like data for testing."""
    print(f"[INFO] Generating simulated data: {n_cells} cells × {n_genes} genes")

    np.random.seed(0)

    # Generate expression matrix
    X = np.random.poisson(lam=5, size=(n_cells, n_genes)).astype(np.float32)

    # Create AnnData
    adata = sc.AnnData(X=X)
    adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
    adata.var_names = [f"gene_{i}" for i in range(n_genes)]

    # Add cell type labels
    adata.obs["cell_type"] = np.random.choice([f"type_{i}" for i in range(n_types)], n_cells)

    # Add spatial coordinates
    spatial_coords = np.random.uniform(0, 1000, size=(n_cells, 2))
    adata.obsm["spatial"] = spatial_coords
    adata.obs["x_centroid"] = spatial_coords[:, 0]
    adata.obs["y_centroid"] = spatial_coords[:, 1]

    # Add library metadata
    adata.obs["library_id"] = "simulated"
    adata.layers["counts"] = adata.X.copy()
    adata.uns["spatial"] = {
        "simulated": {
            "metadata": {"coordinate_system": "xenium"},
            "scalefactors": {"spot_diameter_fullres": 1.0},
        }
    }

    print(f"[INFO] Generated: {adata.shape[0]} cells × {adata.shape[1]} genes")

    if output_file:
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(output_file, compression="gzip")
        print(f"[INFO] Saved to {output_file}")

    return adata


def main():
    """Main function with argument parsing."""
    parser = argparse.ArgumentParser(
        description="Xenium Preprocessing - Paper Optimized (Fixed Parameters)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic preprocessing with fixed optimal parameters
  python xenium_preprocessing_paper_optimized.py --input data.h5ad --output out.h5ad

  # Generate and preprocess simulated data
  python xenium_preprocessing_paper_optimized.py --simulate --output simulated.h5ad

  # From Xenium directory
  python xenium_preprocessing_paper_optimized.py --xenium-dir /path/to/xenium --output out.h5ad
        """,
    )

    parser.add_argument("--input", help="Input H5AD file")
    parser.add_argument("--xenium-dir", help="Input Xenium directory")
    parser.add_argument("--output", default="preprocessed_paper_optimized.h5ad", help="Output H5AD file")
    parser.add_argument("--output-dir", default="./output", help="Output directory")
    parser.add_argument("--library-id", default="xenium_sample", help="Library ID")
    parser.add_argument("--simulate", action="store_true", help="Generate simulated data")
    parser.add_argument("--sim-cells", type=int, default=5000, help="Simulated cell count")
    parser.add_argument("--sim-genes", type=int, default=2000, help="Simulated gene count")

    args = parser.parse_args()

    # Validate input
    if not (args.input or args.xenium_dir or args.simulate):
        parser.error("Provide --input, --xenium-dir, or --simulate")

    # Generate simulated data if requested
    if args.simulate:
        adata = generate_simulated_data(
            n_cells=args.sim_cells,
            n_genes=args.sim_genes,
            output_file=args.output_dir + "/simulated_raw.h5ad",
        )
        input_file = args.output_dir + "/simulated_raw.h5ad"
    else:
        input_file = args.input or args.xenium_dir

    # Initialize and run preprocessing
    preprocessor = XeniumPreprocessorPaperOptimized(
        input_file=input_file if (args.input or args.simulate) else None,
        output_dir=args.output_dir,
        library_id=args.library_id,
    )

    # Load data
    if args.xenium_dir:
        preprocessor.load_xenium_directory(args.xenium_dir)
    elif args.input or args.simulate:
        preprocessor.load_data()
    else:
        parser.error("No data source provided")

    # Run pipeline
    preprocessor.run_pipeline(output_h5ad=args.output)

    print(f"\n✅ Preprocessing complete! Output saved to {args.output}")


if __name__ == "__main__":
    main()
