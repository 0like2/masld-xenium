import os
import logging
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
from matplotlib.lines import Line2D
import scipy.sparse as sp
from scipy.sparse import eye, vstack, spmatrix, csr_matrix
from scipy.ndimage import zoom
from sklearn.cluster import MiniBatchKMeans
from sklearn.preprocessing import OneHotEncoder, normalize
from sklearn.neighbors import kneighbors_graph, radius_neighbors_graph
from typing import Any, List, Literal, Optional, Union
import warnings

# Configure Logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# ==============================================================================
# Visualization Helpers
# ==============================================================================

def plot_scalebar(ax, x_pos, y_pos, length_value, text='1mm', color='w', edge_color='k', linewidth=3, fontsize=15):
    """
    Adds a scalebar to a matplotlib axis.
    """
    ax.plot([x_pos, x_pos + length_value], [y_pos, y_pos], color=edge_color, linewidth=linewidth+2)
    ax.plot([x_pos, x_pos + length_value], [y_pos, y_pos], color=color, linewidth=linewidth)
    ax.text(x_pos + length_value/2, y_pos - 50, text, color=color, 
            ha='center', va='bottom', fontsize=fontsize, weight='bold',
            path_effects=[PathEffects.withStroke(linewidth=2, foreground=edge_color)])

# ==============================================================================
# Points2Regions Implementation
# ==============================================================================

COLORS_P2R = [
    [0.9019607843137255, 0.09803921568627451, 0.29411764705882354],
    [0.23529411764705882, 0.7058823529411765, 0.29411764705882354],
    [1.0, 0.8823529411764706, 0.09803921568627451],
    [0.2627450980392157, 0.38823529411764707, 0.8470588235294118],
    [0.9607843137254902, 0.5098039215686274, 0.19215686274509805],
    [0.5686274509803921, 0.11764705882352941, 0.7058823529411765],
    [0.27450980392156865, 0.9411764705882353, 0.9411764705882353],
    [0.9411764705882353, 0.19607843137254902, 0.9019607843137255],
    [0.7372549019607844, 0.9647058823529412, 0.047058823529411764],
    [0.9803921568627451, 0.7450980392156863, 0.7450980392156863],
    [0.0, 0.5019607843137255, 0.5019607843137255],
    [0.9019607843137255, 0.7450980392156863, 1.0],
    [0.6039215686274509, 0.38823529411764707, 0.1411764705882353],
    [1.0, 0.9803921568627451, 0.7843137254901961],
    [0.5019607843137255, 0.0, 0.0],
    [0.6666666666666666, 1.0, 0.7647058823529411],
    [0.5019607843137255, 0.5019607843137255, 0.0],
    [1.0, 0.8470588235294118, 0.6941176470588235],
    [0.0, 0.0, 0.4588235294117647],
    [0.5019607843137255, 0.5019607843137255, 0.5019607843137255],
    [1.0, 1.0, 1.0],
    [0.0, 0.0, 0.0],
]
COLORS_P2R = [[int(255 * v) for v in RGB] for RGB in COLORS_P2R]

# Helper Functions for P2R
def connectivity_matrix(xy: np.ndarray, method="knn", k: int = 5, r: Optional[float] = None, include_self: bool = False) -> sp.spmatrix:
    if method == "knn":
        A = kneighbors_graph(xy, k, include_self=include_self).astype('bool')
    else:
        A = radius_neighbors_graph(xy, r, include_self=include_self).astype('bool')
    return A

def attribute_matrix(cat: np.ndarray, unique_cat: Union[np.ndarray, Literal["auto"]] = "auto", return_encoder: bool = False):
    X = np.array(cat).reshape((-1, 1))
    if not isinstance(unique_cat, str):
        unique_cat_list = [np.array(unique_cat)]
    elif unique_cat == "auto":
        unique_cat_list = "auto"
    else:
        raise ValueError("unique_cat must be a numpy array or the string 'auto'.")
    encoder = OneHotEncoder(categories=unique_cat_list, sparse_output=True, handle_unknown="ignore")
    encoder.fit(X)
    y = encoder.transform(X)
    categories = list(encoder.categories_[0])
    if return_encoder:
        return y, categories, encoder
    return y, categories

def spatial_binning_matrix(xy: np.ndarray, box_width: float, return_grid_props: bool = False):
    mi, ma = xy.min(axis=0, keepdims=True), xy.max(axis=0, keepdims=True)
    xys = xy - mi
    grid = (ma - mi).flatten()
    bin_ids = (xys // box_width).astype("int")
    bin_ids = tuple(x for x in bin_ids.T)
    size = (grid // box_width + 1).astype("int")
    size = tuple(x for x in size)
    linear_ind = np.ravel_multi_index(bin_ids, size)
    bin_matrix, linear_unique_bin_ids = attribute_matrix(linear_ind)
    bin_matrix = bin_matrix.T

    if return_grid_props:
        sub_unique_bin_ids = np.unravel_index(linear_unique_bin_ids, size)
        grid_props = dict(
            grid_coords=sub_unique_bin_ids,
            grid_size=size,
            grid_offset=mi.flatten(),
            grid_scale=1.0/box_width
        )
        return bin_matrix, grid_props
    return bin_matrix

def kde_per_label(xy: np.ndarray, features: sp.spmatrix, sigma: float, return_neighbors: bool = False):
    logging.debug("Computing connectivity matrix for KDE...")
    adj = connectivity_matrix(xy, method="radius", r=2.0 * sigma, include_self=True)
    row, col = adj.nonzero()
    d2 = (xy[row,0] - xy[col,0])**2 + (xy[row,1] - xy[col,1])**2
    d2 = np.exp(-d2 / (2 * sigma * sigma))
    aff = sp.csr_matrix((d2, (row, col)), shape=adj.shape, dtype='float32')
    if not return_neighbors:
        return aff @ features
    else:
        return aff @ features, adj

def create_features(xy: np.ndarray, labels: np.ndarray, unique_labels:np.ndarray, sigma: float, bin_width: Union[float, str, None], min_genes_per_bin: int):
    if isinstance(bin_width, str):
        if bin_width == 'auto':
            bin_width = sigma/3.0

    grid_props = {}
    if bin_width is not None:
        B, grid_props = spatial_binning_matrix(xy, box_width=bin_width, return_grid_props=True)
    else:
        B = eye(len(xy))
    B = B.astype('float32')

    if bin_width is not None:
        x = grid_props['grid_coords'][0] / grid_props['grid_scale'] + grid_props['grid_offset'][0]
        y = grid_props['grid_coords'][1] / grid_props['grid_scale'] + grid_props['grid_offset'][1]
        xy = np.vstack((x,y)).T

    attributes, _ = attribute_matrix(labels, unique_labels)
    attributes = attributes.astype('bool')

    features = kde_per_label(xy, B @ attributes, sigma)

    bin_size = features.sum(axis=1).A.flatten()
    good_bins = bin_size >= min_genes_per_bin
    norms_r = 1.0 / features.sum(axis=1)
    norms_r[np.isinf(norms_r)] = .0
    features = features.multiply(norms_r).tocsr()

    return dict(
        features=features,
        xy_bin=xy,
        norms=norms_r.A.flatten(),
        grid_props=grid_props,
        good_bins=good_bins,
        back_map=B.T.nonzero()[1]
    )

def predict_clusters(kmeans_model, features: spmatrix, good_bins: np.ndarray, back_map: np.ndarray):
    clusters = np.zeros(features.shape[0], dtype='int') - 1
    if np.sum(good_bins) > 0:
        clusters[good_bins] = kmeans_model.predict(features[good_bins,:])
    return clusters[back_map], clusters

# Main Points2Regions Function
def points2regions(
        xy: np.ndarray, 
        gene_labels: np.ndarray, 
        sigma: float, 
        n_clusters: int,
        bin_width: Union[float, str, None] = 'auto', 
        min_genes_per_bin:int = 1, 
        groupids: Union[np.ndarray,None] = None, 
        convert_to_geojson: bool = False, seed:int=42, 
        region_name:str="My regions", 
        return_anndata: bool = False
    ) -> Any:
    
    print(f"    [P2R] Running Points2Regions with sigma={sigma}, n_clusters={n_clusters}, bin_width={bin_width}...")
    xy = np.array(xy, dtype="float32")

    if groupids is not None:
        unique_library_id = np.unique(groupids)
        iterdata = [
            (lib_id, (
                xy[groupids==lib_id,:],
                gene_labels[groupids==lib_id]
            )) for lib_id in unique_library_id
        ]
        get_slice = lambda library_id, data: data == library_id
    else:
        iterdata = [('id', (xy, gene_labels))]
        get_slice = lambda library_id, data: np.ones(len(data), dtype='bool')

    unique_genes = np.unique(gene_labels)
    results = {
        library_id : create_features(
            xy_slice,
            labels_slice,
            unique_genes,
            sigma,
            bin_width,
            min_genes_per_bin
        )
        for library_id, (xy_slice, labels_slice) in iterdata
    }

    # Create train features
    features_list = [r['features'][r['good_bins']] for r in results.values()]
    if len(features_list) == 0:
         logging.error("No valid features found for clustering.")
         return None
         
    X_train = vstack(features_list)

    # Train K-Means
    print("    [P2R] Training KMeans...")
    kmeans = MiniBatchKMeans(n_clusters=n_clusters, n_init='auto', random_state=seed)
    kmeans = kmeans.fit(X_train)

    # Predict
    for library_id, result_dict in results.items():
        cluster_per_gene, cluster_per_bin = predict_clusters(kmeans, result_dict['features'], result_dict['good_bins'], result_dict['back_map'])
        results[library_id]['cluster_per_gene'] = cluster_per_gene
        results[library_id]['cluster_per_bin'] = cluster_per_bin

    # Add clusters to dataframe
    clusters_per_gene = np.zeros(len(xy), dtype='int')
    for library_id in results.keys():
        if groupids is not None:
            library_id_slice_ind = get_slice(library_id, groupids)
        else:
            library_id_slice_ind = get_slice(library_id, xy)
        clusters_per_gene[library_id_slice_ind] = results[library_id]['cluster_per_gene']

    if return_anndata:
        print("    [P2R] Creating output AnnData...")
        import anndata
        # Get position of bins
        valid_xy_list = [r['xy_bin'][r['good_bins']] for r in results.values()]
        if not valid_xy_list:
             return None
        xy_bin = np.vstack(valid_xy_list)
        
        # Get labels of bins
        labels_bin = np.hstack([r['cluster_per_bin'][r['good_bins']] for r in results.values()])

        obs = {}
        obs['points2regions'] = labels_bin
        if len(results) > 1:
            obs['groupid'] =  np.hstack([[id]*len(r['cluster_per_bin']) for id, r in results.items()])

        # Multiply back features with the norm
        norms = 1.0 / np.hstack([r['norms'][r['good_bins']] for r in results.values()])
        norms[np.isinf(norms)] = 0
        norms = norms.reshape((-1,1))

        adata_p2r = anndata.AnnData(
            X=X_train.multiply(norms).tocsc(),
            obsm={'spatial' : xy_bin},
            obs=obs,
            var=pd.DataFrame(index=unique_genes)
        )
        adata_p2r.obs['points2regions'] = adata_p2r.obs['points2regions'].astype('category')
        if len(results) > 1:
            adata_p2r.obs['groupid'] = adata_p2r.obs['groupid'].astype('category')

        # Store spot-level assignments in uns
        reads = {}
        reads['x'] = xy[:,0]
        reads['y'] = xy[:,1]
        reads['labels'] = gene_labels
        reads['points2regions'] = clusters_per_gene
        if groupids is not None:
             reads['groupid'] = groupids
             
        reads_df = pd.DataFrame(reads)
        reads_df['labels'] = reads_df['labels'].astype('category')
        reads_df['points2regions'] = reads_df['points2regions'].astype('category')
        if groupids is not None:
             reads_df['groupid'] = reads_df['groupid'].astype('category')
             
        adata_p2r.uns['reads'] = reads_df
        return adata_p2r

    return clusters_per_gene

# ==============================================================================
# Functional Utilities
# ==============================================================================

def ssam_to_anndata(analysis, gene_names):
    """
    Converts SSAM analysis Local Maxima results into an AnnData object.
    Requires analysis.local_maxima to be computed.
    """
    print("    [SSAM] Converting SSAM Local Maxima to AnnData...")
    import anndata
    
    # Check if local maxima exists
    if hasattr(analysis, "local_maxima") and analysis.local_maxima is not None:
        centroids = analysis.local_maxima 
        # local_maxima is typically a DataFrame with 'x', 'y' or numpy array of coords
        
        if isinstance(centroids, pd.DataFrame):
            xy = centroids[['x', 'y']].values
        else:
            xy = centroids[:, :2]
            
        print(f"      - Found {len(xy)} potential cells.")
        
        # Check integrity
        if not hasattr(analysis, 'dataset') or analysis.dataset is None:
             logging.error("SSAM Analysis object missing dataset/kde_matrix.")
             return None
             
        # Sample KDE at centroids
        # Coordinates in local_maxima are in pixel space of the KDE map
        
        x_idx = np.clip(np.round(xy[:, 0]).astype(int), 0, analysis.dataset.width - 1)
        y_idx = np.clip(np.round(xy[:, 1]).astype(int), 0, analysis.dataset.height - 1)
        
        n_cells = len(xy)
        n_genes = len(gene_names)
        X = np.zeros((n_cells, n_genes))
        
        # We iterate genes to fetch slices
        for i, gene in enumerate(gene_names):
            # Access KDE matrix
            # Note: analysis.dataset[gene] returns the 2D density map for that gene
            density_map = analysis.dataset[gene] 
            X[:, i] = density_map[x_idx, y_idx]
            
        # Create AnnData
        obs = pd.DataFrame(index=[f"cell_{i}" for i in range(n_cells)])
        
        # Store centroids (optionally scaled to microns if we knew resolution offset, but pixels ok for now)
        obs['x_centroid_px'] = xy[:, 0]
        obs['y_centroid_px'] = xy[:, 1]
        
        adata_ssam = anndata.AnnData(X=X, obs=obs)
        adata_ssam.var_names = gene_names
        adata_ssam.obsm['spatial'] = xy 
        
        return adata_ssam
    else:
        logging.warning("    [SSAM] No local maxima found. Cannot create AnnData.")
        return None

# ==============================================================================
# Pipeline Execution Functions
# ==============================================================================

def run_points2regions_analysis(adata, output_dir, sample_tag, config):
    """
    Runs Points2Regions analysis on Xenium data.
    """
    # Create Subdirectory
    step_dir = os.path.join(output_dir, f"{sample_tag}_step2_segmentation_free", "points2regions")
    os.makedirs(step_dir, exist_ok=True)
    
    print(f"\n[Step 2-1] Segmentation-Free: Points2Regions Analysis")
    print(f"    - Output Directory: {step_dir}")
    
    p2r_config = config.get("segmentation_free", {}).get("points2regions", {})
    sigma = p2r_config.get("sigma", 2.0)
    n_clusters = p2r_config.get("n_clusters", 10)
    bin_width = p2r_config.get("bin_width", "auto")
    min_genes_per_bin = p2r_config.get("min_genes_per_bin", 5)
    
    if 'spots' not in adata.uns:
        logging.warning("    - 'spots' (transcripts) not found in adata.uns. Skipping.")
        return

    spots = adata.uns['spots']
    # Check column names
    if 'feature_name' in spots.columns:
        gene_labels = spots['feature_name'].values
    else:
        logging.error("    - 'feature_name' column missing in spots.")
        return

    if 'x_location' in spots.columns and 'y_location' in spots.columns:
        xy = spots[['x_location', 'y_location']].values
    else:
        logging.error("    - 'x_location'/'y_location' columns missing in spots.")
        return

    # Run Points2Regions
    adata_p2r = points2regions(
        xy, 
        gene_labels, 
        sigma=sigma, 
        n_clusters=n_clusters, 
        bin_width=bin_width, 
        min_genes_per_bin=min_genes_per_bin,
        return_anndata=True
    )

    if adata_p2r is not None:
        output_file = os.path.join(step_dir, f"{sample_tag}_step2_points2regions.h5ad")
        adata_p2r.write(output_file)
        print(f"    - Saved Points2Regions result to: {output_file}")
        
        # Plotting
        plt.figure(figsize=(10, 10))
        sc.pl.spatial(adata_p2r, color="points2regions", spot_size=10, show=False, title=f"Points2Regions (Sigma={sigma}, K={n_clusters})")
        plot_file = os.path.join(step_dir, f"{sample_tag}_step2_points2regions_map.png")
        plt.savefig(plot_file, bbox_inches='tight')
        plt.close()
        print(f"    - Saved P2R Map to: {plot_file}")
    else:
        print("    - Points2Regions returned None (failed).")

def run_ssam_analysis(adata, output_dir, sample_tag, config):
    # Create Subdirectory
    step_dir = os.path.join(output_dir, f"{sample_tag}_step2_segmentation_free", "ssam")
    os.makedirs(step_dir, exist_ok=True)
    
    print(f"\n[Step 2-2] Segmentation-Free: SSAM Analysis")
    print(f"    - Output Directory: {step_dir}")
    
    ssam_config = config.get("segmentation_free", {}).get("ssam", {})
    vf_res = ssam_config.get("ssam_vf_um_per_px", 2.0)
    kde_bw = ssam_config.get("kde_bandwidth_um", 2.5)
    norm_thres = ssam_config.get("norm_thres", 5)
    exp_thres = ssam_config.get("exp_thres", 0.2)
    min_dist = ssam_config.get("min_dist", 7)
    
    try:
        import ssam
        print(f"    - SSAM library found version {ssam.__version__ if hasattr(ssam, '__version__') else 'unknown'}")
        
        spots = None
        if 'spots' in adata.uns:
            spots = adata.uns['spots']
        else:
            print("    - 'spots' dataframe not found in adata.uns. Skipping SSAM.")
            return

        if 'x_location' not in spots.columns or 'y_location' not in spots.columns or 'feature_name' not in spots.columns:
             logging.error("    - standard columns (x_location, y_location, feature_name) missing in spots dataframe.")
             return
             
        # Determine dimensions
        x_min, x_max = spots['x_location'].min(), spots['x_location'].max()
        y_min, y_max = spots['y_location'].min(), spots['y_location'].max()
        width = int(np.ceil((x_max - x_min) / vf_res))
        height = int(np.ceil((y_max - y_min) / vf_res))
        
        print(f"    - Configuring SSAM: Res={vf_res}, KDE_BW={kde_bw}, Size={width}x{height}")
        
        # 1. Create SSAM Dataset
        genes = spots['feature_name'].unique()
        mrna_loci = []
        grouped = spots.groupby('feature_name')
        for g in genes:
            if g in grouped.groups:
                g_spots = grouped.get_group(g)
                coords = g_spots[['x_location', 'y_location']].values
                coords[:, 0] = (coords[:, 0] - x_min) / vf_res
                coords[:, 1] = (coords[:, 1] - y_min) / vf_res
                mrna_loci.append(coords)
            else:
                mrna_loci.append(np.empty((0, 2)))

        ds = ssam.SSAMDataset(list(genes), mrna_loci, width, height)
        
        # 2. Compute KDE
        print(f"    - Running SSAM Analysis (Fast KDE)...")
        analysis = ssam.SSAMAnalysis(ds, save_dir=step_dir)
        bw_px = kde_bw / vf_res
        analysis.run_fast_kde(bandwidth=bw_px, use_mmap=False, re_run=True)
        
        # 3. Find Local Maxima
        print("    - Finding Local Maxima / Cell Types...")
        # Note: If de-novo, we typically look for local maxima in total density or per-gene?
        analysis.find_local_maxima(min_dist=min_dist, norm_thres=norm_thres, exp_thres=exp_thres)
        
        # 4. Visualization (Enhanced)
        print("    - Generating SSAM Visualizations...")
        
        # A. Cell Type Map (Rich)
        try:
            plt.figure(figsize=(25, 30))
            # Use 'nipy_spectral' or similar if de-novo, or colors=None for ssam default
            ssam.plot.cell_types_map(analysis, rotate=3, colors=None)
            
            # Formatting
            plt.xlim(0, None) # Auto
            plt.ylim(None, 0) # Flip Y for image coords usually? Or match Xenium.
            
            # Scalebar (custom)
            # 1mm = 1000um. At vf=2.0 um/px, 1000um = 500 px.
            scale_len_px = 1000.0 / vf_res
            plot_scalebar(plt.gca(), 100, 100, scale_len_px, text='1mm', color='w', edge_color='k')
            
            plt.axis('off')
            
            plot_file = os.path.join(step_dir, f"{sample_tag}_step2_ssam_denovo_map.png")
            plt.savefig(plot_file, dpi=300, bbox_inches='tight')
            # Save PDF too
            plt.savefig(plot_file.replace('.png', '.pdf'), dpi=300, bbox_inches='tight')
            plt.close()
            print(f"    - Saved De-Novo Map to: {plot_file}")
            
        except Exception as e_plot:
            print(f"      [WARNING] Enhanced plotting failed ({e_plot}). Trying basic plot.")
            try:
                ssam.plot.plot_total_density(analysis)
                plt.savefig(os.path.join(step_dir, f"{sample_tag}_step2_ssam_density.png"))
                plt.close()
            except:
                pass

        # 5. Convert to AnnData & Downstream Analysis
        adata_ssam = ssam_to_anndata(analysis, list(genes))
        
        if adata_ssam is not None:
            print("    - Running Downstream Scanpy Analysis on SSAM results...")
            
            # Save Raw object
            adata_ssam.write(os.path.join(step_dir, f"{sample_tag}_step2_ssam_cells.h5ad"))
            
            # Normalize & Log (Simple for visualization)
            sc.pp.normalize_total(adata_ssam, target_sum=1e4)
            sc.pp.log1p(adata_ssam)
            
            # Highest Expr Genes
            sc.pl.highest_expr_genes(adata_ssam, n_top=20, show=False)
            plt.title("Highest Expressed Genes (SSAM Cells)")
            plt.savefig(os.path.join(step_dir, f"{sample_tag}_step2_ssam_highest_expr.png"), bbox_inches='tight')
            plt.close()
            
            # PCA & UMAP
            if adata_ssam.n_vars > 0 and adata_ssam.n_obs > 10:
                print("      > Computing PCA...")
                sc.tl.pca(adata_ssam, svd_solver='arpack')
                
                print("      > Computing Neighbors & UMAP...")
                sc.pp.neighbors(adata_ssam, n_neighbors=15, n_pcs=min(20, adata_ssam.n_vars-1))
                sc.tl.umap(adata_ssam)
                
                # Clustering (Leiden) for coloring UMAP
                sc.tl.leiden(adata_ssam, resolution=0.5)
                
                # Plot PCA
                sc.pl.pca(adata_ssam, color=['leiden'], show=False)
                plt.savefig(os.path.join(step_dir, f"{sample_tag}_step2_ssam_pca.png"), bbox_inches='tight')
                plt.close()
                
                # Plot UMAP
                sc.pl.umap(adata_ssam, color=['leiden'], title="UMAP (SSAM Cells)", show=False)
                plt.savefig(os.path.join(step_dir, f"{sample_tag}_step2_ssam_umap.png"), bbox_inches='tight')
                plt.close()
                
                # Save processed
                adata_ssam.write(os.path.join(step_dir, f"{sample_tag}_step2_ssam_cells_processed.h5ad"))
                print(f"    - Saved processed SSAM AnnData & Plots.")
            else:
                print("      [WARNING] Not enough cells/genes for PCA/UMAP.")

    except ImportError:
        print("    - 'ssam' library not installed. Skipping SSAM analysis.")
    except Exception as e:
        logging.error(f"    - SSAM analysis failed: {e}")
        import traceback
        traceback.print_exc()

def run_overlaps_analysis(adata, output_dir, sample_tag, config):
    # Create Subdirectory
    step_dir = os.path.join(output_dir, f"{sample_tag}_step2_segmentation_free", "overlaps")
    os.makedirs(step_dir, exist_ok=True)
    
    print(f"  - Running Overlaps analysis for {sample_tag}...")
    
    try:
        import ovrlpy as ovrlp
        import matplotlib.pyplot as plt
        
        print("    - ovrlpy library found.")
        
        if 'X_pca' not in adata.obsm:
             sc.tl.pca(adata, svd_solver='arpack')
        
        spots = adata.uns.get('spots', None)
        
        if spots is not None and 'z_location' in spots.columns:
             # Prepare DF
             df = spots.rename(columns={'x_location': 'x', 'y_location': 'y', 'z_location': 'z', 'feature_name': 'gene'})
             
             z_median = df['z'].median()
             df_top = df[df['z'] < z_median]
             df_bot = df[df['z'] > z_median]
             
             x_max = df['x'].max()
             y_max = df['y'].max()
             
             distance = None
             genes = list(adata.var_names)
             
             for i, g in enumerate(genes):
                 if g not in df['gene'].values: continue
                 
                 try:
                     hist_top_g = ovrlp.create_histogram(df_top, genes=[g], x_max=x_max, y_max=y_max, KDE_bandwidth=1.0)
                     hist_bot_g = ovrlp.create_histogram(df_bot, genes=[g], x_max=x_max, y_max=y_max, KDE_bandwidth=1.0)
                     
                     if distance is None:
                         distance = np.zeros_like(hist_top_g)

                     diff = np.abs(hist_top_g - hist_bot_g)
                     distance += diff
                     
                 except Exception:
                     continue

             if distance is not None:
                 try:
                     distance = ovrlp.gaussian_filter(distance, sigma=1)
                 except:
                     pass
                     
                 plt.figure(figsize=(10, 10))
                 plt.imshow(distance, cmap='inferno')
                 plt.title(f"Overlaps Incoherence Map")
                 plt.axis('off')
                 plot_file = os.path.join(step_dir, f"{sample_tag}_step2_overlaps_map.png")
                 plt.savefig(plot_file, dpi=300)
                 plt.close()
                 print(f"    - Saved Overlaps Map to {plot_file}")
             else:
                 print("    - No valid overlap histograms generated.")
        else:
             print("    - 'spots' or 'z_location' missing. Skipping Overlaps.")

    except ImportError:
        print("    - 'ovrlpy' library not installed. Skipping Overlaps analysis.")
    except Exception as e:
        logging.error(f"    - Overlaps analysis failed: {e}")

def run_step2(config):
    """
    Main function for Step 2: Segmentation Free Analysis
    """
    print("============================================================")
    print("Step 2: Segmentation Free Analysis")
    print("============================================================")
    
    input_path = config["output_dir"] # output_dir for Step 2
    sample_tag = config["sample_tag"]
    
    # Load Data (Step 0 or Step 1 Output)
    if 'previous_step_adata_path' in config and config['previous_step_adata_path']:
        file_path = config['previous_step_adata_path']
        print(f"  > Loading input from explicit path: {file_path}")
    else:
        # Fallback to current or parent
        # We try to find Step 0 or Step 1 data
        # Assume standard structure if flat:
        file_path = os.path.join(input_path, f"{sample_tag}.h5ad")
        
        if not os.path.exists(file_path):
             # Try parent directory structure (Step 1 first as it has QC)
             parent_dir = os.path.dirname(input_path)
             step1_input = os.path.join(parent_dir, "step1_exploration", f"{sample_tag}_step1_exploration.h5ad")
             if os.path.exists(step1_input):
                 file_path = step1_input
             else:
                 # Try step0
                 step0_input = os.path.join(parent_dir, "step0_formatting", f"{sample_tag}.h5ad")
                 if os.path.exists(step0_input):
                     file_path = step0_input

    if not os.path.exists(file_path):
        # Fallback to step1 if exists? usually we utilize step0 basic file which has spots
        logging.error(f"Input file not found at {file_path}. Cannot run Step 2.")
        return

    print(f"[Step 2-0] Loading Data: {file_path}")
    
    try:
        adata = sc.read_h5ad(file_path)
    except Exception as e:
        logging.error(f"Error loading AnnData: {e}")
        return
    
    # Verify spots in uns
    if 'spots' not in adata.uns:
         print("    [WARNING] 'spots' key missing in adata.uns. Analysis requiring transcripts will fail.")
    
    seg_config = config.get("segmentation_free", {})
    
    if seg_config.get("run_points2regions", False):
        run_points2regions_analysis(adata, input_path, sample_tag, config)
    else:
        print("[Step 2-1] Points2Regions skipped.")

    if seg_config.get("run_ssam", False):
        run_ssam_analysis(adata, input_path, sample_tag, config)
    
    if seg_config.get("run_overlaps", False):
        run_overlaps_analysis(adata, input_path, sample_tag, config)

    print("\n=== Step 2 Analysis Complete ===")

if __name__ == "__main__":
    import yaml
    print("Xenium Pipeline Step 2: Segmentation Free Analysis (Standalone)")
    config_path = os.path.join(os.path.dirname(__file__), 'config.yaml')
    if os.path.exists(config_path):
        with open(config_path) as f:
            config = yaml.safe_load(f)
        run_step2(config)
