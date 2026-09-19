#!/usr/bin/env python3
"""Processes ESMC protein embeddings, clusters them, and outputs species UMAP panels.

Outputs:
    Saves individual species UMAP plots, merged panel image, and cluster CSVs.
"""
import argparse
import datetime
from pathlib import Path
import resource
import time
import anndata as ad
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
import scanpy as sc
from tqdm import tqdm

###################################### CONFIG
DEFAULT_OUTDIR = "data/scanpy/tmp/"
N_PCS = 50
N_NEIGHBORS = 15
LEIDEN_RES = 1.0
PANEL_SIZE_INCHES = 2.5
PANEL_DPI = 300
COLOR_BG = "#E5E7EB"
COLORMAP_SPEC = LinearSegmentedColormap.from_list("teal_copper", ["#68d1de", "#d08c66"])
PT_SIZE_BG = 0.8
PT_SIZE_FG = 2.2
ALPHA_BG = 0.20
ALPHA_HIGHLIGHT = 0.85

###################################### MAIN
def main():
    ###################################### CLI_ARGS
    parser = argparse.ArgumentParser(description="ESMC Proteome Embedding Clustering")
    parser.add_argument("--npzs", type=str, required=True, help="Comma-separated paths to .npz files")
    parser.add_argument("--outdir", type=str, default=DEFAULT_OUTDIR, help="Base output directory")
    parser.add_argument("--n_pcs", type=int, default=N_PCS, help="Number of PCs (0 to bypass PCA)")
    parser.add_argument("--n_neighbors", type=int, default=N_NEIGHBORS, help="Number of neighbors for KNN")
    parser.add_argument("--resolution", type=float, default=LEIDEN_RES, help="Leiden cluster resolution")
    args = parser.parse_args()
    npz_files = [p.strip() for p in args.npzs.split(",")]
    target_dir = Path(args.outdir) / "_".join([Path(f).stem for f in npz_files])
    target_dir.mkdir(parents=True, exist_ok=True)
    ###################################### LOAD_DATA
    start_time = time.time()
    adatas = []
    for f in tqdm(npz_files):
        data = np.load(f)
        k_emb = next(k for k in data.files if 'embeddings' in k)
        k_ids = next(k for k in data.files if 'ids' in k)
        adata = ad.AnnData(X=data[k_emb].astype(np.float32))
        adata.obs['protein_id'] = data[k_ids]
        adata.obs['species'] = Path(f).stem.replace('_reps', '')
        adatas.append(adata)
    adata = ad.concat(adatas, join="outer") if len(adatas) > 1 else adatas[0]
    adata.obs_names_make_unique()
    ###################################### PREPROCESSING
    norms = np.linalg.norm(adata.X, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    adata.X = adata.X / norms
    ###################################### CLUSTERING
    if args.n_pcs > 0:
        sc.tl.pca(adata, n_comps=args.n_pcs)
        sc.pp.neighbors(adata, n_neighbors=args.n_neighbors, n_pcs=args.n_pcs, metric='cosine')
    else:
        sc.pp.neighbors(adata, n_neighbors=args.n_neighbors, use_rep='X', metric='cosine')
    sc.tl.leiden(adata, resolution=args.resolution, key_added='cluster', flavor='igraph')
    sc.tl.umap(adata)
    ###################################### PLOTTING
    species_list = adata.obs['species'].unique()
    umap_coords = adata.obsm['X_umap']
    cluster_sp_frac = pd.crosstab(adata.obs['cluster'], adata.obs['species'], normalize='index')
    for sp in tqdm(species_list, desc="Saving individual plots"):
        fig, ax = plt.subplots(figsize=(PANEL_SIZE_INCHES, PANEL_SIZE_INCHES), dpi=PANEL_DPI)
        is_sp = (adata.obs['species'] == sp).values
        spec_scores = adata.obs['cluster'].map(cluster_sp_frac[sp]).values
        ax.scatter(umap_coords[~is_sp, 0], umap_coords[~is_sp, 1], c=COLOR_BG, s=PT_SIZE_BG, alpha=ALPHA_BG, rasterized=True, linewidths=0)
        fg_idx = np.where(is_sp)[0]
        sort_order = fg_idx[np.argsort(spec_scores[is_sp])]
        ax.scatter(umap_coords[sort_order, 0], umap_coords[sort_order, 1], c=spec_scores[sort_order], cmap=COLORMAP_SPEC, vmin=0, vmax=1, s=PT_SIZE_FG, alpha=ALPHA_HIGHLIGHT, rasterized=True, linewidths=0)
        ax.set_title(f"{sp} (n={is_sp.sum():,})", fontsize=8, fontweight='bold', pad=3)
        ax.set_xlabel("UMAP 1", fontsize=7); ax.set_ylabel("UMAP 2", fontsize=7)
        ax.tick_params(labelsize=6)
        ax.set_aspect('equal', 'datalim')
        plt.tight_layout()
        plt.savefig(target_dir / f"umap_{sp}.png", dpi=PANEL_DPI, bbox_inches='tight')
        plt.close(fig)
    ###################################### MERGE_PANEL
    n_sp = len(species_list)
    n_cols = min(n_sp, 4)
    n_rows = (n_sp + n_cols - 1) // n_cols
    fig = plt.figure(figsize=(PANEL_SIZE_INCHES * n_cols, PANEL_SIZE_INCHES * n_rows + 0.5), dpi=PANEL_DPI)
    gs = fig.add_gridspec(n_rows + 1, n_cols, height_ratios=[1] * n_rows + [0.06], hspace=0.30, wspace=0.25)
    axes = [fig.add_subplot(gs[r, c]) for r in range(n_rows) for c in range(n_cols)]
    for idx, sp in enumerate(species_list):
        ax = axes[idx]
        is_sp = (adata.obs['species'] == sp).values
        spec_scores = adata.obs['cluster'].map(cluster_sp_frac[sp]).values
        ax.scatter(umap_coords[~is_sp, 0], umap_coords[~is_sp, 1], c=COLOR_BG, s=PT_SIZE_BG, alpha=ALPHA_BG, rasterized=True, linewidths=0)
        fg_idx = np.where(is_sp)[0]
        sort_order = fg_idx[np.argsort(spec_scores[is_sp])]
        ax.scatter(umap_coords[sort_order, 0], umap_coords[sort_order, 1], c=spec_scores[sort_order], cmap=COLORMAP_SPEC, vmin=0, vmax=1, s=PT_SIZE_FG, alpha=ALPHA_HIGHLIGHT, rasterized=True, linewidths=0)
        ax.set_title(f"{sp} (n={is_sp.sum():,})", fontsize=8, fontweight='bold', pad=3)
        ax.set_xlabel("UMAP 1", fontsize=7, labelpad=1); ax.set_ylabel("UMAP 2", fontsize=7, labelpad=1)
        ax.tick_params(labelsize=6)
        ax.set_aspect('equal', 'datalim')
    for idx in range(n_sp, len(axes)): fig.delaxes(axes[idx])
    cax = fig.add_subplot(gs[-1, :])
    sm = plt.cm.ScalarMappable(cmap=COLORMAP_SPEC, norm=plt.Normalize(vmin=0, vmax=1))
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cax, orientation='horizontal')
    cbar.set_label('Cluster Species Fraction', fontsize=8)
    cbar.ax.tick_params(labelsize=6)
    plt.savefig(target_dir / "umap_species_panel.png", dpi=PANEL_DPI, bbox_inches='tight')
    plt.close(fig)
    ###################################### SAVE_RESULTS
    adata.obs[['protein_id', 'species', 'cluster']].to_csv(target_dir / "protein_clusters.csv", index=False)
    if len(npz_files) > 1:
        enrichment = cluster_sp_frac * 100
        enrichment.to_csv(target_dir / "species_cluster_enrichment.csv")
    ###################################### LOGGING
    max_mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    exec_time = time.time() - start_time
    date_str = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[PYTHON-INFO] {date_str} {exec_time:.2f}s {max_mem:.2f}MB")

if __name__ == "__main__":
    main()