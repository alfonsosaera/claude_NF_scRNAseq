#!/usr/bin/env python3
"""
BBKNN integration script for batch-corrected single-cell RNA-seq analysis.
Merges multiple per-sample AnnData objects and performs BBKNN batch correction.
"""

import sys
import json
import argparse
import logging
import traceback
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import anndata

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Set Scanpy settings
try:
    sc.settings.verbosity = 3
    sc.settings.set_figure_params(dpi=100, facecolor='white')
except AttributeError:
    # Older versions of scanpy may not have settings attribute
    pass


def load_bbknn_config(config_file: str) -> Dict:
    """Load BBKNN configuration from JSON file."""
    try:
        with open(config_file, 'r') as f:
            config = json.load(f)
        logger.info(f"Loaded BBKNN config from {config_file}")
        return config
    except Exception as e:
        logger.error(f"Failed to load config file: {e}")
        raise


def load_anndata_files(input_files: List[str]) -> List[anndata.AnnData]:
    """Load multiple AnnData files."""
    adata_list = []
    for file_path in input_files:
        try:
            adata = anndata.read_h5ad(file_path)
            logger.info(f"Loaded {Path(file_path).name}: shape={adata.shape}")
            adata_list.append(adata)
        except Exception as e:
            logger.error(f"Failed to load {file_path}: {e}")
            raise

    return adata_list


def merge_anndata(adata_list: List[anndata.AnnData]) -> anndata.AnnData:
    """Merge multiple AnnData objects."""
    try:
        # Concatenate all samples
        adata = anndata.concat(adata_list, axis=0, join='outer', label='sample_id', keys=None)
        logger.info(f"Merged {len(adata_list)} samples: shape={adata.shape}")

        # Fill missing values in sparse matrix if using outer join
        if hasattr(adata.X, 'toarray'):
            adata.X.data = np.nan_to_num(adata.X.data)
            adata.X.eliminate_zeros()

        return adata

    except Exception as e:
        logger.error(f"Failed to merge AnnData objects: {e}")
        raise


def ensure_batch_metadata(adata: anndata.AnnData, batch_key: str) -> anndata.AnnData:
    """Ensure batch metadata exists in obs."""
    if batch_key not in adata.obs.columns:
        if 'sample_id' in adata.obs.columns:
            logger.info(f"Using 'sample_id' as batch key ('{batch_key}' not found)")
            adata.obs[batch_key] = adata.obs['sample_id']
        else:
            logger.warning(f"Neither '{batch_key}' nor 'sample_id' found in obs")
            logger.info("Assigning batch=0 to all cells")
            adata.obs[batch_key] = '0'

    logger.info(f"Using batch key: '{batch_key}'")
    logger.info(f"Batches: {adata.obs[batch_key].unique().tolist()}")

    return adata


def clean_nans(adata: anndata.AnnData) -> anndata.AnnData:
    """Remove genes with NaN values that can cause issues with HVG selection."""
    # Check for NaN values in the expression matrix
    if hasattr(adata.X, 'data'):  # sparse matrix
        nan_genes = []
        for i in range(adata.n_vars):
            if np.isnan(adata.X[:, i].data).any():
                nan_genes.append(i)
    else:  # dense matrix
        nan_mask = np.isnan(adata.X).any(axis=0)
        nan_genes = np.where(nan_mask)[0].tolist()

    if nan_genes:
        logger.info(f"Found {len(nan_genes)} genes with NaN values, removing them")
        adata = adata[:, ~adata.var_names.isin(adata.var_names[nan_genes])].copy()

    return adata


def recompute_hvg(adata: anndata.AnnData, config: Dict) -> anndata.AnnData:
    """Recompute HVGs on merged data."""
    n_hvgs = config.get('n_hvgs', 3000)
    hvg_flavor = config.get('hvg_flavor', 'seurat')

    # Clean NaN values first to avoid HVG selection errors
    adata = clean_nans(adata)

    # Use variable genes - seurat flavor works well with log-normalized data
    sc.pp.highly_variable_genes(adata, flavor=hvg_flavor, n_top_genes=n_hvgs)
    logger.info(f"Selected {n_hvgs} HVGs on merged data using flavor '{hvg_flavor}'")

    return adata


def scale_and_pca_merged(adata: anndata.AnnData, config: Dict) -> anndata.AnnData:
    """Scale data and run PCA on merged data."""
    # Scale
    sc.pp.scale(adata, max_value=10)
    logger.info("Scaled merged data")

    # PCA
    n_pcs = config.get('n_pcs', 50)
    sc.tl.pca(adata, n_comps=n_pcs)
    logger.info(f"Computed PCA with {n_pcs} components")

    return adata


def run_bbknn(
    adata: anndata.AnnData,
    config: Dict,
    batch_key: str
) -> anndata.AnnData:
    """Run BBKNN batch correction."""
    try:
        neighbors_within_batch = config.get('neighbors_within_batch', 3)
        n_pcs = config.get('n_pcs', 50)

        sc.external.pp.bbknn(
            adata,
            batch_key=batch_key,
            neighbors_within_batch=neighbors_within_batch,
            n_pcs=n_pcs
        )
        logger.info(f"BBKNN computed with neighbors_within_batch={neighbors_within_batch}")

    except Exception as e:
        logger.error(f"BBKNN failed: {e}")
        logger.warning("Falling back to standard kNN without batch correction")
        sc.pp.neighbors(adata, n_neighbors=15, n_pcs=n_pcs)

    return adata


def compute_umap(adata: anndata.AnnData, config: Dict) -> anndata.AnnData:
    """Compute UMAP embedding."""
    umap_min_dist = config.get('umap_min_dist', 0.1)
    umap_spread = config.get('umap_spread', 1.0)
    umap_metric = config.get('umap_metric', 'correlation')

    try:
        # Try with metric parameter (older scanpy versions)
        sc.tl.umap(adata, min_dist=umap_min_dist, spread=umap_spread, metric=umap_metric)
        logger.info(f"Computed UMAP with min_dist={umap_min_dist}, spread={umap_spread}, metric={umap_metric}")
    except TypeError as e:
        # Fallback if metric parameter not supported (newer versions)
        if 'metric' in str(e):
            logger.warning(f"UMAP metric parameter not supported, computing without it: {e}")
            sc.tl.umap(adata, min_dist=umap_min_dist, spread=umap_spread)
            logger.info(f"Computed UMAP with min_dist={umap_min_dist}, spread={umap_spread}")
        else:
            raise

    return adata


def leiden_clustering(adata: anndata.AnnData, config: Dict) -> anndata.AnnData:
    """Run Leiden clustering."""
    resolution = config.get('leiden_resolution', 1.0)

    sc.tl.leiden(adata, resolution=resolution, key_added='leiden')
    logger.info(f"Leiden clustering with resolution={resolution}")

    n_clusters = len(adata.obs['leiden'].unique())
    logger.info(f"Found {n_clusters} clusters")

    return adata


def generate_umap_plots(
    adata: anndata.AnnData,
    batch_key: str,
    plots_dir: str
) -> None:
    """Generate UMAP plots colored by different features."""
    plots_dir = Path(plots_dir)
    plots_dir.mkdir(parents=True, exist_ok=True)

    try:
        # UMAP by batch
        fig, ax = plt.subplots(figsize=(10, 8))
        scatter = ax.scatter(
            adata.obsm['X_umap'][:, 0],
            adata.obsm['X_umap'][:, 1],
            c=pd.Categorical(adata.obs[batch_key]).codes,
            cmap='tab20',
            alpha=0.7,
            s=30
        )
        ax.set_xlabel('UMAP1')
        ax.set_ylabel('UMAP2')
        ax.set_title(f'UMAP colored by {batch_key}')
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label(batch_key)
        plt.tight_layout()
        output_file = plots_dir / f"umap_by_{batch_key}.png"
        plt.savefig(output_file, dpi=100, bbox_inches='tight')
        plt.close()
        logger.info(f"Saved UMAP plot to {output_file}")

        # UMAP by leiden clustering
        if 'leiden' in adata.obs:
            fig, ax = plt.subplots(figsize=(10, 8))
            scatter = ax.scatter(
                adata.obsm['X_umap'][:, 0],
                adata.obsm['X_umap'][:, 1],
                c=pd.Categorical(adata.obs['leiden']).codes,
                cmap='tab20',
                alpha=0.7,
                s=30
            )
            ax.set_xlabel('UMAP1')
            ax.set_ylabel('UMAP2')
            ax.set_title('UMAP colored by Leiden cluster')
            cbar = plt.colorbar(scatter, ax=ax)
            cbar.set_label('Leiden cluster')
            plt.tight_layout()
            output_file = plots_dir / "umap_by_leiden.png"
            plt.savefig(output_file, dpi=100, bbox_inches='tight')
            plt.close()
            logger.info(f"Saved UMAP plot to {output_file}")

        # UMAP by sample_id if available
        if 'sample_id' in adata.obs:
            fig, ax = plt.subplots(figsize=(10, 8))
            scatter = ax.scatter(
                adata.obsm['X_umap'][:, 0],
                adata.obsm['X_umap'][:, 1],
                c=pd.Categorical(adata.obs['sample_id']).codes,
                cmap='tab20',
                alpha=0.7,
                s=30
            )
            ax.set_xlabel('UMAP1')
            ax.set_ylabel('UMAP2')
            ax.set_title('UMAP colored by sample_id')
            cbar = plt.colorbar(scatter, ax=ax)
            cbar.set_label('sample_id')
            plt.tight_layout()
            output_file = plots_dir / "umap_by_sample_id.png"
            plt.savefig(output_file, dpi=100, bbox_inches='tight')
            plt.close()
            logger.info(f"Saved UMAP plot to {output_file}")

        # UMAP by other metadata columns
        metadata_cols = [col for col in adata.obs.columns
                        if col not in ['batch', 'sample_id', 'leiden', 'doublet_score', 'doublet_class']
                        and adata.obs[col].dtype == 'object']

        for col in metadata_cols[:3]:  # Limit to first 3 to avoid too many plots
            fig, ax = plt.subplots(figsize=(10, 8))
            scatter = ax.scatter(
                adata.obsm['X_umap'][:, 0],
                adata.obsm['X_umap'][:, 1],
                c=pd.Categorical(adata.obs[col]).codes,
                cmap='tab20',
                alpha=0.7,
                s=30
            )
            ax.set_xlabel('UMAP1')
            ax.set_ylabel('UMAP2')
            ax.set_title(f'UMAP colored by {col}')
            cbar = plt.colorbar(scatter, ax=ax)
            cbar.set_label(col)
            plt.tight_layout()
            output_file = plots_dir / f"umap_by_{col}.png"
            plt.savefig(output_file, dpi=100, bbox_inches='tight')
            plt.close()
            logger.info(f"Saved UMAP plot to {output_file}")

    except Exception as e:
        logger.error(f"Error generating UMAP plots: {e}")
        raise


def main():
    parser = argparse.ArgumentParser(
        description='BBKNN integration for single-cell RNA-seq'
    )
    parser.add_argument('--input-files', nargs='+', required=True, help='Input H5AD files')
    parser.add_argument('--output', required=True, help='Output merged H5AD file')
    parser.add_argument('--config', required=True, help='BBKNN config JSON file')
    parser.add_argument('--plots-dir', default='.', help='Output directory for plots')

    args = parser.parse_args()

    try:
        # Load configuration
        config = load_bbknn_config(args.config)

        # Load and merge data
        adata_list = load_anndata_files(args.input_files)
        adata = merge_anndata(adata_list)

        # Ensure batch metadata
        batch_key = config.get('batch_key', 'batch')
        adata = ensure_batch_metadata(adata, batch_key)

        # Recompute HVGs on merged data
        adata = recompute_hvg(adata, config)

        # Scale and PCA
        adata = scale_and_pca_merged(adata, config)

        # BBKNN
        adata = run_bbknn(adata, config, batch_key)

        # UMAP
        adata = compute_umap(adata, config)

        # Leiden clustering
        adata = leiden_clustering(adata, config)

        # Generate plots
        generate_umap_plots(adata, batch_key, args.plots_dir)

        # Save output
        adata.write_h5ad(args.output, compression='gzip')
        logger.info(f"Saved integrated AnnData to {args.output}")

        logger.info("Integration completed successfully")
        sys.exit(0)

    except Exception as e:
        logger.error(f"Fatal error: {e}")
        logger.debug(traceback.format_exc())
        sys.exit(1)


if __name__ == '__main__':
    main()
