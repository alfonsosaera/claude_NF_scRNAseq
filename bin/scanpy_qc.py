#!/usr/bin/env python3
"""
Scanpy QC and filtering script for single-cell RNA-seq data.
Performs QC metrics calculation, cell/gene filtering, normalization, and HVG selection.
"""

import sys
import json
import argparse
import logging
import traceback
from pathlib import Path
from typing import Dict, Optional

import numpy as np
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
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=100, facecolor='white')


def load_qc_config(config_file: str) -> Dict:
    """Load QC configuration from JSON file."""
    try:
        with open(config_file, 'r') as f:
            config = json.load(f)
        logger.info(f"Loaded QC config from {config_file}")
        return config
    except Exception as e:
        logger.error(f"Failed to load config file: {e}")
        raise


def load_anndata(input_file: str) -> anndata.AnnData:
    """Load AnnData object from H5AD file."""
    try:
        adata = sc.read_h5ad(input_file)
        logger.info(f"Loaded AnnData from {input_file}: shape={adata.shape}")
        return adata
    except Exception as e:
        logger.error(f"Failed to load AnnData file: {e}")
        raise


def add_metadata(adata: anndata.AnnData, sample_id: str, metadata: Dict) -> anndata.AnnData:
    """Add sample metadata to AnnData obs."""
    adata.obs['sample_id'] = sample_id

    for key, value in metadata.items():
        if key != 'sample_id':
            adata.obs[key] = value

    logger.info(f"Added metadata to obs: {list(metadata.keys())}")
    return adata


def flag_mitochondrial_genes(adata: anndata.AnnData, mt_prefix: str) -> anndata.AnnData:
    """Flag mitochondrial genes."""
    adata.var['mt'] = adata.var_names.str.startswith(mt_prefix)
    n_mt = adata.var['mt'].sum()
    logger.info(f"Flagged {n_mt} mitochondrial genes with prefix '{mt_prefix}'")
    return adata


def calculate_qc_metrics(adata: anndata.AnnData) -> anndata.AnnData:
    """Calculate QC metrics."""
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True, log1p=False)
    logger.info("Calculated QC metrics")
    return adata


def generate_prefilter_plots(
    adata: anndata.AnnData,
    sample_id: str,
    plots_dir: str
) -> None:
    """Generate pre-filtering QC plots."""
    plots_dir = Path(plots_dir)
    plots_dir.mkdir(parents=True, exist_ok=True)

    try:
        # Histogram plots
        fig, axes = plt.subplots(1, 3, figsize=(15, 4))
        fig.suptitle(f"{sample_id} - Pre-filtering QC Metrics")

        axes[0].hist(adata.obs['n_genes_by_counts'], bins=50, color='steelblue', edgecolor='black')
        axes[0].set_xlabel('Genes per Cell')
        axes[0].set_ylabel('Count')
        axes[0].set_title('n_genes_by_counts')

        axes[1].hist(adata.obs['total_counts'], bins=50, color='steelblue', edgecolor='black')
        axes[1].set_xlabel('Total Counts')
        axes[1].set_ylabel('Count')
        axes[1].set_title('total_counts')
        axes[1].set_yscale('log')

        axes[2].hist(adata.obs['pct_counts_mt'], bins=50, color='steelblue', edgecolor='black')
        axes[2].set_xlabel('% Mitochondrial Counts')
        axes[2].set_ylabel('Count')
        axes[2].set_title('pct_counts_mt')

        plt.tight_layout()
        output_file = plots_dir / f"{sample_id}_prefilter_histograms.png"
        plt.savefig(output_file, dpi=100, bbox_inches='tight')
        plt.close()
        logger.info(f"Saved prefilter histograms to {output_file}")

        # Violin plots
        fig, axes = plt.subplots(1, 3, figsize=(15, 4))
        fig.suptitle(f"{sample_id} - Pre-filtering QC Violin Plots")

        metrics = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']
        for i, metric in enumerate(metrics):
            parts = axes[i].violinplot([adata.obs[metric].values], positions=[0], widths=0.7)
            axes[i].set_ylabel(metric)
            axes[i].set_title(metric)
            axes[i].set_xticks([])

        plt.tight_layout()
        output_file = plots_dir / f"{sample_id}_prefilter_violins.png"
        plt.savefig(output_file, dpi=100, bbox_inches='tight')
        plt.close()
        logger.info(f"Saved prefilter violins to {output_file}")

        # Scatter plot
        fig, ax = plt.subplots(figsize=(8, 6))
        scatter = ax.scatter(
            adata.obs['n_genes_by_counts'],
            adata.obs['total_counts'],
            c=adata.obs['pct_counts_mt'],
            cmap='viridis',
            alpha=0.6,
            s=20
        )
        ax.set_xlabel('Genes per Cell')
        ax.set_ylabel('Total Counts')
        ax.set_title(f"{sample_id} - Total Counts vs Genes (colored by % MT)")
        ax.set_yscale('log')
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label('% Mitochondrial Counts')
        plt.tight_layout()
        output_file = plots_dir / f"{sample_id}_prefilter_scatter.png"
        plt.savefig(output_file, dpi=100, bbox_inches='tight')
        plt.close()
        logger.info(f"Saved prefilter scatter to {output_file}")

    except Exception as e:
        logger.error(f"Error generating prefilter plots: {e}")
        raise


def apply_filters(
    adata: anndata.AnnData,
    config: Dict
) -> anndata.AnnData:
    """Apply cell and gene filters based on configuration."""
    filters = config.get('filters', {})
    n_cells_before = adata.n_obs
    n_genes_before = adata.n_vars

    # Cell filters
    min_genes = filters.get('min_genes_per_cell', 200)
    max_genes = filters.get('max_genes_per_cell')
    min_counts = filters.get('min_counts_per_cell', 500)
    max_counts = filters.get('max_counts_per_cell')
    max_pct_mt = filters.get('max_pct_mt', 20)

    # Apply cell filters
    adata = adata[adata.obs['n_genes_by_counts'] >= min_genes].copy()
    if max_genes is not None:
        adata = adata[adata.obs['n_genes_by_counts'] <= max_genes].copy()

    adata = adata[adata.obs['total_counts'] >= min_counts].copy()
    if max_counts is not None:
        adata = adata[adata.obs['total_counts'] <= max_counts].copy()

    adata = adata[adata.obs['pct_counts_mt'] <= max_pct_mt].copy()

    # Gene filters
    min_cells = filters.get('min_cells_per_gene', 3)
    sc.pp.filter_genes(adata, min_cells=min_cells)

    n_cells_after = adata.n_obs
    n_genes_after = adata.n_vars

    logger.info(f"Cells: {n_cells_before} → {n_cells_after} ({100*n_cells_after/n_cells_before:.1f}%)")
    logger.info(f"Genes: {n_genes_before} → {n_genes_after} ({100*n_genes_after/n_genes_before:.1f}%)")

    return adata


def run_doublet_detection(adata: anndata.AnnData) -> anndata.AnnData:
    """Run doublet detection using scDblFinder via rpy2."""
    try:
        import rpy2.robjects as ro
        from rpy2.robjects.packages import importr
        from rpy2.robjects import pandas2ri

        # Enable automatic conversion
        pandas2ri.activate()

        # Import scDblFinder
        scdblfinder = importr('scDblFinder')

        logger.info("Running doublet detection with scDblFinder...")

        # Convert to dense matrix for R
        X = np.array(adata.X.toarray()) if hasattr(adata.X, 'toarray') else np.array(adata.X)

        # Create R dataframe
        ro.r('library(SingleCellExperiment)')
        ro.r.assign('X', X.T)  # Transpose: genes x cells
        ro.r.assign('gene_names', adata.var_names.tolist())
        ro.r.assign('cell_names', adata.obs_names.tolist())

        # Create SingleCellExperiment object
        ro.r('''
            sce <- SingleCellExperiment(assays=list(counts=X))
            rownames(sce) <- gene_names
            colnames(sce) <- cell_names
        ''')

        # Run scDblFinder
        ro.r('sce <- scDblFinder(sce)')

        # Extract results
        doublet_class = ro.r('sce$scDblFinder.class')
        doublet_score = ro.r('sce$scDblFinder.score')

        # Add to AnnData
        adata.obs['doublet_class'] = list(doublet_class)
        adata.obs['doublet_score'] = np.array(doublet_score)

        n_doublets = (adata.obs['doublet_class'] == 'doublet').sum()
        logger.info(f"Identified {n_doublets} doublets ({100*n_doublets/len(adata):.1f}%)")

    except ImportError:
        logger.warning("rpy2 or scDblFinder not available, skipping doublet detection")
        logger.info("To use doublet detection, install: pip install rpy2")
        logger.info("And in R: install.packages('scDblFinder')")
        adata.obs['doublet_score'] = np.nan
        adata.obs['doublet_class'] = 'unknown'

    except Exception as e:
        logger.warning(f"Doublet detection failed: {e}")
        logger.debug(traceback.format_exc())
        adata.obs['doublet_score'] = np.nan
        adata.obs['doublet_class'] = 'unknown'

    return adata


def normalize_and_hvg(adata: anndata.AnnData, config: Dict) -> anndata.AnnData:
    """Normalize data and select highly variable genes."""
    normalization = config.get('normalization', {})
    target_sum = normalization.get('target_sum', 1e6)

    # Normalize
    sc.pp.normalize_total(adata, target_sum=target_sum)
    logger.info(f"Normalized total counts to {target_sum}")

    sc.pp.log1p(adata)
    logger.info("Log-transformed data")

    # HVG selection
    n_hvgs = config.get('n_hvgs', 3000)
    hvg_flavor = config.get('hvg_flavor', 'seurat_v3')

    sc.pp.highly_variable_genes(adata, flavor=hvg_flavor, n_top_genes=n_hvgs)
    logger.info(f"Selected {n_hvgs} HVGs using flavor '{hvg_flavor}'")

    return adata


def scale_and_pca(adata: anndata.AnnData, config: Dict) -> anndata.AnnData:
    """Scale data and run PCA."""
    # Scale
    sc.pp.scale(adata, max_value=10)
    logger.info("Scaled data")

    # PCA
    n_pcs = config.get('pca_n_comps', 50)
    sc.tl.pca(adata, n_comps=n_pcs)
    logger.info(f"Computed PCA with {n_pcs} components")

    return adata


def generate_postfilter_plots(
    adata: anndata.AnnData,
    sample_id: str,
    plots_dir: str
) -> None:
    """Generate post-filtering QC plots."""
    plots_dir = Path(plots_dir)
    plots_dir.mkdir(parents=True, exist_ok=True)

    try:
        # Histogram plots after filtering
        fig, axes = plt.subplots(1, 3, figsize=(15, 4))
        fig.suptitle(f"{sample_id} - Post-filtering QC Metrics")

        axes[0].hist(adata.obs['n_genes_by_counts'], bins=50, color='steelblue', edgecolor='black')
        axes[0].set_xlabel('Genes per Cell')
        axes[0].set_ylabel('Count')
        axes[0].set_title('n_genes_by_counts')

        axes[1].hist(adata.obs['total_counts'], bins=50, color='steelblue', edgecolor='black')
        axes[1].set_xlabel('Total Counts')
        axes[1].set_ylabel('Count')
        axes[1].set_title('total_counts')
        axes[1].set_yscale('log')

        axes[2].hist(adata.obs['pct_counts_mt'], bins=50, color='steelblue', edgecolor='black')
        axes[2].set_xlabel('% Mitochondrial Counts')
        axes[2].set_ylabel('Count')
        axes[2].set_title('pct_counts_mt')

        plt.tight_layout()
        output_file = plots_dir / f"{sample_id}_postfilter_histograms.png"
        plt.savefig(output_file, dpi=100, bbox_inches='tight')
        plt.close()
        logger.info(f"Saved postfilter histograms to {output_file}")

        # PCA variance explained
        if 'pca' in adata.obsm:
            fig, ax = plt.subplots(figsize=(10, 6))
            variance_ratio = adata.uns['pca']['variance_ratio']
            cumsum_var = np.cumsum(variance_ratio)
            ax.plot(range(1, len(variance_ratio) + 1), cumsum_var, 'b-', marker='o', markersize=4)
            ax.set_xlabel('PC')
            ax.set_ylabel('Cumulative Variance Explained')
            ax.set_title(f"{sample_id} - PCA Variance Explained")
            ax.grid(True, alpha=0.3)
            plt.tight_layout()
            output_file = plots_dir / f"{sample_id}_pca_variance.png"
            plt.savefig(output_file, dpi=100, bbox_inches='tight')
            plt.close()
            logger.info(f"Saved PCA variance plot to {output_file}")

    except Exception as e:
        logger.error(f"Error generating postfilter plots: {e}")
        raise


def write_stats(
    adata: anndata.AnnData,
    sample_id: str,
    stats_file: str,
    config: Dict
) -> None:
    """Write QC statistics to file."""
    try:
        with open(stats_file, 'w') as f:
            f.write(f"QC Statistics for {sample_id}\n")
            f.write("=" * 50 + "\n\n")

            f.write(f"Number of cells: {adata.n_obs}\n")
            f.write(f"Number of genes: {adata.n_vars}\n\n")

            f.write("Cell metrics:\n")
            f.write(f"  Mean genes per cell: {adata.obs['n_genes_by_counts'].mean():.1f}\n")
            f.write(f"  Median genes per cell: {adata.obs['n_genes_by_counts'].median():.1f}\n")
            f.write(f"  Mean total counts: {adata.obs['total_counts'].mean():.1f}\n")
            f.write(f"  Median total counts: {adata.obs['total_counts'].median():.1f}\n")
            f.write(f"  Mean % MT: {adata.obs['pct_counts_mt'].mean():.2f}%\n")
            f.write(f"  Median % MT: {adata.obs['pct_counts_mt'].median():.2f}%\n\n")

            if 'doublet_score' in adata.obs:
                n_doublets = (adata.obs['doublet_class'] == 'doublet').sum()
                f.write(f"Doublets detected: {n_doublets} ({100*n_doublets/len(adata):.1f}%)\n\n")

            f.write("Gene metrics:\n")
            f.write(f"  Mean counts per gene: {adata.var['total_counts'].mean():.1f}\n")
            f.write(f"  Median counts per gene: {adata.var['total_counts'].median():.1f}\n")
            f.write(f"  Genes with >0 cells: {(adata.var['n_cells_by_counts'] > 0).sum()}\n\n")

            if 'highly_variable' in adata.var:
                n_hvgs = adata.var['highly_variable'].sum()
                f.write(f"Highly variable genes: {n_hvgs}\n\n")

            if 'pca' in adata.obsm:
                variance_ratio = adata.uns['pca']['variance_ratio']
                cumsum_var = np.cumsum(variance_ratio)
                f.write(f"PCA cumulative variance (PC1-10): {cumsum_var[9]:.3f}\n")
                f.write(f"PCA cumulative variance (PC1-50): {cumsum_var[-1]:.3f}\n")

        logger.info(f"Wrote statistics to {stats_file}")

    except Exception as e:
        logger.error(f"Error writing statistics: {e}")
        raise


def main():
    parser = argparse.ArgumentParser(
        description='Scanpy QC and filtering for single-cell RNA-seq'
    )
    parser.add_argument('--input', required=True, help='Input H5AD file')
    parser.add_argument('--output', required=True, help='Output H5AD file')
    parser.add_argument('--sample-id', required=True, help='Sample ID')
    parser.add_argument('--config', required=True, help='QC config JSON file')
    parser.add_argument('--metadata', default='{}', help='Metadata JSON string')
    parser.add_argument('--stats-file', help='Output statistics file')
    parser.add_argument('--plots-dir', default='.', help='Output directory for plots')

    args = parser.parse_args()

    try:
        # Load configuration and data
        config = load_qc_config(args.config)
        adata = load_anndata(args.input)

        # Parse metadata
        metadata = json.loads(args.metadata)

        # Add metadata
        adata = add_metadata(adata, args.sample_id, metadata)

        # Flag mitochondrial genes
        mt_prefix = config.get('mt_prefix', 'MT-')
        adata = flag_mitochondrial_genes(adata, mt_prefix)

        # Calculate QC metrics
        adata = calculate_qc_metrics(adata)

        # Generate pre-filter plots
        generate_prefilter_plots(adata, args.sample_id, args.plots_dir)

        # Apply filters
        adata = apply_filters(adata, config)

        # Doublet detection
        if config.get('run_doublets', True):
            adata = run_doublet_detection(adata)

        # Normalize and select HVGs
        adata = normalize_and_hvg(adata, config)

        # Scale and run PCA
        adata = scale_and_pca(adata, config)

        # Generate post-filter plots
        generate_postfilter_plots(adata, args.sample_id, args.plots_dir)

        # Write statistics
        if args.stats_file:
            write_stats(adata, args.sample_id, args.stats_file, config)

        # Save output
        adata.write_h5ad(args.output, compression='gzip')
        logger.info(f"Saved filtered AnnData to {args.output}")

        logger.info("QC processing completed successfully")
        sys.exit(0)

    except Exception as e:
        logger.error(f"Fatal error: {e}")
        logger.debug(traceback.format_exc())
        sys.exit(1)


if __name__ == '__main__':
    main()
