#MELD_pipeline.py
import os
import psutil
import torch
import time
import platform
from datetime import datetime
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import mudata as md
import graphtools as gt
import meld
import phate
import magic
import scprep
import cmocean
import sklearn
import scipy
import seaborn as sns
import matplotlib.pyplot as plt
from typing import Optional
from utils import setup_logger, log_system_info, log_loaded_packages
from pathlib import Path
from scipy.stats import chi2_contingency
import tempfile
import gc

# ================== CPU Limitation ==================
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
torch.set_num_threads(1)

# Bind the current process to CPU core 18
p = psutil.Process(os.getpid())
p.cpu_affinity([18])
print(f"Bound to CPU core: {p.cpu_affinity()}")

# ================== Cache Directories ==================
cache_dir = Path("/Data/test/xiaoyingxue/")
cache_dir.mkdir(parents=True, exist_ok=True)
os.environ["TMPDIR"] = f"{cache_dir}/tmp"
os.environ["JOBLIB_TEMP_FOLDER"] = f"{cache_dir}/joblib"
os.environ["TRANSFORMERS_CACHE"] = f"{cache_dir}/transformers"
os.environ["TORCH_HOME"] = f"{cache_dir}/torch"
os.environ["XDG_CACHE_HOME"] = f"{cache_dir}"

print("Current temporary directory:", tempfile.gettempdir())

# ================== Base Paths ==================
base_dir = Path("/Data/test/xiaoyingxue/").resolve()
data_dir = base_dir / "simulation_data"  # Root directory for preprocessed data
log_dir = base_dir / "logs"
results_dir = base_dir / "da_results"

log_dir.mkdir(parents=True, exist_ok=True)
results_dir.mkdir(parents=True, exist_ok=True)

# ================== Recursive File Reading ==================
def get_all_h5ad_files(data_dir):
    """Recursively search for all .h5ad files in subdirectories"""
    return list(data_dir.rglob("*.h5ad"))

# ================== MELD Core Function ==================
def runMELD(adata, k, sample_col, label_col, replicate_col=None, beta=20):
    """Run MELD analysis"""
    samplem = pd.DataFrame(index=pd.Series(adata.obs[sample_col]).unique())
    samplem.loc[:, label_col] = (
        adata.obs[[sample_col, label_col]]
        .groupby(by=sample_col)
        .aggregate(lambda x: x.iloc[0])
    )
    if replicate_col is not None:
        samplem.loc[:, replicate_col] = (
            adata.obs[[sample_col, replicate_col]]
            .groupby(by=sample_col)
            .aggregate(lambda x: x.iloc[0])
        )
    adata.uns['samplem'] = samplem

    if adata.n_vars <= 50:
        G = gt.Graph(adata.X, knn=k, use_pygsp=True)
    else:
        if 'X_pca' not in adata.obsm:
            sc.tl.pca(adata, n_comps=50)
        G = gt.Graph(adata.obsm['X_pca'], knn=k, use_pygsp=True)

    meld_op = meld.MELD(beta=beta)

    if replicate_col:
        sample_densities = meld_op.fit_transform(G, sample_labels=adata.obs[sample_col])
        sample_likelihoods = sample_densities.div(sample_densities.sum(axis=1), axis=0)
        sample_likelihoods.columns = sample_likelihoods.columns.map(samplem[label_col].to_dict())
        likelihoods_per_cond = sample_likelihoods.groupby(axis=1, level=0).mean()
    else:
        sample_densities = meld_op.fit_transform(G, sample_labels=adata.obs[label_col])
        likelihoods_per_cond = sample_densities.div(sample_densities.sum(axis=1), axis=0)

    adata.obsm['meld_likelihoods'] = likelihoods_per_cond
    return likelihoods_per_cond

# ================== Single Dataset Processing ==================
def get_dataset_name(data_file_path: Path) -> str:
    if data_file_path.suffix == ".h5ad":
        return data_file_path.stem
    else:
        raise ValueError("Only .h5ad files are supported as input datasets")

# ================== Single Dataset Processing ==================
def process_single_dataset(data_file: Path, celltype_col: str = 'celltype'):
    dataset_name = get_dataset_name(data_file)
    print(f"Processing dataset: {dataset_name}")

    output_dir = results_dir / dataset_name
    output_dir.mkdir(parents=True, exist_ok=True)

    log_path = log_dir / f"{dataset_name}_meld.log"
    figures_dir = base_dir / "figures" / dataset_name
    figures_dir.mkdir(parents=True, exist_ok=True)

    # Read data, avoid virtual files (those starting with ._)
    if data_file.suffix == ".h5ad" and not data_file.stem.startswith('._'):
        adata = ad.read_h5ad(data_file)
    else:
        print(f"Skipping unnecessary file: {data_file}")
        return dataset_name, 0, 0

    if celltype_col not in adata.obs.columns:
        print(f"❌ Dataset {dataset_name} does not have column {celltype_col}, skipping")
        return dataset_name, 0, 0

    # Log initialization
    logger = setup_logger(log_path)
    log_system_info(logger)
    log_loaded_packages(logger)

    # Memory & Time tracking
    process = psutil.Process(os.getpid())
    start_time = time.time()

    # Memory cleanup and recording
    gc.collect()  # Trigger garbage collection to clean up unused objects
    start_mem = process.memory_info().rss / 1024 ** 2  # Get initial memory usage
    logger.info(f"Starting MELD, initial memory: {start_mem:.2f} MB")

    # Run MELD
    likelihoods = runMELD(
        adata,
        k=7,
        sample_col='sample',
        label_col='condition',
        replicate_col=None,
        beta=20
    )

    # Record end time and memory usage
    end_time = time.time()
    gc.collect()  # Trigger garbage collection after MELD run
    end_mem = process.memory_info().rss / 1024 ** 2  # Get memory usage at the end
    runtime = end_time - start_time
    mem_used = max(end_mem - start_mem, 0)  # Ensure memory usage is not negative

    # Log recording
    logger.info(f"MELD method runtime: {runtime:.2f} seconds")
    logger.info(f"End memory: {end_mem:.2f} MB, memory consumed: {mem_used:.2f} MB")

    # Save predicted labels
    adata.obs['meld_pred_condition'] = likelihoods.idxmax(axis=1)

    # Generate count and proportion tables
    celltype_cond_counts = adata.obs.groupby([celltype_col, 'meld_pred_condition']).size().unstack(fill_value=0)
    celltype_cond_props = celltype_cond_counts.div(celltype_cond_counts.sum(axis=1), axis=0)
    adata.uns['celltype_cond_counts'] = celltype_cond_counts
    adata.uns['celltype_cond_props'] = celltype_cond_props

    # Chi-squared test
    chi2_results = {}
    for celltype, row in celltype_cond_counts.iterrows():
        chi2, p, dof, expected = chi2_contingency([row.values, celltype_cond_counts.sum().values])
        chi2_results[celltype] = p
    chi2_df = pd.DataFrame.from_dict(chi2_results, orient='index', columns=['p_value'])
    chi2_df['adj_p'] = chi2_df['p_value'] * len(chi2_df)
    adata.uns['chi2_results'] = chi2_df

    # Save CSV files
    method_name = "MELD"
    props_csv_path = output_dir / f"{dataset_name}_{method_name}_celltype_cond_props.csv"
    chi2_csv_path = output_dir / f"{dataset_name}_{method_name}_p_value_results.csv"
    celltype_cond_props.to_csv(props_csv_path)
    chi2_df.to_csv(chi2_csv_path)
    logger.info(f"Proportion table saved to: {props_csv_path}")
    logger.info(f"Chi-squared test results saved to: {chi2_csv_path}")

    # Save heatmap
    plt.figure(figsize=(16, 9))
    sns.heatmap(celltype_cond_props, annot=True, fmt=".2f", cmap="YlGnBu", cbar_kws={'label': 'Proportion'})
    plt.title(f"Cell type proportions across conditions ({method_name})", fontsize=14)
    plt.ylabel("Cell Type", fontsize=12)
    plt.xlabel("Condition", fontsize=12)
    plt.tight_layout()
    heatmap_pdf_path = figures_dir / f"{dataset_name}_{method_name}_celltype_condition_heatmap.pdf"
    plt.savefig(heatmap_pdf_path, format='pdf', dpi=300)
    plt.close()
    logger.info(f"Heatmap saved to: {heatmap_pdf_path}")

    # Save h5mu file
    mdata = md.MuData({"meld": adata})
    mdata.write_h5mu(output_dir / f"{dataset_name}_meld.h5mu")
    logger.info(f"h5mu file saved to {output_dir / f'{dataset_name}_meld.h5mu'}")

    return dataset_name, runtime, mem_used

# ================== Main Program ==================
if __name__ == "__main__":
    summary_list = []

    # Get all .h5ad files
    h5ad_files = get_all_h5ad_files(data_dir)
    print(f"Found {len(h5ad_files)} .h5ad files: {h5ad_files}")

    for data_file in h5ad_files:
        try:
            dataset_name, runtime, mem_used = process_single_dataset(data_file, celltype_col='celltype')
            summary_list.append({
                "dataset": dataset_name,
                "runtime_sec": runtime,
                "memory_MB": mem_used
            })
        except Exception as e:
            print(f"❌ Failed to process dataset {data_file.stem}: {e}")

    # Save summary table
    if summary_list:
        summary_df = pd.DataFrame(summary_list)
        summary_path = results_dir / "meld_summary.csv"
        summary_df.to_csv(summary_path, index=False)
        print(f"✅ All datasets processed, summary table saved to: {summary_path}")
