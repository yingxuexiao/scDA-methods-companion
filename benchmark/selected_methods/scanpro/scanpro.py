import sys
import pandas as pd
from scanpro import scanpro
from scanpro.utils import convert_counts_to_df
from scanpy import AnnData
import matplotlib.pyplot as plt
from typing import Optional, List, Tuple, Dict
from python_utils import setup_logger, log_system_info, log_loaded_packages, run_analysis_and_log, save_results
from pathlib import Path
import os
import time
import psutil
import torch
import platform
import logging
import numpy as np
import importlib.metadata
import anndata as ad
from anndata import AnnData
from datetime import datetime
from matplotlib.figure import Figure
import torch

# ================== CPU Limits ==================
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
torch.set_num_threads(1)

# Bind current process to CPU 0
p = psutil.Process(os.getpid())
p.cpu_affinity([17])
print(f"Bound to CPU core: {p.cpu_affinity()}")

base_dir = Path('/data1_Literature/xiao_DA/da_project').resolve()
data_dir = base_dir / "preprocessed_data"

# ================== Get all .h5ad files ==================
def get_all_h5ad_files(data_dir: Path):
    """Recursively search for all .h5ad files in subdirectories"""
    return list(data_dir.rglob("*.h5ad"))

# ================== Execute Scanpro method ==================
def run_scanpro_method(
    adata,
    clusters_col: str,
    conds_col: str,
    samples_col: Optional[str] = None,
    covariates: Optional[List[str]] = None,
    transform: str = "logit",
    n_sims: int = 100,
    n_reps: int = 8,
    plot: bool = True,
    specific_clusters: Optional[List[str]] = None,
    n_cols: int = 4,
) -> Tuple[pd.DataFrame, Optional[Dict[str, plt.Figure]]]:

    df = adata.obs.copy()

    print("Running Scanpro...")
    if samples_col:
        results = scanpro(
            df,
            samples_col=samples_col,
            clusters_col=clusters_col,
            conds_col=conds_col,
            covariates=covariates,
            transform=transform,
        )
    else:
        results = scanpro(
            df,
            clusters_col=clusters_col,
            conds_col=conds_col,
            covariates=covariates,
            transform=transform,
            n_sims=n_sims,
            n_reps=n_reps,
        )

    result_df = results.results
    figures = {}

    if plot:
        try:
            fig_props = results.props
            if isinstance(fig_props, plt.Figure):
                figures["props"] = fig_props
        except Exception as e:
            print(f"Warning: props image failed: {e}")

        try:
            fig_strip = results.plot(kind="boxplot")
            if fig_strip:
                figures["boxplot"] = fig_strip
        except Exception as e:
            print(f"Warning: boxplot failed: {e}")

    return result_df, figures

# ================== Process a single dataset ==================
def get_dataset_name(data_file_path: Path) -> str:
    if data_file_path.suffix == ".h5ad":
        return data_file_path.stem
    else:
        raise ValueError("Only .h5ad files are supported as input datasets")

def process_single_dataset(data_file: Path):
    dataset_name = get_dataset_name(data_file)
    print(f"Processing dataset: {dataset_name}")

    output_dir = base_dir / "da_results" / dataset_name
    output_dir.mkdir(parents=True, exist_ok=True)

    log_dir = base_dir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / f"{dataset_name}_scanpro.log"

    figures_dir = base_dir / "figures" / dataset_name
    figures_dir.mkdir(parents=True, exist_ok=True)

    # Read data
    adata = ad.read_h5ad(data_file)

    # Log initialization
    logger = setup_logger(log_path)
    log_system_info(logger)
    log_loaded_packages(logger)

    # === Execute method + memory tracking ===
    process = psutil.Process(os.getpid())
    start_mem = process.memory_info().rss / 1024 ** 2  # MB
    start_time = time.time()
    logger.info(f"Starting Scanpro, initial memory: {start_mem:.2f} MB")

    da_results, figures = run_scanpro_method(
        adata,
        clusters_col="celltype",
        conds_col="condition",
        samples_col="sample",
        covariates=None,
        transform="logit",
        n_sims=100,
        n_reps=8,
        plot=True
    )

    end_time = time.time()
    end_mem = process.memory_info().rss / 1024 ** 2  # MB
    logger.info(f"Scanpro method runtime: {end_time - start_time:.2f} seconds")
    logger.info(f"End memory: {end_mem:.2f} MB, memory used: {end_mem - start_mem:.2f} MB")

    # === Save results ===
    results_path = output_dir / "scanpro_results.csv"
    da_results.to_csv(results_path)
    logger.info(f"Results saved to: {results_path}")

    for name, fig in figures.items():
        fig_path = figures_dir / f"{name}.png"
        fig.savefig(fig_path)
        plt.close(fig)  # Release memory
        logger.info(f"Image saved to: {fig_path}")

    # Return time and memory for summary table
    return dataset_name, end_time - start_time, end_mem - start_mem

# ================== Main Program ==================
if __name__ == "__main__":
    summary_list = []

    # Get all .h5ad files
    h5ad_files = get_all_h5ad_files(data_dir)
    print(f"Found {len(h5ad_files)} .h5ad files: {h5ad_files}")

    for data_file in h5ad_files:
        try:
            dataset_name, runtime, mem_used = process_single_dataset(data_file)
            summary_list.append({
                "dataset": dataset_name,
                "runtime_sec": runtime,
                "memory_MB": mem_used
            })
        except Exception as e:
            print(f"❌ Failed to process dataset {data_file.stem}: {e}")

    # === Save summary table ===
    summary_df = pd.DataFrame(summary_list)
    summary_path = base_dir / "da_results" / "scanpro_summary.csv"
    summary_df.to_csv(summary_path, index=False)
    print(f"✅ All datasets processed, summary table saved to: {summary_path}")
