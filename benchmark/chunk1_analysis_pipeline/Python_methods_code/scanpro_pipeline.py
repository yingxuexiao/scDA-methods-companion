#scanpro_pipeline:
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
import gc

# ================== CPU Limitations ==================
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
torch.set_num_threads(1)

# Bind the current process to CPU core 17
p = psutil.Process(os.getpid())
p.cpu_affinity([17])
print(f"CPU core bound: {p.cpu_affinity()}")

base_dir = Path("/Data/test/xiaoyingxue/").resolve()
data_dir = base_dir / "simulation_data"

# ================== Memory Management Functions ==================
def cleanup_memory():
    """Clean up memory thoroughly to ensure each dataset starts from a clean state"""
    # Python garbage collection
    gc.collect()
    
    # Clear numpy cache
    try:
        np._globals._NoValue = None
    except:
        pass
    
    # Clear PyTorch cache (if used)
    try:
        torch.cuda.empty_cache() if torch.cuda.is_available() else None
    except:
        pass
    
    # Clear matplotlib cache
    plt.close('all')
    
    # Wait for memory to stabilize
    time.sleep(1.0)

def get_stable_memory_usage(process, samples=3, interval=0.3):
    """Get stable memory usage measurement"""
    measurements = []
    for _ in range(samples):
        measurements.append(process.memory_info().rss / 1024 ** 2)
        time.sleep(interval)
    # Use median to reduce the impact of outliers
    return np.median(measurements)

# ================== Get All .h5ad Files ==================
def get_all_h5ad_files(data_dir: Path):
    """Recursively search for all .h5ad files in subfolders"""
    return list(data_dir.rglob("*.h5ad"))

# ================== Run Scanpro Method ==================
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
            print(f"Warning: Failed to generate props plot: {e}")

        try:
            fig_strip = results.plot(kind="boxplot")
            if fig_strip:
                figures["boxplot"] = fig_strip
        except Exception as e:
            print(f"Warning: Failed to generate boxplot: {e}")

    return result_df, figures

# ================== Process Single Dataset ==================
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

    # Thoroughly clean memory before processing
    cleanup_memory()
    
    # Get process object
    process = psutil.Process(os.getpid())
    
    # Get stable starting memory
    start_mem = get_stable_memory_usage(process, samples=3, interval=0.3)
    
    # Record start time
    start_time = time.time()
    
    # Log initialization (placed after memory measurement to avoid interference)
    logger = setup_logger(log_path)
    log_system_info(logger)
    log_loaded_packages(logger)
    logger.info(f"Initial memory: {start_mem:.2f} MB")
    
    # Load the dataset
    adata = ad.read_h5ad(data_file)
    logger.info(f"Data loaded, number of cells: {adata.n_obs}, number of genes: {adata.n_vars}")
    
    logger.info(f"Running Scanpro...")

    # Run Scanpro method
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

    # Clean up memory after processing
    cleanup_memory()
    
    # Get stable ending memory
    end_mem = get_stable_memory_usage(process, samples=3, interval=0.3)
    
    # Calculate total runtime
    end_time = time.time()
    runtime = end_time - start_time
    
    # Calculate memory usage (ensure non-negative)
    mem_used = max(end_mem - start_mem, 0)
    
    # If memory change is small, treat it as noise
    if mem_used < 2.0:  # Changes smaller than 2MB are considered noise
        mem_used = 0
    
    # Log results
    logger.info(f"Scanpro method runtime: {runtime:.2f} seconds")
    logger.info(f"Ending memory: {end_mem:.2f} MB")
    logger.info(f"Memory used: {mem_used:.2f} MB")

    # === Save Results ===
    results_path = output_dir / "scanpro_results.csv"
    da_results.to_csv(results_path)
    logger.info(f"Results saved to: {results_path}")

    # Save figures
    for name, fig in figures.items():
        fig_path = figures_dir / f"{name}.png"
        fig.savefig(fig_path, dpi=300, bbox_inches='tight')
        plt.close(fig)  # Release memory
        logger.info(f"Figure saved to: {fig_path}")

    # Clean up adata object
    del adata
    
    # Clean up memory again for the next dataset
    cleanup_memory()
    
    # Return time and memory
    return dataset_name, runtime, mem_used

# ================== Main Program ==================
if __name__ == "__main__":
    # Clean up memory at the start of the program
    cleanup_memory()
    
    summary_list = []

    # Get all .h5ad files
    h5ad_files = get_all_h5ad_files(data_dir)
    print(f"Found {len(h5ad_files)} .h5ad files")
    
    # Sort datasets by size and process them
    h5ad_files_info = []
    for file in h5ad_files:
        try:
            file_size_mb = os.path.getsize(file) / (1024 * 1024)
            h5ad_files_info.append((file, file_size_mb))
        except:
            h5ad_files_info.append((file, 0))
    
    # Sort by file size (from small to large)
    h5ad_files_info.sort(key=lambda x: x[1])
    h5ad_files = [file for file, _ in h5ad_files_info]
    
    print("Dataset processing order (from small to large):")
    for i, (file, size_mb) in enumerate(h5ad_files_info):
        print(f"  {i+1}. {file.stem}: {size_mb:.1f} MB")
    
    print("=" * 50)
    
    for data_file in h5ad_files:
        try:
            print(f"\nProcessing dataset: {data_file.stem}")
            dataset_name, runtime, mem_used = process_single_dataset(data_file)
            summary_list.append({
                "dataset": dataset_name,
                "runtime_sec": round(runtime, 4),
                "memory_MB": round(mem_used, 4)
            })
            print(f"✓ Completed - Time: {runtime:.2f} seconds, Memory: {mem_used:.2f} MB")
        except Exception as e:
            print(f"❌ Failed to process dataset {data_file.stem}: {e}")
            summary_list.append({
                "dataset": data_file.stem,
                "runtime_sec": None,
                "memory_MB": None
            })
        
        # Add pause between datasets
        print("-" * 50)
        time.sleep(1.0)

    # === Save Summary Table ===
    if summary_list:
        summary_df = pd.DataFrame(summary_list)
        summary_path = base_dir / "da_results" / "scanpro_summary.csv"
        summary_df.to_csv(summary_path, index=False)
        print(f"\n✅ All datasets processed, summary table saved to: {summary_path}")
        
        # Print summary statistics
        successful_runs = summary_df[summary_df['memory_MB'].notna()]
        if len(successful_runs) > 0:
            print("\n" + "="*50)
            print("Run Summary Statistics:")
            print("="*50)
            print(f"Number of successfully processed datasets: {len(successful_runs)}")
            print(f"Average runtime: {successful_runs['runtime_sec'].mean():.2f} seconds")
            print(f"Average memory usage: {successful_runs['memory_MB'].mean():.2f} MB")
            print(f"Maximum memory usage: {successful_runs['memory_MB'].max():.2f} MB")
            print(f"Minimum memory usage: {successful_runs['memory_MB'].min():.2f} MB")
            print(f"Total runtime: {successful_runs['runtime_sec'].sum():.2f} seconds")
            print("="*50)
    else:
        print("\n⚠️ No successfully processed datasets")
