#scCODA_pipeline.py

import sys
import os
import time
import pickle as pkl
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import pertpy as pt
import anndata as ad
from pathlib import Path
from python_utils import setup_logger, log_system_info, log_loaded_packages
from typing import Tuple
import psutil
import logging
import mudata as mu
import traceback
import torch
import gc
import numpy as np

# ================== CPU Limitations ==================
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
torch.set_num_threads(1)

# Bind the current process to CPU core 0
p = psutil.Process(os.getpid())
p.cpu_affinity([0])
print(f"CPU core bound: {p.cpu_affinity()}")

# Clear all handlers for the root logger
root_logger = logging.getLogger()
if root_logger.hasHandlers():
    root_logger.handlers.clear()
root_logger.setLevel(logging.WARNING)  # Only show warnings and above to reduce distractions

# ====== Global Summary List ======
summary_records = []

# ====== Base Paths ======
base_dir = Path("/Data/test/xiaoyingxue").resolve()
data_dir = base_dir / "simulation_data"

# ================== Memory Cleanup Functions ==================
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

# ====== Dataset Processing Function ======
def get_dataset_name(data_file_path: Path) -> str:
    if data_file_path.suffix == ".h5ad":
        return data_file_path.stem
    else:
        raise ValueError("Only .h5ad files are supported as input datasets")

# ====== SCCODA Analysis Function ======
def run_sccoda_analysis(
    adata,
    celltype_col="celltype",
    sample_col="sample",
    condition_col="condition"
):
    try:
        print(f"[INFO] Reading dataset, number of cells: {adata.n_obs}, number of genes: {adata.n_vars}")

        # Check number and types of cell types
        num_cell_types = adata.obs[celltype_col].nunique()
        print(f"[INFO] Found {num_cell_types} cell types")
        print(f"[INFO] Cell types: {adata.obs[celltype_col].unique()}")

        # Check relationship between sample and condition
        relation = pd.crosstab(adata.obs[sample_col], adata.obs[condition_col])
        print("[INFO] Sample-condition relationship:")
        print(relation)

        # Initialize scCODA model
        sccoda_model = pt.tl.Sccoda()
        print("[INFO] Initialized scCODA model")

        # Generate sample-level data
        sccoda_data = sccoda_model.load(
            adata,
            type="cell_level",
            generate_sample_level=True,
            cell_type_identifier=celltype_col,
            sample_identifier=sample_col,
            covariate_obs=[condition_col],
        )
        print("[INFO] Generated sample-level data")

        # Extract unique conditions
        unique_conditions = sccoda_data["coda"].obs[condition_col].unique().tolist()
        print(f"[INFO] Detected conditions: {unique_conditions}")

        subset_key = "coda_salm"
        sccoda_data.mod[subset_key] = sccoda_data["coda"][
            sccoda_data["coda"].obs[condition_col].isin(unique_conditions)
        ].copy()
        print("[INFO] Subset data preparation complete")

        # Automatically select reference cell type (with minimum coefficient of variation)
        X = sccoda_data["coda"].X
        cv = X.std(axis=0) / X.mean(axis=0)
        ref_cell_type = sccoda_data["coda"].var.index[cv.argmin()]
        print(f"[INFO] Automatically selected reference cell type: {ref_cell_type}")

        # Prepare the model
        sccoda_data = sccoda_model.prepare(
            sccoda_data,
            modality_key=subset_key,
            formula=condition_col,
            reference_cell_type=ref_cell_type,
        )
        print("[INFO] Model preparation complete")

        # Run MCMC sampling
        print("[INFO] Starting MCMC sampling...")
        mcmc_start = time.time()
        sccoda_model.run_nuts(sccoda_data, modality_key=subset_key)
        mcmc_end = time.time()
        print(f"[INFO] MCMC sampling complete, took {mcmc_end - mcmc_start:.2f} seconds")

        # Default FDR 0.05
        sccoda_model.summary(sccoda_data, modality_key=subset_key)
        sccoda_model.credible_effects(sccoda_data, modality_key=subset_key)

        # Save default results copy
        sccoda_data.mod["coda_salm_fdr_0.05"] = sccoda_data.mod[subset_key].copy()

        # Change FDR threshold to 0.2
        sccoda_model.set_fdr(sccoda_data, modality_key=subset_key, est_fdr=0.2)
        sccoda_model.summary(sccoda_data, modality_key=subset_key)

        return sccoda_data

    except Exception as e:
        print("[ERROR] scCODA analysis failed!")
        traceback.print_exc()
        return None

# ====== Single Dataset Processing Function ======
def process_single_dataset(data_file: Path):
    dataset_name = get_dataset_name(data_file)
    print(f"[INFO] Starting processing dataset: {dataset_name}")

    output_dir = base_dir / "da_results" / dataset_name
    output_dir.mkdir(parents=True, exist_ok=True)

    log_dir = base_dir / "log"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / f"{dataset_name}_sccoda.log"

    # Thoroughly clean memory before reading data
    cleanup_memory()
    
    process = psutil.Process(os.getpid())
    
    # Get stable starting memory (multiple samples)
    start_mem = get_stable_memory_usage(process, samples=3, interval=0.3)
    
    # Record start time
    start_time = time.time()
    
    # Load data
    adata = ad.read_h5ad(data_file)

    # Log initialization
    logger = setup_logger(log_path)
    log_system_info(logger)
    log_loaded_packages(logger)
    logger.info(f"Initial memory: {start_mem:.2f} MB")

    status = "success"
    try:
        sccoda_data = run_sccoda_analysis(adata)
        if sccoda_data is None:
            status = "failed"
    except Exception as e:
        logger.error(f"scCODA analysis failed: {e}")
        status = "failed"
        sccoda_data = None

    # Clean up memory after running
    cleanup_memory()
    
    # Get stable ending memory
    end_mem = get_stable_memory_usage(process, samples=3, interval=0.3)
    
    # Calculate total time and memory usage
    elapsed = time.time() - start_time
    mem_used = max(end_mem - start_mem, 0)  # Ensure non-negative
    
    # If memory change is very small, treat as noise
    if mem_used < 2.0:  # Changes smaller than 2MB are considered noise
        mem_used = 0
    
    logger.info(f"Ending memory: {end_mem:.2f} MB")
    logger.info(f"Memory used: {mem_used:.2f} MB")
    logger.info(f"Elapsed time: {elapsed:.2f} seconds")

    # Update summary table
    summary_records.append({
        "dataset": dataset_name,
        "elapsed_time_sec": round(elapsed, 4) if status == "success" else None,
        "memory_MB": round(mem_used, 4) if status == "success" else None,
        "status": status
    })

    # Save results
    if sccoda_data is not None:
        sccoda_data.write_h5mu(output_dir / f"{dataset_name}_sccoda.h5mu")
        logger.info(f"Results saved to: {output_dir / f'{dataset_name}_sccoda.h5mu'}")
    
    # Finally, clean up for the next dataset
    cleanup_memory()

# ====== Main Program Entry ======
if __name__ == "__main__":
    # Clean memory at the start of the program
    cleanup_memory()
    
    # Get all .h5ad files
    h5ad_files = get_all_h5ad_files(data_dir)
    print(f"Found {len(h5ad_files)} .h5ad files")
    
    # Sort datasets by size (optional), to observe memory usage patterns
    h5ad_files_info = []
    for file in h5ad_files:
        try:
            # Quickly check file size
            file_size_mb = os.path.getsize(file) / (1024 * 1024)
            h5ad_files_info.append((file, file_size_mb))
        except:
            h5ad_files_info.append((file, 0))
    
    # Sort by file size (smallest to largest)
    h5ad_files_info.sort(key=lambda x: x[1])
    h5ad_files = [file for file, _ in h5ad_files_info]
    
    print(f"Datasets sorted by size (smallest to largest):")
    for i, (file, size_mb) in enumerate(h5ad_files_info):
        print(f"  {i+1}. {file.stem}: {size_mb:.1f} MB")
    
    for data_file in h5ad_files:
        try:
            process_single_dataset(data_file)
        except Exception as e:
            print(f"❌ Failed to process dataset {data_file.stem}: {e}")
            summary_records.append({
                "dataset": data_file.stem,
                "elapsed_time_sec": None,
                "memory_MB": None,
                "status": "failed"
            })
        
        # Add short pause between datasets
        print("-" * 50)
        time.sleep(1.0)

    # Save summary table
    summary_df = pd.DataFrame(summary_records)
    summary_path = base_dir / "da_results" / "1216_sccoda_summary.csv"
    summary_df.to_csv(summary_path, index=False)
    print(f"[INFO] All datasets' summary table saved to: {summary_path}")
    
    # Print summary statistics
    print("\n" + "="*50)
    print("Memory Usage Statistics:")
    print("="*50)
    successful_runs = summary_df[summary_df['status'] == 'success']
    if len(successful_runs) > 0:
        print(f"Number of successfully run datasets: {len(successful_runs)}")
        print(f"Average memory usage: {successful_runs['memory_MB'].mean():.2f} MB")
        print(f"Maximum memory usage: {successful_runs['memory_MB'].max():.2f} MB")
        print(f"Minimum memory usage: {successful_runs['memory_MB'].min():.2f} MB")
        print(f"Total runtime: {successful_runs['elapsed_time_sec'].sum():.2f} seconds")
    else:
        print("No successfully run datasets")
