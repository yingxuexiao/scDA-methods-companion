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

# ================== CPU LIMITATIONS ==================
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
torch.set_num_threads(1)

# Bind the current process to CPU 0
p = psutil.Process(os.getpid())
p.cpu_affinity([0])
print(f"CPU core bound: {p.cpu_affinity()}")

# Clear all handlers from the root logger
root_logger = logging.getLogger()
if root_logger.hasHandlers():
    root_logger.handlers.clear()
root_logger.setLevel(logging.WARNING)  # Only show warnings and above, reduce distractions

# ====== Global summary list ======
summary_records = []

# ====== Base path ======
base_dir = Path('/data1_Literature/xiao_DA/da_project').resolve()
data_dir = base_dir / "preprocessed_data"

# ================== Get all .h5ad files ==================
def get_all_h5ad_files(data_dir: Path):
    """Recursively search for all .h5ad files in subfolders"""
    return list(data_dir.rglob("*.h5ad"))

# ====== Dataset processing function ======
def get_dataset_name(data_file_path: Path) -> str:
    if data_file_path.suffix == ".h5ad":
        return data_file_path.stem
    else:
        raise ValueError("Only .h5ad files are supported as input datasets")

# ====== SCCODA analysis function ======
def run_sccoda_analysis(
    adata,
    celltype_col="celltype",
    sample_col="sample",
    condition_col="condition"
):
    process = psutil.Process(os.getpid())
    start_mem = process.memory_info().rss / 1024 ** 2  # MB
    start_time = time.time()

    try:
        print(f"[INFO] Reading dataset, number of cells: {adata.n_obs}, number of genes: {adata.n_vars}")

        # Check the number of cell types and their categories
        num_cell_types = adata.obs[celltype_col].nunique()
        print(f"[INFO] There are {num_cell_types} cell types")
        print(f"[INFO] Example cell types: {adata.obs[celltype_col].unique()}")

        # Check the relationship between sample and condition
        relation = pd.crosstab(adata.obs[sample_col], adata.obs[condition_col])
        print("[INFO] Sample-condition relationship:")
        print(relation)

        # Initialize scCODA model
        sccoda_model = pt.tl.Sccoda()
        print("[INFO] scCODA model initialization complete")

        # Generate sample-level data
        sccoda_data = sccoda_model.load(
            adata,
            type="cell_level",
            generate_sample_level=True,
            cell_type_identifier=celltype_col,
            sample_identifier=sample_col,
            covariate_obs=[condition_col],
        )
        print("[INFO] Sample-level data generation complete")

        # Extract all conditions
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

        # Save a copy of the default results
        sccoda_data.mod["coda_salm_fdr_0.05"] = sccoda_data.mod[subset_key].copy()

        # Change FDR threshold to 0.2
        sccoda_model.set_fdr(sccoda_data, modality_key=subset_key, est_fdr=0.2)
        sccoda_model.summary(sccoda_data, modality_key=subset_key)

        return sccoda_data

    except Exception as e:
        print("[ERROR] SCCODA run failed!")
        traceback.print_exc()
        return None

# ====== Single dataset processing function ======
def process_single_dataset(data_file: Path):
    dataset_name = get_dataset_name(data_file)
    print(f"[INFO] Processing dataset: {dataset_name}")

    output_dir = base_dir / "da_results" / dataset_name
    output_dir.mkdir(parents=True, exist_ok=True)

    log_dir = base_dir / "log"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / f"{dataset_name}_sccoda.log"

    # Read data
    adata = ad.read_h5ad(data_file)

    # Log initialization
    logger = setup_logger(log_path)
    log_system_info(logger)
    log_loaded_packages(logger)

    # Record initial memory and time
    process = psutil.Process(os.getpid())
    start_mem = process.memory_info().rss / 1024 ** 2  # MB
    start_time = time.time()

    status = "success"
    try:
        sccoda_data = run_sccoda_analysis(adata)
        if sccoda_data is None:
            status = "failed"
    except Exception as e:
        logger.error(f"SCCODA run failed: {e}")
        status = "failed"
        sccoda_data = None

    # Record total time and memory usage
    elapsed = time.time() - start_time
    end_mem = process.memory_info().rss / 1024 ** 2
    mem_used = end_mem - start_mem

    # Update summary table (handled in the outer layer)
    summary_records.append({
        "dataset": dataset_name,
        "elapsed_time_sec": elapsed if status=="success" else None,
        "memory_MB": mem_used if status=="success" else None,
        "status": status
    })

    # Save results
    if sccoda_data is not None:
        sccoda_data.write_h5mu(output_dir / f"{dataset_name}_sccoda.h5mu")

# ====== Main program entry ======
if __name__ == "__main__":
    # Get all .h5ad files
    h5ad_files = get_all_h5ad_files(data_dir)
    print(f"Found {len(h5ad_files)} .h5ad files: {h5ad_files}")

    for data_file in h5ad_files:
        try:
            process_single_dataset(data_file)
        except Exception as e:
            print(f"❌ Processing dataset {data_file.stem} failed: {e}")
            summary_records.append({
                "dataset": data_file.stem,
                "elapsed_time_sec": None,
                "memory_MB": None,
                "status": "failed"
            })

    # Save summary table
    summary_df = pd.DataFrame(summary_records)
    summary_path = base_dir / "da_results" / "sccoda_summary.csv"
    summary_df.to_csv(summary_path, index=False)
    print(f"[INFO] All dataset run summary saved to: {summary_path}")
