import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import anndata as ad
import cna
np.random.seed(0)
import time
from pathlib import Path
import logging
from python_utils import setup_logger, log_system_info, log_loaded_packages
import mudata as mu
import psutil
import os
import torch

# ================== CPU  ==================
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
torch.set_num_threads(1)

# Bind the current process to CPU 19
p = psutil.Process(os.getpid())
p.cpu_affinity([19])
print(f"CPU : {p.cpu_affinity()}")
# ====== 基本路径设置 ======
base_dir = Path('/data1_Literature/xiao_DA/da_project').resolve()
data_dir = base_dir / "preprocessed_data"
output_root = base_dir / "da_results"
log_root = base_dir / "log"
figures_root = base_dir / "figures"

# ====== Global summary list ======
summary_records = []

# ====== Utility function ======
def get_dataset_name(data_file_path: Path) -> str:
    if data_file_path.suffix == ".h5ad":
        return data_file_path.stem
    else:
        raise ValueError("Only.h5ad files are supported as input datasets")

def add_one_hot_encoding(adata, condition_col='condition', sample_col='sample', batch_col='batch'):
    meta = adata.obs.copy()
    condition_dummies = pd.get_dummies(meta[condition_col], prefix=condition_col)
    sample_dummies = pd.get_dummies(meta[sample_col], prefix=sample_col)
    dfs_to_concat = [adata.obs, condition_dummies, sample_dummies]

    if batch_col in meta.columns and meta[batch_col].notna().any():
        batch_dummies = pd.get_dummies(meta[batch_col], prefix=batch_col)
        dfs_to_concat.append(batch_dummies)
        print(f"Batch column '{batch_col}' detected and encoded.")
    else:
        print(f"No batch column '{batch_col}' detected or it's empty. Skip batch encoding.")

    adata.obs = pd.concat(dfs_to_concat, axis=1)
    return adata

def prepare_samplem_from_adata(adata, condition_prefix='condition_', batch_prefix='batch_', sample_col='sample'):
    obs = adata.obs.copy()
    # condition
    condition_cols = [col for col in obs.columns if col.startswith(condition_prefix)]
    if not condition_cols:
        if 'condition' in obs.columns:
            dummies = pd.get_dummies(obs['condition'], prefix=condition_prefix.rstrip('_'))
            obs = pd.concat([obs, dummies], axis=1)
            condition_cols = dummies.columns.tolist()
        else:
            raise ValueError("The condition column was not found, and there is no encoded condition_column")
    # batch
    batch_cols = [col for col in obs.columns if col.startswith(batch_prefix)]
    if not batch_cols and 'batch' in obs.columns:
        dummies_batch = pd.get_dummies(obs['batch'], prefix=batch_prefix.rstrip('_'))
        obs = pd.concat([obs, dummies_batch], axis=1)
        batch_cols = dummies_batch.columns.tolist()
    adata.obs = obs
    cols_for_samplem = condition_cols + batch_cols if batch_cols else condition_cols
    samplem = cna.ut.obs_to_sample(adata, cols_for_samplem, sample_col)
    return samplem

def preprocess_neighbors(adata, n_neighbors=30, use_rep='X_pca'):
    sc.pp.neighbors(adata, use_rep=use_rep, n_neighbors=n_neighbors)

def run_all_condition_associations(adata, samplem, sample_col='sample', prefix='condition_'):
    condition_cols = [col for col in samplem.columns if col.startswith(prefix)]
    results = {}
    batch_cols = [col for col in samplem.columns if col.startswith('batch')]
    batch_arg = samplem[batch_cols] if batch_cols else None
    if batch_arg is not None:
        print("Batch information detected, will be used in association.")
    else:
        print("No batch information detected, running without batch correction.")

    for cond in condition_cols:
        p = cna.tl.association(
            adata,
            samplem[cond],
            sample_col,
            covs=None,
            batches=batch_arg,
            allow_low_sample_size=True,
            key_added=cond + '_coef'
        )
        results[cond] = p
        print(f"{cond} global association p-value: {p}")
    return results

def summarize_condition_coefs(adata, celltype_col='celltype', coef_suffix='_coef', fdr_suffix='_coef_fdr', agg_func_coef='mean', agg_func_fdr='mean'):
    coef_cols = [col for col in adata.obs.columns if col.endswith(coef_suffix)]
    fdr_cols = [col for col in adata.obs.columns if col.endswith(fdr_suffix)]
    rows = []
    for celltype in adata.obs[celltype_col].unique():
        sub_df = adata.obs[adata.obs[celltype_col] == celltype]
        for coef_col, fdr_col in zip(coef_cols, fdr_cols):
            condition = coef_col.replace(coef_suffix, '')
            coef_vals = sub_df[coef_col].values
            fdr_vals = sub_df[fdr_col].values
            mean_coef = getattr(pd.Series(coef_vals), agg_func_coef)() if isinstance(agg_func_coef, str) else agg_func_coef(coef_vals)
            mean_fdr = getattr(pd.Series(fdr_vals), agg_func_fdr)() if isinstance(agg_func_fdr, str) else agg_func_fdr(fdr_vals)
            rows.append({'celltype': celltype, 'condition': condition, 'mean_coef': mean_coef, 'mean_fdr': mean_fdr})
    summary_df = pd.DataFrame(rows)
    return summary_df

def save_results_to_h5mu(adata, summary_df, output_path):
    summary_adata = ad.AnnData()
    summary_adata.uns['summary_table'] = summary_df
    mdata = mu.MuData({'main': adata, 'summary': summary_adata})
    mdata.write_h5mu(output_path)
    print(f"Save the result to {output_path}")

# ====== Dataset processing function ======
def process_single_dataset(data_file: Path):
    dataset_name = get_dataset_name(data_file)
    print(f"Start processing the dataset: {dataset_name}")

    output_dir = output_root / dataset_name
    output_dir.mkdir(parents=True, exist_ok=True)

    log_dir = log_root
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / f"{dataset_name}_cna_python.log"

    figures_dir = figures_root / dataset_name
    figures_dir.mkdir(parents=True, exist_ok=True)

    # Read the.h5ad dataset and exclude virtual files
    if data_file.suffix == ".h5ad" and not data_file.stem.startswith('._'):
        adata = ad.read_h5ad(data_file)
    else:
        print(f"Skip the files that do not need to be processed: {data_file}")
        return

    logger = setup_logger(log_path)
    log_system_info(logger)
    log_loaded_packages(logger)
    
    process = psutil.Process(os.getpid())
    start_time = time.time()
    start_mem = process.memory_info().rss / 1024 ** 2  # MB
    logger.info(f"Start running cna and the initial memory: {start_mem:.2f} MB")

    try:
        adata = add_one_hot_encoding(adata)
        samplem = prepare_samplem_from_adata(adata)
        batch_cols = [col for col in samplem.columns if col.startswith('batch_')]
        batches = samplem[batch_cols] if batch_cols else None

        preprocess_neighbors(adata)
        run_all_condition_associations(adata, samplem)
        summary_df = summarize_condition_coefs(adata)
        save_results_to_h5mu(adata, summary_df, output_dir / f"{dataset_name}_cna_python.h5mu")

        end_time = time.time()
        end_mem = process.memory_info().rss / 1024 ** 2
        mem_used = end_mem - start_mem

        summary_records.append({
            "dataset": dataset_name,
            "elapsed_time_sec": end_time - start_time,
            "memory_MB": mem_used,
            "status": "success"
        })

        logger.info(f"Memory at the end of operation: {end_mem:.2f} MB，Consume memory: {mem_used:.2f} MB")
        logger.info(f"CNA Method running time: {end_time - start_time:.2f} 秒")
        logger.info(f"{dataset_name} Processing completed\n")

    except Exception as e:
        logger.error(f"CNA Failed to run: {e}")
        summary_records.append({
            "dataset": dataset_name,
            "elapsed_time_sec": None,
            "memory_MB": None,
            "status": "failed"
        })
        print(f"❌ Process the dataset {dataset_name} failed: {e}")

# ====== Main program entry ======
def main():
    # Traverse the data folder and read all.h5ad files
    for sub_dir in data_dir.glob("*/"):  
        for data_file in sub_dir.glob("*.h5ad"):  
            process_single_dataset(data_file)

    # Save the summary table
    summary_df = pd.DataFrame(summary_records)
    summary_path = output_root / "cna_python_summary.csv"
    summary_df.to_csv(summary_path, index=False)
    print(f"The summary table of all dataset operations has been saved to: {summary_path}")

if __name__ == "__main__":
    main()
