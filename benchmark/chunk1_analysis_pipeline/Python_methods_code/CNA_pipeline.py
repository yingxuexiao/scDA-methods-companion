#CNA_pipeline:
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
import gc

# ================== CPU Limitations ==================
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
torch.set_num_threads(1)

# Bind the current process to CPU core 19
p = psutil.Process(os.getpid())
p.cpu_affinity([19])
print(f"CPU core bound: {p.cpu_affinity()}")

# ====== Basic Path Configuration ======
base_dir = Path("/Data/test/xiaoyingxue").resolve()
data_dir = base_dir / "simulation_data"
output_root = base_dir / "da_results"
log_root = base_dir / "log"
figures_root = base_dir / "figures"

# ====== Global Summary List ======
summary_records = []

# ================== Memory Management Functions ==================
def cleanup_memory():
    """Clean up memory thoroughly, ensuring each dataset starts from a clean state"""
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
    """Get stable memory usage measurements"""
    measurements = []
    for _ in range(samples):
        measurements.append(process.memory_info().rss / 1024 ** 2)
        time.sleep(interval)
    # Use median to reduce the impact of outliers
    return np.median(measurements)

# ====== Utility Functions ======
def get_dataset_name(data_file_path: Path) -> str:
    if data_file_path.suffix == ".h5ad":
        return data_file_path.stem
    else:
        raise ValueError("Only .h5ad files are supported as input datasets")

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
    # Handle condition columns
    condition_cols = [col for col in obs.columns if col.startswith(condition_prefix)]
    if not condition_cols:
        if 'condition' in obs.columns:
            dummies = pd.get_dummies(obs['condition'], prefix=condition_prefix.rstrip('_'))
            obs = pd.concat([obs, dummies], axis=1)
            condition_cols = dummies.columns.tolist()
        else:
            raise ValueError("No 'condition' column found and no encoded 'condition_' columns")
    # Handle batch columns
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
    print(f"Results saved to {output_path}")

# ====== Dataset Processing Functions ======
def process_single_dataset(data_file: Path):
    dataset_name = get_dataset_name(data_file)
    print(f"Processing dataset: {dataset_name}")
    print("-" * 50)

    output_dir = output_root / dataset_name
    output_dir.mkdir(parents=True, exist_ok=True)

    log_dir = log_root
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / f"{dataset_name}_cna_python.log"

    figures_dir = figures_root / dataset_name
    figures_dir.mkdir(parents=True, exist_ok=True)

    # Thoroughly clean memory before processing
    cleanup_memory()
    
    # Get process object
    process = psutil.Process(os.getpid())
    
    # Get stable starting memory
    start_mem = get_stable_memory_usage(process, samples=3, interval=0.3)
    
    # Record start time
    start_time = time.time()
    
    # Read .h5ad dataset, excluding virtual files
    if data_file.suffix == ".h5ad" and not data_file.stem.startswith('._'):
        adata = ad.read_h5ad(data_file)
        print(f"Data loaded, number of cells: {adata.n_obs}, number of genes: {adata.n_vars}")
    else:
        print(f"Skipping file: {data_file}")
        return

    # Logger initialization (placed after memory measurement to avoid interference)
    logger = setup_logger(log_path)
    log_system_info(logger)
    log_loaded_packages(logger)
    logger.info(f"Dataset: {dataset_name}")
    logger.info(f"Initial memory usage: {start_mem:.2f} MB")
    logger.info(f"Cells: {adata.n_obs}, Genes: {adata.n_vars}")

    try:
        logger.info("Starting CNA analysis...")

        # 1. Add one-hot encoding
        print("Step 1: Adding one-hot encoding...")
        logger.info("Step 1: Adding one-hot encoding")
        adata = add_one_hot_encoding(adata)
        
        # 2. Prepare sample matrix
        print("Step 2: Preparing sample matrix...")
        logger.info("Step 2: Preparing sample matrix")
        samplem = prepare_samplem_from_adata(adata)
        batch_cols = [col for col in samplem.columns if col.startswith('batch_')]
        batches = samplem[batch_cols] if batch_cols else None
        logger.info(f"Sample matrix shape: {samplem.shape}")
        logger.info(f"Condition columns: {[col for col in samplem.columns if col.startswith('condition_')]}")
        if batch_cols:
            logger.info(f"Batch columns: {batch_cols}")
        
        # 3. Preprocess neighbors
        print("Step 3: Calculating neighbors...")
        logger.info("Step 3: Calculating neighbors")
        preprocess_neighbors(adata)
        
        # 4. Run condition association analysis
        print("Step 4: Running condition association analysis...")
        logger.info("Step 4: Running condition association analysis")
        run_all_condition_associations(adata, samplem)
        
        # 5. Summarize results
        print("Step 5: Summarizing results...")
        logger.info("Step 5: Summarizing results")
        summary_df = summarize_condition_coefs(adata)
        logger.info(f"Summary table shape: {summary_df.shape}")
        
        # 6. Save results
        print("Step 6: Saving results...")
        logger.info("Step 6: Saving results")
        save_results_to_h5mu(adata, summary_df, output_dir / f"{dataset_name}_cna_python.h5mu")
        
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
        logger.info(f"CNA method runtime: {runtime:.2f} seconds")
        logger.info(f"Ending memory usage: {end_mem:.2f} MB")
        logger.info(f"Memory consumed: {mem_used:.2f} MB")
        logger.info(f"{dataset_name} processing completed")

        summary_records.append({
            "dataset": dataset_name,
            "elapsed_time_sec": round(runtime, 4),
            "memory_MB": round(mem_used, 4),
            "status": "success"
        })

        print(f"✓ Completed - Time: {runtime:.2f} seconds, Memory: {mem_used:.2f} MB")

    except Exception as e:
        logger.error(f"CNA analysis failed: {e}")
        summary_records.append({
            "dataset": dataset_name,
            "elapsed_time_sec": None,
            "memory_MB": None,
            "status": "failed"
        })
        print(f"❌ Failed to process dataset {dataset_name}: {e}")
    
    # Clean up objects and memory for the next dataset
    try:
        del adata
    except:
        pass
    
    cleanup_memory()

# ====== Main Program Entry ======
def main():
    # Clean up memory at the start of the program
    cleanup_memory()
    
    print("=" * 60)
    print("Optimized CNA Method Version")
    print("=" * 60)
    
    # Collect all h5ad files
    h5ad_files = []
    for sub_dir in data_dir.glob("*/"):  # Process all subdirectories
        for data_file in sub_dir.glob("*.h5ad"):  # Only process .h5ad files
            h5ad_files.append(data_file)
    
    print(f"Found {len(h5ad_files)} .h5ad files")
    
    if not h5ad_files:
        print("❌ No .h5ad files found, please check the data directory")
        return
    
    # Sort by file size (small to large)
    h5ad_files_info = []
    for file in h5ad_files:
        try:
            file_size_mb = os.path.getsize(file) / (1024 * 1024)
            h5ad_files_info.append((file, file_size_mb))
        except:
            h5ad_files_info.append((file, 0))
    
    # Sort by file size (small to large)
    h5ad_files_info.sort(key=lambda x: x[1])
    h5ad_files = [file for file, _ in h5ad_files_info]
    
    print("\nDataset processing order (from small to large):")
    for i, (file, size_mb) in enumerate(h5ad_files_info):
        print(f"  {i+1}. {file.stem}: {size_mb:.1f} MB")
    
    print("=" * 60)
    
    # Process each dataset
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
        
        # Add a pause between datasets
        print("-" * 50)
        time.sleep(1.0)

    # Save summary table
    if summary_records:
        summary_df = pd.DataFrame(summary_records)
        summary_path = output_root / "cna_python_summary.csv"
        summary_df.to_csv(summary_path, index=False)
        print(f"\n✅ All datasets processed, summary table saved to: {summary_path}")
        
        # Print summary statistics
        successful_runs = summary_df[summary_df['status'] == 'success']
        if len(successful_runs) > 0:
            print("\n" + "="*60)
            print("CNA Method Summary Statistics:")
            print("="*60)
            print(f"Number of successfully processed datasets: {len(successful_runs)}")
            print(f"Number of failed datasets: {len(summary_df) - len(successful_runs)}")
            print(f"Average runtime: {successful_runs['elapsed_time_sec'].mean():.2f} seconds")
            print(f"Average memory usage: {successful_runs['memory_MB'].mean():.2f} MB")
            print(f"Maximum memory usage: {successful_runs['memory_MB'].max():.2f} MB")
            print(f"Minimum memory usage: {successful_runs['memory_MB'].min():.2f} MB")
            print(f"Total runtime: {successful_runs['elapsed_time_sec'].sum():.2f} seconds")
            print("="*60)
            
            # Show memory usage ranking
            print("\nMemory usage ranking (from high to low):")
            sorted_results = successful_runs.sort_values('memory_MB', ascending=False)
            for i, (_, row) in enumerate(sorted_results.iterrows()):
                print(f"  {i+1}. {row['dataset']}: {row['memory_MB']:.2f} MB")
        else:
            print("\n⚠️ No datasets were processed successfully")
    else:
        print("\n⚠️ No datasets were processed")

if __name__ == "__main__":
    main()
