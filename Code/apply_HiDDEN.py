import anndata as ad
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import found
from found.adapters import Pipeline
from found import methods as m
RANDOM_STATE = 42
adata = ad.read_h5ad("Data/human_prefrontal_cortex_HBCC_Cohort_BBB_cell_types_Schizophrenia_CaseControl_preprocessed.h5ad")

# %% [markdown]P
# next, we run the standard HiDDEN pipeline to classify affected cells on:
# normal bone marrow, smoldering multiple myeloma, and multiple myeloma patients

# %%
found.set_seed(RANDOM_STATE)  # set a fixed seed for replicability
algo = Pipeline.from_proc_ad("X_pca", m.log_reg, m.kmeans_bin)
p_hat, labs = found.find(adata, "Schizophrenia", "No", algo, k=30, 
                         adata=adata, regopt_maxiter=300, regopt_solver="newton-cg")
