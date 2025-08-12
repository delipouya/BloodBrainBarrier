import anndata as ad
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import found
from found.adapters import Pipeline
from found.tune import NaiveMinScoreTuner, score_phatdiff, score_deg
from found import methods as m
import seaborn as sns
import numpy as np
from matplotlib_venn import venn2
import pickle

RANDOM_STATE = 42
found.set_seed(RANDOM_STATE)  # set a fixed seed for replicability

#adata = ad.read_h5ad("Data/human_prefrontal_cortex_HBCC_Cohort_BBB_cell_types_Schizophrenia_CaseControl_preprocessed.h5ad")
#adata = ad.read_h5ad("/home/delaram/BloodBrainBarrier/Data/human_prefrontal_cortex_HBCC_Cohort_BBB_cell_types_Schizophrenia_extended.h5ad")
adata = ad.read_h5ad("/home/delaram/BloodBrainBarrier/Data/human_prefrontal_cortex_HBCC_Cohort_BBB_cell_types_Schizophrenia_extended_normalized.h5ad")

###

algo = Pipeline.from_proc_ad("X_pca", m.log_reg, m.kmeans_bin) 
#algo = Pipeline(m.run_lognorm_pca, m.log_reg, m.kmeans_bin, True)

#algo = Pipeline.from_proc_ad("X_pca", m.log_reg, m.kmeans_bin) ## run this if already processed
#p_hat, labs = found.find(adata, cond_col="Schizophrenia", control_val="No",
#                         algo=algo, k=30,  adata=adata, regopt_maxiter=300,
#                           regopt_solver="newton-cg")

optimal_k, res_dict = found.findt(
    adata,
    cond_col="Schizophrenia",
    control_val="No",
    tuner=NaiveMinScoreTuner(score_phatdiff, range(5, 31)),
    X=adata.X,  
    adata=adata,
    regopt_maxiter=300,
    regopt_solver="newton-cg"
)

for test_k in res_dict.keys():
    adata.obs['phat'] = res_dict[test_k][0]  # p_hat
    adata.obs['labs'] = res_dict[test_k][1]  # labs
    df = pd.DataFrame(adata.obs)

    data_for_plot = [df['phat'][df['Schizophrenia'] == cat] for cat in df['Schizophrenia'].unique()]
    categories = df['Schizophrenia'].unique()
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.violinplot(data_for_plot, showmeans=True, showmedians=True)
    ax.set_xticks(np.arange(1, len(categories) + 1))
    ax.set_xticklabels(categories, fontsize=16)
    ## set y ticks font size
    ax.tick_params(axis='y', labelsize=16)
    ax.set_title(f'Violin Plot of p-hat by Schizophrenia=Yes/No (K={test_k})', fontsize=16)
    ax.set_xlabel('Schizophrenia Status', fontsize=16)
    ax.set_ylabel('p-hat', fontsize=16)
    plt.ylim(0, 1) 
    ax.grid(True, linestyle='--', alpha=0.7)
    plt.show()

k = 29
adata.obs['phat'] = res_dict[k][0]  # p_hat
adata.obs['labs'] = res_dict[k][1]  # labs
df = pd.DataFrame(adata.obs)
print(adata.obs.cell_type.value_counts())
### violin plot of p-hat over cell types
plt.figure(figsize=(10, 6))
sns.violinplot(x='cell_type', y='phat', hue='cell_type', data=df)
plt.title(f'Grouped Violin Plot of p-hat by Cell Type - K={k}', fontsize=20)
plt.xlabel('Cell Type', fontsize=20)
plt.ylabel('p-hat', fontsize=20)
plt.xticks(rotation=45, ha='right', fontsize=20)
plt.yticks(fontsize=18)
plt.ylim(0, 1)
plt.grid(True, linestyle='--', alpha=0.7)
plt.show()

adata.obs['phat'] = res_dict[optimal_k][0]  # p_hat
adata.obs['labs'] = res_dict[optimal_k][1]  # labs


###########################################################################
sc.pl.umap(adata, color=['phat'], frameon=False, size=20, title='P-hat UMAP', wspace=0.4, hspace=0.4)
sc.pl.umap(adata, color=['Schizophrenia'], frameon=False, size=20, 
           title='Schizophrenia (original labels) UMAP', wspace=0.4, hspace=0.4)
sc.pl.umap(adata, color=['labs'], frameon=False, size=20, title='Corrected Labels UMAP', wspace=0.4, hspace=0.4)
###########################################################################
