import anndata as ad
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import found
from found.adapters import Pipeline
from found.tune import NaiveMinScoreTuner, score_phatdiff, score_deg
from found import methods as m
RANDOM_STATE = 42
adata = ad.read_h5ad("Data/human_prefrontal_cortex_HBCC_Cohort_BBB_cell_types_Schizophrenia_CaseControl_preprocessed.h5ad")

# next, we run the standard HiDDEN pipeline to classify affected cells on:
# normal bone marrow, smoldering multiple myeloma, and multiple myeloma patients


found.set_seed(RANDOM_STATE)  # set a fixed seed for replicability
algo = Pipeline.from_proc_ad("X_pca", m.log_reg, m.kmeans_bin)
p_hat, labs = found.find(adata, "Schizophrenia", "No", algo, k=30, 
                         adata=adata, regopt_maxiter=300, regopt_solver="newton-cg")

p_hat, labs = found.findt(adata, "Schizophrenia", "No", algo, 
                          tuner = NaiveMinScoreTuner(score_phatdiff, range(10, 30)), 
                          adata=adata, regopt_maxiter=300, regopt_solver="newton-cg")



import anndata as ad
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import found
import seaborn as sns
from found.adapters import Pipeline
import numpy as np
from found import methods as m
from matplotlib_venn import venn2

RANDOM_STATE = 42

#### import processed data
adata = ad.read_h5ad("Data/human_prefrontal_cortex_HBCC_Cohort_BBB_cell_types_Schizophrenia_CaseControl_preprocessed.h5ad")

found.set_seed(RANDOM_STATE)  # set a fixed seed for replicability
algo = Pipeline.from_proc_ad("X_pca", m.log_reg, m.kmeans_bin)
### control="No"
p_hat, labs = found.find(adata, cond_col="Schizophrenia", control_val="No",
                         algo=algo, k=30,  adata=adata, regopt_maxiter=300,
                           regopt_solver="newton-cg")

adata.obs['phat'] = p_hat
adata.obs['labs'] = labs
df = pd.DataFrame(adata.obs)

print(adata.obs.cell_type.value_counts())

adata_sub_dict = {}
for cell_type in adata.obs.cell_type.unique():
    print(f"Cell type: {cell_type}")
   
    #### subsetting data to include endothelial cells only
    adata_cell_type = adata[adata.obs.cell_type==cell_type,]
    p_hat, labs = found.find(adata_cell_type, cond_col="Schizophrenia", 
                             control_val="No", algo=algo, k=30,  
                             adata=adata_cell_type, regopt_maxiter=300, 
                             regopt_solver="newton-cg")

    adata_cell_type.obs['phat'] = p_hat
    adata_cell_type.obs['labs'] = labs
    adata_sub_dict[cell_type] = adata_cell_type

    df_cell_type = pd.DataFrame(adata_cell_type.obs)
    #### compare the p-hat for endothelial cells with the full dataset using scatter plot   
    plt.figure(figsize=(10, 6))
    df_full = df[df['cell_type'] == cell_type]
    plt.scatter(df_full['phat'], df_cell_type['phat'], alpha=0.5, color='blue',
                edgecolor='k', s=50)
    plt.title(f'Comparison of p-hat for {cell_type} vs Full Dataset')
    plt.xlabel('p-hat (Full Dataset)')
    plt.ylabel(f'p-hat ({cell_type})')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.plot([0, 1], [0, 1], color='red', linestyle='--', linewidth=1)
    plt.show()

    ### create a venn diagram to compare labs from full dataset and cell type specific dataset
    cell_type_labs = labs.tolist()
    full_labs = df_full['labs'].tolist()

    label_df = pd.DataFrame({
        'cell_type_labs': cell_type_labs,
        'full_labs': full_labs
    })
    label_df['shared'] = label_df['cell_type_labs'].isin(label_df['full_labs'])
    label_df['unique_full'] = ~label_df['full_labs'].isin(label_df['cell_type_labs'])
    label_df['unique_cell_type'] = ~label_df['cell_type_labs'].isin(label_df['full_labs'])
    shared_count = label_df['shared'].sum()
    unique_full_count = label_df['unique_full'].sum()
    unique_cell_type_count = label_df['unique_cell_type'].sum()
    print(f"Shared labels: {shared_count}, Unique to full dataset: {unique_full_count}, Unique to {cell_type}: {unique_cell_type_count}")



keys = ['astrocyte', 'endothelial cell', 'pericyte']
i = 0
df = adata_sub_dict[keys[i]].obs
#df = pd.DataFrame(adata_sub_dict['astrocyte'].obs)
data_for_plot = [df['phat'][df['Schizophrenia'] == cat] for cat in df['Schizophrenia'].unique()]
categories = df['Schizophrenia'].unique()
fig, ax = plt.subplots(figsize=(8, 6))
ax.violinplot(data_for_plot, showmeans=True, showmedians=True)
ax.set_xticks(np.arange(1, len(categories) + 1))
ax.set_xticklabels(categories)
ax.set_title('Violin Plot of p-hat by Schizophrenia=Yes/No')
ax.set_xlabel('Schizophrenia Status')
ax.set_ylabel('p-hat')
ax.grid(True, linestyle='--', alpha=0.7)
plt.show()


data_for_plot = [df['phat'][df['Schizophrenia'] == cat] for cat in df['Schizophrenia'].unique()]
categories = df['Schizophrenia'].unique()
fig, ax = plt.subplots(figsize=(8, 6))
ax.violinplot(data_for_plot, showmeans=True, showmedians=True)
ax.set_xticks(np.arange(1, len(categories) + 1))
ax.set_xticklabels(categories)
ax.set_title('Violin Plot of p-hat by Schizophrenia=Yes/No')
ax.set_xlabel('Schizophrenia Status')
ax.set_ylabel('p-hat')
ax.grid(True, linestyle='--', alpha=0.7)
plt.show()

plt.figure(figsize=(10, 6))
sns.violinplot(x='cell_type', y='phat', hue='Schizophrenia', data=df, palette='viridis', split=True)
plt.title('Grouped Violin Plot of p-hat by Cell Type and Schizophrenia Status')
plt.xlabel('Cell Type')
plt.ylabel('p-hat')
plt.legend(title='Schizophrenia Status', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True, linestyle='--', alpha=0.7)
plt.show()


####################################################################################################
######################### p-hat distribution in various Schizophrenia samples #########################
plt.figure(figsize=(20, 12))
sns.violinplot(x='disease', y='phat', data=df, split=False) #palette='viridis',
plt.title('Violin Plot of p-hat by Disease Status', fontsize=32)
plt.xlabel('Disease Status', fontsize=28)
plt.ylabel('p-hat', fontsize=28)
plt.xticks(rotation=45, ha='right', fontsize=28)
plt.yticks(fontsize=26)
plt.legend(title='Disease', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True, linestyle='--', alpha=0.7)
plt.show()


plt.figure(figsize=(24, 12))
sns.violinplot(x='disease', y='phat', hue='labs', data=df, split=False) #palette='viridis',
plt.title('Violin Plot of p-hat by Disease Status', fontsize=32)
plt.xlabel('Disease Status', fontsize=28)
plt.ylabel('p-hat', fontsize=28)
plt.xticks(rotation=45, ha='right', fontsize=28)
plt.yticks(fontsize=28)
plt.legend(title='Disease', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True, linestyle='--', alpha=0.7)
plt.show()
###########################################################################

####################################################################################################
######################### cell type specific plots - corrected and uncorrected #########################
plt.figure(figsize=(24, 12))
sns.violinplot(x='cell_type', y='phat', hue='labs', data=df, split=False) #palette='viridis',
plt.title('Violin Plot of p-hat by Cell Type (After label correction)', fontsize=32)
plt.xlabel('Cell Type', fontsize=28)
plt.ylabel('p-hat', fontsize=28)
plt.xticks(rotation=45, ha='right', fontsize=28)
plt.yticks(fontsize=28)
plt.legend(title='Disease', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True, linestyle='--', alpha=0.7)
plt.show()


plt.figure(figsize=(24, 12))
sns.violinplot(x='cell_type', y='phat', hue='Schizophrenia', data=df, split=False) #palette='viridis',
plt.title('Violin Plot of p-hat by Cell Type (Before label correction)', fontsize=32)
plt.xlabel('Cell Type', fontsize=28)
plt.ylabel('p-hat', fontsize=28)
plt.xticks(rotation=45, ha='right', fontsize=28)
plt.yticks(fontsize=28)
plt.legend(title='Disease', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True, linestyle='--', alpha=0.7)
plt.show()
###########################################################################
sc.pl.umap(adata, color=['phat'], frameon=False, size=20, title='P-hat UMAP', wspace=0.4, hspace=0.4)

sc.pl.umap(adata, color=['Schizophrenia'], frameon=False, size=20, title='Schizophrenia (original labels) UMAP', wspace=0.4, hspace=0.4)
sc.pl.umap(adata, color=['labs'], frameon=False, size=20, title='Corrected Labels UMAP', wspace=0.4, hspace=0.4)
