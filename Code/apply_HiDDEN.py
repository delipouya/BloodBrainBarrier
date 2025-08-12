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

algo = Pipeline.from_proc_ad("X_pca", m.log_reg, m.kmeans_bin) 
#algo = Pipeline(m.run_lognorm_pca, m.log_reg, m.kmeans_bin, True)
#algo = Pipeline.from_proc_ad("X_pca", m.log_reg, m.kmeans_bin) ## run this if already processed

vis_k_choice = False  # set to True if you want to visualize the k choice
adata_sub_dict = {}

for cell_type in adata.obs.cell_type.unique():
    print(f"Cell type: {cell_type}")
   
    #### subsetting data to include endothelial cells only
    adata_cell_type = adata[adata.obs.cell_type==cell_type,]
    #p_hat, labs = found.find(adata_cell_type, cond_col="Schizophrenia", 
    #                         control_val="No", algo=algo, k=30,  
    #                         adata=adata_cell_type, regopt_maxiter=300, 
    #                         regopt_solver="newton-cg")
    optimal_k, res_dict = found.findt(
        adata_cell_type,
        cond_col="Schizophrenia",
        control_val="No",
        tuner=NaiveMinScoreTuner(score_phatdiff, [2, 5, 10, 15, 20, 25, 30]), #range(5, 31)
        X=adata_cell_type.X,  
        adata=adata_cell_type,
        regopt_maxiter=300,
        regopt_solver="newton-cg"
    )

    print('optimal K is', optimal_k)

    if vis_k_choice:
        for test_k in res_dict.keys():
            adata_cell_type.obs['phat'] = res_dict[test_k][0]  # p_hat
            adata_cell_type.obs['labs'] = res_dict[test_k][1]  # labs            
            df_vis = pd.DataFrame(adata_cell_type.obs)
            data_for_plot = [df_vis['phat'][df_vis['Schizophrenia'] == cat] for cat in df_vis['Schizophrenia'].unique()]
            categories = df_vis['Schizophrenia'].unique()
            fig, ax = plt.subplots(figsize=(8, 6))
            ax.violinplot(data_for_plot, showmeans=True, showmedians=True)
            ax.set_xticks(np.arange(1, len(categories) + 1))
            ax.set_xticklabels(categories, fontsize=16)
            ax.set_title(f'Violin Plot of p-hat by Schizophrenia=Yes/No (K={test_k})', fontsize=16)
            ax.set_xlabel('Schizophrenia Status', fontsize=16)
            ax.set_ylabel('p-hat', fontsize=16)
            plt.ylim(0, 1) 
            ax.grid(True, linestyle='--', alpha=0.7)
            plt.show()

    ### adding the optimal-K based results to the object
    adata_cell_type.obs['phat'] = res_dict[optimal_k][0]
    adata_cell_type.obs['labs'] = res_dict[optimal_k][1]
    adata_sub_dict[cell_type] = adata_cell_type


# Save the dictionary to a binary file
with open('/home/delaram/BloodBrainBarrier/Data/human_prefrontal_cortex_HBCC_Cohort_BBB_adata_sub_dict.pkl', 'wb') as f:
    pickle.dump(adata_sub_dict, f)

# Load the dictionary back from the file
with open('/home/delaram/BloodBrainBarrier/Data/human_prefrontal_cortex_HBCC_Cohort_BBB_adata_sub_dict.pkl', 'rb') as f:
    adata_sub_dict = pickle.load(f)


######################################################################
df = pd.DataFrame(adata.obs)
for adata_cell_type in adata_sub_dict.values():
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

    ### compare labs from full dataset and cell type specific dataset
    cell_type_labs = adata_cell_type.obs['labs'].tolist()
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



for adata_cell_type in adata_sub_dict.values():
    ### if a cell status was as control, keep as "control_0", 
    # if status was "disease" before, and updated label is 'Yes' (labs), label as "case"
    # if status was "disease" before, but updated label is 'No' (labs), label as "control_1"
    
    obs = adata_cell_type.obs 
    obs['sub_label'] = np.select(
        [
            obs['status'].eq('control'),
            obs['status'].eq('disease') & obs['labs'].eq('Yes'),
            obs['status'].eq('disease') & obs['labs'].eq('No'),
        ],
        ['control_0', 'case', 'control_1'],
        default='unknown'
    )

    # make it an ordered categorical for consistent plotting
    obs['sub_label'] = pd.Categorical(obs['sub_label'],
                                    categories=['control_0', 'control_1', 'case'],
                                    ordered=True)
    adata_cell_type.obs = obs

######################################################################
############### Cell type specific p-hat evaluation ###############
######################################################################
keys = ['astrocyte', 'endothelial cell', 'pericyte']
keys = ['astrocyte', 'endothelial cell', 'microglial cell', 'perivascular macrophage', 
       'pericyte', 'smooth muscle cell', 'vascular leptomeningeal cell']

for i in range(len(keys)):
    print('cell type: ', keys[i])
    df = adata_sub_dict[keys[i]].obs

    ###### Violin plot for p-hat by Schizophrenia status ######
    plt.figure(figsize=(16, 12))
    sns.violinplot(x='Schizophrenia', y='phat', hue='Schizophrenia', data=df, split=False) #palette='viridis',
    plt.title(f'{keys[i]} p-hat (Before label correction)', fontsize=32)
    plt.xlabel('Schizophrenia Status', fontsize=28)
    plt.ylabel('p-hat', fontsize=28)
    plt.xticks(rotation=45, ha='right', fontsize=28)
    plt.yticks(fontsize=28)
    plt.ylim(0, 1)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.show()

for i in range(len(keys)):
    print('cell type: ', keys[i])
    df = adata_sub_dict[keys[i]].obs
    plt.figure(figsize=(12, 10))
    sns.violinplot(x='sub_label', y='phat', hue='sub_label', data=df, split=False) #palette='viridis',
    plt.title(f'{keys[i]} p-hat (after label correction)', fontsize=32)
    plt.ylabel('p-hat', fontsize=28)
    plt.xticks(rotation=45, ha='right', fontsize=28)
    plt.yticks(fontsize=28)
    plt.ylim(0, 1)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.show()


###########################################################################
### add astrocyte metadata to the adata object 
df = adata_sub_dict[keys[i]].obs

adata_tmp = adata.copy()
df_full = adata_tmp.obs
### merge df with the df_full based on barcodekey - put the rest of the column as NA
df_merged = df_full.merge(df[['phat', 'labs']], left_index=True, right_index=True, how='left')
adata_tmp.obs['phat_celltype'] = df_merged['phat_y']
adata_tmp.obs['labs_celltype'] = df_merged['labs_y']
sc.pl.umap(adata_tmp, color=['phat_celltype'], frameon=False, size=20, title='P-hat UMAP', wspace=0.4, hspace=0.4)
sc.pl.umap(adata_tmp, color=['labs_celltype'], frameon=False, size=20, title='Corrected Labels UMAP', wspace=0.4, hspace=0.4)


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


plt.figure(figsize=(27, 12))
sns.violinplot(x='disease', y='phat', hue='labs', data=df, split=False) #palette='viridis',
plt.title('Violin Plot of p-hat by Disease Status', fontsize=32)
plt.xlabel('Disease Status', fontsize=28)
plt.ylabel('p-hat', fontsize=28)
plt.xticks(rotation=45, ha='right', fontsize=31)
plt.yticks(fontsize=28)
plt.legend(title='Label', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True, linestyle='--', alpha=0.7)
plt.show()


df2 = df[df['status'].isin(['disease'])]
plt.figure(figsize=(180, 70))
sns.violinplot(x='donor_id', y='phat', hue='disease', data=df2, split=False) #palette='viridis',
plt.title('Violin Plot of p-hat by Donor ID', fontsize=32)
plt.xlabel('Donor ID', fontsize=50)
plt.ylabel('p-hat', fontsize=50)
plt.xticks(rotation=45, ha='right', fontsize=40)
plt.yticks(fontsize=40)
plt.legend(title='Label', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True, linestyle='--', alpha=0.7)
plt.show()
