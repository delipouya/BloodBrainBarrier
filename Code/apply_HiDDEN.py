import anndata as ad
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import found
from found.adapters import Pipeline
from found.tune import NaiveMinScoreTuner, score_phatdiff, score_deg, score_deg2
from found import methods as m
import seaborn as sns
import numpy as np
from matplotlib_venn import venn2
import pickle
import math
from scipy.sparse import isspmatrix_csr, isspmatrix_csc
from scipy.sparse import csr_array, csc_array

RANDOM_STATE = 42
found.set_seed(RANDOM_STATE)  # set a fixed seed for replicability

##############################
### add UMAP coordinates form adata_full_merged to adata
adata_full_normed = ad.read_h5ad("/home/delaram/BloodBrainBarrier/Data/human_prefrontal_cortex_HBCC_Cohort_BBB_cell_types_Schizophrenia_extended_normalized.h5ad")
UMAP_coordinates = pd.DataFrame(adata_full_normed.obsm['X_umap'])
UMAP_coordinates.index = adata_full_normed.obs.index
##############################

#adata = ad.read_h5ad("Data/human_prefrontal_cortex_HBCC_Cohort_BBB_cell_types_Schizophrenia_CaseControl_preprocessed.h5ad")
adata = ad.read_h5ad("/home/delaram/BloodBrainBarrier/Data/human_prefrontal_cortex_HBCC_Cohort_BBB_cell_types_Schizophrenia_extended.h5ad")

algo = Pipeline(m.run_lognorm_pca, m.log_reg, m.kmeans_bin, True)
#algo = Pipeline.from_proc_ad("X_pca", m.log_reg, m.kmeans_bin) ## run this if already processed

vis_k_choice = True  # set to True if you want to visualize the k choice
adata_sub_dict = {}
de_optimizer = False

optimal_K_dict = {}
for cell_type in adata.obs.cell_type.unique():
    print(f"Cell type: {cell_type}")
    #### subsetting data to include endothelial cells only
    adata_cell_type = adata[adata.obs.cell_type==cell_type,]
    X_in = adata_cell_type.X
    if de_optimizer:
        if isspmatrix_csr(X_in):
            X_in = csr_array(X_in)
        elif isspmatrix_csc(X_in):
            X_in = csc_array(X_in)

    #p_hat, labs = found.HiDDENt(adata_cell_type, cond_col="Schizophrenia", 
    #                         control_val="No", algo=algo, k=30,  
    #                         adata=adata_cell_type, regopt_maxiter=300, 
    #                         regopt_solver="newton-cg")
   
    optimal_k, res_dict = found.HiDDENt(
        adata_cell_type,
        cond_col="Schizophrenia",
        control_val="No",
        #tuner=NaiveMinScoreTuner(score_phatdiff, [2, 5, 10, 20, 25, 30]), #range(5, 31)
        tuner=NaiveMinScoreTuner(score_phatdiff, k_range=range(2, 31)),#[2, 5, 10, 20, 25, 30]
        X=X_in, #adata_cell_type.X,  
        algo=algo,
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
    adata_cell_type.obs['opt_k'] = optimal_k
    optimal_K_dict[cell_type] = optimal_k

    ### add phat and labs for various k as phat_{k value} and labs_{k_value}
    for test_k in res_dict.keys():
        adata_cell_type.obs[f'phat_{test_k}'] = res_dict[test_k][0]
        adata_cell_type.obs[f'labs_{test_k}'] = res_dict[test_k][1]
    
    adata_sub_dict[cell_type] = adata_cell_type
    print(f"Processed {cell_type} with optimal K={optimal_k}")
    print('-------------------------------------------------------')



# Save the dictionary to a binary file
with open('/home/delaram/BloodBrainBarrier/Data/human_prefrontal_cortex_HBCC_Cohort_BBB_adata_sub_dict.pkl', 'wb') as f:
    pickle.dump(adata_sub_dict, f)

# Load the dictionary back from the file
with open('/home/delaram/BloodBrainBarrier/Data/human_prefrontal_cortex_HBCC_Cohort_BBB_adata_sub_dict.pkl', 'rb') as f:
    adata_sub_dict = pickle.load(f)

for cell_type in adata_sub_dict.keys():
    adata_cell_type = adata_sub_dict[cell_type]
    print(cell_type, ' ', adata_cell_type.shape[0])

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
    adata_cell_type.obs = obs

inconsistent_relabel_dict = {}

######################################################################
### for each adata_cell_type, evaluate how many labs_2, labs_5, labs_10, labs_20, labs_25, labs_30 ('Yes', 'No') lists are consistent with each other given the exact order
for cell_type in adata_sub_dict.keys():

    adata_cell_type = adata_sub_dict[cell_type]
    lab_columns = [f'labs_{k}' for k in [2, 5, 10, 20, 25, 30]]
    lab_values = pd.DataFrame(adata_cell_type.obs[lab_columns].values)
    lab_values.columns = lab_columns

    ### add a column to the dataframe based on whether values in each line are the same
    K_vs_K30 = {}
    for k in [2, 5, 10, 20, 25, 30]:
        lab_values['notequal_{}_30'.format(k)] = lab_values['labs_{}'.format(k)] != lab_values['labs_30']
        adata_cell_type.obs[f'notequal_{k}_30'] = lab_values['notequal_{}_30'.format(k)].values

        print(adata_cell_type.obs['status'][adata_cell_type.obs[f'notequal_{k}_30']].value_counts())
        print(adata_cell_type.obs['sub_label'][adata_cell_type.obs[f'notequal_{k}_30']].value_counts())
        print(f"Cell type: {cell_type}, Not Equal labels (labs_{k}, labs_30): {lab_values['notequal_{}_30'.format(k)].sum()/ lab_values.shape[0] * 100:.2f}%")
        
        K_vs_K30[f"(labs_{k}, labs_30)"] = round(lab_values['notequal_{}_30'.format(k)].sum()/ lab_values.shape[0] * 100, 3)
    print('--------------------------')
    inconsistent_relabel_dict[cell_type] = K_vs_K30

### visualize the stats on inconsistency
for cell_type in inconsistent_relabel_dict.keys():
    print(cell_type, ' ', inconsistent_relabel_dict[cell_type])
    plt.figure(figsize=(10, 6))
    sns.barplot(x=list(inconsistent_relabel_dict[cell_type].keys()), y=list(inconsistent_relabel_dict[cell_type].values()))
    plt.title(f'Inconsistent Labeling for {cell_type}', fontsize=20)
    plt.xlabel('Label Comparison', fontsize=18)
    plt.ylabel('Percentage (%)', fontsize=18)
    plt.xticks(rotation=45, fontsize=16)
    plt.yticks(fontsize=16)
    plt.ylim(0, 100)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.show()

######################################################################

for adata_cell_type in adata_sub_dict.values():
    obs = adata_cell_type.obs 
    for k in [2, 5, 10, 20, 25, 30]:
        obs['sub_label_{}'.format(k)] = np.select(
            [
                obs['status'].eq('control'),
                obs['status'].eq('disease') & obs['labs_{}'.format(k)].eq('Yes'),
                obs['status'].eq('disease') & obs['labs_{}'.format(k)].eq('No'),
            ],
        ['control_0', 'case', 'control_1'],
        default='unknown'
    )
    adata_cell_type.obs = obs

for cell_type in adata_sub_dict.keys():
    adata_cell_type = adata_sub_dict[cell_type]
    lab_columns = [f'sub_label_{k}' for k in [2, 5, 10, 20, 25, 30]] #
    lab_values = pd.DataFrame(adata_cell_type.obs[lab_columns].values)
    lab_values.columns = lab_columns
    ### add a column to the dataframe based on whether values in each line are the same
    K_vs_K30 = {}
    for k in [2, 5, 10, 20, 25, 30]:
        lab_values['sub_label_notequal_{}_30'.format(k)] = lab_values['sub_label_{}'.format(k)] != lab_values['sub_label_30']
        adata_cell_type.obs[f'sub_label_notequal_{k}_30'] = lab_values['sub_label_notequal_{}_30'.format(k)].values
        sub_label_notequal_control_1 = adata_cell_type.obs[f'sub_label_notequal_{k}_30'] & (adata_cell_type.obs['sub_label_30'] == 'control_1')
        print(f"Cell type: {cell_type}, Not Equal labels (labs_{k}, labs_30): {sum(sub_label_notequal_control_1)/ sum(lab_values['sub_label_30']=='control_1') * 100:.2f}%")
        
        K_vs_K30[f"(labs_{k}, labs_30)"] = round(sum(sub_label_notequal_control_1)/ sum(lab_values['sub_label_30']=='control_1') * 100, 3)
    print('--------------------------')
    inconsistent_relabel_dict[cell_type] = K_vs_K30

for cell_type in inconsistent_relabel_dict.keys():
    print(cell_type, ' ', inconsistent_relabel_dict[cell_type])
    plt.figure(figsize=(10, 6))
    sns.barplot(x=list(inconsistent_relabel_dict[cell_type].keys()), y=list(inconsistent_relabel_dict[cell_type].values()))
    plt.title(f'Inconsistent Labeling for {cell_type}', fontsize=20)
    plt.xlabel('Label Comparison', fontsize=18)
    plt.ylabel('Percentage (%)', fontsize=18)
    plt.xticks(rotation=45, fontsize=16)
    plt.yticks(fontsize=16)
    plt.ylim(0, 100)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.show()

### the total number of relabeling in each K seems to the same within each cell type
for cell_type in inconsistent_relabel_dict.keys():
    adata_cell_type = adata_sub_dict[cell_type]
    print(cell_type)
    for k in [2, 5, 10, 20, 25, 30]:
        print(k, ' : ', sum(adata_cell_type.obs['sub_label_30']=='control_1'))
    print('--------------------------')



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
    print(adata_sub_dict[keys[i]].shape[0])

for i in range(len(keys)):
    print('cell type: ', keys[i])
    df = adata_sub_dict[keys[i]].obs
    k_name = optimal_K_dict[keys[i]]
    plt.figure(figsize=(12, 10))
    sns.violinplot(x='sub_label', y='phat', hue='sub_label', 
                   data=df, split=False, cut=0) #palette='viridis',
    plt.title(f'{keys[i]} p-hat (after label correction - opt K: {k_name})', fontsize=32)
    plt.ylabel('p-hat', fontsize=28)
    plt.xticks(rotation=45, ha='right', fontsize=28)
    plt.yticks(fontsize=28)
    plt.ylim(0, 1)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.show()

    '''
    for k in [2, 5, 10, 20, 25, 30]:
        print('K =', k)
        plt.figure(figsize=(12, 10))
        sns.violinplot(x='sub_label', y='phat_'+str(k), cut=0,
                       hue='sub_label', data=df, split=False) #palette='viridis',
        plt.title(f'{keys[i]} p-hat (after label correction) - K: {k}', fontsize=32)
        plt.ylabel('p-hat', fontsize=28)
        plt.xticks(rotation=45, ha='right', fontsize=28)
        plt.yticks(fontsize=28)
        plt.ylim(0, 1)
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.show()
    '''

################# Donor based plotting

for i in range(len(keys)):
    print('cell type: ', keys[i])
    df = adata_sub_dict[keys[i]].obs
    df_normal = df[df['sub_label'] == 'control_0']
    # donors per plot
    batch_size = 50  # increase as you like
    donors_all = sorted(df_normal['donor_id'].unique())
    n_batches = math.ceil(len(donors_all) / batch_size)

    for b in range(n_batches):
        donor_batch = donors_all[b*batch_size:(b+1)*batch_size]
        df_batch = df_normal[df_normal['donor_id'].isin(donor_batch)].copy()
        # Make donor_id a categorical with ONLY this batch's donors
        df_batch['donor_id'] = pd.Categorical(df_batch['donor_id'],
                                            categories=donor_batch,
                                            ordered=True)
        # (sometimes pandas keeps hidden categories; drop any)
        df_batch['donor_id'] = df_batch['donor_id'].cat.remove_unused_categories()
        fig, ax = plt.subplots(figsize=(16, 8))
        sns.violinplot(data=df_batch, x='donor_id', y='phat',
                    order=donor_batch, cut=0, scale='width', ax=ax)

        ax.set_title(f"{keys[i]} p-hat (after label correction) - batch {b+1}", fontsize=25)
        ax.set_ylabel('p-hat', fontsize=24)
        ax.set_xlabel('donor_id', fontsize=24)
        ax.set_ylim(0, 1)
        ### set fontsize of y ticks
        ax.tick_params(axis='y', labelsize=23)
        ax.set_xticklabels(donor_batch, rotation=45, ha='right', fontsize=17)
        # y-grid only (prevents the forest of vertical lines)
        ax.yaxis.grid(True, linestyle='--', alpha=0.6)
        ax.xaxis.grid(False)

        plt.tight_layout()
        plt.show()
    print('------------------------------')


for i in range(len(keys)):
    print('cell type: ', keys[i])
    df = adata_sub_dict[keys[i]].obs
    df_normal = df[df['sub_label'].isin(['control_1', 'case'])]
    # donors per plot
    batch_size = 30  # increase as you like
    donors_all = sorted(df_normal['donor_id'].unique())
    n_batches = math.ceil(len(donors_all) / batch_size)

    for b in range(n_batches):
        donor_batch = donors_all[b*batch_size:(b+1)*batch_size]
        df_batch = df_normal[df_normal['donor_id'].isin(donor_batch)].copy()
        # Make donor_id a categorical with ONLY this batch's donors
        df_batch['donor_id'] = pd.Categorical(df_batch['donor_id'],
                                            categories=donor_batch,
                                            ordered=True)
        # (sometimes pandas keeps hidden categories; drop any)
        df_batch['donor_id'] = df_batch['donor_id'].cat.remove_unused_categories()
        fig, ax = plt.subplots(figsize=(16, 8))
        sns.violinplot(data=df_batch, x='donor_id', y='phat', hue='sub_label',
                    order=donor_batch, cut=0, scale='width', ax=ax)

        ax.set_title(f"{keys[i]} p-hat (after label correction) - batch {b+1}", fontsize=25)
        ax.set_ylabel('p-hat', fontsize=24)
        ax.set_xlabel('donor_id', fontsize=24)
        ax.set_ylim(0, 1)
        ### set fontsize of y ticks
        ax.tick_params(axis='y', labelsize=23)
        ax.set_xticklabels(donor_batch, rotation=45, ha='right', fontsize=17)
        # y-grid only (prevents the forest of vertical lines)
        ax.yaxis.grid(True, linestyle='--', alpha=0.6)
        ax.xaxis.grid(False)

        plt.tight_layout()
        plt.show()
    print('------------------------------')


###########################################################################
### add astrocyte metadata to the adata object 
for i in range(len(keys)):
    df = adata_sub_dict[keys[i]].obs
    adata_tmp = adata.copy()
    df_full = adata_tmp.obs
    ### merge df with the df_full based on barcodekey - put the rest of the column as NA
    df_merged = df_full.merge(df[['phat', 'labs', 'sub_label']], left_index=True, right_index=True, how='left')
    adata_tmp.obs['phat_celltype'] = df_merged['phat']
    adata_tmp.obs['labs_celltype'] = df_merged['labs']
    adata_tmp.obs['sub_label_celltype'] = df_merged['sub_label']

    ### subset adata_tmp to only include cells in adata_full_normed
    adata_tmp = adata_tmp[adata_tmp.obs.index.isin(adata_full_normed.obs.index), :]
    adata_tmp.obsm['X_umap'] = UMAP_coordinates.loc[adata_tmp.obs.index].values
    ### plot umap using X_umap_2
    sc.pl.umap(adata_tmp, color=['phat_celltype'], frameon=False, size=20, 
               title=f'{keys[i]} P-hat UMAP', wspace=0.4, hspace=0.4)
    sc.pl.umap(adata_tmp, color=['sub_label_celltype'], frameon=False, 
               size=20, title=f'{keys[i]} Corrected Labels UMAP', wspace=0.4, hspace=0.4)



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
