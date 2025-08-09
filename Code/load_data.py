import anndata as ad
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc

adata = ad.read_h5ad("/home/delaram/BloodBrainBarrier/Data/human_prefrontal_cortex_HBCC_Cohort.h5ad")
print(adata.shape)
metadata = pd.DataFrame(adata.obs)
metadata.head()

### save/read metadata
#metadata.to_csv("Data/human_prefrontal_cortex_HBCC_Cohort_metadata.csv", index=False)
# Load metadata from CSV
#metadata = pd.read_csv("../Data/human_prefrontal_cortex_HBCC_Cohort_metadata.csv")
#metadata = pd.read_csv("../Data/human_prefrontal_cortex_HBCC_Cohort_Schizophrenia_metadata.csv")
metadata_immune_cells = metadata[metadata['cell_type'].isin(['B cell', 
                                                             'plasma cell', 
                                                             'T cell',
                                                              'natural killer cell'])]
print(sum(metadata_immune_cells.disease=='normal'))
### exclude disease=='normal
metadata_immune_cells = metadata_immune_cells[metadata_immune_cells['disease'] != 'normal']

plt.figure(figsize=(13, 10)) # (25, 18) (25, 15) for larger plots
metadata_immune_cells['disease'].value_counts().plot(kind='bar')
plt.title(f"Counts for immune cell types", fontsize=26)
plt.xlabel('disease by immune cell', fontsize=26)
plt.ylabel("Count", fontsize=26)
plt.xticks(rotation=45, ha='right', fontsize=26)
plt.yticks(fontsize=26)
plt.tight_layout()
plt.show()

#print(metadata[column_name].value_counts())
#adata = adata[adata.obs[column_name] == 'Yes']
#adata_subset.write("Data/human_prefrontal_cortex_HBCC_Cohort_Schizophrenia.h5ad")
#adata = ad.read_h5ad("Data/human_prefrontal_cortex_HBCC_Cohort_Schizophrenia.h5ad")

#'tissue_type'=tissue! al
columns_to_check = ['class', 'subclass',  'suspension_type',
    'genetic_ancestry', 'assay', 'sex', 'tissue']
columns_to_check = ['disease', 'cell_type']

# Make a bar plot for each column
for column_name in columns_to_check:
    if column_name in adata.obs.columns:
        plt.figure(figsize=(25, 18)) # (25, 15) for larger plots
        adata.obs[column_name].value_counts().plot(kind='bar')
        plt.title(f"Counts for {column_name}")
        plt.xlabel(column_name, fontsize=26)
        plt.ylabel("Count", fontsize=26)
        plt.xticks(rotation=45, ha='right', fontsize=26)
        plt.yticks(fontsize=26)
        plt.tight_layout()
        plt.show()

## count unique 'donor_id'
unique_donors = adata.obs['donor_id'].nunique()
print(f"Unique donor IDs: {unique_donors}")

### make sample-level metadata for donor_id and sex, genetic_ancestry, disease and then make bar plots for each
sample_metadata = adata.obs[['donor_id', 'sex', 'genetic_ancestry', 'disease']]
sample_metadata = sample_metadata.drop_duplicates()

for column_name in sample_metadata.columns:
    plt.figure(figsize=(25, 18)) #(10, 6)
    sample_metadata[column_name].value_counts().plot(kind='bar')
    plt.title(f"Counts for {column_name}")
    plt.xlabel(column_name, fontsize=26)
    plt.ylabel("Count", fontsize=26)
    plt.xticks(rotation=45, ha='right', fontsize=26)
    plt.yticks(fontsize=26)
    plt.tight_layout()
    plt.show()

##########################################################
BBB_cell_types = ['astrocyte', 'endothelial cell', 'pericyte']
BBB_cell_types = ['astrocyte', 'endothelial cell', 'pericyte', 
                  'microglial cell', 
                  'vascular leptomeningeal cell', 
                  'perivascular macrophage', 
                  'smooth muscle cell']
## subset to only BBB cell types
adata_BBB = adata[adata.obs['cell_type'].isin(BBB_cell_types)]
print(adata.shape)
print(adata_BBB.shape)

columns_to_check = ['class', 'subclass', 'genetic_ancestry', 'sex']

# Make a bar plot for each column
for column_name in columns_to_check:
    if column_name in adata_BBB.obs.columns:
        plt.figure(figsize=(12, 10)) # (25, 18) for larger plots
        adata_BBB.obs[column_name].value_counts().plot(kind='bar')
        #plt.title(f"Counts for {column_name}")
        plt.xlabel(column_name, fontsize=26)
        #plt.ylabel("Count", fontsize=26)
        plt.xticks(rotation=45, ha='right', fontsize=30)
        plt.yticks(fontsize=28)
        plt.tight_layout()
        plt.show()

adata_BBB.obs['status'] = adata_BBB.obs['disease'].apply(lambda x: 'disease' if x != 'normal' else 'control')
column_name = 'Schizophrenia'
print(adata_BBB.obs[column_name].value_counts())
condition = (adata_BBB.obs[column_name] == 'Yes') | (adata_BBB.obs['disease'] == 'normal')
adata_BBB_SZ_CR = adata_BBB[condition]
print(adata_BBB_SZ_CR.obs['disease'].value_counts())
### save AnnData file using ad
adata_BBB_SZ_CR.write("/home/delaram/BloodBrainBarrier/human_prefrontal_cortex_HBCC_Cohort_BBB_cell_types_Schizophrenia_extended.h5ad")

####### disease distribution #######
### exclude normal ones
adata_temp = adata_BBB_SZ_CR[adata_BBB_SZ_CR.obs['disease'] != 'normal']
plt.figure(figsize=(22, 15)) # (25, 18) for larger plots
adata_temp.obs['disease'].value_counts().plot(kind='bar')
#plt.title(f"Counts for {column_name}")
plt.xlabel('disease', fontsize=26)
#plt.ylabel("Count", fontsize=26)
plt.xticks(rotation=45, ha='right', fontsize=30)
plt.yticks(fontsize=28)
plt.tight_layout()
plt.show()

###################################
df = pd.DataFrame(adata_BBB_SZ_CR.obs)
# Count cells per cell_type per status
counts = df.groupby(["cell_type", "status"]).size().unstack(fill_value=0)
# Plot as grouped bars (side-by-side)
counts.plot(kind="bar", stacked=False)
plt.ylabel("Number of cells")
plt.xlabel("Cell type")
plt.xticks(rotation=45, ha='right')
plt.legend(title="Status")
plt.tight_layout()
plt.show()

unique_donors = adata_BBB_SZ_CR.obs['donor_id'].nunique()
print(f"Unique donor IDs: {unique_donors}")

### make sample-level metadata for donor_id and sex, genetic_ancestry, disease and then make bar plots for each
sample_metadata = adata_BBB_SZ_CR.obs[['donor_id', 'sex', 'genetic_ancestry', 'disease', 'status']]
sample_metadata = sample_metadata.drop_duplicates()


df = pd.DataFrame(sample_metadata)
counts = df.groupby(["sex", "status"]).size().unstack(fill_value=0)

plt.figure(figsize=(3, 4)) #(10, 6)
counts.plot(kind="bar", stacked=False)
plt.title('')
plt.xticks(rotation=45, ha='right', fontsize=15)
plt.yticks(fontsize=15)
plt.legend(title="Status")
plt.tight_layout()
plt.show()

counts = df.groupby(["genetic_ancestry", "status"]).size().unstack(fill_value=0)
plt.figure(figsize=(12, 6)) #(10, 6)
counts.plot(kind="bar", stacked=False)
plt.title('')
plt.xticks(rotation=45, ha='right', fontsize=15)
plt.yticks(fontsize=15)
plt.legend(title="Status")
plt.tight_layout()
plt.show()

column_name = 'status'
plt.figure(figsize=(9, 9)) #(10, 6)
sample_metadata[column_name].value_counts().plot(kind='bar')
plt.title(f"Counts for {column_name}")
plt.xlabel(column_name, fontsize=26)
plt.ylabel("Count", fontsize=26)
plt.xticks(rotation=45, ha='right', fontsize=26)
plt.yticks(fontsize=26)
plt.tight_layout()
plt.show()

sample_metadata_tmp = sample_metadata[sample_metadata['disease'] != 'normal']
## remove normal from xtciks
temp = sample_metadata_tmp[column_name].value_counts()
### remove normal
temp = temp[temp.index != 'normal']
plt.figure(figsize=(12, 6)) #(10, 6)
temp.plot(kind='bar')
plt.xlabel(column_name, fontsize=26)
plt.ylabel("Count", fontsize=26)
plt.xticks(rotation=45, ha='right', fontsize=26)
plt.yticks(fontsize=26)
plt.tight_layout()
plt.show()

#### save subset data
#adata_sub.write("Data/human_prefrontal_cortex_HBCC_Cohort_BBB_cell_types.h5ad")
####################################################

####
columns_to_check = ['class', 'subclass',  'suspension_type',
    'genetic_ancestry', 'assay', 'sex', 'tissue']
columns_to_check = ['disease', 'cell_type']

# Make a bar plot for each column
for column_name in columns_to_check:
    if column_name in adata_sub2.obs.columns:
        plt.figure(figsize=(14, 10)) # (25, 18) for larger plots
        adata_sub2.obs[column_name].value_counts().plot(kind='bar')
        plt.title(f"Counts for {column_name}")
        plt.xlabel(column_name, fontsize=26)
        plt.ylabel("Count", fontsize=26)
        plt.xticks(rotation=45, ha='right', fontsize=28)
        plt.yticks(fontsize=26)
        plt.tight_layout()
        plt.show()
unique_donors = adata_sub2.obs['donor_id'].nunique()
print(f"Unique donor IDs: {unique_donors}")

### make sample-level metadata for donor_id and sex, genetic_ancestry, disease and then make bar plots for each
sample_metadata = adata_sub2.obs[['donor_id', 'sex', 'genetic_ancestry', 'disease']]
sample_metadata = sample_metadata.drop_duplicates()

for column_name in sample_metadata.columns:
    plt.figure(figsize=(16, 12)) #(10, 6)
    sample_metadata[column_name].value_counts().plot(kind='bar')
    plt.title(f"Counts for {column_name}")
    plt.xlabel(column_name, fontsize=26)
    plt.ylabel("Count", fontsize=26)
    plt.xticks(rotation=45, ha='right', fontsize=26)
    plt.yticks(fontsize=26)
    plt.tight_layout()
    plt.show()

##### save subset data
#adata_sub2.write("Data/human_prefrontal_cortex_HBCC_Cohort_BBB_cell_types_Schizophrenia_CaseControl.h5ad")
########################################################

### print the head of the data X
print(adata_sub2.X.shape)
print(adata_sub2.X[:5])  # X is count not normalized
sc.pp.filter_genes(adata_sub2, min_cells=3)
sc.pp.normalize_total(adata_sub2, target_sum=1e4)
sc.pp.log1p(adata_sub2)
sc.pp.highly_variable_genes(adata_sub2)#min_mean=0.0125, max_mean=3, min_disp=0.5
sc.pp.pca(adata_sub2, n_comps=50, svd_solver='arpack')
print(adata_sub2.obsm['X_pca'][:5])  # Print the first 5 rows of the PCA data

sc.pp.neighbors(adata_sub2, n_neighbors=10, use_rep='X_pca')
sc.tl.umap(adata_sub2)
sc.pl.umap(adata_sub2, color=['cell_type'], frameon=False, size=20, title='Cell Type UMAP', wspace=0.4, hspace=0.4)
sc.pl.umap(adata_sub2, color=['disease'], frameon=False, size=20, title='Disease UMAP', wspace=0.4, hspace=0.4)
#sc.pl.umap(adata_sub2, color=['donor_id'], frameon=False, size=20, title='Donor ID UMAP', wspace=0.4, hspace=0.4)


columns_to_check = ['class', 'subclass','genetic_ancestry', 'assay', 'sex', 'disease', 'cell_type']
for column_name in columns_to_check:
    if column_name in adata_sub2.obs.columns:
        sc.pl.umap(adata_sub2, color=[column_name], 
        frameon=False, size=20, title=f'{column_name} UMAP', wspace=0.4, hspace=0.4)
print(adata_sub2.shape)

#### save subset2 preprocessed data
adata_sub2.write("Data/human_prefrontal_cortex_HBCC_Cohort_BBB_cell_types_Schizophrenia_CaseControl_preprocessed.h5ad")


columns_to_check = ['class', 'subclass','genetic_ancestry', 'assay', 'sex', 'disease', 'cell_type']
for column_name in columns_to_check:
    if column_name in adata.obs.columns:
        sc.pl.umap(adata, color=[column_name], 
        frameon=False, size=20, title=f'{column_name} UMAP', wspace=0.4, hspace=0.4)
print(adata.shape)


adata = ad.read_h5ad("../Data/human_prefrontal_cortex_HBCC_Cohort_Schizophrenia_preprocessed.h5ad")
adata.obs['immune_cell'] = adata.obs['cell_type']
### if the cell_type is in not the list of immune cells, ste it to 'other'
immune_cell_types = ['B cell', 'plasma cell', 'T cell', 'natural killer cell']
adata.obs['immune_cell'] = adata.obs['cell_type'].apply(lambda x: x if x in immune_cell_types else 'other')
custom_palette = {
    'B cell': '#1f77b4',
    'plasma cell': '#ff7f0e',
    'T cell': '#2ca02c',
    'natural killer cell': '#d62728',
    'astrocyte': '#9467bd',
    'endothelial cell': '#8c564b',
    'pericyte': '#e377c2',
    #'other': '#7f7f7f'  # Grey for other cell types
    'other': '#D3D3D3'
}
sc.pl.umap(adata, color=['immune_cell'], frameon=False, 
           size=20, title='Immune cells Schizophrenia UMAP',
             palette=custom_palette, 
           wspace=0.4, hspace=0.4, alpha=0.5)
