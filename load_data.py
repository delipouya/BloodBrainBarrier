import anndata as ad
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc

adata = ad.read_h5ad("Data/human_prefrontal_cortex_HBCC_Cohort.h5ad")
print(adata.shape)
metadata = pd.DataFrame(adata.obs)
metadata.head()

### save metadata
#metadata.to_csv("Data/human_prefrontal_cortex_HBCC_Cohort_metadata.csv")
#column_name = 'Schizophrenia'
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
## subset to only BBB cell types
adata_sub = adata[adata.obs['cell_type'].isin(BBB_cell_types)]
print(adata.shape)
print(adata_sub.shape)

columns_to_check = ['class', 'subclass',  'suspension_type',
    'genetic_ancestry', 'assay', 'sex', 'tissue']
columns_to_check = ['disease', 'cell_type']

# Make a bar plot for each column
for column_name in columns_to_check:
    if column_name in adata_sub.obs.columns:
        plt.figure(figsize=(12, 10)) # (25, 18) for larger plots
        adata_sub.obs[column_name].value_counts().plot(kind='bar')
        plt.title(f"Counts for {column_name}")
        plt.xlabel(column_name, fontsize=26)
        plt.ylabel("Count", fontsize=26)
        plt.xticks(rotation=45, ha='right', fontsize=28)
        plt.yticks(fontsize=26)
        plt.tight_layout()
        plt.show()
unique_donors = adata_sub.obs['donor_id'].nunique()
print(f"Unique donor IDs: {unique_donors}")

### make sample-level metadata for donor_id and sex, genetic_ancestry, disease and then make bar plots for each
sample_metadata = adata_sub.obs[['donor_id', 'sex', 'genetic_ancestry', 'disease']]
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


#### save subset data
#adata_sub.write("Data/human_prefrontal_cortex_HBCC_Cohort_BBB_cell_types.h5ad")
####################################################
column_name = 'Schizophrenia'
print(adata_sub.obs[column_name].value_counts())
condition = (adata_sub.obs[column_name] == 'Yes') | (adata_sub.obs['disease'] == 'normal')
adata_sub2 = adata_sub[condition]
print(adata_sub2.obs['disease'].value_counts())
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
