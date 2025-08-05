import anndata as ad
import pandas as pd
import matplotlib as plt

adata = ad.read_h5ad("Data/human_prefrontal_cortex_HBCC_Cohort.h5ad")
print(adata.shape)
metadata = pd.DataFrame(adata.obs)
### save metadata
metadata.to_csv("Data/human_prefrontal_cortex_HBCC_Cohort_metadata.csv")
column_name = 'Schizophrenia'
print(metadata[column_name].value_counts())
adata_subset = adata[adata.obs[column_name] == 'Yes']


#adata_subset.write("Data/human_prefrontal_cortex_HBCC_Cohort_Schizophrenia.h5ad")
adata_sub = ad.read_h5ad("Data/human_prefrontal_cortex_HBCC_Cohort_Schizophrenia.h5ad")
print(adata_sub.shape)
adata_sub.obs
metadata = pd.DataFrame(adata_sub.obs)

for column_name in adata_sub.obs.columns:
    print(f"Counts for {column_name}:")
    if adata_sub.obs[column_name].dtype == 'object':
        print(adata_sub.obs[column_name].value_counts())
    else:
        print("This column is not categorical or has no value counts.")

#'tissue_type'=tissue! al
columns_to_check = [
    'class', 'subclass',  'suspension_type',
    'genetic_ancestry', 'disease_ontology_term_id',  'assay', , 'sex', 'tissue'
]
columns_to_check = ['disease', 'cell_type']
# Make a bar plot for each column
for column_name in columns_to_check:
    if column_name in adata_sub.obs.columns:
        plt.figure(figsize=(25, 15))  # corrected from plt.pyplot.figure to plt.figure
        adata_sub.obs[column_name].value_counts().plot(kind='bar')
        plt.title(f"Counts for {column_name}")
        plt.xlabel(column_name, fontsize=16)
        plt.ylabel("Count", fontsize=16)
        plt.xticks(rotation=45, ha='right', fontsize=13)
        plt.yticks(fontsize=16)
        plt.tight_layout()
        plt.show()
## count unique 'donor_id'
unique_donors = adata_sub.obs['donor_id'].nunique()
print(f"Unique donor IDs: {unique_donors}")

### make sample-level metadata for donor_id and sex, genetic_ancestry, disease and then make bar plots for each
sample_metadata = adata_sub.obs[['donor_id', 'sex', 'genetic_ancestry', 'disease']]
sample_metadata = sample_metadata.drop_duplicates()
for column_name in sample_metadata.columns:
    plt.figure(figsize=(10, 10)) #(10, 6)
    sample_metadata[column_name].value_counts().plot(kind='bar')
    plt.title(f"Counts for {column_name}")
    plt.xlabel(column_name, fontsize=18)
    plt.ylabel("Count", fontsize=16)
    plt.xticks(rotation=45, ha='right', fontsize=13)
    plt.yticks(fontsize=16)
    plt.tight_layout()
    plt.show()

### draw the UMAP colored by cell type
import scanpy as sc
sc.pp.neighbors(adata_sub, n_neighbors=10, use_rep='X_pca')
sc.tl.umap(adata_sub)
sc.pl.umap(adata_sub, color=['cell_type'], frameon=False, size=20, title='Cell Type UMAP', legend_loc='on data', wspace=0.4, hspace=0.4, fontsize=16)
sc.pl.umap(adata_sub, color=['disease'], frameon=False, size=20, title='Disease UMAP', legend_loc='on data', wspace=0.4, hspace=0.4, fontsize=16)
sc.pl.umap(adata_sub, color=['donor_id'], frameon=False, size=20, title='Donor ID UMAP', legend_loc='on data', wspace=0.4, hspace=0.4, fontsize=16)
