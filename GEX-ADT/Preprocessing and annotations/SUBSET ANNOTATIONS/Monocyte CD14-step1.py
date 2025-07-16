import scanpy as sc

ct = 'Monocyte CD14'
adata = sc.read_h5ad(f'{ct}.h5ad')

sc.tl.leiden(adata, key_added='subset_leiden_060', resolution = 0.6)
sc.tl.umap(adata)

adata.obs['new_subset_annotations'] = adata.obs['subset_leiden_060'].map({
    '0': 'Monocyte CD14',
    '1': 'Monocyte CD14',
    '2': 'Monocyte CD14',
    '3': 'Monocyte CD14 ISG',
    '4': 'LowQC',
    '5': 'Monocyte CD14 ISG',
})


adata.write_h5ad(f'{ct}-1.h5ad')