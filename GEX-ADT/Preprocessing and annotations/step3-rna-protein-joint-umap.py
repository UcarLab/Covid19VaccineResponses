import muon as mu
import scanpy as sc
import anndata
import os

print(os.getcwd())

mudata = mu.read_h5mu('mudata-qc-1.h5mu')

rna = mudata['rna'].copy()
adt = mudata['adt'].copy()

rna = rna.raw.to_adata()
adt = adt.raw.to_adata()

rna.var['is_protein'] = False
rna.var['is_rna'] = True

adt.var['is_protein'] = True
adt.var['is_rna'] = False

combined = anndata.concat([rna, adt], axis = 1)

combined.obs[rna.obs.columns] = rna.obs.copy()

combined.raw = combined.copy()

combined = combined[:, combined.var.highly_variable]

sc.pp.scale(combined, max_value=10)
sc.tl.pca(combined, svd_solver='arpack')

sc.external.pp.bbknn(combined, batch_key='pool', n_pcs=50) 
sc.tl.umap(combined)
sc.tl.leiden(combined, key_added = 'leiden_060', resolution=0.6)
sc.tl.leiden(combined, key_added = 'leiden_080', resolution=0.8)
sc.tl.leiden(combined, key_added = 'leiden_100', resolution=1.0)

combined.write('adata-rna-adt-joint-1.h5ad')

