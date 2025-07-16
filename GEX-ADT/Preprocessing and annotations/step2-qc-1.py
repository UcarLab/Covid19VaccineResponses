import muon as mu
import scanpy as sc
import pandas as pd
import numpy as np
import anndata
import os


print(os.getcwd())

mudata = mu.read_h5mu('combined_raw_mudata.h5mu')

######################################### RNA

rna = mudata['rna'].copy()

rna.var['mt'] = rna.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(rna, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

rna.obs[mudata.obs.columns] = mudata.obs

mu.pp.filter_var(rna, 'n_cells_by_counts', lambda x: x >= 3)
mu.pp.filter_obs(rna, 'n_genes_by_counts', lambda x: (x >= 200) & (x < 3000))

sc.pp.normalize_total(rna)
sc.pp.log1p(rna)

sc.pp.highly_variable_genes(rna, batch_key='library')

rna.raw = rna.copy()

sc.pp.scale(rna, max_value=10)

rna = rna[:, rna.var.highly_variable]

sc.tl.pca(rna, svd_solver='arpack')
sc.external.pp.bbknn(rna, batch_key='library', n_pcs=30)
sc.tl.leiden(rna, key_added = 'leiden_060', resolution=0.6)
sc.tl.leiden(rna, key_added = 'leiden_100', resolution=1.0)
sc.tl.umap(rna)

mudata.mod['rna'] = rna.copy()

print('RNA done.')


######################################### ADT

adt = mudata.mod['adt'].copy()

adt.obs[mudata.obs.columns] = mudata.obs

mu.prot.pp.clr(adt)
sc.pp.normalize_total(adt)

# Hack to exclude not-in-common proteins from pca
adt.var['highly_variable'] = (~adt.var['is_control'] & adt.var['is_in_common'])

# Save normalized but not scaled
adt.raw = adt.copy()

sc.pp.scale(adt, max_value=10)

sc.tl.pca(adt, svd_solver='arpack')
sc.external.pp.bbknn(adt, batch_key='library', n_pcs=20)
sc.tl.umap(adt)
sc.tl.leiden(adt, key_added = 'leiden_060', resolution=0.6)
sc.tl.leiden(adt, key_added = 'leiden_100', resolution=1.0)


mudata.mod['adt'] = adt

print('Protein done.')

######################################### Save

mudata.update()

mu.pp.intersect_obs(mudata)

mudata.write_h5mu('mudata-qc-1.h5mu')



