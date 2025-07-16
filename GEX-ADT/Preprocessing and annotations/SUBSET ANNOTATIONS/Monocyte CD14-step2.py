import my_plotting
import markers
import importlib
import utils
import covax_constants
import CovidVAX_specific_utils
import GeneModules


import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy import stats
import numpy as np

import warnings
warnings.simplefilter("ignore")


import scanpy as sc


ct = 'Monocyte CD14'
adata = sc.read_h5ad(f'{ct}-1.h5ad')

adata = adata[adata.obs.new_subset_annotations != 'LowQC']

adata_isg = adata[adata.obs.new_subset_annotations == 'Monocyte CD14 ISG'].copy()

adata_isg = adata_isg.raw.to_adata()

adata_isg.raw = adata_isg.copy()

adata_isg = adata_isg[:, adata_isg.var.is_rna]

sc.pp.highly_variable_genes(adata_isg, batch_key='library')

adata_isg = adata_isg[:, adata_isg.var.highly_variable]

sc.pp.scale(adata_isg, max_value=10)

sc.tl.pca(adata_isg, svd_solver='arpack')

sc.external.pp.bbknn(adata_isg, batch_key='library', n_pcs=20)

sc.tl.umap(adata_isg)

sc.tl.leiden(adata_isg, resolution = 0.4, key_added='leiden_040')

sc.tl.leiden(adata_isg, resolution = 0.7, key_added='leiden_040_1', restrict_to=('leiden_040', ['0']))


adata_isg.obs['new_subset_annotations_1'] = adata_isg.obs.leiden_040_1.map({
    '0,0': 'Monocyte CD14 ISG 2',
    '0,1': 'Monocyte CD14 ISG 2',
    '0,2': 'Monocyte CD14 ISG 2',
    '0,3': 'Monocyte CD14 ISG 2',
    '0,4': 'Monocyte CD14 ISG 2',
    '0,5': 'Monocyte CD14 ISG 2',
    '0,6': 'Monocyte CD14 ISG Inflamm.',
    '1': 'Monocyte CD14 ISG 1',
    '2': 'Doublet'
})


sc.pl.umap(adata_isg, color= ['new_subset_annotations_1', 'IL1B', 'ISG15', 'IFITM3'])
plt.savefig('adata_isg_umap.png')




adata_nonisg = adata[adata.obs.new_subset_annotations == 'Monocyte CD14'].copy()

adata_nonisg = adata_nonisg.raw.to_adata()

adata_nonisg.raw = adata_nonisg.copy()

adata_nonisg = adata_nonisg[:, adata_nonisg.var.is_rna]

sc.pp.highly_variable_genes(adata_nonisg, batch_key='library')

adata_nonisg = adata_nonisg[:, adata_nonisg.var.highly_variable]

sc.pp.scale(adata_nonisg, max_value=10)

sc.tl.pca(adata_nonisg, svd_solver='arpack')

sc.external.pp.bbknn(adata_nonisg, batch_key='library', n_pcs=15)

sc.tl.umap(adata_nonisg)

sc.tl.leiden(adata_nonisg, resolution = .5, key_added = 'leiden_050')

adata_nonisg.obs['new_subset_annotations_1'] = adata_nonisg.obs.leiden_050.map({
    '0': 'Monocyte CD14',
    '1': 'Monocyte CD14',
    '2': 'Monocyte CD14',
    '3': 'Monocyte CD14',
    '4': 'Monocyte CD14 ISG (Low)',
    '5': 'Doublet'
})

sc.pl.umap(adata_nonisg, color = ['leiden_050', 'new_subset_annotations_1', 'PID1', 'XAF1', 'n_genes_by_counts'])
plt.savefig('adata_nonisg_umap.png')


pbmc = sc.read_h5ad('adata-rna-adt-joint-1.h5ad')

pbmc.obs['new_subset_annotations_2'] = pbmc.obs['new_subset_annotations_1'].astype(str)
pbmc.obs.loc[adata_isg.obs_names, 'new_subset_annotations_2'] = adata_isg.obs.new_subset_annotations_1
pbmc.obs.loc[adata_nonisg.obs_names, 'new_subset_annotations_2'] = adata_nonisg.obs.new_subset_annotations_1

monocyte = pbmc[pbmc.obs.new_subset_annotations_2.str.startswith('Monocyte')].copy()


CD14 = monocyte[monocyte.obs.new_subset_annotations_2.str.startswith('Monocyte CD14')].copy()

adata_c = CD14.copy()

adata_c = adata_c.raw.to_adata()

adata_c.raw = adata_c.copy()

adata_c = adata_c[:, adata_c.var.is_rna]

sc.pp.highly_variable_genes(adata_c, batch_key='library')

adata_c = adata_c[:, adata_c.var.highly_variable]

sc.pp.scale(adata_c, max_value=10)

sc.tl.pca(adata_c, svd_solver='arpack')

sc.external.pp.bbknn(adata_c, batch_key='library', n_pcs=25)

sc.tl.umap(adata_c)

sc.external.pp.harmony_integrate(adata_c, 'library')

pcs = [10, 15, 25]
perplexities = [20, 30, 50, 100, 200, 500, 800, 1000]


for n_pcs in pcs:
    for perplexity in perplexities:
        sc.tl.tsne(adata_c, use_rep='X_pca_harmony', n_pcs = n_pcs, perplexity=perplexity, learning_rate=1000)
        adata_c.obsm[f'X_tsne_{n_pcs}_{perplexity}'] = adata_c.obsm['X_tsne'].copy()
        print(f'{n_pcs}, {perplexity} done.')
        

adata_c.write_h5ad('Monocyte CD14-2.h5ad')

















