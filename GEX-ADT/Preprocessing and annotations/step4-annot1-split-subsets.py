import scanpy as sc

import pandas as pd
import numpy as np
import os
from os import path
from glob import glob

import anndata

adata = sc.read_h5ad('adata-rna-adt-joint-1.h5ad')

adata.obs['pbmc_level_annotations'] = adata.obs['leiden_080'].map({
    '0': 'T CD4 Naive',
    '1': 'T CD4 Memory',
    '2': 'Monocyte CD14',
    '3': 'T CD8 Memory',
    '4': 'T CD8 Naive',
    '5': 'NK',
    '6': 'B Naive',
    '7': 'Monocyte CD14',
    '8': 'B Memory',
    '9': 'T GD',
    '10': 'Monocyte CD16',
    '11': 'T CD4 Treg',
    '12': 'T Mait',
    '13': 'DC',
    '14': 'PC',
    '15': 'HSC'
})

######## Save the subsets

for ct in adata.obs['pbmc_level_annotations'].unique():
    subset = adata[adata.obs['pbmc_level_annotations'] == ct].copy()
    path = f'{ct}.h5ad'
    subset.write_h5ad(path)


