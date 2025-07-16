import pandas as pd
from scipy import stats
import scanpy as sc
import anndata
from typing import Union, Literal, Optional, Iterable, List
from statsmodels.stats.multitest import multipletests
import warnings
import numpy as np

def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier


def donor_level_diff_analysis(adata, features, groupby, order = None, use_raw=True, pairs = [], layer=None, statistical_test = stats.ttest_rel, pvalue_adj_method = 'fdr_bh'):
    from statsmodels.stats.multitest import multipletests
    # This assumes paired samples, does not work properly for comparing different people!
    all_people_here = adata.obs['person'].unique()
    df = sc.get.obs_df(adata, [*features, 'person', groupby], use_raw=use_raw, layer=layer)
    
    has_negative_values = df[[*features]].lt(0).any(axis=None)
    if has_negative_values:
        warnings.warn('There are negative values in the matrix. This method assumes non-negative values!', UserWarning)
        
    df = df.groupby([groupby, 'person'], observed=False).mean().reset_index().dropna()

    results = {}   
    
    for a, b in pairs:
        if a not in df[groupby].cat.categories:
            raise ValueError(f'{a} is not in the {groupby} column.')
        if b not in df[groupby].cat.categories:
            raise ValueError(f'{b} is not in the {groupby} column.')
        
        df_a = df[df[groupby] == a].set_index('person').sort_index()
        df_b = df[df[groupby] == b].set_index('person').sort_index()
    
        genes1 = df_a[list(features)]
        genes1 = genes1.reindex(all_people_here, fill_value=0)

        genes2 = df_b[list(features)]
        genes2 = genes2.reindex(all_people_here, fill_value=0)
        
        res = statistical_test(genes1, genes2, axis = 0)
        fc = (genes1.mean(axis = 0))/(genes2.mean(axis = 0))

        stat_fc_inconsistent = (res.statistic * (fc-1) < 0).any()
        
        if stat_fc_inconsistent:
            warnings.warn(f'The statistic and fold change are inconsistent in {a} vs {b}', UserWarning)


        result_df = pd.DataFrame({'statistic' : res.statistic, 'pvalue' : res.pvalue, 'foldchange' : fc}, index = features)
        nan_features = result_df.isna().any(axis = 'columns')

        result_df['pvalue_adj'] = np.nan
        result_df.loc[~nan_features, 'pvalue_adj'] = multipletests(result_df.loc[~nan_features, 'pvalue'], method=pvalue_adj_method)[1]

        results[f'{a} vs {b}'] = result_df
    
    return pd.concat(results, axis=1)


def filter_genes(adata: anndata.AnnData, 
                 genes: Union[Literal['all', 'all_reasonable'], Iterable[str]]='all_reasonable', 
                 perc_exp_min: int = 15,
                 ubiquitous_exp_max_perc: int = 85,
                 min_avg_expr: int = 0.1, 
                 groupby: Optional[str] = None, 
                 use_raw: bool = True,
                 layer: Optional[str]=None,
                ):
    
    adata = adata.raw.to_adata() if use_raw else adata

    if type(genes) == list:
        genes = pd.Index(genes)
    if type(genes) ==str:
        if genes == 'all':
            genes = adata.var_names
        elif genes == 'all_reasonable':
            genes = adata.var_names
            mask = ~(genes.str.startswith('MT-') | 
                     genes.str.startswith('LINC') |
                     genes.str.startswith('RPL') |
                     genes.str.startswith('RPS') |
                     genes.str.contains('.', regex = False) |
                     genes.str.contains('_', regex = False))
            genes = genes[mask]
        else:
            raise ValueError('Invalid value for the parameter `genes`')
    feats = [*genes, groupby] if groupby else genes
    df = sc.get.obs_df(adata, list(feats), layer = layer)

    if groupby:
        means = df.groupby(groupby, observed=False).mean().max()
    else:
        means = df.mean()

    df[genes] = df[genes] > 0
    perc_expr_overall = df[genes].mean() * 100
    if groupby:
        perc_expr = df.groupby(groupby, observed=False).mean().max() * 100
    else:
        perc_expr = perc_expr_overall


    mask1 = means > min_avg_expr
    mask2 = perc_expr > perc_exp_min
    mask3 = perc_expr_overall < ubiquitous_exp_max_perc

    the_stats = pd.concat({'Percentage_Expressed' : perc_expr, 'Mean_GE' : means}, axis=1)
    return genes[mask1 & mask2 & mask3], the_stats

def calculate_pvalues_between_visits_vectorized(vaccine, visit1, visit2, adata, features, use_raw=True, statistical_test = stats.ttest_rel):
    _adata = adata[(adata.obs.Vaccine==vaccine) & adata.obs.Visit.isin({visit1, visit2})]
    df = sc.get.obs_df(_adata, [*features, 'Vaccine', 'person', 'Visit'], use_raw=use_raw)
    
    df = df.groupby(['Vaccine', 'person', 'Visit'], observed=False).mean().reset_index().dropna()
    
    df1 = df[df.Visit == visit1].set_index('person').sort_index()
    df2 = df[df.Visit == visit2].set_index('person').sort_index()

    genes1 = df1[list(features)]
    genes2 = df2[list(features)]
    
    res = statistical_test(genes1, genes2, axis = 0)
    
    return pd.DataFrame(res, index = ['statistic', 'p-value'], columns=features)

def calculate_pvalues_between_visits(vaccine, visit1, visit2, adata, feature, use_raw=True, statistical_test = stats.ttest_rel):
    
    df = sc.get.obs_df(adata, [feature, 'Vaccine', 'person', 'Visit'], use_raw=use_raw)
    
    df = df.groupby(['Vaccine', 'person', 'Visit'], observed=False).mean().reset_index().dropna()
    
    df = pd.melt(df, id_vars=['Vaccine', 'person', 'Visit'], value_vars=[feature])

    df1 = df.query('(Vaccine == @vaccine) and (Visit == @visit1)').set_index('person')
    df2 = df.query('(Vaccine == @vaccine) and (Visit == @visit2)').set_index('person')

    res = statistical_test(df1.value.values, df2.value.values)
  
    return res.statistic, res.pvalue


def make_combined_column(adata, columns_list, sep = ' '):
    new_col = adata.obs[columns_list[0]].astype(str)
    for col_name in columns_list[1:]:
        new_col = new_col + sep + adata.obs[col_name].astype(str)

    df = adata.obs[columns_list]
    categories = []
    for col in df:
        if df[col].dtype.name == 'category':
            categories.append(df[col].cat.categories)
        else:
            categories.append(df[col].astype(str).astype('category').cat.categories)
            
    combined_categories = pd.MultiIndex.from_product(categories, names=columns_list)
    combined_categories = combined_categories.map(lambda x: sep.join(x))
    return pd.Categorical(new_col, categories=combined_categories)


def clean_visits(adata):
    pfizer_i = (adata.obs['Vaccine'] == 'Pfizer') & (adata.obs['visit'].isin({'V1', 'V2', 'V4', 'V5', 'V6'}))
    moderna_i = (adata.obs['Vaccine'] == 'Moderna') & (adata.obs['visit'].isin({'V1', 'V2', 'V4', 'V5', 'V6'}))
    jnj_i = (adata.obs['Vaccine'] == 'J&J') & (adata.obs['visit'].isin({'V1', 'V2', 'V3'}))
    valid_i = pfizer_i | moderna_i | jnj_i
    return adata[valid_i].copy()

def make_Visits_col(adata):
    adata.obs['Visit'] = adata.obs['visit'].map({
        'V1': 'Baseline 1',
        'V2': 'Vac1 D1',
        'V3': 'Vac1 D7',
        'V4': 'Baseline 2',
        'V5': 'Vac2 D1',
        'V6': 'Vac2 D7',
        'V7': 'Day70', # From the start
        'V8': 'Day180', # From the start
    })

def kite_output_as_adata(folder, postprocessed=False):
    import scanpy as sc
    if not postprocessed:
        adt = sc.read_mtx(f'{folder}/featurecounts.mtx')
        features = pd.read_csv(f'{folder}/featurecounts.genes.txt', header=None)[0]
        barcodes = pd.read_csv(f'{folder}/featurecounts.barcodes.txt', header=None)[0]
        adt.var_names = [f[3:-6] if "IgG" not in f else f for f in features]
        adt.obs_names = barcodes
        return adt
    else:
        adt = sc.read_mtx(f'{folder}/matrix.mtx.gz')
        adt = adt.transpose()
        features = pd.read_csv(f'{folder}/features.tsv.gz', header=None)[0]
        barcodes = pd.read_csv(f'{folder}/barcodes.tsv.gz', header=None)[0]
        adt.var_names = [f[3:-6] if "IgG" not in f else f for f in features]
        adt.obs_names = barcodes
        return adt

