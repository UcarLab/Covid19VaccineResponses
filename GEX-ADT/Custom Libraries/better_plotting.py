import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import pandas as pd
from matplotlib.patches import Patch
from scipy import stats

def regression_plot(df, x, y, ax=None, legend_loc = 'upper left', hue=None, ci = None, annot_xy = (.3, .9)):
    
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.get_figure() # Did no test
    ax = sns.regplot(data=df.reset_index(), x=x, y=y, ci=ci, color = 'k', scatter=False, ax = ax)
    sns.scatterplot(df, x = x, y = y, hue = hue, ax = ax)
    
    res = stats.linregress(df[x], df[y])
    a_x, a_y = annot_xy
    ax.text(s = f'pval:{round(res.pvalue, 5)}', x = a_x, y = a_y, transform=ax.transAxes)
    ax.text(s = f'rval:{round(res.rvalue, 5)}', x = a_x, y = a_y - 0.05, transform=ax.transAxes)
    ax.spines[['right', 'top']].set_visible(False);
    return fig

def mean_expr_longitudinal_box_plot(adata, feature, visit_col = 'Visit', donor_col = 'person', use_raw=True, ax = None, legend_loc = None,
                   x_tick_rot = 45, figsize = (8, 4), pval_kwargs = {}, connect_dots=False, hue='person', palette=None, layer = None):

    df = sc.get.obs_df(adata, [feature, donor_col, visit_col], use_raw=use_raw, layer = layer)
    df = df.groupby([donor_col, visit_col], observed=True).mean().reset_index().dropna()
    df = pd.melt(df, id_vars=[donor_col, visit_col], value_vars=[feature])

    fig, ax = plt.subplots(figsize = figsize)

    df[donor_col] = df[donor_col].astype(str)

    df = df.sort_values(by = donor_col)
    sns.boxplot(df, x=visit_col, y='value', color = 'Gray', fill=False, showfliers=False, saturation = 1, ax = ax)
    sns.stripplot(df, x=visit_col, y="value", hue=donor_col, palette=palette, s = 9, jitter= .3, ax=ax)
    ax.tick_params(axis='x', rotation=90)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax.set_ylabel(f'{feature} - Mean Expression')
    ax.set_xlabel(visit_col)

    fig.suptitle(f'Mean {feature} Expression Per Donor/Visit', fontsize = 16);
    return fig

def cell_proportions_longitudinal_box_plot_2(adata, celltype_name, annots_col, figsize = (8, 4), donor_col = 'person', visit_col = 'Visit', palette=None):
    df = adata.obs.loc[adata.obs[annots_col] == celltype_name, [donor_col, visit_col]].value_counts()
    df = 100*df/adata.obs[[donor_col, visit_col]].value_counts()
    df = df.fillna(0)
    df = df.reset_index()

    
    fig, ax = plt.subplots(figsize = figsize)

    df[donor_col] = df[donor_col].astype(str)

    df = df.sort_values(by = donor_col)
    
    sns.boxplot(df, x=visit_col, y='count', hue=visit_col, fill=True, showfliers=False, dodge=False, ax = ax, width=0.7, saturation = 1, palette=palette)
    strip = sns.stripplot(df, x=visit_col, y='count', jitter=False, color='k', linewidth=1, size=3, legend=True, ax=ax)
    
    ax.tick_params(axis='x', rotation=90)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    ax.set_ylabel('Cell Percentages')
    ax.set_xlabel(visit_col)

    fig.suptitle(f'{celltype_name} Percentages Per Donor/Visit', fontsize = 16);
    return fig


def cell_proportions_longitudinal_box_plot(adata, celltype_name, annots_col, figsize = (8, 4), donor_col = 'person', visit_col = 'Visit', palette=None):
    df = adata.obs.loc[adata.obs[annots_col] == celltype_name, [donor_col, visit_col]].value_counts()
    df = 100*df/adata.obs[[donor_col, visit_col]].value_counts()
    df = df.fillna(0)
    df = df.reset_index()

    
    fig, ax = plt.subplots(figsize = figsize)

    df[donor_col] = df[donor_col].astype(str)

    df = df.sort_values(by = donor_col)
    sns.boxplot(df, x=visit_col, y='count', color = 'Gray', fill=False, showfliers=False, saturation = 1, ax = ax)
    sns.stripplot(df, x=visit_col, y="count", hue=donor_col, palette=palette, s = 9, jitter= .3, ax=ax)
    ax.tick_params(axis='x', rotation=90)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax.set_ylabel('Cell Percentages')
    ax.set_xlabel(visit_col)

    fig.suptitle(f'{celltype_name} Percentages Per Donor/Visit', fontsize = 16);
    return fig


def compare_groups(adata, feature, groupby, order = None, use_raw=True, ax = None, legend_loc = None,
                   x_tick_rot = 45, figsize = (8, 4), pval_kwargs = {}, connect_dots=False, hue='person', donor_col = 'person', palette=None, jitter=True, layer = None, dodge=False, groupby_observed=True):

    
    df = sc.get.obs_df(adata, [feature, donor_col, groupby], use_raw=use_raw, layer = layer)
    df = df.groupby([donor_col, groupby], observed=groupby_observed).mean().reset_index().dropna()
    df = pd.melt(df, id_vars=[donor_col, groupby], value_vars=[feature])

    if hue and hue!=donor_col:
        hues = sc.get.obs_df(adata, [donor_col, hue])
        hues = hues.drop_duplicates().set_index(donor_col).squeeze().to_dict()
        df[hue] = df[donor_col].map(hues)

    if ax is None:
        fig, ax = plt.subplots(figsize = figsize)
    else:
        fig = ax.get_figure()
    y_max = df['value'].max()
    y_min = df['value'].min()

    df[donor_col] = df[donor_col].astype(str)
    # df[groupby] = df[groupby].cat.remove_unused_categories()
    sns.boxplot(df, x=groupby, y='value', color = 'Gray', fill=False, showfliers=False, saturation = 1, order = order, ax = ax)
    sns.stripplot(df, x=groupby, y="value", hue=hue, palette=palette, s = 9, jitter=jitter, order = order, dodge=dodge, ax=ax)
    ax.tick_params(axis='x', rotation=x_tick_rot)
    ax.set_ylim(y_min*0.90, y_max*1.05)    
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    fig.suptitle(f'Average Expression of {feature}', fontsize = 16);

    # Draw lines connecting the points for each person
    if connect_dots:
        for person, group in df.groupby(donor_col):
            group = group.sort_values(by=groupby)
            ax.plot(group[groupby], group['value'], marker='o', markersize=3, linewidth=0.8, alpha=0.7)

    return fig