import numpy as np
# import muon as mu
import scanpy as sc
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches

import seaborn as sns
import pandas as pd

import radarPlot
import scipy.stats as stats
import covax_constants
from datetime import datetime
import os

from IPython.display import display, HTML
import io
import base64



def cell_numbers_longitudinal_box_plot(adata_vc, restrict_to_subset=None, annots_col = None, title = None, figsize = (8, 4)):
    if restrict_to_subset is None:
        i = adata_vc.obs_names
        df_ = adata_vc.obs.loc[i, ['Vaccine', 'person', 'Visit']].value_counts()
        
    elif annots_col is not None:
        i = adata_vc.obs[annots_col] == restrict_to_subset
        df_ = adata_vc.obs.loc[i, ['Vaccine', 'person', 'Visit']].value_counts()
        # Not a perfect solution for not-removing-zero-counts
        all_combinations = adata_vc.obs.loc[:, ['Vaccine', 'person', 'Visit']].value_counts().index
        df_ = df_.reindex(all_combinations, fill_value=0)
    else:
        raise ValueError('Either pass both restrict_to_subset and annots_col or neither of them.')
    
    
    df_ = df_.reset_index()
    
    fig, axs = plt.subplots(ncols = 3, nrows = 1, width_ratios=[3, 5, 5], figsize = figsize)
    y_max = df_['count'].max()
    for ax, vax in zip(axs, adata_vc.obs.Vaccine.cat.categories):
        df_vax = df_.query('Vaccine == @vax')
        df_vax['person'] = df_vax['person'].astype(str)
        df_vax['Visit'] = df_vax['Visit'].cat.remove_unused_categories()
        df_vax = df_vax.sort_values(by = 'person')
        sns.boxplot(df_vax, x='Visit', y='count', color = 'Gray', fill=False, showfliers=False, saturation = 1, ax = ax)
        sns.stripplot(df_vax, x='Visit', y="count", hue='person', palette='tab10', marker = covax_constants.vax_markers[vax], s = 9, ax=ax)
        ax.tick_params(axis='x', rotation=90)
        ax.set_ylim(0, y_max*1.05)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
    
    axs[0].set_title('J&J', color = covax_constants.vax_colors['J&J'], fontsize= 18)
    axs[1].set_title('Moderna', color = covax_constants.vax_colors['Moderna'], fontsize= 18)
    axs[2].set_title('Pfizer', color = covax_constants.vax_colors['Pfizer'], fontsize= 18)
    axs[0].set_ylabel('Cell Counts')
    axs[0].set_xlabel('')
    axs[2].set_xlabel('')
    axs[1].set_xlabel('Visit')
    for ax in axs[1:]:
        ax.set_ylabel('')
        ax.set_yticks([])
        ax.spines['left'].set_visible(False)
    
    
    # Merge legends into one
    handles, labels = [], []
    for ax in axs:
        handles_, labels_ = ax.get_legend_handles_labels()
        handles.extend(handles_)
        labels.extend(labels_)
        ax.get_legend().remove()
    
    # Create the merged legend
    fig.legend(handles, labels, loc='center left', bbox_to_anchor=(.91, 0.5), frameon = False, title = 'Donors')
    
    fig.subplots_adjust(wspace=0.1, top=0.8)
    title = 'Counts Per Donor/Visit' if title is None else title
    fig.suptitle(title, fontsize = 16);
    return fig


def compare_groups(adata, feature, groupby, order = None, use_raw=True, ax = None, legend_loc = None, x_tick_rot = 45, figsize = (8, 4), pval_kwargs = {}, connect_dots=False, hue='Vaccine', palette=covax_constants.vax_colors,):
    df = sc.get.obs_df(adata, [feature, 'Vaccine', 'person', groupby], use_raw=use_raw)
    
    df = df.groupby(['person', 'Vaccine', groupby]).mean().reset_index().dropna()
    
    df = pd.melt(df, id_vars=['person', 'Vaccine', groupby], value_vars=[feature])


    if ax is None:
        fig, ax = plt.subplots(figsize = figsize)
    else:
        fig = ax.get_figure()
    y_max = df['value'].max()
    y_min = df['value'].min()

    df['person'] = df['person'].astype(str)
    # df[groupby] = df[groupby].cat.remove_unused_categories()
    sns.boxplot(df, x=groupby, y='value', color = 'Gray', fill=False, showfliers=False, saturation = 1, order = order, ax = ax)
    sns.stripplot(df, x=groupby, y="value", hue=hue, palette=palette, s = 9, jitter=False, order = order, ax=ax)
    ax.tick_params(axis='x', rotation=x_tick_rot)
    ax.set_ylim(y_min*0.90, y_max*1.05)    
    fig.suptitle(f'Average Expression of {feature}', fontsize = 16);

    # Draw lines connecting the points for each person
    if connect_dots:
        for vax, vax_df in  df.groupby('Vaccine'):
            for person, group in vax_df.groupby('person'):
                group = group.sort_values(by=groupby)
                ax.plot(group[groupby], group['value'], marker='o', markersize=3, linewidth=0.8, color=covax_constants.vax_colors[vax], alpha=0.7)

    if pval_kwargs:
        stat_test = pval_kwargs.get("test", stats.ttest_rel)
        pairs = pval_kwargs['pairs']
        for a, b in pairs:
            a = df.loc[df[groupby] == a, 'value']
            b = df.loc[df[groupby] == b, 'value']
            res = stat_test(a, b)
            display(res)

    return fig

def remove_all_texts_from(figure):
    # Remove text elements from each axes
    for ax in figure.get_axes():
        # Remove axis titles and labels
        ax.set_title('')
        ax.set_xlabel('')
        ax.set_ylabel('')
        
        # Remove tick labels
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        
        # Remove text objects in the axes
        for text in ax.texts:
            text.remove()
        
        # Remove tick parameters
        ax.tick_params(axis='both', which='both', labelleft=False, labelright=False, labeltop=False, labelbottom=False)

    # Remove figure-wide text elements
    for text in figure.texts:
        text.remove()

def savefig_png_eps(fig, fname, use_timestamp=True, make_folder=False, show=True):
    
    if use_timestamp:
        timestamp = datetime.now().strftime("%Y_%m_%d_%H_%M")
        fname = f"{fname}_{timestamp}"

    if make_folder:
        try:
            os.mkdir(fname)
        except FileExistsError:
            pass
        path_png = f'{fname}/{fname}.png'
        path_eps = f'{fname}/{fname}.eps'
    else:
        path_png = f'{fname}.png'
        path_eps = f'{fname}.eps'

    fig.savefig(path_png,  bbox_inches='tight')
    if show:
        plt.show()
    remove_all_texts_from(fig)
    fig.savefig(path_eps)
    plt.close(fig)

def average_expression_level_longitudinal_box_plot_2(adata, feature, use_raw=True, figsize = (8, 4), width_ratios=[3, 5, 5]):
    df = sc.get.obs_df(adata, [feature, 'Vaccine', 'person', 'Visit'], use_raw=use_raw)
    
    df = df.groupby(['Vaccine', 'person', 'Visit']).mean().reset_index().dropna()
    
    df = pd.melt(df, id_vars=['Vaccine', 'person', 'Visit'], value_vars=[feature])

    fig, axs = plt.subplots(ncols = 3, nrows = 1, width_ratios=width_ratios, figsize = figsize)
    y_max = df['value'].max()
    y_min = df['value'].min()
    for ax, vax in zip(axs, adata.obs.Vaccine.cat.categories):
        df_vax = df.query('Vaccine == @vax')
        df_vax['person'] = df_vax['person'].astype(str)
        df_vax['Visit'] = df_vax['Visit'].cat.remove_unused_categories()
        df_vax = df_vax.sort_values(by = 'person')
        sns.boxplot(df_vax, x='Visit', y='value', hue='Visit', fill=True, palette=covax_constants.vax_visit_colors[vax].values(), showfliers=False, dodge=False, ax = ax, width=0.7, saturation = 1)
        # strip = sns.stripplot(df_vax, x='Visit', y='value', linewidth=1, size=0, legend=True, ax = ax)
        strip = sns.stripplot(data=df_vax, x='Visit', y='value', jitter=False, color='k', linewidth=1, size=3, legend=True, ax=ax)
       
        # Draw lines connecting the points for each person
        for person, group in df_vax.groupby('person'):
            group = group.sort_values(by='Visit')
            ax.plot(group['Visit'], group['value'], marker='o', markersize=3, linewidth=0.8, color='gray', alpha=0.7)


        ax.tick_params(axis='x', rotation=90)
        ax.set_ylim(y_min*0.90, y_max*1.05)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_title(vax, color = covax_constants.vax_colors[vax], fontsize= 18)
        
    axs[0].set_ylabel('Average Expression')
    axs[0].set_xlabel('')
    axs[2].set_xlabel('')
    axs[1].set_xlabel('Visit')
    for ax in axs[1:]:
        ax.set_ylabel('')
        ax.set_yticks([])
        ax.spines['left'].set_visible(False)
        
    fig.subplots_adjust(wspace=0.1, top=0.8)
    fig.suptitle(f'Average Expression of {feature}', fontsize = 16);
    return fig

def cell_proportions_longitudinal_box_plot_2(pbmc_vc, celltype_name, annots_col = 'pbmc_level_annotations_1', figsize = (8, 4), width_ratios=[3, 5, 5]):
    df = pbmc_vc.obs.loc[pbmc_vc.obs[annots_col] == celltype_name, ['Vaccine', 'person', 'Visit']].value_counts()
    df = 100*df/pbmc_vc.obs[['Vaccine', 'person', 'Visit']].value_counts()
    df = df.reindex(covax_constants.long_index, fill_value=0)
    df = df.reset_index()

    fig, axs = plt.subplots(ncols = 3, nrows = 1, width_ratios=width_ratios, figsize = figsize)
    y_max = df['count'].max()
    y_min = df['count'].min()
    
    # for ax, vax in zip(axs, pbmc_vc.obs.Vaccine.cat.categories):
    for ax, (vax, df_vax) in zip(axs, df.groupby('Vaccine', observed=True)):
        df_vax.fillna({'count': 0}, inplace = True)
        y_min = min(y_min, df_vax['count'].min())
        # df_vax['person'] = df_vax['person'].astype(str)
        # df_vax['Visit'] = df_vax['Visit'].cat.remove_unused_categories()
        df_vax = df_vax.sort_values(by = 'person')
        sns.boxplot(df_vax, x='Visit', y='count', hue='Visit', fill=True, palette=covax_constants.vax_visit_colors[vax].values(), showfliers=False, dodge=False, ax = ax, width=0.7, saturation = 1)
        strip = sns.stripplot(data=df_vax, x='Visit', y='count', jitter=False, color='k', linewidth=1, size=3, legend=True, ax=ax)
       
        # Draw lines connecting the points for each person
        for person, group in df_vax.groupby('person'):
            group['Visit'] = pd.Categorical(group['Visit'], categories = pbmc_vc.obs['Visit'].cat.categories)
            group['Visit'] = group['Visit'].cat.remove_unused_categories()
            group = group.sort_values(by='Visit')
            
            ax.plot(group['Visit'], group['count'], marker='o', markersize=3, linewidth=0.8, color='gray', alpha=0.7)


        ax.tick_params(axis='x', rotation=90)
        ax.set_ylim(y_min*0.90, y_max*1.05)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_title(vax, color = covax_constants.vax_colors[vax], fontsize= 18)
        
    axs[0].set_ylabel('Cell Percentages')
    axs[0].set_xlabel('')
    axs[2].set_xlabel('')
    axs[1].set_xlabel('Visit')
    for ax in axs[1:]:
        ax.set_ylabel('')
        ax.set_yticks([])
        ax.spines['left'].set_visible(False)
        
    fig.subplots_adjust(wspace=0.1, top=0.8)
    fig.suptitle(f'{celltype_name} Percentages Per Donor/Visit', fontsize = 16);
    return fig

def average_expression_level_longitudinal_box_plot(adata, feature, use_raw=True, figsize = (8, 4), width_ratios=[3, 5, 5]):
    df = sc.get.obs_df(adata, [feature, 'Vaccine', 'person', 'Visit'], use_raw=use_raw)
    
    df = df.groupby(['Vaccine', 'person', 'Visit']).mean().reset_index().dropna()
    
    df = pd.melt(df, id_vars=['Vaccine', 'person', 'Visit'], value_vars=[feature])

    fig, axs = plt.subplots(ncols = 3, nrows = 1, width_ratios=width_ratios, figsize = figsize)
    y_max = df['value'].max()
    y_min = df['value'].min()
    for ax, vax in zip(axs, adata.obs.Vaccine.cat.categories):
        df_vax = df.query('Vaccine == @vax')
        df_vax['person'] = df_vax['person'].astype(str)
        df_vax['Visit'] = df_vax['Visit'].cat.remove_unused_categories()
        df_vax = df_vax.sort_values(by = 'person')
        sns.boxplot(df_vax, x='Visit', y='value', color = 'Gray', fill=False, showfliers=False, saturation = 1, ax = ax)
        sns.stripplot(df_vax, x='Visit', y="value", hue='person', palette='tab10', marker = covax_constants.vax_markers[vax], s = 9, ax=ax)
        ax.tick_params(axis='x', rotation=90)
        ax.set_ylim(y_min*0.90, y_max*1.05)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    axs[0].set_title('J&J', color = covax_constants.vax_colors['J&J'], fontsize= 18)
    axs[1].set_title('Moderna', color = covax_constants.vax_colors['Moderna'], fontsize= 18)
    axs[2].set_title('Pfizer', color = covax_constants.vax_colors['Pfizer'], fontsize= 18)
    axs[0].set_ylabel('Average Expression')
    axs[0].set_xlabel('')
    axs[2].set_xlabel('')
    axs[1].set_xlabel('Visit')
    for ax in axs[1:]:
        ax.set_ylabel('')
        ax.set_yticks([])
        ax.spines['left'].set_visible(False)
        


    # Merge legends into one
    handles, labels = [], []
    for ax in axs:
        handles_, labels_ = ax.get_legend_handles_labels()
        handles.extend(handles_)
        labels.extend(labels_)
        if ax.get_legend():
            ax.get_legend().remove()
    
    # Create the merged legend
    fig.legend(handles, labels, loc='center left', bbox_to_anchor=(.91, 0.5), frameon = False, title = 'Donors')
    
    fig.subplots_adjust(wspace=0.1, top=0.8)
    fig.suptitle(f'Average Expression of {feature}', fontsize = 16);
    return fig


def cell_proportions_longitudinal_box_plot(pbmc_vc, celltype_name, annots_col, figsize = (8, 4)):
    df_ = pbmc_vc.obs.loc[pbmc_vc.obs[annots_col] == celltype_name, ['Vaccine', 'person', 'Visit']].value_counts()
    df_ = 100*df_/pbmc_vc.obs[['Vaccine', 'person', 'Visit']].value_counts()
    df_ = df_.reindex(covax_constants.long_index, fill_value=0)
    df_ = df_.reset_index()
    
    fig, axs = plt.subplots(ncols = 3, nrows = 1, width_ratios=[3, 5, 5], figsize = figsize)
    y_max = df_['count'].max()
    for ax, vax in zip(axs, pbmc_vc.obs.Vaccine.cat.categories):
        df_vax = df_.query('Vaccine == @vax')
        df_vax['person'] = df_vax['person'].astype(str)
        # df_vax['Visit'] = df_vax['Visit'].cat.remove_unused_categories()
        df_vax = df_vax.sort_values(by = 'person')
        df_vax = df_vax.fillna(0)
        sns.boxplot(df_vax, x='Visit', y='count', color = 'Gray', fill=False, showfliers=False, saturation = 1, ax = ax)
        sns.stripplot(df_vax, x='Visit', y="count", hue='person', palette='tab10', marker = covax_constants.vax_markers[vax], s = 9, ax=ax)
        ax.tick_params(axis='x', rotation=90)
        ax.set_ylim(0, y_max*1.05)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    
    axs[0].set_title('J&J', color = covax_constants.vax_colors['J&J'], fontsize= 18)
    axs[1].set_title('Moderna', color = covax_constants.vax_colors['Moderna'], fontsize= 18)
    axs[2].set_title('Pfizer', color = covax_constants.vax_colors['Pfizer'], fontsize= 18)
    axs[0].set_ylabel('Cell Percentages')
    axs[0].set_xlabel('')
    axs[2].set_xlabel('')
    axs[1].set_xlabel('Visit')
    for ax in axs[1:]:
        ax.set_ylabel('')
        ax.set_yticks([])
        ax.spines['left'].set_visible(False)
    
    
    # Merge legends into one
    handles, labels = [], []
    for ax in axs:
        handles_, labels_ = ax.get_legend_handles_labels()
        handles.extend(handles_)
        labels.extend(labels_)
        ax.get_legend().remove()
    
    # Create the merged legend
    fig.legend(handles, labels, loc='center left', bbox_to_anchor=(.91, 0.5), frameon = False, title = 'Donors')
    
    fig.subplots_adjust(wspace=0.1, top=0.8)
    fig.suptitle(f'{celltype_name} Percentages Per Donor/Visit', fontsize = 16);
    return fig



def donor_level_vaccine_comparison_box_plot(adata, feature, the_test = None, exclude_persons = []):
    def vaccines_pairwise_stats(df, the_test, feature,):
        r1 = the_test(df.loc[df.Vaccine=='J&J', feature], df.loc[df.Vaccine=='Moderna', feature])
        r2 = the_test(df.loc[df.Vaccine=='J&J', feature], df.loc[df.Vaccine=='Pfizer', feature])
        r3 = the_test(df.loc[df.Vaccine=='Moderna', feature], df.loc[df.Vaccine=='Pfizer', feature])
        r4 = the_test(df.loc[df.Vaccine=='J&J', feature], df.loc[df.Vaccine!='J&J', feature])
        return {'J&J vs Moderna' : r1, 'J&J vs Pfizer' : r2, 'Moderna vs Pfizer' : r3, 'J&J vs Pfizer + Moderna' : r4}
    df = sc.get.obs_df(adata, [feature, 'Vaccine', 'person'], use_raw=True)
    df = df.groupby(['Vaccine', 'person']).mean().reset_index().dropna()
    df = df[~df.person.isin(exclude_persons)]
    df['Person'] = df['Vaccine'].astype(str) + ' ' + df['person'].astype(str)
    df_mRNA = df[df['Vaccine']!='J&J'].copy()
    df_mRNA['Vaccine'] = 'Pfizer + Moderna'
    df_mRNA['Person'] = '-'
    df = pd.concat([df, df_mRNA], axis=0)
    res_dict = vaccines_pairwise_stats(df, the_test = the_test, feature = feature)
    ax = sns.boxplot(df, x = 'Vaccine', y = feature, showfliers=False)
    sns.stripplot(df, x = 'Vaccine', y = feature, hue = 'Person', palette='tab20', ax = ax)
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5));
    ax.set_title(feature)
    return ax, res_dict




    
