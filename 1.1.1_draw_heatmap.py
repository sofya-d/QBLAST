print('Importing modules')
import sys
import os
import pandas as pd
import numpy as np
from math import log
import seaborn as sns
from matplotlib import pyplot as plt


tsv_dirpath = sys.argv[1]
out_dirpath = sys.argv[2]


print('Defying functions')
def tsv_w_path(path):
    files = os.listdir(path)
    tsvs = [file for file in files if '.tsv' in file]
    tsvs.sort()
    tsvs_w_path = [f'{path}/{tsv}' for tsv in tsvs]
    return tsvs_w_path

def htmp_df_row(tsv_filepath: str):
    df = pd.read_csv(tsv_filepath, sep = '\t')
    assembly = tsv_filepath.split('/')[-1].split('.')[0]
    assembly_col = pd.DataFrame({'assembly': [assembly]})
    counts_col = pd.DataFrame(df.iloc[:, 4:].sum(axis = 1))
    counts_row = counts_col.transpose()
    counts_row = counts_row.set_axis(list(df['query']), axis = 1)
    counts_row = counts_row.set_axis([0], axis = 0)
    counts_row = pd.concat([assembly_col, counts_row], axis = 1)
    return counts_row

def heatmap_df_to_mx(heatmap_df):
    heatmap_mx = heatmap_df.iloc[:, 1:]
    heatmap_mx = heatmap_mx.set_axis(list(heatmap_df['assembly']), axis = 0)
    return heatmap_mx

def find_max_value(mx):
    row_max = mx.max(axis = 1)
    max_ = row_max.max()
    return max_

def log_scale_mx(mx, log_base):
    if log_base == 2:
        scaled_mx = (np.log2(heatmap_mx)).replace(-np.inf, np.nan)
    else:
        scaled_mx = (np.log10(heatmap_mx)).replace(-np.inf, np.nan)
    return scaled_mx

def make_legend_ticks_labels(max_value, log_base):
    max_pow = round(log(max_value, log_base))
    legend_ticks = []
    legend_tick_labels = []
    max_x = 2 if log_base == 2 else 5
    for power in range(max_pow + 1):
        for x in [1, max_x]:
            if log_base == 2:
                tick = log(x**power, log_base)
            else:
                tick = log(x*10**power, log_base)
            if tick not in legend_ticks:
                round_exponed = round(log_base**tick)
                legend_ticks.append(tick)
                legend_tick_labels.append(str(round_exponed))
    return legend_ticks, legend_tick_labels

def plot_heatmap(heatmap_mx, title, name, legend_ticks, legend_tick_labels, plot_dirpath):
    heatmap = sns.heatmap(heatmap_mx, cmap=sns.cubehelix_palette(start=.5, rot=-.75, as_cmap=True), 
                          linewidths = 0.4, cbar_kws={'label': 'accession number', 
                                                      'ticks': legend_ticks})
    heatmap.figure.set_size_inches(12, 10)
    heatmap.set(title = title, xlabel = 'query', ylabel = 'assembly')
    heatmap.set_facecolor('xkcd:light grey')
    heatmap.set_yticks(np.arange(0.5, heatmap_mx.shape[0])) # <--- set the ticks first
    heatmap.set_yticklabels(list(heatmap_mx.index)) # <--- set the ticks labels
    heatmap.collections[0].colorbar.set_ticklabels(legend_tick_labels)
    fig = heatmap.get_figure()
    fig.savefig(f'{out_dirpath}/{name}.svg', bbox_inches='tight')
    plt.close()


print('Making count matrix')
tsv_filepaths = tsv_w_path(tsv_dirpath)

tsv_dir = tsv_dirpath.split('/')[-2]
e_val_thresh = tsv_dir.split('_')[0]
len_thresh = tsv_dir.split('_')[1]
name = 'heatmap'

heatmap_df = pd.DataFrame()
for tsv_filepath in tsv_filepaths:
    counts_row = htmp_df_row(tsv_filepath)
    heatmap_df = pd.concat([heatmap_df, counts_row], ignore_index = True)

print('Writing count matrix as .tsv')
heatmap_df.to_csv(f'{out_dirpath}/{name}.tsv', sep = '\t', index = False)

heatmap_mx = heatmap_df_to_mx(heatmap_df)
max_value = find_max_value(heatmap_mx)

log_base = 2 if max_value <= 1000 else 10

print('Log-scaling count matrix')
heatmap_mx = log_scale_mx(heatmap_mx, log_base)
legend_ticks, legend_tick_labels = make_legend_ticks_labels(max_value, log_base)

title = f'BLAST accessions. e-value thresh. = {e_val_thresh}, length thresh. = {len_thresh}%'

print('Drawing and saving heatmap')
plot_heatmap(heatmap_mx, title, name, legend_ticks, legend_tick_labels, out_dirpath)