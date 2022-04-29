# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib import patches
import wot
from wot import *
import os
os.getcwd()
os.chdir("D:/Jupyter/scRNA")

## Load input data
FULL_DS_PATH = 'data/ExprMatrix.h5ad'
VAR_DS_PATH = 'data/ExprMatrix.var.genes.h5ad'
CELL_DAYS_PATH = 'data/cell_days.txt'
CELL_SETS_PATH = 'data/cell_sets.gmt'
COORDS_PATH = 'data/fle_coords.txt'

## Read expression matrix, cell days, 2-d coordinates, and cell sets
coord_df = pd.read_csv(COORDS_PATH, index_col='id', sep='\t')
days_df = pd.read_csv(CELL_DAYS_PATH, index_col='id', sep='\t')
adata = wot.io.read_dataset(FULL_DS_PATH, obs=[days_df,coord_df])
unique_days = adata.obs['day'].unique()
unique_days = unique_days[np.isnan(unique_days) == False]
cell_sets = wot.io.read_sets(CELL_SETS_PATH)

#%% 1. Data visualization
# Visualize cells in 2D
figure1 = plt.figure(figsize=(6, 4))
plt.axis('off')
plt.tight_layout()
plt.scatter(adata.obs['x'], adata.obs['y'],c=adata.obs['day'], 
            s=4, marker=',', edgecolors='none', alpha=0.8)
cb = plt.colorbar()
cb.ax.set_title('Day')

# Visualize cells group by phrases
adata_obs = adata.obs
adata_obs.loc[adata_obs.index.str.contains('Dox')==True, 'phrase'] = 'Dox'
adata_obs.loc[adata_obs.index.str.contains('2i')==True, 'phrase'] = '2i'
adata_obs.loc[adata_obs.index.str.contains('serum')==True, 'phrase'] = 'serum'

adata_ph1 = adata_obs.loc[adata_obs['phrase'] == 'Dox']
adata_ph2 = adata_obs.loc[adata_obs['phrase'] == '2i']
adata_ph3 = adata_obs.loc[adata_obs['phrase'] == 'serum']

figure2 = plt.figure(figsize=(6, 4))
plt.axis('off')
plt.tight_layout()
plt.scatter(coord_df['x'], coord_df['y'], c='lightgray',
            s=4, marker=',', edgecolors='none', alpha=0.6)
plt.scatter(adata_ph3['x'], adata_ph3['y'],
            s=4, marker=',', edgecolors='none', c='tab:red', alpha=0.6, label='Serum')
plt.scatter(adata_ph1['x'], adata_ph1['y'],
            s=4, marker=',', edgecolors='none', c='tab:gray', alpha=0.6, label='Dox')
plt.scatter(adata_ph2['x'], adata_ph2['y'],
            s=4, marker=',', edgecolors='none', c='tab:blue', alpha=0.6, label='2i')
plt.legend(fontsize=15)
plt.show()

# Visualize cells group by conditions
adata_c1 = adata_obs.loc[adata_obs['phrase'] != 'serum']
adata_c2 = adata_obs.loc[adata_obs['phrase'] != '2i']

figure3 = plt.figure(figsize=(6, 4))
plt.title('2i condition', fontsize=14)
plt.axis('off')
plt.tight_layout()
plt.scatter(coord_df['x'], coord_df['y'], c='lightgray',
            s=4, marker=',', edgecolors='none', alpha=0.9)
plt.scatter(adata_c1['x'], adata_c1['y'], c=adata_c1['day'],
            s=4, marker=',', edgecolors='none', alpha=0.6)
cb = plt.colorbar()
cb.ax.set_title('Day')

figure4 = plt.figure(figsize=(6, 4))
plt.title('Serum condition', fontsize=14)
plt.axis('off')
plt.tight_layout()
plt.scatter(coord_df['x'], coord_df['y'], c='lightgray',
            s=4, marker=',', edgecolors='none', alpha=0.9)
plt.scatter(adata_c2['x'], adata_c2['y'], c=adata_c2['day'],
            s=4, marker=',', edgecolors='none', alpha=0.6)
cb = plt.colorbar()
cb.ax.set_title('Day')


# Visualize major cell sets
for name in cell_sets.var.index:
    print(name)
    cell_set = cell_sets[:, name]
    cell_set_coords = cell_set[cell_set.X>0].obs.join(coord_df).join(days_df)
    figure = plt.figure(figsize=(5, 5))
    plt.axis('off')
    plt.tight_layout()
    plt.title(name, fontsize=16)
    plt.scatter(coord_df['x'], coord_df['y'], c='lightgray',
                s=4, marker=',', edgecolors='none', alpha=0.8)
    plt.scatter(cell_set_coords['x'], cell_set_coords['y'], c=cell_set_coords['day'],
                s=4, marker=',', edgecolors='none', vmin=unique_days[0], vmax=unique_days[len(unique_days)-1])
    cb = plt.colorbar()
    cb.ax.set_title('Day')
    # plt.savefig('{}.png'.format(name))


#%% 2. Trajectory analysis
# Compute the long-term coupling
ot_model = wot.ot.OTModel(adata, epsilon = 0.05, lambda1 = 1, lambda2 = 50, growth_iters = 3) 
ot_model.compute_all_transport_maps(tmap_out='tmaps/serum')
tmap_model = wot.tmap.TransportMapModel.from_directory('tmaps/serum')

for key in list(cell_sets.keys()):
    if len(cell_sets.get(key))<4000:
        cell_sets.pop(key)
populations = tmap_model.population_from_cell_sets(cell_sets, at_time=18)
trajectory_ds = tmap_model.trajectories(populations)

# Load embedding coordinates
coord_df = pd.read_csv(COORDS_PATH, sep='\t', index_col=0)
nbins = 500
xrange = coord_df['x'].min(), coord_df['x'].max()
yrange = coord_df['y'].min(), coord_df['y'].max()
coord_df['x'] = np.floor(np.interp(coord_df['x'], [xrange[0], xrange[1]], [0, nbins - 1])).astype(int)
coord_df['y'] = np.floor(np.interp(coord_df['y'], [yrange[0], yrange[1]], [0, nbins - 1])).astype(int)
trajectory_ds.obs = trajectory_ds.obs.join(coord_df)

# Visualize trajectories
for name in trajectory_ds.var.index:
    figure = plt.figure(figsize=(6, 5))
    plt.axis('off')
    plt.tight_layout()
    plt.title('{} ancestors'.format(name), fontsize=15)
    plt.scatter(coord_df['x'], coord_df['y'], c='lightgray',
                s=4, marker=',', edgecolors='none', alpha=0.8)
    binned_df = trajectory_ds.obs.copy()
    binned_df['values'] = trajectory_ds[:, name].X
    binned_df.loc[binned_df['values']<binned_df['values'].quantile(0.6), 'values'] = np.nan
    binned_df = binned_df.dropna(axis=0,how='any')
    binned_df = binned_df.groupby(['x', 'y'], as_index=False).sum()
    plt.scatter(binned_df['x'], binned_df['y'], c=binned_df['values'],
                s=6, marker=',', edgecolors='none', vmax=binned_df['values'].quantile(0.95))
    plt.colorbar().ax.set_title('Trajectory')


#%% 3. Shared ancestors
adata_var = wot.io.read_dataset(VAR_DS_PATH)

# Divergence type 1: Total variation
divergence_df = wot.tmap.trajectory_divergence(adata_var, trajectory_ds, distance_metric='total_variation')
divergence_df['name'] = divergence_df['name1'].str.split('/').str.get(0) + ' vs. ' + divergence_df['name2'].str.split('/').str.get(0)

plt.figure(figsize=(8, 8))
plt.xlabel("Day", fontsize=15)
plt.ylabel("Total variation", fontsize=15)
for p, d in divergence_df.groupby('name'):
    plt.plot(d['day2'], d['distance'], '-o', label=p)
plt.tick_params(labelsize=15)
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=15)
plt.show()

# Divergence type: 2 Wasserstein distance
divergence_df = wot.tmap.trajectory_divergence(adata_var, trajectory_ds, distance_metric='emd')
divergence_df['name'] = divergence_df['name1'].str.split('/').str.get(0) + ' vs. ' + divergence_df['name2'].str.split('/').str.get(0)
# divergence_df.to_csv('data/divergence_df.csv', index=False)
# divergence_df = pd.read_csv('data/divergence_df.csv')

plt.figure(figsize=(8, 8))
plt.xlabel("Day", fontsize=15)
plt.ylabel("Wasserstein distance", fontsize=15)
for p, d in divergence_df.groupby('name'):
    plt.plot(d['day2'], d['distance'], '-o', label=p)
plt.tick_params(labelsize=15)
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=15)
plt.show()


#%% 4. Transition table
start_time = 12
end_time = 18
start_populations = tmap_model.population_from_cell_sets(cell_sets, at_time=start_time)
end_populations = tmap_model.population_from_cell_sets(cell_sets, at_time=end_time)
transition_table = tmap_model.transition_table(start_populations, end_populations)

# Visualize transition table as heat map
fig, ax = plt.subplots(figsize=(6, 6))
im = ax.imshow(transition_table.X)
ax.set_xticks(np.arange(len(transition_table.var_names)))
ax.set_yticks(np.arange(len(transition_table.obs_names)))
ax.set_xticklabels(transition_table.var_names)
ax.set_yticklabels(transition_table.obs_names)

plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

for i in range(transition_table.shape[0]):
    for j in range(transition_table.shape[1]):
        text = ax.text(j, i, '{:.2f}'.format(transition_table.X[i, j]),
                       ha="center", va="center", color="w")
fig.tight_layout()
plt.title('From day {} to day {}'.format(start_time, end_time), fontsize=15)


#%% 5. Validation via Interpolation
SERUM_CELL_IDS_PATH = 'data/serum_cell_ids.txt'
adata = wot.io.read_dataset(VAR_DS_PATH, obs=[CELL_DAYS_PATH], obs_filter=SERUM_CELL_IDS_PATH)
ot_model = wot.ot.OTModel(adata, growth_rate_field='g2', growth_iters = 1)
triplets_stats = wot.ot.compute_validation_summary(ot_model, day_triplets=[(0.5, 1, 1.5), (1.5, 2, 2.5), (2.5, 3, 3.5), (3.5, 4, 4.5),
                                                                           (4.5, 5, 5.5), (5.5, 6, 6.5), (6.5, 7, 7.5), (7.5, 8, 8.5),
                                                                           (8.5, 9, 9.5), (9.5, 10, 10.5), (10.5, 11, 11.5), (11.5, 12, 12.5),
                                                                           (12.5, 13, 13.5), (13.5, 14, 14.5), (14.5, 15, 15.5), (15.5, 16, 16.5),
                                                                           (16.5, 17, 17.5)])
# triplets_stats.to_csv('data/triplets_stats.csv', index=False)
# triplets_stats = pd.read_csv('data/triplets_stats.csv')
triplets_stats = triplets_stats.reset_index()

ot_validation_legend = {
    'I': ["red", "between interpolated and real"],
    'F': ["#4daf4a", "between first and real"],
    'L': ["#984ea3", "between last and real"],
    'R': ["#377eb8", "between random and real"]
}

plt.figure(figsize=(8, 8))
plt.title("Validation", fontsize=15)
plt.xlabel("Day", fontsize=15)
plt.ylabel("Wasserstein distance", fontsize=15)
plt.tick_params(labelsize=15)
legend = {}

for p, d in triplets_stats.groupby('name'):
    if p not in ot_validation_legend.keys():
        continue
    t = np.asarray(d['interval_mid'])
    m = np.asarray(d['mean'])
    s = np.asarray(d['std'])
    legend[p] = ot_validation_legend[p]
    plt.plot(t, m, '-o', color=ot_validation_legend[p][0])
    plt.fill_between(t, m - s, m + s, color=ot_validation_legend[p][0], alpha=0.2)
patch_list = [patches.Patch(color=c, label=l) for c, l in legend.values()]
plt.legend(handles=patch_list, fontsize=15)
