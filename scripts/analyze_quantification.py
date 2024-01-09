""" Analyze normalized quantification and create plots."""

from collections import defaultdict
import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import upsetplot

import utils

# Infiles and params
QUANT = snakemake.input.quant
UNIQUES = snakemake.input.unique_ids
CONDITIONS = snakemake.params.conditions
ALLOWED_NA = snakemake.params.allowed_na

# Outfiles
UPSETPLOT = snakemake.output.upsetplot
DIST_PLOT = snakemake.output.quant_plot
FOMONET_ID = snakemake.output.fomonet_id
QUANT_PG = snakemake.output.quant_pg

# Load list of FOMOnet predictions ids
with open(UNIQUES, 'r', encoding='utf-8') as f:
    unique_ids = [line.strip() for line in f.readlines()]

# MSstats normalized quant
pg_dict = utils.extract_msstats_quant(QUANT, CONDITIONS, ALLOWED_NA)

is_unique_fomonet = utils.get_fomonet_pg(pg_dict, unique_ids)
is_unique_cond = utils.get_condition_unique_pg(pg_dict)
mpl.rcParams['axes.spines.right'] = False
mpl.rcParams['axes.spines.top'] = False


# ------------------------------------------------------------------------------
# Make quantification distribution plot
fig, ax = plt.subplots(4, 3, figsize=(10, 7),
                            gridspec_kw={'height_ratios': [1, .5, 0.2, 0.15]},
                            sharex=True, sharey='row')

for count, cond in enumerate(CONDITIONS):
    # Determine bin size
    log2_quants = pg_dict[cond]['quants']
    min_bin = int(np.min(log2_quants))
    max_bin = int(np.max(log2_quants)) + 1
    bins = np.arange(min_bin, max_bin, 0.5)

    # Quant distribution for all PG identified
    data = log2_quants
    ax[0][count].hist(data, 
                      bins=bins,
                      color = '#9E8399',
                      edgecolor='#DFC8DB',
                      linewidth=.3)
    ax[0, count].text(x=np.mean(ax[0][0].get_xlim())*1.15, 
                      y=ax[0][0].get_ylim()[1]*.8,
                      s=f'PG = {log2_quants.shape[0]}')

    # Quant distribution for PG unique to condition
    data = log2_quants[is_unique_cond[cond]]
    ax[1][count].hist(data, 
                      bins=bins,
                      color = '#9E8399', 
                      edgecolor='#DFC8DB',
                      linewidth=.3)
    ax[1][count].text(x=np.mean(ax[1][0].get_xlim())*1.15,
                      y=ax[1][0].get_ylim()[1]*.8,
                      s=f'PG = {data.shape[0]}')

    # Quant distribution for PG unique to fomonet but not to condition
    data = log2_quants[is_unique_fomonet[cond]]
    ax[2][count].hist(data,
                      bins=bins,
                      color = '#4A6F8F',
                      edgecolor='#9CAEBE',
                      linewidth=.3)
    ax[2][count].text(x=np.mean(ax[2][0].get_xlim())*1.15,
                      y=ax[2][0].get_ylim()[1]*.8,
                      s=f'PG = {data.shape[0]}')

    # Quant distribution for PG unique to condition and fomonet
    data = log2_quants[is_unique_cond[cond]*is_unique_fomonet[cond]]
    ax[3][count].hist(data,
                      bins=bins, 
                      color='#4A6F8F', 
                      edgecolor='#9CAEBE',
                      linewidth=.3)
    ax[3][count].text(x=np.mean(ax[3][0].get_xlim())*1.15,
                      y=ax[3][0].get_ylim()[1]*.8,
                      s=f'PG = {data.shape[0]}')

    # background axes object for plotting the vertical line
    gs = fig.add_gridspec(3, 3, hspace=0)
    ax_com =  fig.add_subplot(gs[:, count], sharex = ax[0][count])
    # set background color to transparent and turn off the frame
    ax_com.patch.set_alpha(0)
    ax_com.axis("off")
    # plot the vertical line
    ax_com.text(np.mean(ax[0][0].get_xlim()), 1.05, f'{cond}')
    ax_com.axvline(np.mean(log2_quants), color='#5B5D5F', linestyle='--')


# Set y ticks
ax[1][0].set_yticks(np.arange(0, ax[1][0].get_ylim()[1], 50))
ax[2][0].set_yticks(np.arange(0, ax[2][0].get_ylim()[1], 5))
ax[3][0].set_yticks(np.arange(0, 7, 3))

# Row labels
last_ax = len(CONDITIONS)-1
for c, text in enumerate(['All PG', 'All PG unique \nto condition',
                          'All FOMOnet \nPG', 'FOMOnet PG \nunique to \ncondition']): 
    ax[c][last_ax].text(38, np.mean(ax[c][last_ax].get_ylim()), text, ha='center')

# Add axis labels and title
fig.text(0.5, 0.98, f'PG quantification distribution ({ALLOWED_NA} NA allowed)',
         ha='center', fontdict={'fontsize': 14})
fig.text(0.05, 0.5,'Number of PG', rotation=90, va='center', fontdict={'fontsize': 14})
fig.text(0.5, 0.03, 'Log2 transformed mean quantification', ha='center', fontdict={'fontsize': 14})

fig.savefig(DIST_PLOT, format='svg')


# ------------------------------------------------------------------------------
# Make upsetplot
extract_pg = lambda cond: pg_dict[cond]['names'][is_unique_fomonet[cond]]
# Plot the results        
fig = plt.figure(figsize=(10, 10))
upsetplot.plot(upsetplot.from_contents({'Contrôle':extract_pg('ctr'),
                                        'Ischémie':extract_pg('hyp'), 
                                        'Reperfusion':extract_pg('rep')}),
               fig=fig)
fig.savefig(UPSETPLOT, format='svg')

# ------------------------------------------------------------------------------
# Get identified predictions informations

# Make lists of all identified FOMOnet predictions
id_request = []
quants_dict = defaultdict(dict)
prots = []
for cond in CONDITIONS:
    # Make clean quantification dict with FOMOnet predictions
    fomonet_pg = pg_dict[cond]['names'][is_unique_fomonet[cond]]
    fomonet_quants = pg_dict[cond]['quants'][is_unique_fomonet[cond]]
    quants_dict[cond] = {pg: quant
                        for pg, quant in zip(fomonet_pg, fomonet_quants)}

    # Id list for Biomart query
    id_request.extend([prot for pg in fomonet_pg for prot in pg.split(';')])

# Write id list to pickle to make query in another rule
with open(FOMONET_ID, 'wb') as f:
    pickle.dump(id_request, f)

# Write quantification dict to pickle
with open(QUANT_PG, 'wb') as f:
    pickle.dump(quants_dict, f)
