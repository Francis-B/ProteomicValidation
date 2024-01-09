# Create plot for each conditions
for cond in CONDITIONS:
    # Bool arrays to select unique pg to condition
    is_unique_cond = np.array([pg in pg_unique[cond] for pg in pg_identified[cond]])
    is_unique_fomonet = np.array([pg in all_unique_fomonet for pg in pg_identified[cond]])

    # Extract log2 transformed mean quantifications 
    log2_quants_unique_cond = log2_quants[cond][is_unique_cond]
    log2_quants_unique_fomonet = log2_quants[cond][is_unique_cond*is_unique_fomonet]

    # Determine bin size
    min_bin = int(np.min(log2_quants[cond]))
    max_bin = int(np.max(log2_quants[cond])) + 1
    bins = np.arange(min_bin, max_bin, 0.2)

    # Initialize figure
    fig, ax = plt.subplots(3, 1, figsize=(10, 7),
                            gridspec_kw={'height_ratios': [1, .2, 0.075]})
    
    # Quant distribution for all PG identified
    x, _, _ = ax[0].hist(log2_quants[cond], bins=bins, color = '#9E8399', edgecolor='#DFC8DB')
    ax[0].text(max_bin-2, x.max()*.9, f'PG = {log2_quants[cond].shape[0]}')
    ax[0].set_xticks([])
    
    # Quant distribution for PG unique to condition
    x, _, _ = ax[1].hist(log2_quants_unique_cond, bins=bins, color = '#9E8399', edgecolor='#DFC8DB')
    ax[1].text(max_bin-2, x.max()*.9, f'PG = {log2_quants_unique_cond.shape[0]}')
    ax[1].set_xticks([])
    
    # Quant distribution for PG unique to condition and fomonet
    x, y, _ = ax[2].hist(log2_quants_unique_fomonet, bins=bins, color = '#9E8399', edgecolor='#DFC8DB')
    ax[2].set_yticks(np.arange(0, x.max()+1, 2))
    ax[2].text(max_bin-2, x.max()*.9, f'PG = {log2_quants_unique_fomonet.shape[0]}')
    
    # Add axis labels and title
    fig.text(0.5, 0.95, f'PG quantification distribution for {cond} condition', ha='center', fontdict={'fontsize': 14})
    fig.text(0.05, 0.5,'Number of protein groups (PG)', rotation=90, va='center', fontdict={'fontsize': 14})
    fig.text(0.5, 0.03, 'Log2 transformed mean quantification', ha='center', fontdict={'fontsize': 14})
    
    for a, t in zip(ax, ['A', 'B', 'C']):
        a.set_title(t, x=1)
        
    # background axes object for plotting the vertical line
    gs = fig.add_gridspec(3, 1, hspace=0)
    ax =  fig.add_subplot(gs[:, :], sharex = ax[0])
    # set background color to transparent and turn off the frame
    ax.patch.set_alpha(0)
    ax.axis("off")
    # plot the vertical line
    ax.axvline(np.mean(log2_quants[cond]), color='#2A5E83', linestyle='--')
    
    # fig.savefig(DIST_PLOT(cond), format='svg')