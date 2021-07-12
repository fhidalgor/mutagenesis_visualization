def plot_neighboreffect(self, offset=1, output_file: Union[None, str, Path] = None, **kwargs):
    """
   DEPRECATED.

   Generate a miniheatmap plot telling you the effect of having a residue in front or behind.
   It corrects for the effect of that amino acid on the rest of the population.

   Parameters
   ----------
   self : object from class *Screen*

   offset : int, default 1
      if you want to study effects of a residue when is behind or in front of another residue.
      offset of 1 means that you evaluate the effect of following residue n+1 on n. On a "MTEY..." sequence,
      you would look at the effect of T on M, E on T, Y on E, etc.. and then group by residue (n+1).
      offset of -1 means that you look at the previous residue (n-1 on n).

   output_file : str, default None
      If you want to export the generated graph, add the path and name of the file.
      Example: 'path/filename.png' or 'path/filename.svg'.

   **kwargs : other keyword arguments
       return_plot_object : boolean, default False
           If true, will return plotting objects (ie. fig, ax).

   Returns
   ----------
   fig, ax, cb : matplotlib figure and subplots
       Needs to have return_plot_object==True. By default they do
       not get returned.

   """
    # load font parameters
    code_kwargs._font_parameters()

    # update kwargs
    temp_kwargs = copy.deepcopy(code_kwargs.kwargs())
    temp_kwargs.update(kwargs)
    if '*' in temp_kwargs['neworder_aminoacids']:
        temp_kwargs['neworder_aminoacids'].remove('*')

    # do offset, no stop codons
    df = _normalize_neighboreffect(self, offset, temp_kwargs['neworder_aminoacids'])

    # Plot
    fig, ax, cb = _plot_miniheatmap(df, output_file, temp_kwargs)

    # return matplotlib object
    if temp_kwargs['return_plot_object']:
        return fig, ax, cb

    # show figure
    if temp_kwargs['show']:
        plt.show()


def _plot_miniheatmap(df, output_file, temp_kwargs):
    """
   Aux plot that will do the heavy lifting to plot a miniheatmap.

   Parameters
   ------------
   df: pandas dataframe
       dataframe with the data to plot.

   output_file : str, default None
       If you want to export the generated graph, add the path and name of the file.
       Example: 'path/filename.png' or 'path/filename.svg'.

   temp_kwargs : kwargs

   Returns
   --------
   fig, ax, cb : matplotlib figure and subplots

   """
    # declare figure and subplots
    coeff = len(df.columns) / 19 * 1.05
    fig = plt.figure(figsize=(2.5 * coeff, 2.5))
    gs = gridspec.GridSpec(nrows=1, ncols=1)
    ax = plt.subplot(gs[0])

    # main heatmap
    heatmap = ax.pcolor(
        df.to_numpy(),
        vmin=temp_kwargs['colorbar_scale'][0],
        vmax=temp_kwargs['colorbar_scale'][1],
        cmap=temp_kwargs['colormap'],
        edgecolors='k',
        linewidths=0.2,
        color='darkgrey'
    )

    # ____________axes manipulation____________________________________________
    # put the major ticks at the middle of each cell
    ax.set_xticks(np.arange(df.shape[1]) + 0.5, minor=False)
    ax.set_yticks(np.arange(df.shape[0]) + 0.5, minor=False)

    # position of axis labels
    ax.tick_params('x', direction='out', pad=-2.5)
    ax.tick_params('y', direction='out', pad=0.4)

    # want a more natural, table-like display
    ax.invert_yaxis()
    ax.xaxis.tick_top()

    # remove ticks
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')

    # so labels of x and y do not show up and my labels show up instead
    ax.set_xticklabels(list(df.columns), fontsize=6.5, fontname="Arial", color='k', minor=False)
    ax.set_yticklabels(temp_kwargs['neworder_aminoacids'], fontsize=6.5, fontname="Arial", color='k', minor=False)

    # align the labels of the y axis
    for ylabel in ax.get_yticklabels():
        ylabel.set_horizontalalignment('center')

    # _____________________________________________________________________________

    # for color bar format
    cb = plt.colorbar(
        heatmap,
        fraction=0.025,
        pad=0.05,
        aspect=5,
        ticks=[
            temp_kwargs['colorbar_scale'][0],
            np.mean(temp_kwargs['colorbar_scale']), temp_kwargs['colorbar_scale'][1]
        ],
        orientation='vertical'
    )
    cb.ax.set_yticklabels(cb.ax.get_yticklabels(), fontsize=7, fontname="Arial", color='k')
    cb.update_ticks()
    plt.text(
        len(df.columns) + 2,
        7.8,
        r'$\langle∆E^x_i\rangle_x$',
        horizontalalignment='center',
        fontsize=7,
        fontname="Arial",
        color='k'
    )

    # for putting title on graph
    plt.title(temp_kwargs['title'], horizontalalignment='center', fontname="Arial", fontsize=10, pad=10)
    plt.ylabel('Amino Acid Substitution', fontsize=10, labelpad=-1)

    # save file
    code_utils._save_work(fig, output_file, temp_kwargs)

    # return matplotlib object
    return fig, ax, cb


def _normalize_neighboreffect(self, offset, neworder):
    '''
   For every residue, subtract the average effect of a substitution
   Returns a normalized dataframe
   '''
    aalist = list('ACDEFGHIKLMNPQRSTVWY')
    # Add offset sequence to df
    df = _transform_dataset_offset(self, offset, False)

    # calculate mean effect using condensed heatmap
    mean = _condense_heatmap(self.dataframe, aalist)

    df_normalized = pd.DataFrame()
    for aa in aalist:
        # Choose the neighbors of an aa
        aa_neighbors = df.loc[df['Sequence'] == aa]
        # Do the mean substitution of amino acids that are repeated
        aa_neighbors = aa_neighbors.groupby(['Sequence_old', 'Aminoacid'], as_index=False).mean()
        # Make into table
        aa_neighbors_pivoted = aa_neighbors.pivot_table(values='Score', index='Aminoacid', columns='Sequence_old')
        aa_neighbors_pivoted.reset_index(drop=True, inplace=True)
        # Get the mean of the amino acids that appear in the aa_neighbors subset
        mean_neighbors = mean[list(aa_neighbors_pivoted.columns)]
        # Subtract average effect and do mean
        df_normalized[aa] = (aa_neighbors_pivoted - mean_neighbors).mean(axis=1)

    # Sort by aa
    df_normalized = df_normalized[neworder]
    # Sort in y axis desired order
    df_normalized = _sort_yaxis_aminoacids(df_normalized, neworder, aalist)
    return df_normalized


def _sort_yaxis_aminoacids(df, neworder, oldorder=list('ACDEFGHIKLMNPQRSTVWY')):
    # Sort in y axis desired order
    df['Aminoacid_new'] = oldorder
    df['Aminoacid_new'] = pd.Categorical(df['Aminoacid_new'], neworder)
    df.sort_values(by=['Aminoacid_new'], inplace=True)
    df.drop(['Aminoacid_new'], inplace=True, axis=1)

    return df
