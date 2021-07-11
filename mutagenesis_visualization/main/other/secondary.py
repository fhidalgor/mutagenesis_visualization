def plot_secondary(self, output_file: Union[None, str, Path] = None, **kwargs):
    """
    Generates a bar plot of data sorted by secondary elements (alpha helices
    and beta sheets).

    Parameters
    -----------
    self : object from class *Screen*

    output_file : str, default None
        If you want to export the generated graph, add the path and name of the file.
        Example: 'path/filename.png' or 'path/filename.svg'.

    **kwargs : other keyword arguments
        return_plot_object : boolean, default False
            If true, will return plotting objects (ie. fig, ax).

    Returns
    ----------
    fig, ax : matplotlib figure and subplots
        Needs to have return_plot_object==True. By default they do
        not get returned.

    """
    # Load parameters
    code_kwargs._parameters()

    # Update kwargs
    temp_kwargs = copy.deepcopy(code_kwargs.kwargs())
    temp_kwargs.update(kwargs)
    temp_kwargs['figsize'] = kwargs.get('figsize', (2.5, 2))
    temp_kwargs['yscale'] = kwargs.get('yscale', (-2, 1))

    # Get data
    df = _calculate_secondary(self.dataframe, self.secondary_dup)

    # Color
    df['Color'] = df.apply(
        code_utils._color_data,
        axis=1,
        args=(temp_kwargs['color_gof'], temp_kwargs['color_lof'])
    )

    # Make figure
    fig, ax = plt.subplots(figsize=temp_kwargs['figsize'])
    ticks = np.arange(0, len(df))  # label locations
    width = 0.5
    labels = df['Secondary']

    # Plot figure
    ax.bar(
        ticks,
        df['Score'],
        width,
        color=df['Color'],
        ec='k',
    )

    # graph parameters
    ax.set_xticks(ticks)
    ax.set_xticklabels(
        labels,
        fontsize=9,
        fontname="Arial",
        color='k',
        minor=False,
        rotation=0
    )
    ax.set_ylabel(
        r'$∆E^i_x$',
        fontsize=10,
        fontname="Arial",
        color='k',
        labelpad=12,
        rotation=0
    )
    ax.set_ylim(temp_kwargs['yscale'])
    plt.title(
        temp_kwargs['title'],
        horizontalalignment='center',
        fontname="Arial",
        fontsize=10,
        pad=5
    )

    # save file
    code_utils._save_work(fig, output_file, temp_kwargs)

    # return matplotlib object
    if temp_kwargs['return_plot_object']:
        return fig, ax

    if temp_kwargs['show']:
        plt.show()


def _calculate_secondary(df, secondary):
    '''
    Returns copy
    '''
    df = df.copy()
    df.insert(4, 'Secondary', secondary)
    df = df.groupby(['Secondary'], as_index=False, sort=False).mean()
    df = df[df['Secondary'].str.startswith(('β', 'α'))]
    df = df.drop(['Position'], axis=1)
    return df
