def plot_cumulative(self, mode='all', output_file: Union[None, str, Path] = None, **kwargs):
    """
    Generates a cumulative plot of the enrichment scores by position.

    Parameters
    -----------
    self : object from class *Screen*

    mode : str, default 'all'
        Options are 'all','SNV' and 'nonSNV'.

    output_file : str, default None
        If you want to export the generated graph, add the path and name of the file.
        Example: 'path/filename.png' or 'path/filename.svg'.

    **kwargs : other keyword arguments
        return_plot_object : boolean, default False
            If true, will return plotting objects (ie. fig, ax).

    Returns
    --------
    fig, ax : matplotlib figure and subplots
        Needs to have return_plot_object==True. By default they do
        not get returned.

    """

    # update kwargs
    temp_kwargs = copy.deepcopy(code_kwargs.kwargs())
    temp_kwargs.update(kwargs)
    temp_kwargs['figsize'] = kwargs.get('figsize', (3, 2))
    temp_kwargs['tick_spacing'] = kwargs.get('tick_spacing', 20)

    # import parameters
    code_kwargs._parameters()

    # create figure
    fig, ax = plt.subplots(figsize=temp_kwargs['figsize'])

    # Get data filtered
    df = _filter_bySNV(self, mode)
    cumsum = df.cumsum(skipna=False)['Score']
    plt.plot(df['Position'], cumsum / list(cumsum)[-1], color='red', lw=2)

    # y label
    y_label = 'Cumulative LoF'
    if list(cumsum)[-1] > 0:
        y_label = 'Cumulative GoF'

    # Graph limits
    plt.xlim(self.dataframe['Position'].min(), self.dataframe['Position'].max() + 1)
    plt.ylim(0, 1.1)

    ax.xaxis.set_major_locator(ticker.MultipleLocator(temp_kwargs['tick_spacing']))  # Plt ticks

    # Axis labels
    plt.title(temp_kwargs['title'], fontsize=12, fontname='Arial', color='k')
    plt.ylabel(y_label, fontsize=12, fontname="Arial", color='k', labelpad=5)
    plt.xlabel('Position', fontsize=12, fontname="Arial", color='k', labelpad=0)

    # x=y line
    plt.plot([0, df['Position'].max()], [0, 1], color='silver', lw=2, linestyle='--')

    # save file
    code_utils._save_work(fig, output_file, temp_kwargs)

    # return matplotlib object
    if temp_kwargs['return_plot_object']:
        return fig, ax

    # show plt figure
    if temp_kwargs['show']:
        plt.show()


def _filter_bySNV(self, mode):

    # Select all, SNV, nonSNV
    if mode == 'all':
        df = self.dataframe
    elif mode == 'SNV':
        df = self.dataframe_SNV
    elif mode == 'nonSNV':
        df = self.dataframe_nonSNV
    df = df.groupby(by='Position', as_index=False).mean()
    return df
