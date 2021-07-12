def plot_box(binned_x, y, output_file: Union[None, str, Path] = None, **kwargs):
    """
    Generates a boxplot. Data needs to be binned prior before using this
    function.

    Parameters
    -----------
    x, binned_y : arrays
        Contain the data is going to plot

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

    # Make figure
    fig, ax = plt.subplots(figsize=temp_kwargs['figsize'])

    # Plot data
    ax = sns.boxplot(binned_x, y, color='white', fliersize=2)

    plt.setp(ax.artists, edgecolor='k', facecolor='w')
    plt.setp(ax.lines, color='k')

    # graph parameters
    plt.title(temp_kwargs['title'], fontsize=10, fontname='Arial', color='k', pad=8)
    plt.ylabel(temp_kwargs['y_label'], fontsize=10, fontname="Arial", color='k', labelpad=0)
    plt.xlabel(temp_kwargs['x_label'], fontsize=10, fontname="Arial", color='k')

    # axes limits
    plt.xlim(temp_kwargs['xscale'])
    plt.ylim(temp_kwargs['yscale'])
    plt.grid()

    # save file
    code_utils._save_work(fig, output_file, temp_kwargs)

    # return matplotlib object
    if temp_kwargs['return_plot_object']:
        return fig, ax

    # show plt figure
    if temp_kwargs['show']:
        plt.show()
