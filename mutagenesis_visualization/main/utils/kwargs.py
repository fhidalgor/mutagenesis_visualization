"""
This module contains standard input kwarg parameters.
"""
from typing import Dict, Any
from matplotlib.colors import LinearSegmentedColormap


def generate_default_kwargs() -> Dict[str, Any]:
    """
    Kwargs used in the methods and some other functions. Not all kwargs work on
    each method, read the individual description. Don't call this function
    on its own, use the parameters within the plotting methods.

    Example: mut.heatmap(colormap=colormap of interest)

    Parameters
    -----------
    colormap : cmap, default custom bluewhitered
        Used for heatmaps. You can use your own colormap or the ones provided by
        matplotlib. Example colormap = copy.copy((plt.cm.get_cmap('Blues_r')))

    colorbar_scale: tuple, default [-1, 1]
        Scale min and max used in heatmaps and correlation heatmaps.

    color: str, default 'k'
        Color used for the kernel plot line and the histogram.

    title : str, default 'Title'
        Title of plot.

    x_label : str, default 'x_label'
        Label of x axis.

    y_label : str, default 'y_label'
        Label of y axis.

    xscale: tuple, default (None, None)
        MinMax of x axis.

    yscale: tuple, default (None, None)
        MinMax of y axis.

    tick_spacing: int, default 1
        Space of axis ticks. Used for scatter and cumulative plots.

    outputfilepath : str, default ''
        Path where file will be exported to.

    outputfilename : str, default ''
        Name of the exported file.

    dpi : int, default 600
        Dots Per Inch in the created image.

    neworder_aminoacids: list, default list('DEKHRGNQASTPCVYMILFW*')
        Order of amino acids to display in heatmaps. Used for heatmaps.

    gof: int, default 1
        Cutoff of the enrichment score to classify a mutation as gain of function.
        Used on pymol and 3D methods.

    lof: int, default -1
        Cutoff of the enrichment score to classify a mutation as loss of funtion.
        Used on pymol and 3D methods.

    color_gof : str, default 'red'
        Color to color mutations above the gof cutoff.
        Used in pymol, 3D and mean methods.

    color_lof : str, default 'blue'
        Color to color mutations below the lof cutoff.
        Used in pymol, 3D and mean methods.

    cartoon_colors: list, default ['lightgreen', 'lavender', 'k']
        Colors used for secondary structure cartoon. Used for heatmap, mean and mean_count plots.

    text_labels: str, default 'None'
        Text labels that you can add to mean and mean_count plots. You will need to specify the coordinates.

    show: boolean, default True
        Whether to execute plt.show() or not on a matplotlib object.

    close: boolean, default False
        Whether to execute plt.close() or not on a matplotlib object.

    random_state : int, default 554
        Random state used for PCA function.

    bins : int or str, default 'auto'.
        Number of bins for the histogram. By default it will
        automatically decide the number of bins.

    return_plot_object : boolean, default False
        If true, will return plotting object.

    figsize_x : int

    figsize_y : int


    Returns
    --------
    default_kwargs : dict
        Dictionary with the default kwargs.
    """

    default_kwargs: Dict[str, Any] = {
        'figsize': None,
        'figsize_x': None,
        'figsize_y': None,
        'colormap': _generate_colormap(),
        'colorbar_scale': (-1, 1),
        'color': 'black',
        'kernel_color_replicates': None,
        'title': 'Title',
        'x_label': 'x_label',
        'y_label': 'y_label',
        'z_label': 'z_label',
        'xscale': (None, None),
        'yscale': (None, None),
        'tick_spacing': 1,
        'dpi': 600,
        'aminoacids': None,
        'neworder_aminoacids': None,
        'gof': 1,
        'lof': -1,
        'color_gof': 'red',
        'color_lof': 'blue',
        'cartoon_colors': ['lightgreen', 'lavender', 'k'],
        'text_labels': None,
        'show': True,
        'close': False,
        'random_state': 554,
        'bins': "auto",
        'return_plot_object': False,
        'grid': False,
        'metric': "",
        'title_fontsize': 12,
        'x_label_fontsize': 10,
        'y_label_fontsize': 10,
    }

    return default_kwargs


def _generate_colormap() -> LinearSegmentedColormap:
    """
    Generates standard blue, white, red color map.
    """
    cmap: LinearSegmentedColormap = LinearSegmentedColormap(
        'BlueWhiteRed',
        {
            'red': ((0.0, 0.0, 0.0), (0.15, 0.0, 0.0), (0.475, 1.0, 1), (0.525, 1.0, 1),
                    (0.85, 1.0, 1.0), (1.0, .8, 1)), 'green':
            ((0.0, 0.0, .0), (0.15, 0.5, 0.5), (0.475, 1.0, 1), (0.525, 1.0, 1), (0.85, 0.0, 0.0),
             (1.0, 0.0, 0.0)), 'blue': ((0.0, .5, .5), (0.15, 1, 1), (0.475, 1.0, 1),
                                        (0.525, 1.0, 1), (0.85, 0.0, 0.0), (1.0, 0.0, 0.0))
        },
    )
    return cmap
