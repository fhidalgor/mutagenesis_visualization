#!/usr/bin/env python
# coding: utf-8

# # Import Modules

# In[1]:


from matplotlib.colors import LinearSegmentedColormap
from matplotlib import rcParams


# # Kwargs

# In[2]:


def kwargs():
    """
    Kwargs used in the methods and some other functions. Not all kwargs work on
    each method, read the individual description.

    Parameters
    -----------
    colormap : cmap, default custom bluewhitered
        Used for heatmaps. You can use your own colormap or the ones provided by
        matplotlib. Example colormap = copy.copy((plt.cm.get_cmap('Blues_r')))

    colorbar_scale: list, default [-1, 1]
        Scale min and max used in heatmaps and correlation heatmaps.

    color: str, default 'k'
        Color used for the kernel plot line.

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

    random_state : int, default 554
        Random state used for PCA function.

    bins : int, default 50
        Number of bins used for kernel and histograms.
    
    return_plot_object : boolean, default False
        If true, will return plotting object.
    
    Returns
    --------
    default_kwargs : dict
        Dictionary with the default kwargs.

    """
    default_kwargs = {
        'colormap': _generatecolormap(),
        'colorbar_scale': [-1, 1],
        'color': 'k',
        'title': 'Title',
        'x_label': 'x_label',
        'y_label': 'y_label',
        'z_label': 'y_label',
        'xscale': (None, None),
        'yscale': (None, None),
        'tick_spacing': 1,
        'dpi': 600,
        'aminoacids': list('ACDEFGHIKLMNPQRSTVWY*'),
        'neworder_aminoacids': list('DEKHRGNQASTPCVYMILFW*'),
        'gof': 1,
        'lof': -1,
        'color_gof': 'red',
        'color_lof': 'blue',
        'cartoon_colors': ['lightgreen', 'lavender', 'k'],
        'text_labels': None,
        'show': True,
        'random_state': 554,
        'bins': 50,
        'return_plot_object':False,
    }

    return default_kwargs


def _generatecolormap():
    cmap = LinearSegmentedColormap(
        'BlueWhiteRed',
        {
            'red': ((0.0, 0.0, 0.0), (0.15, 0.0, 0.0), (0.475, 1.0, 1),
                    (0.525, 1.0, 1), (0.85, 1.0, 1.0), (1.0, .8, 1)), 'green':
            ((0.0, 0.0, .0), (0.15, 0.5, 0.5), (0.475, 1.0, 1), (0.525, 1.0, 1),
             (0.85, 0.0, 0.0), (1.0, 0.0, 0.0)), 'blue':
            ((0.0, .5, .5), (0.15, 1, 1), (0.475, 1.0, 1), (0.525, 1.0, 1),
             (0.85, 0.0, 0.0), (1.0, 0.0, 0.0))
        },
    )
    return cmap


# In[4]:


def _parameters():
    # normal font
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['Arial']

    # math font
    rcParams['mathtext.fontset'] = 'custom'
    rcParams['mathtext.rm'] = 'Arial'
    rcParams['svg.fonttype'] = 'none'

    # add grid
    rcParams['grid.color'] = 'silver'
    rcParams['grid.linestyle'] = '--'
    rcParams['grid.linewidth'] = 1
    rcParams['lines.dashed_pattern'] = [5, 10]
    rcParams['axes.axisbelow'] = True
    # Parameters for all graphs
    rcParams['xtick.labelsize'] = 9
    rcParams['ytick.labelsize'] = 9
    return


def _font_parameters():
    # math font
    rcParams['mathtext.fontset'] = 'custom'
    rcParams['mathtext.rm'] = 'Arial'
    rcParams['svg.fonttype'] = 'none'
    return

