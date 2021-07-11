"""
This module contains standard input default rcParams parameters.
"""
from matplotlib import rcParams

def parameters() -> None:
    """
    Default rcParams.
    """
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

def font_parameters() -> None:
    """
    Default math font rcParams.
    """
    # math font
    rcParams['mathtext.fontset'] = 'custom'
    rcParams['mathtext.rm'] = 'Arial'
    rcParams['svg.fonttype'] = 'none'
