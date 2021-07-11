"""
This module contains the plotly 3D scatter plot.
"""
import copy
from pathlib import Path
from typing import Any, Dict, Union, List
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.io as pio
import plotly.graph_objects as go

from mutagenesis_visualization.main.default_parameters.kwargs import default_kwargs
from mutagenesis_visualization.main.plotly.plotly_utils import (
    save_html,
    color_3d_scatter,
    parse_pdb_coordinates,
    update_axes,
    update_layout,
    matplotlib_to_plotly,
)


def plot_scatter_3D_plotly(
    self,
    mode='mean',
    pdb_path=None,
    df_coordinates=None,
    position_correction=0,
    chain='A',
    squared=False,
    output_html: Union[None, str, Path] = None,
    **kwargs: Dict[str, Any],
):
    """
    Generates a 3-D scatter plot of the x,y,z coordinates of the C-alpha atoms
    of the residues, color coded by the enrichment scores. PDBs may have atoms
    missing, you should fix the PDB before using this method. Use matplotlib
    for interactive plot.

    Parameters
    -----------
    self : object from class *Screen*
        **kwargs : other keyword arguments.

    mode : str, default 'mean'
        Specify what enrichment scores to use. If mode = 'mean', it will use the mean of
        each position to classify the residues. If mode = 'A', it will use the Alanine substitution profile.
        Can be used for each amino acid. Use the one-letter code and upper case.

    pdb : str, default None
        User should specify the path PDB chain.

    df_coordinates: pandas dataframe, default None
        If no pdb is included, the user must pass the 3-D coordinates of the residues to plot.
        In here you have more flexibility and you can select other atoms besides the C-alpha.

    position_correction : int, default 0
        If the pdb structure has a different numbering of positions than you dataset,
        you can correct for that. If your start_position = 2, but in the PDB that same residue
        is at position 20, position_correction needs to be set at 18.

    chain : str, default 'A'
        Chain of the PDB file to get the coordinates and SASA from.

    squared : boolean, False
        If this parameter is True, the algorithm will center the data, and plot the square value of the
        distance.

    output_html : str, default None
        If you want to export the generated graph into html, add the path and name of the file.
        Example: 'path/filename.html'.

    **kwargs : other keyword arguments

    Returns
    ---------
    fig : plotly object

    """

    # update kwargs
    temp_kwargs = copy.deepcopy(default_kwargs())
    temp_kwargs.update(kwargs)
    temp_kwargs['x_label'] = kwargs.get('x_label', 'x (Å)')
    temp_kwargs['y_label'] = kwargs.get('y_label', 'y (Å)')
    temp_kwargs['z_label'] = kwargs.get('z_label', 'z (Å)')

    # Get Scores and colors
    df = color_3d_scatter(self.dataframe, mode, temp_kwargs['lof'], temp_kwargs['gof'])

    # If coordinates is not an input, get it from the pdb
    if df_coordinates is None:
        df_coordinates = parse_pdb_coordinates(self, pdb_path, position_correction, chain)

    # Plot figure
    x, y, z = ('x', 'y', 'z')
    if squared is True:
        x, y, z = ('x_cent', 'y_cent', 'z_cent')
    fig = px.scatter_3d(
        df_coordinates,
        x=x,
        y=y,
        z=z,
        color=df['Score'],
        color_continuous_scale=matplotlib_to_plotly(temp_kwargs['colormap'], ),
        range_color=temp_kwargs['colorbar_scale']
    )

    # update axes
    fig = update_axes(fig, temp_kwargs)

    # for the clickable part
    fig.update_traces(hovertext=df['Position'], hovertemplate='Position: %{hovertext}')
    # title
    fig = update_layout(fig, temp_kwargs)

    # save fig to html
    save_html(fig, output_html)

    # return plotly object
    if temp_kwargs['return_plot_object']:
        return fig

    # show only if asked
    if temp_kwargs['show']:
        fig.show(config={'displayModeBar': False})
