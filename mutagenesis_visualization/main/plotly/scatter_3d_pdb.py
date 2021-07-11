"""
This module contains the plotly 3D scatter plot where you can add PDB
properties.
"""
import copy
from pathlib import Path
from typing import Any, Dict, Union, List
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.io as pio
import plotly.graph_objects as go
from Bio.PDB import PDBParser

from mutagenesis_visualization.main.default_parameters.kwargs import default_kwargs
from mutagenesis_visualization.main.plotly.plotly_utils import (
    save_html,
    color_3d_scatter,
    parse_pdb_coordinates,
    update_axes,
    update_layout,
    matplotlib_to_plotly,
)


def plot_scatter_3D_pdbprop_plotly(
    self,
    plot: List[str] = ['Distance', 'SASA', 'log B-factor'],
    mode: str = 'mean',
    pdb_path=None,
    custom=None,
    position_correction: int = 0,
    chain: str = 'A',
    output_df: bool = False,
    output_html: Union[None, str, Path] = None,
    **kwargs: Dict[str, Any],
):
    """
    Generates a 3-D scatter plot of different properties obtained from
    the PDB. PDBs may have atoms missing, you should fix the PDB before
    using this method. We recommend you use matplotlib for interactive plot.

    Parameters
    -----------
    self : object from class *Screen*
        **kwargs : other keyword arguments.

    plot : list, default ['Distance', 'SASA', 'log B-factor']
        List of 3 elements to plot. Other options are 'Score' and Custom.
        If custom, add the label to the third element of the list ie.
        ['Distance', 'SASA', 'Conservation'].

    mode : str, default 'mean'
        Specify what enrichment scores to use. If mode = 'mean', it will
        use the mean of each position to classify the residues.
        If mode = 'A', it will use the Alanine substitution profile.
        Can be used for each amino acid. Use the one-letter code and
        upper case.

    pdb_path : str, default None
        User should specify the path PDB.

    custom : list or dataframe or np.array, default None
        If you want to add a custom dataset to plot, use custom. On the
        parameter plot, the 3rd item of the list will be the label for
        your custom dataset.

    df_color : pandas dataframe, default None
        The color of each residue can also be included. You must label
        that label column.

    color_by_score : boolean, default True
        If set to False, the points in the scatter will not be colored
        based on the enrichment score.

    position_correction : int, default 0
        If the pdb structure has a different numbering of positions than
        you dataset, you can correct for that. If your start_position = 2,
        but in the PDB that same residue is at position 20,
        position_correction needs to be set at 18.

    chain : str, default 'A'
        Chain of the PDB file to get the coordinates and SASA from.

    output_df : boolean, default False
        If true, this method will return the dataframe with the data.
        Set return_plot_object for this to work.

    output_html : str, default None
        If you want to export the generated graph into html, add the
        path and name of the file.
        Example: 'path/filename.html'.

    **kwargs : other keyword arguments

    Returns
    ---------
    fig : plotly object

    df_items : pandas dataframe
        Contains the plotted data. Needs to have output_df set to true.

    """
    # Update kwargs
    temp_kwargs = copy.deepcopy(default_kwargs())
    temp_kwargs.update(kwargs)
    temp_kwargs['x_label'] = kwargs.get('x_label', plot[0])
    temp_kwargs['y_label'] = kwargs.get('y_label', plot[1])
    temp_kwargs['z_label'] = kwargs.get('z_label', plot[2])

    # Get Scores and colors
    df_scores = color_3d_scatter(self.dataframe, mode, temp_kwargs['lof'], temp_kwargs['gof'])

    # If coordinates is not an input, get it from the pdb
    df_items = parse_pdb_coordinates(self, pdb_path, position_correction, chain, sasa=True)

    # Add scores
    df_items['Score'] = list(df_scores['Score'])

    # Custom data
    if custom is not None:
        df_items[plot[2]] = custom

        # Plot figure
    fig = px.scatter_3d(
        df_items,
        x=plot[0],
        y=plot[1],
        z=plot[2],
        color='Score',
        color_continuous_scale=matplotlib_to_plotly(temp_kwargs['colormap'], ),
        range_color=temp_kwargs['colorbar_scale'],
    )

    # update axes
    fig = update_axes(fig, temp_kwargs)

    # for the clickable part
    fig.update_traces(hovertext=df_items['Position'], hovertemplate='Position: %{hovertext}')
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

    if output_df:
        return df_items, df_scores
