"""
This module contains the plotly rank plot.
"""
import numpy as np
import pandas as pd
import copy
from os import path
import os
from pathlib import Path
from typing import Union
import plotly.express as px
import plotly.io as pio
import plotly.graph_objects as go
from Bio.PDB import PDBParser
import freesasa

from mutagenesis_visualization.main.default_parameters.kwargs import default_kwargs
from mutagenesis_visualization.main.plotly.plotly_utils import (
    save_html,
    color_3d_scatter,
    parse_pdb_coordinates,
    update_axes,
    update_layout,
    matplotlib_to_plotly,
)


def plot_rank_plotly(self, mode='pointmutant', outdf=False, output_html: Union[None, str, Path] = None, **kwargs):
    '''
    Generate a plotly rank plot so every mutation/residue is sorted based
    on enrichment score.

    Parameters
    ----------
    self : object from class *Screen*

    mode : str, default 'pointmutant'.
        Alternative set to "mean" for the mean of each position.

    outdf : boolean, default False
        If set to true, will return the df with the rank of mutations.

    output_html : str, default None
        If you want to export the generated graph into html, add the path and name of the file.
        Example: 'path/filename.html'.

    **kwargs : other keyword arguments

    Returns
    ----------
    df : Pandas dataframe
    fig : plotly object

    '''
    # update kwargs
    temp_kwargs = copy.deepcopy(default_kwargs())
    temp_kwargs.update(kwargs)
    temp_kwargs['figsize'] = kwargs.get('figsize', (4, 3))
    temp_kwargs['x_label'] = kwargs.get('x_label', 'Rank')
    temp_kwargs['y_label'] = kwargs.get('y_label', r'$∆E^i_x$')

    # Sort by enrichment scores
    df = self.dataframe.sort_values(by=['Score']).copy()

    # Chose mode:
    if mode == 'mean':
        df = df.groupby(by=['Position'], as_index=False).mean()
        df.sort_values(by=['Score'], inplace=True)
        df['Variant'] = df['Position']

    # Create figure
    fig = px.scatter(x=np.arange(len(df), 0, -1), y=df['Score'], text=df['Variant'])

    # Style
    pio.templates.default = "plotly_white"

    # Axes https://plotly.com/python/axes/
    fig.update_traces(mode="markers", hovertemplate='Position: %{x}<br>Score: %{y}<br>Variant: %{text}<extra></extra>')
    fig.update_xaxes(
        title_text=temp_kwargs['x_label'], showline=True, linewidth=2, linecolor='black', ticks="outside", mirror=True
    )
    fig.update_yaxes(
        title_text=temp_kwargs['y_label'], showline=True, linewidth=2, linecolor='black', ticks="outside", mirror=True
    )

    # Layout and title parameters https://plotly.com/python/figure-labels/
    fig.update_layout(
        width=temp_kwargs['figsize'][0] * 120,
        height=temp_kwargs['figsize'][1] * 120,
        font=dict(family="Arial, monospace", size=12, color="black"),
        title={'text': temp_kwargs['title'], 'xanchor': 'center', 'yanchor': 'top', 'x': 0.5}
    )

    # save fig to html
    _save_html(fig, output_html)

    # return plotly object
    if temp_kwargs['return_plot_object']:
        return fig

    if temp_kwargs['show']:
        fig.show(config={'displayModeBar': False})

    # return dataframe
    if outdf:
        return df


# save file html
def _save_html(fig, output_html):
    '''save to html'''
    if output_html:
        fig.write_html(str(Path(output_html)))
