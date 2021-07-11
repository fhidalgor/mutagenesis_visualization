"""
This module contains the plotly mean enrichment plot.
"""
import copy
from os import path
import os
from pathlib import Path
from typing import Any, Dict, Union, List
import numpy as np
import pandas as pd
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

def plot_mean_plotly(
    self,
    mode='mean',
    output_html: Union[None, str, Path] = None,
    **kwargs: Dict[str, Any],
):
    '''
    Generate a plotly mean plot.

    Parameters
    ----------
    self : object from class *Screen*

    mode : str, default 'mean'
        Specify what enrichment scores to show. If mode = 'mean', it will show the mean of
        each position. If mode = 'A', it will show the alanine substitution profile. Can be
        used for each amino acid. Use the one-letter code and upper case.

    output_html : str, default None
        If you want to export the generated graph into html, add the path and name of the file.
        Example: 'path/filename.html'.

    **kwargs : other keyword arguments

    Returns
    ----------
    fig : plotly object

    '''
    # update kwargs
    temp_kwargs = copy.deepcopy(default_kwargs())
    temp_kwargs.update(kwargs)
    temp_kwargs['figsize'] = kwargs.get('figsize', (4.5, 3))
    temp_kwargs['x_label'] = kwargs.get('x_label', 'Position')
    temp_kwargs['y_label'] = kwargs.get('y_label', 'Enrichment score')
    temp_kwargs['title'] = kwargs.get('title', 'Filtered by: {}'.format(mode))
    temp_kwargs['color_gof'] = kwargs.get('color_gof', '#FD3216')
    temp_kwargs['color_lof'] = kwargs.get('color_lof', '#6A76FC')

    # Chose mode:
    df = _select_grouping(self, mode)

    # Calculate colors
    df['Color'] = df.apply(code_utils._color_data, axis=1, args=(temp_kwargs['color_gof'], temp_kwargs['color_lof']))

    # Create figure
    #fig = px.bar(data_frame=df, x='Position', y='Score', color='Color')
    #px.bar was switching colors when the first value of Score was negative

    fig = go.Figure(data=[go.Bar(
        x=df['Position'],
        y=df['Score'],
        marker_color=df['Color'],
        marker_line_width=0,
    )])

    # Style
    pio.templates.default = "plotly_white"

    # UPDATE AXES
    fig.update_xaxes(
        title_text=temp_kwargs['x_label'],
        showline=True,
        linewidth=2,
        linecolor='black',
        ticks="outside",
        mirror=True,
    )
    fig.update_yaxes(
        title_text=temp_kwargs['y_label'],
        showline=True,
        linewidth=2,
        linecolor='black',
        ticks="outside",
        mirror=True,
    )

    # Color and width of bars
    #fig.update_traces(marker_line_width=0, )

    # Layout and title parameters https://plotly.com/python/figure-labels/
    fig.update_layout(
        width=temp_kwargs['figsize'][0] * 120,
        height=temp_kwargs['figsize'][1] * 120,
        showlegend=False,
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
