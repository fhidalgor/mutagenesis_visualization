"""
This module contains the plotly histogram plot.
"""
import numpy as np
import pandas as pd
import copy
from os import path
import os
from pathlib import Path
from typing import Any, Dict, Union, List
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

def plot_histogram_plotly(
    self,
    mode='pointmutant',
    output_html: Union[None, str, Path] = None,
    **kwargs: Dict[str, Any],
):
    '''
    Generate a plotly histogram plot.

    Parameters
    ----------
    self : object from class *Screen*

    mode : str, default 'pointmutant'.
        Alternative set to "mean" for the mean of each position.

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
    temp_kwargs['figsize'] = kwargs.get('figsize', (4, 3))
    temp_kwargs['x_label'] = kwargs.get('x_label', 'Enrichment score')
    temp_kwargs['y_label'] = kwargs.get('y_label', 'Probability density')
    temp_kwargs['title'] = kwargs.get('title', 'Histogram')

    # Copy dataframe
    df = self.dataframe.copy()

    #fig = px.histogram(data_frame=df, x='Score', histnorm='probability density')
    # px.histogram doesnt plot 2 traces

    # Create figure
    fig = go.Figure()
    fig.add_trace(
        go.Histogram(
            x=df['Score'],
            histnorm='probability density',
            name='Population',
            marker_color='black',
        )
    )

    # Create second histogram
    if mode != 'pointmutant':
        df = _select_grouping(self, mode)
        # Add second trace
        fig.add_trace(
            go.Histogram(
                x=df['Score'],
                histnorm='probability density',
                name=mode.capitalize(),
                marker_color='green',
            ),
        )
        # Overlay both histograms
        fig.update_layout(
            barmode='overlay',
            legend=dict(orientation="h", y=(1 + (1.5 / temp_kwargs['figsize'][1]**2)), bgcolor='rgba(0,0,0,0)'),
        )
        # Reduce opacity to see both histograms
        fig.update_traces(opacity=0.75)
        # Title
        temp_kwargs['title'] = 'Histogram filtered by: {}'.format(mode)

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
        automargin=True,
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
