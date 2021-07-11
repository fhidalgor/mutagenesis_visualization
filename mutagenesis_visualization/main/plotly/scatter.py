"""
This module contains the plotly scatter plot.
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

def plot_scatter_plotly(
    self,
    obj2,
    mode='pointmutant',
    show_results=False,
    output_html: Union[None, str, Path] = None,
    **kwargs: Dict[str, Any],
):
    '''
    Generate a scatter plot between object and a second object of the
    same class.

    Parameters
    ----------
    self : object from class *Screen*

    obj2 : object from class *Screen* to do the scatter with

    mode : str, default 'pointmutant'.
        Alternative set to "mean" for the mean of each position.

    show_results : boolean, default False
        If set to true, will export the details of the linear fit.

    output_html : str, default None
        If you want to export the generated graph into html, add
        the path and name of the file. Example: 'path/filename.html'.

    **kwargs : other keyword arguments

    Returns
    ----------
    fig : plotly object

    '''

    # update kwargs
    temp_kwargs = copy.deepcopy(default_kwargs())
    temp_kwargs.update(kwargs)
    temp_kwargs['figsize'] = kwargs.get('figsize', (4, 3))

    # Chose mode:
    if mode == 'pointmutant':
        df = code_utils._process_bypointmutant(self, obj2)
    elif mode == 'mean':
        df = code_utils._process_meanresidue(self, obj2)
        df['Variant'] = df['Position']
    # raise error if mode is not "mean" or "pointmutant"

    # Style
    pio.templates.default = "plotly_white"

    # create figure
    fig = px.scatter(x=df['dataset_1'], y=df['dataset_2'], trendline="ols", trendline_color_override="red")

    # Titles
    # hide text labels
    fig.update_traces(
        hovertext=df['Variant'], hovertemplate='Score_x: %{x}<br>Score: %{y}<br>Score_y: %{hovertext}<extra></extra>'
    )
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

    if show_results:
        px.get_trendline_results(fig).px_fit_results.iloc[0].summary()

    # save fig to html
    _save_html(fig, output_html)

    # return plotly object
    if temp_kwargs['return_plot_object']:
        return fig

    if temp_kwargs['show']:
        fig.show(config={'displayModeBar': False})
