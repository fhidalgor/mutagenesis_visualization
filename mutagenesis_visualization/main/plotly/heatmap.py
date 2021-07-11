"""
This module contains the plotly heatmap plot.
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

def plot_heatmap_plotly(
    self,
    output_html: Union[None, str, Path] = None,
    **kwargs: Dict[str, Any],
):
    '''
    Generate a plotly histogram plot.

    Parameters
    ----------
    self : object from class *Screen*

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
    temp_kwargs['figsize'] = kwargs.get('figsize', (8, 3))
    temp_kwargs['x_label'] = kwargs.get('x_label', '')
    temp_kwargs['y_label'] = kwargs.get('y_label', '')

    # sort data by rows in specified order by user
    df = code_utils._df_rearrange(
        code_utils._add_SNV_boolean(self.dataframe_stopcodons.copy()),
        temp_kwargs['neworder_aminoacids'],
        values='Score_NaN',
        show_snv=False
    )

    # get labels for texthover and reindex
    text_hover = self.dataframe_stopcodons.pivot(values='Variant', index='Aminoacid', columns='Position')
    text_hover = text_hover.reindex(list(df.index))

    # Create figure
    fig = go.Figure(data=go.Heatmap(
        z = np.around(df.to_numpy(),2),
        x = list(df.columns),
        y = list(df.index),
        zmin=temp_kwargs['colorbar_scale'][0],
        zmax=temp_kwargs['colorbar_scale'][1],
        colorscale=_matplotlib_to_plotly(temp_kwargs['colormap']),
        text = text_hover,
        colorbar=dict( # modify colorbar properties
            len=0.8,
            thickness=10,
            outlinewidth=2,
            outlinecolor='rgb(0,0,0)',
            showticklabels=True,
            xpad=0,
        ),
    ))

    fig.update_traces(hovertemplate='Aminoacid substitution: %{text}<br>Enrichment score: %{z}<extra></extra>')

    # Style
    pio.templates.default = "plotly_white"

    # UPDATE AXES
    fig.update_xaxes(
        title_text=temp_kwargs['x_label'],
        showline=True,
        linewidth=2,
        linecolor='black',
        ticks=None,
        mirror=True,
        side="top",
        dtick=1,  #frequency of ticks
        tickangle=0,
        tickvals=list(df.columns),
        ticktext=list(self.sequence),
        tickfont=dict(size=8, color='black')
    )
    fig.update_yaxes(
        title_text=temp_kwargs['y_label'],
        showline=True,
        linewidth=2,
        linecolor='black',
        ticks=None,
        mirror=True,
        dtick=1,  #frequency of ticks
        autorange="reversed",
        #tickvals = list(df.columns),
        #ticktext= ['Healthy', 'Healthy', 'Moderate', 'Diseased', 'Diseased'],
        tickfont=dict(size=8, color='black'),
    )

    # Layout and title parameters https://plotly.com/python/figure-labels/
    fig.update_layout(
        width=14 * len(df.columns) / 165 * 90,
        height=2.65 * 120,
        font=dict(family="Arial, monospace", size=12, color="black"),
        title={'text': temp_kwargs['title'], 'xanchor': 'center', 'yanchor': 'top', 'x': 0.5},
    )

    # save fig to html
    _save_html(fig, output_html)

    # return plotly object
    if temp_kwargs['return_plot_object']:
        return fig

    if temp_kwargs['show']:
        fig.show(config={'displayModeBar': False})
