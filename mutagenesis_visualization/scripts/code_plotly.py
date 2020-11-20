#!/usr/bin/env python
# coding: utf-8

# # Import Modules

# In[1]:


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

# kwargs and parameters
try:
    import import_notebook
except ModuleNotFoundError:
    pass

import code_kwargs
import code_utils


# # Plot Functions

# ## Plotly plots

# ### Rank

# In[2]:


def plot_rank_plotly(self, mode='pointmutant', outdf=False,
                     output_file: Union[None, str, Path] = None,
                     **kwargs):
    '''
    Generate a plotlu rank plot so every mutation/residue is sorted based on enrichment score.

    Parameters
    ----------
    self : object from class *Screen*z

    mode : str, default 'pointmutant'. 
        Alternative set to "mean" for the mean of each position

    outdf : boolean, default False
        If set to true, will return the df with the rank of mutations

    output_file : str, default None
        If you want to export the generated graph, add the path and name of the file.
        Example: 'path/filename.png' or 'path/filename.svg'. 

    **kwargs : other keyword arguments

    Returns
    ----------
    Pandas dataframe
    '''
    # update kwargs
    temp_kwargs = copy.deepcopy(code_kwargs.kwargs())
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
    fig = px.scatter(x=np.arange(len(df), 0, -1), y=df['Score'],
                     text=df['Variant'])

    # Style
    pio.templates.default = "plotly_white"

    # Axes https://plotly.com/python/axes/
    fig.update_traces(mode="markers",
                      hovertemplate='Position: %{x}<br>Score: %{y}<br>Variant: %{text}<extra></extra>')
    fig.update_xaxes(title_text=temp_kwargs['x_label'], showline=True,
                     linewidth=2, linecolor='black', ticks="outside", mirror=True)
    fig.update_yaxes(title_text=temp_kwargs['y_label'], showline=True,
                     linewidth=2, linecolor='black', ticks="outside", mirror=True)

    # Layout and title parameters https://plotly.com/python/figure-labels/
    fig.update_layout(width=temp_kwargs['figsize'][0]*120,
                      height=temp_kwargs['figsize'][1]*120,
                      font=dict(family="Arial, monospace", size=12,
                                color="black"),
                      title={'text': temp_kwargs['title'],
                             'xanchor': 'center', 'yanchor': 'top', 'x': 0.5})

    if temp_kwargs['show']:
        fig.show()

    if outdf:
        return df


# ### Scatter

# In[3]:


def plot_scatter_plotly(self, obj2, mode='pointmutant',
                        output_file: Union[None, str, Path] = None,
                        show_results=False, **kwargs):
    '''
    Generate a scatter plot between object and a second object of the same class.

    Parameters
    ----------
    self : object from class *Screen*

    obj2 : object from class *Screen* to do the scatter with

    mode : str, default 'pointmutant'. 
        Alternative set to "mean" for the mean of each position

    output_file : str, default None
        If you want to export the generated graph, add the path and name of the file.
        Example: 'path/filename.png' or 'path/filename.svg'. 
    
    show_results : boolean, default False
        If set to true, will export the details of the linear fit.
    
    **kwargs : other keyword arguments

    Returns
    ----------
    None.
    '''

    # update kwargs
    temp_kwargs = copy.deepcopy(code_kwargs.kwargs())
    temp_kwargs.update(kwargs)
    temp_kwargs['figsize'] = kwargs.get('figsize', (4, 3))

    # Chose mode:
    if mode == 'pointmutant':
        df = code_utils._process_bypointmutant(self, obj2)
    else:
        df = code_utils._process_meanresidue(self, obj2)
        df['Variant'] = df['Position']

    # Style
    pio.templates.default = "plotly_white"

    # create figure
    fig = px.scatter(x=df['dataset_1'], y=df['dataset_2'],
                     trendline="ols", trendline_color_override="red")

    # Titles
    # hide text labels
    fig.update_traces(
        hovertext=df['Variant'], hovertemplate='Position: %{x}<br>Score: %{y}<br>Variant: %{hovertext}<extra></extra>')
    fig.update_xaxes(title_text=temp_kwargs['x_label'], showline=True,
                     linewidth=2, linecolor='black', ticks="outside", mirror=True)
    fig.update_yaxes(title_text=temp_kwargs['y_label'], showline=True,
                     linewidth=2, linecolor='black', ticks="outside", mirror=True)

    # Layout and title parameters https://plotly.com/python/figure-labels/
    fig.update_layout(width=temp_kwargs['figsize'][0]*120,
                      height=temp_kwargs['figsize'][1]*120,
                      font=dict(family="Arial, monospace", size=12,
                                color="black"),
                      title={'text': temp_kwargs['title'],
                             'xanchor': 'center', 'yanchor': 'top', 'x': 0.5})

    if temp_kwargs['show']:
        fig.show()
    
    if show_results:
        return px.get_trendline_results(fig).px_fit_results.iloc[0].summary()
        

