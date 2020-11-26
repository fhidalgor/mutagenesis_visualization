#!/usr/bin/env python
# coding: utf-8

# # Import Modules

# In[2]:


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

# local modules
try:
    import mutagenesis_visualization.main.scripts.code_kwargs as code_kwargs
    import mutagenesis_visualization.main.scripts.code_utils as code_utils
    from mutagenesis_visualization.main.scripts.code_3D import (
        _color_3D_scatter, _parse_pdbcoordinates
    )
except ModuleNotFoundError:
    import import_notebook
    import code_kwargs
    import code_utils
    from code_3D import (_color_3D_scatter, _parse_pdbcoordinates)


# # Plotly plots

# ## Rank

# In[9]:


def plot_rank_plotly(
    self,
    mode='pointmutant',
    outdf=False,
    return_plotly_object=False,
    output_html: Union[None, str, Path] = None,
    **kwargs
):
    '''
    Generate a plotlu rank plot so every mutation/residue is sorted based on enrichment score.

    Parameters
    ----------
    self : object from class *Screen*z

    mode : str, default 'pointmutant'. 
        Alternative set to "mean" for the mean of each position

    outdf : boolean, default False
        If set to true, will return the df with the rank of mutations

    return_plotly_object : boolean, default False
        If true, will return plotly object.
        
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
    fig = px.scatter(
        x=np.arange(len(df), 0, -1), y=df['Score'], text=df['Variant']
    )

    # Style
    pio.templates.default = "plotly_white"

    # Axes https://plotly.com/python/axes/
    fig.update_traces(
        mode="markers",
        hovertemplate=
        'Position: %{x}<br>Score: %{y}<br>Variant: %{text}<extra></extra>'
    )
    fig.update_xaxes(
        title_text=temp_kwargs['x_label'],
        showline=True,
        linewidth=2,
        linecolor='black',
        ticks="outside",
        mirror=True
    )
    fig.update_yaxes(
        title_text=temp_kwargs['y_label'],
        showline=True,
        linewidth=2,
        linecolor='black',
        ticks="outside",
        mirror=True
    )

    # Layout and title parameters https://plotly.com/python/figure-labels/
    fig.update_layout(
        width=temp_kwargs['figsize'][0] * 120,
        height=temp_kwargs['figsize'][1] * 120,
        font=dict(family="Arial, monospace", size=12, color="black"),
        title={
            'text': temp_kwargs['title'], 'xanchor': 'center', 'yanchor': 'top',
            'x': 0.5
        }
    )

    if temp_kwargs['show']:
        fig.show()
    
    # save fig to html
    _save_html(fig, output_html)
    
    # return plotly object
    if return_plotly_object:
        return fig
    
    # return dataframe
    if outdf:
        return df
    
    
# save file html
def _save_html(fig, output_html):
    '''save to html'''
    if output_html:
        fig.write_html(str(Path(output_html)))


# ## Scatter

# In[ ]:


def plot_scatter_plotly(
    self,
    obj2,
    mode='pointmutant',
    show_results=False,
    return_plotly_object=False,
    output_html: Union[None, str, Path] = None,
    **kwargs
):
    '''
    Generate a scatter plot between object and a second object of the same class.

    Parameters
    ----------
    self : object from class *Screen*

    obj2 : object from class *Screen* to do the scatter with

    mode : str, default 'pointmutant'. 
        Alternative set to "mean" for the mean of each position.
    
    show_results : boolean, default False
        If set to true, will export the details of the linear fit.
    
    return_plotly_object : boolean, default False
        If true, will return plotly object.
        
    output_html : str, default None
        If you want to export the generated graph into html, add the path and name of the file.
        Example: 'path/filename.html'.
        
    **kwargs : other keyword arguments

    Returns
    ----------
    fig : plotly object

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
    fig = px.scatter(
        x=df['dataset_1'],
        y=df['dataset_2'],
        trendline="ols",
        trendline_color_override="red"
    )

    # Titles
    # hide text labels
    fig.update_traces(
        hovertext=df['Variant'],
        hovertemplate=
        'Position: %{x}<br>Score: %{y}<br>Variant: %{hovertext}<extra></extra>'
    )
    fig.update_xaxes(
        title_text=temp_kwargs['x_label'],
        showline=True,
        linewidth=2,
        linecolor='black',
        ticks="outside",
        mirror=True
    )
    fig.update_yaxes(
        title_text=temp_kwargs['y_label'],
        showline=True,
        linewidth=2,
        linecolor='black',
        ticks="outside",
        mirror=True
    )

    # Layout and title parameters https://plotly.com/python/figure-labels/
    fig.update_layout(
        width=temp_kwargs['figsize'][0] * 120,
        height=temp_kwargs['figsize'][1] * 120,
        font=dict(family="Arial, monospace", size=12, color="black"),
        title={
            'text': temp_kwargs['title'], 'xanchor': 'center', 'yanchor': 'top',
            'x': 0.5
        }
    )

    if temp_kwargs['show']:
        fig.show()
        
    if show_results:
        px.get_trendline_results(fig).px_fit_results.iloc[0].summary()
    
    # save fig to html
    _save_html(fig, output_html)
    
    # return plotly object
    if return_plotly_object:
        return fig
    


# ## 3D

# In[ ]:


def plot_scatter_3D_plotly(
    self,
    mode='mean',
    pdb_path=None,
    df_coordinates=None,
    df_color=None,
    position_correction=0,
    chain='A',
    squared=False,
    return_plotly_object=False,
    output_html: Union[None, str, Path] = None,
    **kwargs
):
    '''
    Generates a 3-D scatter plot of the x,y,z coordinates of the C-alpha atoms of the residues, 
    color coded by the enrichment scores. PDBs may have atoms missing, 
    you should fix the PDB before using this method. Use matplotlib for interactive plot.

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

    df_color : pandas dataframe, default None     
        The color of each residue can also be included. You must label that label column.

    position_correction : int, default 0
        If the pdb structure has a different numbering of positions than you dataset,
        you can correct for that. If your start_position = 2, but in the PDB that same residue
        is at position 20, position_correction needs to be set at 18.

    chain : str, default 'A'
        Chain of the PDB file to get the coordinates and SASA from.

    squared : boolean, False
        If this parameter is True, the algorithm will center the data, and plot the square value of the 
        distance.
    
    return_plotly_object : boolean, default False
        If true, will return plotly object.
        
    output_html : str, default None
        If you want to export the generated graph into html, add the path and name of the file.
        Example: 'path/filename.html'.
        
    **kwargs : other keyword arguments

    Returns
    ---------
    fig : plotly object
    '''

    # update kwargs
    temp_kwargs = copy.deepcopy(code_kwargs.kwargs())
    temp_kwargs.update(kwargs)

    # Get Scores and colors
    if df_color is None:
        df = _color_3D_scatter(
            self.dataframe, mode, temp_kwargs['lof'], temp_kwargs['gof']
        )

    # If coordinates is not an input, get it from the pdb
    if df_coordinates is None:
        df_coordinates = _parse_pdbcoordinates(
            self, pdb_path, position_correction, chain
        )

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
        color_continuous_scale='RdBu_r',
        range_color=temp_kwargs['colorbar_scale']
    )

    # Layout and title parameters https://plotly.com/python/figure-labels/
    fig.update_layout(
        scene=dict(
            xaxis=dict(
                title_text=temp_kwargs['x_label'],
                showline=True,
                linewidth=2,
                linecolor='black',
                ticks="outside",
                mirror=True,
                backgroundcolor='#9467bd',
            ),
            yaxis=dict(
                title_text=temp_kwargs['y_label'],
                showline=True,
                linewidth=2,
                linecolor='black',
                ticks="outside",
                mirror=True,
                backgroundcolor='#9467bd',
            ),
            zaxis=dict(
                title_text=temp_kwargs['z_label'],
                showline=True,
                linewidth=2,
                linecolor='black',
                ticks="outside",
                mirror=True,
                backgroundcolor='#9467bd',  #change the color of axis
            )
        )
    )
    # for the clickable part
    fig.update_traces(
        hovertext=df['Position'], hovertemplate='Position: %{hovertext}'
    )
    # title
    fig.update_layout(
        font=dict(family="Arial, monospace", size=12, color="black"),
        title={
            'text': temp_kwargs['title'], 'xanchor': 'center', 'yanchor': 'top',
            'x': 0.5
        },
        coloraxis_colorbar=dict(title="Enrichment Score") # colorbar
    )

    # show only if asked
    if temp_kwargs['show']:
        fig.show()
        
    # save fig to html
    _save_html(fig, output_html)
    
    # return plotly object
    if return_plotly_object:
        return fig

