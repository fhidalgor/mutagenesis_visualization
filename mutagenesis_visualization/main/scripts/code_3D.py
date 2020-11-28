#!/usr/bin/env python
# coding: utf-8

# # Import Modules

# In[2]:


import numpy as np
import pandas as pd
import copy
import matplotlib.pyplot as plt
import freesasa
from os import path
from pathlib import Path
from typing import Union
from Bio.PDB import PDBParser

# Local imports
try:
    import mutagenesis_visualization.main.scripts.code_kwargs as code_kwargs
    import mutagenesis_visualization.main.scripts.code_utils as code_utils
except ModuleNotFoundError:
    import import_notebook
    import code_kwargs
    import code_utils


# # Plot Functions

# ## 3D plot

# ### 3D Scatter

# In[3]:


def plot_scatter_3D(
    self,
    mode='mean',
    pdb_path=None,
    df_coordinates=None,
    df_color=None,
    position_correction=0,
    chain='A',
    squared=False,
    output_file: Union[None, str, Path] = None,
    **kwargs
):
    """
    DEPRECATED.

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

    squared : booleand, False
        If this parameter is True, the algorithm will center the data, and plot the square value of the
        distance.

    output_file : str, default None
        If you want to export the generated graph, add the path and name of the file.
        Example: 'path/filename.png' or 'path/filename.svg'.

    **kwargs : other keyword arguments
        return_plot_object : boolean, default False
            If true, will return plotting objects (ie. fig, ax).
        gof : int, default is 1
                 cutoff for determining gain of function mutations based on mutagenesis data.
        lof : int, default is -1
            cutoff for determining loss of function mutations based on mutagenesis data.

    Returns
    ---------
    fig, ax : matplotlib figure and subplots
        Needs to have return_plot_object==True. By default they do
        not get returned.
        
    """

    # Load parameters
    code_kwargs._parameters()

    # Update kwargs
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
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    if squared is False:
        ax.scatter(
            df_coordinates['x'],
            df_coordinates['y'],
            df_coordinates['z'],
            c=df['Color']
        )
    else:
        ax.scatter(
            df_coordinates['x_cent'],
            df_coordinates['y_cent'],
            df_coordinates['z_cent'],
            c=df['Color']
        )

    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')

    # save file
    code_utils._save_work(fig, output_file, temp_kwargs)
    
    # return matplotlib object
    if temp_kwargs['return_plot_object']:
        return fig, ax
    
    if temp_kwargs['show']:
        plt.show()


def _color_3D_scatter(df, mode, lof, gof):
    '''
    Color the data points by enrichment scores.
    
    Parameters
    -----------
    df : pandas dataframe
        The input is a dataframe that has colum with ['Position', 'Aminoacid', 'Score'].
    
    mode : str
        Specify what enrichment scores to use. If mode = 'mean', it will use the mean of 
        each position to classify the residues. If mode = 'A', it will use the Alanine substitution profile. 
        Can be used for each amino acid. Use the one-letter code and upper case.

    gof : int, default is 1
        cutoff for determining gain of function mutations based on mutagenesis data.
    
    lof : int, default is -1
        cutoff for determining loss of function mutations based on mutagenesis data.
    
    Returns
    ---------
    df_grouped: pandas dataframe
        New dataframe with added column of ['Color'] and the ['Score'] values of the 
        mode you chose.
        
    '''

    # Copy df
    df_grouped = df.copy()

    # Select grouping.
    if mode == 'mean':
        df_grouped = df_grouped.groupby(['Position'], as_index=False).mean()
    else:
        df_grouped = df_grouped.loc[df_grouped['Aminoacid'] == mode]

    # Select colors based on Score values
    df_grouped['Color'] = 'green'
    df_grouped.loc[df_grouped['Score'] < lof, 'Color'] = 'blue'
    df_grouped.loc[df_grouped['Score'] > gof, 'Color'] = 'red'

    return df_grouped


def _centroid(df):
    '''
    Find center of x,y,z using centroid method. 
    
    The input is a dataframe with columns x, y, z.
    Returns the center of each of the three dimensions
    
    '''
    return df['x'].mean(), df['y'].mean(), df['z'].mean()


def _parse_pdbcoordinates(
    self, pdb_path, position_correction, chain, sasa=False
):
    """
    parse coordinate of CA atoms.

    Will also return the bfactor and SASA using freesasa. If PDB is
    missing atoms, it can handle it.

    """

    # Get structure from PDB
    structure = PDBParser().get_structure('pdb', pdb_path)

    coordinates = []
    commands = []
    bfactors = []
    positions_worked = []  # positions present in pdb

    # Iterate over each CA atom and geet coordinates
    for i in np.arange(self.start_position + position_correction,
                       self.end_position + position_correction):
        # first check if atom exists
        try:
            structure[0][chain][int(i)].has_id("CA")
            # Get atom from pdb and geet coordinates
            atom = list(structure[0][chain][int(i)]["CA"].get_vector()) + [i]
            coordinates.append(atom)
            # Get SASA command for each residue and bfactor
            residue = "s{}, chain {} and resi {}".format(str(i), chain, str(i))
            commands.append(residue)
            bfactor = (structure[0][chain][int(i)]["CA"].get_bfactor())
            bfactors.append(bfactor)
            positions_worked.append(i)
        except:
            print("residue {} not found".format(str(i)))
            coordinates.append([np.nan, np.nan, np.nan, i])

    # Convert to df
    df_coordinates = pd.DataFrame(
        columns=['x', 'y', 'z', 'Position'], data=coordinates
    )

    # Center data
    x, y, z = _centroid(df_coordinates)
    df_coordinates['x_cent'] = (df_coordinates['x'] - x).abs()**2
    df_coordinates['y_cent'] = (df_coordinates['y'] - y).abs()**2
    df_coordinates['z_cent'] = (df_coordinates['z'] - z).abs()**2
    df_coordinates['Distance'] = df_coordinates['x_cent'] +         df_coordinates['y_cent']+df_coordinates['z_cent']

    # Add sasa values
    if sasa:
        # Get structure for SASA
        structure_sasa = freesasa.Structure(pdb_path)
        result = freesasa.calc(structure_sasa)
        # Calculate sasa
        sasa = freesasa.selectArea(commands, structure_sasa, result)
        df_sasa = pd.DataFrame(columns=['SASA'], data=sasa.values())
        df_sasa['B-factor'] = bfactors
        df_sasa['Position'] = positions_worked

        # Merge
        df_coordinates = df_coordinates.merge(
            df_sasa, how='outer', on='Position'
        )

    return df_coordinates


# ### 3D Scatter Second version

# In[6]:


def plot_scatter_3D_pdbprop(
    self,
    plot=['Distance', 'SASA', 'B-factor'],
    mode='mean',
    pdb_path=None,
    custom=None,
    axis_scale=["linear", "linear", "linear"],
    df_color=None,
    color_by_score=True,
    position_correction=0,
    chain='A',
    output_df=False,
    output_file: Union[None, str, Path] = None,
    **kwargs
):
    """
    DEPRECATED.

    Generates a 3-D scatter plot of different properties obtained from the PDB.
    PDBs may have atoms missing, you should fix the PDB before using this
    method. We recommend you use matplotlib for interactive plot.

    Parameters
    -----------
    self : object from class *Screen*
        **kwargs : other keyword arguments.

    plot : list, default ['Distance', 'SASA', 'B-factor']
        List of 3 elements to plot. Other options are 'Score' and Custom. If custom, add the
        label to the third element of the list ie ['Distance', 'SASA', 'Conservation'].

    mode : str, default 'mean'
        Specify what enrichment scores to use. If mode = 'mean', it will use the mean of
        each position to classify the residues. If mode = 'A', it will use the Alanine substitution profile. Can be
        used for each amino acid. Use the one-letter code and upper case.

    pdb_path : str, default None
        User should specify the path PDB.

    custom : list or dataframe or np.array, default None
        If you want to add a custom dataset to plot, use custom. On the parameter
        plot, the 3rd item of the list will be the label for your custom dataset.

    axis_scale : list, default ["linear", "linear", "linear"]
        Check matplotlib.axes.Axes.set_xscale documentation for more information.
        The axis scale type to apply. Some options are {"linear", "log", "symlog", "logit", ...}.

    df_color : pandas dataframe, default None
        The color of each residue can also be included. You must label that label column.

    color_by_score : boolean, default True
        If set to False, the points in the scatter will not be colored based on the enrichment score.

    position_correction : int, default 0
        If the pdb structure has a different numbering of positions than you dataset,
        you can correct for that. If your start_position = 2, but in the PDB that same residue
        is at position 20, position_correction needs to be set at 18.

    chain : str, default 'A'
        Chain of the PDB file to get the coordinates and SASA from.

    output_df : boolean, default False
        If true, this method will return the dataframe with the data.

    output_file : str, default None
        If you want to export the generated graph, add the path and name of the file.
        Example: 'path/filename.png' or 'path/filename.svg'.

    **kwargs : other keyword arguments
        return_plot_object : boolean, default False
            If true, will return plotting objects (ie. fig, ax).
        gof : int, default is 1
                 cutoff for determining gain of function mutations based on mutagenesis data.
        lof : int, default is -1
            cutoff for determining loss of function mutations based on mutagenesis data.

    Returns
    ---------
    fig, ax : matplotlib figure and subplots
        Needs to have return_plot_object==True. By default they do
        not get returned.
        
    df_items : pandas dataframe
        Contains the plotted data. Needs to have output_df set to true.

    """
    # Load parameters
    code_kwargs._parameters()

    # Update kwargs
    temp_kwargs = copy.deepcopy(code_kwargs.kwargs())
    temp_kwargs.update(kwargs)
    temp_kwargs['x_label'] = kwargs.get('x_label', plot[0])
    temp_kwargs['y_label'] = kwargs.get('y_label', plot[1])
    temp_kwargs['z_label'] = kwargs.get('z_label', plot[2])

    # Get Scores and colors
    df_scores = _color_3D_scatter(
        self.dataframe, mode, temp_kwargs['lof'], temp_kwargs['gof']
    )

    # If coordinates is not an input, get it from the pdb
    df_items = _parse_pdbcoordinates(
        self, pdb_path, position_correction, chain, sasa=True
    )

    # Add scores
    df_items['Score'] = list(df_scores['Score'])

    # Custom data
    if custom is not None:
        df_items[plot[2]] = custom

    # Plot figure
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    if df_color is None and color_by_score is True:
        c = df_scores['Color']
    elif df_color is None and color_by_score is False:
        c = 'k'
    else:
        c = df_scores['Color']

    ax.scatter(df_items[plot[0]], df_items[plot[1]], df_items[plot[2]], c=c)

    # axis labels
    ax.set_xlabel(temp_kwargs['x_label'])
    ax.set_ylabel(temp_kwargs['y_label'])
    ax.set_zlabel(temp_kwargs['z_label'])

    # axis scales
    ax.set_xscale(axis_scale[0])
    ax.set_yscale(axis_scale[1])
    ax.set_zscale(axis_scale[2])

    # save file
    code_utils._save_work(fig, output_file, temp_kwargs)
    
    # return matplotlib object
    if temp_kwargs['return_plot_object']:
        return fig, ax
    
    if temp_kwargs['show']:
        plt.show()

    if output_df:
        return df_items, df_scores

