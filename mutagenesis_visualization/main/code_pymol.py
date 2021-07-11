#!/usr/bin/env python
# coding: utf-8

# # Import Modules

# In[1]:


import numpy as np
import pandas as pd
import copy
import matplotlib.pyplot as plt
import freesasa
from os import path
from pathlib import Path
from typing import Union

# local modules
try:
    import mutagenesis_visualization.main.scripts.code_kwargs as code_kwargs
    import mutagenesis_visualization.main.scripts.code_utils as code_utils
except ModuleNotFoundError:
    import import_notebook
    import code_kwargs
    import code_utils
    
try:
    from ipymol import viewer as pymol
except ModuleNotFoundError:
    pass


# # Plot Functions

# ## Map into Pymol

# In[ ]:


def plot_pymol(
    self,
    pdb,
    mode='mean',
    residues=None,
    position_correction=0,
    quit=False,
    output_file: Union[None, str, Path] = None,
    **kwargs
):
    """
    Color pymol structure residues. User can specify the residues to color, or
    can use the mutagenesis data. Activating mutations will be colored red and
    loss of function blue. Neutral mutations in green. Only works if pymol is
    your $PATH as pymol or you can start PyMOL in server mode. Uses the ipymol
    package, which needs to be installed from Github $pip install
    git+https://github.com/cxhernandez/ipymol , not from pypi (not updated
    there).

    Parameters
    ----------
    pdb : str
        User should specify the PDB chain in the following format 4G0N_A.
        If you have internet connection, Pymol will download the pdb. Otherwise,
        include the path were your PDB is stored locally.

    mode : str, default 'mean'
        Specify what enrichment scores to use. If mode = 'mean', it will use the mean of
        each position to classify the residues. If mode = 'A', it will use the Alanine substitution profile. 
        Can be used for each amino acid. Use the one-letter code and upper case.

    residues : list , optional
        If user decides to pass custom arguments, use the following format
        residues = ['1,2,3,4-10','12-15,23,24,35','48,49,50,52-60'] which are [blue,red,green].

    position_correction : int, default 0
        If the pdb structure has a different numbering of positions than you dataset,
        you can correct for that. If your start_position = 2, but in the PDB that same residue
        is at position 20, position_correction needs to be set at 18.

    quit : boolean, default False
        if quit, close pymol after executing code.

    output_file : str, default None
        If you want to export the generated graph, add the path and name of the file.
        Example: 'path/filename.png' or 'path/filename.svg'.

    **kwargs : other keyword arguments
        gof : int, default is 1
             cutoff for determining gain of function mutations based on mutagenesis data.
        lof : int, default is -1
             cutoff for determining loss of function mutations based on mutagenesis data.
        color_gof : str, default 'red'
            Choose color to color positions with an enrichment score > gof.
        color_lof : str, default 'neptunium'
            Choose color to color positions with an enrichment score < lof.


    Returns
    ----------
    Open pymol session with a fetched pdb structure where the residues are colored according to the enrichment scores.

    """
    # update kwargs
    temp_kwargs = copy.deepcopy(code_kwargs.kwargs())
    temp_kwargs.update(kwargs)
    temp_kwargs['color_lof'] = kwargs.get('color_lof', 'neptunium')

    # Calculate residues only if they are not given by the user
    if residues is None:
        residues = _pymol_fitness(
            self.dataframe.copy(), temp_kwargs['gof'], temp_kwargs['lof'], mode,
            position_correction
        )

    # Start Pymol
    if not pymol._process_is_running():
        pymol.start()

    # Fetch structure. If pdb contains a "/", it will assume it is stored locally
    if '/' in pdb:
        pymol.load(pdb)
        pdb = (path.basename(pdb)).partition('.')[
            0]  # Extract filename from pdb and then extract pdb code
    else:
        pymol.fetch(pdb)

    # Hide everything
    pymol.do('hide everything')

    # Selection names
    blue = pdb + '_blue'
    red = pdb + '_red'
    white = pdb + '_white'

    # Do selections
    pymol.select(blue, 'resi ' + residues[0])
    pymol.select(red, 'resi ' + residues[1])
    pymol.select(white, 'resi ' + residues[2])

    # Representation parameters
    pymol.show_as('cartoon', pdb)
    pymol.set('cartoon_color', temp_kwargs['color_lof'], blue)
    pymol.set('cartoon_color', temp_kwargs['color_gof'], red)
    pymol.set('cartoon_color', 'chlorine', white)
    pymol.bg_color('white')
    pymol.remove('solvent')

    # light parameters
    _light_parameters()

    # deselect everything
    pymol.deselect()

    if quit:
        pymol.quit()
    return


# Convert fitness scores into pymol residues


def _pymol_fitness(df, gof, lof, mode, position_correction):
    '''You input the dataframe. Removes stop codons. 
    Returns the positions that are going to be colored blue,red and white'''

    # Select grouping
    if mode == 'mean':
        df_grouped = df.groupby(['Position'], as_index=False).mean()
    else:
        df_grouped = df.loc[df['Aminoacid'] == mode]

    # Color of mutations
    blue_mutations = df_grouped[df_grouped['Score'] < lof]
    red_mutations = df_grouped[df_grouped['Score'] > gof]
    white_mutations = df_grouped[
        df_grouped['Score'].between(lof, gof, inclusive=True)]

    # Pymol Format
    blue_pymol = _array_to_pymol(
        blue_mutations['Position'] + position_correction
    )
    red_pymol = _array_to_pymol(red_mutations['Position'] + position_correction)
    white_pymol = _array_to_pymol(
        white_mutations['Position'] + position_correction
    )

    residues = [blue_pymol, red_pymol, white_pymol]

    # If one group does not have any position, color position 0. Otherwise it gives an error
    for i, residue in enumerate(residues):
        if residue == '':
            residues[i] = '0'

    return residues


def _array_to_pymol(array):
    '''Input an array with positions of aminoacids, return it in pymol format'''
    pymol = ''
    for aminoacid in array:
        pymol += str(aminoacid) + '+'

    # delete last '+'
    pymol = pymol[:-1]
    return pymol


def _light_parameters():
    '''Group the light and ray parameters for pymol figures'''
    # Light parameters
    pymol.set('antialias', '3')
    pymol.set('ambient', '0.15')
    pymol.set('spec_count', '5')
    pymol.set('shininess', '50')
    pymol.set('specular', '0')
    pymol.set('light_count', '4')
    pymol.set('direct', '0.45')
    pymol.set('reflect', '0.5')
    pymol.set('opaque_background', 'off')
    pymol.set('dash_gap', 0.5)
    pymol.set('dash_radius', 0.1)

    # Stick parameters
    pymol.set('stick_radius', '0.2')
    pymol.set('sphere_scale', '0.2')
    pymol.set('sphere_quality', '4')
    return

