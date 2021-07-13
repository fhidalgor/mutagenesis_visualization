from pandas.core.frame import DataFrame
from typing import List

try:
    from ipymol import viewer as pymol
except ModuleNotFoundError:
    pass

def pymol_fitness(df_input: DataFrame, gof, lof, mode, position_correction):
    """
    You input the dataframe. Removes stop codons.
    Returns the positions that are going to be colored blue,red and white."""

    # Select grouping
    if mode.lower() == 'mean':
        df_grouped = df_input.groupby(['Position'], as_index=False).mean()
    else:
        df_grouped = df_input.loc[df_input['Aminoacid'] == mode]

    # Color of mutations
    blue_mutations = df_grouped[df_grouped['Score'] < lof]
    red_mutations = df_grouped[df_grouped['Score'] > gof]
    white_mutations = df_grouped[df_grouped['Score'].between(lof, gof, inclusive=True)]

    # Pymol Format
    blue_pymol = _array_to_pymol(blue_mutations['Position'] + position_correction)
    red_pymol = _array_to_pymol(red_mutations['Position'] + position_correction)
    white_pymol = _array_to_pymol(white_mutations['Position'] + position_correction)

    residues = [blue_pymol, red_pymol, white_pymol]

    # If one group does not have any position, color position 0. Otherwise it gives an error
    for i, residue in enumerate(residues):
        if residue == '':
            residues[i] = '0'

    return residues

def _array_to_pymol(aminoacid_positions: List[str]):
    """
    Input an array with positions of aminoacids, return it in pymol format.
    """
    pymol = ''
    for aminoacid in aminoacid_positions:
        pymol += str(aminoacid) + '+'

    # delete last '+'
    return pymol[:-1]


def light_parameters():
    """
    Group the light and ray parameters for pymol figures.
    """
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
