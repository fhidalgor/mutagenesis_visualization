"""
This module contains utils for the plotly figures.
"""

from typing import Dict, Any, Tuple
import numpy as np
from Bio.PDB import PDBParser
import freesasa
from pandas.core.frame import DataFrame
from plotly.graph_objects import Figure


def select_grouping(dataframe: DataFrame, mode: str) -> DataFrame:
    """
    Choose the subset of substitutions based on mode input.
    For example, if mode=='A', then return data for Alanine.
    """
    # Select grouping
    if mode.upper() == 'MEAN':
        return dataframe.groupby('Position', as_index=False).mean()
    return dataframe.loc[dataframe['Aminoacid'] == mode.upper()].copy()


def color_3d_scatter(df_input: DataFrame, mode: str, lof: float, gof: float) -> DataFrame:
    """
    Color the data points by enrichment scores.

    Parameters
    -----------
    df : pandas dataframe
        The input is a dataframe that has colum with
        ['Position', 'Aminoacid', 'Score'].

    mode : str
        Specify what enrichment scores to use. If mode = 'mean', it will
        use the mean of each position to classify the residues. If
        mode = 'A', it will use the Alanine substitution profile. Can be
        used for each amino acid. Use the one-letter code and upper case.

    gof : int, default is 1
        cutoff for determining gain of function mutations based on
        mutagenesis data.

    lof : int, default is -1
        cutoff for determining loss of function mutations based on
        mutagenesis data.

    Returns
    ---------
    df_grouped: pandas dataframe
        New dataframe with added column of ['Color'] and the ['Score']
        values of the mode you chose.
    """

    # Copy df
    df_grouped: DataFrame = df_input.copy()

    # Select grouping.
    if mode.lower() == 'mean':
        df_grouped = df_grouped.groupby(['Position'], as_index=False).mean()
    else:
        df_grouped = df_grouped.loc[df_grouped['Aminoacid'] == mode]

    # Select colors based on Score values
    df_grouped['Color'] = 'green'
    df_grouped.loc[df_grouped['Score'] < lof, 'Color'] = 'blue'
    df_grouped.loc[df_grouped['Score'] > gof, 'Color'] = 'red'
    return df_grouped


def centroid(df_input: DataFrame) -> Tuple:
    """
    Find center of x,y,z using centroid method.
    The input is a dataframe with columns x, y, z.
    Returns the center of each of the three dimensions
    """
    return df_input['x'].mean(), df_input['y'].mean(), df_input['z'].mean()


def parse_pdb_coordinates(
    pdb_path: str,
    start_position: int,
    end_position: int,
    position_correction: int,
    chain: str,
    sasa: bool = False
) -> DataFrame:
    """
    Parse coordinate of CA atoms. Will also return the bfactor and SASA using freesasa.
    If PDB is missing atoms, it can handle it.
    """

    # Get structure from PDB
    structure = PDBParser().get_structure('pdb', pdb_path)

    coordinates = []
    commands = []
    bfactors = []
    positions_worked = []  # positions present in pdb

    # Iterate over each CA atom and geet coordinates
    for i in np.arange(start_position + position_correction, end_position + position_correction):
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
            bfactors.append(np.log10(bfactor))
            positions_worked.append(i)
        except:
            print("residue {} not found".format(str(i)))
            coordinates.append([np.nan, np.nan, np.nan, i])

    # Convert to df
    df_coordinates = DataFrame(columns=['x', 'y', 'z', 'Position'], data=coordinates)

    # Center data
    x, y, z = centroid(df_coordinates)
    df_coordinates['x_cent'] = (df_coordinates['x'] - x).abs()**2
    df_coordinates['y_cent'] = (df_coordinates['y'] - y).abs()**2
    df_coordinates['z_cent'] = (df_coordinates['z'] - z).abs()**2
    df_coordinates[
        'Distance'] = df_coordinates['x_cent'] + df_coordinates['y_cent'] + df_coordinates['z_cent']

    # Add sasa values
    if sasa:
        # Get structure for SASA
        structure_sasa = freesasa.Structure(pdb_path)
        result = freesasa.calc(structure_sasa)
        # Calculate sasa
        sasa_area = freesasa.selectArea(commands, structure_sasa, result)
        df_sasa: DataFrame = DataFrame(columns=['SASA'], data=sasa_area.values())
        df_sasa['log B-factor'] = bfactors
        df_sasa['Position'] = positions_worked

        # Merge
        df_coordinates = df_coordinates.merge(df_sasa, how='outer', on='Position')

    return df_coordinates


def matplotlib_to_plotly(cmap: Any, pl_entries: int = 255) -> list:
    """
    Convert a matplotlib colorscale into plotly rgb format.
    """
    h = 1.0 / (pl_entries - 1)
    pl_colorscale = []

    for k in range(pl_entries):
        colors = list(map(np.uint8, np.array(cmap(k * h)[: 3]) * 255))
        pl_colorscale.append([k * h, 'rgb' + str((colors[0], colors[1], colors[2]))])

    return pl_colorscale


def update_layout(fig: Figure, temp_kwargs: Dict[str, Any]) -> None:
    """
    Update layout of plotly figure.
    """
    fig.update_layout(
        font=dict(family="Arial, monospace", size=12, color="black"),
        title={
            'text': temp_kwargs['title'], 'xanchor': 'center', 'yanchor': 'top',
            'x': 0.5
        },
        coloraxis_colorbar=dict( # modify colorbar properties
        title = 'Fitness',
        len=0.4,
        thickness=15,
        outlinewidth=2,
        outlinecolor='rgb(0,0,0)',
        showticklabels=True,
        )
    )


def update_axes(fig: Figure, temp_kwargs: Dict[str, Any]) -> None:
    """
    Separeated this portion of the code because it is clumsy. It changes
    the axes looks.
    """
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
