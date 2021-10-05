"""
Auxiliar functions used in the heatmaps.
"""
from typing import List, Tuple, Union
from collections import Counter
import numpy as np
from numpy import typing as npt
from pandas.core.frame import DataFrame
from scipy.cluster import hierarchy
from matplotlib import patches
from matplotlib import pyplot as plt
from matplotlib import gridspec


def labels(start_position: int = 1) -> Tuple[List[str], List[Union[str, int]]]:
    """
    Write heatmap labels
    """
    # residue label and color
    number_sequencelabels: List[str] = list([
        'b' if index in np.arange(10 - (start_position % 10), 1000, 10) else 'k'
        for index in range(0, 1000)
    ])
    color_sequencelabels: List[Union[str, int]] = list([
        index + start_position if index in np.arange(10 - (start_position % 10), 1000, 10) else ''
        for index in range(0, 1000)
    ])
    return number_sequencelabels, color_sequencelabels


def hierarchical_sort(df_input: DataFrame) -> npt.NDArray:
    """
    Sorts columns of dataset using hierarchical clustering and returns
    order to rearrange dataset.
    """

    # replaces NaN values with 0
    df_nan_replaced = df_input.fillna(0)

    # give sorted column order
    hierarchy_z = hierarchy.ward(df_nan_replaced.T)
    new_order: npt.NDArray = hierarchy.leaves_list(hierarchy_z)
    return new_order


def generate_cartoon(
    secondary: list,
    start_position: int,
    gs_object: gridspec.GridSpec,
    n_row: int,
    colors: List[str],
    bottom_space: float = 0,
    fig_inches: float = 13.91,
    show_labels: bool = True,
) -> plt.Axes:
    """
    Generates cartoon for heatmap.
    """
    # Create subplot
    cartoon = plt.subplot(gs_object[n_row, 0])

    # Generate coordinates of labels
    secondary_labels = list(Counter(secondary).keys())
    secondary_lengths = list(Counter(secondary).values())
    cumsum = secondary_lengths[:-1]
    cumsum.insert(0, start_position)

    # Create cartoon
    for label, length, cum in zip(secondary_labels, secondary_lengths, np.cumsum(cumsum)):
        if 'β' in label:
            loopstructure = _loop(cum, length, color=colors[2])
            cartoon.add_patch(loopstructure)
            sheetstructure = _sheet(cum, length, colors[0])
            cartoon.add_patch(sheetstructure)
            x_label = cum + length - 3.5
            if length > 2 and show_labels:  # If beta sheet is too small, label does not fit
                if length == 3:
                    cartoon.text((x_label + 0.6),
                                 -0.25,
                                 label,
                                 fontweight='normal',
                                 size=8.5 * fig_inches / 13.91,
                                 multialignment='right')
                else:
                    cartoon.text((x_label),
                                 -0.25,
                                 label,
                                 fontweight='normal',
                                 size=8.5 * fig_inches / 13.91,
                                 multialignment='right')
        elif 'α' in label:
            helixstructure = _helix(cum, length, colors[1])
            cartoon.add_patch(helixstructure)
            x_label = cum + length / 2 - 1
            if length > 2 and show_labels:
                cartoon.text((x_label),
                             -0.3,
                             label,
                             fontweight='normal',
                             size=9 * fig_inches / 14,
                             multialignment='center')
        elif 'L' in label:
            loopstructure = _loop(cum, length, colors[2])
            cartoon.add_patch(loopstructure)

    # format of secondary cartoon
    cartoon.xaxis.set_ticks_position('none')
    cartoon.yaxis.set_ticks_position('none')
    cartoon.axis('off')

    # size
    cartoon.set_xlim(start_position - 0.1, len(secondary) + start_position + 0.2)
    cartoon.set_ylim(-2, 2.5)

    # adjust proximity to heatmap
    box = cartoon.get_position()
    box.y0 = box.y0 - bottom_space
    box.y1 = box.y1 - bottom_space
    cartoon.set_position(box)


def _sheet(starting_aa: float, length_aa: int, color: str = 'lightgreen') -> patches.FancyArrow:
    return patches.FancyArrow(
        starting_aa,
        0.25,
        length_aa,
        0,
        width=2,
        length_includes_head=True,
        head_width=4,
        head_length=3,
        shape='full',
        overhang=0,
        head_starts_at_zero=False,
        ec='k',
        fc=color
    )


def _helix(starting_aa: float, length_aa: int, color: str = 'lavender') -> plt.Rectangle:
    """
    Produces matplotlib.Rectangle of length specified by length_aa.
    Used for the helix lines.
    """
    return plt.Rectangle((starting_aa, -0.85), length_aa, 2.2, fc=color, ec='k')


def _loop(starting_aa: float, length_aa: int, color: str = 'k') -> plt.Rectangle:
    """
    Produces matplotlib.Rectangle of length specified by length_aa.
    Used for the loop lines.
    """
    return plt.Rectangle((starting_aa, 0), length_aa, 0.5, fc=color)


def add_border_self_substitution(
    ax_object, sequence: str, aminoacids_order: Union[str, List[str]], color: str, lw: float
) -> None:
    """
    Add a border to each square that belongs to a self-substitution
    """
    if isinstance(aminoacids_order, list):
        aminoacids_order = "".join(aminoacids_order)

    for i, aminoacid in enumerate(sequence):
        # Check if aminoacid in the sequence is in the positions
        if aminoacid in aminoacids_order:
            ax_object.add_patch(
                patches.Rectangle((i, aminoacids_order.index(aminoacid)),
                                  1,
                                  1,
                                  edgecolor=color,
                                  fill=False,
                                  lw=lw)
            )
