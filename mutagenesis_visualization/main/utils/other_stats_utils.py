"""
Utilities used in the other stats methods.
"""
from typing import List, Tuple

from copy import deepcopy
import numpy as np
from numpy import typing as npt
from pandas import merge, Categorical
from pandas.core.frame import DataFrame
from sklearn.metrics import roc_auc_score, roc_curve

from mutagenesis_visualization.main.utils.pandas_functions import return_common_elements


def select_grouping(df_input: DataFrame, mode: str) -> DataFrame:
    """
    Choose the subset of substitutions based on mode input.
    For example, if mode=='A', then return data for Alanine.
    """

    # Select grouping
    if mode.upper() == 'POINTMUTANT':
        return df_input
    if mode.upper() == 'MEAN':
        return df_input.groupby('Position', as_index=False).mean()
    return df_input.loc[df_input['Aminoacid'] == mode].copy()


def roc_auc(df: DataFrame) -> Tuple:
    """
    Calculate roc rates and auc.

    The input is a dataframe that contains [Variants,Class,Score]
    """
    fpr, tpr, thresholds = roc_curve(df['Class'], df['Score'], drop_intermediate=True)
    auc = roc_auc_score(df['Class'], df['Score'])
    return fpr, tpr, auc, thresholds


def merge_class_variants(df_score: DataFrame, df_class: DataFrame, mode: str) -> DataFrame:
    """
    Merge the input dataframe containing the class (true score) for
    variants and the enrichment scores
    """

    if mode.upper() == 'POINTMUTANT':
        # Cut other data
        df_class = df_class[['Variant', 'Class']].copy()
        # Merge DMS with true score dataset
        df_merged = merge(df_score, df_class, on=['Variant'], how='left')
    else:
        # Cut other data
        df_class = df_class[['Position', 'Class']].copy()
        df_merged = merge(df_score, df_class, on=['Position'], how='left')

    # Drop rows with Nan values
    df_merged.dropna(inplace=True)

    return df_merged


def condense_heatmap(df_input: DataFrame, new_order: List[str]) -> DataFrame:
    """
    Converts the np.array with stored enrichment scores into the condensed heatmap
    """
    df_input = df_input.copy()
    df_input.drop(['Position'], axis=1, inplace=True)

    # Group by sequence and aminoacid, and then pivot table
    df_grouped = df_input.groupby(['Sequence', 'Aminoacid'], sort=False).mean()
    df_pivoted = df_grouped.pivot_table(values='Score', index='Aminoacid', columns='Sequence')
    df_pivoted.reset_index(drop=False, inplace=True)

    # Sort in y axis desired order
    df_pivoted['Aminoacid'] = Categorical(df_pivoted['Aminoacid'], new_order)
    df_pivoted = df_pivoted.sort_values(by=['Aminoacid'])

    # Sort in x axis desired order
    x_order = return_common_elements(new_order, list(df_pivoted.columns))

    # Drop amino acid column
    data_dropped = df_pivoted.drop(['Aminoacid'], axis=1)

    return data_dropped[x_order]


def _offset_sequence(
    dataset: npt.NDArray, sequence: str, start_position: int, position_offset: int
) -> str:
    """
    Internal function that offsets the input sequence.

    Returns
    --------
    string containing trimmed sequence    """
    # Deep copy sequence
    sequence = deepcopy(sequence)

    # truncate sequence
    if position_offset > 0:
        sequence = sequence + 'X' * np.absolute(position_offset)
        trimmedsequence = sequence[start_position - 1 + position_offset : len(dataset[0]) +
                                   start_position - 1 + position_offset]
    else:
        sequence = 'X' * (np.absolute(position_offset)) + sequence
        trimmedsequence = sequence[start_position - 1 : len(dataset[0]) + start_position - 1]

    return trimmedsequence


def transform_dataset_offset(
    dataset: npt.NDArray,
    dataframe: DataFrame,
    dataframe_stopcodons: DataFrame,
    sequence_raw: str,
    start_position: int,
    position_offset: int,
    stopcodons: bool,
) -> DataFrame:
    """
    Generate a dataframe with the sequence position_offset.
    """
    # Add position_offset sequence
    offset_sequence = _offset_sequence(dataset, sequence_raw, start_position, position_offset)
    df_output = dataframe_stopcodons.copy() if stopcodons is True else dataframe.copy()

    # Copy old sequence
    df_output['Sequence_old'] = df_output['Sequence']
    # Count amino acids
    aa_number = len(set(df_output['Aminoacid']))
    # Generate new position_offset sequence
    df_output['Sequence'] = np.ravel([[aa] * aa_number for aa in offset_sequence])

    # Drop rows with X
    df_output.drop(df_output.index[df_output['Sequence'] == 'X'], inplace=True)
    return df_output


def normalize_neighbor_effect(
    df_input: DataFrame, df_condensed_heatmap: DataFrame, neworder_aminoacids: List[str],
    aminoacids: List[str]
) -> DataFrame:
    """
    For every residue, subtract the average effect of a substitution
    Returns a normalized dataframe
    """

    df_output: DataFrame = DataFrame()
    for aa in aminoacids:
        # Choose the neighbors of an aa
        aa_neighbors = df_input.loc[df_input['Sequence'] == aa]
        # Do the mean substitution of amino acids that are repeated
        aa_neighbors = aa_neighbors.groupby(['Sequence_old', 'Aminoacid'], as_index=False).mean()
        # Make into table
        aa_neighbors_pivoted = aa_neighbors.pivot_table(
            values='Score', index='Aminoacid', columns='Sequence_old'
        )
        aa_neighbors_pivoted.reset_index(drop=True, inplace=True)
        # Get the mean of the amino acids that appear in the aa_neighbors subset
        mean_neighbors = df_condensed_heatmap[list(aa_neighbors_pivoted.columns)]
        # Subtract average effect and do mean
        df_output[aa] = (aa_neighbors_pivoted - mean_neighbors).mean(axis=1)

    # Sort by aa
    df_output = df_output[neworder_aminoacids]
    # Sort in y axis desired order
    df_output = _sort_yaxis_aminoacids(df_output, neworder_aminoacids, aminoacids)
    return df_output


def _sort_yaxis_aminoacids(
    df_input: DataFrame, neworder_aminoacids: list, old_order: list
) -> DataFrame:
    # Sort in y axis desired order
    df_input['Aminoacid_new'] = old_order
    df_input['Aminoacid_new'] = Categorical(df_input['Aminoacid_new'], neworder_aminoacids)
    df_input.sort_values(by=['Aminoacid_new'], inplace=True)
    df_input.drop(['Aminoacid_new'], inplace=True, axis=1)
    return df_input
