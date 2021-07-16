"""
This module contains utils to manipulate dataframes.
"""
from typing import Any, Tuple
import numpy as np
import pandas as pd
import itertools


def transform_dataset(
    dataset,
    sequence: str,
    aminoacids: str,
    start_position: int,
    fillna: float,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Internal function that constructs a dataframe from user inputs
    Returns dataframe containing [Position, Sequence, Aminoacid, Variant, Score]
    """

    # make a dataframe
    df_input: pd.DataFrame = pd.DataFrame()

    # Define Columns
    df_input['Sequence'] = np.ravel([[aa] * len(aminoacids) for aa in sequence])

    # Create column with position label
    df_input['Position'] = np.ravel(
        [[i] * len(aminoacids) for i in range(start_position,
                                              len(dataset[0]) + start_position)]
    )
    df_input['Aminoacid'] = aminoacids * len(dataset[0])
    df_input['Variant'
             ] = df_input['Sequence'] + df_input['Position'].astype(str) + df_input['Aminoacid']
    df_input['Score'] = np.ravel(dataset.T)
    df_input['Score_NaN'] = np.ravel(dataset.T)

    # Eliminate NaNs
    df_input['Score'].fillna(fillna, inplace=True)

    # Eliminate stop codons
    df_clean = df_input[df_input['Aminoacid'] != '*'].copy()

    return df_input, df_clean


def transform_sequence(dataset: Any, sequence: str, start_position: int) -> str:
    """
    Internal function that trims the input sequence.
    """

    # truncate sequence
    return sequence[start_position - 1 : len(dataset[0]) + start_position - 1]


def transform_secondary(
    dataset: Any,
    secondary: list,
    start_position: int,
    aminoacids: str,
) -> Tuple[list, list]:
    """
    Internal function that trims the input secondary structure. Returns
    list containing trimmed secondary structure (20 times each element).
    """

    # Convert lists of lists to list
    secondary_list = list(itertools.chain.from_iterable(secondary))

    # Truncate list
    trimmedsecondary = secondary_list[start_position - 1 : len(dataset[0]) + start_position - 1]

    # Multiply each element by number of aminoacids. not use stop codon
    aminoacids = list(np.copy(aminoacids))
    if '*' in aminoacids:
        aminoacids.remove('*')
    secondary_dup = [
        x for item in trimmedsecondary for x in itertools.repeat(item, len(aminoacids))
    ]

    return trimmedsecondary, secondary_dup


def df_rearrange(
    df_input: pd.DataFrame,
    new_order: str,
    values: str = 'Score',
    show_snv: bool = False,
) -> pd.DataFrame:
    """
    Convert a df_input into a numpy array for mutagenesis data. Allows the
    option of keeping NaN scores. Returns copy
    """
    df_copy = df_input.copy()

    # If only SNVs, turn rest to NaN
    if show_snv is True:
        df_copy.loc[df_copy['SNV?'] == False, values] = np.nan

    df_pivoted = df_copy.pivot_table(
        values=values, index='Aminoacid', columns=['Position'], dropna=False
    )
    return df_pivoted.reindex(index=list(new_order))


def common_elements_list(a: list, b: list) -> list:
    """
    Return common elements of two lists.
    """
    return [value for value in a if value in b]


def _transpose(df_input: pd.DataFrame, values: str = 'Score') -> pd.DataFrame:
    """
    Convert a df_input into a numpy array for mutagenesis data. Returns copy
    """
    return df_input.pivot_table(values=values, index='Aminoacid', columns=['Position']).T


def select_aa(df_input: pd.DataFrame, selection, values: str = 'Score') -> pd.DataFrame:
    """returns copy"""
    df_input = _transpose(df_input.copy(), values)
    return df_input[selection].T


def parse_pivot(
    df_imported: pd.DataFrame, col_variant='variant', col_data='DMS', fill_value=np.nan
):
    """
    Parses a dataframe that contains saturation mutagenesis data in the Variant/Scores format.

    Parameters
    -----------
    df_imported : pandas dataframe
        Dataframe with the data imported with pd.read_excel.

    col_variant : str, default 'variant'
        Name of the column that contains the variants (ie T31A).

    col_data : str, default 'DMS'
        Name of the column that contains the saturation mutagenesis scores.

    fill_value : float, default np.nan
        What number to replace values that are omitted. It is possible that your
        dataset does not have a wt value.

    Returns
    --------
    df_pivoted : pandas dataframe
        Dataframe that has been pivoted. Values are the saturation mutagenesis data. Columns are
        the amino acid substitutions. Rows are the positions of the protein.

    sequence : list
        List of the amino acids that form the protein sequence.
    """

    # Copy
    df_input = df_imported.copy()

    # Extract position and amino acids that are being mutated
    df_input['Position'] = df_input[col_variant].str.extract(r'(\d+)').astype(int)
    df_input['Original'] = df_input[col_variant].str[0 : 1]
    df_input['Substitution'] = df_input[col_variant].str[-1 :]

    # Get sequence
    sequence = list(
        df_input.groupby(by=['Position', 'Original'], as_index=False,
                         group_keys=False).sum()['Original']
    )

    # Pivot
    df_pivoted = df_input.pivot_table(
        index='Substitution',
        columns='Position',
        values=col_data,
        fill_value=fill_value,
        dropna=False
    )

    return df_pivoted, sequence


def color_data(row, color_gof: str, color_lof: str) -> str:
    if row["Score"] > 0:
        return color_gof
    return color_lof


def process_mean_residue(
    dataframe_1: pd.DataFrame,
    dataframe_2: pd.DataFrame,
) -> pd.DataFrame:
    """
    Given two dataframes, it groups by position and truncates the
    longer one. It also drops nan values. Returns joined dataframe
    that contains the Scores, the position and the score of
    d1 - score of d2.
    """

    # truncate so both datasets have same length and delete stop codons
    dataset_1 = dataframe_1.groupby(['Position'], as_index=False).mean()
    dataset_2 = dataframe_2.groupby(['Position'], as_index=False).mean()
    min_length: int = min(len(dataset_1), len(dataset_2))

    # convert to dataframe and eliminate Nans
    df_ouput: pd.DataFrame = pd.DataFrame()
    df_ouput['dataset_1'] = list(dataset_1['Score'])[0 : min_length]
    df_ouput['dataset_2'] = list(dataset_2['Score'])[0 : min_length]
    df_ouput['Position'] = list(dataset_1['Position'])[0 : min_length]
    df_ouput['d1 - d2'] = df_ouput['dataset_1'] - df_ouput['dataset_2']
    df_ouput.dropna(how='any', inplace=True)
    return df_ouput


def process_by_pointmutant(
    dataframe_1: pd.DataFrame,
    dataframe_2: pd.DataFrame,
) -> pd.DataFrame:
    """
    Given two dataframes, it truncates the longer one. It also drops nan values.
    Returns joined dataframe that contains the Scores and the Variants.
    """
    # truncate so both datasets have same length and delete stop codons
    min_length: int = min(len(dataframe_1), len(dataframe_2))
    df_ouput: pd.DataFrame = pd.DataFrame()
    df_ouput['dataset_1'] = list(dataframe_1['Score_NaN'])[: min_length]
    df_ouput['dataset_2'] = list(dataframe_2['Score_NaN'])[: min_length]
    df_ouput['Variant'] = list(dataframe_1['Variant'])[: min_length]

    # eliminate Nans
    df_ouput.dropna(how='any', inplace=True)
    return df_ouput
