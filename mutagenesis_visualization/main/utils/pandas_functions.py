"""
This module contains utils to manipulate dataframes.
"""
import numpy as np
import pandas as pd
import itertools

def _transform_dataset(dataset, sequence, aminoacids, start_position, fillna):
    '''
    Internal function that constructs a dataframe from user inputs

    Parameters
    -----------
    dataset, sequence, aminoacids, start_position,fillna

    Returns
    --------
    Dataframe containing [Position, Sequence, Aminoacid, Variant, Score]
    '''

    # make a dataframe
    df = pd.DataFrame()

    # Define Columns
    df['Sequence'] = np.ravel([[aa] * len(aminoacids) for aa in sequence])

    # Create column with position label
    df['Position'] = np.ravel(
        [[i] * len(aminoacids)
         for i in range(start_position,
                        len(dataset[0]) + start_position)]
    )
    df['Aminoacid'] = aminoacids * len(dataset[0])
    df['Variant'] = df['Sequence'] + df['Position'].astype(str
                                                           ) + df['Aminoacid']
    df['Score'] = np.ravel(dataset.T)
    df['Score_NaN'] = np.ravel(dataset.T)

    # Eliminate NaNs
    df['Score'].fillna(fillna, inplace=True)

    # Eliminate stop codons
    df_clean = df[df['Aminoacid'] != '*'].copy()

    return df, df_clean


def _transform_sequence(dataset, sequence, start_position):
    '''
    Internal function that trims the input sequence

    Parameters
    -----------
    dataset, sequence, start_position

    Returns
    --------
    string containing trimmed sequence
    '''

    # truncate sequence
    trimmedsequence = sequence[start_position - 1:len(dataset[0]) +
                               start_position - 1]

    return trimmedsequence


def _transform_secondary(dataset, secondary, start_position, aminoacids):
    '''
    Internal function that trims the input secondary structure

    Parameters
    -----------
    dataset, sequence, start_position

    Returns
    --------
    list containing trimmed secondary structure (20 times each element)
    '''

    # Convert lists of lists to list
    secondary_list = list(itertools.chain.from_iterable(secondary))

    # Truncate list
    trimmedsecondary = secondary_list[start_position - 1:len(dataset[0]) +
                                      start_position - 1]

    # Multiply each element by number of aminoacids. not use stop codon
    aminoacids = list(np.copy(aminoacids))
    if '*' in aminoacids:
        aminoacids.remove('*')
    secondary_dup = [
        x for item in trimmedsecondary
        for x in itertools.repeat(item, len(aminoacids))
    ]

    return trimmedsecondary, secondary_dup

def df_rearrange(df, new_order, values='Score', show_snv=False):
    '''
    convert a df into a numpy array for mutagenesis data.
    Allows the option of keeping NaN scores

    Returns copy
    '''
    dfcopy = df.copy()

    # If only SNVs, turn rest to NaN
    if show_snv is True:
        dfcopy.loc[dfcopy['SNV?'] == False, values] = np.nan

    df_pivoted = dfcopy.pivot_table(
        values=values, index='Aminoacid', columns=['Position'], dropna=False
    )
    df_reindexed = df_pivoted.reindex(index=list(new_order))

    return df_reindexed


def _common(a, b):
    '''
    return common elements of two lists
    '''
    c = [value for value in a if value in b]
    return c


def _transpose(df, values='Score'):
    '''
    convert a df into a numpy array for mutagenesis data

    Returns copy
    '''
    df = df.pivot_table(
        values=values, index='Aminoacid', columns=['Position']
    ).T
    return df


def select_aa(df, selection, values='Score') -> pd.DataFrame:
    '''returns copy'''
    df = _transpose(df.copy(), values)
    return df[selection].T

def parse_pivot(
    df_imported, col_variant='variant', col_data='DMS', fill_value=np.nan
):
    '''
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
    '''

    # Copy
    df = df_imported.copy()

    # Extract position and amino acids that are being mutated
    df['Position'] = df[col_variant].str.extract(r'(\d+)').astype(int)
    df['Original'] = df[col_variant].str[0:1]
    df['Substitution'] = df[col_variant].str[-1:]

    # Get sequence
    sequence = list(
        df.groupby(
            by=['Position', 'Original'], as_index=False, group_keys=False
        ).sum()['Original']
    )

    # Pivot
    df_pivoted = df.pivot_table(
        index='Substitution',
        columns='Position',
        values=col_data,
        fill_value=fill_value,
        dropna=False
    )

    return df_pivoted, sequence
