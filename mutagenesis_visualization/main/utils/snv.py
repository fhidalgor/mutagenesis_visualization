"""
This module contains functions that manipulate SNV-related variants.
"""
from typing import List, Dict
import pandas as pd
import itertools
from collections import defaultdict

def select_nonsnv(df_input: pd.DataFrame) -> pd.DataFrame:
    '''
    Generate a dataframe that contains the non-SNV variants and the enrichment score

    Parameters
    -----------
    df : pd.DataFrame

    Returns
    --------
    Dataframe containing a column of variants that are non-SNV, and the Score.
    '''
    # Dataframe with SNV
    df_snv: pd.DataFrame = select_snv(df_input)

    # Merge and eliminate duplicates. Keep Non-SNV
    df_nonsnv: pd.DataFrame = pd.concat([df_snv, df_input], sort=False)[[
        'Position', 'Variant', 'Score', 'Score_NaN'
    ]]
    df_nonsnv.drop_duplicates(subset='Variant', keep=False, inplace=True)

    return df_nonsnv


def select_snv(df: pd.DataFrame) -> pd.DataFrame:
    '''
    Select for SNV variants in DSM dataset

    Parameters
    -----------
    df : pandas dataframe containing DSM data

    Returns
    --------
    Modified dataframe('Variant','Score') where 'SNV?'== True. Returns copy
    '''

    # Use _add_SNV_boolean funciton
    df = add_snv_boolean(df.copy())

    # Select SNV? == True only
    df = df[df['SNV?'] == True].copy()

    # Select columns of interest
    df = df[['Position', 'Variant', 'Score', 'Score_NaN']].copy()

    # Reset index
    df.reset_index(drop=True, inplace=True)

    return df


def _aminoacids_snv(aa1: str, aa2: str, codontable, same_aa_SNV: bool=True) -> bool:
    '''
    Determine if two amino acids are snv (one base difference)

    Parameters
    -----------
    aa1 : str
    aa2 : str
    codontable : dict (did not want to generate each time I run the function)
    same_aa_SNV : boolean, default True
        If True, it will consider the same amino acid to be SNV of itself

    Returns
    --------
    boolean, True/False
    '''
    # Check if aa1 is aa2
    if not (same_aa_SNV) and (aa1.upper() == aa2.upper()):
        return False

    # Convert amino acids to codons
    codons1 = codontable[aa1.upper()]
    codons2 = codontable[aa2.upper()]

    # Generate a list of combination pairs between all codons in aa1 and aa2
    codon_combinations = list(itertools.product(codons1, codons2))

    # If one pair of combinations is a SNV, then return True
    for combination in codon_combinations:
        if _codons_pointmutants(combination[0], combination[1]) == True:
            return True
    return False


def add_snv_boolean(df: pd.DataFrame) -> pd.DataFrame:
    '''
    Add a column to dataframe indication if the variant is a SNV or not

    Parameters
    -----------
    df : pandas dataframe containing DSM data

    Returns
    --------
    Modified dataframe. Returns copy
    '''

    # Generate dictionary with aa and codon translation
    codontable: dict = _dict_codontoaa()

    # Add column with True/False input
    df['SNV?'] = df.apply(
        lambda x: _aminoacids_snv(x['Sequence'], x['Aminoacid'], codontable),
        axis=1
    )

    return df


def _codons_pointmutants(codon1: str, codon2: str, same_codon_SNV: bool=False) -> bool:
    '''
    Determine if two codons are SNV. Returns a boolean.
    If the codon is the same, will return False.
    Not case sensitive.

    Parameters
    -----------
    codon1 : str
    codon2 : str
    same_codon_SNV : boolean, default False
        If True, it will consider the same codon to be SNV of itself

    Returns
    --------
    boolean, True/False
    '''

    # Check if codons are the same
    if same_codon_SNV and codon1.upper() == codon2.upper():
        return True

    counter_occurrences = 0
    for index, base1 in enumerate(codon1.upper()):
        base2 = list(codon2.upper())[index]
        if base1 == base2:
            counter_occurrences = counter_occurrences + 1
    if counter_occurrences == 2:
        return True
    return False


def _are_pointmutants(aa: str, seqbase: str) -> bool:
    '''
    converts the amino acid to all possible degenerate codons and then checks if they are point mutants

    Parameters
    -----------
    aa: str
    seqbase: str

    Returns
    --------
    Boolean
    '''
    codontoaadict: dict = _dict_codontoaa()
    for codon in codontoaadict[aa]:
        if _codons_pointmutants(seqbase, codon):
            return True
    return False


def _are_pointmutants_list(aa: str, seqbase_list: List[str]) -> List[bool]:
    '''
    converts the amino acid to all possible degenerate codons and then checks if they are point mutants
    Same as _are_pointmutants but in list format

    Parameters
    -----------
    aa: str
    seqbase_list: list of str

    Returns
    --------
    List of Boolean
    '''
    pointmutants_list = []

    for seqbase in seqbase_list:
        pointmutants_list.append(_are_pointmutants(aa, seqbase))
    return pointmutants_list


def _dict_codontoaa() -> Dict[str, str]:
    '''
    Generates a dictionary with all amino acids and all possible codons.
    aa is the aminoacid of the mutation and seqbase is the original codon of the wtsequence
    '''
    bases = ['T', 'C', 'A', 'G']
    codons = [a + b + c for a in bases for b in bases for c in bases]
    aminoacids = list(
        'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
    )

    # dictionary with more than one value for each key
    codontoaadict = defaultdict(list)
    for codon, aminoacid in zip(codons, aminoacids):
        codontoaadict[aminoacid].append(codon)
    return codontoaadict


def _aatocodons(aminoacid):
    '''
    Inputs an aminoacid, returns all codons. Used dict_codontoaa()

    Parameters
    -----------
    aminoacid : str

    Returns
    --------
    List with all the codons that code for that amino acid
    '''

    # Dictionary with all codons and aa
    codontoaadict = _dict_codontoaa()

    # Codons for that amino acid
    codons = codontoaadict[aminoacid]

    return codons


def _aatocodons_df(df: pd.DataFrame, namecolumn: str) -> pd.DataFrame:
    '''
    Inputs a dataframe with a column of amino acids, returns all syn for each amino acidcodons.
    Used dict_codontoaa() and _aatocodons.

    Parameters
    -----------
    df : pandas dataframe
    namecolumn : str
        Name of the column containing the amino acids.

    Returns
    --------
    Dataframe with a column containing all the codons that code for that amino acid. Returns copy
    '''
    # Copy df
    df = df.copy()

    # Calculate each possible codon for every amino acid
    df['Codons_' + namecolumn] = df.apply(
        lambda x: _aatocodons(x[namecolumn]), axis=1
    )

    return df
