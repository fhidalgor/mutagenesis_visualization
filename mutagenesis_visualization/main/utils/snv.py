"""
This module contains functions that manipulate SNV-related variants.
"""
from itertools import product
from typing import List, Dict
from collections import defaultdict
from Bio.Seq import Seq
from pandas.core.frame import DataFrame
from pandas import concat


def select_nonsnv(df_input: DataFrame) -> DataFrame:
    """
    Generate a dataframe that contains the non-SNV variants and the enrichment score

    Parameters
    -----------
    df_input : DataFrame

    Returns
    --------
    Dataframe containing a column of variants that are non-SNV, and the Score.
    """
    # Dataframe with SNV
    df_snv: DataFrame = select_snv(df_input)

    # Merge and eliminate duplicates. Keep Non-SNV
    df_nonsnv: DataFrame = concat([df_snv, df_input],
                                  sort=False)[["Position", "Variant", "Score", "Score_NaN"]]
    df_nonsnv.drop_duplicates(subset="Variant", keep=False, inplace=True)

    return df_nonsnv


def select_snv(df_input: DataFrame) -> DataFrame:
    """
    Select for SNV variants in DSM dataset

    Parameters
    -----------
    df_input : pandas dataframe containing DSM data

    Returns
    --------
    Modified dataframe("Variant","Score") where "SNV?"== True. Returns copy
    """

    # Use _add_SNV_boolean funciton
    df_input = add_snv_boolean(df_input.copy())

    # Select SNV? == True only
    df_input = df_input[df_input["SNV?"] == True].copy()  # pylint: disable=singleton-comparison

    # Select columns of interest
    df_input = df_input[["Position", "Variant", "Score", "Score_NaN"]].copy()

    # Reset index
    df_input.reset_index(drop=True, inplace=True)

    return df_input


def _aminoacids_snv(
    aa1: str, aa2: str, codon_table: Dict[str, List[str]], same_aa_snv: bool = True
) -> bool:
    """
    Determine if two amino acids are snv (one base difference)

    Parameters
    -----------
    aa1 : str
    aa2 : str
    codon_table : dict (did not want to generate each time I run the function)
    same_aa_snv : boolean, default True
        If True, it will consider the same amino acid to be SNV of itself

    Returns
    --------
    boolean, True/False
    """
    # Check if aa1 is aa2
    if not (same_aa_snv) and (aa1.upper() == aa2.upper()):
        return False

    # Convert amino acids to codons
    codons1 = codon_table[aa1.upper()]
    codons2 = codon_table[aa2.upper()]

    # Generate a list of combination pairs between all codons in aa1 and aa2
    codon_combinations = list(product(codons1, codons2))

    # If one pair of combinations is a SNV, then return True
    for combination in codon_combinations:
        if _codons_pointmutants(combination[0], combination[1]) is True:
            return True
    return False


def add_snv_boolean(df_input: DataFrame, column_sequence: str = "Sequence", column_aminoacid: str = "Aminoacid") -> DataFrame:
    """
    Add a column to dataframe indication if the variant is a SNV or not

    Parameters
    -----------
    df_input : pandas dataframe containing DSM data.
    column_sequence: the column that contains the original aa.
    column_aminoacid: the column that contains the substitution.

    Returns
    --------
    Modified dataframe. Returns copy
    """

    # Generate dictionary with aa and codon translation
    codon_table: Dict[str, List[str]] = _dict_codon_to_aa()

    # Add column with True/False input
    df_input["SNV?"] = df_input.apply(
        lambda x: _aminoacids_snv(x[column_sequence], x[column_aminoacid], codon_table), axis=1
    )

    return df_input


def _codons_pointmutants(codon1: str, codon2: str, same_codon_snv: bool = False) -> bool:
    """
    Determine if two codons are SNV. Returns a boolean.
    If the codon is the same, will return False.
    Not case sensitive.

    Parameters
    -----------
    codon1 : str
    codon2 : str
    same_codon_snv : boolean, default False
        If True, it will consider the same codon to be SNV of itself

    Returns
    --------
    boolean, True/False
    """

    # Check if codons are the same
    if same_codon_snv and codon1.upper() == codon2.upper():
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
    """
    converts the amino acid to all possible degenerate codons and then checks if they are point mutants

    Parameters
    -----------
    aa: str
    seqbase: str

    Returns
    --------
    Boolean
    """
    codon_to_aa_dict: dict = _dict_codon_to_aa()
    for codon in codon_to_aa_dict[aa]:
        if _codons_pointmutants(seqbase, codon):
            return True
    return False


def _are_pointmutants_list(aa: str, seqbase_list: List[str]) -> List[bool]:
    """
    converts the amino acid to all possible degenerate codons and then checks if they are point mutants
    Same as _are_pointmutants but in list format

    Parameters
    -----------
    aa: str
    seqbase_list: list of str

    Returns
    --------
    List of Boolean
    """
    pointmutants_list = []

    for seqbase in seqbase_list:
        pointmutants_list.append(_are_pointmutants(aa, seqbase))
    return pointmutants_list


def _dict_codon_to_aa() -> Dict[str, List[str]]:
    """
    Generates a dictionary with all amino acids and all possible codons.
    aa is the aminoacid of the mutation and seqbase is the original codon of the wtsequence
    """
    bases: List[str] = ["T", "C", "A", "G"]
    codons: List[str] = [a + b + c for a in bases for b in bases for c in bases]
    aminoacids: List[str] = list("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")

    # dictionary with more than one value for each key
    codon_to_aa_dict: Dict[str, List[str]] = defaultdict(list)
    for codon, aminoacid in zip(codons, aminoacids):
        codon_to_aa_dict[aminoacid].append(codon)
    return codon_to_aa_dict


def _aa_to_codons(aminoacid: str) -> List[str]:
    """
    Inputs an aminoacid, returns all codons. Used dict_codon_to_aa()

    Parameters
    -----------
    aminoacid : str

    Returns
    --------
    List with all the codons that code for that amino acid
    """

    # Dictionary with all codons and aa
    codon_to_aa_dict: Dict[str, List[str]] = _dict_codon_to_aa()

    # Codons for that amino acid
    codons = codon_to_aa_dict[aminoacid]

    return codons


def _aa_to_codons_df(df_input: DataFrame, namecolumn: str) -> DataFrame:
    """
    Inputs a dataframe with a column of amino acids, returns all syn for each amino acidcodons.
    Used dict_codon_to_aa() and _aa_to_codons.

    Parameters
    -----------
    df_input : pandas dataframe
    namecolumn : str
        Name of the column containing the amino acids.

    Returns
    --------
    Dataframe with a column containing all the codons that code for that amino acid. Returns copy
    """
    # Copy df_input
    df_input = df_input.copy()

    # Calculate each possible codon for every amino acid
    df_input["Codons_" + namecolumn] = df_input.apply(
        lambda x: _aa_to_codons(x[namecolumn]), axis=1
    )

    return df_input


def translate_codons(df_input: DataFrame) -> List[str]:
    """
    Translate the index of the df_input from codons to AA.
    """
    return [str(Seq(codon).translate()) for codon in list(df_input.index)]


def is_dna(df_input: DataFrame) -> bool:
    """
    Check if the index of the dataframe are the DNA codons.
    """
    aminoacids: str = "DEFHIKLMNPQRSVWY*"
    for aa in aminoacids:
        if aa in "".join(list(df_input.index)):
            return False
    return True
