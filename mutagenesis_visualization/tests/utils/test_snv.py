"""
This module tests the utils.snv folder.
"""
from typing import List, Dict
from pandas.core.frame import DataFrame

from mutagenesis_visualization.main.utils.snv import (
    _aminoacids_snv, _dict_codon_to_aa, _codons_pointmutants, _are_pointmutants,
    _are_pointmutants_list, translate_codons, is_dna
)


def test_aminoacids_snv() -> None:
    """
    Testing full capabilities of function.
    """

    # Create dict with codons
    codon_table: Dict[str, List[str]] = _dict_codon_to_aa()

    # Create input arguments
    pairs: List[List[str]] = [['F', 'L'], ['I', 'M'], ['T', 'A'], ['S', 'R'], ['F', 'P'],
                              ['I', 'G'], ['T', 'L'], ['S', 'H'], ['A', 'A']]

    # Expected answer
    expected_answers: List[bool] = [True] * 4 + [False] * 5

    # Calculate answer and assert
    for pair, expected_answer in zip(pairs, expected_answers):
        calculated_answer = _aminoacids_snv(pair[0], pair[1], codon_table, same_aa_snv=False)
        assert expected_answer == calculated_answer, 'Error in determining snv'

    # Now change the same_aa_snv parameter
    assert _aminoacids_snv(
        'A', 'A', codon_table, same_aa_snv=True
    ) is True, 'Error in determining snv when the two amino acids are the same'


def test_codons_snv() -> None:
    """
    Testing full capabilities of function.
    """

    # Create input arguments
    pairs: List[List[str]] = [['AAA', 'AAT'], ['ACA', 'ACT'], ['ATA', 'TTA'], ['CCC', 'CCT'],
                              ['AAA', 'ACC'], ['CAA', 'CCC'], ['ATA', 'TTT'], ['CCC', 'AAA'],
                              ['AAA', 'AAA']]

    # Expected answer
    expected_answers: List[bool] = [True] * 4 + [False] * 5

    # Calculate answer and assert
    for pair, expected_answer in zip(pairs, expected_answers):
        calculated_answer = _codons_pointmutants(pair[0], pair[1], same_codon_snv=False)
        assert expected_answer == calculated_answer, 'Error in determining snv'

    # Now change the same_aa_snv parameter
    assert _codons_pointmutants( # pylint: disable=singleton-comparison
        'CAA', 'CAA', same_codon_snv=True
    ) == True, 'Error in determining snv when the two codons are the same'


def test_are_pointmutants() -> None:
    """Test."""
    assert _are_pointmutants('M', 'ATG') is False
    assert _are_pointmutants('A', 'ATG') is False
    assert _are_pointmutants('F', 'TTA') is True


def test_are_pointmutants_list() -> None:
    """Test."""
    assert _are_pointmutants_list('A', ['ATG', 'ATT', 'ACT']) == [False, False, True]


def test_translate_codons() -> None:
    """
    Testing full capabilities of function.
    """
    list_codons: List[str] = [
        'AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 'AGC', 'AGG', 'AGT', 'ATA',
        'ATC', 'ATG', 'ATT', 'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC',
        'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT', 'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG',
        'GCT', 'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT', 'TAA', 'TAC', 'TAG', 'TAT',
        'TCA', 'TCC', 'TCG', 'TCT', 'TGA', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT'
    ]
    list_aminoacids: List[str] = [
        'K', 'N', 'K', 'N', 'T', 'T', 'T', 'T', 'R', 'S', 'R', 'S', 'I', 'I', 'M', 'I', 'Q', 'H',
        'Q', 'H', 'P', 'P', 'P', 'P', 'R', 'R', 'R', 'R', 'L', 'L', 'L', 'L', 'E', 'D', 'E', 'D',
        'A', 'A', 'A', 'A', 'G', 'G', 'G', 'G', 'V', 'V', 'V', 'V', '*', 'Y', '*', 'Y', 'S', 'S',
        'S', 'S', '*', 'C', 'W', 'C', 'L', 'F', 'L', 'F'
    ]
    df_input: DataFrame = DataFrame(index=list_codons)
    assert (
        translate_codons(df_input) == list_aminoacids
    ), 'error when translating the codons of the dataframe index'


def test_is_dna() -> None:
    """
    Testing full capabilities of function.
    """
    df_1: DataFrame = DataFrame(index=['A', 'C', 'T', 'G', 'P', 'L'])
    df_2: DataFrame = DataFrame(index=['ATA', 'CAT', 'TGG', 'TGT'])
    assert (is_dna(df_1) is False), 'error determining if the index of the dataframe contains DNA'
    assert (is_dna(df_2) is True), 'error determining if the index of the dataframe contains DNA'
