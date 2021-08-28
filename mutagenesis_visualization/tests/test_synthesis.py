"""
Test dna synthesis classes and functions.
"""

from typing import List, Tuple
from pandas.core.frame import DataFrame
from mutagenesis_visualization.main.classes.create_variants import CreateVariants
from mutagenesis_visualization.main.classes.generate_primers import (
    GeneratePrimers, _primer_design, _create_primers_list
)


def test_generate_primers() -> None:
    """
    Test the class GeneratePrimers.
    """
    dna: str = 'GGCAATGCGcccccaATGaaaaaaTAAaaACGGGGTTTTaaa'
    start: str = 'GGCAATGCGccccca'
    end: str = 'aaACGGGGTTTTaaa'
    generate_primers: GeneratePrimers = GeneratePrimers(dna, start, end)

    #create a temporary directory using the context manager
    df_primers: DataFrame = generate_primers(codon='NNS', length_primer=18, melting_temp=60)

    assert len(df_primers) == 9, 'error the primer length does not match'
    assert df_primers.iloc[
        1, 1] == 'NNSGCGCCCCCAATGAAAAAATAAAAACGGG', 'error the first primer is not mutated at Met'
    assert df_primers.iloc[
        8, 1
    ] == 'CGCCCCCAATGAAAAAANNSAAACGGGGTTTTAA', 'error the last primer is not mutated at the stop codon'


def test_primer_design() -> None:
    """
    Test the primer design function.
    """
    dna: str = 'GGCAATGCGcccccaATGaaaaaaTAAaaACGGGGTTTTaaa'  # 42 codons
    primers0: Tuple[str, str] = _primer_design(
        dna, codon='NNS', codon_position=0, length_primer=None, melting_temp=60
    )
    primers1: Tuple[str, str] = _primer_design(
        dna, codon='NNS', codon_position=1, length_primer=None, melting_temp=60
    )
    primers2: Tuple[str, str] = _primer_design(
        dna, codon='NNS', codon_position=2, length_primer=None, melting_temp=60
    )
    primers3: Tuple[str, str] = _primer_design(
        dna, codon='NNS', codon_position=3, length_primer=None, melting_temp=60
    )
    primers21: Tuple[str, str] = _primer_design(
        dna, codon='NNS', codon_position=21, length_primer=18, melting_temp=None
    )
    assert isinstance(primers3, tuple), 'error the output is not a tuple'

    assert primers0 == (
        'NNSAATGCGcccccaATGaaaaaaTAAaaACGG', 'CCGTttTTAttttttCATtgggggCGCATTSNN'
    ), 'error the output is not the forward and reverse primers at codon position 0'
    assert primers1 == (
        'NNSATGCGcccccaATGaaaaaaTAAaaACGGG', 'CCCGTttTTAttttttCATtgggggCGCATSNN'
    ), 'error the output is not the forward and reverse primers at codon position 1'
    assert primers2 == (
        'NNSTGCGcccccaATGaaaaaaTAAaaACGGG', 'CCCGTttTTAttttttCATtgggggCGCASNN'
    ), 'error the output is not the forward and reverse primers at codon position 2'
    assert primers3 == (
        'NNSGCGcccccaATGaaaaaaTAAaaACGGG', 'CCCGTttTTAttttttCATtgggggCGCSNN'
    ), 'error the output is not the forward and reverse primers at codon position 3'

    assert primers21 == (
        'AATGCGcccccaATGaaaNNSTAAaaACGGGGTTTTaaa', 'tttAAAACCCCGTttTTASNNtttCATtgggggCGCATT'
    ), 'error the codon is not replaced after position 21'


def test_create_primers_list() -> None:
    """
    Test the function _create_primers_list
    """
    dna: str = 'GGCAATGCGcccccaATGaaaaaaTAAaaACGGGGTTTTaaa'

    primers_list: Tuple[List[str], List[str]] = _create_primers_list(
        dna, start_codon=3, end_codon=6, codon='NNS', length_primer=15, melting_temp=60
    )
    assert isinstance(primers_list[0], list), 'error the type is not a list'


def test_create_variants() -> None:
    """
    Test the create variants function.
    """
    dna: str = 'AAGAAGAAG'
    codon_list1: List[str] = ['GAA', 'TAA', 'AAA']
    create_variants: CreateVariants = CreateVariants()
    variants1: DataFrame = create_variants(dna, codon_list1)
    assert variants1.shape == (10, 1), 'error there are an incorrect number of variants'
    assert variants1.iloc[0, 0] == 'AAGAAGAAG', 'error the first output is not wild type'
    assert variants1.iloc[
        1, 0] == 'GAAAAGAAG', 'error the second output is not replaced with the first codon'
    assert variants1.iloc[
        2, 0] == 'TAAAAGAAG', 'error the third output is not replaced with the second codon'

    #create a temporary directory using the context manager
    codon_list2: List[str] = ['CAA', 'CCC', 'TTT', 'AAA']
    variants2: DataFrame = create_variants(dna, codon_list2)

    assert variants2.iloc[
        1, 0] == 'CAAAAGAAG', 'error the second output is not replaced with the first codon'
    assert variants2.shape == (13, 1), 'error there are an incorrect number of variants'
