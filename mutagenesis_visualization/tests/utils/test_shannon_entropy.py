"""
This module tests the utils.
"""
from typing import List, Dict
import pandas as pd
import numpy as np
from pandas.core.frame import DataFrame

#from mutagenesis_visualization.main.utils.shannon_entropy import (par)

def test_parseMSA() -> None:
    """
    Test of parse_msa.
    """
    # Load file
    try:
        location = os.path.dirname(os.path.realpath(__file__))
        my_file = os.path.join(location, '../../data/for_tests', "msa.fasta")
    except NameError:
        my_file = os.path.join('../../data/for_tests', "msa.fasta")
    # Read MSA
    msa, seq_lengths, index = _parseMSA(my_file, "fasta", 0)

    assert seq_lengths == [20]


def test_shannon_entropy_list_msa() -> None:
    # Load file
    try:
        location = os.path.dirname(os.path.realpath(__file__))
        my_file = os.path.join(location, '../../data/for_tests', "msa.fasta")
    except NameError:
        my_file = os.path.join('../../data/for_tests', "msa.fasta")

    # Read MSA
    msa, seq_lengths, index = _parseMSA(my_file, "fasta", 0)

    # Calculate entropy
    shannon = _shannon_entropy_list_msa(msa)

    assert _shannon_entropy_list_msa(msa) == [-0.0]*20
