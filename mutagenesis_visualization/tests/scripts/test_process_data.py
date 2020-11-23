#!/usr/bin/env python
# coding: utf-8

# In[1]:


try:
    from mutagenesis_visualization.main.scripts.code_process_data import (
        count_reads
    )
except ModuleNotFoundError:
    import import_notebook
    import os
    directory = os.getcwd()
    new_directory = directory.replace('tests', 'main')
    os.chdir(new_directory)
    from code_process_data import (count_reads)
    os.chdir(directory)

from pathlib import Path
import pandas as pd
from pandas.testing import assert_frame_equal
import os


# In[10]:


"""
To run these tests run pytest from the root directory. "Test CCT" is the one
that fails.

Other bits and bobs:
    In _enumerate_variants you don't use firstwtseq so you can remove that. And
    if that goes you can also not pass dna_sequence.
"""


def test_count_reads():
    # Need to remove firstwtseq from _enumerate_variants

    def mut_assert_df_equal(left: pd.DataFrame, right: pd.DataFrame):
        assert_frame_equal(
            left,
            right,
            check_dtype=False,
            check_index_type=False,
            check_column_type=False,
            check_frame_type=False,
            check_names=False,
            check_like=True,
        )

    # File location
    # Use relative file import to access the data folder
    try:
        location = os.path.dirname(os.path.realpath(__file__))
        my_file = os.path.join(location, '../data/for_tests', "short.fastq")
    except NameError:
        my_file = os.path.join('../data/for_tests', "short.fastq")

    # Create dataframe
    codon_list = [
        'AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 'AGC',
        'AGG', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT', 'CAA', 'CAC', 'CAG', 'CAT',
        'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC', 'CGG', 'CGT', 'CTA', 'CTC',
        'CTG', 'CTT', 'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT',
        'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT', 'TAA', 'TAC',
        'TAG', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT', 'TGA', 'TGC', 'TGG', 'TGT',
        'TTA', 'TTC', 'TTG', 'TTT'
    ]
    index = pd.Index(codon_list)
    column_counts = pd.Index([2])
    column_wt = pd.Index(["Position", "Codon", "Aminoacid", "Counts"])
    values_cct = values_atc = [0] * 23 + [1] + [0] * 40

    # Test ATG
    expected_atc_counts = pd.DataFrame(
        values_atc, index=index, columns=column_counts
    )
    expected_atc_wt = pd.DataFrame([], columns=column_wt)
    atg_counts, atg_wt = count_reads("atg", my_file, codon_list)
    # return atg_counts
    mut_assert_df_equal(atg_counts, expected_atc_counts)
    mut_assert_df_equal(atg_wt, expected_atc_wt)

    # Test CCT
    expected_cct_counts = pd.DataFrame(
        values_cct, index=index, columns=column_counts
    )
    expected_cct_wt = pd.DataFrame([[2, 'CCA', 'P', 0], [2, 'CCC', 'P', 0],
                                    [2, 'CCG', 'P', 0]],
                                   index=[20, 21, 22],
                                   columns=column_wt)
    cct_counts, cct_wt = count_reads("cCt", my_file, codon_list)
    mut_assert_df_equal(cct_counts, expected_cct_counts)
    mut_assert_df_equal(cct_wt, expected_cct_wt)

    # Test CCT when not in codon list
    index = pd.Index([
        "GCC",
        "GCG",
        "TGC",
        "GAC",
        "GAG",
        "TTC",
        "GGC",
        "GGG",
        "CAC",
        "ATC",
        "AAG",
        "CTC",
        "CTG",
        "TTG",
        "ATG",
        "AAC",
        "CCC",
        "CCG",
        "CAG",
        "CGC",
        "CGG",
        "AGG",
        "TCC",
        "TCG",
        "AGC",
        "ACC",
        "ACG",
        "GTC",
        "GTG",
        "TGG",
        "TAC",
        "TAG",
    ])
    values_cct = [0] * 32
    expected_cct_counts = pd.DataFrame(
        values_cct, index=index, columns=column_counts
    )
    expected_cct_wt = pd.DataFrame([[2, 'CCC', 'P', 0], [2, 'CCG', 'P', 0]],
                                   index=[16, 17],
                                   columns=column_wt)
    cct_counts, cct_wt = count_reads("cCt", my_file, 'NNS')
    mut_assert_df_equal(cct_counts, expected_cct_counts)
    mut_assert_df_equal(cct_wt, expected_cct_wt)

