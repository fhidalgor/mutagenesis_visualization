#!/usr/bin/env python
# coding: utf-8

# In[17]:


from pathlib import Path
import pandas as pd
from pandas.testing import assert_frame_equal
import os
from itertools import product
from random import randint, random
import logging

log: logging.Logger = logging.getLogger('test_process_data')


try:
    from mutagenesis_visualization.main.scripts.code_process_data import (
        count_reads, calculate_enrichment
    )
except ModuleNotFoundError:
    import import_notebook
    import os
    directory = os.getcwd()
    new_directory = directory.replace('tests', 'main')
    os.chdir(new_directory)
    from code_process_data import count_reads, calculate_enrichment
    os.chdir(directory)


# In[3]:


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
        my_file = os.path.join(location, '../../data/for_tests', "short.fastq")
    except NameError:
        my_file = os.path.join('../../data/for_tests', "short.fastq")

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
        "GCC", "GCG", "TGC", "GAC", "GAG", "TTC", "GGC", "GGG", "CAC",
        "ATC", "AAG", "CTC", "CTG", "TTG", "ATG", "AAC", "CCC", "CCG",
        "CAG", "CGC", "CGG", "AGG", "TCC", "TCG", "AGC", "ACC", "ACG",
        "GTC", "GTG", "TGG", "TAC", "TAG",
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


# In[23]:


def test_calculate_enrichment():
    # Read counts from file (could be txt, csv, xlsx, etc...)
    prefix = "mutagenesis_visualization/"
    # prefix = "../../"
    df_counts_pre = pd.read_excel(
        prefix + 'data/hrasGAPGEF_counts.xlsx',
        'R1_before',
        skiprows=1,
        index_col='Codons',
        usecols='E:FN',
        nrows=32
    )

    df_counts_sel = pd.read_excel(
        prefix + 'data/hrasGAPGEF_counts.xlsx',
        'R1_after',
        skiprows=1,
        index_col='Codons',
        usecols='E:FN',
        nrows=32
    )

    # Ras parameters to create an object

    # Order of amino acids (from count_reads)
    aminoacids_NNS = list('AACDEFGGHIKLLLMNPPQRRRSSSTTVVWY*')

    # TODO: do 0 and then a random number from 100 - 1000
    stopcodon = ["True", "False"]
    min_counts = [0, randint(100, 1000)]
    mpop = [0.01, (random() + 0.01) * 10]  # mpop 0 causes an error
    common_args = [stopcodon, min_counts, mpop]

    # "wt" requires pre_wt
    zeroing_compatible = ["population"]
    # "zscore" zeroing broken at "elif zeroing == "zscore""
    zeroing_other = ["kernel", "counts"]
    how = ["median", "mean", "mode"]
    std_scale = [0.1, randint(0, 100)]

    log.info(f"{min_counts=}")
    log.info(f"{mpop=}")
    log.info(f"{std_scale=}")

    args_how_scale = product(
        zeroing_compatible, how, [True], std_scale, *common_args
    )
    args_how_no_scale = product(
        zeroing_compatible, how, [False], *common_args
    )
    args_no_how_scale = product(
        zeroing_other, [True], std_scale, *common_args
    )
    args_no_how_no_scale = product(
        zeroing_other, [False], *common_args
    )

    for args in args_how_scale:
        print(args)
        zeroing, how, norm_std, std_scale, *common_args = args
        stopcodon, min_counts, mpop = common_args
        frequencies = calculate_enrichment(
            df_counts_pre.iloc[:, :54],
            df_counts_sel.iloc[:, :54],

            zeroing=zeroing,
            how=how,
            norm_std=norm_std,
            std_scale=std_scale,

            aminoacids=aminoacids_NNS,
            stopcodon=stopcodon,
            min_counts=min_counts,
            min_countswt=100,
            mpop=mpop,
            mwt=2,
            infinite=3
        )
    for args in args_no_how_scale:
        print(args)
        zeroing, how, norm_std, *common_args = args
        stopcodon, min_counts, mpop = common_args
        frequencies = calculate_enrichment(
            df_counts_pre.iloc[:, :54],
            df_counts_sel.iloc[:, :54],

            zeroing=zeroing,
            how=how,
            norm_std=norm_std,

            aminoacids=aminoacids_NNS,
            stopcodon=stopcodon,
            min_counts=min_counts,
            min_countswt=100,
            mpop=mpop,
            mwt=2,
            infinite=3
        )
    for args in args_how_no_scale:
        print(args)
        zeroing, norm_std, std_scale, *common_args = args
        stopcodon, min_counts, mpop = common_args
        frequencies = calculate_enrichment(
            df_counts_pre.iloc[:, :54],
            df_counts_sel.iloc[:, :54],

            zeroing=zeroing,
            norm_std=norm_std,
            std_scale=std_scale,

            aminoacids=aminoacids_NNS,
            stopcodon=stopcodon,
            min_counts=min_counts,
            min_countswt=100,
            mpop=mpop,
            mwt=2,
            infinite=3
        )
    for args in args_how_scale:
        print(args)
        zeroing, norm_std, *common_args = args
        stopcodon, min_counts, mpop = common_args
        frequencies = calculate_enrichment(
            df_counts_pre.iloc[:, :54],
            df_counts_sel.iloc[:, :54],

            zeroing=zeroing,
            norm_std=norm_std,

            aminoacids=aminoacids_NNS,
            stopcodon=stopcodon,
            min_counts=min_counts,
            min_countswt=100,
            mpop=mpop,
            mwt=2,
            infinite=3
        )


# In[21]:


test_calculate_enrichment()


# In[ ]:




