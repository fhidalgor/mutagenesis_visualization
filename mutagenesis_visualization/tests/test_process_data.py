"""
This module tests the process data methods and functions.
"""
from pathlib import Path
import pandas as pd
from pandas.testing import assert_frame_equal
import os
from itertools import product
from random import randint, random
import logging
import tempfile
from collections import OrderedDict
import numpy as np

log: logging.Logger = logging.getLogger('test_process_data')

try:
    from mutagenesis_visualization.main.scripts.code_process_data import (
        count_reads, calculate_enrichment, assemble_sublibraries,
        _initialize_ordereddict, msa_enrichment
    )
except ModuleNotFoundError:
    import import_notebook
    import os
    directory = os.getcwd()
    new_directory = directory.replace('tests', 'main')
    os.chdir(new_directory)
    from code_process_data import count_reads, calculate_enrichment, assemble_sublibraries, _initialize_ordereddict, msa_enrichment
    os.chdir(directory)


# In[ ]:


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

    # Now check with NNK
    # Create a temporary directory using the context manager
    with tempfile.TemporaryDirectory() as tmpdirname:
        cct_counts, cct_wt, info = count_reads(
            dna_sequence="cCt",
            input_file=my_file,
            codon_list='NNK',
            output_file=tmpdirname + '/counts.xlsx',
            full=True,
        )

    values_cct = [0] * 17 + [1] + [0] * 14
    index = [
        'GCG', 'GCT', 'TGT', 'GAT', 'GAG', 'TTT', 'GGG', 'GGT', 'CAT', 'ATT',
        'AAG', 'CTG', 'CTT', 'TTG', 'ATG', 'AAT', 'CCG', 'CCT', 'CAG', 'AGG',
        'CGG', 'CGT', 'AGT', 'TCG', 'TCT', 'ACG', 'ACT', 'GTG', 'GTT', 'TGG',
        'TAT', 'TAG'
    ]
    expected_cct_counts = pd.DataFrame(
        values_cct, index=index, columns=column_counts
    )
    expected_cct_wt = pd.DataFrame([[2, 'CCG', 'P', 0]],
                                   index=[16],
                                   columns=column_wt)
    mut_assert_df_equal(cct_counts, expected_cct_counts)
    mut_assert_df_equal(cct_wt, expected_cct_wt)


# In[ ]:


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
    args_how_no_scale = product(zeroing_compatible, how, [False], *common_args)
    args_no_how_scale = product(zeroing_other, [True], std_scale, *common_args)
    args_no_how_no_scale = product(zeroing_other, [False], *common_args)

    for args in args_how_scale:
        #print(args)
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
        #print(args)
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
        #print(args)
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
        #print(args)
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


# In[ ]:


def test_assemble_sublibraries():
    # There aren't actually very many arguments to test here.  Once you remove
    # - all arguments that are just forwarded to calculate_enrichment
    # - the filename and excel sheet arguments
    # - treat the columns_wt arguments as either there or not
    # You're left with testing columns, nrows_pop, and columns_wt as a bool

    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    pairs = [a + b for a, b in product(alphabet, alphabet)]
    all_columns = list(alphabet[5:]) + pairs[:144]
    col_lists = [partition_list(all_columns, i) for i in range(2, 6)]

    nrows_list = [randint(1, 31) for _ in range(3)] + [32]
    aminos = list(reversed("AACDEFGGHIKLLLMNPPQRRRSSSTTVVWY*"))
    amino_list = [list(reversed(aminos[:rows])) for rows in nrows_list]

    col_wt_list = [['A', 'B', 'C']]

    args = product(col_lists, col_wt_list, zip(nrows_list, amino_list))

    filename = "mutagenesis_visualization/data/hrasGAPGEF_counts.xlsx"
    # filename = "../../data/hrasGAPGEF_counts.xlsx"
    sheet_pre = "R1_before"
    sheet_post = "R1_after"
    columns_wt = ['A', 'B', 'C']
    nrows_pop = 32
    nrows_wt = [50, 37, 57]

    for columns, columns_wt, nrows_aminos in args:
        nrows_pop, aminos = nrows_aminos
        #print(f"{columns=}\t{columns_wt=}\t{nrows_pop=}")
        df = assemble_sublibraries(
            excel_path=filename,
            sheet_pre=sheet_pre,
            sheet_post=sheet_post,
            columns=columns,
            nrows_pop=nrows_pop,
            nrows_wt=nrows_wt,
            columns_wt=columns_wt,
            aminoacids=aminos,
            output_file=None
        )


def partition_list(array, num_partitions):
    """Partition array randomly where each partition has at least one item."""
    if num_partitions < 2:
        return [f"{array[0]}:{array[-1]}"]
    partition_idxs = []
    while len(partition_idxs) < num_partitions - 1:
        num = randint(0, len(array) - 1)
        if num not in partition_idxs:
            tmp_parts = partition_idxs.copy()
            tmp_parts.append(num)
            tmp_parts.sort()
            idx = tmp_parts.index(num)
            if idx != 0:
                if tmp_parts[idx] - tmp_parts[idx - 1] < 1:
                    continue
            if idx != len(tmp_parts) - 1:
                if tmp_parts[idx + 1] - tmp_parts[idx] < 1:
                    continue
            partition_idxs.append(num)
    partition_idxs.sort()
    parts = [f"{array[0]}:{array[partition_idxs[0] - 1]}"]
    for start, end in zip(partition_idxs, partition_idxs[1:]):
        parts.append(f"{array[start]}:{array[end - 1]}")
    parts.append(f"{array[partition_idxs[-1]]}:{array[-1]}")
    return parts


# In[ ]:


def test_initialize_ordereddict():
    variants = _initialize_ordereddict(['ACC', 'CCT'])
    assert (isinstance(variants, OrderedDict))


# # Shannon tests

# In[2]:


def test_msa_enrichment():

    # File location
    # Use relative file import to access the data folder
    try:
        location = os.path.dirname(os.path.realpath(__file__))
        my_file = os.path.join(location, '../../data/for_tests', "msa.fasta")
    except NameError:
        my_file = os.path.join('../../data/for_tests', "msa.fasta")

    # Create fake data
    class test_obj:
        dataframe = pd.DataFrame(
            index=[
                'Q', 'V', 'W', 'L', 'I', 'M', 'K', 'C', 'N', 'S', 'Y', 'D', 'F',
                'G', 'R', 'E', 'H', 'P', 'T', 'A'
            ]
        )
        dataframe['Position'] = np.arange(0,20)
        dataframe['Score'] = [0]*20
        dataframe['Aminoacid'] = list('ACDEFGHIKLMNPQRSTVWY')
    # Create test object
    test = test_obj()
    # Get df
    df, df_freq = msa_enrichment(test, my_file, 0)
    assert df['Shannon'].sum() == 0
