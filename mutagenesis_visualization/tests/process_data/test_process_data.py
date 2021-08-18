"""
This module tests the process data methods and functions.
"""
from pathlib import Path
import pandas as pd
from pandas.core.frame import DataFrame
from pandas.testing import assert_frame_equal
import os
from itertools import product
from random import randint, random
import logging
import tempfile
from collections import OrderedDict
import numpy as np

log: logging.Logger = logging.getLogger('test_process_data')


def test_assemble_sublibraries() -> None:
    # There aren't actually very many arguments to test here.  Once you remove
    # - all arguments that are just forwarded to calculate_enrichment
    # - the filename and excel sheet arguments
    # - treat the columns_wt arguments as either there or not
    # You're left with testing columns, nrows_pop, and columns_wt as a bool

    alphabet: str = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
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

def test_initialize_ordereddict() -> None:
    variants = _initialize_ordereddict(['ACC', 'CCT'])
    assert (isinstance(variants, OrderedDict))


def test_msa_enrichment() -> None:

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
