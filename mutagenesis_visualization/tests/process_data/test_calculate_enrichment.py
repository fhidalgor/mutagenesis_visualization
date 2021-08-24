"""
This module tests the process data methods and functions.
"""
from typing import List
from itertools import product
from random import randint, random
import logging
import pandas as pd
from pandas.core.frame import DataFrame

from mutagenesis_visualization.main.process_data.calculate_enrichment import calculate_enrichment

log: logging.Logger = logging.getLogger('test_calculate_enrichment.py')


def test_calculate_enrichment() -> None:
    # Read counts from file (could be txt, csv, xlsx, etc...)
    prefix: str = "mutagenesis_visualization/"
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
    aminoacids_NNS: List[str] = list('AACDEFGGHIKLLLMNPPQRRRSSSTTVVWY*')

    # TODO: do 0 and then a random number from 100 - 1000
    stopcodon: List[str] = ["True", "False"]
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

    args_how_scale = product(zeroing_compatible, how, [True], std_scale, *common_args)
    args_how_no_scale = product(zeroing_compatible, how, [False], *common_args)
    args_no_how_scale = product(zeroing_other, [True], std_scale, *common_args)
    args_no_how_no_scale = product(zeroing_other, [False], *common_args)

    for args in args_how_scale:
        #print(args)
        zeroing, how, norm_std, std_scale, *common_args = args
        stopcodon, min_counts, mpop = common_args
        frequencies = calculate_enrichment(
            df_counts_pre.iloc[:, : 54],
            df_counts_sel.iloc[:, : 54],
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
            df_counts_pre.iloc[:, : 54],
            df_counts_sel.iloc[:, : 54],
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
            df_counts_pre.iloc[:, : 54],
            df_counts_sel.iloc[:, : 54],
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
            df_counts_pre.iloc[:, : 54],
            df_counts_sel.iloc[:, : 54],
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
