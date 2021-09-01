"""
This module tests the process data methods and functions.
"""
from typing import List
from itertools import product
from random import randint, random
import logging
from pandas import read_excel

from mutagenesis_visualization.main.process_data.calculate_enrichment import calculate_enrichment
from mutagenesis_visualization.main.utils.data_paths import HRAS_GAPGEF_COUNTS

log: logging.Logger = logging.getLogger('test_calculate_enrichment.py')


def test_calculate_enrichment() -> None:
    """
    Test the calculate_enrichment function.
    """
    # Read counts from file (could be txt, csv, xlsx, etc...)
    df_counts_pre = read_excel(
        HRAS_GAPGEF_COUNTS, 'R1_before', skiprows=1, index_col='Codons', usecols='E:FN', nrows=32
    )

    df_counts_sel = read_excel(
        HRAS_GAPGEF_COUNTS, 'R1_after', skiprows=1, index_col='Codons', usecols='E:FN', nrows=32
    )

    # Ras parameters to create an object

    # Order of amino acids (from count_reads)
    aminoacids_nns: List[str] = list('AACDEFGGHIKLLLMNPPQRRRSSSTTVVWY*')

    # do 0 and then a random number from 100 - 1000
    stopcodon = ["True", "False"]
    min_counts = [0, randint(100, 1000)]
    mpop = [0.01, (random() + 0.01) * 10]  # mpop 0 causes an error
    common_args = [["True", "False"], min_counts, mpop]

    # "wt" requires pre_wt
    zeroing_compatible = ["population"]
    # "zscore" zeroing broken at "elif zeroing == "zscore""
    zeroing_other = ["kernel", "counts"]
    zeroing_metric = ["median", "mean", "mode"]
    std_scale = [0.1, randint(0, 100)]

    log.info(f"{min_counts=}")
    log.info(f"{mpop=}")
    log.info(f"{std_scale=}")

    args_how_scale = product(zeroing_compatible, zeroing_metric, std_scale, *common_args)
    args_how_no_scale = product(zeroing_compatible, zeroing_metric, *common_args)
    args_no_how_scale = product(zeroing_other, std_scale, *common_args)
    args_no_how_no_scale = product(zeroing_other, *common_args)

    for args in args_how_scale:
        zeroing_method, zeroing_metric, std_scale, *common_args = args
        stopcodon, min_counts, mad_filtering = common_args
        _ = calculate_enrichment(
            aminoacids=aminoacids_nns,
            pre_lib=df_counts_pre.iloc[:, : 54],
            post_lib=df_counts_sel.iloc[:, : 54],
            zeroing_method=zeroing_method,
            zeroing_metric=zeroing_metric,
            std_scale=std_scale,
            stopcodon=stopcodon,
            min_counts=min_counts,
            min_countswt=100,
            mad_filtering=mad_filtering,
            mwt=2,
            infinite=3
        )
    for args in args_no_how_scale:
        #print(args)
        zeroing_method, zeroing_metric, *common_args = args
        stopcodon, min_counts, mad_filtering = common_args
        _ = calculate_enrichment(
            aminoacids=aminoacids_nns,
            pre_lib=df_counts_pre.iloc[:, : 54],
            post_lib=df_counts_sel.iloc[:, : 54],
            zeroing_method=zeroing_method,
            zeroing_metric=zeroing_metric,
            stopcodon=stopcodon,
            min_counts=min_counts,
            min_countswt=100,
            mad_filtering=mad_filtering,
            mwt=2,
            infinite=3
        )
    for args in args_how_no_scale:
        #print(args)
        zeroing_method, zeroing_metric, *common_args = args
        stopcodon, min_counts, mad_filtering = common_args
        _ = calculate_enrichment(
            aminoacids=aminoacids_nns,
            pre_lib=df_counts_pre.iloc[:, : 54],
            post_lib=df_counts_sel.iloc[:, : 54],
            zeroing_method=zeroing_method,
            std_scale=0.2,
            stopcodon=stopcodon,
            min_counts=min_counts,
            min_countswt=100,
            mwt=2,
            infinite=3
        )
    for args in args_how_scale:
        #print(args)
        zeroing_method, zeroing_metric, std_scale, *common_args = args
        stopcodon, min_counts, mad_filtering = common_args
        _ = calculate_enrichment(
            aminoacids=aminoacids_nns,
            pre_lib=df_counts_pre.iloc[:, : 54],
            post_lib=df_counts_sel.iloc[:, : 54],
            zeroing_method=zeroing_method,
            stopcodon=stopcodon,
            min_counts=min_counts,
            min_countswt=100,
            mad_filtering=mad_filtering,
            mwt=2,
            infinite=3
        )
