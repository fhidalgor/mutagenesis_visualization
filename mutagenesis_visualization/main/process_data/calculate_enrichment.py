"""
This module contais the enrichment calculations from dna counts.
"""
from typing import Literal, Union, List, Optional, Tuple
from pathlib import Path
import numpy as np
from numpy import typing as npt
from pandas.core.frame import DataFrame
from scipy.stats import zscore
from mutagenesis_visualization.main.process_data.process_data_utils import (
    stopcodon_correction, array_to_df_enrichments, group_by_aa, filter_by_mad, nan_mode,
    kernel_correction, replace_inf
)

ZEROING_METHODS = Literal['none', 'zscore', 'counts', 'wt', 'wt synonymous', 'kernel', 'population']  # pylint: disable=invalid-name
ZEROING_METRICS = Literal['mean', 'mode', 'median']  # pylint: disable=invalid-name
STOPCODON: bool = True
MIN_COUNTS: int = 25
MIN_COUNTSWT: int = 100
STD_SCALE: float = 0.2
MAD_FILTERING: int = 2
MWT: float = 2
INFINITE: float = 3


def calculate_enrichment(
    aminoacids: Union[List[str], str],
    pre_lib: Union[str, DataFrame, npt.NDArray],
    post_lib: Union[str, DataFrame, npt.NDArray],
    pre_wt: Union[str, None, npt.NDArray] = None,
    post_wt: Union[str, None, npt.NDArray] = None,
    zeroing_method: ZEROING_METHODS = 'population',
    zeroing_metric: ZEROING_METRICS = 'median',
    stopcodon: bool = STOPCODON,
    min_counts: int = MIN_COUNTS,
    min_countswt: int = MIN_COUNTSWT,
    std_scale: Optional[float] = STD_SCALE,
    mad_filtering: int = MAD_FILTERING,
    mwt: float = MWT,
    infinite: float = INFINITE,
    output_file: Union[None, str, Path] = None
) -> npt.NDArray:
    """
    Determine the enrichment scores of a selection experiment, where there
    is a preselected population (input) and a selected population (output).

    Parameters
    -----------
    aminoacids : list, str
        Index of aminoacids (in order). Stop codon needs to be '*'.

    pre_lib : str, pandas dataframe or np.array
        Can be filepath and name of the exported txt file, dataframe or
        np.array.

    post_lib : str, pandas dataframe or np.array
        Can be filepath and name of the exported txt file, dataframe or
        np.array.

    pre_wt : str, or np.array, optional
        Str with filepath and name of the exported txt file or np.array.

    post_wt : str, or np.array, optional
        Str with filepath and name of the exported txt file or np.array.

    zeroing_method : str, default 'population'
        Method to normalize the data.
        Can also use 'none', 'zscore', 'counts', 'wt' or 'kernel'.
        If 'wt' is used 'pre_wt' must not be set to None.

    zeroing_metric : str, default 'median'
        Metric to zero the data. Only works if zeroing_method='population'
        or 'wt'. Can also be set to 'mean' or 'mode'.

    stopcodon : boolean, default False
        Use the enrichment score stop codons as a metric to determine
        the minimum enrichment score.

    min_counts : int, default 25
        If mutant has less than the min_counts, it will be replaced by
        np.nan.

    min_countswt : int, default 100
        If synonymous wild-type mutant has less than the min_counts, it
        will be replaced by np.nan.

    std_scale : float, default 0.2
        Factor by which the population is scaled. Set to None if you don't
        want to scale the data.

    mad_filtering : int, default 2
        Will apply MAD (median absolute deviation) filtering to data.

    mwt : int, default 2
        When MAD filtering is applied, mad_filtering is the number of
        medians away a data point must be to be discarded. mwt is only
        used when the population of wild-type alleles is the reference
        for data zeroing_method.

    infinite : int, default 3
        It will replace +infinite values with +3 and -infinite with -3.

    output_file : str, default None
        If you want to export the generated files, add the path and name.
        Example: 'path/filename.txt'. File will be save as a txt, csv, xlsx file.

    Returns
    --------
    zeroed : ndarray
        A np.array containing the enrichment scores.
    """

    # Convert to numpy if libraries are in dataframe format.
    # If input is a filepath, then load the txt files
    if isinstance(pre_lib, DataFrame):
        pre_lib = pre_lib.to_numpy()
    elif isinstance(pre_lib, str):
        pre_lib = np.loadtxt(pre_lib)
    if isinstance(post_lib, DataFrame):
        post_lib = post_lib.to_numpy()
    elif isinstance(post_lib, str):
        post_lib = np.loadtxt(post_lib)

    # Same thing for wt allele files
    if isinstance(pre_wt, str):
        pre_wt = np.loadtxt(pre_wt)
    if isinstance(post_wt, str):
        post_wt = np.loadtxt(post_wt)

    if isinstance(aminoacids, str):
        aminoacids = list(aminoacids)

    # Convert to df
    pre_lib_df: DataFrame = array_to_df_enrichments(pre_lib, aminoacids)
    post_lib_df: DataFrame = array_to_df_enrichments(post_lib, aminoacids)

    # Locate stop codons
    input_stopcodon: str = ''
    output_stopcodon: str = ''
    if stopcodon:
        input_stopcodon = pre_lib_df.loc['*'].astype(float)
        output_stopcodon = post_lib_df.loc['*'].astype(float)

    # Log10 of the counts for library and wt alleles
    log10_counts, output_lib_np = get_log_enrichment(
        pre_lib_df, post_lib_df, input_stopcodon, output_stopcodon, min_counts, stopcodon, infinite
    )
    # Group by amino acid
    log10_counts_grouped: DataFrame = group_by_aa(DataFrame(log10_counts), aminoacids)
    log10_counts_mad = np.ravel(np.array(log10_counts_grouped))
    if mad_filtering:
        log10_counts_mad = filter_by_mad(log10_counts_mad, mad_filtering)

    # Wt counts
    log10_wtcounts = None
    if pre_wt is not None:
        log10_wtcounts, _ = get_log_enrichment(
            pre_wt, post_wt, input_stopcodon, output_stopcodon, min_countswt, stopcodon, infinite
        )
        # MAD filtering
        # If set to m=1, if tosses out about 50% of the values. the mean barely changes though
        if mad_filtering:
            log10_wtcounts = filter_by_mad(log10_wtcounts, mwt)

    zeroed = zero_data(
        ratio_counts=np.log10(np.nansum(output_lib_np)/np.nansum(pre_lib)),
        log10_counts_grouped=log10_counts_grouped,
        log10_counts_mad=log10_counts_mad,
        log10_wtcounts=log10_wtcounts,
        zeroing_method=zeroing_method,
        zeroing_metric=zeroing_metric,
        std_scale=std_scale
    )

    if output_file:
        Path(output_file).parent.mkdir(parents=True, exist_ok=True)
        if Path(output_file).suffix == ".txt":
            np.savetxt(Path(output_file), zeroed, fmt='%i', delimiter='\t')
        elif Path(output_file).suffix == ".csv":
            df_output: DataFrame = DataFrame(zeroed)
            df_output.to_csv(Path(output_file), index=False, header=False)
        elif Path(output_file).suffix == ".xlsx":
            df_output = DataFrame(zeroed)
            df_output.to_excel(
                Path(output_file), sheet_name='Enrichment_Scores', index=False, header=False
            )

    return zeroed


def get_log_enrichment(
    input_lib: DataFrame,
    output_lib: DataFrame,
    input_stopcodon: DataFrame,
    output_stopcodon: DataFrame,
    min_counts: int,
    stopcodon: bool,
    infinite: float,
) -> Tuple[npt.NDArray, npt.NDArray]:
    """
    Calculate log10 enrichment scores from input and output counts.
    """
    # Copy data and replace low counts by np.nan
    input_lib_np: npt.NDArray = np.copy(input_lib.astype(float))
    output_lib_np: npt.NDArray = np.copy(output_lib.astype(float))
    input_lib_np[input_lib_np < min_counts] = np.nan

    # Stop codon correction
    if stopcodon:
        output_lib_np = stopcodon_correction(
            input_lib_np, output_lib_np, np.array(input_stopcodon), np.array(output_stopcodon)
        )

    # log10 of library and replace infinite values. This will potentially divide by zero.
    with np.errstate(divide='ignore'):
        counts_log10_ratio: npt.NDArray = replace_inf(
            np.log10(output_lib_np / input_lib_np), infinite
        )

    return counts_log10_ratio, output_lib_np


def zero_data(
    ratio_counts: float, log10_counts_grouped: DataFrame,
    log10_counts_mad: npt.NDArray, log10_wtcounts: Optional[npt.NDArray], zeroing_method: str,
    zeroing_metric: str, std_scale: float
) -> npt.NDArray:
    """
    Zeroes the data according to the input parameters.
    """
    # Statistics population using mad filtered data
    mean_pop = np.nanmean(log10_counts_mad)
    median_pop = np.nanmedian(log10_counts_mad)
    std_pop = np.nanstd(log10_counts_mad)
    mode_pop = nan_mode(log10_counts_mad)

    # Zero data, select case
    if zeroing_method == 'wt' or zeroing_method == 'wt synonymous' and log10_wtcounts is not None:
        mean_wt = np.nanmean(log10_wtcounts)
        median_wt = np.nanmedian(log10_wtcounts)
        std_wt = np.nanstd(log10_wtcounts)
        if len(log10_wtcounts) > 1:
            mode_wt = nan_mode(log10_wtcounts)
        else:  #case for only 1 wt
            mode_wt = log10_wtcounts
        if zeroing_metric == 'mean':
            zeroed = log10_counts_grouped - mean_wt
        elif zeroing_metric == 'median':
            zeroed = log10_counts_grouped - median_wt
        elif zeroing_metric == 'mode':
            zeroed = log10_counts_grouped - mode_wt
        if std_scale:
            zeroed = zeroed * std_scale / 2 / std_wt
    elif zeroing_method == 'population':
        if zeroing_metric == 'mean':
            zeroed = log10_counts_grouped - mean_pop
        elif zeroing_metric == 'median':
            zeroed = log10_counts_grouped - median_pop
        elif zeroing_metric == 'mode':
            zeroed = log10_counts_grouped - mode_pop
        if std_scale:
            zeroed = zeroed * std_scale / std_pop
    elif zeroing_method == 'counts':
        zeroed = log10_counts_grouped - ratio_counts
        if std_scale:
            zeroed = zeroed * std_scale / std_pop
    elif zeroing_method == 'kernel':
        zeroed_0, kernel_std = kernel_correction(log10_counts_grouped, cutoff=1)
        zeroed, kernel_std = kernel_correction(zeroed_0, cutoff=1)
        if std_scale:
            zeroed = zeroed * std_scale / kernel_std
    elif zeroing_method == 'zscore':
        zeroed = zscore(log10_counts_grouped, nan_policy='omit')
    elif zeroing_method == 'none':
        zeroed = log10_counts_grouped
    else:
        raise ValueError('Wrong zeroed parameter.')

    return np.array(zeroed)
