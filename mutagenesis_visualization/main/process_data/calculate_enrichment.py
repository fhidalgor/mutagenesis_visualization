"""

"""
import numpy as np
import pandas as pd
import copy
from Bio import SeqIO
from Bio.Seq import Seq
from collections import OrderedDict
from os import path
from pathlib import Path
from typing import Union
from scipy import stats
from logomaker import alignment_to_matrix

def _array_to_df_enrichments(lib, aminoacids) -> pd.DataFrame:
    '''
    Aux function to transform array in df with index of amino acids.
    '''

    df = pd.DataFrame(index=aminoacids, data=lib)
    return df.astype(float)

def calculate_enrichment(
    pre_lib,
    post_lib,
    pre_wt=None,
    post_wt=None,
    aminoacids=list('AACDEFGGHIKLLLMNPPQRRRSSSTTVVWY*'),
    zeroing='population',
    how='median',
    norm_std=True,
    stopcodon=False,
    min_counts=25,
    min_countswt=100,
    std_scale=0.2,
    mpop=2,
    mad_filtering=True,
    mwt=2,
    infinite=3,
    output_file: Union[None, str, Path] = None
):
    """
    Determine the enrichment scores of a selection experiment, where there is a
    preselected population (input) and a selected population (output).

    Parameters
    -----------
    pre_lib : str, pandas dataframe or np.array
        Can be filepath and name of the exported txt file, dataframe or np.array.

    post_lib : str, pandas dataframe or np.array
        Can be filepath and name of the exported txt file, dataframe or np.array.

    pre_wt : str, or np.array, optional
        Str with filepath and name of the exported txt file or np.array.

    post_wt : str, or np.array, optional
        Str with filepath and name of the exported txt file or np.array.

    aminoacids : list, default ('AACDEFGGHIKLLLMNPPQRRRSSSTTVVWY*')
        Index of aminoacids (in order). Stop codon needs to be '*'.

    zeroing : str, default 'population'
        Method to normalize the data.
        Can also use 'none', 'zscore', 'counts', 'wt' or 'kernel'. If 'wt' is used 'pre_wt' must not be set to None.

    how : str, default 'median'
        Metric to zero the data. Only works if zeroing='population' or 'wt'.
        Can also be set to 'mean' or 'mode'.

    norm_std : boolean, default True
        If norm_std is set to True, it will scale the data.

    stopcodon : boolean, default False
        Use the enrichment score stop codons as a metric to determine the minimum enrichment score.

    min_counts : int, default 25
        If mutant has less than the min_counts, it will be replaced by np.nan.

    min_countswt : int, default 100
        If synonymous wild-type mutant has less than the min_counts, it will be replaced by np.nan.

    std_scale : float, default 0.2
        Factor by which the population is scaled. Only works if norm_std is set to True.

    mpop : int, default 2
        When using the median absolute deviation (MAD) filtering, mpop is the number of medians away
        a data point must be to be discarded.

    mad_filtering : boolean, default True
        Will apply MAD filtering to data.

    mwt : int, default 2
        When MAD filtering, mpop is the number of medians away a data point must be to
        be discarded. The difference with mpop is that mwt is only used when the population of wild-type
        alleles is the reference for data zeroing.

    infinite : int, default 3
        It will replace +infinite values with +3 and -infinite with -3.

    output_file : str, default None
        If you want to export the generated files, add the path and name of the file without suffix.
        Example: 'path/filename'. File will be save as a txt file.

    Returns
    --------
    zeroed : ndarray
        A np.array containing the enrichment scores.

    """

    # Convert to numpy if libraries are in dataframe format.
    # If input is a filepath, then load the txt files
    if type(pre_lib) is pd.DataFrame:
        pre_lib = pre_lib.to_numpy()
    elif type(pre_lib) is str:
        pre_lib = np.loadtxt(pre_lib)
    if type(post_lib) is pd.DataFrame:
        post_lib = post_lib.to_numpy()
    elif type(post_lib) is str:
        post_lib = np.loadtxt(post_lib)

    # Same thing for wt allele files
    if type(pre_wt) is str:
        pre_wt = np.loadtxt(pre_wt)
    if type(post_wt) is str:
        post_wt = np.loadtxt(post_wt)

    # Convert to df
    pre_lib = _array_to_df_enrichments(pre_lib, aminoacids)
    post_lib = _array_to_df_enrichments(post_lib, aminoacids)

    # Locate stop codons
    if stopcodon:
        input_stopcodon = pre_lib.loc['*'].astype(float)
        output_stopcodon = post_lib.loc['*'].astype(float)
    else:
        input_stopcodon = ''
        output_stopcodon = ''

    # Log10 of the counts for library and wt alleles
    log10_counts = _get_enrichment(
        pre_lib, post_lib, input_stopcodon, output_stopcodon, min_counts,
        stopcodon, infinite
    )
    # Group by amino acid
    df = pd.DataFrame(data=log10_counts)
    log10_counts_grouped = _group_byaa(df, aminoacids)

    # MAD filtering
    if mad_filtering:
        log10_counts_mad = _MAD_filtering(
            np.ravel(np.array(log10_counts_grouped)), mpop
        )
    else:
        log10_counts_mad = np.ravel(np.array(log10_counts_grouped))

    # Statistics population using mad filtered data
    mean_pop = np.nanmean(log10_counts_mad)
    median_pop = np.nanmedian(log10_counts_mad)
    std_pop = np.nanstd(log10_counts_mad)
    mode_pop = _nanmode(log10_counts_mad)

    # Wt counts
    if pre_wt is not None:
        log10_wtcounts = _get_enrichment(
            pre_wt, post_wt, input_stopcodon, output_stopcodon, min_countswt,
            stopcodon, infinite
        )
        # MAD filtering
        # If set to m=1, if tosses out about 50% of the values. the mean barely changes though
        if mad_filtering:
            log10_wtcounts = _MAD_filtering(log10_wtcounts, mwt)

        mean_wt = np.nanmean(log10_wtcounts)
        median_wt = np.nanmedian(log10_wtcounts)
        std_wt = np.nanstd(log10_wtcounts)
        if len(log10_wtcounts)>1:
            mode_wt = _nanmode(log10_wtcounts)
        else: #case for only 1 wt
            mode_wt = log10_wtcounts

    # Zero data, select case
    if zeroing == 'wt':
        if how == 'mean':
            zeroed = log10_counts_grouped - mean_wt
        elif how == 'median':
            zeroed = log10_counts_grouped - median_wt
        elif how == 'mode':
            zeroed = log10_counts_grouped - mode_wt
        elif norm_std == True:
            zeroed = zeroed * std_scale / 2 / std_wt
    elif zeroing == 'population':
        if how == 'mean':
            zeroed = log10_counts_grouped - mean_pop
        elif how == 'median':
            zeroed = log10_counts_grouped - median_pop
        elif how == 'mode':
            zeroed = log10_counts_grouped - mode_pop
        elif norm_std == True:
            zeroed = zeroed * std_scale / std_pop
    elif zeroing == 'counts':
        # Get the ratio of counts
        ratio = np.log10(post_lib.sum().sum()/pre_lib.sum().sum())
        zeroed = log10_counts_grouped + ratio
        if norm_std == True:
            zeroed = zeroed * std_scale / std_pop
    elif zeroing == 'kernel':
        zeroed_0, kernel_std = _kernel_correction(
            log10_counts_grouped, aminoacids
        )
        zeroed, kernel_std = _kernel_correction(zeroed_0, aminoacids, cutoff=1)
        if norm_std is True:
            zeroed = zeroed * std_scale / kernel_std
    elif zeroing == 'zscore':
        zeroed = stats.zscore(log10_counts_grouped, nan_policy='omit')
    elif zeroing == 'none':
        zeroed = log10_counts_grouped
    else:
        raise ValueError('Wrong zeroed parameter')
    # Export files
    if output_file:
        np.savetxt(Path(output_file), zeroed, fmt='%i', delimiter='\t')

    return zeroed
