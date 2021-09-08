"""
Use calculate enrichment in batch mode.
"""
from typing import Union, List, Optional
from pathlib import Path
from numpy import typing as npt
from numpy import savetxt, array, ravel, nanstd
from pandas.core.frame import DataFrame
from pandas import concat
from mutagenesis_visualization.main.process_data.calculate_enrichment import (
    calculate_enrichment, ZEROING_METRICS, ZEROING_METHODS, STD_SCALE, STOPCODON, MAD_FILTERING,
    MIN_COUNTS, MIN_COUNTSWT, MWT, INFINITE
)
from mutagenesis_visualization.main.process_data.process_data_utils import filter_by_mad


def batch_calculate_enrichment(
    aminoacids: Union[List[str], str],
    list_pre_lib: List[Union[str, DataFrame, npt.NDArray]],
    list_sel_lib: List[Union[str, DataFrame, npt.NDArray]],
    list_pre_wt: List[Union[str, None, npt.NDArray]] = None,
    list_sel_wt: List[Union[str, None, npt.NDArray]] = None,
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
) -> DataFrame:
    """
    Uses calculate_enrichment in batch to process multiple sublibraries,
    and concatenates them.

    Parameters
    -----------
    aminoacids : list
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
        Metric to zero the data. Only works if zeroing_method='population' or 'wt'.
        Can also be set to 'mean' or 'mode'.

    std_scale : float, default 0.2
        Factor by which the population is scaled. Set to None if you don't
        want to scale the data.

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
        Factor by which the population is scaled. Only works if norm_std
        is set to True.

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
    df_output : DataFrame
        A dataframe containing the enrichment scores assembled.    """

    enrichment_lib: List[DataFrame] = []

    for pre, sel, pre_wt, sel_wt in zip(list_pre_lib, list_sel_lib, list_pre_wt, list_sel_wt):
        enrichment_log10 = calculate_enrichment(
            aminoacids=aminoacids,
            pre_lib=pre,
            post_lib=sel,
            pre_wt=pre_wt,
            post_wt=sel_wt,
            zeroing_method=zeroing_method,
            zeroing_metric=zeroing_metric,
            stopcodon=stopcodon,
            min_counts=min_counts,
            min_countswt=min_countswt,
            std_scale=std_scale,
            mad_filtering=mad_filtering,
            mwt=mwt,
            infinite=infinite,
            output_file=None
        )
        enrichment_lib.append(DataFrame(enrichment_log10))

    # Concatenate sublibraries
    df_output: DataFrame = concat(enrichment_lib, ignore_index=True, axis=1)

    # Perform the std scaling
    log10_counts_mad = ravel(array(df_output))
    if mad_filtering:
        log10_counts_mad = filter_by_mad(log10_counts_mad, mad_filtering)
    std_pop = nanstd(log10_counts_mad)
    df_output = df_output * std_scale / std_pop

    if output_file:
        Path(output_file).parent.mkdir(parents=True, exist_ok=True)
        if Path(output_file).suffix == ".txt":
            savetxt(Path(output_file), df_output.to_numpy(), fmt='%i', delimiter='\t')
        elif Path(output_file).suffix == ".csv":
            df_output.to_csv(Path(output_file), index=False, header=False)
        elif Path(output_file).suffix == ".xlsx":
            df_output.to_excel(
                Path(output_file), sheet_name='Enrichment_Scores', index=False, header=False
            )
    return df_output
