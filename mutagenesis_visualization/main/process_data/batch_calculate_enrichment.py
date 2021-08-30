"""
Use calculate enrichment in batch mode.
"""
from typing import Union, List, Literal
from pathlib import Path
from numpy import typing as npt
from numpy import savetxt
from pandas.core.frame import DataFrame
from pandas import concat
from mutagenesis_visualization.main.process_data.calculate_enrichment import calculate_enrichment


def batch_calculate_enrichment(
    aminoacids: Union[List[str], str],
    list_pre_lib: List[Union[str, DataFrame, npt.NDArray]],
    list_sel_lib: List[Union[str, DataFrame, npt.NDArray]],
    list_pre_wt: List[Union[str, None, npt.NDArray]] = None,
    list_sel_wt: List[Union[str, None, npt.NDArray]] = None,
    zeroing_method: Literal['none', 'zscore', 'counts', 'wt', 'kernel',
                            'population'] = 'population',
    zeroing_metric: Literal['mean', 'mean', 'mode', 'median'] = 'median',
    norm_std: bool = True,
    stopcodon: bool = False,
    min_counts: int = 25,
    min_countswt: int = 100,
    std_scale: float = 0.2,
    mpop: float = 2,
    mad_filtering: bool = True,
    mwt: float = 2,
    infinite: float = 3,
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

    norm_std : boolean, default True
        If norm_std is set to True, it will scale the data.

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

    mpop : int, default 2
        When using the median absolute deviation (MAD) filtering, mpop is
        the number of medians away a data point must be to be discarded.

    mad_filtering : boolean, default True
        Will apply MAD filtering to data.

    mwt : int, default 2
        When MAD filtering, mpop is the number of medians away a data
        point must be to be discarded. The difference with mpop is that
        mwt is only used when the population of wild-type alleles is
        the reference for data zeroing_method.

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
            aminoacids, pre, sel, pre_wt, sel_wt, zeroing_method, zeroing_metric, norm_std,
            stopcodon, min_counts, min_countswt, std_scale, mpop, mad_filtering, mwt, infinite, None
        )
        enrichment_lib.append(DataFrame(enrichment_log10))

    # Concatenate sublibraries
    df_output: DataFrame = concat(enrichment_lib, ignore_index=True, axis=1)

    if output_file:
        if Path(output_file).suffix == ".txt":
            savetxt(Path(output_file), df_output.to_numpy(), fmt='%i', delimiter='\t')
        elif Path(output_file).suffix == ".csv":
            df_output.to_csv(Path(output_file), index=False, header=False)
        elif Path(output_file).suffix == ".xlsx":
            df_output.to_excel(
                Path(output_file), sheet_name='Enrichment_Scores', index=False, header=False
            )
    return df_output
