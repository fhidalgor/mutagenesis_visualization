"""
This module will handle the input of datasets (including replicates)
to the class Screen.
"""

from typing import List, Union
from dataclasses import dataclass
from copy import deepcopy
import numpy as np
from numpy import typing as npt
from pandas import DataFrame
from mutagenesis_visualization.main.utils.snv import select_snv, select_nonsnv


@dataclass
class DataframesHolder:
    """
    This dataclass will handle the massaging of input datasets.
    """
    df_notstopcodons: List[DataFrame]
    df_stopcodons: List[DataFrame]

    @property
    def df_snv(self) -> List[DataFrame]:
        """
        Select SNV variants only.
        """
        return [select_snv(df) for df in self.df_notstopcodons]

    @property
    def df_nonsnv(self) -> List[DataFrame]:
        """
        Select nonSNV variants only.
        """
        return [select_nonsnv(df) for df in self.df_notstopcodons]

    @property
    def df_wildtype_scores(self) -> List[DataFrame]:
        """
        Filter wild-type variants.
        """
        return [
            df.loc[df["Sequence"] == df["Aminoacid"]].drop_duplicates("Score_NaN")
            for df in self.df_notstopcodons
        ]

    @property
    def df_notstopcodons_mean(self) -> DataFrame:
        return [df.groupby(['Position'], as_index=False).mean() for df in self.df_notstopcodons]

    def df_notstopcodons_limit_score(self, min_score: float, max_score: float) -> List[DataFrame]:
        """
        Make values above a threshold to be the threshold.
        using df.where was not respecting the NaN values.
        """
        df_list: List[DataFrame] = deepcopy(self.df_notstopcodons)
        for df in df_list:
            df.loc[df["Score"] < min_score, "Score"] = min_score
            df.loc[df["Score"] > max_score, "Score"] = max_score
            df.loc[df["Score_NaN"] < min_score, "Score_NaN"] = min_score
            df.loc[df["Score_NaN"] > max_score, "Score_NaN"] = max_score
        return df_list

def handle_input_datasets(
    datasets: Union[npt.NDArray, DataFrame, List[Union[npt.NDArray, DataFrame]]]
) -> List[npt.NDArray]:
    """
    Transform the input into a list of np.arrays. The mean will be the
    last element (-1).
    """
    if isinstance(datasets, DataFrame):
        output_dataset: List[npt.NDArray] = [np.array(datasets)]
    elif isinstance(datasets, np.ndarray):
        output_dataset = [datasets]
    elif isinstance(datasets, list):
        output_dataset = [np.array(replicate) for replicate in datasets]
        output_dataset.append(np.nanmean(output_dataset, 0))
    return output_dataset
