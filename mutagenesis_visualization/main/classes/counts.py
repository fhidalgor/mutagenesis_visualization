"""
This module contains the class Counts.
"""
from typing import List, Optional, Union
import numpy as np
from pandas.core.frame import DataFrame

from mutagenesis_visualization.main.utils.snv import is_dna, translate_codons
from mutagenesis_visualization.main.bar_graphs.library_representation import LibraryRepresentation
from mutagenesis_visualization.main.bar_graphs.mean_counts import MeanCounts


class Counts:
    """
    *Counts* represents the output of reading a fastq file.

    Parameters
    ----------
    dataframes : dataframe, list dataframes
        2D matrix containing the counts per codon.
        Columns will contain the amino acid substitutions, rows will
        contain the counts for each residue in the protein sequence.
        If multiple replicates, pass items in a list.

    start_position : int, default None
        First position in the protein sequence that will be used for the
        first column of the array. If a protein has been mutated only
        from residue 100-150, then if start_position = 100, the algorithm
        will trim the first 99 amino acids in the input sequence. The last
        residue will be calculated based on the length of the input array.
        We have set the default value to 2 because normally the Methionine
        in position 1 is not mutated.

    aminoacids : list, default None
        List of aminoacids (in order). Stop codon needs to be '*'.
        If none, it will use the index of the dataframe


    Methods
    --------
    mean_counts
    library_representation
    """
    def __init__(
        self,
        dataframes: Union[DataFrame, List[DataFrame]],
        start_position: Optional[int] = None,
        aminoacids: Optional[List[str]] = None,
    ):
        """
        Start.
        """
        if isinstance(dataframes, DataFrame):
            self.dataframes: List[DataFrame] = [dataframes]
        else:
            self.dataframes = dataframes
            df_mean: DataFrame = DataFrame(
                np.nanmean([np.array(df) for df in self.dataframes]),
                columns=self.dataframes[0].columns
            )
            self.dataframes.append(df_mean)

        self.positions: List[int] = list(self.dataframes[0].columns)
        self.start_position: int = self.positions[0]
        if start_position:
            self.start_position = start_position
            self.positions = np.arange(
                start_position,
                len(self.dataframes[0].columns) + start_position
            )

        # if aminoacids is none, use the index of the dataframe
        if aminoacids:
            self.aminoacids: List[str] = aminoacids
        elif is_dna(self.dataframes[0]):
            self.aminoacids = translate_codons(self.dataframes[0])
        else:
            self.aminoacids = list(self.dataframes[0].index)

        assert len(self.aminoacids) > 0

        self.mean_counts = MeanCounts(
            dataframes_raw=self.dataframes, positions=self.positions, aminoacids=self.aminoacids
        )
        self.library_representation = LibraryRepresentation(
            dataframes_raw=self.dataframes,
            aminoacids=self.aminoacids,
            positions=self.positions,
        )
