"""
This module contains the class Counts.
"""
from typing import List, Optional
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
    dataframe : pandas dataframe
        Dataframe containing the counts per codon.

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

    """
    def __init__(
        self,
        dataframe: DataFrame,
        start_position: Optional[int] = None,
        aminoacids: Optional[List[str]] = None,
    ):
        self.dataframe: DataFrame = dataframe
        self.positions: List[int] = list(self.dataframe.columns)
        self.start_position: int = self.positions[0]
        if start_position:
            self.start_position = start_position
            self.positions = np.arange(start_position, len(self.dataframe.columns) + 1)

        # if aminoacids is none, use the index of the dataframe
        if aminoacids:
            self.aminoacids: List[str] = aminoacids
        elif is_dna(self.dataframe):
            self.aminoacids = translate_codons(self.dataframe)
        else:
            self.aminoacids = list(self.dataframe.index)

        self.mean_counts = MeanCounts(dataframe=self.dataframe, positions=self.positions, aminoacids = self.aminoacids)
        self.library_representation = LibraryRepresentation(
            dataframe=self.dataframe,
            aminoacids=self.aminoacids,
            positions=self.positions,
        )
