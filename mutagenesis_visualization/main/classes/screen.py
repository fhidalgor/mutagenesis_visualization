"""
This module contains the class Screen, which groups the plotting classes.
"""
from typing import List, Optional, Union
import logging
from numpy import typing as npt
from numpy import delete
from pandas import DataFrame
from mutagenesis_visualization.main.bar_graphs.enrichment_bar import EnrichmentBar
from mutagenesis_visualization.main.bar_graphs.differential import Differential
from mutagenesis_visualization.main.bar_graphs.position_bar import PositionBar
from mutagenesis_visualization.main.bar_graphs.secondary import Secondary
from mutagenesis_visualization.main.kernel.kernel import Kernel
from mutagenesis_visualization.main.kernel.histogram import Histogram
from mutagenesis_visualization.main.kernel.multiple_kernels import MultipleKernel
from mutagenesis_visualization.main.kernel.sequence_differences import SequenceDifferences
from mutagenesis_visualization.main.heatmaps.heatmap import Heatmap
from mutagenesis_visualization.main.heatmaps.heatmap_columns import HeatmapColumns
from mutagenesis_visualization.main.heatmaps.heatmap_rows import HeatmapRows
from mutagenesis_visualization.main.heatmaps.miniheatmap import Miniheatmap
from mutagenesis_visualization.main.other_stats.rank import Rank
from mutagenesis_visualization.main.other_stats.cumulative import Cumulative
from mutagenesis_visualization.main.other_stats.roc_analysis import ROC
from mutagenesis_visualization.main.pca_analysis.correlation import Correlation
from mutagenesis_visualization.main.pca_analysis.individual_correlation import IndividualCorrelation
from mutagenesis_visualization.main.pca_analysis.pca import PCA
from mutagenesis_visualization.main.plotly.differential import DifferentialP
from mutagenesis_visualization.main.plotly.enrichment_bar import EnrichmentBarP
from mutagenesis_visualization.main.plotly.heatmap import HeatmapP
from mutagenesis_visualization.main.plotly.histogram import HistogramP
from mutagenesis_visualization.main.plotly.rank import RankP
from mutagenesis_visualization.main.plotly.scatter_3d_pdb import Scatter3DPDB
from mutagenesis_visualization.main.plotly.scatter_3d import Scatter3D
from mutagenesis_visualization.main.plotly.scatter import ScatterP
from mutagenesis_visualization.main.pymol.pymol import Pymol
from mutagenesis_visualization.main.scatter.scatter import Scatter
from mutagenesis_visualization.main.scatter.scatter_replicates import ScatterReplicates
from mutagenesis_visualization.main.utils.pandas_functions import (
    transform_dataset, transform_sequence, transform_secondary
)
from mutagenesis_visualization.main.utils.replicates_screen_input import (
    handle_input_datasets, DataframesHolder
)


class Screen:
    """
    *Screen* represents a mutagenesis experiment. If you are doing deep scan
    mutagenesis, then every amino acid in the protein has been mutated to
    every possible amino acid. For example, if there was a leucine at
    position 2, then this leucine would be mutated to the other 19 naturally
    occurring amino acids. However, you can also use the package if you
    only have a handful of amino acid substitutions.


    Parameters
    -----------
    datasets : array, list arrays, dataframe, list dataframes
        2D matrix containing the enrichment scores of the point mutants.
        Columns will contain the amino acid substitutions, rows will
        contain the enrichment for each residue in the protein sequence.
        If multiple replicates, pass items in a list.

    sequence : str
        Protein sequence in 1 letter code format.

    aminoacids : list, default list('ACDEFGHIKLMNPQRSTVWY*')
        Amino acid substitutions (rows). Submit in the same order that
        is used for the array.

    start_position : int, default 2
        First position in the protein sequence that will be used for the
        first column of the array. If a protein has been mutated only
        from residue 100-150, then if start_position = 100, the algorithm
        will trim the first 99 amino acids in the input sequence. The last
        residue will be calculated based on the length of the input array.
        We have set the default value to 2 because normally the Methionine
        in position 1 is not mutated.

    delete_position : int, List[int], default None
        Can delete positions (columns) in the dataset. For example, if you
        set start_position = 2 and delete_position = 122, you will be deleting
        the column 120 of the input dataset. The sequence parameter won't
        delete anything, so if you plan on deleting a few columns in your
        dataset, adjust the input sequence and secondary list.

    fillna : float, default 0
        How to replace NaN values.

    secondary : list, optional
        This parameter is used to group the data by secondary structure.
        The format is the name of the secondary structure multiplied by
        the residue length of that motif.
        Example : [['β1']*(8),['L1']*(7),['α1']*(9),...,].

    replicates: list[np.array, Dataframe], optional
        If you have multiple replicates for that experiment, pass them
        in the same format as dataset.

    Attributes
    ------------
    dataframe : pandas dataframe
        Contains the enrichment scores, position, sequence.

    Other attributes are same as input parameters: dataset, aminoacids,
    start_position, roc_df, secondary
    """
    def __init__(
        self,
        datasets: Union[npt.NDArray, DataFrame, List[Union[npt.NDArray, DataFrame]]],
        sequence: str,
        aminoacids: List[str],
        start_position: int = 2,
        delete_position: Union[int, List[int], None]= None,
        fillna: float = 0,
        secondary: Optional[List[List[str]]] = None,
    ):
        if isinstance(datasets,list):
            self.nreplicates = len(datasets)
        else:
            self.nreplicates = 1
        self.datasets: List[npt.NDArray] = handle_input_datasets(datasets)


        assert len(sequence) >= len(self.datasets[0][0]), "Input sequence is not long enough."

        if isinstance(aminoacids, list):
            self.aminoacids: List[str] = aminoacids
        elif isinstance(aminoacids, str):
            self.aminoacids = list(aminoacids)

        self.start_position: int = start_position
        self.end_position: int = len(self.datasets[0][0]) + start_position
        self.sequence_raw: str = sequence  # unchanged protein sequence
        self.sequence: str = transform_sequence(self.datasets[0], sequence, self.start_position)

        df_notstopcodons: List[DataFrame] = []
        df_stopcodons: List[DataFrame] = []
        for i, dataset in enumerate(self.datasets):
            if delete_position:
                self.datasets[i] = delete(dataset, delete_position, 1)

            df_output, df_not_stopcodons = transform_dataset(
                self.datasets[i], self.sequence, self.aminoacids, self.start_position, fillna
            )
            df_notstopcodons.append(df_not_stopcodons)
            df_stopcodons.append(df_output)
        self.dataframes: DataframesHolder = DataframesHolder(df_notstopcodons, df_stopcodons)

        # Optional parameters
        self.secondary_dup: Optional[List[str]] = None
        self.secondary: Optional[List[str]] = None
        if secondary:
            self.secondary, self.secondary_dup = transform_secondary(
                self.datasets[-1], secondary, self.start_position, self.aminoacids
            )

        # Only do this if there were replicates
        if len(self.datasets) > 1:
            self.scatter_replicates = ScatterReplicates(
                dataframes=self.dataframes,
                aminoacids=self.aminoacids,
            )

        # Assert messages

        # kernel
        self.kernel: Kernel = Kernel(
            dataframes=self.dataframes,
            aminoacids=self.aminoacids,
        )

        self.histogram: Histogram = Histogram(
            dataframes=self.dataframes,
            aminoacids=self.aminoacids,
        )
        self.multiple_kernel: MultipleKernel = MultipleKernel(
            dataframes=self.dataframes, aminoacids=self.aminoacids
        )
        self.sequence_differences: SequenceDifferences = SequenceDifferences(
            dataframes=self.dataframes, aminoacids=self.aminoacids
        )

        # heatmaps
        self.heatmap: Heatmap = Heatmap(
            dataframes=self.dataframes,
            sequence=self.sequence,
            start_position=self.start_position,
            secondary=self.secondary,
            aminoacids=self.aminoacids
        )

        self.heatmap_rows: HeatmapRows = HeatmapRows(
            dataframes=self.dataframes,
            sequence=self.sequence,
            start_position=self.start_position,
            aminoacids=self.aminoacids
        )
        self.heatmap_columns: HeatmapColumns = HeatmapColumns(
            dataframes=self.dataframes,
            sequence=self.sequence,
            start_position=self.start_position,
            aminoacids=self.aminoacids
        )
        self.miniheatmap: Miniheatmap = Miniheatmap(
            dataframes=self.dataframes,
            sequence_raw=self.sequence_raw,
            start_position=self.start_position,
            datasets=self.datasets,
            aminoacids=self.aminoacids
        )

        # bar
        self.enrichment_bar: EnrichmentBar = EnrichmentBar(
            dataframes=self.dataframes,
            start_position=self.start_position,
            aminoacids=self.aminoacids,
            secondary=self.secondary
        )
        self.differential: Differential = Differential(
            dataframes=self.dataframes,
            start_position=self.start_position,
            aminoacids=self.aminoacids,
            secondary=self.secondary
        )
        self.position_bar: PositionBar = PositionBar(
            dataframes=self.dataframes, aminoacids=self.aminoacids
        )
        self.secondary_mean: Secondary = Secondary(
            dataframes=self.dataframes,
            secondary_dup=self.secondary_dup,
            aminoacids=self.aminoacids
        )

        self.scatter: Scatter = Scatter(dataframes=self.dataframes, aminoacids=self.aminoacids)

        self.correlation: Correlation = Correlation(
            dataframes=self.dataframes,
            start_position=self.start_position,
            aminoacids=self.aminoacids
        )
        self.individual_correlation: IndividualCorrelation = IndividualCorrelation(
            dataframes=self.dataframes, aminoacids=self.aminoacids
        )
        self.pca: PCA = PCA(
            dataframes=self.dataframes,
            secondary_dup=self.secondary_dup,
            aminoacids=self.aminoacids
        )

        # other stats
        self.cumulative: Cumulative = Cumulative(
            dataframes=self.dataframes, aminoacids=self.aminoacids
        )
        self.rank: Rank = Rank(dataframes=self.dataframes, aminoacids=self.aminoacids)
        self.roc: ROC = ROC(dataframes=self.dataframes, aminoacids=self.aminoacids)

        # plotly
        self.plotly_differential: DifferentialP = DifferentialP(
            dataframes=self.dataframes,
            start_position=self.start_position,
            aminoacids=self.aminoacids
        )
        self.plotly_enrichment_bar: EnrichmentBarP = EnrichmentBarP(
            dataframes=self.dataframes, aminoacids=self.aminoacids
        )
        self.plotly_heatmap: HeatmapP = HeatmapP(
            sequence=self.sequence, dataframes=self.dataframes, aminoacids=self.aminoacids
        )
        self.plotly_histogram: HistogramP = HistogramP(
            dataframes=self.dataframes, aminoacids=self.aminoacids
        )
        self.plotly_rank: RankP = RankP(dataframes=self.dataframes, aminoacids=self.aminoacids)
        self.plotly_scatter: ScatterP = ScatterP(
            dataframes=self.dataframes, aminoacids=self.aminoacids
        )
        self.plotly_scatter_3d: Scatter3D = Scatter3D(
            dataframes=self.dataframes,
            start_position=self.start_position,
            end_position=self.end_position,
            aminoacids=self.aminoacids
        )
        self.plotly_scatter_3d_pdbprop: Scatter3DPDB = Scatter3DPDB(
            dataframes=self.dataframes,
            start_position=self.start_position,
            end_position=self.end_position,
            aminoacids=self.aminoacids
        )
        # pymol
        try:
            self.pymol: Pymol = Pymol(dataframes=self.dataframes)
        except Exception:  # pylint: disable=broad-except
            logging.info(
                "Pymol capability not loaded. Check the documentation for more information."
            )
