"""
This module contains the class Screen, which groups the plotting classes.
"""
# Regular libraries
from typing import List, Optional
import numpy as np
import pandas as pd

from mutagenesis_visualization.main.kernel.kernel import Kernel

class Screen:
    '''
    *Screen* represents a mutagenesis experiment. If you are doing deep scan
    mutagenesis, then every amino acid in the protein has been mutated to
    every possible amino acid. For example, if there was a leucine at
    position 2, then this leucine would be mutated to the other 19 naturally
    occurring amino acids. However, you can also use the package if you
    only have a handful of amino acid substitutions.


    Parameters
    -----------
    dataset : array
        2D matrix containing the enrichment scores of the point mutants.
        Columns will contain the amino acid substitutions, rows will
        contain the enrichment for each residue in the protein sequence.

    sequence : str
        Protein sequence (columns) in 1 letter code format.

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

    secondary : list, optional
        This parameter is used to group the data by secondary structure.
        The format is the name of the secondary structure multiplied by
        the residue length of that motif.
        Example : [['β1']*(8),['L1']*(7),['α1']*(9),...,].

    roc_df: Pandas dataframe, optional
        A dataframe that contains a column of variants labeled 'Variant'
        with a column labeled 'Class' containing the true class of that
        mutation. This can be used to compare enrichment scores to some
        label (such as pathogenicity as found in a Cancer database)
        using ROC AUC.

    fillna : float, default 0
        How to replace NaN values.


    Attributes
    ------------
    dataframe : pandas dataframe
        Contains the enrichment scores, position, sequence.

    Other attributes are same as input parameters: dataset, aminoacids,
    start_position, roc_df, secondary

    '''
    def __init__(
        self,
        dataset,
        sequence: str,
        aminoacids: List[str]=list('ACDEFGHIKLMNPQRSTVWY*'),
        start_position: Optional[int]=2,
        fillna: Optional[float]=0,
        secondary: Optional[list]=None,
        roc_df: Optional[pd.DataFrame]=None
    ):
        # Instances
        self.dataset = np.array(dataset)
        self.aminoacids = aminoacids
        self.start_position = start_position
        self.end_position = len(self.dataset[0]) + start_position
        self.sequence_raw = ''.join(sequence)  # why I am doing this?
        self.sequence = _transform_sequence(
            self.dataset, self.sequence_raw, self.start_position
        )
        self.dataframe_stopcodons, self.dataframe = _transform_dataset(
            self.dataset, self.sequence, self.aminoacids, self.start_position,
            fillna
        )
        self.dataframe_SNV = _select_SNV(self.dataframe)
        self.dataframe_nonSNV = _select_nonSNV(self.dataframe)

        # Optional parameters
        self.roc_df = roc_df
        self.secondary = secondary
        if self.secondary is not None:
            self.secondary, self.secondary_dup = _transform_secondary(
                self.dataset, self.secondary, self.start_position,
                self.aminoacids
            )

        # Assert messages
        assert len(sequence) >= len(
            self.dataset[0]
        ), 'Input sequence is not long enough'


        # code_kernel
        self.kernel = Kernel(dataset = self.dataframe['Score_NaN'])
        self.histogram = Histogram()

        # code_heatmaps
        heatmap = plot_heatmap
        heatmap_rows = plot_heatmap_rows
        heatmap_columns = plot_heatmap_columns

        # code_bar
        mean = plot_mean
        differential = plot_meandifferential
        position = plot_position
        meancounts = plot_meancounts
        library_representation = plot_library_representation

        # code_scatter
        scatter = plot_scatter

        # code_PCA
        correlation = plot_correlation
        individual_correlation = plot_individual_correlation
        group_correlation = plot_group_correlation
        pca = plot_pca

        # code_Other
        rank = plot_rank
        miniheatmap = plot_miniheatmap
        neighboreffect = plot_neighboreffect
        secondary_mean = plot_secondary
        roc = plot_roc
        cumulative = plot_cumulative

        # code_plotly
        rank_plotly = plot_rank_plotly
        scatter_plotly = plot_scatter_plotly
        histogram_plotly = plot_histogram_plotly
        mean_plotly = plot_mean_plotly
        scatter_3D_plotly = plot_scatter_3D_plotly
        scatter_3D_pdbprop_plotly = plot_scatter_3D_pdbprop_plotly
        heatmap_plotly = plot_heatmap_plotly


        # pymol
        try:
            pymol = plot_pymol
        except:
            pass
