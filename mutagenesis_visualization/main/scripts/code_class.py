#!/usr/bin/env python
# coding: utf-8

# # Import Modules

# In[1]:


# Regular libraries
import numpy as np
import itertools

# Local imports
try: # Work with .py files
    from mutagenesis_visualization.main.scripts.code_kernel import (
        plot_kernel,
        plot_hist,
    )
    from mutagenesis_visualization.main.scripts.code_heatmaps import (
        plot_heatmap,
        plot_heatmap_rows,
        plot_heatmap_columns,
    )
    from mutagenesis_visualization.main.scripts.code_bar import (
        plot_mean,
        plot_meandifferential,
        plot_position,
        plot_meancounts,
        plot_library_representation,
    )
    from mutagenesis_visualization.main.scripts.code_scatter import (
        plot_scatter,
    )
    from mutagenesis_visualization.main.scripts.code_PCA import (
        plot_correlation,
        plot_individual_correlation,
        plot_group_correlation,
        plot_pca,
    )
    from mutagenesis_visualization.main.scripts.code_other import (
        plot_rank,
        plot_miniheatmap,
        plot_neighboreffect,
        plot_secondary,
        plot_roc,
        plot_cumulative,
    )

    from mutagenesis_visualization.main.scripts.code_plotly import (
        plot_mean_plotly,
        plot_heatmap_plotly,
        plot_histogram_plotly,
        plot_rank_plotly,
        plot_scatter_plotly,
        plot_scatter_3D_plotly,
        plot_scatter_3D_pdbprop_plotly,
    )
    from mutagenesis_visualization.main.scripts.code_utils import (
        _is_DNA, _translate_codons, _transform_sequence, _transform_dataset,
        _select_SNV, _select_nonSNV, _transform_secondary
    )
except ModuleNotFoundError: # Work on jupyter notebook
    import import_notebook
    from code_kernel import (plot_kernel, plot_hist)
    from code_heatmaps import (
        plot_heatmap,
        plot_heatmap_rows,
        plot_heatmap_columns,
    )
    from code_bar import (
        plot_mean,
        plot_meandifferential,
        plot_position,
        plot_meancounts,
        plot_library_representation,
    )
    from code_scatter import (plot_scatter)
    from code_PCA import (
        plot_correlation,
        plot_individual_correlation,
        plot_group_correlation,
        plot_pca,
    )
    from code_other import (
        plot_rank,
        plot_miniheatmap,
        plot_neighboreffect,
        plot_secondary,
        plot_roc,
        plot_cumulative,
    )

    from code_plotly import (
        plot_mean_plotly,
        plot_heatmap_plotly,
        plot_histogram_plotly,
        plot_rank_plotly,
        plot_scatter_plotly,
        plot_scatter_3D_plotly,
        plot_scatter_3D_pdbprop_plotly,
    )
    from code_utils import (
        _is_DNA,
        _translate_codons,
        _transform_sequence,
        _transform_dataset,
        _select_SNV,
        _select_nonSNV,
        _transform_secondary,
    )

try:
    from mutagenesis_visualization.main.scripts.code_pymol import (plot_pymol)
except ModuleNotFoundError:
    from code_pymol import (plot_pymol)
except ModuleNotFoundError:
    pass


# # Define Classes

# In[ ]:


class Counts:
    '''
    *Counts* represents the output of reading a fastq file.

    Parameters
    ----------
    df : pandas dataframe
        Dataframe containing the counts per codon.

    start_position : int, default None
        First position in the protein sequence that will be used for the first column of the
        array. If a protein has been mutated only from residue 100-150, then if start_position = 100,
        the algorithm will trim the first 99 amino acids in the input sequence. The last 
        residue will be calculated based on the length of the input array. We have set the default value to 2
        because normally the Methionine in position 1 is not mutated.

    aminoacids : list, default None
        List of aminoacids (in order). Stop codon needs to be '*'.
        If none, it will use the index of the dataframe   

    '''
    def __init__(self, df, start_position=None, aminoacids=None):
        self.dataframe = df

        if start_position:
            self.start_position = start_position
            self.positions = np.arange(start_position, len(df.columns) + 1)
        else:  # if none, use the columns of the dataframe
            self.positions = list(df.columns)
            self.start_position = self.positions[0]

        if aminoacids:
            self.aminoacids = aminoacids
        else:  # if aminoacids is none, use the index of the dataframe
            if _is_DNA(df):
                self.aminoacids = _translate_codons(df)
            else:
                self.aminoacids = list(df.index)
        # bar plots

    mean_counts = plot_meancounts
    library_representation = plot_library_representation


# In[ ]:


class Screen:
    '''
    *Screen* represents a mutagenesis experiment. If you are doing deep scan 
    mutagenesis, then every amino acid in the protein has been mutated to every possible 
    amino acid. For example, if there was a leucine at position 2, then this leucine would
    be mutated to the other 19 naturally occurring amino acids. However, you can also use the 
    package if you only have a handful of amino acid substitutions.


    Parameters
    -----------
    dataset : array
        2D matrix containing the enrichment scores of the point mutants. Columns will contain the
        amino acid substitutions, rows will contain the enrichment for each residue in the protein sequence.

    sequence : str
        Protein sequence (columns) in 1 letter code format.

    aminoacids : list, default list('ACDEFGHIKLMNPQRSTVWY*')
        Amino acid substitutions (rows). Submit in the same order that is used for the array.

    start_position : int, default 2
        First position in the protein sequence that will be used for the first column of the
        array. If a protein has been mutated only from residue 100-150, then if start_position = 100,
        the algorithm will trim the first 99 amino acids in the input sequence. The last 
        residue will be calculated based on the length of the input array. We have set the default value to 2
        because normally the Methionine in position 1 is not mutated.

    secondary : list, optional
        This parameter is used to group the data by secondary structure. The format is 
        the name of the secondary structure multiplied by the residue length of that motif.
        Example : [['β1']*(8),['L1']*(7),['α1']*(9),...,].

    roc_df: Pandas dataframe, optional
        A dataframe that contains a column of variants labeled 'Variant' with a column labeled 'Class'
        containing the true class of that mutation. This can be used to compare enrichment scores to some label (such as 
        pathogenicity as found in a Cancer database) using ROC AUC.

    fillna : int, default 0
        How to replace NaN values.


    Attributes
    ------------
    dataframe : pandas dataframe
        Contains the enrichment scores, position, sequence.

    Other attributes are same as input parameters: dataset, aminoacids, start_position, roc_df, secondary

    '''
    def __init__(
        self,
        dataset,
        sequence,
        aminoacids=list('ACDEFGHIKLMNPQRSTVWY*'),
        start_position=2,
        fillna=0,
        secondary=None,
        roc_df=None
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

    # Methods (Associated functions)

    # code_kernel
    kernel = plot_kernel
    histogram = plot_hist

    # code_heatmaps
    heatmap = plot_heatmap
    heatmap_rows = plot_heatmap_rows
    heatmap_columns = plot_heatmap_columns

    # code_bar
    mean = plot_mean
    differential = plot_meandifferential
    position = plot_position

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

