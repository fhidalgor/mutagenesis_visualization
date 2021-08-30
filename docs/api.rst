Classes
****************

*CreateVariants* class
=========================
.. autoclass:: mutagenesis_visualization.main.classes.create_variants.CreateVariants
    :members: __call__

*Counts* class
=========================
.. autoclass:: mutagenesis_visualization.main.classes.counts.Counts

*LibraryRepresentation* pyplot class
---------------------------------------
.. autoclass:: mutagenesis_visualization.main.bar_graphs.library_representation.LibraryRepresentation
    :members: __call__

*MeanCounts* pyplot class
-----------------------------
.. autoclass:: mutagenesis_visualization.main.bar_graphs.mean_counts.MeanCounts
    :members: __call__

*GeneratePrimers* class
=========================
.. autoclass:: mutagenesis_visualization.main.classes.generate_primers.GeneratePrimers
    :members: __call__

*Screen* class
=========================
.. autoclass:: mutagenesis_visualization.main.classes.screen.Screen


The following classes are integrated into *Screen*, thus, you only have to use the __call__ method.


*EnrichmentBar* pyplot class
-----------------------------
.. autoclass:: mutagenesis_visualization.main.bar_graphs.enrichment_bar.EnrichmentBar
    :members: __call__

*Differential* pyplot class
-----------------------------
.. autoclass:: mutagenesis_visualization.main.bar_graphs.differential.Differential
    :members: __call__

*PositionBar* pyplot class
-----------------------------
.. autoclass:: mutagenesis_visualization.main.bar_graphs.position_bar.PositionBar
    :members: __call__

*Secondary* pyplot class
-----------------------------
.. autoclass:: mutagenesis_visualization.main.bar_graphs.secondary.Secondary
    :members: __call__

*Kernel* pyplot class
-----------------------------
.. autoclass:: mutagenesis_visualization.main.kernel.kernel.Kernel
    :members: __call__

*Histogram* pyplot class
-----------------------------
.. autoclass:: mutagenesis_visualization.main.kernel.histogram.Histogram
    :members: __call__

*MultipleKernel* pyplot class
-----------------------------
.. autoclass:: mutagenesis_visualization.main.kernel.multiple_kernels.MultipleKernel

*Heatmap* pyplot class
-----------------------------
.. autoclass:: mutagenesis_visualization.main.heatmaps.heatmap.Heatmap
    :members: __call__

*HeatmapColumns* pyplot class
-----------------------------
.. autoclass:: mutagenesis_visualization.main.heatmaps.heatmap_columns.HeatmapColumns
    :members: __call__

*HeatmapRows* pyplot class
-----------------------------
.. autoclass:: mutagenesis_visualization.main.heatmaps.heatmap_rows.HeatmapRows

*Miniheatmap* pyplot class
-----------------------------
.. autoclass:: mutagenesis_visualization.main.heatmaps.miniheatmap.Miniheatmap
    :members: __call__

*Rank* pyplot class
-----------------------------
.. autoclass:: mutagenesis_visualization.main.other_stats.rank.Rank
    :members: __call__

*Cumulative* pyplot class
-----------------------------
.. autoclass:: mutagenesis_visualization.main.other_stats.cumulative.Cumulative
    :members: __call__

*ROC* pyplot class
-----------------------------
.. autoclass:: mutagenesis_visualization.main.other_stats.roc_analysis.ROC

*Correlation* pyplot class
-----------------------------
.. autoclass:: mutagenesis_visualization.main.pca_analysis.correlation.Correlation
    :members: __call__

*IndividualCorrelation* pyplot class
------------------------------------------
.. autoclass:: mutagenesis_visualization.main.pca_analysis.individual_correlation.IndividualCorrelation
    :members: __call__

*PCA* pyplot class
-----------------------------
.. autoclass:: mutagenesis_visualization.main.pca_analysis.pca.PCA
    :members: __call__

*DifferentialP* plotly class
-----------------------------
.. autoclass:: mutagenesis_visualization.main.plotly.differential.DifferentialP
    :members: __call__

*EnrichmentBarP* plotly class
-----------------------------
.. autoclass:: mutagenesis_visualization.main.plotly.enrichment_bar.EnrichmentBarP

*HeatmapP* plotly class
-----------------------------
.. autoclass:: mutagenesis_visualization.main.plotly.heatmap.HeatmapP
    :members: __call__

*HistogramP* plotly class
-----------------------------
.. autoclass:: mutagenesis_visualization.main.plotly.histogram.HistogramP
    :members: __call__

*RankP* plotly class
-----------------------------
.. autoclass:: mutagenesis_visualization.main.plotly.rank.RankP
    :members: __call__

*Scatter3DPDB* plotly class
-----------------------------
.. autoclass:: mutagenesis_visualization.main.plotly.scatter_3d_pdb.Scatter3DPDB
    :members: __call__

*Scatter3D* plotly class
-----------------------------
.. autoclass:: mutagenesis_visualization.main.plotly.scatter_3d.Scatter3D
    :members: __call__

*ScatterP* plotly class
-----------------------------
.. autoclass:: mutagenesis_visualization.main.plotly.scatter.ScatterP
    :members: __call__

*Scatter* pyplot class
-----------------------------
.. autoclass:: mutagenesis_visualization.main.scatter.scatter.Scatter
    :members: __call__

*ScatterReplicates* pyplot class
----------------------------------
.. autoclass:: mutagenesis_visualization.main.scatter.scatter_replicates.ScatterReplicates
    :members: __call__

Functions
****************


.. autofunction:: mutagenesis_visualization.calculate_enrichment

.. autofunction:: mutagenesis_visualization.count_reads

.. autofunction:: mutagenesis_visualization.count_fastq

.. autofunction:: mutagenesis_visualization.load_demo_datasets

.. autofunction:: mutagenesis_visualization.run_demo

The following function *generate_default_kwargs* is not called by the user as a function. It contains the kwargs that
are parameters of the Screen methods.

.. autofunction:: mutagenesis_visualization.main.utils.kwargs.generate_default_kwargs
