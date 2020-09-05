==============================
Mutagenesis Visualization
==============================

Overview
**********

Description
=============

``Mutagenesis_visualization`` is a Python package aimed to generate publication-quality figures for saturation mutagenesis datasets. The package will take as an input fastq files, count the reads and calculate the enrichment scores for each variant in the population. Then it will create a stunning heatmap as well as perform further statistical analysis on the data.

Key Features
=============

	1. Calculate enrichment scores from fastq files, allowing for different ways of data processing and normalization.
	2. Produce publication-quality heatmaps from enrichment scores.
	3. Numerous visualization plots to slice the data.
	4. Principal component analysis (PCA) and receiver operating characteristic (ROC) curve tools.
	5. Calculates correlation between evolutionary sequence alignment and enrichment scores.
	6. Map enrichment scores on a Pymol session.
	

Attribution and Contact
========================

This package was written by Frank Hidalgo, as part of his PhD work. Should you have any questons or concerns, please feel free to get in touch at fhidalgoruiz<at>berkeley<dot>edu


.. toctree::
   :maxdepth: 1
   :Caption: Getting Started
	
	Installation <install>
     Implementation <implementation>
	Processing <mut_processing>
	Plotting <mut_plotting>  

