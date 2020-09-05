=====================================
Welcome to Mutagenesis Visualization!
=====================================

Overview
==========

Description
-------------

``mutagenesis_visualization`` is a Python package aimed to generate publication-quality figures for saturation mutagenesis datasets. The package will take as an input fastq files, count the reads and calculate the enrichment scores for each variant in the population. Then it will create a stunning heatmap as well as perform further statistical analysis on the data.

Key Features
-------------

	1. Calculate enrichment scores from fastq files, allowing for different ways of data processing and normalization.
	2. Produce publication-quality heatmaps from enrichment scores.
	3. Numerous visualization plots to slice the data.
	4. Principal component analysis (PCA) and receiver operating characteristic (ROC) curve tools.
	5. Calculates correlation between evolutionary sequence alignment and enrichment scores.
	6. Map enrichment scores on a Pymol session.
	

Attribution and Contact
-------------------------

This package was written by Frank Hidalgo, as part of his PhD work. Should you have any questons or concerns, please feel free to get in touch at fhidalgoruiz<at>berkeley<dot>edu

Getting Started
===============

In this section, you will find how to install the package and a manual with the different classes and functions in this package.

.. toctree::
   :maxdepth: 1
   :hidden:
   :Caption: Getting Started

	Installation <install>
	Implementation <implementation>

Examples
===============

In this section, we will walk the user through the different functions and methods that can be implemented. We will start from a fastq file, we will process the data and we will do each type of plot.

.. toctree::
   :maxdepth: 1
   :hidden:
   :Caption: Examples
	
	Processing <mut_processing>
	Plotting <mut_plotting>





