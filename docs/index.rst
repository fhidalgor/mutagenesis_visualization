=====================================
Welcome to Mutagenesis Visualization!
=====================================

Overview
==========

Description
-------------

``mutagenesis_visualization`` is a Python package aimed to generate publication-quality figures for saturation mutagenesis datasets. The package will take as an input fastq files and count the DNA reads. The log10 ratio of the counts of each variant from the pre-selected and post-selection samples will yield the enrichment scores. Once the enrichment scores are calculated, the package will create a series of plots such as a stunning heatmap, as well as perform a statistical analysis on the data.

Key Features
-------------

	1. Calculate enrichment scores from fastq files, allowing for different ways of data processing and normalization.
	2. Produce publication-quality heatmaps from enrichment scores.
	3. Numerous visualization plots to analyze the data.
	4. Principal component analysis (PCA) and receiver operating characteristic (ROC) curve tools.
	5. Calculates correlation between evolutionary sequence alignment and enrichment scores.
	6. Map enrichment scores on a Pymol session. Extract properties from PDB files (ie. SASA) and compare it to the enrichment scores.
	

Attribution and Contact
-------------------------

This package was written by Frank Hidalgo, as part of his PhD work. Should you have any questons or concerns, please feel free to get in touch at fhidalgoruiz@berkeley.edu


Getting Started
===============

In this chapter, you will find how to install the package (:ref:`installation`) and a manual (:ref:`implementation`) with the class, methods and functions in this package.

.. toctree::
   :maxdepth: 1
   :hidden:
   :Caption: Getting Started

	Installation <install>
	Implementation <implementation>

Examples
===============

In this chapter, we will walk the user through the different functions and methods ``mutagenesis_visualization`` offers. We will start from a fastq file, we will process the data (:ref:`processing`) and we will do each type of plot (:ref:`plotting`). :ref:`Normalizing` shows the different options of data normalization that the package allows for. :ref:`more examples`) shows uses other datasets to showcase the different options that the software gives you. The jupyter notebooks used to generate the examples can be found on `GitHub <https://github.com/fhidalgor/mutagenesis_visualization>`_ and are named ``mut_processing.ipynb``, ``mut_normalizing.ipynb``, ``mut_plotting.ipynb`` and ``mut_moreexamples.ipynb``. 

.. toctree::
   :maxdepth: 1
   :hidden:
   :Caption: Examples
	
	Processing <mut_processing>
	Normalizing <mut_normalizing>
	Plotting <mut_plotting>
	More examples <mut_moreexamples>




