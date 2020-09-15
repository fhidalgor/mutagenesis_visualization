=====================================
Welcome to Mutagenesis Visualization!
=====================================

Overview
==========

Description
-------------

``mutagenesis_visualization`` is a Python package aimed to generate publication-quality figures for saturation mutagenesis datasets. While the software offers tools to calculate enrichment scores from fastq files, the main focus of this package is to perform the statistical analysis and and the visualization steps of your pipeline. Unlike other available python packages, ``mutagenesis_visualization`` has a user-centered design. The customer will convert their input dataset into an object, and from this object each plot can be called by typing the command ``object.method``. A vast number of examples is provided in the documentation so the user can see hot to tune the parameters to suit their own requirements.

Key Features
-------------

	1. Calculate enrichment scores from fastq files, allowing for different ways of data processing and normalization.
	2. Produce publication-quality heatmaps from enrichment scores as well as a wide range of visualization plots.
	3. Principal component analysis (PCA), hierarchical clusterint and receiver operating characteristic (ROC) curve tools.
	4. Map enrichment scores effortlessly onto a PDB structure using Pymol. Structural properties such as SASA, B-factor or atom coordinates can be extracted from the PDB and visualized using a built-in method.

	
Attribution and Contact
-------------------------

This package was written by Frank Hidalgo, as part of his PhD work. Should you have any questons or concerns, please feel free to get in touch at fhidalgoruiz@berkeley.edu


Getting Started
===============

In this chapter, you will find how to install the package (:ref:`installation`) and the :ref:`API` containing the class, methods and functions in this package.

.. toctree::
   :maxdepth: 1
   :hidden:
   :Caption: Getting Started

	Installation <install>
	API <api>

Examples
===============

In this chapter, we will walk the user through the different functions and methods ``mutagenesis_visualization`` offers. We will start from a fastq file, we will process the data (:ref:`processing`) and we will do each type of plot (:ref:`plotting`). :ref:`Normalizing` shows the different options of data normalization that the package allows for. :ref:`more examples` uses other datasets to showcase the different options that the software gives you. The jupyter notebooks used to generate the examples can be found on `GitHub <https://github.com/fhidalgor/mutagenesis_visualization>`_ and are named ``mut_processing.ipynb``, ``mut_normalizing.ipynb``, ``mut_plotting.ipynb`` and ``mut_moreexamples.ipynb``. 

.. toctree::
   :maxdepth: 1
   :hidden:
   :Caption: Examples
	
	Processing <mut_processing>
	Normalizing <mut_normalizing>
	Plotting <mut_plotting>
	More examples <mut_moreexamples>




