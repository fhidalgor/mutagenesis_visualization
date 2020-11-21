=====================================
Welcome to Mutagenesis Visualization!
=====================================

Overview
==========

.. toctree::
   :hidden:
   :Caption: Overview

	Overview <overview>
	

Description
-------------

``mutagenesis_visualization`` is a Python package aimed to generate publication-quality figures for saturation mutagenesis datasets.

.. image:: images/exported_images/heatmap_intro_v2.png
   :width: 350px
   :align: center

The package main focus is to perform the statistical analysis and visualization steps of your pipeline, but it additionally offers tools to calculate enrichment scores from fastq files.

Unlike other available python packages, we have developed a user-centered API which does not require prior experience with Python nor statistics. The documentation provides multiple examples of how to perform each step. As the user, you will be guided to input your dataset and the protein sequence. From here, the software *prend le contrôle*, and will produce a wide range of stunning and detailed plots. 

.. image:: images/exported_images/user_experience_v2.png
   :width: 250px
   :align: center

Key Features
-------------

	- Calculate enrichment scores from fastq files, allowing for different ways of data processing and normalization.
	- Produce publication-quality heatmaps from enrichment scores as well as a wide range of visualization plots.
	- Principal component analysis (PCA), hierarchical clustering and receiver operating characteristic (ROC) curve tools.
	- Map enrichment scores effortlessly onto a PDB structure using Pymol. Structural properties such as SASA, B-factor or atom coordinates can be extracted from the PDB and visualized using a built-in method.


Getting Started
===============

In this chapter, you will find how to install the package (:ref:`installation guide`) and how to rapidly test that the software is up and running (:ref:`quick demo`). You will also find a workflow.

.. toctree::
   :maxdepth: 1
   :hidden:
   :Caption: Getting Started

	Installation guide <install>

API Description
===============

In here, you will find the :ref:`class`, methods and :ref:`functions` used in this API.

.. toctree::
   :maxdepth: 2
   :hidden:
   :Caption: API Description

	API <api>
	
Tutorial
===============

In this chapter, we will walk the user through the different functions and methods of this Python library. We will start with :ref:`Design DNA libraries` by seeing how to generate the primers to synthesize the DNA library, or the input FASTA file containing all possible site-saturation sequences that companies like Twist Bioscience need in order to synthesize the library for you. Then, from a FASTQ file, we will process the data (:ref:`Processing DNA reads`) and we will do each type of plot (:ref:`Creating plots`). :ref:`Normalizing datasets` shows the different options of data normalization that the package allows for. :ref:`Holistic testing` uses other datasets to showcase the different options that the software gives you. The jupyter notebooks used to generate the examples can be found on `GitHub <https://github.com/fhidalgor/mutagenesis_visualization>`_ and are named ``doc_library.ipynb``, ``doc_processing.ipynb``, ``doc_normalizing.ipynb``, ``doc_plotting.ipynb``, ``doc_plotly.ipynb`` and ``doc_other_datasets.ipynb``. 

.. toctree::
   :maxdepth: 2
   :hidden:
   :Caption: Tutorial
	
	Design DNA libraries <doc_libraries>	
	Processing DNA reads <doc_processing>
	Normalizing datasets <doc_normalizing>
	Creating plots <doc_plotting>
	Using Plotly <doc_plotly>
	Holistic testing <doc_other_datasets>

About Us
===============

Get to know more about the :ref:`Frank Hidalgo`, :ref:`Sage Templeton`, :ref:`Joanne Wang`, and :ref:`Che Olavarria Gallegos`.

.. toctree::
   :maxdepth: 2
   :hidden:
   :Caption: About Us
   
	Bio <aboutme>
