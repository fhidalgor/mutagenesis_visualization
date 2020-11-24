# Mutagenesis Visualization
[![Documentation Status](https://readthedocs.org/projects/mutagenesis-visualization/badge/?version=latest)](https://mutagenesis-visualization.readthedocs.io/en/latest/?badge=latest)
[![Build Status](https://travis-ci.org/fhidalgor/mutagenesis_visualization.svg?branch=master)](https://travis-ci.org/fhidalgor/mutagenesis_visualization)
[![PyPI version](https://badge.fury.io/py/mutagenesis-visualization.svg)](https://badge.fury.io/py/mutagenesis-visualization)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-brightgreen.svg)](https://www.gnu.org/licenses/gpl-3.0)

Mutagenesis_visualization is a Python package aimed to generate publication-quality figures for saturation mutagenesis datasets.

The package main focus is to perform the statistical analysis and visualization steps of your pipeline, but it additionally offers tools to calculate enrichment scores from fastq files.

## Key Features

- Calculate enrichment scores from FASTQ files, allowing for different ways of data processing and normalization.
- Produce publication-quality heatmaps from enrichment scores as well as a wide range of visualization plots.
- Principal component analysis (PCA), hierarchical clustering and receiver operating characteristic (ROC) curve tools.
- Map enrichment scores effortlessly onto a PDB structure using Pymol. Structural properties such as SASA, B-factor or atom coordinates can be extracted from the PDB and visualized using a built-in method.

## Installation

Mutagenesis Visualization can be installed from PyPI by executing:

```
pip install mutagenesis_visualization
```

If you prefer to install from Github, use:

```
pip install git+https://github.com/fhidalgor/mutagenesis_visualization
```

## Documentation

You can find the documentation at `readthedocs <https://mutagenesis-visualization.readthedocs.io/en/latest/>`_


