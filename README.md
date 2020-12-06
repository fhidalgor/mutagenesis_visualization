# Mutagenesis Visualization
[![Maintenance](https://img.shields.io/badge/maintained%3F-yes-brightgreen.svg)](https://github.com/fhidalgor/mutagenesis_visualization/graphs/commit-activity)
[![Documentation Status](https://readthedocs.org/projects/mutagenesis-visualization/badge/?version=latest)](https://mutagenesis-visualization.readthedocs.io/en/latest/?badge=latest)
[![Build Status](https://travis-ci.org/fhidalgor/mutagenesis_visualization.svg?branch=master)](https://travis-ci.org/fhidalgor/mutagenesis_visualization)
[![codecov](https://codecov.io/gh/fhidalgor/mutagenesis_visualization/branch/master/graph/badge.svg?token=QEAHI2DQDE)](https://codecov.io/gh/fhidalgor/mutagenesis_visualization)
[![PyPI version](https://badge.fury.io/py/mutagenesis-visualization.svg)](https://badge.fury.io/py/mutagenesis-visualization)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-brightgreen.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/fhidalgor/mutagenesis_visualization/HEAD?filepath=mutagenesis_visualization%2Ftutorial%2F)

## Overview
Mutagenesis_visualization is a Python package aimed to generate publication-quality figures for site-saturation mutagenesis datasets.

The package main focus is to perform the statistical analysis and visualization steps of your pipeline, but it additionally offers tools to calculate enrichment scores from FASTQ files.

## Key Features

- Calculate enrichment scores from FASTQ files, allowing for different ways of data processing and normalization.
- Produce publication-quality heatmaps from enrichment scores as well as a wide range of visualization plots.
- Principal component analysis (PCA), hierarchical clustering and receiver operating characteristic (ROC) curve tools.
- Map enrichment scores effortlessly onto a PDB structure using Pymol. Structural properties such as SASA, B-factor or atom coordinates can be extracted from the PDB and visualized using a built-in method.

## Workflow

![Workflow](/docs/_static/workflow_v3.png)

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

You can find the documentation [here](https://mutagenesis-visualization.readthedocs.io/en/latest/).

## Tutorial

There are 7 jupyter notebooks in the folder mutagenesis_visualization/tutorial that go through the basics on how to use the software. You can play with them online via [mybinder](https://mybinder.org/v2/gh/fhidalgor/mutagenesis_visualization/HEAD?filepath=mutagenesis_visualization%2Ftutorial%2F) without having to download anything.

