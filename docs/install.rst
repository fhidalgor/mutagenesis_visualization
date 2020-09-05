.. role:: py(code)
      :language: python

.. role:: bash(code)
      :language: bash
	  
	  
Installation
***************

Installing with pip
====================

``Mutagenesis_visualization`` is compatible with Python 3.6. The code is available on `GitHub <https://github.com/fhidalgor/mutagenesis_visualization>`_ under a GNU GENERAL PUBLIC LICENSE. The package be installed from `PyPI <XXXX>`_ using the ``pip`` package manager by executing the following at the commandline:

.. code-block:: bash

     pip install mutagenesis_visualization

If you are working on a jupyter notebook, use the following command:

.. code-block:: python

    import sys
    !{sys.executable} -m pip install mutagenesis_visualization

Alternative installation
=========================

You may decide to download the jupyter notebook called ``mutagenesis_visualization`` which contains all the functions used in this package, and do some modifications of your own. If you do so, there is an easy way to use that same notebook without having to convert it to a py file first. The way to do that is to have the script ``Import_notebook`` on the same folder as your current notebook and the notebook you are trying to import. Note that you will need to manually install the required dependencies.

.. code-block:: python

    import Import_notebook
    import mutagenesis_visualization as mut	

Quick Start
=============
For a quick demonstration of mutagenesis_visualization, execute the following within Python:

.. code-block:: python

	import mutagenesis_visualization as mut
	mut.demo()

This command will load the mutagenesis_visualization package and create a heatmap plot of an H-Ras saturation mutagenesis dataset.

.. image:: ../example/exported_images/hras_fullheatmap.png

    
Dependencies
==============

Required Dependencies
-----------------------
- `numpy <http://numpy.org/>`_ (version 1.18.5 or later)

- `matplotlib <http://matplotlib.org/>`_ (version 3.2.2 or later)

- `seaborn <https://seaborn.pydata.org/>`_ (version 0.10.1 or later)

- `pandas <http://pandas.pydata.org/>`_ (version 1.0.5 or later)

- `scipy <http://www.scipy.org/scipylib/index.html>`_ (version 1.5.0 or later)

- `scikit-learn <http://scikit-learn.org/stable/>`_ (version 0.23.1 or later)

- `copy <https://docs.python.org/2/library/copy.html>`_ 

- `itertools <https://docs.python.org/3/library/itertools.html>`_ (version 8.4.0 or later)

- `biopython <https://pypi.org/project/biopython/>`_ (version 1.77 or later)

- `collections <https://docs.python.org/2/library/collections.html>`_ (version 1.2.1 or later)

Optional dependencies
---------------------
- `ipymol <https://github.com/cxhernandez/ipymol>`_ (version 0.5 or later)

- `logomaker <https://logomaker.readthedocs.io>`_ (version 0.8 or later)

- `adjustText <https://pypi.org/project/adjustText/>`_ (version 0.7.3 or later)

- `Shannon <https://pypi.org/project/shannon/>`_ (version 1.0.0 or later)

A way to ensure Pymol is on the same as Python is to install ``Pymol`` using the following command:

.. code:: ipython3

	conda install -c schrodinger pymol-bundle