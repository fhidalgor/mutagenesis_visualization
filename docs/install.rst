.. role:: py(code)
      :language: python

.. role:: bash(code)
      :language: bash
	  
	  
Installation
***************

Installing with pip
====================

``Mutagenesis_visualization`` is compatible with Python 3.6. The code is available on `GitHub <https://github.com/fhidalgor/mutagenesis_visualization>`_ under a GNU GENERAL PUBLIC LICENSE. The package can be installed from `PyPI <XXXX>`_ using the ``pip`` package manager by executing the following at the command line:

.. code-block:: bash

     pip install mutagenesis_visualization

If you are working on a jupyter notebook, use the following command:

.. code:: ipython3

    import sys
    !{sys.executable} -m pip install mutagenesis_visualization

Alternative installation
=========================

You may decide to download the jupyter notebook called ``mutagenesis_visualization.ipynb`` also found on `GitHub <https://github.com/fhidalgor/mutagenesis_visualization>`, which contains the code used in this package, and do some modifications of your own. If you do so, there is an easy way to use that same notebook without having to convert it everyt time to a .py file. The way to do that is to download the script ``Import_notebook.py`` and place it in the same folder where you have the ``mutagenesis_visualization.ipynb`` notebook and your current notebook. Should you go through this route, you will need to manually install the required :ref:`dependencies`.

.. code:: ipython3

    import Import_notebook
    import mutagenesis_visualization as mut	

Quick Start
=============
Now that you have installed ``mutagenesis_visualization``, execute the following within Python to evaluate whether it is working propertly:

.. code:: ipython3

	import mutagenesis_visualization as mut
	mut.demo()

This command will load the ``mutagenesis_visualization`` package, create a ``Screen.object`` with sample data, call the ``object.heatmap`` method and show a heatmap plot of the sample data.

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

- `Shannon`_ (version 1.0.0)

You may have already installed ``Pymol``. However, if it is not on the same path as Python, there will not be communication between the two. An easy way to circumvent the problem is to reinstall ``Pymol`` using the following command:

.. code:: ipython3

	conda install -c schrodinger pymol-bundle