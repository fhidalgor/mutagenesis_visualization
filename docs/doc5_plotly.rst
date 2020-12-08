Visualizing with plotly
=======================

In this section we will use Plotly to make interactive plots. Please let
us know if you have suggestions for new figures that could be made with
Plotly.

.. code:: ipython3

    import numpy as np
    import pandas as pd
    
    try:
        import mutagenesis_visualization as mut
    except ModuleNotFoundError:  # This step is only for when I run the notebooks locally
        import sys
        sys.path.append('../../')
        import mutagenesis_visualization as mut

Create objects.

.. code:: ipython3

    # Load enrichment scores. This is how you would load them from a local file.
    hras_enrichment_GAPGEF = np.genfromtxt(
        '../data/HRas166_GAPGEF.csv', delimiter=','
    )
    
    hras_enrichment_RBD = np.genfromtxt('../data/HRas166_RBD.csv', delimiter=',')
    
    # Define protein sequence
    hras_sequence = 'MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEY'\
                    + 'SAMRDQYMRTGEGFLCVFAINNTKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDLAARTVES'\
                    + 'RQAQDLARSYGIPYIETSAKTRQGVEDAFYTLVREIRQHKLRKLNPPDESGPG'
    
    # Order of amino acid substitutions in the hras_enrichment dataset
    aminoacids = list('ACDEFGHIKLMNPQRSTVWY*')
    
    # First residue of the hras_enrichment dataset. Because 1-Met was not mutated, the dataset starts at residue 2
    start_position = 2
    
    # Define secondary structure
    secondary = [['L0'], ['β1'] * (9 - 1), ['L1'] * (15 - 9), ['α1'] * (25 - 15),
                 ['L2'] * (36 - 25), ['β2'] * (46 - 36), ['L3'] * (48 - 46),
                 ['β3'] * (58 - 48), ['L4'] * (64 - 58), ['α2'] * (74 - 64),
                 ['L5'] * (76 - 74), ['β4'] * (83 - 76), ['L6'] * (86 - 83),
                 ['α3'] * (103 - 86), ['L7'] * (110 - 103), ['β5'] * (116 - 110),
                 ['L8'] * (126 - 116), ['α4'] * (137 - 126), ['L9'] * (140 - 137),
                 ['β6'] * (143 - 140), ['L10'] * (151 - 143), ['α5'] * (172 - 151),
                 ['L11'] * (190 - 172)]
    
    # Substitute Nan values with 0
    fillna = 0
    
    # Create objects
    hras_GAPGEF = mut.Screen(
        hras_enrichment_GAPGEF, hras_sequence, aminoacids, start_position, fillna,
        secondary
    )
    hras_RBD = mut.Screen(
        hras_enrichment_RBD, hras_sequence, aminoacids, start_position, fillna,
        secondary
    )

Heatmap
-------

Plot an interactive heatmap. Hopover individual pixels to get their
characteristics.

.. code:: ipython3

    hras_RBD.heatmap_plotly(
        title='H-Ras Heatmap',
        output_html='../../docs/html/hras_heatmap.html',
        figsize=(6, 2.5),
    )

.. raw:: html
    :file: html/hras_heatmap.html

Mean
----

Analogous function to ``object.mean`` but rendered using plotly. Will
plot the mean enrichment score for every position on a bar chart. It
will be colored blue for loss of function and red for gain of function.
Additionally, setting the parameter ``mode`` to an amino acid (using the
one letter code) will plot the enrichment for that particular amino acid
along the protein. In this example, we are showing the mean enrichment
scores (top) and an alanine scan (bottom)

.. code:: ipython3

    hras_RBD.mean_plotly(
        title='Mean',
        output_html='../../docs/html/hras_mean.html',
        figsize=(6, 2.5),
    )
    
    hras_RBD.mean_plotly(
        title='A scan',
        mode='A',
        output_html='../../docs/html/hras_mean_A.html',
        figsize=(6, 2.5),
    )

.. raw:: html
    :file: html/hras_mean.html
.. raw:: html
    :file: html/hras_mean_A.html

Histogram
---------

Plot a histogram.

.. code:: ipython3

    hras_RBD.histogram_plotly(
        title='Histogram',
        output_html='../../docs/html/hras_histogram.html',
        figsize=(3, 2.5),
    )

.. raw:: html
    :file: html/hras_histogram.html

Rank
----

Methods reviewed in this section:
    - :meth:`mutagenesis_visualization.Screen.rank_plotly`


Create an interactive rank figure that displays each mutant. You can
export to an html file by giving a path to the variable ``output_html``.

.. code:: ipython3

    hras_RBD.rank_plotly(
        title='Rank of pointmutants',
        output_html='../../docs/html/hras_rankpointmutants.html'
    )

.. raw:: html
    :file: html/hras_rankpointmutants.html

Now display the rank of the positional mean.

.. code:: ipython3

    hras_RBD.rank_plotly(
        mode='mean',
        title='Rank of positions',
        output_html='../../docs/html/hras_rankposition.html',
    )

.. raw:: html
    :file: html/hras_rankposition.html

The following property is applicable to any of the plotly figures. If
you set ``return_plotly_object=True``, you will be able to get the
plotly object and edit it.

.. code:: ipython3

    # Obtain the object
    rank_plotly = hras_RBD.rank_plotly(
        mode='mean',
        title='Rank of positions',
        return_plotly_object=True,
    )
    
    # Then edit rank_plotly

Scatter
-------

Methods reviewed in this section:
    - :meth:`mutagenesis_visualization.Screen.scatter_plotly`


If you have two datasets, you can create a scatter plot. The advantage
of using plotly over matplotlib is that you can visually check each data
point by putting the mouse pointer on top.

.. code:: ipython3

    hras_RBD.scatter_plotly(
        hras_GAPGEF,
        show_results=False,
        title='Scatter Point Mutants',
        x_label='hras_RBD',
        y_label='hras_GAPGEF',
        output_html='../../docs/html/hras_scatterpointmutants.html',
    )

.. raw:: html
    :file: html/hras_scatterpointmutants.html

Now we just look at the positional average.

.. code:: ipython3

    hras_RBD.scatter_plotly(
        hras_GAPGEF,
        mode='mean',
        title='Scatter Positional Average',
        x_label='hras_RBD',
        y_label='hras_GAPGEF',
        output_html='../../docs/html/hras_scatterposition.html',
    )

.. raw:: html
    :file: html/hras_scatterposition.html

3D scatter plot
---------------

Methods reviewed in this section:
    - :meth:`mutagenesis_visualization.Screen.scatter_3D_plotly`


If there is an available PDB structure, you can input it and the
software will plot a 3d plot of the C-alpha atoms, colored by their
enrichment score.

The method ``object.scatter_3D_plotly`` will take as an input either a
PDB file (``pdb_path=/path/to/file``) or the x,y,z coordinates
(``df_coordinates``).

.. code:: ipython3

    hras_RBD.scatter_3D_plotly(
        mode='mean',
        pdb_path='../data/5p21.pdb',
        title='Scatter 3D',
        squared=False,
        x_label='x',
        y_label='y',
        z_label='z',
        output_html='../../docs/html/hras_3dscatter.html',
    )

.. raw:: html
    :file: html/hras_3dscatter.html

By setting up mode=‘V’, we can evaluate the impact of valine
substitutions. Mode can be set up to any residue. In this example,
residues in the core are tolerant to valine substitutions.

.. code:: ipython3

    hras_RBD.scatter_3D_plotly(
        mode='V',
        pdb_path='../data/5p21.pdb',
        title='Scatter 3D - Valine substitution',
        squared=False,
        x_label='x',
        y_label='y',
        z_label='z',
        output_html='../../docs/html/hras_3dvalsubstitution.html',
    )

.. raw:: html
    :file: html/hras_3dvalsubstitution.html

When we set mode=‘D’, the core of the protein turns completely blue.

.. code:: ipython3

    hras_RBD.scatter_3D_plotly(
        mode='D',
        pdb_path='../data/5p21.pdb',
        title='Scatter 3D - Aspartate substitution',
        squared=False,
        x_label='x',
        y_label='y',
        z_label='z',
        output_html='../../docs/html/hras_3daspsubstitution.html',
    )

.. raw:: html
    :file: html/hras_3daspsubstitution.html

By setting squared = True, we plot the distance to the center of the
protein of each residue. In this example, we see that residues in the
core of the protein are blue, indicating a sensitivity to mutations.

.. code:: ipython3

    hras_RBD.scatter_3D_plotly(
        mode='mean',
        pdb_path='../data/5p21.pdb',
        title='Scatter 3D - Distance to center',
        squared=True,
        x_label='x',
        y_label='y',
        z_label='z',
        output_html='../../docs/html/hras_3ddistcenter.html',
    )

.. raw:: html
    :file: html/hras_3ddistcenter.html

PDB properties
--------------

From the PDB, properties such as B-factor or SASA can be extracted.
Using plotly we allow the user to have a 3-D scatter plot colored by the
enrichment scores. You can additionally include other properties to
include such as the conservation scores using the parameter ``custom``.

.. code:: ipython3

    # Calculate conservation score from MSA
    path = '../data/Ras_family_trimmed.fasta'
    
    # Calculate shannon scores
    df_shannon, df_freq = mut.msa_enrichment(
        hras_RBD, path, start_position=1, threshold=0.1
    )
    
    # Plot 3-D SASA, log B-factor and Shannon Entropy
    hras_RBD.scatter_3D_pdbprop_plotly(
        plot=['SASA', 'log B-factor', 'Shannon'],
        custom=df_shannon['Shannon'],
        pdb_path='../data/5p21.pdb',
        title='Scatter 3D - PDB properties',
        output_html='../../docs/html/hras_3d_pdbprop.html',
    )

.. raw:: html
    :file: html/hras_3d_pdbprop.html
