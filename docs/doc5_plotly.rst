Visualizing with plotly
=======================

In this section we will use Plotly to make interactive plots. Please let
us know if you have suggestions for new figures that could be made with
Plotly.

.. code:: ipython3

    %matplotlib inline
    import numpy as np
    from mutagenesis_visualization.main.demo.demo_objects import DemoObjects
    from mutagenesis_visualization.main.utils.data_paths import PDB_5P21
    
    DEMO_OBJECTS:DemoObjects = DemoObjects()
    hras_rbd = DEMO_OBJECTS.hras_rbd
    hras_gapgef = DEMO_OBJECTS.hras_gapgef

Heatmap
-------

Plot an interactive heatmap. Hopover individual pixels to get their
characteristics.

.. code:: ipython3

    hras_rbd.plotly_heatmap(
        title='H-Ras Heatmap',
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

    hras_rbd.plotly_enrichment_bar(
        title='Mean',
        figsize=(6, 2.5),
    )
    
    hras_rbd.plotly_enrichment_bar(
        title='A scan',
        mode='A',
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

    hras_rbd.plotly_histogram(
        title='Histogram',
        figsize=(3, 2.5),
    )

.. raw:: html
    :file: html/hras_histogram.html

Rank
----

Methods reviewed in this section:
    - :meth:`mutagenesis_visualization.Screen.plotly_rank`


Create an interactive rank figure that displays each mutant. The default
mode is set to pointmutant to provide the ranking on the mutation level.
You can download the plot as a png file by clicking the camera icon
which appears on the far left when our cursor is over the plot. You can
export to an html file by giving a path to the variable ``output_html``.

.. code:: ipython3

    hras_rbd.plotly_rank(
        title='Rank of pointmutants',
    )

.. raw:: html
    :file: html/hras_rankpointmutants.html

Now display the rank of the positional mean.

.. code:: ipython3

    hras_rbd.plotly_rank(mode='mean',title='Rank of positions')

.. raw:: html
    :file: html/hras_rankposition.html

Scatter
-------

Methods reviewed in this section:
    - :meth:`mutagenesis_visualization.Screen.plotly_scatter`


If you have two datasets, you can create a scatter plot. The advantage
of using plotly over matplotlib is that you can visually check each data
point by hovoring your cursor over a point. By setting show_results =
True, the OLS regression results will also be printed as output. The
mode = ‘pointmutant’ is default which shows a comparison as mutation by
mutation.

.. code:: ipython3

    hras_rbd.plotly_scatter(
        hras_gapgef,
        show_results=False,
        title='Scatter Point Mutants',
        x_label='hras_rbd',
        y_label='hras_gapgef',
    )

.. raw:: html
    :file: html/hras_scatterpointmutants.html

Now we just look at the positional average.

.. code:: ipython3

    hras_rbd.plotly_scatter(
        hras_gapgef,
        mode='mean',
        title='Scatter Positional Average',
        x_label='hras_rbd',
        y_label='hras_gapgef',
    )

.. raw:: html
    :file: html/hras_scatterposition.html

3D scatter plot
---------------

Methods reviewed in this section:
    - :meth:`mutagenesis_visualization.Screen.plotly_scatter_3d`


If there is an available PDB structure, you can input it and the
software will plot a 3d plot of the C-alpha atoms, colored by their
enrichment score.

The method ``object.plotly_scatter_3d`` will take as an input either a
PDB file (``pdb_path=/path/to/file``) or the x,y,z coordinates
(``df_coordinates``).

.. code:: ipython3

    hras_rbd.plotly_scatter_3d(
        mode='mean',
        pdb_path=PDB_5P21,
        title='Scatter 3D',
        squared=False,
        x_label='x',
        y_label='y',
        z_label='z',
    )

.. raw:: html
    :file: html/hras_3dscatter.html

By setting up mode=‘V’, we can evaluate the impact of valine
substitutions. Mode can be set up to any residue. In this example,
residues in the core are tolerant to valine substitutions.

.. code:: ipython3

    hras_rbd.plotly_scatter_3d(
        mode='V',
        pdb_path=PDB_5P21,
        title='Scatter 3D - Valine substitution',
        squared=False,
        x_label='x',
        y_label='y',
        z_label='z',
    )

.. raw:: html
    :file: html/hras_3dvalsubstitution.html

When we set mode=‘D’, the core of the protein turns completely blue.

.. code:: ipython3

    hras_rbd.plotly_scatter_3d(
        mode='D',
        pdb_path=PDB_5P21,
        title='Scatter 3D - Aspartate substitution',
        squared=False,
        x_label='x',
        y_label='y',
        z_label='z',
    )

.. raw:: html
    :file: html/hras_3daspsubstitution.html

By setting squared = True, we plot the distance to the center of the
protein of each residue. In this example, we see that residues in the
core of the protein are blue, indicating a sensitivity to mutations.

.. code:: ipython3

    hras_rbd.plotly_scatter_3d(
        mode='mean',
        pdb_path=PDB_5P21,
        title='Scatter 3D - Distance to center',
        squared=True,
        x_label='x',
        y_label='y',
        z_label='z',
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

    # Plot 3-D SASA, log B-factor and Shannon Entropy
    hras_rbd.plotly_scatter_3d_pdbprop(
        plot = ['Distance', 'SASA', 'log B-factor'],
        pdb_path=PDB_5P21,
        title='Scatter 3D - PDB properties',
    )

.. raw:: html
    :file: html/hras_3d_pdbprop.html
