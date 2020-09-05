#!/usr/bin/env python
# coding: utf-8

# # Plotting

# This notebook shows how to use the mutagenesis_visualization package. The plotting functions can be used regardless of how you process your data. For the examples, we are using two datasets that are derived from Pradeep's [#Pradeep2017]_legacy.

# ## Import module

# In[ ]:


# Import Module
import Import_notebook
import mutagenesis_visualization as mut
import numpy as np
import pandas as pd


# ## Create object of class Screen

# In order to create plots, the data needs to be passed as an object of the class ``Screen``. The first step is to load the enrichment score data. Once the data is loaded, other parameters such as the protein sequence or the amino acid substitutions need to be defined for the object to be created. Adding the secondary structure is optional, but without it some plots won't work. In this example, we are importing two datasets and creating two objects named ``hras_GAPGEF`` and ``hras_RBD``.

# In[ ]:


# Load enrichment scores
hras_enrichment_GAPGEF = np.genfromtxt('Exported/HRas166_GAPGEF.csv', delimiter=',')
hras_enrichment_RBD = np.genfromtxt('Exported/HRas166_RBD.csv', delimiter=',')

# Define protein sequence
hras_sequence = 'MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDLAARTVESRQAQDLARSYGIPYIETSAKTRQGVEDAFYTLVREIRQHKLRKLNPPDESGPG'

# Order of amino acid substitutions in the hras_enrichment dataset
aminoacids = list('ACDEFGHIKLMNPQRSTVWY*')

# First residue of the hras_enrichment dataset. Because 1-Met was not mutated, the dataset starts at residue 2
start_position = 2

# Define secondary structure
secondary = [['L0'], ['β1']*(9-1), ['L1']*(15-9), ['α1']*(25-15), ['L2']*(36-25), ['β2']*(46-36), ['L3']*(48-46), ['β3']*(58-48), ['L4'] * (64-58), ['α2']*(74-64), ['L5']*(76-74), ['β4']*(83-76), [
    'L6']*(86-83), ['α3']*(103-86), ['L7']*(110-103), ['β5']*(116-110), ['L8']*(126-116), ['α4']*(137-126), ['L9']*(140-137), ['β6']*(143-140), ['L10']*(151-143), ['α5']*(172-151), ['L11']*(190-172)]

# Substitute Nan values with 0
fillna = 0

# Create objects
hras_GAPGEF = mut.Screen(hras_enrichment_GAPGEF, hras_sequence,
                         aminoacids, start_position, fillna, secondary)
hras_RBD = mut.Screen(hras_enrichment_RBD, hras_sequence,
                         aminoacids, start_position, fillna, secondary)

# Parameters to save output images, will be the same for each plot
outputfilepath = 'mutagenesis_visualization_package/example/exported_images/'
outputformat = 'png'
savefile = True


# ## Heatmaps

# Once the object ``hras_RBD`` is created, we will plot a heatmap of the enrichment scores using the attribute ``object.heatmap``.

# In[ ]:


# Create full heatmap
hras_RBD.heatmap(title='H-Ras 2-166', outputfilename='hras_fullheatmap', outputfilepath=outputfilepath,
                    show_cartoon=True, outputformat=outputformat, savefile=savefile)

.. image:: ../example/exported_images/hras_fullheatmap.png
# If you set the parameter ``show_snv=True``, the algorithm will color green every mutation that is not a single nucleotide variant of the wild-type protein. You will notice how many mutations are not accessible through a nucleotide change.

# In[ ]:


# Create full heatmap
hras_RBD.heatmap(title='H-Ras 2-166', outputfilename='hras_fullheatmap_snv', outputfilepath=outputfilepath,
                    show_cartoon=True, outputformat=outputformat, show_snv=True, savefile=savefile)

.. image:: ../example/exported_images/hras_fullheatmap_snv.png
# From the full heatmap, we can slice into just a few rows (a few amino acid substitution profiles). For that, we will use the attribute ``object.heatmap_selection``. Note that we need to specify which amino acids to show.

# In[ ]:


# Create heatmap of selected aminoacid substitutions
hras_RBD.heatmap_selection(title='H-Ras 2-166', outputfilename='hras_selectionheatmap', outputfilepath=outputfilepath,
                              selection=['E', 'Q', 'A', 'P', 'V', 'Y'], outputformat=outputformat, savefile=savefile)

.. image:: ../example/exported_images/hras_selectionheatmap.png
# Instead of displaying only a few rows, we can display only a few columns. For that, we will use the attribute ``object.heatmap_subset``.

# In[ ]:


# Create a heatmap of a subset region in the protein
hras_RBD.heatmap_subset(segment=[20, 40], outputfilename='hras_subsetheatmap', outputfilepath=outputfilepath,
                           outputformat=outputformat, savefile=savefile)

.. image:: ../example/exported_images/hras_subsetheatmap.png
   :width: 200px
   :align: center
# A summarized heatmap can also be generated. It is useful to evaluate global trends in the data. The command to use is ``object.miniheatmap``.

# In[ ]:


# Condensed heatmap
hras_RBD.miniheatmap(title = 'Wt residue H-Ras', outputfilename='hras_miniheatmap', 
             outputfilepath=outputfilepath, outputformat=outputformat, savefile=savefile)

.. image:: ../example/exported_images/hras_miniheatmap.png
   :width: 250px
   :align: center
# ## Histogram, scatter and more

# There are different tools to analyze the data. The package can plot the kernel density estimation (``object.kernel``). There is the option to fit other functions to the data (see Implementation for more). You could also only plot the histogram (``object.histogram``). For the histograms, we can select to plot only the single nucleotide variants (SNVs) or the non-SNVs. In the example, it actually changes the shape of the population. Non-SNVs are more sensitive to mutations than SNVs because there is a higher proportion of non-conservative amino acid replacements. 

# In[ ]:


# Plot histograms and a PDF
hras_RBD.kernel(histogram=True, title='H-Ras 2-166', xscale=[-2, 1], outputfilename='hras_kde',
                   outputfilepath=outputfilepath, outputformat=outputformat, savefile=savefile)
hras_RBD.histogram(population='SNV', title='H-Ras 2-166 SNV', xscale=[-2, 1], outputfilename='hras_histsnv',
                      outputfilepath=outputfilepath, outputformat=outputformat, savefile=savefile)
hras_RBD.histogram(population='nonSNV', title='H-Ras 2-166 non-SNV', xscale=[-2, 1], outputfilename='hras_histnonsnv',
                      outputfilepath=outputfilepath, outputformat=outputformat, savefile=savefile)

.. image:: ../example/exported_images/hras_kde.png
   :width: 240px

.. image:: ../example/exported_images/hras_histsnv.png
   :width: 200px

.. image:: ../example/exported_images/hras_histnonsnv.png
   :width: 200px

# If you have multiple objects, you can create a scatter plot using ``object.scatter``. Or if you have multiple replicates of the same experiment, that would be a way to compare them. We give the option to do the comparison at a mutation by mutation level, or at a position level.

# In[ ]:


# Plot a scatter plot of each mutation
hras_RBD.scatter(hras_GAPGEF, title='Individual mutations', mode='pointmutant',
                 xscale=(-2.5, 1.5), yscale=(-2.5, 1.5),
                    x_label='H-Ras Unregulated', y_label='H-Ras Regulated', outputfilename='hras_scatter',
                    outputfilepath=outputfilepath, outputformat=outputformat, savefile=savefile)

# Plot a scatter plot of the mean position
hras_RBD.scatter(hras_GAPGEF, title='Positional average', mode='mean', xscale=(-2, 1), yscale=(-2, 1),
                    x_label='H-Ras Unregulated', y_label='H-Ras Regulated', outputfilename='hras_scatter_mean',
                    outputfilepath=outputfilepath, outputformat=outputformat, savefile=savefile)

.. image:: ../example/exported_images/hras_scatter.png
   :width: 200px

.. image:: ../example/exported_images/hras_scatter_mean.png
   :width: 200px

# The attribute ``object.rank`` sorts each mutation (or position) by the enrichment score. If ``outdf=True``, it will return a dataframe with the mutations ranked.

# In[ ]:


# Rank plot
hras_RBD.rank(mode='pointmutant', outdf=True, title = 'Rank of mutations', outputfilename='hras_rank',
                    outputfilepath=outputfilepath, outputformat=outputformat, savefile=savefile)

.. image:: ../example/exported_images/hras_rank.png
   :width: 400px
   
.. image:: ../example/exported_images/hras_ranktable.png
   :width: 200px

# The attribute ``object.cumulative`` draws a cumulative plot that sums the mean enrichment score of every position. This  plot is useful to determine if the sensitivity to mutations is constant throughout the protein or not. In the example, we see that the cumulative function follows the x=y line, suggestion a homogeneous mutational tolerance. 
# 

# In[ ]:


# Cumulative plot
hras_RBD.cumulative(mode = 'all', title = 'Cumulative Score', outputfilename='hras_cumulative',
                    outputfilepath=outputfilepath, outputformat=outputformat, savefile=savefile)

.. image:: ../example/exported_images/hras_cumulative.png
   :width: 300px
   :align: center
# ## Bar and line charts

# The attribute ``object.mean`` will plot the mean enrichment score for every position on a bar chart. It will be colored blue for loss of function and red for gain of function. Additionally, setting the parameter ``mode`` to an amino acid will plot the enrichment for that particular amino acid along the protein. 

# In[ ]:


# Plot a bar graph with the mean enrichment score
hras_RBD.mean(figsize=[6, 2.5], mode='mean', show_cartoon=True, yscale=[-2, 0.5], outputfilename='hras_bar_mean',
                 title='', outputfilepath=outputfilepath, outputformat=outputformat, savefile=savefile)

# Plot a bar graph with the alanine enrichment score
hras_RBD.mean(figsize=[6, 2.5], mode='A', show_cartoon=True, yscale=[-2, 0.5], outputfilename='hras_bar_alanine',
                 title='', outputfilepath=outputfilepath, outputformat=outputformat, savefile=savefile)

.. image:: ../example/exported_images/hras_bar_mean.png
   :width: 500px
   :align: center
   
.. image:: ../example/exported_images/hras_bar_alanine.png
   :width: 500px
   :align: center
# The mean differential effect between the two H-Ras datasets is displayed (``object.differential``). This plot is useful to compare either orthologs/paralogs or the same protein with different effectors, and determine which areas of the protein have a different sensitivity to mutations. 

# In[ ]:


# Plot the difference between H-Ras unregulated and H-Ras regulated datasets
# The subtraction is hras_RBD - hrasGAPGEF
hras_RBD.differential(hras_GAPGEF, figsize=[6, 2.5], show_cartoon=True, yscale=[-1, 1], outputfilename='hras_diffenrichment',
                         title='', outputfilepath=outputfilepath, outputformat=outputformat, savefile=savefile)

.. image:: ../example/exported_images/hras_diffenrichment.png
   :width: 500px
   :align: center
# You can check the individual mutational profile of a residue by using ``object.position``.

# In[ ]:


# Create plot for position 117
hras_RBD.position(position = 117, yscale = (-1.5, 0.8), figsize = (3.5,2), title = 'Position 117', outputfilename='hras_position117',
                  outputfilepath=outputfilepath, outputformat=outputformat, savefile=savefile)

.. image:: ../example/exported_images/hras_position117.png
   :width: 350px
   :align: center
# If you added the secondary structure as an attribute of the object, you can plot the mean enrichment score for each alpha and beta motif in the protein (``object.secondary_mean``).

# In[ ]:


# Graph bar of the mean of each secondary motif
hras_RBD.secondary_mean(yscale=[-1, 0], figsize=[3, 2], title='Mean of secondary motifs',
                           outputfilename='hras_secondary', outputfilepath=outputfilepath,
                           outputformat=outputformat, savefile=savefile)

.. image:: ../example/exported_images/hras_secondary.png
   :width: 300px
   :align: center
# ## Correlation, PCA and ROC AUC

# The correlation of amino acid substitution profiles can be calculated for each amino acid and graphed using ``object.correlation``. In the example we observe that polar amino acids have high correlation between themselves but low correlation with hydrophobic amino acids.

# In[ ]:


# Correlation between amino acids
hras_RBD.correlation(colorbar_scale=[0.5,1],title = 'Correlation', outputfilename='hras_correlation', 
             outputfilepath=outputfilepath, outputformat=outputformat, savefile=savefile)

.. image:: ../example/exported_images/hras_correlation.png
   :width: 250px
   :align: center
# If you were to do a single amino acid scan (ie. alanine scan), how would that predict the rest of the amino acid mutational profiles? That can be determined using ``object.meancorrelation``.

# In[ ]:


# Explained variability by amino acid
hras_RBD.meancorrelation(yscale=[0,0.6], title = 'Explained variability by amino acid', outputfilename='hras_variability', 
             outputfilepath=outputfilepath, outputformat=outputformat, savefile=savefile)

.. image:: ../example/exported_images/hras_variability.png
   :width: 300px
   :align: center
# Grouping amino acids improves the predictive power. ``object.representatives`` lets you manually group amino acids. The algorithm picks one amino acid per group and evaluates the predictive power of the subset. Such operation will be done for every possible combination. In the example, 8 amino acids explain 0.75 of the data. The ``logomaker`` package [#Tareen2019]_will show for each group which is the most represented amino acid in of the subset that has an R2 value greater than the cutoff that you have set using the parameter ``r2``. Such plot will let you see if there is any preference for a particular amino acid within a group. This tool can be useful when ordering new libraries, since you can save some cost by ordering less mutants.

# In[ ]:


# Get list of all combinations and their associated R2 value
df_r2 = hras_RBD.representatives(r2=0.75, groups=['DE', 'HKR', 'QN', 'CST', 'AG', 'ILMV', 'WYF', 'P'],
                                 output=False, title='', outputfilename='hras_logo',
                                 outputfilepath=outputfilepath, outputformat=outputformat, savefile=savefile)

# Only show the top 5
df_r2.sort_values(by = 'R2', ascending=False).head()

.. image:: ../example/exported_images/hras_logo.png
   :align: center

.. image:: ../example/exported_images/hras_meanrepresentatives_rank.png
   :width: 200px
   :align: center

# The package can perform principal component analysis (PCA) using the attribute ``object.pca``. The parameter ``mode`` can be set to ``aminoacid``, in which will cluster amino acids based on their similarity, ``individual`` in which will do the same for each individual residue and ``secondary``, in which will cluster for each motif. By default, the first two dimensions will be plotted (0 and 1 in Python notation), but that can be changed by ``dimensions`` parameter.

# In[ ]:


# PCA by amino acid substitution
hras_RBD.pca(title='', dimensions=[0, 1], figsize=(2, 2), adjustlabels=True, outputfilename='hras_pcaaminoacid', 
             outputfilepath=outputfilepath, outputformat=outputformat, savefile=savefile)

# PCA by secondary structure motif
hras_RBD.pca(title='', mode='secondary', dimensions=[0, 1], figsize=(2, 2), 
             adjustlabels=True, outputfilename='hras_pcasecondary',
                    outputfilepath=outputfilepath, outputformat=outputformat, savefile=savefile)

# PCA by each individual residue. Don't set adjustlabels = True unless really big figsize
hras_RBD.pca(title='', mode='individual', dimensions=[0, 1], figsize=(5, 5), 
             adjustlabels=False, outputfilename='hras_pcaindividual',
                    outputfilepath=outputfilepath, outputformat=outputformat, savefile=savefile)
             

.. image:: ../example/exported_images/hras_pcaaminoacid.png
   :width: 200px

.. image:: ../example/exported_images/hras_pcasecondary.png
   :width: 200px

.. image:: ../example/exported_images/hras_pcaindividual.png
   :width: 250px
# Another type of plot that can be done is a receiver operating characteristic (ROC) curve for classification. You will use the attribute ``object.roc`` and as an input you will pass a dataframe that contains the label for each variant. In this example, we are using it to evaluate whether we can use evolutionary conservation data to predict the mutational tolerance of the protein. The area under the curve (AUC) is above 0.5, implying that there is a small relationship between enrichment/conservation.

# In[ ]:


# Calculate conservation score from MSA
path = 'Other/2020_pfam/Ras_family_trimmed.fasta'
df_shannon, df_freq = mut.msa_enrichment(hras_RBD, path, start_position=1, threshold=0.1)

# Plot ROC curve
hras_RBD.roc(df_freq[['Variant', 'Class']], title='MSA predictive power', outputfilename='hras_roc',
                    outputfilepath=outputfilepath, outputformat=outputformat, savefile=savefile)

.. image:: ../example/exported_images/hras_roc.png
   :width: 250px
   :align: center
# The package also allows to do a boxplot using the function ``plot_box``. Note that the data needs to be binned separately. In this example, we have used it to study if the Shannon entropy is related to the mutational tolerance. Although there is high variability, on average residues with a lower enrichment score are more conserved.

# In[ ]:


# Bin data
binned_shannon = (2*df_shannon['Shannon']).round(0)/2

# Plot box plot.
mut.plot_box(binned_x = binned_shannon, y = df_shannon['Score'], title='Shannon vs Enrichment',
             x_label = 'Shannon Entropy', y_label=r'$∆E^i_x$', outputfilename='hras_shannon',
                    outputfilepath=outputfilepath, outputformat=outputformat, savefile=savefile)

.. image:: ../example/exported_images/hras_shannon.png
   :width: 300px
   :align: center
# ## Pymol

# The data can be graphed on a Pymol session using ``object.pymol``. The parameter ``pdb`` will fetch that pdb. Note that the chain to fetch needs to be specified (see example). Red for gain of function and blue for loss of function. ``mode`` lets you specifiy whether to plot an individual amino acid profile (left - Leucine, right - Aspartate) or the mean.

# In[ ]:


# Start pymol and color residues. Cut offs are set with gof and lof parameters.
hras_RBD.pymol(pdb='5p21_A', mode = 'mean', gof=0.2, lof=-0.5)

# Now check the mutational profile of Leucine (left image)
hras_RBD.pymol(pdb='5p21_A', mode = 'L', gof=0.2, lof=-0.5)

# Now check the mutational profile of Aspartate (right image)
hras_RBD.pymol(pdb='5p21_A', mode = 'D', gof=0.2, lof=-0.5)

.. image:: ../example/exported_images/hras_pymol_combLD.png
   :align: center

# ## Reference
.. [#Pradeep2017] Bandaru et al. (2017). Deconstruction of the Ras switching cycle through saturation mutagenesis. `DOI: 10.7554/eLife.27810  <https://elifesciences.org/articles/27810>`_

.. [#Tareen2019] Tareen A, Kinney JB (2019). Logomaker: beautiful sequence logos in Python. `bioRxiv DOI:10.1101/635029. <https://www.biorxiv.org/content/10.1101/635029v1>`_