Processing DNA reads
====================

This section will teach you how to use the built-in data muting
functions. If you already have your own muting pipeline built, you can
skip this section and go to the plotting examples.

Import module
-------------

.. code:: ipython3

    import numpy as np
    import pandas as pd
    try:
        import mutagenesis_visualization as mut
    except ModuleNotFoundError:  # This step is only for when I run the notebooks locally
        import sys
        sys.path.append('../../')
        import mutagenesis_visualization as mut

Count DNA reads from fastq file
-------------------------------

Site saturation mutagenesis
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Methods and functions reviewed in this section:
    - :meth:`mutagenesis_visualization.Screen.meancounts`
    - :func:`mutagenesis_visualization.count_reads`

After sequencing your DNA library, using other packages you will
assemble the forward and reverse reads and trim the flanking bases. That
will produce a trimmed fastq file that contains the DNA reads. This is
where ``mutagenesis_visualization`` kicks in. The following function
``count_reads`` will read your trimmed fastq file and count the number
of times a DNA sequence is present. You will have to pass as inputs a
``dna_sequence`` and a ``codon_list`` with the codons that were used to
make the point mutant library. If ``savefile=True`` , it will export the
results to txt files. Below there is a prettified example of the output
file.

.. code:: ipython3

    # H-Ras dna sequence
    hras_dnasequence = 'acggaatataagctggtggtggtgggcgccggcggtgtgggcaagagtgcgctgaccat'\
        + 'ccagctgatccagaaccattttgtggacgaatacgaccccactatagaggattcctaccggaagcaggtgg'\
        + 'tcattgatggggagacgtgcctgttggacatcctg'
    
    # Codons used to make the NNS library. I could also have used 'NNS' and the package will use the NNS codons
    codon_list = [
        "GCC", "GCG", "TGC", "GAC", "GAG", "TTC", "GGC", "GGG", "CAC", "ATC", "AAG",
        "CTC", "CTG", "TTG", "ATG", "AAC", "CCC", "CCG", "CAG", "CGC", "CGG", "AGG",
        "TCC", "TCG", "AGC", "ACC", "ACG", "GTC", "GTG", "TGG", "TAC", "TAG"
    ]
    
    # Input and output files
    input_file = '../data/hras.trimmed.fastq'
    output_file = 'hras_counts.xlsx'
    counts_wt = False
    start_position = 2
    
    # Execute count reads
    df_counts_pre, wt_counts_pre = mut.count_reads(
        hras_dnasequence, input_file, codon_list, counts_wt, start_position,
        output_file
    )

.. image:: images/exported_images/hras_tablecounts.png
   :width: 450px
   :align: center

Create object of class ``Counts``.

.. code:: ipython3

    hras_obj = mut.Counts(df_counts_pre)

Once the reads have been counted, the method ``mean_counts`` can be used
to evaluate the coverage by position. The method
``library_representation`` will tell you the percentage coverage of each
amino acid per position.

.. code:: ipython3

    hras_obj.mean_counts(title='H-Ras mean counts per position')
    
    hras_obj.library_representation(title='H-Ras amino acid coverage')

.. image:: images/exported_images/hras_countspre.png
   :width: 500px
   :align: center
        
.. image:: images/exported_images/hras_countspre_aacoverage.png
   :width: 500px
   :align: center

Custom DNA list
~~~~~~~~~~~~~~~

Use a custom input DNA list. That way it does not matter if you are
using NNS or you have second order mutations. Create a list of variants
on your own, and the software will count the frequency of each of those
variants on the fastq file you provide as an input. In the example non
of the sequences we are specifying are found in the trimmed file, thus
there are 0% of useful reads.

.. code:: ipython3

    # Create your list of variants
    variants = [
        'acggaatataagctggtggtggtgggcgccggcggtgtgggcaagagtgcgctgaccat' +
        'ccagctgatccagaaccattttgtggacgaatacgaccccactatagaggattcctaccggaagcaggtgg' +
        'tcattgatggggagacgtgcctgttggacatcctg',
        'aaaaaatataagctggtggtggtgggcgccggcggtgtgggcaagagtgcgctgaccat' +
        'ccagctgatccagaaccattttgtggacgaatacgaccccactatagaggattcctaccggaagcaggtgg' +
        'tcattgatggggagacgtgcctgttggacatcctg',
        'tttttttataagctggtggtggtgggcgccggcggtgtgggcaagagtgcgctgaccat' +
        'ccagctgatccagaaccattttgtggacgaatacgaccccactatagaggattcctaccggaagcaggtgg' +
        'tcattgatggggagacgtgcctgttggacatcctg'
    ]
    
    # Count DNA variants in the fastq file
    input_file = '../data/hras.trimmed.fastq'
    
    variants, totalreads, usefulreads = mut.count_fastq(variants, input_file)
    
    # Evaluate how many variants in the fastq file were useful
    print(
        '{}/{} useful reads ({}%)'.format(
            str(usefulreads), str(totalreads),
            str(int(usefulreads / totalreads * 100))
        )
    )

Calculate enrichment scores
---------------------------

Methods and functions reviewed in this section:
    - :class:`mutagenesis_visualization.Screen`
    - :meth:`mutagenesis_visualization.Screen.heatmap`
    - :func:`mutagenesis_visualization.calculate_enrichment`

If you are performing a selection experiment, where you sequence your
library before and after selection, you will need to calculate the
enrichment score of each mutant. The function to do so is
``calculate_enrichment``. This function allows for different parameters
to tune how the data is muted and normalized.

In this example, we show two different ways of using ``calculate_enrichment``. Note that the parameters of choice will have a say on the final result. In the example, the tonality of red of the two heatmaps is slightly different. A more detailed explanation of the parameters can be found in :ref:`Normalizing datasets`.

.. code:: ipython3

    # Read counts from file (could be txt, csv, xlsx, etc...)
    df_counts_pre = pd.read_excel(
        '../data/hrasGAPGEF_counts.xlsx',
        'R1_before',
        skiprows=1,
        index_col='Codons',
        usecols='E:FN',
        nrows=32
    )
    
    df_counts_sel = pd.read_excel(
        '../data/hrasGAPGEF_counts.xlsx',
        'R1_after',
        skiprows=1,
        index_col='Codons',
        usecols='E:FN',
        nrows=32
    )

.. code:: ipython3

    # Ras parameters to create an object
    
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

.. code:: ipython3

    # Order of amino acids (from count_reads)
    aminoacids_NNS = list('AACDEFGGHIKLLLMNPPQRRRSSSTTVVWY*')
    
    # Different parameters can be used to calculate the enrichment scores. They are described in the implementation section
    
    # Zeroing using the median of the population, and not using stop codons to correct.
    frequencies = mut.calculate_enrichment(
        df_counts_pre.iloc[:, :54],
        df_counts_sel.iloc[:, :54],
        aminoacids=aminoacids_NNS,
        zeroing='population',
        how='median',
        norm_std=True,
        stopcodon=True,
        min_counts=25,
        min_countswt=100,
        mpop=2,
        mwt=2,
        infinite=3,
        std_scale=0.3
    )
    
    hras_example1 = mut.Screen(
        np.array(frequencies), hras_sequence, aminoacids, start_position, fillna,
        secondary
    )
    
    hras_example1.heatmap(title='Normal distribution zeroing', output_file=None)
    
    # Zeroing using the median of the population, and not using stop codons to correct.
    frequencies = mut.calculate_enrichment(
        df_counts_pre.iloc[:, :54],
        df_counts_sel.iloc[:, :54],
        aminoacids=aminoacids_NNS,
        zeroing='kernel',
        how='median',
        norm_std=True,
        stopcodon=True,
        min_counts=25,
        min_countswt=100,
        mpop=2,
        mwt=2,
        infinite=3,
        std_scale=0.15
    )
    
    hras_example2 = mut.Screen(
        np.array(frequencies), hras_sequence, aminoacids, start_position, fillna,
        secondary
    )
    
    hras_example2.heatmap(title='KDE zeroing', output_file=None)
    
    # Note that the two heatmaps look quite similar but the red tonality is slighly different. That is caused by
    # small differences in zeroing the data.

.. image:: images/exported_images/hras_tableenrichment.png
   :width: 450px
   :align: center

.. image:: images/exported_images/hras_zeronormal.png
   :width: 300px
   :align: center

.. image:: images/exported_images/hras_zerokernel.png
   :width: 300px
   :align: center

Assemble multiple sublibraries
------------------------------

Function reviewed in this section:
    - :func:`mutagenesis_visualization.assemble_avengers`

If you split your library into multiple pools, you can use ``assemble_avengers`` to use ``calculate_enrichment`` in an automated loop and return the assembled dataframe. To use this function, you need to import the data in an excel file in the same format as the provided in Example/hrasGAPGEF_counts.xlsx. Note that the parameters for normalization used in ``calculate_enrichment`` also apply here. See :ref:`Normalizing datasets` for more details.

.. code:: ipython3

    # Sheet that stores input/preselected counts within the Excel file
    sheet_pre = 'R1_before'
    # Sheet that stores output/selected counts within the Excel file
    sheet_post = 'R1_after'
    # Columns of each sublibrary. In this example, there are three pools.
    columns = ['F:BG', 'BH:DK', 'DL:FN']
    # Columns of the wt pools (optional)
    columns_wt = ['A', 'B', 'C']
    # Path were the excel file is stored.
    excel_path = '../data/hrasGAPGEF_counts.xlsx'
    # Parameter for pd.read_excel function
    nrows_pop = 32  # For nrows of the sublibrary
    nrows_wt = [50, 37, 57]  # For ncolumns of each of the three wild-type columns
    skiprows = 1  # Skip one row when reading the columns specified in the list `columns`
    
    # Normalization parameters also need to be specified. In here we
    # are using the default ones.
    
    # Call the function and return a df
    df = mut.assemble_avengers(
        excel_path,
        sheet_pre,
        sheet_post,
        columns,
        nrows_pop,
        nrows_wt,
        columns_wt,
        output_file=None
    )
    
    # The output is a combined dataframe

Combine MSA with enrichment scores
----------------------------------

Function and class reviewed in this section:
    - :class:`mutagenesis_visualization.Screen`
    - :func:`mutagenesis_visualization.msa_enrichment`

Function ``msa_enrichment`` will calculate the frequency of each substitution in an input MSA. The frequency of each substitution will be merged into the enrichment score dataframe. The function also calculates the Shannon entropy for each position in the protein. This function has been used to generate the data that is plotted in box plot and the ROC AUC charts :ref:`Correlation, PCA and ROC AUC`. We will first need to create the object.

.. code:: ipython3

    # Load enrichment scores
    hras_enrichment_RBD = np.genfromtxt('../data/HRas166_RBD.csv', delimiter=',')
    
    # Define protein sequence
    hras_sequence = 'MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDLAARTVESRQAQDLARSYGIPYIETSAKTRQGVEDAFYTLVREIRQHKLRKLNPPDESGPG'
    
    # Create object (more detail about this in plotting examples)
    hras_RBD = mut.Screen(hras_enrichment_RBD, hras_sequence)

Now we can get the frequency of each substituion in the MSA and the
Shannon entropy. You can use the example fasta file by loading
``fasta_dict = mut_fasta()`` and then ``path = fasta_dict['ras']``.

.. code:: ipython3

    # Calculate conservation score from MSA
    #path = '../data/Ras_family_trimmed.fasta'  # local file
    path = mut.demo_fasta()[
        'ras']  # Load example file (only if you are trying to reproduce the plots)
    
    # Calculate msa scores
    df_shannon, df_freq = mut.msa_enrichment(
        hras_RBD, path, start_position=1, threshold=0.1
    )
    
    # In the example, for position 2, in 3.63% of the cases there was an Ala.
    df_freq.head(5)

.. image:: images/exported_images/hras_table_msa.png
   :width: 300px
   :align: center

Note: The Shannon entropy is calculated using a script created by Joe R.
J. Healey from Warwick University. Could not find the script on Github
or Pypi so I included it in the package (shannon.py).
