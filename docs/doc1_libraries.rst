Design DNA libraries
====================

In this section we will generate the primers that are used to do
saturation mutagenesis on proteins (ie. NNS primers).

We will also generate each possible point mutant sequence and export it
to a Fasta file, which can be useful if you use Twist Bioscience to
generate your site saturation library.

.. code:: ipython3

    try:
        import mutagenesis_visualization as mut
    except ModuleNotFoundError:  # This step is only for when I run the notebooks locally
        import sys
        sys.path.append('../../')
        import mutagenesis_visualization as mut

Design primers
--------------

Now we will define the dna sequence, the beginning and end of the
mutable part.

.. code:: ipython3

    # DNA
    dna = 'TGTACAGTAATACAAGGGGTGTTATGGAAAAAATTATGCCGGAAGAAGAATACAGCGAATTTAAAGAACTGATTCTGCAGAAGGAACTGCACGTGGTGTATGCACTGAGCCACGTGTGTGGCCAGGATCGTACCCTGCTGGCCAGTATCTTACTGCGCATCTTTCTGCACGAGAAGCTGGAGAGCCTGTTACTGTGCACACTGAACGATCGCGAGATCAGCATGGAAGATGAAGCCACCACCCTGTTCCGCGCAACAACCCTGGCCAGCACCCTGATGGAGCAGTATATGAAAGCCACCGCCACCCAGTTCGTGCATCATGCCCTGAAAGATAGCATTTTAAAAATTATGGAAAGCAAACAGAGCTGCGAACTGAGCCCGAGCAAGCTGGAGAAAAACGAGGACGTGAACACCAACCTGACCCACCTGCTGAACATTCTGAGCGAACTGGTGGAAAAAATCTTTATGGCAAGCGAAATCCTGCCTCCGACCCTGCGTTACATCTACGGCTGCCTGCAGAAGAGCGTGCAGCATAAATGGCCGACCAATACCACCATGCGCACACGTGTGGTGAGCGGTTTTGTGTTCCTGCGTCTGATCTGCCCGGCAATCCTGAACCCGCGCATGTTCAACATCATTAGCGACAGCCCGAGTCCTATCGCAGCACGTACCCTGATCCTGGTGGCAAAAAGCGTGCAAAATCTGGCCAACCTGGTGGAATTTGGCGCCAAAGAGCCGTACATGGAAGGCGTGAATCCGTTTATCAAAAGTAACAAACATCGCATGATCATGTTCCTGGACGAACTGGGCAACGTTCCGGAACTGCCGGATACAACCGAACATAGTCGCACAGACCTGAGTCGTGACCTGGCCGCCCTGCATGAAATCTGCGTGGCCCATAGCGATGAGCTGCGCACACTGAGCAACGAGCGTGGCGCCCAGCAGCACGTGCTGAAGAAACTGCTGGCCATTACCGAACTGCTGCAACAAAAGCAGAACCAGTACACCAAAACCAACGACGTGCGTtatccgtatgatgtgccggattatgcgTAAccatcacttggctagaggcatc'
    
    # Start of protein. Note 'ATG' codon
    start = ('ATGGAAAAAATTATGCCGGAAGAA')
    
    # The 'tat' codon will be the first codon that is not mutated
    end = ('tatccgtatgatgtgccggattatgcg')

Set all primers to have the same base pair length.

.. code:: ipython3

    df_primers = mut.generate_primers(
        dna,
        start,
        end,
        output_file=None,
        codon='NNS',
        length_primer=15,
        return_df=True
    )

Set all primers to have the same melting temperature.

.. code:: ipython3

    df_primers_tm = mut.generate_primers(
        dna, start, end, output_file=None, codon='NNS', tm=60, return_df=True
    )

If you just want to export the file to excel:

.. code:: ipython3

    mut.generate_primers(
        dna,
        start,
        end,
        output_file='primers.xlsx',
        codon='NNS',
        tm=60,
        return_df=False
    )

.. image:: images/exported_images/primers.png
   :width: 450px
   :align: center

Design site-saturation sequences
--------------------------------

Define dna sequence and the list of codons that we want to use to
generate the mutants.

.. code:: ipython3

    # list of codons we want to use
    codon_list = ["GCC", "GCG", "TGC", "GAC", "GAG", "TTC"]
    # DNA sequence we are going to use as the template
    dna = 'ATGGCCGTGGGGTGTTATGGATGTACAGTAATACAAGGGGTGTTATGGAAAAAATTATGCCGGAAGAAGAATACAGCGAATTTAAAG'

Get a dataframe with the sequences:

.. code:: ipython3

    df = mut.create_variants(dna, codon_list, output_file=None, return_df=True)

If you just want to export the file to fasta:

.. code:: ipython3

    mut.create_variants(dna, codon_list, output_file='sequences.fasta')

.. image:: images/exported_images/fasta.png
   :width: 300px
   :align: center

If you just want to export the file to excel:

.. code:: ipython3

    mut.create_variants(dna, codon_list, output_file='sequences.xlsx')
