#!/usr/bin/env python
# coding: utf-8

# In[2]:


from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from pathlib import Path
from typing import Union
import pandas as pd


# # Generate NNS primers

# In[3]:


def generate_primers(
    dna,
    start,
    end,
    output_file: Union[str, Path],
    codon='NNS',
    length_primer=15,
    tm=None,
    return_df=False
):
    """
    Generate primers for saturation mutagenesis.

    Parameters
    -----------
    dna : string
        DNA sequence containing the protein of study.
        The DNA sequence should also contain at least 15 base pairs before the
        starting ATG and 15 base pairs after the stop codon.

    start : string
        Beginning of the DNA sequence that will be mutageneized.
        For example, if you will start mutating the first methionine, copy a few more dna bases
        so the algorithm can identify it 'ATGACCAGC'.

    end : string
        The algorithm will stop creating primers once it reaches that base.
        For example, if you will stop mutating at the stop codon, copy a few more dna bases ie. 'TAAATGATT'.

    output_file : str
        File where the list of primers will be exported to. Only exports to excel.
        Example: 'path/primers.xlsx'.

    codon : str, default 'NNS'
        Degenerate codon that will be used to create the primers. Check idt's website for a list
        of all mixed bases and letter code (https://www.idtdna.com/pages/products/custom-dna-rna/mixed-bases).
        This parameter should contain 3 letters, although can contain more.

    length_primer: int, default 15
        Number of bases that the primers will have to each side of the mutated codon.
        Total primer length will be 2*length_primer+3.

    tm : int, default None
        Melting temperature in Celsius of the primers. Will override length_primer.
        If none, primers will have a total length of 2*length_primer+3

    return_df : boolean, default False
        If true, will export a dataframe with the primers.

    Returns
    --------
    df : pandas dataframe, optional
        Dataframe containing the primers.

    """
    # Transform to upper case
    dna = dna.upper()

    # Find the first and last codons
    start_codon = dna.find(start.upper())
    end_codon = dna.find(end.upper())

    # loop through DNA and make a list with fp and second list with rp

    label_fp = [
        'fp ' + str(i) for i in range(0, int((end_codon - start_codon) / 3))
    ]
    label_rp = [
        'rp ' + str(i) for i in range(0, int((end_codon - start_codon) / 3))
    ]
    forward_primers, reverse_primers = _create_primers_list(
        dna, start_codon, end_codon, codon, length_primer, tm
    )

    # Create dataframe
    dictionary = {
        'FP_label': label_fp, 'FP_seq': forward_primers, 'RP_label': label_rp,
        'RP_seq': reverse_primers
    }
    df = pd.DataFrame(dictionary)

    # Export dataframe
    if output_file:
        df.to_excel(Path(output_file), sheet_name='Primers', index=False)

    # Return dataframe
    if return_df:
        return df


# In[4]:


def _create_primers_list(dna, start_codon, end_codon, codon, length_primer, tm):
    '''Aux function to create list with fp and list with rp'''
    forward_primers = []
    reverse_primers = []
    for codonposition in range(start_codon, end_codon, 3):
        # Create fp, rp for that position
        fp, rp = _primerdesign(dna, codon, codonposition, length_primer, tm)
        # Append to list
        forward_primers.append(fp)
        reverse_primers.append(rp)
    return forward_primers, reverse_primers


def _reverse_complement(dna):
    '''aux function that uses biopython to calculate the reverse complement of a DNA string.
    Includes mixed-base code. More info in https://biopython.org/docs/1.75/api/Bio.Seq.html'''

    # Needs to be converted to str
    reverse_dna = str(Seq(dna).reverse_complement())

    return reverse_dna


def _primerdesign(dna, codon, codonposition, length_primer, tm):
    '''aux function to design the degenerate primers given a sequence and a codon position. 
    The length of the primer is fixed.

    Parameters
    -----------
    dna : string
        DNA sequence containing the protein of study. 
        The DNA sequence should also contain at least 15 base pairs before the 
        starting ATG and 15 base pairs after the stop codon.

    codon : str
        Degenerate codon that will be used to create the primers. Check idt's website for a list
        of all mixed bases and letter code (https://www.idtdna.com/pages/products/custom-dna-rna/mixed-bases).

    codonposition : int
        Position of the codon  to mutate with respect to the gene. 
        The first codon is 0 and the mutation occurs after the inputed codon number.

    length_primer : int, default 15
        Number of bases that the primers will have to each side of the mutated codon.
        Total primer length will be 2*length_primer+3.
        If using tm to generate primer length, then set length_primer = None. 

    tm : int, default None
        Melting temperature in Celsius of the primers. Will override length_primer.
        If none, primers will have a total length of 2*length_primer+3
        
    Returns
    ---------
    forward_primer, reverse_primer
    '''

    if tm:
        # loop until tm is achieved
        x = meltingT_fp = 6
        while meltingT_fp < tm:
            forward_primer = dna[(codonposition-x):codonposition] +                 codon + dna[(codonposition+3):(codonposition+x)]
            meltingT_fp = mt.Tm_NN(forward_primer)
            x += 1
    else:
        forward_primer = dna[(codonposition-length_primer):codonposition] +             codon + dna[(codonposition+3):(codonposition+length_primer+3)]

    reverse_primer = _reverse_complement(forward_primer)
    return forward_primer, reverse_primer


# # Generate Variants

# In[7]:


def create_variants(
    dna, codon_list, output_file: Union[str, Path], return_df=False
):
    '''
    Generate a list of all point mutants given a dna sequence and a list of codons.

    Parameters
    -----------
    dna : str, 
        Contains the DNA sequence of the allele of reference (usually wild-type).

    codon_list : list or str
        Input a list of the codons that were used to create point mutations. Example: ["GCC", "GCG", "TGC"].
        It is important to know that the order of the codon_list will determine the output order.

    output_file : str
        File where the list of primers will be exported to. Only exports to 'xlsx',
        'fasta', 'txt'.
        Example: 'path/sequences.xlsx'. 

    return_df : boolean, default False
        If true, will export a dataframe with the primers.

    Returns
    --------
    df : pandas dataframe, optional
        Dataframe containing the generated sequences.

    '''
    # Make upper case in case input was lower case
    dna = dna.upper()
    codon_list = [item.upper() for item in codon_list]

    # Generate list of variants
    seq_list = _enumerate_variants_2(dna, codon_list)

    # Make dataframe
    df = pd.DataFrame()
    df['Sequences'] = seq_list

    # Export to excel or fasta
    if output_file:
        if Path(output_file).suffix == '.xlsx':
            df.to_excel(Path(output_file), sheet_name='Variants', index=False)
        elif Path(output_file).suffix == '.fasta' or Path(output_file
                                                          ).suffix == '.txt':
            _list_to_fasta(seq_list, output_file)

    # Return dataframe
    if return_df:
        return df


def _enumerate_variants_2(dna, codon_list):
    '''
    Copy of _enumerate variants with slight changes. Does return the wild-type sequence
    as the first item of the list.
    '''
    # Create list with codons of sequence
    wtSeqList = [dna[i:i + 3] for i in range(0, len(dna), 3)]

    # List of sequences
    seq_list = [dna]

    # Loop over the dna sequence
    for position in range(0, len(wtSeqList)):
        for codon in (codon_list):
            variant = ''.join(wtSeqList[0:position]) +                 codon + ''.join(wtSeqList[position+1:])
            if (variant != dna):
                seq_list.append(variant)
    return seq_list


def _list_to_fasta(seq_list, output_file):
    '''Export list to fasta format'''

    # Open file
    ofile = open(str(output_file), "w")

    # Loop through list and write into file
    for i, seq in enumerate(seq_list):
        line = (">{}\n{}\n").format(str(i), seq)
        ofile.write(line)

    # Close file
    ofile.close()
    return

