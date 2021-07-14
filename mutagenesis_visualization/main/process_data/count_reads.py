"""

"""
import numpy as np
import pandas as pd
import copy
from Bio import SeqIO
from Bio.Seq import Seq
from collections import OrderedDict
from os import path
from pathlib import Path
from typing import Union
from scipy import stats
from logomaker import alignment_to_matrix


def count_reads(
    dna_sequence,
    input_file: Union[str, Path],
    codon_list='NNS',
    counts_wt=True,
    start_position=2,
    output_file: Union[None, str, Path] = None,
    full=False
):
    """
    Process a trimmed fastq file containing DNA reads and returns the counts of
    each DNA sequence specified by the user.

    Parameters
    -----------
    dna_sequence : str,
        Contains the DNA sequence of the allele of reference (usually wild-type).

    input_file : str, default None
        Path and name of the fastq file (full name including suffix ".fastq").

    codon_list : list or str, default 'NNS'
        Input a list of the codons that were used to create point mutations.
        Example: ["GCC", "GCG", "TGC"].
        If the library was built using NNS and NNK codons, it is enough to input
        'NNS' or 'NNK' as a string. It is important to know that the order of
        the codon_list will determine the output order.

    counts_wt : boolean, default True
        If true it will add the counts to the wt allele. If false, it will set it up to np.nan.

    start_position : int, default 2
        First position in the protein sequence that will be used for the first column of the
        array. If a protein has been mutated only from residue 100-150, then if start_position = 100,
        the algorithm will trim the first 99 amino acids in the input sequence. The last
        residue will be calculated based on the length of the input array. We have set the default value to 2
        because normally the Methionine in position 1 is not mutated.

    output_file : str, default None
        If you want to export the generated files, add the path and name of the file without suffix.
        Example: 'path/filename.xlsx'.

    full: bool, optional
        Switch determining nature of return value.
        When it is False (the default) just the reads are
        returned, when True diagnostic information from the
        fastq analysis is also returned.

    Returns
    --------
    df_counts : dataframe
        Dataframe with the counts for each point mutant.

    wt_counts : list
        List of the counts for each for each DNA sequence that codes for the wild-type protein.

    useful_reads : str
        Present only if `full` = True. Contains the useful reads.

    """
    # Assert messages
    assert len(dna_sequence) % 3 == 0, 'The dna_sequence length is not a multiple of 3'

    # Make upper case in case input was lower case
    dna_sequence = dna_sequence.upper()
    if isinstance(codon_list, str):  #check if codon_list is a String
        codon_list = codon_list.upper()
    else:
        codon_list = [item.upper() for item in codon_list]

    # Create list with codons of sequence
    wtSeqList = [dna_sequence[i : i + 3] for i in range(0, len(dna_sequence), 3)]

    # codon_list
    if codon_list == 'NNS':
        codon_list = [
            "GCC", "GCG", "TGC", "GAC", "GAG", "TTC", "GGC", "GGG", "CAC", "ATC", "AAG", "CTC",
            "CTG", "TTG", "ATG", "AAC", "CCC", "CCG", "CAG", "CGC", "CGG", "AGG", "TCC", "TCG",
            "AGC", "ACC", "ACG", "GTC", "GTG", "TGG", "TAC", "TAG"
        ]
    elif codon_list == 'NNK':
        codon_list = [
            'GCG', 'GCT', 'TGT', 'GAT', 'GAG', 'TTT', 'GGG', 'GGT', 'CAT', 'ATT', 'AAG', 'CTG',
            'CTT', 'TTG', 'ATG', 'AAT', 'CCG', 'CCT', 'CAG', 'AGG', 'CGG', 'CGT', 'AGT', 'TCG',
            'TCT', 'ACG', 'ACT', 'GTG', 'GTT', 'TGG', 'TAT', 'TAG'
        ]

    # Enumerate variants
    variants = _enumerate_variants(wtSeqList, codon_list, dna_sequence)

    # Count variant frequency
    variants, totalreads, usefulreads = count_fastq(variants, input_file)

    # Convert to df
    wtProtein = Seq(dna_sequence).translate()
    df = pd.DataFrame()
    df['Position'] = np.ravel(
        [[pos] * len(codon_list)
         for pos in np.arange(start_position,
                              len(wtProtein) + start_position).astype(int)]
    )
    df['Codon'] = codon_list * len(wtProtein)
    df['WTCodon'] = np.ravel([[codon] * len(codon_list) for codon in wtSeqList])
    df['Aminoacid'] = np.ravel([[aa] * len(codon_list) for aa in wtProtein])
    df['SynWT'] = df.apply(lambda x: _are_syn(x['Codon'], x['WTCodon'], _codon_table()), axis=1)
    df['Counts'] = list(variants.values())

    if counts_wt:
        try:  # try is to fix the Bug Che discovered
            df.loc[df['Codon'] == df['WTCodon'], 'Counts'] = variants[dna_sequence]
        except:
            pass
    else:
        df.loc[df['Codon'] == df['WTCodon'], 'Counts'] = np.nan

    # Pivot table and reindex
    df_counts = df.pivot_table(values='Counts', index='Codon', columns=['Position'], dropna=False)
    df_counts = df_counts.reindex(index=codon_list)

    # Get WT counts syn. Added or operator so also chooses WT codon
    df_wt = df.loc[(df['SynWT'] == True) | (df['SynWT'] == 'wt codon')][
        [  # perhaps I need to remove this again
            'Counts'  # removed
        ]]

    # Export files
    if output_file:
        # Create a Pandas Excel writer using XlsxWriter as the engine.
        writer = pd.ExcelWriter(str(output_file), engine='xlsxwriter')

        # export
        df_counts.to_excel(writer, sheet_name='Counts', index=True)
        df_wt.to_excel(writer, sheet_name='WT', index=False)

        # Close the Pandas Excel writer and output the Excel file.
        writer.save()

    if full:
        # Print total reads
        percent_useful = usefulreads / totalreads * 100
        return df_counts, df_wt, f"{usefulreads}/{totalreads} useful reads ({percent_useful:.1f}%)"
    else:
        return df_counts, df_wt
