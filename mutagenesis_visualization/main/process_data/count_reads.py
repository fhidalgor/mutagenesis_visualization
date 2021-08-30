"""
This module contains the count dna reads function.
"""

from typing import List, Union, Tuple
from pathlib import Path
import numpy as np
from pandas import ExcelWriter
from pandas.core.frame import DataFrame
from Bio.Seq import Seq

from mutagenesis_visualization.main.process_data.count_fastq import count_fastq
from mutagenesis_visualization.main.process_data.process_data_utils import (
    enumerate_variants, generate_codon_table, are_syn
)


def count_reads(
    dna_sequence: str,
    input_file: Union[str, Path],
    codon_list: Union[List[str], str] = 'NNS',
    counts_wt: bool = True,
    start_position: int = 2,
    output_file: Union[None, str, Path] = None,
    full: bool = False
) -> Tuple[DataFrame, DataFrame]:
    """
    Process a trimmed fastq file containing DNA reads and returns the
    counts of each DNA sequence specified by the user.

    Parameters
    -----------
    dna_sequence : str,
        Contains the DNA sequence of the allele of reference (usually
        wild-type).

    input_file : str, default None
        Path and name of the fastq file (full name including suffix
        ".fastq").

    codon_list : list or str, default 'NNS'
        Input a list of the codons that were used to create point mutations.
        Example: ["GCC", "GCG", "TGC"].
        If the library was built using NNS and NNK codons, it is enough
        to input 'NNS' or 'NNK' as a string. It is important to know
        that the order of the codon_list will determine the output order.

    counts_wt : boolean, default True
        If true it will add the counts to the wt allele. If false, it
        will set it up to np.nan.

    start_position : int, default 2
        First position in the protein sequence that will be used for the
        first column of the array. If a protein has been mutated only
        from residue 100-150, then if start_position = 100, the algorithm
        will trim the first 99 amino acids in the input sequence. The last
        residue will be calculated based on the length of the input array.
        We have set the default value to 2 because normally the Methionine
        in position 1 is not mutated.

    output_file : str, default None
        If you want to export the generated files, add the path and name
        of the file without suffix.
        Example: 'path/filename.xlsx'.

    full: bool, optional
        Switch determining nature of return value. When it is False
        (the default) just the reads are returned, when True diagnostic
        information from the fastq analysis is also returned.

    Returns
    --------
    df_counts : dataframe
        Dataframe with the counts for each point mutant.

    wt_counts : list
        List of the counts for each for each DNA sequence that codes for
        the wild-type protein.

    useful_reads : str
        Present only if `full` = True. Contains the useful reads.    """
    # Assert messages
    assert len(dna_sequence) % 3 == 0, 'The dna_sequence length is not a multiple of 3'

    # Make upper case in case input was lower case
    dna_sequence = dna_sequence.upper()
    # Create list with codons of sequence
    wtseq_list: List[str] = [dna_sequence[i : i + 3] for i in range(0, len(dna_sequence), 3)]

    if isinstance(codon_list, str):  #check if codon_list is a String
        codon_list = codon_list.upper()
    else:
        codon_list = [item.upper() for item in codon_list]

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
    assert not isinstance(codon_list, str), "Error when using codon list"
    variants = enumerate_variants(wtseq_list, codon_list, dna_sequence)
    assert dna_sequence in list(variants.keys()), "DNA sequence not compatible with codon list."

    # Count variant frequency
    variants, totalreads, usefulreads = count_fastq(list(variants.keys()), input_file)
    # Convert to df_output
    wt_protein: Seq = Seq(dna_sequence).translate()
    df_output: DataFrame = DataFrame()
    df_output['Position'] = np.ravel(
        [[pos] * len(codon_list)
         for pos in np.arange(start_position,
                              len(wt_protein) + start_position).astype(int)]
    )
    df_output['Codon'] = codon_list * len(wt_protein)
    df_output['WTCodon'] = np.ravel([[codon] * len(codon_list) for codon in wtseq_list])
    df_output['Aminoacid'] = np.ravel([[aa] * len(codon_list) for aa in wt_protein])
    df_output['SynWT'] = df_output.apply(
        lambda x: are_syn(x['Codon'], x['WTCodon'], generate_codon_table()), axis=1
    )
    df_output['Counts'] = list(variants.values())

    if counts_wt:
        # try is to fix the Bug Che discovered
        df_output.loc[df_output['Codon'] == df_output['WTCodon'], 'Counts'] = variants[dna_sequence]
    else:
        df_output.loc[df_output['Codon'] == df_output['WTCodon'], 'Counts'] = np.nan

    # Pivot table and reindex
    df_counts = df_output.pivot_table(
        values='Counts', index='Codon', columns=['Position'], dropna=False
    )
    df_counts = df_counts.reindex(index=codon_list)

    # Get WT counts syn. Added or operator so also chooses WT codon
    df_wt = df_output.loc[df_output['SynWT'] == True][['Counts']]  # pylint: disable=singleton-comparison

    # Export files
    if output_file:
        # Create a Pandas Excel writer using XlsxWriter as the engine.
        writer = ExcelWriter(str(output_file), engine='xlsxwriter')  # pylint: disable=abstract-class-instantiated

        # export
        df_counts.to_excel(writer, sheet_name='Counts', index=True)
        df_wt.to_excel(writer, sheet_name='WT', index=False)

        # Close the Pandas Excel writer and output the Excel file.
        writer.save()

    if full:
        percent_useful: float = usefulreads / totalreads * 100
        print(f"{usefulreads}/{totalreads} useful reads ({percent_useful:.1f}%)")
    return df_counts, df_wt
