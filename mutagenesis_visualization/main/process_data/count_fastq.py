"""
This module contains the function to count dna variants in fastq files.
"""
from typing import Tuple, Union, List
from pathlib import Path
import numpy as np
from Bio import SeqIO

from mutagenesis_visualization.main.process_data.process_data_utils import initialize_ordered_dict


def count_fastq(variants: List[str], input_file: Union[str, Path]) -> Tuple[dict, int, int]:
    """
    Count the frequency of variants in the input fastq file.

    Parameters
    -----------
    variants : list

    input_file : str, default None
        Path and name of the fastq file (full name including suffix ".fastq").

    Returns
    --------
    variants : ordered dict
        Same input dictionary by now has the values updated with the counts.
    totalreads : int
        Total number of DNA chains that appear in the fastq file.
    usefulreads : int
        Total number of identified DNA chains. Calculated as the sum of
        all the key values.
    """
    # if variant input is not an ordered dict, convert to ordered dict
    variants_dict = initialize_ordered_dict(variants)

    # iterate over fastq file and count reads
    total_reads: int = 0
    for nuc in SeqIO.parse(str(input_file), "fastq"):
        total_reads += 1
        nucleicsequence = str(nuc.seq)
        if nucleicsequence in variants_dict:
            variants_dict[nucleicsequence] += 1
    usefulreads: int = np.nansum(list(variants_dict.values()))
    return variants_dict, total_reads, usefulreads
