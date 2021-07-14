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


def count_fastq(variants, input_file):
    """
    Count the frequency of variants in the input fastq file.

    Parameters
    -----------
    variants : ordered dict
        Contains each DNA sequence that you want to count from the fastq file.
        If your input is a list of strings, use the auxiliar function _initialize_ordereddict
        to convert it to an ordered dictionary.
        If you input a list, it will convert it to an ordered dict.

    input_file : str, default None
        Path and name of the fastq file (full name including suffix ".fastq").

    Returns
    --------
    variants : ordered dict
        Same input dictionary by now has the values updated with the counts.
    totalreads : int
        Total number of DNA chains that appear in the fastq file.
    usefulreads : int
        Total number of identified DNA chains. Calculated as the sum of all the key values.
    """
    # if variant input is not an ordered dict, convert to ordered dict
    if not (isinstance(variants, OrderedDict)):
        variants = _initialize_ordereddict(variants)

    # iterate over fastq file and count reads
    totalreads = 0
    for nuc in SeqIO.parse(str(input_file), "fastq"):
        totalreads += 1
        nucleicsequence = str(nuc.seq)
        if nucleicsequence in variants:
            variants[nucleicsequence] += 1
    usefulreads = np.nansum(list(variants.values()))
    return variants, totalreads, usefulreads
