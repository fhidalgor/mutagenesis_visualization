from typing import Dict, Union, List
from pandas.core.frame import DataFrame
import numpy as np

from mutagenesis_visualization.main.classes.counts import Counts
from mutagenesis_visualization.main.classes.screen import Screen
from mutagenesis_visualization.main.demo.demo_data import load_demo_datasets
from mutagenesis_visualization.main.process_data.count_reads import count_reads

def _return_hras_counts() -> Counts:
    """
    This method will generate a *Counts* object.
    """
    # H-Ras dna sequence
    hras_dnasequence: str = 'acggaatataagctggtggtggtgggcgccggcggtgtgggcaagagtgcgctgaccat' + 'ccagctgatccagaaccattttgtggacgaatacgaccccactatagaggattcctaccggaagcaggtgg' + 'tcattgatggggagacgtgcctgttggacatcctg'

    # Codons used to make the NNS library. I could also have used 'NNS' and the package will use the NNS codons
    codon_list: List[str] = [
        "GCC", "GCG", "TGC", "GAC", "GAG", "TTC", "GGC", "GGG", "CAC", "ATC", "AAG", "CTC", "CTG",
        "TTG", "ATG", "AAC", "CCC", "CCG", "CAG", "CGC", "CGG", "AGG", "TCC", "TCG", "AGC", "ACC",
        "ACG", "GTC", "GTG", "TGG", "TAC", "TAG"
    ]

    df_counts_pre, _ = count_reads(
        hras_dnasequence,
        "mutagenesis_visualization/data/hras.trimmed.fastq",
        codon_list,
        counts_wt=False,
        start_position=2,
        output_file=
    )

    return Counts(dataframe=df_counts_pre)

_return_hras_counts()
