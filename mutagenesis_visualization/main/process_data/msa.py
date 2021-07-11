"""

"""
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from collections import OrderedDict
from os import path
from pathlib import Path
from typing import Union
from scipy import stats
from logomaker import alignment_to_matrix

def msa_enrichment(self, path, start_position, threshold=0.01):
    """
    Generate a dataframe with the Shannon entropy by residue and the mean
    enrichment score, and a second dataframe with the frequency of each
    substitution and the enrichment score.

    Parameters
    -----------
    self : object from class *Screen*

    path : str
        Path where is located the fasta MSA that will be parsed. That MSA
        needs to have removed any insertions that are not present in the
        target sequence. For example, if a Ras ortholog has an extra
        amino acid at position 123, that needs to be removed from the
        aligment. Otherwise, everything will be shifted by 1 residue.

    start_position : int
        This is the position in the protein sequence of the first
        position in the MSA.

    threshold : float, default 0.01
        The conservation frequency for each amino acid subsitution will
        be binarized, and a threshold between 0-1 needs to be selected.

    Returns
    --------
    df_shannon: pandas dataframe
        Shannon entropy by residue and mean enrichment score by residue.

    df_freq : pandas dataframe
        Frequency of each susbsitution merged to the enrichment score.
    """
    # Read MSA
    msa, seq_lengths, index = code_utils._parseMSA(path, "fasta", 0)

    # Calculate Shannon entropy from alignment
    shannon_entropy = code_utils._shannon_entropy_list_msa(msa)

    # Merge enrichment scores and MSA conservation
    df_freq = _merge_msa_enrichment(
        self, _msa_to_df(msa), start_position, threshold
    )

    # Merge shannon and mean enrichment score
    df_shannon = _merge_shannon_enrichment(
        self, shannon_entropy, start_position
    )

    return df_shannon, df_freq


def _merge_shannon_enrichment(self, shannon_entropy, start_position):

    # Create df with shannon entropy by residue and average enrichment score by residue
    df_shannon = pd.DataFrame()
    df_shannon['Position'] = np.arange(
        start_position,
        len(shannon_entropy) + start_position
    )
    df_shannon['Shannon'] = shannon_entropy

    # group by enrichment scores
    df_enrichment = self.dataframe.groupby(by='Position', as_index=False).mean()

    # Merge Shannon with enrichment scores
    df_shannon = df_shannon.merge(df_enrichment, how='inner', on=['Position'])

    return df_shannon


def _flatten_msa(msa):
    '''Flatten an msa so each sequence is in one string'''
    msa_flattened = []
    for sequence in msa:
        msa_flattened.append(''.join(sequence))
    return msa_flattened


def _msa_to_df(msa, correctionfactor=1):
    '''
    Convert a msa from a fasta file into a df ready to plot with
    logomaker. Returns frequency
    '''
    # Flatten MSA
    msa_flattened = _flatten_msa(msa)

    # Make matrix
    df = alignment_to_matrix(msa_flattened)

    # Reindex
    df.index = np.arange(correctionfactor, len(df) + correctionfactor)

    # Return only common aa
    aminoacids = list('ACDEFGHIKLMNPQRSTVWY')

    # Normalize by the total number of counts
    df_final = df[aminoacids].copy()

    return df_final / df_final.sum(axis=1).max()


def _merge_msa_enrichment(self, df_msa, start_position, threshold):
    '''
    Merges msa conservation of each individual amino acid with the
    enrichment scores.
    '''

    # make a dataframe
    df = pd.DataFrame()

    # Create column with position and aminoacid label
    df['Position'] = np.ravel(
        [[i] * len(df_msa.T)
         for i in range(start_position,
                        len(df_msa) + start_position)]
    )
    df['Aminoacid'] = list(df_msa.columns) * len(df_msa)

    # Add conservation from MSA
    df['Conservation'] = list(df_msa.stack(dropna=False))

    # Merge with enrichment scores
    df_merged = self.dataframe.merge(
        df, how='inner', on=['Position', 'Aminoacid']
    )

    # Copycat conservation
    df_merged['Class'] = df_merged['Conservation']

    # Binarize conservation scores. 0 means not conserved
    df_merged.loc[df_merged['Conservation'] > threshold, 'Class'] = 1
    df_merged.loc[df_merged['Conservation'] <= threshold, 'Class'] = 0
    return df_merged
