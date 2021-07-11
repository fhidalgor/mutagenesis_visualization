"""
This module contains Function to calcuate the Shannon's entropy per
alignment column.
"""
import math

##################################################################
# Function to calcuate the Shannon's entropy per alignment column
# H=-\sum_{i=1}^{M} P_i\,log_2\,P_i (http://imed.med.ucm.es/Tools/svs_help.html)
# Gaps and N's are included in the calculation
##################################################################


def shannon_entropy(list_input: list) -> list:
    """
    Calculate Shannon's Entropy per column of the alignment.
    """
    unique_base = set(list_input)
    M = len(list_input)
    entropy_list = []
    # Number of residues in column
    for base in unique_base:
        n_i = list_input.count(base)  # Number of residues of type i
        P_i = n_i / float(
            M
        )  # n_i(Number of residues of type i) / M(Number of residues in column)
        entropy_i = P_i * (math.log(P_i, 2))
        entropy_list.append(entropy_i)

    sh_entropy = -(sum(entropy_list))

    return sh_entropy


def shannon_entropy_list_msa(alignment) -> list:
    """
    Calculate Shannon Entropy across the whole MSA.
    """
    shannon_entropy_list = []
    for col_no in range(len(list(alignment[0]))):
        list_input = list(alignment[:, col_no])
        shannon_entropy_list.append(shannon_entropy(list_input))

    return shannon_entropy_list
