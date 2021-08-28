"""
This module tests the process data methods and functions.
"""
import tempfile
from typing import List
from pandas import Index
from pandas.core.frame import DataFrame
from pandas.testing import assert_frame_equal

from mutagenesis_visualization.main.process_data.count_reads import count_reads
from mutagenesis_visualization.main.utils.data_paths import TEST_SHORT_FASTQ


def test_count_reads() -> None:
    """
    Method to test function count reads.
    """

    # Need to remove firstwtseq from _enumerate_variants

    def mut_assert_df_equal(left: DataFrame, right: DataFrame) -> None:
        assert_frame_equal(
            left,
            right,
            check_dtype=False,
            check_index_type=False,
            check_column_type=False,
            check_frame_type=False,
            check_names=False,
            check_like=True,
        )

    # Create dataframe
    codon_list: List[str] = [
        'AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 'AGC', 'AGG', 'AGT', 'ATA',
        'ATC', 'ATG', 'ATT', 'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC',
        'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT', 'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG',
        'GCT', 'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT', 'TAA', 'TAC', 'TAG', 'TAT',
        'TCA', 'TCC', 'TCG', 'TCT', 'TGA', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT'
    ]
    index = Index(codon_list)
    column_counts = Index([2])
    column_wt = Index(["Position", "Codon", "Aminoacid", "Counts"])
    values_cct = values_atc = [0] * 23 + [1] + [0] * 40

    # Test ATG
    expected_atc_counts = DataFrame(values_atc, index=index, columns=column_counts)
    expected_atc_wt = DataFrame([], columns=column_wt)
    atg_counts, atg_wt = count_reads("atg", TEST_SHORT_FASTQ, codon_list)
    # return atg_counts
    mut_assert_df_equal(atg_counts, expected_atc_counts)
    mut_assert_df_equal(atg_wt, expected_atc_wt[["Counts"]])

    # Test CCT
    expected_cct_counts = DataFrame(values_cct, index=index, columns=column_counts)
    expected_cct_wt = DataFrame([[2, 'CCA', 'P', 0], [2, 'CCC', 'P', 0], [2, 'CCG', 'P', 0]],
                                index=[20, 21, 22],
                                columns=column_wt)
    cct_counts, cct_wt = count_reads("cCt", TEST_SHORT_FASTQ, codon_list)
    mut_assert_df_equal(cct_counts, expected_cct_counts)
    mut_assert_df_equal(cct_wt, expected_cct_wt[["Counts"]])

    # Test CCT when not in codon list
    index = Index([
        "GCC", "GCG", "TGC", "GAC", "GAG", "TTC", "GGC", "GGG", "CAC", "ATC", "AAG", "CTC", "CTG",
        "TTG", "ATG", "AAC", "CCC", "CCG", "CAG", "CGC", "CGG", "AGG", "TCC", "TCG", "AGC", "ACC",
        "ACG", "GTC", "GTG", "TGG", "TAC", "TAG"
    ])
    values_cct = [0] * 32
    expected_cct_counts = DataFrame(values_cct, index=index, columns=column_counts)
    expected_cct_wt = DataFrame([[2, 'CCC', 'P', 0], [2, 'CCG', 'P', 0]],
                                index=[16, 17],
                                columns=column_wt)

    # Now check with NNK
    # Create a temporary directory using the context manager
    with tempfile.TemporaryDirectory() as tmpdirname:
        cct_counts, cct_wt = count_reads(
            dna_sequence="cCt",
            input_file=TEST_SHORT_FASTQ,
            codon_list='NNK',
            output_file=tmpdirname + '/counts.xlsx',
            full=True,
        )

    values_cct = [0] * 17 + [1] + [0] * 14
    index: List[str] = [
        'GCG', 'GCT', 'TGT', 'GAT', 'GAG', 'TTT', 'GGG', 'GGT', 'CAT', 'ATT', 'AAG', 'CTG', 'CTT',
        'TTG', 'ATG', 'AAT', 'CCG', 'CCT', 'CAG', 'AGG', 'CGG', 'CGT', 'AGT', 'TCG', 'TCT', 'ACG',
        'ACT', 'GTG', 'GTT', 'TGG', 'TAT', 'TAG'
    ]
    expected_cct_counts = DataFrame(values_cct, index=index, columns=column_counts)
    expected_cct_wt = DataFrame([[2, 'CCG', 'P', 0]], index=[16], columns=column_wt)
    mut_assert_df_equal(cct_counts, expected_cct_counts)
    mut_assert_df_equal(cct_wt, expected_cct_wt[["Counts"]])
