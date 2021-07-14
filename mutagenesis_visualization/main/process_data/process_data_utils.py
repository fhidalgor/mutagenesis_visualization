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


def _codon_table():
    codontable = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M', 'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT':
        'T', 'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K', 'AGC': 'S', 'AGT': 'S', 'AGA': 'R',
        'AGG': 'R', 'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L', 'CCA': 'P', 'CCC': 'P', 'CCG':
        'P', 'CCT': 'P', 'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGA': 'R', 'CGC': 'R',
        'CGG': 'R', 'CGT': 'R', 'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V', 'GCA': 'A', 'GCC':
        'A', 'GCG': 'A', 'GCT': 'A', 'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E', 'GGA': 'G',
        'GGC': 'G', 'GGG': 'G', 'GGT': 'G', 'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S', 'TTC':
        'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L', 'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
        'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W'
    }
    return codontable


def _are_syn(codon1, codon2, codontable):
    """Determine if 2 codons are synonymous"""
    if codon1 == codon2:
        return 'wt codon'  # changed from False
    if _translate(codon1, codontable) is not _translate(codon2, codontable):
        return False
    return True


def _translate(seq, codontable):
    """Translate DNA sequence to protein."""
    # I forgot why I made this custom function instead of using a biopython function
    protein = ''
    if len(seq) % 3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i : i + 3]
            protein += codontable[codon]
    return protein


def _enumerate_variants(wtSeqList, codon_list, dna_sequence):
    """Will return an ordered dictionary with variants initialized to 0 counts"""
    # Create ordered dictionary
    variants = OrderedDict()

    # First instance that we see the wt?
    firstwtseq = False

    # Loop over codons
    for position in range(0, len(wtSeqList)):
        for codons in (codon_list):
            variant = ''.join(wtSeqList[0 : position]) + codons + ''.join(wtSeqList[position + 1 :])
            if (variant == dna_sequence):  # Store redundant wild-types
                if firstwtseq:
                    variant = 'wtSeq' + str(position)
                firstwtseq = True
            variants[variant] = 0
    return variants


def _initialize_ordereddict(list_variants):
    """
    Will return an ordered dictionary with variants initialized to 0 counts.
    Here the user specifies the variants as a list.

    This function should be used when you want to use _count_fastq

    """

    # Create normal dictionary
    dictionary = dict(zip(list_variants, np.zeros(len(list_variants))))

    # Create ordered dictionary
    variants = OrderedDict(dictionary)

    return variants


def _get_enrichment(
    input_lib, output_lib, input_stopcodon, output_stopcodon, min_counts, stopcodon, infinite
):
    """Calculate log10 enrichment scores from input and output counts"""
    # Copy data and replace low counts by np.nan
    input_lib = np.copy(input_lib.astype(float))
    output_lib = np.copy(output_lib.astype(float))
    input_lib[input_lib < min_counts] = np.nan

    # Stop codon correction
    if stopcodon:
        output_lib = _stopcodon_correction(input_lib, output_lib, input_stopcodon, output_stopcodon)

    # log10 of library and replace infinite values. This will potentially divide by zero.
    with np.errstate(divide='ignore'):
        counts_log10_ratio = _replace_inf(np.log10(output_lib / input_lib), infinite)

    return counts_log10_ratio


def _stopcodon_correction(input_lib, output_lib, input_stopcodon, output_stopcodon):
    """This aux function will take as an input the counts for pre and post selection (and also for wT subset),
    and will return the corrected output counts"""

    # calculate stop codons frequencies
    frequency_stopcodons = output_stopcodon / input_stopcodon

    # MAD filtering
    frequency_stopcodons_filtered = _MAD_filtering(frequency_stopcodons, m=2)
    median_frequency = np.nanmedian(frequency_stopcodons_filtered)

    # subtract to output counts
    output_lib_corr = output_lib - input_lib * median_frequency

    # eliminate negative values so they wont get turned into np.nan
    output_lib_corr[output_lib_corr < 0] = 0

    return output_lib_corr


def _MAD_filtering(data, m=2):
    """This aux function will take a numpy array, calculate median and MAD,
    and filter the data removing outliers"""

    # turn data into df to do mad calculations
    df = pd.DataFrame(np.array(data), columns=['Data'])
    median = df['Data'].median(axis=0)
    mad = df['Data'].mad(axis=0)
    df['Abs_Dev'] = np.abs(data - median) / mad

    # filter values m times away from median, by default m = 2
    df['Abs_Dev'].mask(df['Abs_Dev'] > m, inplace=True)  # mask values
    df.dropna(how='any', inplace=True)  # eliminte NANs

    return df['Data'].to_numpy()


def _replace_inf(array, infinite):
    """Replace values over a threshold with a min or max value"""
    np.warnings.filterwarnings('ignore')
    array[array == -np.inf] = -infinite
    array[array < -infinite] = -infinite
    array[array == +np.inf] = +infinite
    array[array > +infinite] = +infinite
    return array


def _group_byaa(df, aminoacids):
    """Group different codons that are synonymous"""
    # copy df
    df = df.copy()

    # Set up amino acid column
    df['Aminoacid'] = aminoacids

    # Group by mean
    df = df.groupby(as_index=True, by='Aminoacid', sort=False).mean()
    return df


def _nanmode(data):
    """
    Input is wt log enrichments, and return the mode of the histogram
    (aka the x coordinate at which y is max).

    """

    # Copy data
    data = np.copy(data)
    # Remove NaN values
    data_corrected = data[np.invert(np.isnan(data))]
    # Adjust kernel
    kernel_processed_data = stats.gaussian_kde(data_corrected)
    # Find mode
    indexmax = np.where(
        kernel_processed_data(data_corrected) == kernel_processed_data(data_corrected).max()
    )
    # Return mean in case there are two x values with equal y-axis height
    return data_corrected[indexmax].mean()


# corrects the mutagenesis data and returns the height of the peak
def _kernel_correction(data, aminoacids, cutoff=2):
    """input the library matrix, returns the corrected version. I set to 0 the max of the peak of the normal dist
    ignores stop codons. Not used for dataframes, only numpy arrays"""

    # Get data into right format
    data_corrected, kernel_processed_data = _kernel_datapreparation(data, cutoff)

    # Find max of kernel peak
    indexmax = np.where(
        kernel_processed_data(data_corrected) == kernel_processed_data(data_corrected).max()
    )

    # Normalize the max of peak os it has an x = 0
    data_final = data - data_corrected[indexmax].mean()

    # find std of kernel. It uses the already max peak x=0 normalized data
    data_final_flatten, data_final_kernel_processed_data = _kernel_datapreparation(
        data_final, cutoff
    )
    std = _kernel_std(data_final_flatten, data_final_kernel_processed_data)

    return data_final, std


def _kernel_datapreparation(data, cutoff):
    """
    This function will copy the data, eliminate stop codon, eliminate values lower than -1,
    flatten and eliminate np.nan. Will return the data in that format + the adjusted kernel

    """

    # Eliminate stop codon
    data_corrected = np.array(data.drop('*', errors='ignore').copy())

    # Eliminate values lower than -1
    data_corrected = data_corrected[(data_corrected >= -cutoff) & (data_corrected <= cutoff)]

    # Get rid of np.nan values and convert matrix into 1d matrix
    data_corrected = data_corrected[np.invert(np.isnan(data_corrected))]

    # Adjust gaussian kernel
    kernel_processed_data = stats.gaussian_kde(data_corrected)

    return data_corrected, kernel_processed_data


def _kernel_std(data, kernel):
    """
    Input the library matrix (and wont count stop codon), and will return the std of the normal distribution.
    To calculate the std, it will find the FWHM and divide by 2.355
    https://en.wikipedia.org/wiki/Full_width_at_half_maximum
    The algorithm will give back the min std between both sides of the peak.

    """

    # find ymax and the x value of the max height
    y_max = kernel(data).max()
    index_y_max = np.where(kernel(data) == y_max)
    x_ymax = data[index_y_max].mean()

    # find the two x value of ymax/2. One at each side of the center. l of left and r of right
    # so I can select only the positive side of the distribution
    y_temp = kernel(data)

    # left side
    y_hw_l = (min(y_temp[data < 0], key=lambda x: abs(x - y_max / 2)))
    index_yhw_l = np.where(kernel(data) == y_hw_l)
    x_yhw_l = data[index_yhw_l].mean()

    # right side
    y_hw_r = (min(y_temp[data > 0], key=lambda x: abs(x - y_max / 2)))
    index_yhw_r = np.where(kernel(data) == y_hw_r)
    x_yhw_r = data[index_yhw_r].mean()

    # calculate half width at half maximum
    hwhm_l = abs(x_yhw_l - x_ymax)
    hwhm_r = abs(x_yhw_r - x_ymax)

    # calculate std from fwhm
    std_l = hwhm_l / ((2 * np.log(2))**0.5)
    std_r = hwhm_r / ((2 * np.log(2))**0.5)

    return min(std_l, std_r)


def _array_to_df_enrichments(lib, aminoacids):
    """
    aux function to transform array in df with index of amino acids.

    """

    df = pd.DataFrame(index=aminoacids, data=lib)
    return df.astype(float)


def _read_counts(
    excel_path, sheet_pre, sheet_post, columns, nrows_pop, nrows_wt, columns_wt=None, skiprows=1
):
    """Aux"""
    # Create dictionary with data. Loading 3 replicates, each of them is divided into 3 pools
    list_pre, list_sel, list_pre_wt, list_sel_wt = ([] for i in range(4))

    # Read counts from excel
    replicates = np.arange(0, len(sheet_pre))
    for column, column_wt, nrow_wt, rep in zip(columns, columns_wt, nrows_wt, replicates):
        # Pre counts
        list_pre.append(
            pd.read_excel(
                excel_path, sheet_pre, skiprows=skiprows, usecols=column, nrows=nrows_pop
            )
        )
        # Sel counts
        list_sel.append(
            pd.read_excel(
                excel_path, sheet_post, skiprows=skiprows, usecols=column, nrows=nrows_pop
            )
        )
        if columns_wt is None:
            list_pre_wt.append(None)
            list_sel_wt.append(None)
        else:
            # Pre counts wild-type alleles
            list_pre_wt.append(
                pd.read_excel(excel_path, sheet_pre, usecols=column_wt, nrows=nrow_wt)
            )
            # Sel counts wild-type alleles
            list_sel_wt.append(
                pd.read_excel(excel_path, sheet_post, usecols=column_wt, nrows=nrow_wt)
            )
    return list_pre, list_sel, list_pre_wt, list_sel_wt


def _assemble_list(
    list_pre,
    list_sel,
    list_pre_wt,
    list_sel_wt,
    aminoacids,
    zeroing,
    how,
    norm_std,
    stopcodon,
    min_counts,
    min_countswt,
    std_scale,
    mpop,
    mwt,
    infinite,
    output_file: Union[None, str, Path] = None
):
    """
    gets the output from _read_counts and assembles the sublibraries

    """

    enrichment_lib = []

    for pre, sel, pre_wt, sel_wt in zip(list_pre, list_sel, list_pre_wt, list_sel_wt):
        # log 10
        enrichment_log10 = calculate_enrichment(
            pre, sel, pre_wt, sel_wt, aminoacids, zeroing, how, norm_std, stopcodon, min_counts,
            min_countswt, std_scale, mpop, mwt, infinite
        )
        # Store in list
        enrichment_lib.append(enrichment_log10)

    # Concatenate sublibraries
    df = pd.concat(enrichment_lib, ignore_index=True, axis=1)
    return df
