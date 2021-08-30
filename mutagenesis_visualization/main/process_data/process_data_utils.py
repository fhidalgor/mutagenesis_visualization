"""
Utils for process data.
"""
from typing import Any, Dict, Tuple, Union, List
from collections import OrderedDict
import numpy as np
from numpy import typing as npt
from pandas.core.frame import DataFrame
from scipy.stats import gaussian_kde


def generate_codon_table() -> Dict[str, str]:
    """
    Returns codon table dictionary.
    """
    return {
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


def are_syn(codon1: str, codon2: str, codon_table: Dict[str, str]) -> Union[bool, str]:
    """
    Determine if 2 codons are synonymous.
    """
    if codon1 == codon2:
        return 'wt codon'  # changed from False
    if _translate_dna(codon1, codon_table) is not _translate_dna(codon2, codon_table):
        return False
    return True


def _translate_dna(seq: str, codon_table: Dict[str, str]) -> str:
    """
    Translate DNA sequence to protein.
    """
    protein: str = ''
    if len(seq) % 3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i : i + 3]
            protein += codon_table[codon]
    return protein


def enumerate_variants(wtseq_list: List[str], codon_list: List[str],
                       dna_sequence: str) -> Dict[str, int]:
    """
    Will return an ordered dictionary with variants initialized to 0 counts.
    """
    # Create ordered dictionary
    variants: dict = OrderedDict()

    # First instance that we see the wt?
    firstwtseq: bool = False

    # Loop over codons
    for position in range(0, len(wtseq_list)):
        for codons in codon_list:
            variant = ''.join(wtseq_list[0 : position]
                              ) + codons + ''.join(wtseq_list[position + 1 :])
            if variant == dna_sequence:  # Store redundant wild-types
                if firstwtseq:
                    variant = 'wtSeq' + str(position)
                firstwtseq = True
            variants[variant] = 0
    return variants


def initialize_ordered_dict(list_variants: List[str]) -> dict:
    """
    Will return an ordered dictionary with variants initialized to 0 counts.
    Here the user specifies the variants as a list.

    This function should be used when you want to use _count_fastq    """

    # Create normal dictionary
    dictionary = dict(zip(list_variants, np.zeros(len(list_variants))))

    # Create ordered dictionary
    return OrderedDict(dictionary)


def stopcodon_correction(
    input_lib: npt.NDArray, output_lib: npt.NDArray, input_stopcodon: npt.NDArray,
    output_stopcodon: npt.NDArray
) -> npt.NDArray:
    """
    This aux function will take as an input the counts for pre and post
    selection (and also for wT subset), and will return the corrected
    output counts.
    """

    # calculate stop codons frequencies
    frequency_stopcodons = output_stopcodon / input_stopcodon

    # MAD filtering
    frequency_stopcodons_filtered = filter_by_mad(frequency_stopcodons, m=2)
    median_frequency = np.nanmedian(frequency_stopcodons_filtered)

    # subtract to output counts
    output_lib_corr = output_lib - input_lib * median_frequency

    # eliminate negative values so they wont get turned into np.nan
    output_lib_corr[output_lib_corr < 0] = 0

    return output_lib_corr


def filter_by_mad(data: npt.NDArray, m: float = 2) -> npt.NDArray:
    """
    This aux function will take a numpy array, calculate median and MAD,
    and filter the data removing outliers.
    """

    # turn data into df_output to do mad calculations
    df_output = DataFrame(np.array(data), columns=['Data'])
    median = df_output['Data'].median(axis=0)
    mad = df_output['Data'].mad(axis=0)
    df_output['Abs_Dev'] = np.abs(data - median) / mad

    # filter values m times away from median, by default m = 2
    df_output['Abs_Dev'].mask(df_output['Abs_Dev'] > m, inplace=True)  # mask values
    df_output.dropna(how='any', inplace=True)  # eliminte NANs

    return df_output['Data'].to_numpy()


def replace_inf(array: npt.NDArray, infinite: float) -> npt.NDArray:
    """
    Replace values over a threshold with a min or max value.
    """
    #np.warnings.filterwarnings('ignore')
    array[array == -np.inf] = -infinite
    array[array < -infinite] = -infinite
    array[array == +np.inf] = +infinite
    array[array > +infinite] = +infinite
    return array


def group_by_aa(df_input: DataFrame, aminoacids: List[str]) -> DataFrame:
    """
    Group different codons that are synonymous.
    """
    # copy df
    df_output = df_input.copy()

    # Set up amino acid column
    df_output['Aminoacid'] = aminoacids

    # Group by mean
    df_output = df_output.groupby(as_index=True, by='Aminoacid', sort=False).mean()
    return df_output


def nan_mode(data: npt.NDArray) -> npt.NDArray:
    """
    Input is wt log enrichments, and return the mode of the histogram
    (aka the x coordinate at which y is max).
    """

    # Copy data
    data = np.copy(data)
    # Remove NaN values
    data_corrected = data[np.invert(np.isnan(data))]
    # Adjust kernel
    kernel_processed_data = gaussian_kde(data_corrected)
    # Find mode
    indexmax = np.where(
        kernel_processed_data(data_corrected) == kernel_processed_data(data_corrected).max()
    )
    # Return mean in case there are two x values with equal y-axis height
    return data_corrected[indexmax].mean()


# corrects the mutagenesis data and returns the height of the peak
def kernel_correction(data: Union[DataFrame, npt.NDArray],
                      cutoff: float = 2) -> Tuple[npt.NDArray, float]:
    """
    Input the library matrix, returns the corrected version. I set
    to 0 the max of the peak of the normal dist ignores stop codons.
    Not used for dataframes, only numpy arrays.
    """

    # Get data into right format
    data_corrected, kernel_processed_data = _kernel_data_preparation(data, cutoff)

    # Find max of kernel peak
    indexmax = np.where(
        kernel_processed_data(data_corrected) == kernel_processed_data(data_corrected).max()
    )

    # Normalize the max of peak os it has an x = 0
    data_final = data - data_corrected[indexmax].mean()

    # find std of kernel. It uses the already max peak x=0 normalized data
    data_final_flatten, data_final_kernel_processed_data = _kernel_data_preparation(
        data_final, cutoff
    )
    std = _kernel_std(data_final_flatten, data_final_kernel_processed_data)

    return data_final, std


def _kernel_data_preparation(data: DataFrame, cutoff: float) -> Tuple[npt.NDArray, npt.NDArray]:
    """
    This function will copy the data, eliminate stop codon, eliminate
    values lower than -1, flatten and eliminate np.nan. Will return the
    data in that format + the adjusted kernel.
    """

    # Eliminate stop codon
    data_corrected: npt.NDArray = np.array(data.drop('*', errors='ignore').copy())

    # Eliminate values lower than -1
    data_corrected = data_corrected[(data_corrected >= -cutoff) & (data_corrected <= cutoff)]

    # Get rid of np.nan values and convert matrix into 1d matrix
    data_corrected = data_corrected[np.invert(np.isnan(data_corrected))]

    # Adjust gaussian kernel
    kernel_processed_data = gaussian_kde(data_corrected)

    return data_corrected, kernel_processed_data


def _kernel_std(data: npt.NDArray, kernel: Any) -> float:
    """
    Input the library matrix (and wont count stop codon), and will return
    the std of the normal distribution.
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


def array_to_df_enrichments(lib: npt.NDArray, aminoacids: List[str]) -> DataFrame:
    """
    Aux function to transform array in df with index of amino acids.
    """

    df_output: DataFrame = DataFrame(index=aminoacids, data=lib)
    return df_output.astype(float)


def rearrange_dataframe(
    df_input: DataFrame, values_column: str, index_order: List[str]
) -> DataFrame:
    """
    Convert a df that contains mutagenesis data in column format into a
    matrix format.
    """
    df_pivoted: DataFrame = df_input.pivot_table(
        values=values_column, index='Aminoacid', columns=['Position'], dropna=False
    )
    return df_pivoted.reindex(index=index_order)
