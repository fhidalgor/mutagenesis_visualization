#!/usr/bin/env python
# coding: utf-8

# # Import Modules

# In[1]:


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


# # Data Process Functions

# ## Process trimmed fastq file

# In[7]:


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
    assert len(
        dna_sequence
    ) % 3 == 0, 'The dna_sequence length is not a multiple of 3'

    # Make upper case in case input was lower case
    dna_sequence = dna_sequence.upper()
    if isinstance(codon_list, str):  #check if codon_list is a String
        codon_list = codon_list.upper()
    else:
        codon_list = [item.upper() for item in codon_list]

    # Create list with codons of sequence
    wtSeqList = [dna_sequence[i:i + 3] for i in range(0, len(dna_sequence), 3)]

    # codon_list
    if codon_list == 'NNS':
        codon_list = [
            "GCC", "GCG", "TGC", "GAC", "GAG", "TTC", "GGC", "GGG", "CAC",
            "ATC", "AAG", "CTC", "CTG", "TTG", "ATG", "AAC", "CCC", "CCG",
            "CAG", "CGC", "CGG", "AGG", "TCC", "TCG", "AGC", "ACC", "ACG",
            "GTC", "GTG", "TGG", "TAC", "TAG"
        ]
    elif codon_list == 'NNK':
        codon_list = [
            'GCG', 'GCT', 'TGT', 'GAT', 'GAG', 'TTT', 'GGG', 'GGT', 'CAT',
            'ATT', 'AAG', 'CTG', 'CTT', 'TTG', 'ATG', 'AAT', 'CCG', 'CCT',
            'CAG', 'AGG', 'CGG', 'CGT', 'AGT', 'TCG', 'TCT', 'ACG', 'ACT',
            'GTG', 'GTT', 'TGG', 'TAT', 'TAG'
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
    df['SynWT'] = df.apply(
        lambda x: _are_syn(x['Codon'], x['WTCodon'], _codon_table()), axis=1
    )
    df['Counts'] = list(variants.values())

    if counts_wt:
        try:  # try is to fix the Bug Che discovered
            df.loc[df['Codon'] == df['WTCodon'],
                   'Counts'] = variants[dna_sequence]
        except:
            pass
    else:
        df.loc[df['Codon'] == df['WTCodon'], 'Counts'] = np.nan

    # Pivot table and reindex
    df_counts = df.pivot_table(
        values='Counts', index='Codon', columns=['Position'], dropna=False
    )
    df_counts = df_counts.reindex(index=codon_list)

    # Get WT counts
    df_wt = df.loc[df['SynWT'] == True][[
        'Position', 'Codon', 'Aminoacid', 'Counts'
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


def _codon_table():
    codontable = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M', 'ACA': 'T', 'ACC': 'T',
        'ACG': 'T', 'ACT': 'T', 'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R', 'CTA': 'L', 'CTC': 'L',
        'CTG': 'L', 'CTT': 'L', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGA': 'R', 'CGC': 'R',
        'CGG': 'R', 'CGT': 'R', 'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A', 'GAC': 'D', 'GAT': 'D',
        'GAA': 'E', 'GAG': 'E', 'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S', 'TTC': 'F', 'TTT': 'F',
        'TTA': 'L', 'TTG': 'L', 'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
        'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W'
    }
    return codontable


def _are_syn(codon1, codon2, codontable):
    '''Determine if 2 codons are synonymous'''
    if codon1 == codon2:
        return False
    if _translate(codon1, codontable) is not _translate(codon2, codontable):
        return False
    return True


def _translate(seq, codontable):
    '''Translate DNA sequence to protein.'''
    # I forgot why I made this custom function instead of using a biopython function
    protein = ''
    if len(seq) % 3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            protein += codontable[codon]
    return protein


def _enumerate_variants(wtSeqList, codon_list, dna_sequence):
    '''Will return an ordered dictionary with variants initialized to 0 counts'''
    # Create ordered dictionary
    variants = OrderedDict()

    # First instance that we see the wt?
    firstwtseq = False

    # Loop over codons
    for position in range(0, len(wtSeqList)):
        for codons in (codon_list):
            variant = ''.join(wtSeqList[0:position]) +                 codons + ''.join(wtSeqList[position+1:])
            if (variant == dna_sequence):  # Store redundant wild-types
                if firstwtseq:
                    variant = 'wtSeq' + str(position)
                firstwtseq = True
            variants[variant] = 0
    return variants


def count_fastq(variants, input_file):
    '''
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
    '''
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


def _initialize_ordereddict(list_variants):
    '''
    Will return an ordered dictionary with variants initialized to 0 counts.
    Here the user specifies the variants as a list.

    This function should be used when you want to use _count_fastq
    
    '''

    # Create normal dictionary
    dictionary = dict(zip(list_variants, np.zeros(len(list_variants))))

    # Create ordered dictionary
    variants = OrderedDict(dictionary)

    return variants


# ## Process count files and return enrichment

# In[ ]:


def calculate_enrichment(
    pre_lib,
    post_lib,
    pre_wt=None,
    post_wt=None,
    aminoacids=list('AACDEFGGHIKLLLMNPPQRRRSSSTTVVWY*'),
    zeroing='population',
    how='median',
    norm_std=True,
    stopcodon=False,
    min_counts=25,
    min_countswt=100,
    std_scale=0.2,
    mpop=2,
    mwt=2,
    infinite=3,
    output_file: Union[None, str, Path] = None
):
    """
    Determine the enrichment scores of a selection experiment, where there is a
    preselected population (input) and a selected population (output).

    Parameters
    -----------
    pre_lib : str, pandas dataframe or np.array
        Can be filepath and name of the exported txt file, dataframe or np.array.

    post_lib : str, pandas dataframe or np.array
        Can be filepath and name of the exported txt file, dataframe or np.array.

    pre_wt : str, or np.array, optional
        Str with filepath and name of the exported txt file or np.array.

    post_wt : str, or np.array, optional
        Str with filepath and name of the exported txt file or np.array.

    aminoacids : list, default ('AACDEFGGHIKLLLMNPPQRRRSSSTTVVWY*')
        Index of aminoacids (in order). Stop codon needs to be '*'.

    zeroing : str, default 'population'
        Method to normalize the data.
        Can also use 'zscore, 'counts', wt' or 'kernel'.

    how : str, default 'median'
        Metric to zero the data. Only works if zeroing='population' or 'wt'.
        Can also be set to 'mean' or 'mode'.

    norm_std : boolean, default True
        If norm_std is set to True, it will scale the data.

    stopcodon : boolean, default False
        Use the enrichment score stop codons as a metric to determine the minimum enrichment score.

    min_counts : int, default 25
        If mutant has less than the min_counts, it will be replaced by np.nan.

    min_countswt : int, default 100
        If synonymous wild-type mutant has less than the min_counts, it will be replaced by np.nan.

    std_scale : float, default 0.2
        Factor by which the population is scaled. Only works if norm_std is set to True.

    mpop : int, default 2
        When using the median absolute deviation (MAD) filtering, mpop is the number of medians away
        a data point must be to be discarded.

    mwt : int, default 2
        When MAD filtering, mpop is the number of medians away a data point must be to
        be discarded. The difference with mpop is that mwt is only used when the population of wild-type
        alleles is the reference for data zeroing.

    infinite : int, default 3
        It will replace +infinite values with +3 and -infinite with -3.

    output_file : str, default None
        If you want to export the generated files, add the path and name of the file without suffix.
        Example: 'path/filename'. File will be save as a txt file.

    Returns
    --------
    zeroed : ndarray
        A np.array containing the enrichment scores.

    """

    # Convert to numpy if libraries are in dataframe format.
    # If input is a filepath, then load the txt files
    if type(pre_lib) is pd.DataFrame:
        pre_lib = pre_lib.to_numpy()
    elif type(pre_lib) is str:
        pre_lib = np.loadtxt(pre_lib)
    if type(post_lib) is pd.DataFrame:
        post_lib = post_lib.to_numpy()
    elif type(post_lib) is str:
        post_lib = np.loadtxt(post_lib)

    # Same thing for wt allele files
    if type(pre_wt) is str:
        pre_wt = np.loadtxt(pre_wt)
    if type(post_wt) is str:
        post_wt = np.loadtxt(post_wt)

    # Convert to df
    pre_lib = _array_to_df_enrichments(pre_lib, aminoacids)
    post_lib = _array_to_df_enrichments(post_lib, aminoacids)

    # Locate stop codons
    if stopcodon:
        input_stopcodon = pre_lib.loc['*'].astype(float)
        output_stopcodon = post_lib.loc['*'].astype(float)
    else:
        input_stopcodon = ''
        output_stopcodon = ''

    # Log10 of the counts for library and wt alleles
    log10_counts = _get_enrichment(
        pre_lib, post_lib, input_stopcodon, output_stopcodon, min_counts,
        stopcodon, infinite
    )
    # Group by amino acid
    df = pd.DataFrame(data=log10_counts)
    log10_counts_grouped = _group_byaa(df, aminoacids)

    # MAD filtering
    log10_counts_mad = _MAD_filtering(
        np.ravel(np.array(log10_counts_grouped)), mpop
    )
    mean_pop = np.nanmean(log10_counts_mad)
    median_pop = np.nanmedian(log10_counts_mad)
    std_pop = np.nanstd(log10_counts_mad)
    mode_pop = _nanmode(log10_counts_mad)

    # Wt counts
    if pre_wt is not None:
        log10_wtcounts = _get_enrichment(
            pre_wt, post_wt, input_stopcodon, output_stopcodon, min_countswt,
            stopcodon, infinite
        )
        # MAD filtering
        # If set to m=1, if tosses out about 50% of the values. the mean barely changes though
        log10_wtcounts = _MAD_filtering(log10_wtcounts, mwt)
        mean_wt = np.nanmean(log10_wtcounts)
        median_wt = np.nanmedian(log10_wtcounts)
        std_wt = np.nanstd(log10_wtcounts)
        mode_wt = _nanmode(log10_wtcounts)

    # Zero data, select case
    if zeroing == 'wt':
        if how == 'mean':
            zeroed = log10_counts_grouped - mean_wt
        elif how == 'median':
            zeroed = log10_counts_grouped - median_wt
        elif how == 'mode':
            zeroed = log10_counts_grouped - mode_wt
        elif norm_std == True:
            zeroed = zeroed * std_scale / 2 / std_wt
    elif zeroing == 'population':
        if how == 'mean':
            zeroed = log10_counts_grouped - mean_pop
        elif how == 'median':
            zeroed = log10_counts_grouped - median_pop
        elif how == 'mode':
            zeroed = log10_counts_grouped - mode_pop
        elif norm_std == True:
            zeroed = zeroed * std_scale / std_pop
    elif zeroing == 'counts':
        # Get the ratio of counts
        ratio = np.log10(pre_lib.sum().sum() / post_lib.sum().sum())
        zeroed = log10_counts_grouped + ratio
        if norm_std == True:
            zeroed = zeroed * std_scale / std_pop
    elif zeroing == 'kernel':
        zeroed_0, kernel_std = _kernel_correction(
            log10_counts_grouped, aminoacids
        )
        zeroed, kernel_std = _kernel_correction(zeroed_0, aminoacids, cutoff=1)
        if norm_std is True:
            zeroed = zeroed * std_scale / kernel_std
    elif zeroing == 'zscore':
        zeroed = stats.zscore(zeroed, nan_policy='propagate')

    # Export files
    if output_file:
        np.savetxt(Path(output_file), zeroed, fmt='%i', delimiter='\t')

    return zeroed


# 
# ## Aux functions

# In[ ]:


def _get_enrichment(
    input_lib, output_lib, input_stopcodon, output_stopcodon, min_counts,
    stopcodon, infinite
):
    '''Calculate log10 enrichment scores from input and output counts'''
    # Copy data and replace low counts by np.nan
    input_lib = np.copy(input_lib.astype(float))
    output_lib = np.copy(output_lib.astype(float))
    input_lib[input_lib < min_counts] = np.nan

    # Stop codon correction
    if stopcodon:
        output_lib = _stopcodon_correction(
            input_lib, output_lib, input_stopcodon, output_stopcodon
        )

    # log10 of library and replace infinite values
    counts_log10_ratio = _replace_inf(
        np.log10(output_lib / input_lib), infinite
    )

    return counts_log10_ratio


def _stopcodon_correction(
    input_lib, output_lib, input_stopcodon, output_stopcodon
):
    '''This aux function will take as an input the counts for pre and post selection (and also for wT subset), 
    and will return the corrected output counts'''

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
    '''This aux function will take a numpy array, calculate median and MAD, 
    and filter the data removing outliers'''

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
    '''Replace values over a threshold with a min or max value'''
    np.warnings.filterwarnings('ignore')
    array[array == -np.inf] = -infinite
    array[array < -infinite] = -infinite
    array[array == +np.inf] = +infinite
    array[array > +infinite] = +infinite
    return array


def _group_byaa(df, aminoacids):
    '''Group different codons that are synonymous'''
    # copy df
    df = df.copy()

    # Set up amino acid column
    df['Aminoacid'] = aminoacids

    # Group by mean
    df = df.groupby(as_index=True, by='Aminoacid', sort=False).mean()
    return df


def _nanmode(data):
    '''
    Input is wt log enrichments, and return the mode of the histogram 
    (aka the x coordinate at which y is max).
    
    '''

    # Copy data
    data = np.copy(data)
    # Remove NaN values
    data_corrected = data[np.invert(np.isnan(data))]
    # Adjust kernel
    kernel_processed_data = stats.gaussian_kde(data_corrected)
    # Find mode
    indexmax = np.where(
        kernel_processed_data(data_corrected) ==
        kernel_processed_data(data_corrected).max()
    )
    # Return mean in case there are two x values with equal y-axis height
    return data_corrected[indexmax].mean()


# corrects the mutagenesis data and returns the height of the peak
def _kernel_correction(data, aminoacids, cutoff=2):
    '''input the library matrix, returns the corrected version. I set to 0 the max of the peak of the normal dist
    ignores stop codons. Not used for dataframes, only numpy arrays'''

    # Get data into right format
    data_corrected, kernel_processed_data = _kernel_datapreparation(
        data, cutoff
    )

    # Find max of kernel peak
    indexmax = np.where(
        kernel_processed_data(data_corrected) ==
        kernel_processed_data(data_corrected).max()
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
    '''
    This function will copy the data, eliminate stop codon, eliminate values lower than -1, 
    flatten and eliminate np.nan. Will return the data in that format + the adjusted kernel
    
    '''
    
    # Eliminate stop codon
    data_corrected = np.array(data.drop('*', errors='ignore').copy())

    # Eliminate values lower than -1
    data_corrected = data_corrected[(data_corrected >= -cutoff)
                                    & (data_corrected <= cutoff)]

    # Get rid of np.nan values and convert matrix into 1d matrix
    data_corrected = data_corrected[np.invert(np.isnan(data_corrected))]

    # Adjust gaussian kernel
    kernel_processed_data = stats.gaussian_kde(data_corrected)

    return data_corrected, kernel_processed_data


def _kernel_std(data, kernel):
    '''
    Input the library matrix (and wont count stop codon), and will return the std of the normal distribution.
    To calculate the std, it will find the FWHM and divide by 2.355
    https://en.wikipedia.org/wiki/Full_width_at_half_maximum
    The algorithm will give back the min std between both sides of the peak.
    
    '''

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
    '''
    aux function to transform array in df with index of amino acids.
    
    '''
    
    df = pd.DataFrame(index=aminoacids, data=lib)
    return df.astype(float)


# ## Assemble sublibraries

# In[ ]:


def assemble_avengers(
    excel_path,
    sheet_pre,
    sheet_post,
    columns,
    nrows_pop,
    nrows_wt,
    columns_wt=None,
    skiprows=1,
    aminoacids=list('AACDEFGGHIKLLLMNPPQRRRSSSTTVVWY*'),
    zeroing='population',
    how='median',
    norm_std=True,
    stopcodon=False,
    min_counts=25,
    min_countswt=100,
    std_scale=0.2,
    mpop=2,
    mwt=2,
    infinite=3,
    output_file: Union[None, str, Path] = None
):
    """
    Assembles different sublibraries into one. Uses calculate_enrichments. Can
    only read from excel files that are in the same format as the example
    provided.

    Parameters
    -----------
    excel_path : str
        Location of the excel file to read.

    sheet_pre : str
        Name of the sheet with input (pre-selected) counts.

    sheet_post : str
        Name of the sheet with output (post-selected) counts.

    columns : list
        List of columns for each sublibrary to read from the excel file.

    nrows_pop : int,
        Number of rows to read from the excel.

    nrows_wt : list,
        Contains a list of integers, with the number of rows to read from each wt subset.

    columns_wt : list,
        Contains a list of strings, specifying the excel columns to read for each wt subset.

    skiprows : int, default 1
        Parameter for pd.read_excel. Only works for the main columns, not for wt.

    aminoacids : list, default ('AACDEFGGHIKLLLMNPPQRRRSSSTTVVWY*')
        Index of aminoacids (in order). Stop codon needs to be '*'.

    zeroing : str, default 'population'
        Method to zero the data.
        Can also use 'counts', wt' or 'kernel'.

    how : str, default 'median'
        Metric to zero the data. Only works if zeroing='population' or 'wt'.
        Can also be set to 'mean' or 'mode'.

    norm_std : boolean, default True
        If norm_std is set to True, it will scale the data.

    stopcodon : boolean, default False
        Use the enrichment score stop codons as a metric to determine the minimum enrichment score.

    min_counts : int, default 25
        If mutant has less than the min_counts, it will be replaced by np.nan.

    min_countswt : int, default 100
        If synonymous wild-type mutant has less than the min_counts, it will be replaced by np.nan.

    std_scale : float, default 0.2
        Factor by which the population is scaled. Only works if norm_std is set to True.

    mpop : int, default 2
        When using the median absolute deviation (MAD) filtering, mpop is the number of medians away
        a data point must be to be discarded.

    mwt : int, default 2
        When MAD filtering, mpop is the number of medians away a data point must be to
        be discarded. The difference with mpop is that mwt is only used when the population of wild-type
        alleles is the reference for data zeroing.

    infinite : int, default 3
        It will replace +infinite values with +3 and -infinite with -3.

    output_file : str, default None
        If you want to export the generated files, add the path and name of the file without suffix.
        Example: 'path/filename'. File will be save as a txt file.

    Returns
    --------
    df : Pandas dataframe
        A dataframe that contains the enrichment scores of the assembled sublibraries.

    """

    # Read reads from excel
    list_pre, list_sel, list_pre_wt, list_sel_wt = _read_counts(
        excel_path, sheet_pre, sheet_post, columns, nrows_pop, nrows_wt,
        columns_wt
    )

    # Assemble sublibraries
    df = _assemble_list(
        list_pre, list_sel, list_pre_wt, list_sel_wt, aminoacids, zeroing, how,
        norm_std, stopcodon, min_counts, min_countswt, std_scale, mpop, mwt,
        infinite
    )

    # Export files
    if output_file:
        np.savetxt(Path(output_file), zeroed, fmt='%i', delimiter='\t')

    return df


def _read_counts(
    excel_path,
    sheet_pre,
    sheet_post,
    columns,
    nrows_pop,
    nrows_wt,
    columns_wt=None,
    skiprows=1
):
    '''Aux'''
    # Create dictionary with data. Loading 3 replicates, each of them is divided into 3 pools
    list_pre, list_sel, list_pre_wt, list_sel_wt = ([] for i in range(4))

    # Read counts from excel
    replicates = np.arange(0, len(sheet_pre))
    for column, column_wt, nrow_wt, rep in zip(columns, columns_wt, nrows_wt,
                                               replicates):
        # Pre counts
        list_pre.append(
            pd.read_excel(
                excel_path,
                sheet_pre,
                skiprows=skiprows,
                usecols=column,
                nrows=nrows_pop
            )
        )
        # Sel counts
        list_sel.append(
            pd.read_excel(
                excel_path,
                sheet_post,
                skiprows=skiprows,
                usecols=column,
                nrows=nrows_pop
            )
        )
        if columns_wt is None:
            list_pre_wt.append(None)
            list_sel_wt.append(None)
        else:
            # Pre counts wild-type alleles
            list_pre_wt.append(
                pd.read_excel(
                    excel_path, sheet_pre, usecols=column_wt, nrows=nrow_wt
                )
            )
            # Sel counts wild-type alleles
            list_sel_wt.append(
                pd.read_excel(
                    excel_path, sheet_post, usecols=column_wt, nrows=nrow_wt
                )
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
    '''
    gets the output from _read_counts and assembles the sublibraries
    
    '''

    enrichment_lib = []

    for pre, sel, pre_wt, sel_wt in zip(list_pre, list_sel, list_pre_wt,
                                        list_sel_wt):
        # log 10
        enrichment_log10 = calculate_enrichment(
            pre, sel, pre_wt, sel_wt, aminoacids, zeroing, how, norm_std,
            stopcodon, min_counts, min_countswt, std_scale, mpop, mwt, infinite
        )
        # Store in list
        enrichment_lib.append(enrichment_log10)

    # Concatenate sublibraries
    df = pd.concat(enrichment_lib, ignore_index=True, axis=1)
    return df


# ## Merge enrichment with MSA conservation

# In[ ]:


def msa_enrichment(self, path, start_position, threshold=0.01):
    '''
    Generate a dataframe with the Shannon entropy by residue and the mean enrichment score, and
    a second dataframe with the frequency of each substitution and the enrichment score

    Parameters
    -----------
    self : object from class *Screen*

    path : str
        Path where is located the fasta MSA that will be parsed. That MSA needs to have removed 
        any insertions that are not present in the target sequence. For example, if a Ras ortholog has 
        an extra amino acid at position 123, that needs to be removed from the aligment. Otherwise, everything
        will be shifted by 1 residue.

    start_position : int
        This is the position in the protein sequence of the first position in the MSA.

    threshold : float, default 0.01
        The conservation frequency for each amino acid subsitution will be binarized, and a threshold between 0-1 
        needs to be selected.

    Returns
    --------
    df_shannon: pandas dataframe
        Shannon entropy by residue and mean enrichment score by residue. 

    df_freq : pandas dataframe   
        Frequency of each susbsitution merged to the enrichment score.
    '''
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
    '''Convert a msa from a fasta file into a df ready to plot with logomaker. Returns frequency'''
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

