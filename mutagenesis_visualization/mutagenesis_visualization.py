#!/usr/bin/env python
# coding: utf-8

# # Import Modules

# In[ ]:


from __future__ import unicode_literals
import numpy as np
import seaborn as sns
import pandas as pd
import itertools
import copy
from scipy import stats
from matplotlib import rcParams
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.ticker as ticker
from sklearn import metrics
from sklearn.decomposition import PCA
from adjustText import adjust_text
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.PDB import PDBParser
from collections import defaultdict, OrderedDict, Counter
import freesasa
from os import path
import os
import sys
import logomaker
import adjustText

try:
    import shannon
except ModuleNotFoundError:
    dummy = 0
    
try:
    from ipymol import viewer as pymol
except ModuleNotFoundError:
    dummy = 0


# # Data Process Functions

# ## Process trimmed fastq file

# In[ ]:


def count_reads(dna_sequence, codon_list='NNS', **kwargs):
    '''
    Process a trimmed fastq file containing DNA reads and returns the counts of 
    each DNA sequence specified by the user.

    Parameters
    -----------
    dna_sequence : str, 
        Contains the DNA sequence of the allele of reference (usually wild-type).
    
    codon_list : list or str, default 'NNS'
        Input a list of the codons that were used to create point mutations. Example: ["GCC", "GCG", "TGC"].
        If the library was built using NNS and NNK codons, it is enough to input 'NNS' or 'NNK' as a string. 
        It is important to know that the order of the codon_list will determine the output order.
    
    inputfilepath : str, not optional kwarg
        Path where the input fastq file is stored.
    
    inputfilename : str, not optional kwarg
        Name of the fastq file (full name including ".fastq").
    
    outputfilepath : str, optional kwarg
        Path where the output files will be saved.
    
    outputfilename : str, optional kwarg
        Name of the output files.
    
    savefile : boolean, default False. optional kwarg 
        If set to true, the function will export the two arrays to separate txt files using the
        output filepath and output filename specified by the user.

    Returns
    --------
    df_counts : dataframe 
        Dataframe with the counts for each point mutant.
    
    wt_counts : list
        List of the counts for each for each DNA sequence that codes for the wild-type protein.
    '''
    # update kwargs
    temp_kwargs = copy.deepcopy(default_kwargs)
    temp_kwargs.update(kwargs)

    # Naming files
    trimmedfile = temp_kwargs['inputfilepath']+temp_kwargs['inputfilename']

    # Make upper case in case input was lower case
    dna_sequence = dna_sequence.upper()
    codon_list = [item.upper() for item in codon_list]
    
    # Create list with codons of sequence
    wtSeqList = [dna_sequence[i:i+3] for i in range(0, len(dna_sequence), 3)]

    # codon_list
    if codon_list == 'NNS':
        codon_list = ["GCC", "GCG", "TGC", "GAC", "GAG", "TTC", "GGC", "GGG", "CAC", "ATC", "AAG", "CTC", "CTG", "TTG", "ATG",
                      "AAC", "CCC", "CCG", "CAG", "CGC", "CGG", "AGG", "TCC", "TCG", "AGC", "ACC", "ACG", "GTC", "GTG", "TGG", "TAC", "TAG"]
    elif codon_list == 'NNK':
        codon_list = ['GCG', 'GCT', 'TGT', 'GAT', 'GAG', 'TTT', 'GGG', 'GGT', 'CAT', 'ATT', 'AAG', 'CTG', 'CTT', 'TTG', 'ATG',
                      'AAT', 'CCG', 'CCT', 'CAG', 'AGG', 'CGG', 'CGT', 'AGT', 'TCG', 'TCT', 'ACG', 'ACT', 'GTG', 'GTT', 'TGG', 'TAT', 'TAG']

    # Enumerate variants
    variants = OrderedDict()
    isitthefirstwtseq = False
    for position, codon in enumerate(wtSeqList):
        for position2, codons in enumerate(codon_list):
            variant = ''.join(wtSeqList[0:position]) +                 ''.join(codons) + ''.join(wtSeqList[position+1:])
            if (variant == dna_sequence):
                if isitthefirstwtseq:
                    variant = 'wtSeq' + str(position)
                isitthefirstwtseq = True
            variants[variant] = 0

    # Translate nucleotide sequence and count variant frequency
    totalreads = 0
    for nuc in SeqIO.parse(trimmedfile, "fastq"):
        totalreads += 1
        nucleicsequence = str(nuc.seq)
        if nucleicsequence in variants:
            variants[nucleicsequence] += 1
    usefulreads = np.nansum(list(variants.values()))

    # Convert to df
    wtProtein = Seq(dna_sequence).translate()
    df = pd.DataFrame()
    df['Position'] = np.ravel([[pos]*len(codon_list)
                               for pos in np.arange(1, len(wtProtein)+1).astype(int)])
    df['Codon'] = codon_list*len(wtProtein)
    df['WTCodon'] = np.ravel([[codon]*len(codon_list) for codon in wtSeqList])
    df['Aminoacid'] = np.ravel([[aa]*len(codon_list) for aa in wtProtein])
    codontable = _codon_table()
    df['SynWT'] = df.apply(lambda x: _are_syn(
        x['Codon'], x['WTCodon'], codontable), axis=1)
    df['Counts'] = list(variants.values())
    df.loc[df['Codon'] == df['WTCodon'], 'Counts'] = variants[dna_sequence]

    # Pivot table and reindex
    df_counts = df.pivot_table(values='Counts', index='Codon',
                                     columns=['Position'], dropna=False)
    df_counts = df_counts.reindex(index=codon_list)

    # Get WT counts
    df_wt = df.loc[df['SynWT'] == True]
    wt_counts = list(df_wt['Counts'])
    wt_counts.insert(0, int([variants[dna_sequence]][0]))

    # Export files
    if temp_kwargs['savefile']:
        # Generate file handles
        outputfile_counts = temp_kwargs['outputfilepath'] +             temp_kwargs['outputfilename']+"_counts.txt"
        outputfile_wtcounts = temp_kwargs['outputfilepath'] +             temp_kwargs['outputfilename']+"_wtcounts.txt"
        # Save to txt files
        np.savetxt(outputfile_counts, np.array(df_counts), fmt='%i', delimiter='\t')
        np.savetxt(outputfile_wtcounts, wt_counts, fmt='%i', delimiter='\t')

    # Print total reads
    print('{}/{} useful reads ({}%)'.format(str(usefulreads),
                                            str(totalreads), str(int(usefulreads/totalreads*100))))
    return df_counts, wt_counts


def _codon_table():
    codontable = {'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
                  'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
                  'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
                  'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
                  'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
                  'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
                  'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
                  'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
                  'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
                  'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
                  'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
                  'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
                  'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
                  'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
                  'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
                  'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W'}
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


# ## Process count files and return enrichment

# In[ ]:


def calculate_enrichment(pre_lib, post_lib, pre_wt=None, post_wt=None, aminoacids=list('AACDEFGGHIKLLLMNPPQRRRSSSTTVVWY*'),
                         zeroing='population', how='median', norm_std=True, stopcodon=False, min_counts=25,
                         min_countswt=100, std_scale = 0.2, mpop=2, mwt=2, infinite=3, **kwargs):
    '''
    Determine the enrichment scores of a selection experiment, where there is a preselected population (input)
    and a selected population (output).

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
    
    **kwargs : other keyword arguments
        savefile : boolean, default False. optional kwarg 
            If set to true, the function will export the enrichment scores to a txt file.
    
    Returns
    --------
    zeroed : ndarray
        A np.array containing the enrichment scores. 
    '''
    # update kwargs
    temp_kwargs = copy.deepcopy(default_kwargs)
    temp_kwargs.update(kwargs)

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
    log10_counts = _get_enrichment(pre_lib, post_lib, input_stopcodon,
                                   output_stopcodon, min_counts, stopcodon, infinite)
    # Group by amino acid
    df = pd.DataFrame(data=log10_counts)
    log10_counts_grouped = _group_byaa(df, aminoacids)

    # MAD filtering
    log10_counts_mad = _MAD_filtering(np.ravel(np.array(log10_counts_grouped)), mpop)
    mean_pop = np.nanmean(log10_counts_mad)
    median_pop = np.nanmedian(log10_counts_mad)
    std_pop = np.nanstd(log10_counts_mad)
    mode_pop = _nanmode(log10_counts_mad)

    # Wt counts
    if pre_wt is not None:
        log10_wtcounts = _get_enrichment(pre_wt, post_wt, input_stopcodon,
                                         output_stopcodon, min_countswt, stopcodon, infinite)
        # MAD filtering
        # If set to m=1, if tosses out about 50% of the values. the mean barely changes though
        log10_wtcounts = _MAD_filtering(log10_wtcounts, mwt)
        mean_wt = np.nanmean(log10_wtcounts)
        median_wt = np.nanmedian(log10_wtcounts)
        std_wt = np.nanstd(log10_wtcounts)
        mode_wt = _nanmode(log10_wtcounts)

    # Zero data
    if zeroing == 'wt':
        if how == 'mean':
            zeroed = log10_counts_grouped - mean_wt
        elif how == 'median':
            zeroed = log10_counts_grouped - median_wt
        elif how == 'mode':
            zeroed = log10_counts_grouped - mode_wt
        if norm_std == True:
            zeroed = zeroed*std_scale/2/std_wt
    elif zeroing == 'population':
        if how == 'mean':
            zeroed = log10_counts_grouped - mean_pop
        elif how == 'median':
            zeroed = log10_counts_grouped - median_pop
        elif how == 'mode':
            zeroed = log10_counts_grouped - mode_pop
        if norm_std == True:
            zeroed = zeroed*std_scale/std_pop
    elif zeroing == 'counts':
        # Get the ratio of counts
        ratio = np.log10(pre_lib.sum().sum()/post_lib.sum().sum())
        zeroed = log10_counts_grouped + ratio
        if norm_std == True:
            zeroed = zeroed*std_scale/std_pop
    elif zeroing == 'kernel':
        zeroed_0, kernel_std = _kernel_correction(
            log10_counts_grouped, aminoacids)
        zeroed, kernel_std = _kernel_correction(zeroed_0, aminoacids, cutoff=1)
        if norm_std is True:
            zeroed = zeroed*std_scale/kernel_std

    # Export files
    if temp_kwargs['savefile']:
        np.savetxt(temp_kwargs['outputfilepath']+temp_kwargs['outputfilename'],
                   zeroed, fmt='%i', delimiter='\t')

    return zeroed


# 
# ## Aux functions

# In[ ]:


def _get_enrichment(input_lib, output_lib, input_stopcodon, output_stopcodon,
                    min_counts, stopcodon, infinite):
    '''Calculate log10 enrichment scores from input and output counts'''
    # Copy data and replace low counts by np.nan
    input_lib = np.copy(input_lib.astype(float))
    output_lib = np.copy(output_lib.astype(float))
    input_lib[input_lib < min_counts] = np.nan

    # Stop codon correction
    if stopcodon:
        output_lib = _stopcodon_correction(
            input_lib, output_lib, input_stopcodon, output_stopcodon)

    # log10 of library and replace infinite values
    counts_log10_ratio = _replace_inf(np.log10(output_lib/input_lib), infinite)

    return counts_log10_ratio


def _stopcodon_correction(input_lib, output_lib, input_stopcodon, output_stopcodon):
    '''This aux function will take as an input the counts for pre and post selection (and also for wT subset), 
    and will return the corrected output counts'''

    # calculate stop codons frequencies
    frequency_stopcodons = output_stopcodon/input_stopcodon

    # MAD filtering
    frequency_stopcodons_filtered = _MAD_filtering(frequency_stopcodons, m=2)
    median_frequency = np.nanmedian(frequency_stopcodons_filtered)

    # subtract to output counts
    output_lib_corr = output_lib - input_lib*median_frequency

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
    df['Abs_Dev'] = np.abs(data-median)/mad

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
    df = df.groupby(as_index=True, by='Aminoacid', sort= False).mean()
    return df


def _nanmode(data):
    '''input is wt log enrichments, and return the mode of the histogram 
    (aka the x coordinate at which y is max)'''
    
    # Copy data
    data = np.copy(data)
    # Remove NaN values
    data_corrected = data[np.invert(np.isnan(data))]
    # Adjust kernel
    kernel_processed_data = stats.gaussian_kde(data_corrected)
    # Find mode
    indexmax = np.where(kernel_processed_data(data_corrected)
                        == kernel_processed_data(data_corrected).max())
    # Return mean in case there are two x values with equal y-axis height
    return data_corrected[indexmax].mean()


# corrects the mutagenesis data and returns the height of the peak
def _kernel_correction(data, aminoacids, cutoff=2):
    '''input the library matrix, returns the corrected version. I set to 0 the max of the peak of the normal dist
    ignores stop codons. Not used for dataframes, only numpy arrays'''

    # Get data into right format
    data_corrected, kernel_processed_data = _kernel_datapreparation(data, cutoff) 
                                                                    
    # Find max of kernel peak
    indexmax = np.where(kernel_processed_data(data_corrected) ==
                        kernel_processed_data(data_corrected).max())

    # Normalize the max of peak os it has an x = 0
    data_final = data-data_corrected[indexmax].mean()

    # find std of kernel. It uses the already max peak x=0 normalized data
    data_final_flatten, data_final_kernel_processed_data = _kernel_datapreparation(
        data_final, cutoff)
    std = _kernel_std(data_final_flatten, data_final_kernel_processed_data)

    return data_final, std


def _kernel_datapreparation(data, cutoff):
    '''this function will copy the data, eliminate stop codon, eliminate values lower than -1, 
    flatten and eliminate np.nan. Will return the data in that format + the adjusted kernel'''
    # Eliminate stop codon
    data_corrected = np.array(data.drop('*', errors = 'ignore').copy())
    
    # Eliminate values lower than -1
    data_corrected = data_corrected[(
        data_corrected >= -cutoff) & (data_corrected <= cutoff)]

    # Get rid of np.nan values and convert matrix into 1d matrix
    data_corrected = data_corrected[np.invert(np.isnan(data_corrected))]

    # Adjust gaussian kernel
    kernel_processed_data = stats.gaussian_kde(data_corrected)

    return data_corrected, kernel_processed_data


def _kernel_std(data, kernel):
    '''input the library matrix (and wont count stop codon), and will return the std of the normal distribution.
    To calculate the std, it will find the FWHM and divide by 2.355
    https://en.wikipedia.org/wiki/Full_width_at_half_maximum
    The algorithm will give back the min std between both sides of the peak'''

    # find ymax and the x value of the max height
    y_max = kernel(data).max()
    index_y_max = np.where(kernel(data) == y_max)
    x_ymax = data[index_y_max].mean()

    # find the two x value of ymax/2. One at each side of the center. l of left and r of right
    # so I can select only the positive side of the distribution
    y_temp = kernel(data)

    # left side
    y_hw_l = (min(y_temp[data < 0], key=lambda x: abs(x-y_max/2)))
    index_yhw_l = np.where(kernel(data) == y_hw_l)
    x_yhw_l = data[index_yhw_l].mean()

    # right side
    y_hw_r = (min(y_temp[data > 0], key=lambda x: abs(x-y_max/2)))
    index_yhw_r = np.where(kernel(data) == y_hw_r)
    x_yhw_r = data[index_yhw_r].mean()

    # calculate half width at half maximum
    hwhm_l = abs(x_yhw_l - x_ymax)
    hwhm_r = abs(x_yhw_r - x_ymax)

    # calculate std from fwhm
    std_l = hwhm_l/((2*np.log(2))**0.5)
    std_r = hwhm_r/((2*np.log(2))**0.5)

    return min(std_l, std_r)


def _array_to_df_enrichments(lib, aminoacids):
    '''aux function to transform array in df with index of amino acids'''
    df = pd.DataFrame(index=aminoacids, data=lib)
    return df.astype(float)
                 


# ## Assemble sublibraries

# In[ ]:


def assemble_avengers(excel_path, sheet_pre, sheet_post, columns,
                      nrows_pop, nrows_wt, columns_wt=None, skiprows=1,
                      aminoacids=list('AACDEFGGHIKLLLMNPPQRRRSSSTTVVWY*'),
                      zeroing='population', how='median', norm_std=True, 
                      stopcodon=False, min_counts=25, min_countswt=100, 
                      std_scale=0.2, mpop=2, mwt=2, infinite=3, **kwargs):
    '''
    Assembles different sublibraries into one. Uses calculate_enrichments. 
    Can only read from excel files that are in the same format as the example provided.

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

    **kwargs : other keyword arguments
        savefile : boolean, default False. optional kwarg 
            If set to true, the function will export the enrichment scores to a txt file.

    Returns
    --------
    df : Pandas dataframe
        A dataframe that contains the enrichment scores of the assembled sublibraries.
    '''
    # update kwargs
    temp_kwargs = copy.deepcopy(default_kwargs)
    temp_kwargs.update(kwargs)
    
    # Read reads from excel
    list_pre, list_sel, list_pre_wt, list_sel_wt = _read_counts(excel_path, sheet_pre, sheet_post, columns,
                                                                nrows_pop, nrows_wt, columns_wt)

    # Assemble sublibraries
    df = _assemble_list(list_pre, list_sel, list_pre_wt, list_sel_wt,
                        aminoacids, zeroing, how, norm_std, stopcodon, min_counts,
                        min_countswt, std_scale, mpop, mwt, infinite, **kwargs)
    
        # Export files
    if temp_kwargs['savefile']:
        np.savetxt(temp_kwargs['outputfilepath']+temp_kwargs['outputfilename'],
                   zeroed, fmt='%i', delimiter='\t')

    
    return df


def _read_counts(excel_path, sheet_pre, sheet_post, columns,
                 nrows_pop, nrows_wt, columns_wt=None, skiprows=1):
    '''Aux'''
    # Create dictionary with data. Loading 3 replicates, each of them is divided into 3 pools
    list_pre, list_sel, list_pre_wt, list_sel_wt = ([] for i in range(4))

    # Read counts from excel
    replicates = np.arange(0, len(sheet_pre))
    for column, column_wt, nrow_wt, rep in zip(columns, columns_wt, nrows_wt, replicates):
        # Pre counts
        list_pre.append(pd.read_excel(excel_path, sheet_pre,
                                      skiprows=skiprows, usecols=column, nrows=nrows_pop))
        # Sel counts
        list_sel.append(pd.read_excel(excel_path, sheet_post,
                                      skiprows=skiprows, usecols=column, nrows=nrows_pop))
        if columns_wt is None:
            list_pre_wt.append(None)
            list_sel_wt.append(None)
        else:
            # Pre counts wild-type alleles
            list_pre_wt.append(pd.read_excel(
                excel_path, sheet_pre, usecols=column_wt, nrows=nrow_wt))
            # Sel counts wild-type alleles
            list_sel_wt.append(pd.read_excel(
                excel_path, sheet_post, usecols=column_wt, nrows=nrow_wt))
    return list_pre, list_sel, list_pre_wt, list_sel_wt


def _assemble_list(list_pre, list_sel, list_pre_wt, list_sel_wt, aminoacids, zeroing,
                   how, norm_std, stopcodon, min_counts,
                   min_countswt, std_scale, mpop, mwt, infinite, **kwargs):
    '''gets the output from _read_counts and assembles the sublibraries'''
    
    enrichment_lib = []
    
    for pre, sel, pre_wt, sel_wt in zip(list_pre, list_sel, list_pre_wt, list_sel_wt):
        # log 10
        enrichment_log10 = calculate_enrichment(pre, sel, pre_wt, sel_wt, aminoacids, zeroing,
                                                how, norm_std, stopcodon, min_counts,
                                                min_countswt, std_scale, mpop, mwt, infinite, **kwargs)
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
    self : object from class "Screen"
    
    path : str
        Path where is located the fasta MSA that will be parsed. That MSA needs to have removed 
        any insertions that are not present in the target sequence. For example, if a Ras ortholog has 
        an extra amino acid at position 123, that needs to be removed from the aligment. Otherwise, everything
        will be shifted by 1 residue.
    
    start_position : int
        This is the position in the protein sequence of the first position in the MSA
    
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
    msa, seq_lengths, index = shannon.parseMSA(path,"fasta",0)

    # Calculate Shannon entropy from alignment
    shannon_entropy = shannon.shannon_entropy_list_msa(msa)
    
    # Merge enrichment scores and MSA conservation
    df_freq = _merge_msa_enrichment(self,_msa_to_df(msa),start_position,threshold)
    
    # Merge shannon and mean enrichment score
    df_shannon = _merge_shannon_enrichment(self,shannon_entropy,start_position)
    
    return  df_shannon, df_freq

def _merge_shannon_enrichment(self,shannon_entropy,start_position):
        
    # Create df with shannon entropy by residue and average enrichment score by residue
    df_shannon = pd.DataFrame()
    df_shannon['Position'] = np.arange(start_position, len(shannon_entropy)+start_position)
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

def _msa_to_df(msa,correctionfactor=1):
    '''Convert a msa from a fasta file into a df ready to plot with logomaker. Returns frequency'''
    # Flatten MSA
    msa_flattened = _flatten_msa(msa)

    # Make matrix
    df = logomaker.alignment_to_matrix(msa_flattened)

    # Reindex
    df.index = np.arange(correctionfactor,len(df)+correctionfactor)
    
    # Return only common aa 
    aminoacids = list ('ACDEFGHIKLMNPQRSTVWY')
    
    # Normalize by the total number of counts
    df_final = df[aminoacids].copy()
    
    return df_final/df_final.sum(axis=1).max()


def _merge_msa_enrichment(self, df_msa, start_position, threshold):
    '''merges msa conservation of each individual amino acid with the enrichment scores'''
    
    # make a dataframe
    df = pd.DataFrame()

    # Create column with position and aminoacid label
    df['Position'] = np.ravel([[i]*len(df_msa.T) for i in range(start_position, len(df_msa)+start_position)])
    df['Aminoacid'] = list(df_msa.columns) * len(df_msa)

    # Add conservation from MSA
    df['Conservation'] = list(df_msa.stack(dropna=False))

    # Merge with enrichment scores
    df_merged = self.dataframe.merge(df, how='inner', on=['Position', 'Aminoacid'])

    # Copycat conservation
    df_merged['Class'] = df_merged['Conservation']

    # Binarize conservation scores. 0 means not conserved
    df_merged.loc[df_merged['Conservation']>threshold, 'Class'] = 1
    df_merged.loc[df_merged['Conservation']<=threshold, 'Class'] = 0
    return df_merged
    


# # Plot Functions

# ## Kernel

# In[ ]:


def plot_kernel(self, kernel='gau', kernel_label='KDE', histogram=False,
                fit=None, fit_label='_nolegend_', extra_dist=None,
                extra_dist_label='_nolegend_',  **kwargs):
    '''
    Generate a kernel density plot. If specified it can also draw a histogram. Uses sns.distplot.

    Parameters
    ----------
    self : object from class "Screen"
    
    kernel : str, default gau
        options are ['biw','cos','epa','gau','tri','triw']
    
    kernel_label : str, default '_nolegend_'
    
    histogram : boolean, default False
    
    fit : boolean, optional
        ask sns.distplot to fit a function
    
    fit_label : str, default '_nolegend_'
    
    extra_dist : [x,y], optional
        fit any distribution you want. Input the x and y coordinates
    
    extra_dist_label : str, default '_nolegend_'
    
    **kwargs : other keyword arguments

    Returns
    ----------
    None
    '''
    # update kwargs
    temp_kwargs = copy.deepcopy(default_kwargs)
    temp_kwargs.update(kwargs)
    temp_kwargs['figsize'] = kwargs.get('figsize', (2.5, 2))

    # create figure
    fig = plt.figure(figsize=temp_kwargs['figsize'])

    # import parameters
    _parameters()

    # plot
    ax = sns.distplot(self.dataframe['Score_NaN'], kde=True, hist=histogram, norm_hist=True,
                      kde_kws={'kernel': kernel, 'color': temp_kwargs['color'], 'lw': 2, 'label': kernel_label}, fit=fit,
                      fit_kws={'label': fit_label, 'linestyle': 'dotted', 'color': 'red'})

    # plot extra distribution
    if extra_dist is not None:
        plt.plot(extra_dist[0], extra_dist[1], linewidth=2, linestyle='dotted',
                 color='green', label=extra_dist_label)

    # tune graph
    plt.xlabel(r'$∆E^i_x$', fontsize=10,
               fontname='Arial', color='k', labelpad=0)
    plt.ylabel('Probability density', fontsize=10,
               fontname='Arial', color='k', labelpad=3)
    plt.title(temp_kwargs['title'], fontsize=12, fontname='Arial', color='k')
    plt.xlim(temp_kwargs['xscale'])
    plt.grid()
    ax.legend(loc='best', frameon=False, fontsize=9,
              handlelength=1, handletextpad=0.5)

    # save file
    _savefile(fig, temp_kwargs)

    if temp_kwargs['show']:
        plt.show()

    return


def _parameters():
    # normal font
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['Arial']

    # math font
    rcParams['mathtext.fontset'] = 'custom'
    rcParams['mathtext.rm'] = 'Arial'
    rcParams['svg.fonttype'] = 'none'

    # add grid
    rcParams['grid.color'] = 'silver'
    rcParams['grid.linestyle'] = '--'
    rcParams['grid.linewidth'] = 1
    rcParams['lines.dashed_pattern'] = [5, 10]
    rcParams['axes.axisbelow'] = True
    # Parameters for all graphs
    rcParams['xtick.labelsize'] = 9
    rcParams['ytick.labelsize'] = 9
    return


# ## Multiple Kernel plots

# In[ ]:


def plot_multiplekernel(dict_entries, kernel='gau',
                        colors=['k', 'crimson', 'dodgerblue', 'g', 'silver'], **kwargs):
    '''
    Generate a kernel density plot for multiple objects passed as a dictionary.
    If specified it can also draw a histogram. Uses sns.distplot. Can manage either Screen objects
    or dataframes out of the calculate_enrichments function.

    Parameters
    ----------
    dict_entries : dictionary containing dataframes
        Allows for either putting multiple objects as inputs or to use dataframes 
        that come out of the calculate_enrichments function. If you use an object, 
        you need to say object.dataframe.

    kernel : str, default gau
        options are ['biw','cos','epa','gau','tri','triw'].

    colors : list, default ['k', 'crimson', 'dodgerblue', 'g', 'silver']
        List of the colors (in order of arguments) that the kernels will have.

    **kwargs : other keyword arguments

    Returns
    ----------
    None
    '''

    # update kwargs
    temp_kwargs = copy.deepcopy(default_kwargs)
    temp_kwargs.update(kwargs)
    temp_kwargs['figsize'] = kwargs.get('figsize', (3.5, 2))
    temp_kwargs['xscale'] = kwargs.get('xscale', (-2, 2))

    # hard copy of data
    dict_copy = copy.deepcopy(dict_entries)

    # create figure
    fig = plt.figure(figsize=temp_kwargs['figsize'])

    # import parameters
    _parameters()

    # plot (allows two types of data input)
    for (label, dataset, color) in zip(dict_copy.keys(), dict_copy.values(), colors[0:len(dict_copy)]):
        if 'Score' in dataset.columns:
            # plot objects scores
            sns.distplot(dataset['Score_NaN'], hist=False,
                         kde_kws={"color": color, "lw": 2, "label": label})
        else:
            # get rid of stop codons
            dataset.drop('*', errors='ignore', inplace=True)
            dataset = dataset.stack()
            # plot stacked matrix
            sns.distplot(dataset[~numpy.isnan(dataset)], kde=True, hist=False, 
                         kde_kws={"color": color, "lw": 2, "label": label})

    # tune graph
    plt.xlabel(r'$∆E^i_x$', fontsize=10,
               fontname='Arial', color='k', labelpad=0)
    plt.ylabel('Probability density', fontsize=10,
               fontname='Arial', color='k', labelpad=3)
    plt.title(temp_kwargs['title'], fontsize=12, fontname='Arial', color='k')
    plt.xlim(temp_kwargs['xscale'])
    plt.grid()
    plt.legend(dict_copy.keys(), loc='best', frameon=False, fontsize=9,
               handlelength=1, handletextpad=0.5)

    # save file
    _savefile(fig, temp_kwargs)

    if temp_kwargs['show']:
        plt.show()

    return


# ## Heatmap

# ### Full Heatmap

# In[ ]:


def plot_heatmap(self, nancolor='lime', show_cartoon=False, show_snv = False, **kwargs):
    '''
    Generate a heatmap plot of the enrichment scores.

    Parameters
    ----------
    self : object from class "Screen"
    
    nancolor : str, default 'lime'
        Will color np.nan values with the specified color.
    
    show_carton : boolean, default False
        If true, the plot will display a cartoon with the secondary structure. The user must have added the secondary structure to the object. 
    
    show_snv : boolean, default False
        If true, it will only display mutants that are a single nucleotide variant (SNV) of the wild-type
        protein sequence. The algorithm does not take into account the wild-type DNA allele, so it 
        will include any possible mutant that is one base away.
        
    **kwargs : other keyword arguments
        
    Returns
    ----------
    None    
    '''
    # load font parameters
    _font_parameters()

    # update kwargs
    temp_kwargs = copy.deepcopy(default_kwargs)
    temp_kwargs.update(kwargs)

    # load labels
    temp_kwargs['color_sequencelabels'] = _labels(self.start_position)[0]
    temp_kwargs['number_sequencelabels'] = _labels(self.start_position)[1]

    # sort data in specified order by user
    df = _df_rearrange(_add_SNV_boolean(self.dataframe_stopcodons.copy()),
                            temp_kwargs['neworder_aminoacids'], values='Score_NaN', show_snv = show_snv)

    # declare figure and subplots
    figwidth = 14*len(df.columns)/165

    # Change parameters depending on whether cartoon is on or off
    if show_cartoon:
        figheight = 2.45
        fig = plt.figure(figsize=(figwidth, figheight))
        gs = gridspec.GridSpec(nrows=3, ncols=2, height_ratios=[
                               len(df), 1, 5], width_ratios=[len(df.columns), 1])
    else:
        figheight = 2
        fig = plt.figure(figsize=(figwidth, figheight))
        gs = gridspec.GridSpec(nrows=2, ncols=2, height_ratios=[
                               len(df), 1], width_ratios=[len(df.columns), 1])

    ax = plt.subplot(gs[0, 0])
    averageresidue = plt.subplot(gs[1, 0])
    cbar1 = plt.subplot(gs[0, 1])
    cbar2 = plt.subplot(gs[1, 1])

    # Change color of values that are NaN
    cmap = temp_kwargs['colormap']
    cmap.set_bad(color=nancolor)

    # main heatmap
    heatmap = ax.pcolormesh(df, vmin=temp_kwargs['colorbar_scale'][0], vmax=temp_kwargs['colorbar_scale'][1],
                            cmap=cmap, edgecolors='k', linewidths=0.2, antialiased=True, color='darkgrey')

    # average of residues by positon
    average = [_add_SNV_boolean(self.dataframe.copy()).groupby(by='Position').mean()['Score_NaN']]

    # average by position
    heatmapaverageresidues = averageresidue.pcolormesh(average, vmin=temp_kwargs['colorbar_scale'][0], vmax=temp_kwargs['colorbar_scale'][1],
                        cmap=cmap, edgecolors='k', linewidths=0.2, antialiased=True, color='darkgrey')

    # ____________axes manipulation____________________________________________
    # put the major ticks at the middle of each cell
    ax.set_xticks(np.arange(len(df.columns)) + 0.5, minor=False)
    ax.set_yticks(np.arange(len(df)) + 0.5, minor=False)

    # position of axis labels
    ax.tick_params('x', direction='out', pad=-2.5)
    ax.tick_params('y', direction='out', pad=0.4)

    # make new axes
    ax2 = ax.twiny()
    ax3 = ax.twinx()

    # tune the axes
    ax2.set_xticks(np.arange(len(df.columns)) + 0.5, minor=False)
    ax3.set_yticks(np.arange(len(df)) + 0.5, minor=False)
    ax2.tick_params(direction='out', pad=4)
    ax3.tick_params(direction='out', pad=0.4)
    averageresidue.tick_params(direction='out', pad=-2)
    averageresidue.set_xticks(np.arange(len(df.columns)) + 0.5, minor=False,)
    averageresidue.set_yticks(np.arange(0.5)+0.5)

    # Set the limits of the new axis from the original axis limits
    ax2.set_xlim(ax.get_xlim())
    ax3.set_ylim(ax.get_ylim())

    # want a more natural, table-like display
    ax.invert_yaxis()
    ax.xaxis.tick_top()
    ax3.invert_yaxis()

    # remove ticks
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')
    ax2.yaxis.set_ticks_position('none')
    ax2.xaxis.set_ticks_position('none')
    ax3.yaxis.set_ticks_position('none')
    averageresidue.xaxis.set_ticks_position('none')
    averageresidue.yaxis.set_ticks_position('none')

    # so labels of x and y do not show up and my labels show up instead
    ax.set_xticklabels(list(self.sequence), fontsize=6.5,
                       fontname="Arial", color='k', minor=False)
    ax.set_yticklabels(temp_kwargs['neworder_aminoacids'], fontsize=6,
                       fontname="Arial", color='k', minor=False)
    ax2.set_xticklabels(temp_kwargs['number_sequencelabels'][0:len(df.columns)],
                        fontsize=10, fontname="Arial", color='k', minor=False)
    ax3.set_yticklabels(temp_kwargs['neworder_aminoacids'],
                        fontsize=6, fontname="Arial", color='k', minor=False)
    averageresidue.set_xticklabels(list(self.sequence), fontsize=6.5,
                                   fontname="Arial", color='k', minor=False)
    rowaverage = ''
    averageresidue.set_yticklabels(
        rowaverage, fontsize=6, fontname="Arial", color='k', minor=False)

    # align the labels of the y axis
    for ylabel in ax.get_yticklabels():
        ylabel.set_horizontalalignment('center')
    for ylabel in ax3.get_yticklabels():
        ylabel.set_horizontalalignment('center')

    # for coloring the residues that are 10,20...
    for xtick, color in zip(ax.get_xticklabels(), temp_kwargs['color_sequencelabels']):
        xtick.set_color(color)
    for xtick, color in zip(averageresidue.get_xticklabels(), temp_kwargs['color_sequencelabels']):
        xtick.set_color(color)
    # _____________________________________________________________________________

    # for color bar format
    cbar1.axis('off')
    cbar2.axis('off')
    cb = plt.colorbar(heatmap, fraction=1, pad=0, ax=[cbar1], aspect=5, ticks=[temp_kwargs['colorbar_scale'][0], np.mean(
        temp_kwargs['colorbar_scale']), temp_kwargs['colorbar_scale'][1]], orientation='vertical')
    cb.ax.set_yticklabels(cb.ax.get_yticklabels(),
                          fontsize=6, fontname="Arial", color='k')
    cb.update_ticks()
    plt.text(1.2+10/len(df.columns), 0.7, r'$\langle∆E^x_i\rangle_x$', transform=cbar1.transAxes,
             horizontalalignment='center', fontsize=7, fontname="Arial", color='k')

    gs.update(hspace=0.1, wspace=0.1/len(df.columns)*50)

    # for putting title on graph
    plt.title(temp_kwargs['title'], horizontalalignment='center',
              fontname="Arial", fontsize=12)

    # Cartoon
    if show_cartoon:
        _generate_cartoon(self, gs, 2, temp_kwargs['cartoon_colors'], 0.025)

    # save file
    _savefile(fig, temp_kwargs)
    plt.show()
    return


def _labels(start_position=1):
    # residue label and color
    emptylist = [''] * 1000
    number_sequencelabels = list(['b' if index in np.arange(
        10-(start_position % 10), 1000, 10) else 'k' for index, x in enumerate(emptylist)])
    color_sequencelabels = list([index+start_position if index in np.arange(10 -
                                                                            (start_position % 10), 1000, 10) else '' for index, x in enumerate(emptylist)])
    return number_sequencelabels, color_sequencelabels


def _font_parameters():
    # math font
    rcParams['mathtext.fontset'] = 'custom'
    rcParams['mathtext.rm'] = 'Arial'
    rcParams['svg.fonttype'] = 'none'
    return


def generatecolormap():
    cmap = LinearSegmentedColormap(
        'BlueWhiteRed',
        {
            'red':  ((0.0, 0.0, 0.0),
                     (0.15, 0.0, 0.0),
                     (0.475, 1.0, 1), (0.525, 1.0, 1), (0.85, 1.0, 1.0), (1.0, .8, 1)),
            'green': ((0.0, 0.0, .0),
                      (0.15, 0.5, 0.5), (0.475, 1.0,
                                         1), (0.525, 1.0, 1), (0.85, 0.0, 0.0),
                      (1.0, 0.0, 0.0)),
            'blue': ((0.0, .5, .5),
                     (0.15, 1, 1),
                     (0.475, 1.0, 1), (0.525, 1.0, 1), (0.85, 0.0, 0.0), (1.0, 0.0, 0.0))
        },)
    return cmap


def _generate_cartoon(self, gs, n_row, colors, bottom_space=0,
                      fig_inches=13.91, show_labels=True):
    '''Generates cartoon for heatmap'''
    # Create subplot
    cartoon = plt.subplot(gs[n_row, 0])

    # Generate coordinates of labels
    labels = list(Counter(self.secondary).keys())
    length = list(Counter(self.secondary).values())
    cumsum = length[:-1]
    cumsum.insert(0, self.start_position)
    cumsum = np.cumsum(cumsum)

    # Create cartoon
    for label, length, cum in zip(labels, length, cumsum):
        if 'β' in label:
            loopstructure = _loop(cum, length, color=colors[2])
            cartoon.add_patch(loopstructure)
            sheetstructure = _sheet(cum, length, colors[0])
            cartoon.add_patch(sheetstructure)
            x_label = cum + length - 3.5
            if length > 2 and show_labels:  # If beta sheet is too small, label does not fit
                if length == 3:
                    cartoon.text((x_label+0.6), -0.25, label, name='Arial',
                                 fontweight='normal', size=8.5*fig_inches/13.91, multialignment='right')
                else:
                    cartoon.text((x_label), -0.25, label, name='Arial', fontweight='normal',
                                 size=8.5*fig_inches/13.91, multialignment='right')
        elif 'α' in label:
            helixstructure = _helix(cum, length, colors[1])
            cartoon.add_patch(helixstructure)
            x_label = cum + length/2 - 1
            if length > 2 and show_labels:
                cartoon.text((x_label), -0.3, label, name='Arial', fontweight='normal',
                             size=9*fig_inches/14, multialignment='center')
        elif 'L' in label:
            loopstructure = _loop(cum, length, colors[2])
            cartoon.add_patch(loopstructure)

    # format of secondary cartoon
    cartoon.xaxis.set_ticks_position('none')
    cartoon.yaxis.set_ticks_position('none')
    cartoon.axis('off')

    # size
    cartoon.set_xlim(self.start_position-0.1,
                     len(self.secondary)+self.start_position+0.2)
    cartoon.set_ylim(-2, 2.5)

    # adjust proximity to heatmap
    box = cartoon.get_position()
    box.y0 = box.y0-bottom_space
    box.y1 = box.y1-bottom_space
    cartoon.set_position(box)

    return


def _sheet(starting_aa, length_aa, color='lightgreen'):
    dx = length_aa
    sheetstructure = patches.FancyArrow(starting_aa, 0.25, dx, 0, width=2, length_includes_head=True,
                                        head_width=4, head_length=3, shape='full', overhang=0, head_starts_at_zero=False, ec='k', fc=color)
    return sheetstructure


def _helix(starting_aa, length_aa, color='lavender'):
    dx = length_aa  # so i can overlap tip
    helixstructure = plt.Rectangle(
        (starting_aa, -0.85), dx, 2.2, fc=color, ec='k')
    return helixstructure


def _loop(starting_aa, length_aa, color='k'):
    dx = length_aa
    loopstructure = plt.Rectangle((starting_aa, 0), dx, 0.5, fc=color)
    return loopstructure


# ### Grouped Heatmap

# In[ ]:


def plot_heatmap_rows(self, selection=['E', 'Q', 'A', 'P', 'V', 'Y'],
                           nancolor='lime', **kwargs):
    '''
    Generate a heatmap plot enrichment scores of selected aminoacids. 

    Parameters
    ----------
    self : object from class "Screen"
    
    selection : list of aa to show, default ['E','Q','A','P','V','Y']. 
    
    nancolor : str, default 'lime'
        Will color np.nan values with the specified color.
    
    **kwargs : other keyword arguments

    Returns
    ----------
    None.
    '''
    # load font parameters
    _font_parameters()

    # update kwargs
    temp_kwargs = copy.deepcopy(default_kwargs)
    temp_kwargs.update(kwargs)

    # load labels
    temp_kwargs['color_sequencelabels'] = _labels(self.start_position)[0]
    temp_kwargs['number_sequencelabels'] = _labels(self.start_position)[1]

    # Add group and pivot df
    df = _select_aa(self.dataframe_stopcodons, selection, values='Score_NaN')
    dataset = df.to_numpy()

    # The size can be changed. I found it empirically
    figwidth = 14*len(dataset[0])/165
    figheight = 2/21*len(selection)
    fig = plt.figure(figsize=(figwidth, figheight))
    gs = gridspec.GridSpec(nrows=1, ncols=2, width_ratios=[len(dataset[0]), 1])
    ax = plt.subplot(gs[0, 0])
    cbar1 = plt.subplot(gs[0, 1])

    # Change color of values that are NaN
    cmap = temp_kwargs['colormap']
    cmap.set_bad(color=nancolor)

    # main heatmap
    heatmap = ax.pcolormesh(dataset, vmin=temp_kwargs['colorbar_scale'][0], vmax=temp_kwargs['colorbar_scale'][1],
                            cmap=cmap, edgecolors='k', linewidths=0.2, antialiased=True, color='darkgrey')

    # put the major ticks at the middle of each cell
    ax.set_xticks(np.arange(dataset.shape[1]) + 0.5, minor=False)
    ax.set_yticks(np.arange(dataset.shape[0]) + 0.5, minor=False)

    # position of axis labels
    ax.tick_params('x', direction='out', pad=-2.5)
    ax.tick_params('y', direction='out', pad=0.4)

    # second axis
    ax2 = ax.twiny()
    ax2.set_xticks(np.arange(dataset.shape[1]) + 0.5, minor=False)
    ax2.tick_params(direction='out', pad=4)

    # Set the limits of the new axis from the original axis limits
    ax2.set_xlim(ax.get_xlim())

    # want a more natural, table-like display
    ax.invert_yaxis()
    ax.xaxis.tick_top()

    # so labels of x and y do not show up and my labels show up instead
    ax.set_xticklabels(list(self.sequence), fontsize=6.5,
                       fontname="Arial", color='k', minor=False)
    ax.set_yticklabels(list(df.T.columns), fontsize=6,
                       fontname="Arial", color='k', minor=False)
    ax2.set_xticklabels(temp_kwargs['number_sequencelabels'][0:len(dataset[0])],
                        fontsize=10, fontname="Arial", color='k', minor=False)

    # align the labels of the y axis
    for ylabel in ax.get_yticklabels():
        ylabel.set_horizontalalignment('center')

    # for coloring the residues that are 10,20...
    for xtick, color in zip(ax.get_xticklabels(), temp_kwargs['color_sequencelabels']):
        xtick.set_color(color)

    # for putting title on graph
    plt.title(temp_kwargs['title'], horizontalalignment='center',
              fontname="Arial", fontsize=12)

    # for color bar format
    cbar1.axis('off')
    cb = plt.colorbar(heatmap, fraction=1, pad=0, ax=[cbar1], aspect=5, ticks=[temp_kwargs['colorbar_scale'][0], np.mean(
        temp_kwargs['colorbar_scale']), temp_kwargs['colorbar_scale'][1]], orientation='vertical')
    cb.ax.set_yticklabels(cb.ax.get_yticklabels(),
                          fontsize=8, fontname="Arial", color='k')
    cb.update_ticks()
    gs.update(hspace=0.1, wspace=0.1/len(dataset[0])*50)

    # remove ticks
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')
    ax2.yaxis.set_ticks_position('none')
    ax2.xaxis.set_ticks_position('none')

    # save file
    _savefile(fig, temp_kwargs)

    if temp_kwargs['show']:
        plt.show()
    return


# ### Subset Heatmap

# In[ ]:


def plot_heatmap_columns(self, segment, ylabel_color='k', nancolor='lime', **kwargs):
    '''
    Generate a heatmap plot enrichment scores but only plots a selected segment.

    Parameters
    ----------
    self : object from class "Screen"
    
    segment : list
        Segment is typed as [20,40] and includes both residues 20 and 40.
    
    ylabel_color : str, default 'k'
        Choose white if you don't want amino acid y axis label.
    
    nancolor : str, default 'lime'
        Will color np.nan values with the specified color.
    
    **kwargs : other keyword arguments

    Returns
    ----------
    None    
    '''

    # load font parameters
    _font_parameters()

    # update kwargs
    temp_kwargs = copy.deepcopy(default_kwargs)
    temp_kwargs.update(kwargs)

    # load labels
    temp_kwargs['color_sequencelabels'] = _labels(self.start_position)[0]
    temp_kwargs['number_sequencelabels'] = _labels(self.start_position)[1]

    # sort data in specified order by user
    df_whole = _df_rearrange(self.dataframe_stopcodons,
                            temp_kwargs['neworder_aminoacids'], values='Score_NaN')
    
    # select subset
    c0 = segment[0]-self.start_position
    c1 = segment[1]-self.start_position+1
    df = df_whole.iloc[:,c0:c1]

    # the size can be changed
    figwidth = 2*len(df.columns)/22
    figheight = 2
    fig = plt.figure(figsize=(figwidth, figheight))
    gs = gridspec.GridSpec(nrows=1, ncols=1)
    # needed to set autoscale off to avoid missalignment
    ax = plt.subplot(gs[0])

    # Change color of values that are NaN
    cmap = temp_kwargs['colormap']
    cmap.set_bad(color=nancolor)

    # main heatmap
    heatmap = ax.pcolormesh(df, vmin=temp_kwargs['colorbar_scale'][0], vmax=temp_kwargs['colorbar_scale'][1],
                            cmap=cmap, edgecolors='k', linewidths=0.2, antialiased=True, color='darkgrey')

    # put the major ticks at the middle of each cell
    ax.set_xticks(np.arange(len(df.columns)) + 0.5, minor=False,)
    ax.set_yticks(np.arange(len(df)) + 0.5, minor=False)

    # position of axis labels
    ax.tick_params('x', direction='out', pad=-2.5)
    ax.tick_params('y', direction='out', pad=0.4)

    # second axis
    ax2 = ax.twiny()
    ax2.set_xticks(np.arange(len(df.columns)) + 0.5, minor=False)
    ax2.tick_params(direction='out', pad=4)

    # Set the limits of the new axis from the original axis limits
    ax2.set_xlim(ax.get_xlim())

    # want a more natural, table-like display
    ax.invert_yaxis()
    ax.xaxis.tick_top()

    # so labels of x and y do not show up and my labels show up instead
    ax.set_xticklabels(list(self.sequence)[segment[0]-self.start_position:segment[1] -
                                           self.start_position+1], fontsize=6.5, fontname="Arial", color='k', minor=False)
    ax.set_yticklabels(temp_kwargs['neworder_aminoacids'], fontsize=6,
                       fontname="Arial", color=ylabel_color, minor=False)

    ax2_label = (segment[1]-segment[0]+1)*['']
    ax2_label[0] = segment[0]
    ax2_label[-1] = segment[1]
    ax2.set_xticklabels(ax2_label, fontsize=7,
                        fontname="Arial", color='k', minor=False)

    # align the labels of the y axis
    for ylabel in ax.get_yticklabels():
        ylabel.set_horizontalalignment('center')

    # remove ticks
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')
    ax2.yaxis.set_ticks_position('none')
    ax2.xaxis.set_ticks_position('none')

    # save file
    _savefile(fig, temp_kwargs)

    if temp_kwargs['show']:
        plt.show()
    return 


# ## Mean Plots

# ### Bar graph Enrichment

# In[ ]:


def plot_mean(self, mode='mean', show_cartoon=False, **kwargs):
    '''
    Plot in a bargraph the mean enrichment for each residue of the protein. Red for gain of function, blue for loss of function

    Parameters
    ----------
    self : object from class "Screen"
    
    mode : str, default 'mean'
        Specify what enrichment scores to show. If mode = 'mean', it will show the mean of 
        each position. If mode = 'A', it will show the alanine substitution profile. Can be 
        used for each amino acid. Use the one-letter code and upper case.
        
    show_carton : boolean, default False
        If true, the plot will display a cartoon with the secondary structure. The user must have added the secondary structure to the object. 
    
    **kwargs : other keyword arguments

    Returns
    ----------
    None.
    '''
    # update kwargs
    temp_kwargs = copy.deepcopy(default_kwargs)
    temp_kwargs.update(kwargs)
    temp_kwargs['figsize'] = kwargs.get('figsize', (3, 2.5))
    temp_kwargs['yscale'] = kwargs.get('yscale', (-1, 1))
    temp_kwargs['y_label'] = kwargs.get('y_label', r'$∆E^i_x$')

    # load parameters
    parameters_mean()

    # Select grouping
    if mode == 'mean':
        df = self.dataframe.groupby('Position', as_index=False).mean()
    else:
        df = self.dataframe.loc[self.dataframe['Aminoacid']==mode].copy()
    
    df['Color'] = df.apply(color_data, axis=1)

    # make figure
    if show_cartoon:
        fig = plt.figure(figsize=temp_kwargs['figsize'])
        gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[5, 1])
        ax = plt.subplot(gs[0])
    else:
        fig, ax = plt.subplots(figsize=temp_kwargs['figsize'])
    width = 1.2

    # Color based on values
    ax.bar(df['Position'], df['Score'], width, color=df['Color'], snap=False)

    # axes parameters
    ax.set_ylim(temp_kwargs['yscale'])
    ax.set_ylabel(temp_kwargs['y_label'], fontsize=10,
                  fontname="Arial", color='k', labelpad=10, rotation=0)
    ax.set_xticks(np.arange(self.start_position,
                            len(df)+self.start_position, 20))
    ax.set_xlabel('Residue', fontsize=10,
                  fontname="Arial", color='k', labelpad=4)
    ax.set_xlim(self.start_position-0.1, len(df)+self.start_position-1+0.1)
    plt.title(temp_kwargs['title'], fontsize=12, fontname='Arial', color='k')

    # cartoon
    if show_cartoon:
        _generate_cartoon(self, gs, 1, temp_kwargs['cartoon_colors'],
                          bottom_space=-0.78, show_labels=False)
    # Put text labels
    _inputtext(temp_kwargs['text_labels'])

    # save file
    _savefile(fig, temp_kwargs)
    if temp_kwargs['show']:
        plt.show()
    return


def color_data(row):
    if row['Score'] > 0:
        return 'red'
    else:
        return 'blue'


def parameters_mean():
    # normal font
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['Arial']

    # math font
    rcParams['mathtext.fontset'] = 'custom'
    rcParams['mathtext.rm'] = 'Arial'
    rcParams['svg.fonttype'] = 'none'

    # Parameters for all graphs
    rcParams['xtick.labelsize'] = 9
    rcParams['ytick.labelsize'] = 9

    return


# ### Compare two proteins

# In[ ]:


def plot_meandifferential (self, obj2, show_cartoon=False,**kwargs):
    '''
    Plot the mean positional difference between two experiments

    Parameters
    ----------
    self : object from class "Screen"
    
    obj2 : another Screen object to compare with 
    
    show_carton : boolean, default False
        If true, the plot will display a cartoon with the secondary structure. The user must have added the secondary structure to the object. 
    
    **kwargs : other keyword arguments

    Returns
    ----------
    None.
    '''
    # update kwargs
    temp_kwargs = copy.deepcopy(default_kwargs)
    temp_kwargs.update(kwargs)
    temp_kwargs['figsize'] = kwargs.get('figsize', (3,2.5))
    temp_kwargs['yscale'] = kwargs.get('yscale', (-1,1))
    temp_kwargs['y_label'] = kwargs.get('y_label', r'Mean Differential $∆E^i_x$')

    # load parameters
    parameters_mean()
    
    # make pandas
    df = _process_meanresidue(self,obj2)
    
    # make cartoon
    if show_cartoon:
        fig = plt.figure(figsize=temp_kwargs['figsize'])
        gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[5, 1])
        ax = plt.subplot(gs[0])
    else:
        fig, ax = plt.subplots(figsize=temp_kwargs['figsize']) 

    # plot
    ax.plot(df['Position'], df['d1 - d2'], color='k')

    # axes parameters
    ax.set_ylim(temp_kwargs['yscale'])
    ax.set_ylabel(temp_kwargs['y_label'], fontsize=10, fontname="Arial", 
                  color='k', labelpad=-5, rotation=90)
    ax.set_xticks(np.arange(self.start_position, len(df)+self.start_position, 20))
    ax.set_xlabel('Residue', fontsize=10, fontname="Arial", color='k', labelpad=4)
    ax.set_xlim(self.start_position-0.1, len(df)+self.start_position-1+0.1)
    ax.set_title(temp_kwargs['title'], fontsize=12, fontname='Arial', color='k')

    # cartoon
    if show_cartoon:
        obj = obj2
        if len(self.dataframe) < len(obj2.dataframe):
            obj = self
        _generate_cartoon(obj,gs,1,temp_kwargs['cartoon_colors'],
                            bottom_space=-0.78, show_labels=False)
    # save file
    _savefile(fig,temp_kwargs)

    if temp_kwargs['show']: plt.show()
    return


# ### Bar graph Counts

# In[ ]:


def plot_meancounts (self, positions, counts, show_cartoon=False, **kwargs):
    '''
    Plot in a bargraph the mean counts for each residue of the protein.
    
    Parameters
    ----------
    self : object from class "Screen"
    
    positions : list, x coordinates
    
    counts : list, y coordinates
    
    show_carton : boolean, default False
        If true, the plot will display a cartoon with the secondary structure. The user must have added the secondary structure to the object. 
    
    **kwargs : other keyword arguments

    Returns
    ----------
    None.
    '''
    # update kwargs
    temp_kwargs = copy.deepcopy(default_kwargs)
    temp_kwargs.update(kwargs)
    temp_kwargs['figsize'] = kwargs.get('figsize', (3,2.5))
    temp_kwargs['yscale'] = kwargs.get('yscale', (0,5))
    temp_kwargs['y_label'] = kwargs.get('y_label', r'$Log_{10}$ mean counts')

    # load parameters
    parameters_mean()
        
    # make figure
    if show_cartoon:
        fig = plt.figure(figsize=temp_kwargs['figsize'])
        gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[5, 1])
        ax = plt.subplot(gs[0])
    else:
        fig, ax = plt.subplots(figsize=temp_kwargs['figsize']) 
    width = 0.8

    # Color based on values
    ax.bar(positions, np.log10(counts), width, color='red', snap=False)

    # axes parameters
    ax.set_ylim(temp_kwargs['yscale'])
    ax.set_ylabel(temp_kwargs['y_label'], fontsize=10, fontname="Arial", color='k', labelpad=0, rotation=90)
    ax.set_xticks(np.arange(self.start_position, len(self.dataset[0])+self.start_position, 20))
    ax.set_xlabel('Residue', fontsize=10, fontname="Arial", color='k', labelpad=4)
    ax.set_xlim(self.start_position-0.1, len(self.dataset[0])+self.start_position-1+0.1)
    plt.title(temp_kwargs['title'], fontsize=12, fontname='Arial', color='k')
    
    # cartoon
    if show_cartoon:
        _generate_cartoon(self,gs,1,temp_kwargs['cartoon_colors'],
                            bottom_space=-0.78, show_labels=False)
    # Put text labels
    _inputtext(temp_kwargs['text_labels'])
    
    # save file
    _savefile(fig,temp_kwargs)
    if temp_kwargs['show']: plt.show()
    return

def _inputtext(text_entries):
    '''the user can input text as a variable by manually giving the coordinates'''
    if text_entries:
        for entry in text_entries:
            plt.text(entry[0],entry[1],entry[2])
    return


# ### Positional

# In[ ]:


def plot_position(self, position, **kwargs):
    '''
    Choose a position and plot in a bargraph the enrichment score for each substitution.
    Red for gain of function, blue for loss of function.

    Parameters
    ----------
    self : object from class "Screen"
    
    position : int
        number of residue of the protein to display.
    
    **kwargs : other keyword arguments

    Returns
    ----------
    None.
    '''
    # update kwargs
    temp_kwargs = copy.deepcopy(default_kwargs)
    temp_kwargs.update(kwargs)
    temp_kwargs['figsize'] = kwargs.get('figsize', (3.5, 2))
    temp_kwargs['yscale'] = kwargs.get('yscale', (-1, 1))
    temp_kwargs['y_label'] = kwargs.get('y_label', r'$∆E^i_x$')

    # load parameters
    parameters_mean()

    # Select position
    df = self.dataframe.loc[self.dataframe['Position']==position].copy()
    
    # Color
    df['Color'] = df.apply(color_data, axis=1)

    # make figure
    fig, ax = plt.subplots(figsize=temp_kwargs['figsize'])
    width = 0.5

    # Color based on values
    ax.bar(df['Aminoacid'], df['Score'], width, color=df['Color'], ec='k')

    # axes parameters
    ax.set_ylim(temp_kwargs['yscale'])
    ax.set_ylabel(temp_kwargs['y_label'], fontsize=10,
                  fontname="Arial", color='k', labelpad=10, rotation=0)

    ax.set_xlabel('Residue', fontsize=10,
                  fontname="Arial", color='k', labelpad=4)
    plt.title(temp_kwargs['title'], fontsize=12, fontname='Arial', color='k')

    # save file
    _savefile(fig, temp_kwargs)
    
    if temp_kwargs['show']:
        plt.show()
    return


# ## Scatter

# In[ ]:


def plot_scatter(self, obj2, mode='pointmutant', **kwargs):
    '''
    Generate a scatter plot between object and a second object of the same class.

    Parameters
    ----------
    self : object from class "Screen"
    
    obj2 : object from class "Screen" to do the scatter with
    
    mode : str, default 'pointmutant'. 
        Alternative set to "mean" for the mean of each position
    
    **kwargs : other keyword arguments

    Returns
    ----------
    None.
    '''

    # update kwargs
    temp_kwargs = copy.deepcopy(default_kwargs)
    temp_kwargs.update(kwargs)
    temp_kwargs['figsize'] = kwargs.get('figsize', (2, 2))
    temp_kwargs['xscale'] = kwargs.get('xscale', (-2, 2))
    temp_kwargs['yscale'] = kwargs.get('yscale', (-2, 2))

    # Chose mode:
    if mode == 'pointmutant':
        df = _process_bypointmutant(self, obj2)
    else:
        df = _process_meanresidue(self, obj2)

    # create figure
    fig, ax = plt.subplots(figsize=temp_kwargs['figsize'])

    # import parameters
    _parameters()

    # Scatter data points
    plt.scatter(df['dataset_1'], df['dataset_2'], c='k', s=8,
                alpha=0.5, rasterized=True, label='_nolegend_')

    # Titles
    plt.title(temp_kwargs['title'], fontsize=12,
              fontname='Arial', color='k', pad=8)
    plt.ylabel(temp_kwargs['y_label'], fontsize=10,
               fontname="Arial", color='k', labelpad=0)
    plt.xlabel(temp_kwargs['x_label'], fontsize=10,
               fontname="Arial", color='k')

    # correlation and R2
    slope, intercept, r_value, p_value, std_err = stats.linregress(
        df['dataset_1'], df['dataset_2'])
    R2 = str(round(r_value**2, 2))
    legend_label = "$R^2$ = {}".format(R2)
    # fit and graph line
    fit = np.polyfit(df['dataset_1'], df['dataset_2'], 1)
    plt.plot(np.unique(df['dataset_1']), np.poly1d(fit)(
        np.unique(df['dataset_1'])), color='r', linewidth=1, label=legend_label)
    plt.grid()

    # other graph parameters
    plt.xlim(temp_kwargs['xscale'])
    plt.ylim(temp_kwargs['yscale'])
    ax.xaxis.set_major_locator(
        ticker.MultipleLocator(temp_kwargs['tick_spacing']))
    ax.yaxis.set_major_locator(
        ticker.MultipleLocator(temp_kwargs['tick_spacing']))
    plt.gca().set_aspect('equal', adjustable='box')
    plt.draw()

    # Legend
    plt.legend(loc='upper left', handlelength=0,
               handletextpad=0, frameon=False, fontsize=10)

    # save file
    _savefile(fig, temp_kwargs)
    
    if temp_kwargs['show']:
        plt.show()
        
    return


def _process_bypointmutant(self, obj):
    # truncate so both datasets have same length and delete stop codons
    minlength = min(len(self.dataframe), len(obj.dataframe))
    df = pd.DataFrame()
    df['dataset_1'] = list(self.dataframe['Score_NaN'])[:minlength]
    df['dataset_2'] = list(obj.dataframe['Score_NaN'])[:minlength]

    # eliminate Nans
    df.dropna(how='any', inplace=True)
    return df


def _process_meanresidue(self, obj):
    # truncate so both datasets have same length and delete stop codons
    dataset_1 = self.dataframe.groupby(['Position'], as_index=False).mean()
    dataset_2 = obj.dataframe.groupby(['Position'], as_index=False).mean()
    minlength = min(len(dataset_1), len(dataset_2))

    # convert to dataframe and eliminate Nans
    df = pd.DataFrame()
    df['dataset_1'] = list(dataset_1['Score'])[0:minlength]
    df['dataset_2'] = list(dataset_2['Score'])[0:minlength]
    df['Position'] = list(dataset_1['Position'])[0:minlength]
    df['d1 - d2'] = df['dataset_1'] - df['dataset_2']
    df.dropna(how='any', inplace=True)
    return df


# ## Rank 

# In[ ]:


def plot_rank(self, mode='pointmutant', outdf=False, **kwargs):
    '''
    Generate a rank plot so every mutation/residue is sorted based on enrichment score.

    Parameters
    ----------
    self : object from class "Screen"
        
    mode : str, default 'pointmutant'. 
        Alternative set to "mean" for the mean of each position
    
    outdf : boolean, default False
        If set to true, will return the df with the rank of mutations
        
    **kwargs : other keyword arguments

    Returns
    ----------
    Pandas dataframe
    '''

    # update kwargs
    temp_kwargs = copy.deepcopy(default_kwargs)
    temp_kwargs.update(kwargs)
    temp_kwargs['figsize'] = kwargs.get('figsize', (4, 2))
    temp_kwargs['x_label'] = kwargs.get('x_label', 'Rank')
    temp_kwargs['y_label'] = kwargs.get('y_label', r'$∆E^i_x$')
    
    # Sort by enrichment scores
    df = self.dataframe.sort_values(by=['Score']).copy()
    
    # Chose mode:
    if mode == 'mean':
        df = df.groupby(by=['Position'],as_index=False).mean()
        df.sort_values(by=['Score'], inplace=True)

    # create figure
    fig, ax = plt.subplots(figsize=temp_kwargs['figsize'])

    # import parameters
    _parameters()

    # Scatter data points
    plt.scatter(np.arange(len(df),0,-1), df['Score'], c='k', s=1)

    # Titles
    plt.title(temp_kwargs['title'], fontsize=12,
              fontname='Arial', color='k', pad=8)
    # Labels
    plt.ylabel(temp_kwargs['y_label'], fontsize=10,
               fontname="Arial", color='k', labelpad=0)
    plt.xlabel(temp_kwargs['x_label'], fontsize=10,
               fontname="Arial", color='k')
    
    # other graph parameters
    plt.xlim(temp_kwargs['xscale'])
    plt.ylim(temp_kwargs['yscale'])
    
    # save file
    _savefile(fig, temp_kwargs)
    
    if temp_kwargs['show']:
        plt.show()
    if outdf:
        return df


# ## SNV

# ### Plot Histogram

# In[ ]:


def plot_hist(self, population='All', loc='upper left', **kwargs):
    '''
    Generate a histogram plot. Can plot single nucleotide variants (SNVs) or non-SNVs only

    Parameters
    ----------
    population : str, default 'All'. 
        Other options are 'SNV' and 'nonSNV'.
    
    loc : str, default 'upper left'. 
        Position of the legend.
    
    **kwargs : other keyword arguments
        bins : int, default 50. 
            Number of bins for the histogram.
    Returns
    ----------
    None.

    '''
    # update kwargs
    temp_kwargs = copy.deepcopy(default_kwargs)
    temp_kwargs.update(kwargs)
    temp_kwargs['figsize'] = kwargs.get('figsize', (2, 2))
    temp_kwargs['yscale'] = kwargs.get('yscale', (0, 2))
    temp_kwargs['xscale'] = kwargs.get('xscale', (-2, 2))

    # Select case input data
    df = self.dataframe['Score_NaN']
    if population == 'SNV':
        df = self.dataframe_SNV['Score_NaN']
    elif population == 'nonSNV':
        df = self.dataframe_nonSNV['Score_NaN']

    # create figure
    fig = plt.figure(figsize=temp_kwargs['figsize'])

    # Import parameters
    _parameters()

    # plot figure
    plt.hist(df, density=True, bins=temp_kwargs['bins'], color='k')

    # axes labels and title
    plt.xlabel(r'$∆E^i_x$' if temp_kwargs['x_label'] == 'x_label' else temp_kwargs['x_label'],
               fontsize=10, fontname="Arial", color='k', labelpad=0)
    plt.ylabel('Probability density', fontsize=10,
               fontname="Arial", color='k', labelpad=3)
    plt.title(temp_kwargs['title'], fontsize=10, fontname='Arial', color='k')

    # axes limits. spacer will be 1 or the
    plt.xlim(temp_kwargs['xscale'])
    plt.xticks(np.arange(temp_kwargs['xscale'][0], temp_kwargs['xscale']
                         [1]+temp_kwargs['tick_spacing'], temp_kwargs['tick_spacing']))
    plt.ylim(temp_kwargs['yscale'])
    plt.grid()

    # save file
    _savefile(fig, temp_kwargs)
    if temp_kwargs['show']:
        plt.show()
    return


# ### Internal SNV

# In[ ]:


def _select_nonSNV(df):
    '''
    Generate a dataframe that contains the non-SNV variants and the enrichment score

    Parameters
    -----------
    df : pd.dataframe

    Returns
    --------
    Dataframe containing a column of variants that are non-SNV, and the Score.
    '''
    # Dataframe with SNV
    SNV = _select_SNV(df)

    # Merge and eliminate duplicates. Keep Non-SNV
    NonSNV = pd.concat([SNV, df], sort=False)[['Position', 'Variant', 'Score']]
    NonSNV.drop_duplicates(subset='Variant', keep=False, inplace=True)

    return NonSNV


def _select_SNV(df):
    '''
    Select for SNV variants in DSM dataset

    Parameters
    -----------
    df : pandas dataframe containing DSM data

    Returns
    --------
    Modified dataframe('Variant','Score') where 'SNV?'== True. Returns copy
    '''

    # Use _add_SNV_boolean funciton
    df = _add_SNV_boolean(df.copy())

    # Select SNV? == True only
    df = df[df['SNV?'] == True].copy()

    # Select columns of interest
    df = df[['Position', 'Variant', 'Score']].copy()

    # Reset index
    df.reset_index(drop=True, inplace=True)

    return df


def _aminoacids_snv(aa1, aa2, codontable):
    '''
    Determine if two amino acids are snv (one base difference)

    Parameters
    -----------
    aa1 : str
    aa2 : str
    codontable : dict (did not want to generate each time I run the function)

    Returns
    --------
    boolean, True/False
    '''
    # Convert amino acids to codons
    codons1 = codontable[aa1]
    codons2 = codontable[aa2]

    # Generate a list of combination pairs between all codons in aa1 and aa2
    codon_combinations = list(itertools.product(codons1, codons2))

    # If one pair of combinations is a SNV, then return True
    for combination in codon_combinations:
        if _codons_pointmutants(combination[0], combination[1]) == True:
            return True
    return False


def _add_SNV_boolean(df):
    '''
    Add a column to dataframe indication if the variant is a SNV or not

    Parameters
    -----------
    df : pandas dataframe containing DSM data

    Returns
    --------
    Modified dataframe. Returns copy
    '''

    # Generate dictionary with aa and codon translation
    codontable = _dict_codontoaa()

    # Add column with True/False input
    df['SNV?'] = df.apply(lambda x: _aminoacids_snv(
        x['Sequence'], x['Aminoacid'], codontable), axis=1)

    return df


def _codons_pointmutants(codon1, codon2):
    '''
    Determine if two codons are SNV. Returns a boolean

    Parameters
    -----------
    codon1 : str
    codon2 : str

    Returns
    --------
    boolean, True/False
    '''
    counter_occurrences = 0
    for index, base1 in enumerate(codon1):
        base2 = list(codon2)[index]
        if base1 == base2:
            counter_occurrences = counter_occurrences+1
    if counter_occurrences > 1:
        return True
    return False


def _are_pointmutants(aa, seqbase):
    '''
    converts the amino acid to all possible degenerate codons and then checks if they are point mutants

    Parameters
    -----------
    aa: str
    seqbase: str    

    Returns 
    --------
    Boolean
    '''
    codontoaadict = _dict_codontoaa()
    pointmutants = False
    for codon in _codontoaadict[aa]:
        if _codons_pointmutants(seqbase, codon):
            pointmutants = True
    return pointmutants


def _are_pointmutants_list(aa, seqbase_list):
    '''
    converts the amino acid to all possible degenerate codons and then checks if they are point mutants
    Same as _are_pointmutants but in list format

    Parameters
    -----------
    aa: str
    seqbase_list: list of str    

    Returns 
    --------
    List of Boolean
    '''
    pointmutants_list = []

    for seqbase in seqbase_list:
        pointmutants_list.append(_are_pointmutants(aa, seqbase))
    return pointmutants_list


def _dict_codontoaa():
    '''
    Generates a dictionary with all amino acids and all possible codons.
    aa is the aminoacid of the mutation and seqbase is the original codon of the wtsequence
    '''
    bases = ['T', 'C', 'A', 'G']
    codons = [a+b+c for a in bases for b in bases for c in bases]
    aminoacids = list(
        'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG')

    # dictionary with more than one value for each key
    codontoaadict = defaultdict(list)
    for codon, aminoacid in zip(codons, aminoacids):
        codontoaadict[aminoacid].append(codon)
    return codontoaadict


def _aatocodons(aminoacid):
    '''
    Inputs an aminoacid, returns all codons. Used dict_codontoaa()

    Parameters
    -----------
    aminoacid : str

    Returns
    --------
    List with all the codons that code for that amino acid
    '''

    # Dictionary with all codons and aa
    codontoaadict = _dict_codontoaa()

    # Codons for that amino acid
    codons = codontoaadict[aminoacid]

    return codons


def _aatocodons_df(df, namecolumn):
    '''
    Inputs a dataframe with a column of amino acids, returns all syn for each amino acidcodons. 
    Used dict_codontoaa() and _aatocodons

    Parameters
    -----------
    df : pandas dataframe
    namecolumn : str
        name of the column containing the amino acids

    Returns
    --------
    dataframe with a column containing all the codons that code for that amino acid. Returns copy
    '''
    # Copy df
    df = df.copy()

    # Calculate each possible codon for every amino acid
    df['Codons_' +
        namecolumn] = df.apply(lambda x: _aatocodons(x[namecolumn]), axis=1)

    return df


# ## Miniheatmap

# ### Mean substitution heatmap

# In[ ]:


def plot_miniheatmap(self, offset=0, **kwargs):
    '''
    Generate a miniheatmap plot enrichment scores of mutagenesis selection assays.

    Parameters
    ----------
    self : object from class "Screen"
    
    offset : int, default 0
        if you want to study effects of a residue when is behind or in front of another residue.
        offset of 1 means that you evaluate the effect of following residue n+1 on n.
        offset of -1 means that you look at the previous residue (n-1 on n).
    
    **kwargs : other keyword arguments
    
    Returns
    ----------
    None.
    '''

    # load font parameters
    _font_parameters()

    # update kwargs
    temp_kwargs = copy.deepcopy(default_kwargs)
    temp_kwargs.update(kwargs)

    # load labels
    temp_kwargs['color_sequencelabels'] = _labels(self.start_position)[0]
    temp_kwargs['number_sequencelabels'] = _labels(self.start_position)[1]

    # do offset if appropriate
    dataframe_stopcodons = _transform_dataset_offset(self, offset)

    # calculate condensed heatmap
    dataset = _condense_heatmap(
        dataframe_stopcodons, temp_kwargs['neworder_aminoacids'])

    _plot_miniheatmap(dataset, temp_kwargs)

    return


def _condense_heatmap(df, new_order):
    '''
    Converts the np.array with stored enrichment scores into the condensed heatmap
    '''
    # Convert dataset to df
    df = df.copy()
    df.drop(['Position'], axis=1, inplace=True)

    # Group by sequence and aminoacid, and then pivot table
    df_grouped = df.groupby(['Sequence', 'Aminoacid'], sort = False).mean()
    df_pivoted = df_grouped.pivot_table(values='Score',
                                        index='Aminoacid',  columns='Sequence')
    df_pivoted.reset_index(drop=False, inplace=True)

    # Sort in y axis desired order
    df_pivoted['Aminoacid'] = pd.Categorical(
        df_pivoted['Aminoacid'], new_order)
    df_pivoted = df_pivoted.sort_values(by=['Aminoacid'])

    # Sort in x axis desired order
    x_order = _common(new_order, list(df_pivoted.columns))

    # Drop amino acid column
    data_dropped = df_pivoted.drop(['Aminoacid'], axis=1)

    return data_dropped[x_order]


def _offset_sequence(dataset, sequence, start_position, offset):
    '''
    Internal function that offsets the input sequence

    Parameters
    -----------
    dataset, sequence, start_position, offset

    Returns
    --------
    string containing trimmed sequence
    '''
    # Deep copy sequence
    sequence = copy.deepcopy(sequence)

    # truncate sequence
    if offset > 0:
        sequence = sequence+'X'*np.absolute(offset)
        trimmedsequence = sequence[start_position-1 +
                                   offset:len(dataset[0])+start_position-1+offset]
    else:
        sequence = 'X'*(np.absolute(offset))+sequence
        trimmedsequence = sequence[start_position -
                                   1:len(dataset[0])+start_position-1]

    return trimmedsequence


def _transform_dataset_offset(self, offset, stopcodons=True):
    '''
    Generate a dataframe with the sequence offset. Reutilizes _transform_dataset
    '''
    # Add offset sequence
    offset_sequence = _offset_sequence(self.dataset, self.sequence_raw,
                                       self.start_position, offset)
    df = self.dataframe_stopcodons.copy() if stopcodons is True else self.dataframe.copy()

    # Copy old sequence
    df['Sequence_old'] = df['Sequence']
    # Count amino acids
    aa_number = len(set(df['Aminoacid']))
    # Generate new offset sequence
    df['Sequence'] = np.ravel([[aa]*aa_number for aa in offset_sequence])

    # Drop rows with X
    df.drop(df.index[df['Sequence'] == 'X'], inplace=True)

    return df


# ### Neighbor residues

# In[ ]:


def plot_neighboreffect (self, offset=1, **kwargs):
   '''
   Generate a miniheatmap plot telling you the effect of having a residue in front or behind.
   It corrects for the effect of that amino acid on the rest of the population.

   Parameters
   ----------
   self : object from class "Screen"
   
   offset : int, default 1
       if you want to study effects of a residue when is behind or in front of another residue.
       offset of 1 means that you evaluate the effect of following residue n+1 on n. On a "MTEY..." sequence,
       you would look at the effect of T on M, E on T, Y on E, etc.. and then group by residue (n+1).
       offset of -1 means that you look at the previous residue (n-1 on n).
   
   **kwargs : other keyword arguments

   Returns
   ----------
   None.
   '''
   # load font parameters
   _font_parameters()

   # update kwargs
   temp_kwargs = copy.deepcopy(default_kwargs)
   temp_kwargs.update(kwargs)
   if '*' in temp_kwargs['neworder_aminoacids']: temp_kwargs['neworder_aminoacids'].remove('*')
   
   # do offset, no stop codons
   df = _normalize_neighboreffect(self,offset,temp_kwargs['neworder_aminoacids'])
   
   # Plot
   _plot_miniheatmap(df,temp_kwargs)
   
   return
   
def _plot_miniheatmap(df,temp_kwargs):
   # declare figure and subplots
   coeff = len(df.columns)/19*1.05
   fig = plt.figure(figsize=(2.5*coeff, 2.5))
   gs = gridspec.GridSpec(nrows=1, ncols=1)
   ax = plt.subplot(gs[0])

   # main heatmap
   heatmap = ax.pcolor(df.to_numpy(), vmin=temp_kwargs['colorbar_scale'][0], vmax=temp_kwargs['colorbar_scale'][1],
                       cmap=temp_kwargs['colormap'], edgecolors='k', linewidths=0.2, color='darkgrey')

   # ____________axes manipulation____________________________________________
   # put the major ticks at the middle of each cell
   ax.set_xticks(np.arange(df.shape[1]) + 0.5, minor=False)
   ax.set_yticks(np.arange(df.shape[0]) + 0.5, minor=False)

   # position of axis labels
   ax.tick_params('x', direction='out', pad=-2.5)
   ax.tick_params('y', direction='out', pad=0.4)

   # want a more natural, table-like display
   ax.invert_yaxis()
   ax.xaxis.tick_top()

   # remove ticks
   ax.xaxis.set_ticks_position('none')
   ax.yaxis.set_ticks_position('none')

   # so labels of x and y do not show up and my labels show up instead
   ax.set_xticklabels(list(df.columns), fontsize=6.5,
                      fontname="Arial", color='k', minor=False)
   ax.set_yticklabels(temp_kwargs['neworder_aminoacids'],
                      fontsize=6.5, fontname="Arial", color='k', minor=False)

   # align the labels of the y axis
   for ylabel in ax.get_yticklabels():
       ylabel.set_horizontalalignment('center')

   # _____________________________________________________________________________

   # for color bar format
   cb = plt.colorbar(heatmap, fraction=0.025, pad=0.05, aspect=5, ticks=[temp_kwargs['colorbar_scale'][0], np.mean(
       temp_kwargs['colorbar_scale']), temp_kwargs['colorbar_scale'][1]], orientation='vertical')
   cb.ax.set_yticklabels(cb.ax.get_yticklabels(), fontsize=7, fontname="Arial", color='k')
   cb.update_ticks()
   plt.text(len(df.columns)+2, 7.8, r'$\langle∆E^x_i\rangle_x$', horizontalalignment='center',
            fontsize=7, fontname="Arial", color='k')

   # for putting title on graph
   plt.title(temp_kwargs['title'], horizontalalignment='center',
             fontname="Arial", fontsize=10, pad=10)
   plt.ylabel('Amino Acid Substitution', fontsize=10, labelpad=-1)

   # save file
   _savefile(fig,temp_kwargs)
   
   if temp_kwargs['show']: plt.show()
   return

def _normalize_neighboreffect(self,offset,neworder):
   '''
   For every residue, subtract the average effect of a substitution
   Returns a normalized dataframe
   '''
   aalist = list('ACDEFGHIKLMNPQRSTVWY')
   # Add offset sequence to df
   df = _transform_dataset_offset(self,offset,False)
   
   # calculate mean effect using condensed heatmap
   mean = _condense_heatmap(self.dataframe, aalist)
   
   df_normalized = pd.DataFrame()
   for aa in aalist:
       # Choose the neighbors of an aa
       aa_neighbors = df.loc[df['Sequence']==aa]
       # Do the mean substitution of amino acids that are repeated
       aa_neighbors = aa_neighbors.groupby(['Sequence_old','Aminoacid'],as_index=False).mean()
       # Make into table
       aa_neighbors_pivoted = aa_neighbors.pivot_table(values='Score', index='Aminoacid',  columns='Sequence_old')
       aa_neighbors_pivoted.reset_index(drop=True, inplace=True)
       # Get the mean of the amino acids that appear in the aa_neighbors subset
       mean_neighbors = mean[list(aa_neighbors_pivoted.columns)]
       # Subtract average effect and do mean
       df_normalized[aa] = (aa_neighbors_pivoted - mean_neighbors).mean(axis=1)
   
   # Sort by aa
   df_normalized = df_normalized[neworder]
   # Sort in y axis desired order
   df_normalized = _sort_yaxis_aminoacids(df_normalized,neworder,aalist)
   return df_normalized

def _sort_yaxis_aminoacids(df,neworder,oldorder=list('ACDEFGHIKLMNPQRSTVWY')):
   # Sort in y axis desired order
   df['Aminoacid_new'] = oldorder
   df['Aminoacid_new'] = pd.Categorical(df['Aminoacid_new'], neworder)
   df.sort_values(by=['Aminoacid_new'],inplace=True)
   df.drop(['Aminoacid_new'], inplace=True, axis=1)
   
   return df


# ## Correlation

# ### Heatmap correlation

# In[ ]:


def plot_correlation(self, **kwargs):
    '''
    Generate a correlation of each amino acid

    Parameters
    ----------
    **kwargs : other keyword arguments

    Returns
    ----------
    None.
    '''

    # load font parameters
    _font_parameters()

    # update kwargs
    temp_kwargs = copy.deepcopy(default_kwargs)
    temp_kwargs.update(kwargs)

    # load labels
    temp_kwargs['color_sequencelabels'] = _labels(self.start_position)[0]
    temp_kwargs['number_sequencelabels'] = _labels(self.start_position)[1]

    # calculate correlation heatmap
    dataset = _calculate_correlation(
        self.dataframe_stopcodons, temp_kwargs['neworder_aminoacids'])

    # declare figure and subplots
    coeff = len(dataset.columns)/19*1.05
    fig = plt.figure(figsize=(2.5*coeff, 2.5))
    gs = gridspec.GridSpec(nrows=1, ncols=1)
    ax = plt.subplot(gs[0])

    # main heatmap
    heatmap = ax.pcolor(dataset.corr(), vmin=temp_kwargs['colorbar_scale'][0], vmax=temp_kwargs['colorbar_scale'][1],
                        cmap='Greys', edgecolors='k', linewidths=0.2, color='darkgrey')

    # ____________axes manipulation____________________________________________
    # put the major ticks at the middle of each cell
    ax.set_xticks(np.arange(dataset.shape[1]) + 0.5, minor=False)
    ax.set_yticks(np.arange(dataset.shape[0]) + 0.5, minor=False)

    # position of axis labels
    ax.tick_params('x', direction='out', pad=-2.5)
    ax.tick_params('y', direction='out', pad=0.4)

    # want a more natural, table-like display
    ax.invert_yaxis()
    ax.xaxis.tick_top()

    # remove ticks
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')

    # so labels of x and y do not show up and my labels show up instead
    ax.set_xticklabels(list(dataset.columns), fontsize=6.5,
                       fontname="Arial", color='k', minor=False)
    ax.set_yticklabels(temp_kwargs['neworder_aminoacids'],
                       fontsize=6.5, fontname="Arial", color='k', minor=False)

    # align the labels of the y axis
    for ylabel in ax.get_yticklabels():
        ylabel.set_horizontalalignment('center')

    # _____________________________________________________________________________

    # for color bar format
    cb = plt.colorbar(heatmap, fraction=0.025, pad=0.05, aspect=5, ticks=[temp_kwargs['colorbar_scale'][0], temp_kwargs['colorbar_scale'][1]],
                      orientation='vertical')
    cb.ax.set_yticklabels(cb.ax.get_yticklabels(),
                          fontsize=7, fontname="Arial", color='k')
    cb.update_ticks()
    plt.text(len(dataset.columns)+1.2*coeff, len(dataset.columns)/2.5, 'R',
             horizontalalignment='center', fontsize=7, fontname="Arial", color='k')

    # for putting title on graph
    plt.title(temp_kwargs['title'], horizontalalignment='center',
              fontname="Arial", fontsize=10, pad=10)

    # save file
    _savefile(fig, temp_kwargs)

    if temp_kwargs['show']:
        plt.show()
    return


def _calculate_correlation(df, order_aminoacids):

    dataset = df.copy()
    dataset = dataset.pivot_table(
        values='Score', index='Position',  columns='Aminoacid')
    dataset = dataset.corr()
    dataset = dataset.reindex(index=order_aminoacids)[order_aminoacids]

    return dataset


def _calculate_correlation_byresidue(df):

    dataset = df.copy()
    dataset = dataset.pivot_table(
        values='Score', index='Position',  columns='Aminoacid')
    dataset = dataset.T.corr()

    return dataset


# ### Individual correlation

# In[ ]:


def plot_individual_correlation(self, **kwargs):
    '''
    Genereates a bar plot of the correlation of each amino acid mutational 
    profile (row of the heatmap) with the rest of amino acids (rows)
    
    Parameters
    -----------
    self : object from class "Screen"
    
    **kwargs : other keyword arguments

    Returns
    ----------
    None.

    '''
    # Load parameters
    _parameters()

    # Update kwargs
    temp_kwargs = copy.deepcopy(default_kwargs)
    temp_kwargs.update(kwargs)
    temp_kwargs['figsize'] = kwargs.get('figsize', (3.5, 2))
    temp_kwargs['yscale'] = kwargs.get('yscale', (0, 1))

    # Get data
    if '*' in temp_kwargs['neworder_aminoacids']:
        temp_kwargs['neworder_aminoacids'].remove('*')
    df = _calculate_correlation(
        self.dataframe, temp_kwargs['neworder_aminoacids']).mean()**2

    # Make figure
    fig, ax = plt.subplots(figsize=temp_kwargs['figsize'])
    ticks = np.arange(0, len(df))  # label locations
    width = 0.5
    labels = temp_kwargs['neworder_aminoacids']
    # Plot figure
    ax.bar(ticks, df, width, color='blue', ec='k',)

    # graph parameters
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels, fontsize=9, fontname="Arial",
                       color='k', minor=False, rotation=0)
    ax.set_ylabel(r'$R^2$', fontsize=10, fontname="Arial",
                  color='k', labelpad=12, rotation=0)
    ax.set_ylim(temp_kwargs['yscale'])
    plt.title(temp_kwargs['title'], horizontalalignment='center',
              fontname="Arial", fontsize=10, pad=5)

    # save file
    _savefile(fig, temp_kwargs)
    if temp_kwargs['show']:
        plt.show()
    return


# ### Group correlation

# In[ ]:


def plot_group_correlation(self, r2, groups=['DEHKR', 'QN', 'CASTG', 'ILMV', 'WYF'],
                        output=False, **kwargs):
    '''
    Determines which amino acids better represent the heatmap. Requires logomaker package.

    Parameters
    -----------
    self : object from class "Screen"
    
    r2 : float
        cutoff of the r**2 correlation value. Only values above that will be plot at the sequence logo
    
    groups : list, default ['DEHKR','QN','CASTG','ILMV','WYF']
        groups of aa to combine together
    
    output : boolean, default False
    
    **kwargs : other keyword arguments

    Returns
    --------
    Use logomaker to plot the most frequent residues. 
    Optional gives back the different combinations of groups and the R**2 values
    '''

    # Update kwargs
    temp_kwargs = copy.deepcopy(default_kwargs)
    temp_kwargs.update(kwargs)

    # Apply parameters
    _parameters()

    # If there is a stop codon, delete it
    if '*' in temp_kwargs['neworder_aminoacids']:
        temp_kwargs['neworder_aminoacids'].remove('*')
        
    # Get R2 of each combination of amino acid substitutions
    df = _calculate_substitution_correlations(self, temp_kwargs['neworder_aminoacids'], groups)

    # Filter according the the R2 correlation value
    filtered = df.loc[df['R2'] > r2]
    logoplot = logomaker.alignment_to_matrix(list(filtered['Combinations']))

    # create Logo object
    fig = logomaker.Logo(logoplot, font_name='Arial', color_scheme='chemistry', vpad=.1,
                         width=.8, figsize=((len(logoplot)+1)/2.5, 1))

    # style using Logo methods
    fig.style_xticks(anchor=0, spacing=1, rotation=0)

    # No yticks and no xticks (but keep labels)
    plt.yticks([], [])
    fig.ax.tick_params(axis='both', which='both', length=0)

    # style using Axes methods
    fig.ax.set_ylabel('Bits')
    fig.ax.set_xlim([-0.5, len(logoplot)-0.5])

    # for putting title on graph
    plt.title(temp_kwargs['title'], horizontalalignment='center',
              fontname="Arial", fontsize=10, pad=10)

    # save file, cannot save logo file for now
    #_savefile(fig.ax, temp_kwargs)

    if temp_kwargs['show']:
        plt.show()
        
    if output:
        return df


def _calculate_substitution_correlations(self, aminoacids, groups):
    '''if a set of residues was chosen, how well would they represent the entire population'''
    
    # Get correlation values
    corr_values = _calculate_correlation(self.dataframe, aminoacids)**2
    corr_values.reset_index(inplace=True)

    # Get combinations
    replacement_combinations = list(itertools.product(*groups))

    # Retrieve Correlation values
    df = pd.DataFrame()
    df['Aminoacids'] = list(itertools.chain.from_iterable(groups))
    for combination in replacement_combinations:  # Iterate over a combination
        temp_list = []
        
        # Iterate over a group of the combination
        for group, aa_selected in zip(groups, combination):
            for aa_nonselected in group:  # Find correlation values from correlation plot
                if aa_nonselected == aa_selected:
                    temp_list.append(1)
                else:
                    temp_list.append(_find_correlation(
                        aa_selected, aa_nonselected, corr_values))
        df[combination] = temp_list  # Store in df
    return _polishdf(df)


def _polishdf(df):
    df_mean = df.copy()
    df_mean = df.mean().to_frame()
    df_mean.reset_index(drop=False, inplace=True)
    df_mean.rename(columns={0: 'R2'}, inplace=True)
    df_mean['Combinations'] = list(df_mean['index'].apply(lambda x: ''.join(x)))
    df_mean.drop(columns=['index'], inplace=True)
    return df_mean


def _find_correlation(aa1, aa2, corr_values):
    return float(corr_values[aa1].loc[corr_values['Aminoacid'] == aa2])


# ## PCA

# In[ ]:


def plot_pca(self, mode='aminoacid', dimensions=[0, 1], adjustlabels = False, **kwargs):
    '''
    Genereates a plot of two PCA dimensions

    Parameters
    -----------
    self : object from class "Screen"
    
    mode : list, default 'aminoacid'
        Can also do PCA by secondary structure element if set to "secondary" or 
        by individual residue if set to "individual".
    
    dimensions : list, default [0,1]
        Specify which two PCA dimensions to plot. By default PCA1 vs PCA2.
        Max dimension is 5.
    
    adjustlabels : boolean, default False
        If set to true, it will adjust the text labels so there is no overlap. It is convenient to increase
        the size of the figure, otherwise the algorithm will not find a solution. Requires to install adjustText package.
    
    **kwargs : other keyword arguments
        random_state : int, default 554
    Returns
    ----------
    None.

    '''

    # load parameters
    _parameters()

    # update kwargs
    temp_kwargs = copy.deepcopy(default_kwargs)
    temp_kwargs.update(kwargs)
    temp_kwargs['figsize'] = kwargs.get('figsize', (2, 2))

    # calculate correlation heatmap. Choose mode
    dataset = self.dataframe.copy()
    if mode == 'aminoacid':
        if '*' in temp_kwargs['neworder_aminoacids']:
            temp_kwargs['neworder_aminoacids'].remove('*')
        dataset = _calculate_correlation(
            dataset, temp_kwargs['neworder_aminoacids'])
        textlabels = temp_kwargs['neworder_aminoacids']
    elif mode == 'secondary':
        dataset = _calculate_correlation_bysecondary(
            dataset, self.secondary_dup)
        textlabels = list(dataset.columns)
    elif mode == 'individual':
        dataset = _calculate_correlation_byresidue(dataset)
        textlabels = list(dataset.columns)
    
    # plot using plot_clusters
    dimensionstoplot, variance = _calculate_clusters(dataset, dimensions, temp_kwargs['random_state'])

    # x and y
    x = dimensionstoplot.iloc[:, 0]
    y = dimensionstoplot.iloc[:, 1]

    fig, ax = plt.subplots(figsize=temp_kwargs['figsize'])
    ax.scatter(x, y, s=4, c='k')

    # labels
    plt.xlabel('PCA ' + str(dimensions[0]+1) + ': ' + str(int(variance[dimensions[0]]*100))+'%',
               fontsize=10, labelpad=5, fontweight='normal')
    plt.ylabel('PCA ' + str(dimensions[1]+1) + ': ' + str(int(variance[dimensions[1]]*100))+'%',
               fontsize=10, labelpad=-2, fontweight='normal')

    # label of data points
    texts = _auto_text(x, y, textlabels)
    if adjustlabels is True:
        adjust_text(texts, autoalign='xy')
    
    # set title
    plt.title(temp_kwargs['title'], horizontalalignment='center',
              fontname="Arial", fontsize=10, pad=5)

    # save file
    _savefile(fig, temp_kwargs)

    if temp_kwargs['show']:
        plt.show()
    return


def _auto_text(x, y, textlabels):
    '''auto anotates text labels'''
    texts = [plt.annotate(textlabels[i],  # this is the text
                          (x[i], y[i]),  # this is the point to label
                          textcoords="offset points",  # how to position the text
                          xytext=(2, 2),  # distance from text to points (x,y)
                          fontsize=8,
                          ha='center')  # horizontal alignment can be left, right or center
             for i in range(len(textlabels))]
    return texts


def _calculate_clusters(dataset, dimensions, random_state):
    '''input the dataframe that needs to be correlated, the dimensions, and will calculate PCA descomposition. '''

    # call pca model
    pca = PCA(n_components=6, random_state = random_state)

    # fit model to df. use aux function correlation_aminoacids
    model = pca.fit(dataset)

    # create df with PCA data
    df_aa = pd.DataFrame((model.components_).T, columns=[
                         'PCA1', 'PCA2', 'PCA3', 'PCA4', 'PCA5', 'PCA6'])

    # use kmeans to cluster the two dimensions and color
    dimensionstoplot = df_aa.iloc[:, np.r_[dimensions[0], dimensions[1]]]

    return dimensionstoplot, pca.explained_variance_ratio_


def _grouby_secondary(df, secondary):
    '''
    Groups each secondary motif and makes the mean.

    Returns dataframe. Returns copy
    '''
    df = df.copy()
    df.insert(4, 'Secondary', secondary)
    df = df.groupby(['Secondary', 'Aminoacid'], as_index=False).mean()
    df = df.loc[df['Secondary'].str.startswith(('β', 'α'))]
    return df


def _calculate_correlation_bysecondary(df, secondary):
    dataset = _grouby_secondary(df, secondary)
    dataset = dataset.pivot_table(
        values='Score', index='Secondary',  columns='Aminoacid')
    dataset = dataset.T.corr()

    return dataset


# ## Secondary Structure

# In[ ]:


def plot_secondary(self, **kwargs):
    '''
    Genereates a bar plot of data sorted by secondary elements (alpha helices and beta sheets).

    Parameters
    -----------
    self : object from class "Screen"
    
    **kwargs : other keyword arguments

    Returns
    ----------
    None.

    '''
    # Load parameters
    _parameters()

    # Update kwargs
    temp_kwargs = copy.deepcopy(default_kwargs)
    temp_kwargs.update(kwargs)
    temp_kwargs['figsize'] = kwargs.get('figsize', (2.5, 2))
    temp_kwargs['yscale'] = kwargs.get('yscale', (-2, 1))

    # Get data
    df = _calculate_secondary(self.dataframe, self.secondary_dup)

    # Color
    df['Color'] = df.apply(color_data, axis=1)

    # Make figure
    fig, ax = plt.subplots(figsize=temp_kwargs['figsize'])
    ticks = np.arange(0, len(df))  # label locations
    width = 0.5
    labels = df['Secondary']
    
    # Plot figure
    ax.bar(ticks, df['Score'], width, color=df['Color'], ec='k',)

    # graph parameters
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels, fontsize=9, fontname="Arial",
                       color='k', minor=False, rotation=0)
    ax.set_ylabel(r'$∆E^i_x$', fontsize=10, fontname="Arial",
                  color='k', labelpad=12, rotation=0)
    ax.set_ylim(temp_kwargs['yscale'])
    plt.title(temp_kwargs['title'], horizontalalignment='center',
              fontname="Arial", fontsize=10, pad=5)

    # save file
    _savefile(fig, temp_kwargs)
    
    if temp_kwargs['show']:
        plt.show()
    return


def _calculate_secondary(df, secondary):
    '''
    Returns copy
    '''
    df = df.copy()
    df.insert(4, 'Secondary', secondary)
    df = df.groupby(['Secondary'], as_index=False, sort=False).mean()
    df = df[df['Secondary'].str.startswith(('β', 'α'))]
    df = df.drop(['Position'], axis=1)
    return df


# ## ROC AUC

# In[ ]:


def plot_roc(self, df_class=None, **kwargs):
    '''
    Generates ROC AUC plot. It compares enrichment scores to some labels that the user has specified.

    Parameters
    -----------
    self : object from class "Screen"
    
    df_class: Pandas dataframe
        A dataframe that contains a column of variants labeled 'Variant' with a column labeled 'Class'
        containing the true class of that mutation. The true class can also be an input when creating the object.

    Returns
    --------
    None.
    '''
    # Use default class
    if df_class is None:
        df_class = self.roc_df

    # Merge dataframe with classes
    df = _mergeclassvariants(df_class, self.dataframe)

    # Calculate ROC parameters
    fpr, tpr, auc, _ = _rocauc(df)

    # update kwargs
    temp_kwargs = copy.deepcopy(default_kwargs)
    temp_kwargs.update(kwargs)
    temp_kwargs['figsize'] = kwargs.get('figsize', (2.5, 2.5))

    # import parameters
    _parameters()

    # create figure
    fig, ax = plt.subplots(figsize=temp_kwargs['figsize'])
    lw = 2
    plt.plot(fpr, tpr, color='k', lw=lw, label='AUC = %0.2f' % auc)
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')

    # Graph limits
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    tick_spacing = 0.2
    ax.xaxis.set_major_locator(
        ticker.MultipleLocator(tick_spacing))  # Plt ticks
    ax.yaxis.set_major_locator(
        ticker.MultipleLocator(tick_spacing))  # Plt ticks

    # Axis labels
    plt.title(temp_kwargs['title'], fontsize=12, fontname='Arial', color='k')
    plt.ylabel('True Positive Rate', fontsize=12,
               fontname="Arial", color='k', labelpad=0)
    plt.xlabel('False Positive Rate', fontsize=12, fontname="Arial", color='k')

    # Legend
    plt.legend(loc='lower right', handlelength=0,
               handletextpad=0, frameon=False)

    # save file
    _savefile(fig, temp_kwargs)
    if temp_kwargs['show']:
        plt.show()
    return


def _rocauc(df):
    '''
    Calculate roc rates and auc.

    The input is a dataframe that contains [Variants,Class,Score]
    '''
    fpr, tpr, thresholds = metrics.roc_curve(
        df['Class'], df['Score'], drop_intermediate=True)
    auc = metrics.roc_auc_score(df['Class'], df['Score'])
    return fpr, tpr, auc, thresholds


def _mergeclassvariants(df_score, df_class):
    '''
    Merge the input dataframe containing the class (true score) for variants and the enrichment scores
    '''
    # Merge DMS with true score dataset
    df_merged = pd.merge(df_class, df_score, on=['Variant'], how='left')

    # Drop rows with Nan values
    df_merged.dropna(inplace=True)

    return df_merged


def _concattrueposneg(df_tp, df_tn, subset='Variant', keep='first'):
    '''
    Concat a df containing the true positive variants and the true negative variants

    Parameters
    -----------
    df_tp : Dataframe with the true positives
    df_tn : Dataframe with the true negatives
    subset : str, default Variant
    keep : {‘first’, ‘last’, False} 

    Returns
    --------
    None.
    '''
    # Concatenate tp and tn datasets
    df_true = pd.concat([df_tp, df_tn], sort=False)

    # Will keep a variant as true positive if found in both datasets (because could be a mistake in gnomAD)
    df_true.drop_duplicates(subset=subset, keep=keep, inplace=True)

    return df_true


# ## Cumulative

# In[ ]:


def plot_cumulative(self, mode='all', **kwargs):
    '''
    Generates a cumulative plot of the enrichment scores by position. 

    Parameters
    -----------
    self : object from class "Screen"
    
    mode : str, default 'all' 
        Options are 'all','SNV' and 'nonSNV'.
    
    **kwargs : other keyword arguments

    Returns
    --------
    None.
    '''

    # update kwargs
    temp_kwargs = copy.deepcopy(default_kwargs)
    temp_kwargs.update(kwargs)
    temp_kwargs['figsize'] = kwargs.get('figsize', (3, 2))
    temp_kwargs['tick_spacing'] = kwargs.get('tick_spacing', 20)

    # import parameters
    _parameters()

    # create figure
    fig, ax = plt.subplots(figsize=temp_kwargs['figsize'])

    # Get data filtered
    df = _filter(self, mode)
    cumsum = df.cumsum(skipna=False)['Score']
    plt.plot(df['Position'], cumsum/list(cumsum)[-1], color='red', lw=2)

    # y label
    y_label = 'Cumulative LoF'
    if list(cumsum)[-1] > 0:
        y_label = 'Cumulative GoF'

    # Graph limits
    plt.xlim(self.dataframe['Position'].min(),
             self.dataframe['Position'].max()+1)
    plt.ylim(0, 1.1)

    ax.xaxis.set_major_locator(ticker.MultipleLocator(
        temp_kwargs['tick_spacing']))  # Plt ticks

    # Axis labels
    plt.title(temp_kwargs['title'], fontsize=12, fontname='Arial', color='k')
    plt.ylabel(y_label, fontsize=12, fontname="Arial", color='k', labelpad=5)
    plt.xlabel('Position', fontsize=12,
               fontname="Arial", color='k', labelpad=0)

    # x=y line
    plt.plot([0, df['Position'].max()], [0, 1],
             color='silver', lw=2, linestyle='--')

    # save file
    _savefile(fig, temp_kwargs)
    if temp_kwargs['show']:
        plt.show()
    return


def _filter(self, mode):

    # Select all, SNV, nonSNV
    if mode == 'all':
        df = self.dataframe
    elif mode == 'SNV':
        df = self.dataframe_SNV
    elif mode == 'nonSNV':
        df = self.dataframe_nonSNV
    df = df.groupby(by='Position', as_index=False).mean()
    return df


# ## Box Plot

# In[ ]:


def plot_box(binned_x, y, **kwargs):
    '''
    Genereates a boxplot. Data needs to be binned prior before using this function. 

    Parameters
    -----------
    x, binned_y : arrays
        Contain the data is going to plot
        
    **kwargs : other keyword arguments

    Returns
    ----------
    None.

    '''
    # Load parameters
    _parameters()

    # Update kwargs
    temp_kwargs = copy.deepcopy(default_kwargs)
    temp_kwargs.update(kwargs)
    temp_kwargs['figsize'] = kwargs.get('figsize', (2.5, 2))

    # Make figure
    fig, ax = plt.subplots(figsize=temp_kwargs['figsize'])
    
    # Plot data
    ax = sns.boxplot(binned_x, y, color='white', fliersize=2)

    plt.setp(ax.artists, edgecolor='k', facecolor='w')
    plt.setp(ax.lines, color='k')

    # graph parameters
    plt.title(temp_kwargs['title'], fontsize=10,
              fontname='Arial', color='k', pad=8)
    plt.ylabel(temp_kwargs['y_label'], fontsize=10,
               fontname="Arial", color='k', labelpad=0)
    plt.xlabel(temp_kwargs['x_label'], fontsize=10,
               fontname="Arial", color='k')

    # axes limits
    plt.xlim(temp_kwargs['xscale'])
    plt.ylim(temp_kwargs['yscale'])
    plt.grid()

    # save file
    _savefile(fig, temp_kwargs)
    
    if temp_kwargs['show']:
        plt.show()

    return


# ## 3D plot

# ### 3D Scatter

# In[2]:


def plot_scatter_3D(self, mode='mean', pdb_path=None, df_coordinates=None,
                    df_color=None,  position_correction=0, chain='A',
                    squared=False, rotate=False, **kwargs):
    '''
    Generates a 3-D scatter plot of the x,y,z coordinates of the C-alpha atoms of the residues, 
    color coded by the enrichment scores. PDBs may have atoms missing, 
    you should fix the PDB before using this method. Use matplotlib for interactive plot.

    Parameters
    -----------
    self : object from class "Screen"
        **kwargs : other keyword arguments.

    mode : str, default 'mean'
        Specify what enrichment scores to use. If mode = 'mean', it will use the mean of 
        each position to classify the residues. If mode = 'A', it will use the Alanine substitution profile. 
        Can be used for each amino acid. Use the one-letter code and upper case.

    pdb : str, default None
        User should specify the path PDB chain.

    df_coordinates: pandas dataframe, default None
        If no pdb is included, the user must pass the 3-D coordinates of the residues to plot. 
        In here you have more flexibility and you can select other atoms besides the C-alpha.

    df_color : pandas dataframe, default None     
        The color of each residue can also be included. You must label that label column.

    position_correction : int, default 0
        If the pdb structure has a different numbering of positions than you dataset,
        you can correct for that. If your start_position = 2, but in the PDB that same residue
        is at position 20, position_correction needs to be set at 18.

    chain : str, default 'A'
        Chain of the PDB file to get the coordinates and SASA from.

    squared : booleand, False
        If this parameter is True, the algorithm will center the data, and plot the square value of the 
        distance.

    rotate : boolean, False
        If you are using an interactive matplotlib, set up rotate = True so the graph spins.

    **kwargs : other keyword arguments
        gof : int, default is 1
                 cutoff for determining gain of function mutations based on mutagenesis data.
        lof : int, default is -1
            cutoff for determining loss of function mutations based on mutagenesis data.

    Returns
    ---------
    None
    '''

    # Load parameters
    _parameters()

    # Update kwargs
    temp_kwargs = copy.deepcopy(default_kwargs)
    temp_kwargs.update(kwargs)

    # Get Scores and colors
    if df_color is None:
        df = _color_3D_scatter(self.dataframe, mode, temp_kwargs['lof'],
                               temp_kwargs['gof'])

    # If coordinates is not an input, get it from the pdb
    if df_coordinates is None:
        df_coordinates = _parse_pdbcoordinates(
            self, pdb_path, position_correction, chain)

    # Plot figure
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    if squared is False:
        ax.scatter(df_coordinates['x'], df_coordinates['y'],
                   df_coordinates['z'], c=df['Color'])
    else:
        ax.scatter(df_coordinates['x_cent'], df_coordinates['y_cent'],
                   df_coordinates['z_cent'], c=df['Color'])
        

    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')

    # save file
    _savefile(fig, temp_kwargs)

    # rotate if you are using interactive plot
    if rotate is True:
        for angle in range(0, 360):
            ax.view_init(30, angle)
            plt.draw()
            plt.pause(.001)

    if temp_kwargs['show']:
        plt.show()

    return


def _color_3D_scatter(df, mode, lof, gof):
    '''Color the data points by enrichment scores'''
    # Copy df
    df_grouped = df.copy()

    # Select grouping
    if mode == 'mean':
        df_grouped = df_grouped.groupby(['Position'], as_index=False).mean()
    else:
        df_grouped = df_grouped.loc[df_grouped['Aminoacid'] == mode]

    # Select colors
    df_grouped['Color'] = 'green'
    df_grouped.loc[df_grouped['Score'] < lof, 'Color'] = 'blue'
    df_grouped.loc[df_grouped['Score'] > gof, 'Color'] = 'red'

    return df_grouped


def centeroidnp(df):
    '''find center of x,y,z'''
    return df['x'].sum()/len(df['x']), df['y'].sum()/len(df['y']), df['z'].sum()/len(df['z'])


def _parse_pdbcoordinates(self, pdb_path, position_correction, chain, sasa=False):
    '''parse coordinate of CA atoms. Will also return the bfactor and SASA using freesasa.
    If PDB is missing atoms, it can handle it.'''
    
    # Get structure from PDB
    structure = PDBParser().get_structure('pdb', pdb_path)

    coordinates = []
    commands = []
    bfactors = []
    positions_worked =[] # positions present in pdb
    
    # Iterate over each CA atom and geet coordinates
    for i in np.arange(self.start_position+position_correction, self.end_position+position_correction):
        # first check if atom exists
        try:
            structure[0][chain][int(i)].has_id("CA")
            # Get atom from pdb and geet coordinates
            atom = list(structure[0][chain][int(i)]["CA"].get_vector())+[i]
            coordinates.append(atom)
            # Get SASA command for each residue and bfactor
            residue = "s{}, chain {} and resi {}".format(str(i), chain, str(i))
            commands.append(residue)
            bfactor = (structure[0][chain][int(i)]["CA"].get_bfactor())
            bfactors.append(bfactor)
            positions_worked.append(i)
        except:
            print ("residue {} not found".format(str(i)))
            coordinates.append([np.nan, np.nan, np.nan, i])
            
    # Convert to df
    df_coordinates = pd.DataFrame(
        columns=['x', 'y', 'z', 'Position'], data=coordinates)

    # Center data
    x, y, z = centeroidnp(df_coordinates)
    df_coordinates['x_cent'] = (df_coordinates['x']-x).abs()**2
    df_coordinates['y_cent'] = (df_coordinates['y']-y).abs()**2
    df_coordinates['z_cent'] = (df_coordinates['z']-z).abs()**2
    df_coordinates['Distance'] = df_coordinates['x_cent'] +         df_coordinates['y_cent']+df_coordinates['z_cent']

    # Add sasa values
    if sasa:
        # Get structure for SASA
        structure_sasa = freesasa.Structure(pdb_path)
        result = freesasa.calc(structure_sasa)
        # Calculate sasa
        sasa = freesasa.selectArea(commands, structure_sasa, result)
        df_sasa = pd.DataFrame(columns=['SASA'], data=sasa.values())
        df_sasa['B-factor'] = bfactors
        df_sasa['Position'] = positions_worked

        # Merge
        df_coordinates = df_coordinates.merge(df_sasa,how='outer', on='Position')

    return df_coordinates


# ### 3D Scatter Second version

# In[3]:


def plot_scatter_3D_pdbprop(self, plot=['Distance', 'SASA', 'B-factor'],
                            mode='mean', pdb_path=None, custom=None, 
                            axis_scale = ["linear", "linear", "linear"],
                            df_color=None,  color_by_score=True,
                            position_correction=0, chain='A',
                            rotate=False, output_df= False, **kwargs):
    '''
    Generates a 3-D scatter plot of different properties obtained from the PDB. 
    PDBs may have atoms missing, you should fix the PDB before using this
    method. We recommend you use matplotlib for interactive plot. 

    Parameters
    -----------
    self : object from class "Screen"
        **kwargs : other keyword arguments.

    plot : list, default ['Distance', 'SASA', 'B-factor']
        List of 3 elements to plot. Other options are 'Score' and Custom. If custom, add the 
        label to the third element of the list ie ['Distance', 'SASA', 'Conservation']. 

    mode : str, default 'mean'
        Specify what enrichment scores to use. If mode = 'mean', it will use the mean of 
        each position to classify the residues. If mode = 'A', it will use the Alanine substitution profile. Can be 
        used for each amino acid. Use the one-letter code and upper case.

    pdb_path : str, default None
        User should specify the path PDB.
    
    custom : list or dataframe or np.array, default None
        If you want to add a custom dataset to plot, use custom. On the parameter
        plot, the 3rd item of the list will be the label for your custom dataset.
    
    axis_scale : list, default ["linear", "linear", "linear"]
        Check matplotlib.axes.Axes.set_xscale documentation for more information.
        The axis scale type to apply. Some options are {"linear", "log", "symlog", "logit", ...}.
        
    df_color : pandas dataframe, default None     
        The color of each residue can also be included. You must label that label column.
    
    color_by_score : boolean, default True
        If set to False, the points in the scatter will not be colored based on the enrichment score.

    position_correction : int, default 0
        If the pdb structure has a different numbering of positions than you dataset,
        you can correct for that. If your start_position = 2, but in the PDB that same residue
        is at position 20, position_correction needs to be set at 18.

    chain : str, default 'A'
        Chain of the PDB file to get the coordinates and SASA from.

    rotate : boolean, False
        If you are using an interactive matplotlib, set up rotate = True so the graph spins.
    
    output_df : boolean, default False
        If true, this method will return the dataframe with the data.
        
    **kwargs : other keyword arguments
        gof : int, default is 1
                 cutoff for determining gain of function mutations based on mutagenesis data.
        lof : int, default is -1
            cutoff for determining loss of function mutations based on mutagenesis data.

    Returns
    ---------
    df_items : pandas dataframe
        Contains the plotted data. Needs to have output_df set to true.
    '''

    # Load parameters
    _parameters()

    # Update kwargs
    temp_kwargs = copy.deepcopy(default_kwargs)
    temp_kwargs.update(kwargs)
    temp_kwargs['x_label'] = kwargs.get('x_label', plot[0])
    temp_kwargs['y_label'] = kwargs.get('y_label', plot[1])
    temp_kwargs['z_label'] = kwargs.get('z_label', plot[2])
    
    # Get Scores and colors
    df_scores = _color_3D_scatter(self.dataframe, mode, temp_kwargs['lof'],
                           temp_kwargs['gof'])
    
    # If coordinates is not an input, get it from the pdb
    df_items = _parse_pdbcoordinates(self, pdb_path, position_correction,
                                     chain, sasa=True)
    
    # Add scores
    df_items['Score'] = list(df_scores['Score'])
    
    
    # Custom data
    if custom is not None:
        df_items[plot[2]] = custom
        
    # Plot figure
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    if df_color is None and color_by_score is True:
        c = df_scores['Color']
    elif df_color is None and color_by_score is False:
        c = 'k'
    else:
        c = df_scores['Color']
    
    ax.scatter(df_items[plot[0]], df_items[plot[1]], df_items[plot[2]], c=c)

    # axis labels
    ax.set_xlabel(temp_kwargs['x_label'])
    ax.set_ylabel(temp_kwargs['y_label'])
    ax.set_zlabel(temp_kwargs['z_label'])
    
    # axis scales
    ax.set_xscale(axis_scale[0])
    ax.set_yscale(axis_scale[1])
    ax.set_zscale(axis_scale[2])
    
    # save file
    _savefile(fig, temp_kwargs)

    # rotate if you are using interactive plot
    if rotate is True:
        for angle in range(0, 360):
            ax.view_init(30, angle)
            plt.draw()
            plt.pause(.001)

    if temp_kwargs['show']:
        plt.show()
    
    if output_df:
        return df_items, df_scores
    


# ## Map into Pymol

# In[ ]:


def plot_pymol(self, pdb, mode = 'mean', residues=None, position_correction = 0,
               quit=False, **kwargs):
    '''
    Color pymol structure residues. User can specify the residues to color, or can use the mutagenesis data.
    Activating mutations will be colored red and loss of function blue. Neutral mutations in green.
    Only works if pymol is your $PATH as pymol or you can start PyMOL in server mode.
    Uses the ipymol package, which needs to be installed from Github $pip install git+https://github.com/cxhernandez/ipymol , not from pypi (not updated there).
    
    Parameters
    ----------
    pdb : str
        User should specify the PDB chain in the following format 4G0N_A.
        If you have internet connection, Pymol will download the pdb. Otherwise,
        include the path were your PDB is stored locally.
        
    mode : str, default 'mean'
        Specify what enrichment scores to use. If mode = 'mean', it will use the mean of 
        each position to classify the residues. If mode = 'A', it will use the Alanine substitution profile. Can be 
        used for each amino acid. Use the one-letter code and upper case.
        
    residues : list , optional
        If user decides to pass custom arguments, use the following format
        residues = ['1,2,3,4-10','12-15,23,24,35','48,49,50,52-60'] which are [blue,red,green].
    
    position_correction : int, default 0
        If the pdb structure has a different numbering of positions than you dataset,
        you can correct for that. If your start_position = 2, but in the PDB that same residue
        is at position 20, position_correction needs to be set at 18.
    
    quit : boolean, default False
        if quit, close pymol after executing code.
    
    **kwargs : other keyword arguments
         gof : int, default is 1
             cutoff for determining gain of function mutations based on mutagenesis data.
         lof : int, default is -1
             cutoff for determining loss of function mutations based on mutagenesis data.
    Returns
    ----------
    Open pymol session with a fetched pdb structure where the residues are colored according to the enrichment scores.
    '''
    # update kwargs
    temp_kwargs = copy.deepcopy(default_kwargs)
    temp_kwargs.update(kwargs)

    # Calculate residues only if they are not given by the user
    if residues is None:
        residues = _pymol_fitness(self.dataframe.copy(), temp_kwargs['gof'], 
                                  temp_kwargs['lof'], mode, position_correction)

    # Start Pymol
    if not pymol._process_is_running():
        pymol.start()

    # Fetch structure. If pdb contains a "/", it will assume it is stored locally
    if '/' in pdb:
        pymol.load(pdb)
        pdb = (path.basename(pdb)).partition('.')[0] # Extract filename from pdb and then extract pdb code
    else:
        pymol.fetch(pdb)
    
    # Hide everything
    pymol.do('hide everything')

    # Selection names
    blue = pdb + '_blue'
    red = pdb + '_red'
    white = pdb + '_white'

    # Do selections
    pymol.select(blue, 'resi ' + residues[0])
    pymol.select(red, 'resi ' + residues[1])
    pymol.select(white, 'resi ' + residues[2])

    # Representation parameters
    pymol.show_as('cartoon', pdb)
    pymol.set('cartoon_color', 'neptunium', blue)
    pymol.set('cartoon_color', 'red', red)
    pymol.set('cartoon_color', 'chlorine', white)
    pymol.bg_color('white')
    pymol.remove('solvent')

    # light parameters
    _light_parameters()
    
    # deselect everything
    pymol.deselect()

    if quit:
        pymol.quit()
    return

# Convert fitness scores into pymol residues


def _pymol_fitness(df, gof, lof, mode, position_correction):
    '''You input the dataframe. Removes stop codons. 
    Returns the positions that are going to be colored blue,red and white'''
    
    # Select grouping
    if mode == 'mean':
        df_grouped = df.groupby(['Position'], as_index=False).mean()
    else:
        df_grouped = df.loc[df['Aminoacid']==mode]
        
    # Color of mutations
    blue_mutations = df_grouped[df_grouped['Score'] < lof]
    red_mutations = df_grouped[df_grouped['Score'] > gof]
    white_mutations = df_grouped[df_grouped['Score'].between(
        lof, gof, inclusive=True)]

    # Pymol Format
    blue_pymol = _array_to_pymol(blue_mutations['Position']+position_correction)
    red_pymol = _array_to_pymol(red_mutations['Position']+position_correction)
    white_pymol = _array_to_pymol(white_mutations['Position']+position_correction)

    residues = [blue_pymol, red_pymol, white_pymol]
    
    # If one group does not have any position, color position 0. Otherwise it gives an error
    for i, residue in enumerate(residues):
        if residue == '':
            residues[i] = '0'

    return residues


def _array_to_pymol(array):
    '''Input an array with positions of aminoacids, return it in pymol format'''
    pymol = ''
    for aminoacid in array:
        pymol += str(aminoacid)+'+'

    # delete last '+'
    pymol = pymol[:-1]
    return pymol


def _light_parameters():
    '''Group the light and ray parameters for pymol figures'''
    # Light parameters
    pymol.set('antialias', '3')
    pymol.set('ambient', '0.15')
    pymol.set('spec_count', '5')
    pymol.set('shininess', '50')
    pymol.set('specular', '0')
    pymol.set('light_count', '4')
    pymol.set('direct', '0.45')
    pymol.set('reflect', '0.5')
    pymol.set('opaque_background', 'off')
    pymol.set('dash_gap', 0.5)
    pymol.set('dash_radius', 0.1)

    # Stick parameters
    pymol.set('stick_radius', '0.2')
    pymol.set('sphere_scale', '0.2')
    pymol.set('sphere_quality', '4')
    return


# # Internal Functions

# In[ ]:


def _transform_dataset(dataset, sequence, aminoacids, start_position, fillna):
    '''
    Internal function that constructs a dataframe from user inputs

    Parameters
    -----------
    dataset, sequence, aminoacids, start_position,fillna

    Returns
    --------
    Dataframe containing [Position, Sequence, Aminoacid, Variant, Score]
    '''

    # make a dataframe
    df = pd.DataFrame()

    # Define Columns
    df['Sequence'] = np.ravel([[aa]*len(aminoacids) for aa in sequence])

    # Create column with position label
    df['Position'] = np.ravel(
        [[i]*len(aminoacids) for i in range(start_position, len(dataset[0])+start_position)])
    df['Aminoacid'] = aminoacids * len(dataset[0])
    df['Variant'] = df['Sequence']+df['Position'].astype(str)+df['Aminoacid']
    df['Score'] = np.ravel(dataset.T)
    df['Score_NaN'] = np.ravel(dataset.T)

    # Eliminate NaNs
    df['Score'].fillna(fillna, inplace=True)

    # Eliminate stop codons
    df_clean = df[df['Aminoacid'] != '*'].copy()

    return df, df_clean


def _transform_sequence(dataset, sequence, start_position):
    '''
    Internal function that trims the input sequence

    Parameters
    -----------
    dataset, sequence, start_position

    Returns
    --------
    string containing trimmed sequence
    '''

    # truncate sequence
    trimmedsequence = sequence[start_position -
                               1:len(dataset[0])+start_position-1]

    return trimmedsequence


def _transform_secondary(dataset, secondary, start_position, aminoacids):
    '''
    Internal function that trims the input secondary structure

    Parameters
    -----------
    dataset, sequence, start_position

    Returns
    --------
    list containing trimmed secondary structure (20 times each element)
    '''

    # Convert lists of lists to list
    secondary_list = list(itertools.chain.from_iterable(secondary))

    # Truncate list
    trimmedsecondary = secondary_list[start_position -
                                      1:len(dataset[0])+start_position-1]

    # Multiply each element by number of aminoacids. not use stop codon
    aminoacids = list(np.copy(aminoacids))
    if '*' in aminoacids:
        aminoacids.remove('*')
    secondary_dup = [x for item in trimmedsecondary for x in itertools.repeat(
        item, len(aminoacids))]

    return trimmedsecondary, secondary_dup


def _convert_to_df(dataset, sequence, aminoacids, startposition):
    '''
    Convertds np.array with stored enrichment scores into a dataframe
    Makes a copy of data

    Returns dataframe

    '''
    df = pd.DataFrame()
    df['Aminoacid'] = list(aminoacids) * len(dataset[0])
    df['Position'] = np.ravel(
        [[i]*len(aminoacids) for i in range(startposition, len(dataset[0])+startposition)])
    df['Sequence'] = np.ravel([[i]*len(aminoacids)
                               for i in sequence[:len(dataset[0])]])
    df['Score'] = np.copy(dataset.T).ravel()
    return df


def _df_rearrange(df, new_order, values='Score',show_snv = False):
    '''
    convert a df into a numpy array for mutagenesis data. 
    Allows the option of keeping NaN scores

    Returns copy
    '''
    dfcopy = df.copy()
    
    # If only SNVs, turn rest to NaN
    if show_snv is True:
        dfcopy.loc[dfcopy['SNV?']==False, values] = np.nan
    
    df_pivoted = dfcopy.pivot_table(values=values, index='Aminoacid',
                                    columns=['Position'], dropna=False)
    df_reindexed = df_pivoted.reindex(index=list(new_order))

    return df_reindexed


def _common(a, b):
    '''
    return common elements of two lists
    '''
    c = [value for value in a if value in b]
    return c


def _transpose(df, values='Score'):
    '''
    convert a df into a numpy array for mutagenesis data

    Returns copy
    '''
    df = df.pivot_table(values=values, index='Aminoacid',
                        columns=['Position']).T
    return df


def _select_aa(df, selection, values='Score'):
    '''returns copy'''
    df = _transpose(df.copy(), values)

    df = df[selection].T

    return df


def _savefile(fig, temp_kwargs):
    '''Save file function'''
    if temp_kwargs['savefile'] is True:
        filename = temp_kwargs['outputfilepath'] +             temp_kwargs['outputfilename']+"."+temp_kwargs['outputformat']
        fig.savefig(filename, format=temp_kwargs['outputformat'],
                    bbox_inches='tight', dpi=temp_kwargs['dpi'], transparent=True)
    return

def parse_pivot(df_imported, col_variant = 'variant', col_data = 'DMS',
               fill_value = np.nan):
    '''
    Parses a dataframe that contains saturation mutagenesis data in the Variant/Scores format.
    
    Parameters
    -----------
    df_imported : pandas dataframe
        Dataframe with the data imported with pd.read_excel.
    
    col_variant : str, default 'variant'
        Name of the column that contains the variants (ie T31A).
        
    col_data : str, default 'DMS'
        Name of the column that contains the saturation mutagenesis scores.
    
    fill_value : float, default np.nan
        What number to replace values that are omitted. It is possible that your 
        dataset does not have a wt value.
        
    Returns
    --------
    df_pivoted : pandas dataframe
        Dataframe that has been pivoted. Values are the saturation mutagenesis data. Columns are 
        the amino acid substitutions. Rows are the positions of the protein.
    
    sequence : list
        List of the amino acids that form the protein sequence.
    '''
    
    # Copy
    df = df_imported.copy()
    
    # Extract position and amino acids that are being mutated
    df['Position'] = df[col_variant].str.extract('(\d+)').astype(int)
    df['Original'] = df[col_variant].str[0:1]
    df['Substitution'] = df[col_variant].str[-1:]
    
    # Get sequence
    sequence = list(df.groupby(by=['Position', 'Original'], as_index=False, group_keys=False).sum() ['Original'])

    # Pivot
    df_pivoted = df.pivot_table(index='Substitution',columns = 'Position', 
                                values=col_data, fill_value = fill_value, dropna=False)
    
    return df_pivoted, sequence


# # Kwargs

# In[ ]:


def kwargs():
    '''
    Kwargs used in the package. Not all of them work on each function.
    
    Parameters
    -----------
    colormap : cmap, default custom bluewhitered
        Used for heatmaps. You can use your own colormap or the ones provided by 
        matplotlib. Example colormap = copy.copy((plt.cm.get_cmap('Blues_r')))

    colorbar_scale: list, default [-1, 1]
        Scale min and max used in heatmaps and correlation heatmaps.
    
    color: str, default 'k'
        Color used for the ...
        
    title : str, default 'Title'
        Title of plot.
        
    x_label : str, default 'x_label'
        Label of x axis.
        
    y_label : str, default 'y_label'
        Label of y axis.
    
    xscale: tuple, default (None, None)
        MinMax of x axis.
        
    yscale: tuple, default (None, None)
        MinMax of y axis.

    tick_spacing: int, default 1
        Space of axis ticks. Used for scatter and cumulative plots.
        
    inputfilepath : str, default ''
        Path of the input file.
        
    inputfilename : str, default ''
        Name of the input file.
        
    outputfilepath : str, default ''
        Path where file will be exported to.
        
    outputfilename : str, default ''
        Name of the exported file.
        
    outputformat': str, default 'png'
        Fortmat to export matplotlib object.
        
    dpi : int, default 600
        Dots Per Inch in the created image.
        
    neworder_aminoacids: list, default list('DEKHRGNQASTPCVYMILFW*')
        Order of amino acids to display in heatmaps. Used for heatmaps.
    
    savefile: boolean, default False
        If true, will save the matplotlib object into the specified folder.
        
    gof: int, default 1
        Cutoff of the enrichment score to classify a mutation as gain of function.
        Used on pymol function.
        
    lof: int, default -1
        Cutoff of the enrichment score to classify a mutation as loss of funtion.
        Used on pymol function.
    
    color_gof : str, default 'red'
        Color to color mutations above the gof cutoff.
        Used in pymol and mean methods.

    color_lof : str, default 'blue'
        Color to color mutations below the lof cutoff.
        Used in pymol and mean methods.
        
    cartoon_colors: list, default ['lightgreen', 'lavender', 'k']
        Colors used for secondary structure cartoon. Used for heatmap, mean and mean_count plots.
        
    text_labels: str, default 'None'
        Text labels that you can add to mean and mean_count plots. You will need to specify the coordinates.
        
    show: boolean, default True
        Whether to execute plt.show() or not on a matplotlib object.
        
    random_state : int, default 554
        Random state used for PCA function.
      
    '''
    # Do nothing, only so sphinx adds this to the rst file
    return

default_kwargs = {'colormap': generatecolormap(),
                  'colorbar_scale': [-1, 1],
                  'color': 'k',
                  'title': 'Title',
                  'x_label': 'x_label',
                  'y_label': 'y_label',
                  'z_label': 'y_label',
                  'xscale': (None, None),
                  'yscale': (None, None),
                  'tick_spacing': 1,
                  'inputfilepath': '',
                  'inputfilename': '',
                  'outputfilepath': '',
                  'outputfilename': '',
                  'outputformat': 'png',
                  'dpi': 600,
                  'aminoacids': list('ACDEFGHIKLMNPQRSTVWY*'),
                  'neworder_aminoacids': list('DEKHRGNQASTPCVYMILFW*'),
                  'savefile': False,
                  'gof': 1,
                  'lof': -1,
                  'color_gof' : 'red',
                  'color_lof' : 'blue',
                  'cartoon_colors': ['lightgreen', 'lavender', 'k'],
                  'text_labels': None,
                  'show': True,
                  'random_state' : 554,
                  'bins' : 50,
                  }


# # Define Class

# In[ ]:


class Screen:
    '''
    *Screen* represents a saturation mutagenesis experiment, where every amino acid 
    in the protein has been mutated to other amino acids. The mutants are scored based
    on some custom screen.

    Parameters
    -----------
    dataset : array
        2D matrix containing the enrichment scores of the point mutants. Columns will contain the
        amino acid substitutions, rows will contain the enrichment for each residue in the protein sequence.
    
    sequence : str
        Protein sequence (columns) in 1 letter code format.
    
    aminoacids : list, default list('ACDEFGHIKLMNPQRSTVWY*')
        Amino acid substitutions (rows). Submit in the same order that is used for the array.
    
    start_position : int, default 2
        First position in the protein sequence that will be used for the first column of the
        array. If a protein has been mutated only from residue 100-150, then if start_position = 100,
        the algorithm will trim the first 99 amino acids in the input sequence. The last 
        residue will be calculated based on the length of the input array. 
    
    secondary : list, optional
        This parameter is used to group the data by secondary structure. The format is 
        the name of the secondary structure multiplied by the residue length of that motif
        example : [['β1']*(8),['L1']*(7),['α1']*(9),...,].
    
    roc_df: Pandas dataframe, optional
        A dataframe that contains a column of variants labeled 'Variant' with a column labeled 'Class'
        containing the true class of that mutation. This can be used to compare enrichment scores to some label (such as 
        pathogenicity as found in a Cancer database) using ROC AUC.
    
    fillna : int, default 0
        How to replace NaN values.
    
    
    Attributes
    ------------
    dataframe : pandas dataframe
        Contains the enrichment scores, position, sequence.
    
    Other attributes are same as input parameters: dataset, aminoacids, start_position, roc_df, secondary

    '''

    def __init__(self, dataset, sequence, aminoacids=list('ACDEFGHIKLMNPQRSTVWY*'),
                 start_position=2, fillna=0, secondary=None, roc_df=None):
        self.dataset = np.array(dataset)
        self.aminoacids = aminoacids
        self.start_position = start_position
        self.end_position = len(self.dataset[0])+start_position
        self.sequence_raw = ''.join(sequence)
        self.sequence = _transform_sequence(
            self.dataset, self.sequence_raw, self.start_position)
        self.dataframe_stopcodons, self.dataframe = _transform_dataset(
            self.dataset, self.sequence, self.aminoacids, self.start_position, fillna)
        self.dataframe_SNV = _select_SNV(self.dataframe)
        self.dataframe_nonSNV = _select_nonSNV(self.dataframe)

        # Optional parameters
        self.roc_df = roc_df
        self.secondary = secondary
        if self.secondary is not None:
            self.secondary, self.secondary_dup = _transform_secondary(
                self.dataset, self.secondary, self.start_position, self.aminoacids)

    # Methods (Associated functions)
    kernel = plot_kernel
    heatmap = plot_heatmap
    heatmap_rows = plot_heatmap_rows
    heatmap_columns = plot_heatmap_columns
    mean = plot_mean
    meancounts = plot_meancounts
    differential = plot_meandifferential
    position = plot_position
    scatter = plot_scatter
    rank = plot_rank
    histogram = plot_hist
    miniheatmap = plot_miniheatmap
    neighboreffect = plot_neighboreffect
    correlation = plot_correlation
    individual_correlation = plot_individual_correlation
    group_correlation = plot_group_correlation
    pca = plot_pca
    secondary_mean = plot_secondary
    roc = plot_roc
    cumulative = plot_cumulative
    scatter_3D = plot_scatter_3D
    scatter_3D_pdbprop = plot_scatter_3D_pdbprop
    pymol = plot_pymol


# # Demo

# In[ ]:


def demo(figure='heatmap'):
    """
    Performs a demonstration of the mutagenesis_visualization software.

    Parameters
    -----------
    figure : str, default 'heatmap'
        There are 5 example plots that can be displayed to test the package is working on your station.
        The 5 options are 'heatmap', 'miniheatmap', 'mean', 'kernel' and 'pca'. Check the documentation for more information.

    Returns
    -------
    None.
    """
    # Use relative file import to access the data folder
    location = os.path.dirname(os.path.realpath(__file__))
    my_file = os.path.join(location, 'data', 'HRas166_RBD.csv')

    # Load enrichment scores
    hras_enrichment_RBD = np.genfromtxt(my_file, delimiter=',')

    # Define protein sequence
    hras_sequence = 'MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDLAARTVESRQAQDLARSYGIPYIETSAKTRQGVEDAFYTLVREIRQHKLRKLNPPDESGPG'

    # Define secondary structure
    secondary = [['L0'], ['β1']*(9-1), ['L1']*(15-9), ['α1']*(25-15), ['L2']*(36-25), ['β2']*(46-36), ['L3']*(48-46), 
                 ['β3']*(58-48), ['L4'] * (64-58),['α2']*(74-64), ['L5']*(76-74), ['β4']*(83-76), 
                 ['L6']*(86-83), ['α3']*(103-86), ['L7']*(110-103), ['β5']*(116-110), ['L8']*(126-116), ['α4']*(137-126),
                 ['L9']*(140-137), ['β6']*(143-140), ['L10']*(151-143), ['α5']*(172-151), ['L11']*(190-172)]

    # Create object
    hras_RBD = Screen(dataset=hras_enrichment_RBD,
                      sequence=hras_sequence, secondary=secondary)

    if figure == 'heatmap':
        # Create heatmap plot
        hras_RBD.heatmap(title='H-Ras 2-166', show_cartoon=True)
    elif figure == 'miniheatmap':
        # Condensed heatmap
        hras_RBD.miniheatmap(title='Wt residue H-Ras')
    elif figure == 'mean':
        # Mean enrichment by position
        hras_RBD.mean(figsize=[6, 2.5], mode='mean',show_cartoon=True, yscale=[-2, 0.5], title = '')
    elif figure == 'kernel':
        # Plot kernel dist using sns.distplot.
        hras_RBD.kernel(histogram=True, title='H-Ras 2-166', xscale=[-2, 1])
    elif figure == 'pca':
        # PCA by amino acid substitution
        hras_RBD.pca(dimensions=[0, 1], figsize=(2, 2), adjustlabels=True, title = '')
    return
    
def demo_datasets():
    '''
    Loads example datasets so the user can play with it.
    
    Parameters
    -----------
    None
    
    Returns
    --------
    data_dict : dictionary
        Dictionary that contains the datasets used to create the plots on the documentation.
    '''
    
    # Use relative file import to access the data folder
    location = os.path.dirname(os.path.realpath(__file__))
    
    # Create dictionary where to store data
    data_dict = {}
    
    # Retrieve H-Ras dataset and store in dict
    my_file = os.path.join(location, 'data', 'HRas166_RBD.csv')
    hras_enrichment_RBD = np.genfromtxt(my_file, delimiter=',')
    data_dict['array_hras'] = hras_enrichment_RBD
    
    # Beta lactamase data
    my_file = os.path.join(location, 'data', 'df_bla_raw.pkl')
    df_bla_raw = pd.read_pickle(my_file)
    data_dict['df_bla'], sequence_bla = parse_pivot(df_bla_raw, col_data = 'DMS_amp_625_(b)')

    # Sumo
    my_file = os.path.join(location, 'data', 'df_sumo1_raw.pkl')
    df_sumo1_raw = pd.read_pickle(my_file)
    data_dict['df_sumo1'], sequence_sumo1 = parse_pivot(df_sumo1_raw, col_data = 'DMS')
   
    # MAPK1
    my_file = os.path.join(location, 'data', 'df_mapk1_raw.pkl')
    df_mapk1_raw = pd.read_pickle(my_file)
    data_dict['df_mapk1'], sequence_mapk1 = parse_pivot(df_mapk1_raw, col_data = 'DMS_DOX')
    
    #UBE2I
    my_file = os.path.join(location, 'data', 'df_ube2i_raw.pkl')
    df_ube2i_raw = pd.read_pickle(my_file)
    data_dict['df_ube2i'], sequence_ube2i = parse_pivot(df_ube2i_raw, col_data = 'DMS')

    #TAT
    my_file = os.path.join(location, 'data', 'df_tat.pkl')
    data_dict['df_tat'] = pd.read_pickle(my_file)

    #REV
    my_file = os.path.join(location, 'data', 'df_rev.pkl')
    data_dict['df_rev'] = pd.read_pickle(my_file)
    
    # asynuclein
    my_file = os.path.join(location, 'data', 'df_asynuclein.pkl')
    data_dict['df_asynuclein'] = pd.read_pickle(my_file)
    
    # APH
    my_file = os.path.join(location, 'data', 'df_aph.pkl')
    data_dict['df_aph'] = pd.read_pickle(my_file)

    # b11L5
    my_file = os.path.join(location, 'data', 'df_b11L5F_raw.pkl')
    df_b11L5F_raw = pd.read_pickle(my_file)
    data_dict['df_b11L5F'], sequence_b11L5F = parse_pivot(df_b11L5F_raw, col_data = 'relative_tryp_stability_score')
    
    return data_dict


# # Ras Trouble Shooting

# In[ ]:


'''# Load enrichment scores into a np.array
hras_enrichment = np.genfromtxt('../data/HRas166_RBD.csv', delimiter=',')

# Define protein sequence
hras_sequence = 'MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDLAARTVESRQAQDLARSYGIPYIETSAKTRQGVEDAFYTLVREIRQHKLRKLNPPDESGPG'

# Order of amino acid substitutions in the hras_enrichment dataset
aminoacids = list('ACDEFGHIKLMNPQRSTVWY*')

# First residue of the hras_enrichment dataset. Because 1-Met was not mutated, the dataset starts at residue 2
start_position = 2

# Define secondary structure
secondary = [['L0'], ['β1']*(9-1), ['L1']*(15-9), ['α1']*(25-15), ['L2']*(36-25), ['β2']*(46-36),['L3']*(48-46), ['β3']*(58-48), ['L4'] *(64-58), ['α2']*(74-64), ['L5']*(76-74), ['β4']*(83-76), ['L6']*(86-83), ['α3']*(103-86), ['L7']*(110-103), ['β5']*(116-110), ['L8']*(126-116), ['α4']*(137-126), ['L9']*(140-137), ['β6']*(143-140), ['L10']*(151-143), ['α5']*(172-151), ['L11']*(190-172)]

# Create object
hras_object = Screen(hras_enrichment,hras_sequence,aminoacids,start_position,0,secondary)
'''


# In[ ]:


'''path = '../Data/DMS_others.xlsx'
sheet_name='env'
usecols='A:B'

# Read excel file
df = pd.read_excel(path, sheet_name, usecols=usecols)

# Parse
df_env, sequence_env = parse_pivot(df)'''


# In[ ]:


'''# Order of amino acid substitutions in the hras_enrichment dataset
aminoacids = list(df_env.index)
neworder_aminoacids = list('DEKHRGNQASTPCVYMILFW')

# First residue of the hras_enrichment dataset. Because 1-Met was not mutated, the dataset starts at residue 2
start_position = df_env.columns[0]

sequence = ['X']*start_position+sequence_env
# Create objects
env_obj = Screen(df_env, sequence,
                         aminoacids, start_position)'''

