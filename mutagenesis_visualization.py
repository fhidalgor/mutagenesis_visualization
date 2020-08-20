#!/usr/bin/env python
# coding: utf-8

# # Import Modules

# In[1]:


from __future__ import unicode_literals
import numpy as np
import seaborn as sns
import pandas as pd
import copy
from scipy import stats
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.mlab as pl
import matplotlib.ticker as ticker
from sklearn.decomposition import PCA
from adjustText import adjust_text
import itertools
from ipymol import viewer as pymol
from pymol import cmd, stored, math
from collections import defaultdict
from sklearn import metrics
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from collections import OrderedDict,Counter
import logomaker


# # Data Process Functions

# ## Process ILMN reads

# In[170]:


def count_reads(wtsequence, inputfilepath, inputfilename, 
                                 outputfilepath,outputfilename, codonList = 'NNS',save_file=False):
    '''
    Process a trimmed fastq file to return the library counts
    
    Parameters
    -----------
    wtsequence : str, contains the DNA sequence of the wt allele
    inputfilepath : str, path where the fastq file is saved
    inputfilename : str, name of the fastq file (full name including .fastq)
    outputfilepath : str, path where the output files are saved
    outputfilename : str, name that the output files will have
    codonList : list, List of the codons that are mutated. For NNS they are already in the software
    save_file : boolean, default None. Make true if you want to export the files
    
    Returns
    --------
    Dataframe with counts for each point mutant and list of the counts for each synonymous wt allele
    If save_file = True, it will export two txt files containing an array with the counts 
    for each point mutant and an array with the counts for each synonymous wt allele.
    Will also export the dataframe before pivoting???
    '''
    
    # Naming files
    trimmedfile = inputfilepath+inputfilename
    outputfile_counts = outputfilepath+outputfilename+"_counts.txt"
    outputfile_wtcounts = outputfilepath+outputfilename+"_wtcounts.txt"

    # Create list with codons of sequence
    wtsequence = wtsequence.upper()
    wtSeqList = [wtsequence[i:i+3] for i in range(0, len(wtsequence), 3)]
    
    # CodonList
    if codonList == 'NNS':
        codonList = ["GCC","GCG","TGC","GAC","GAG","TTC","GGC","GGG","CAC","ATC","AAG","CTC","CTG","TTG","ATG","AAC","CCC","CCG","CAG","CGC","CGG","AGG","TCC","TCG","AGC","ACC","ACG","GTC","GTG","TGG","TAC","TAG"]
    elif codonList == 'NNK':
        codonList = ['GCG','GCT','TGT','GAT','GAG','TTT','GGG','GGT','CAT','ATT','AAG','CTG','CTT','TTG','ATG','AAT','CCG','CCT','CAG','AGG','CGG','CGT','AGT','TCG','TCT','ACG','ACT','GTG','GTT','TGG','TAT','TAG']
    
    # Enumerate variants
    variants=OrderedDict()
    isitthefirstwtseq = False
    for position, codon in enumerate(wtSeqList):
        for position2,codons in enumerate(codonList): 
            variant = ''.join(wtSeqList[0:position]) + ''.join(codons) + ''.join(wtSeqList[position+1:])
            if (variant == wtsequence):
                if isitthefirstwtseq:
                    variant='wtSeq'+ str(position)
                isitthefirstwtseq=True
            variants[variant] = 0

    # Translate nucleotide sequence and count variant frequency
    totalreads = 0
    for nuc in SeqIO.parse(trimmedfile, "fastq"):
        totalreads +=1
        nucleicsequence = str(nuc.seq)
        if nucleicsequence in variants:
            variants[nucleicsequence] += 1
    usefulreads=np.nansum(list(variants.values()))
    
    # Convert to df
    wtProtein=Seq(wtsequence,generic_dna).translate()
    df = pd.DataFrame()
    df['Position'] =  np.ravel([[pos]*len(codonList) for pos in np.arange(1,len(wtProtein)+1).astype(int)])
    df['Codon'] = codonList*len(wtProtein)
    df['WTCodon'] = np.ravel([[codon]*len(codonList) for codon in wtSeqList])
    df['Aminoacid'] = np.ravel([[aa]*len(codonList) for aa in wtProtein])
    codontable = codon_table()
    df['SynWT'] = df.apply(lambda x: _are_syn(x['Codon'],x['WTCodon'],codontable), axis=1)
    df['Counts'] = list(variants.values())
    df.loc[df['Codon']==df['WTCodon'],'Counts'] = variants[wtsequence]
    
    # Pivot table and reindex
    frequencyMatrix = df.pivot_table(values='Counts', index='Codon',columns=['Position'],dropna=False)
    frequencyMatrix = frequencyMatrix.reindex(index = codonList) 
    
    # Get WT counts
    df_wt = df.loc[df['SynWT']==True]
    wt_counts = list(df_wt['Counts'])
    wt_counts.insert(0,int([variants[wtsequence]][0]))
    
    # Export files
    if save_file:
        np.savetxt(outputfile_counts, np.array(frequencyMatrix), fmt = '%i', delimiter='\t')
        np.savetxt(outputfile_wtcounts, wt_counts, fmt = '%i', delimiter='\t')
    
    # Print total reads
    print ('{}/{} useful reads ({}%)'.format(str(usefulreads),str(totalreads),str(int(usefulreads/totalreads*100))))
    return frequencyMatrix, wt_counts

def codon_table():
    codontable = {'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
            'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 
            'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'}
    return codontable

def _are_syn(codon1, codon2, codontable):
    '''Determine if 2 codons are synonymous'''
    if codon1 == codon2:
        return False
    if _translate(codon1,codontable) is not _translate(codon2,codontable):
        return False
    return True

def _translate(seq, codontable): 
    protein ='' 
    if len(seq)%3 == 0: 
        for i in range(0, len(seq), 3): 
            codon = seq[i:i + 3] 
            protein+= codontable[codon] 
    return protein 


# ## Process counts

# In[3]:


def process_counts(input_lib, output_lib, input_wt=None, output_wt=None, aminoacids=list('AACDEFGGHIKLLLMNPPQRRRSSSTTVVWY*'),
                   zeroing='population', how='median', norm_std=True,stopcodon=False, min_counts=25,
                   min_countswt=100, mpop=2, mwt=2, infinite=3):
    '''
    Calculate enrichment scores of a selection experiment
    
    Returns np.array with the enrichment scores
    '''
    # Convert to numpy if libraries are in dataframe format
    if type(input_lib) is pd.DataFrame: input_lib = input_lib.to_numpy()
    if type(output_lib) is pd.DataFrame: output_lib = output_lib.to_numpy()

    # Locate stop codons
    input_stopcodon = input_lib[-1]
    output_stopcodon = output_lib[-1]

    # Log10 of the counts for library and wt alleles
    log10_counts = _get_enrichment(input_lib, output_lib, input_stopcodon,
                                   output_stopcodon, min_counts, stopcodon, infinite)
    # Group by amino acid
    df = pd.DataFrame(data=log10_counts)
    log10_counts_grouped = _group_byaa(df, aminoacids)

    # MAD filtering
    log10_counts_mad = _MAD_filtering(np.ravel(log10_counts_grouped), mpop)
    mean_pop = np.nanmean(log10_counts_mad)
    median_pop = np.nanmedian(log10_counts_mad)
    std_pop = np.nanstd(log10_counts_mad)
    mode_pop = _nanmode(log10_counts_mad)

    # Wt counts
    if input_wt is not None:
        log10_wtcounts = _get_enrichment(input_wt, output_wt, input_stopcodon,
                                         output_stopcodon, min_countswt, stopcodon, infinite)
        # MAD filtering
        # If set to m=1, if tosses out about 50% of the values. the mean barely changes though
        log10_wtcounts = _MAD_filtering(log10_wtcounts, mwt)
        mean_wt = np.nanmean(log10_wtcounts)
        median_wt = np.nanmedian(log10_wtcounts)
        std_wt = np.nanstd(log10_wtcounts)
        mode_wt = _nanmode(log10_wtcounts)

    # Zero data
    if zeroing is 'wt':
        if how is 'mean':
            zeroed = log10_counts_grouped - mean_wt
        elif how is 'median':
            zeroed = log10_counts_grouped - median_wt
        elif how is 'mode':
            zeroed = log10_counts_grouped - mode_wt
        if norm_std is True: zeroed = zeroed*0.1/std_wt
    elif zeroing is 'population':
        if how is 'mean':
            zeroed = log10_counts_grouped - mean_pop
        elif how is 'median':
            zeroed = log10_counts_grouped - median_pop
        elif how is 'mode':
            zeroed = log10_counts_grouped - mode_pop
        if norm_std is True: zeroed = zeroed*0.2/std_pop
    elif zeroing is 'kernel':
        zeroed_0, kernel_std = _kernel_correction(log10_counts_grouped)
        zeroed, kernel_std = _kernel_correction(zeroed_0, cutoff=1)
        if norm_std is True: zeroed = zeroed*0.2/kernel_std
    return zeroed


# ## Aux functions

# In[4]:


def _get_enrichment(input_lib, output_lib, input_stopcodon, output_stopcodon, 
                   min_counts, stopcodon, infinite):
    '''Calculate log10 enrichment scores from input and output counts'''
    # Copy data and replace low counts by np.nan
    input_lib = np.copy(input_lib.astype(float))
    output_lib = np.copy(output_lib.astype(float))
    input_lib[input_lib < min_counts] = np.nan

    # Stop codon correction
    if stopcodon:
        output_lib = _stopcodon_correction(input_lib, output_lib, input_stopcodon, output_stopcodon)

    # log10 of library and replace infinite values
    counts_log10_ratio = _replace_inf(np.log10(output_lib/input_lib),infinite)

    return counts_log10_ratio

def _stopcodon_correction(input_lib, output_lib,input_stopcodon,output_stopcodon):
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
    df = pd.DataFrame(data, columns=['Data'])
    median = df['Data'].median(axis=0)
    mad = df['Data'].mad(axis=0)
    df['Abs_Dev'] = np.abs(data-median)/mad

    # filter values m times away from median, by default m = 2
    df['Abs_Dev'].mask(df['Abs_Dev'] > m, inplace=True)  # mask values
    df.dropna(how='any', inplace=True)  # eliminte NANs

    return df['Data'].to_numpy()

def _replace_inf(array,infinite):
    '''Replace values over a threshold with a min or max value'''
    np.warnings.filterwarnings('ignore')
    array[array == -np.inf] = -infinite
    array[array < -infinite] = -infinite
    array[array == +np.inf] = +infinite
    array[array > +infinite] = +infinite
    return array

def _group_byaa (df,aminoacids):
    '''Group different codons that are synonymous'''
    # copy df
    df = df.copy()
    
    # Set up amino acid column
    df['Aminoacid'] = aminoacids
    
    # Group by mean
    df = df.groupby(as_index = False, by='Aminoacid').mean()
    
    # Sort amino acids by alphabetical order
    df['Aminoacid'] = pd.Categorical(df['Aminoacid'], list(dict.fromkeys(aminoacids)))
    df = df.sort_values(by=['Aminoacid'])
    
    # Drop amino acids column
    df.drop(columns='Aminoacid', inplace=True)
    return df.to_numpy() 

def _nanmode(data): 
    '''input is wt log enrichments, and return the mode of the histogram 
    (aka the x coordinate at which y is max)'''
    # Copy data
    data=np.copy(data)
    # Remove NaN values
    data_corrected=data[np.invert(np.isnan(data))]
    # Adjust kernel
    kernel_processed_data=stats.gaussian_kde(data_corrected)
    # Find mode
    indexmax = np.where(kernel_processed_data(data_corrected)==kernel_processed_data(data_corrected).max())
    # Return mean in case there are two x values with equal y-axis height
    return data_corrected[indexmax].mean()

def _kernel_correction(data, cutoff=2):  # corrects the mutagenesis data and returns the height of the peak
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
    # Copy data
    data_corrected = np.copy(data[0:20])

    # Eliminate values lower than -1
    data_corrected = data_corrected[(data_corrected >= -cutoff) & (data_corrected <= cutoff)]

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
    y_temp = kernel(data)  # so I can select only the positive side of the distribution

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


# # Plot Functions

# ## Kernel

# In[5]:


def plot_kernel(self, histogram=False, **kwargs):
    '''
    Generate a kernel density plot + histogram (optional).

    Parameters
    ----------
    histogram : boolean, default False
    **kwargs
        color : str, default 'k'
        scale : list, default [-4,2]
        figsize : default (2.5,2)
        title : str, default 'Title'
        start_position : int, default 2
        outputfilepath : str, default ''
        outputfilename : str, default ''
        outputformat : str, default 'png'
        savefile : boolean, default False

    Returns
    ----------
    Exported png file to desired folder
    '''
    # update kwargs
    temp_kwargs = copy.deepcopy(default_kwargs)
    temp_kwargs.update(kwargs)
    temp_kwargs['figsize'] = kwargs.get('figsize', (2.5,2))
    temp_kwargs['xscale'] = kwargs.get('xscale', (-2,2))
    
    # create figure
    fig = plt.figure(figsize=temp_kwargs['figsize'])

    # import parameters
    parameters()

    # plot
    sns.distplot(self.dataframe['Score'], kde=True, hist=histogram, norm_hist=True,
                 kde_kws={'color': temp_kwargs['color'], 'lw': 2})

    # tune graph
    plt.xlabel(r'$∆E^i_x$', fontsize=10,fontname='Arial', color='k', labelpad=0)
    plt.ylabel('Probability density', fontsize=10,fontname='Arial', color='k', labelpad=3)
    plt.title(temp_kwargs['title'], fontsize=12, fontname='Arial', color='k')
    plt.xlim(temp_kwargs['xscale'])
    plt.grid()
    
    # save file
    _savefile(fig,temp_kwargs)

    plt.show()
    return


def parameters():
    # normal font
    mpl.rcParams['font.family'] = 'sans-serif'
    mpl.rcParams['font.sans-serif'] = ['Arial']

    # math font
    mpl.rcParams['mathtext.fontset'] = 'custom'
    mpl.rcParams['mathtext.rm'] = 'Arial'
    mpl.rcParams['svg.fonttype'] = 'none'

    # add grid
    mpl.rcParams['grid.color'] = 'silver'
    mpl.rcParams['grid.linestyle'] = '--'
    mpl.rcParams['grid.linewidth'] = 1
    mpl.rcParams['lines.dashed_pattern'] = [5, 10]
    mpl.rcParams['axes.axisbelow'] = True
    # Parameters for all graphs
    mpl.rcParams['xtick.labelsize'] = 9
    mpl.rcParams['ytick.labelsize'] = 9
    return


# ## Heatmap

# ### Full Heatmap

# In[95]:


def plot_heatmap(self, nancolor='lime', show_cartoon=False, **kwargs):
    '''
    Generate a heatmap plot enrichment scores of mutagenesis selection assays.

    Parameters
    ----------
    **kwargs
        colormap : object, default bluewhitered
        colorbar_scale : list, default [-1,0,1]
            give only 3 values (min,midpoint,max)
        neworder_aminoacids : str, default 'DEKHRGNQASTPCVYMILFW*'
        title : str, default 'Title'
        start_position : int, default 2
        outputfilepath : str, default ''
        outputfilename : str, default ''
        outputformat : str, default 'png'
        savefile : boolean, default False

    Returns
    ----------
    Exported png file to desired folder
    '''
    # load font parameters
    font_parameters()

    # update kwargs
    temp_kwargs = copy.deepcopy(default_kwargs)
    temp_kwargs.update(kwargs)

    # load labels
    temp_kwargs['color_sequencelabels'] = _labels(self.start_position)[0]
    temp_kwargs['number_sequencelabels'] = _labels(self.start_position)[1]

    # sort data in specified order by user
    dataset = _df_rearrange(self.dataframe_stopcodons, temp_kwargs['neworder_aminoacids'],values='Score_NaN').to_numpy()

    # declare figure and subplots
    figwidth = 14*len(dataset[0])/165
    
    # Change parameters depending on whether cartoon is on or off
    if show_cartoon:
        figheight = 2.45
        fig = plt.figure(figsize=(figwidth, figheight))
        gs = gridspec.GridSpec(nrows=3, ncols=2, height_ratios=[
                               len(dataset), 1,5], width_ratios=[len(dataset[0]), 1])
    else:
        figheight = 2
        fig = plt.figure(figsize=(figwidth, figheight))
        gs = gridspec.GridSpec(nrows=2, ncols=2, height_ratios=[
                               len(dataset), 1], width_ratios=[len(dataset[0]), 1])
        

    ax = plt.subplot(gs[0, 0])
    averageresidue = plt.subplot(gs[1, 0])
    cbar1 = plt.subplot(gs[0, 1])
    cbar2 = plt.subplot(gs[1, 1])
    
    # Change color of values that are NaN
    cmap = temp_kwargs['colormap']
    cmap.set_bad(color=nancolor)
    
    # main heatmap
    heatmap = ax.pcolormesh(dataset, vmin=temp_kwargs['colorbar_scale'][0], vmax=temp_kwargs['colorbar_scale'][1],
                        cmap=cmap, edgecolors='k', linewidths=0.2,antialiased=True, color='darkgrey')

    # average heatmap
    heatmapaverageresidues = averageresidue.pcolor([np.nanmean(dataset[0:20], axis=0)], vmin=temp_kwargs['colorbar_scale'][0], vmax=temp_kwargs['colorbar_scale'][1],
                                                   cmap=temp_kwargs['colormap'], edgecolors='k', linewidths=0.2,antialiased=True, color='darkgrey')

    # ____________axes manipulation____________________________________________
    # put the major ticks at the middle of each cell
    ax.set_xticks(np.arange(dataset.shape[1]) + 0.5, minor=False)
    ax.set_yticks(np.arange(dataset.shape[0]) + 0.5, minor=False)

    # position of axis labels
    ax.tick_params('x', direction='out', pad=-2.5)
    ax.tick_params('y', direction='out', pad=0.4)

    # make new axes
    ax2 = ax.twiny()
    ax3 = ax.twinx()

    # tune the axes
    ax2.set_xticks(np.arange(dataset.shape[1]) + 0.5, minor=False)
    ax3.set_yticks(np.arange(dataset.shape[0]) + 0.5, minor=False)
    ax2.tick_params(direction='out', pad=4)
    ax3.tick_params(direction='out', pad=0.4)
    averageresidue.tick_params(direction='out', pad=-2)
    averageresidue.set_xticks(np.arange(dataset.shape[1]) + 0.5, minor=False,)
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
    ax.set_xticklabels(list(self.sequence), fontsize=6.5, fontname="Arial", color='k', minor=False)
    ax.set_yticklabels(temp_kwargs['neworder_aminoacids'], fontsize=6,
                       fontname="Arial", color='k', minor=False)
    ax2.set_xticklabels(temp_kwargs['number_sequencelabels'],
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
    cb.ax.set_yticklabels(cb.ax.get_yticklabels(), fontsize=6, fontname="Arial", color='k')
    cb.update_ticks()
    #plt.text(len(dataset[0])+4, 6.5, r'$\langle∆E^x_i\rangle_x$', horizontalalignment='center',fontsize=7,fontname="Arial",color='k')
    plt.text(1.2+10/len(dataset[0]), 0.7, r'$\langle∆E^x_i\rangle_x$', transform=cbar1.transAxes,
             horizontalalignment='center', fontsize=7, fontname="Arial", color='k')

    gs.update(hspace=0.1, wspace=0.1/len(dataset[0])*50)

    # for putting title on graph
    plt.title(temp_kwargs['title'], horizontalalignment='center', fontname="Arial", fontsize=12)

    # Cartoon
    if show_cartoon:
        _generate_cartoon(self,gs,2,temp_kwargs['cartoon_colors'],0.025)
    
    # save file
    _savefile(fig,temp_kwargs)
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


def font_parameters():
    # math font
    mpl.rcParams['mathtext.fontset'] = 'custom'
    mpl.rcParams['mathtext.rm'] = 'Arial'
    mpl.rcParams['svg.fonttype'] = 'none'
    return


def generatecolormap():
    cmap = LinearSegmentedColormap(
        'BlueWhiteRed',
        {
            'red':  ((0.0, 0.0, 0.0),
                     (0.15, 0.0, 0.0),
                     (0.475, 1.0, 1), (0.525, 1.0, 1), (0.85, 1.0, 1.0), (1.0, .8, 1)),
            'green': ((0.0, 0.0, .0),
                      (0.15, 0.5, 0.5), (0.475, 1.0, 1), (0.525, 1.0, 1), (0.85, 0.0, 0.0),
                      (1.0, 0.0, 0.0)),
            'blue': ((0.0, .5, .5),
                     (0.15, 1, 1),
                     (0.475, 1.0, 1), (0.525, 1.0, 1), (0.85, 0.0, 0.0), (1.0, 0.0, 0.0))
        },)
    return cmap

def _generate_cartoon(self,gs,n_row,colors,bottom_space=0,
                      fig_inches=13.91,show_labels = True):
    '''Generates cartoon for heatmap'''
    # Create subplot
    cartoon = plt.subplot(gs[n_row, 0])
        
    # Generate coordinates of labels
    labels = list(Counter(self.secondary).keys())
    length = list(Counter(self.secondary).values())
    cumsum = length[:-1]
    cumsum.insert(0,self.start_position)
    cumsum = np.cumsum(cumsum)
    
    # Create cartoon
    for label,length,cum in zip(labels,length,cumsum):
        if 'β' in label:
            loopstructure=_loop(cum,length,color=colors[2])
            cartoon.add_patch(loopstructure)
            sheetstructure=_sheet(cum,length,colors[0])
            cartoon.add_patch(sheetstructure)
            x_label = cum + length - 3.5
            if length > 2 and show_labels: # If beta sheet is too small, label does not fit
                if length == 3:
                    cartoon.text((x_label+0.6), -0.25,label,name='Arial',fontweight='normal',size=8.5*fig_inches/13.91,multialignment='right')
                else:
                    cartoon.text((x_label), -0.25,label,name='Arial',fontweight='normal',size=8.5*fig_inches/13.91,multialignment='right')
        elif 'α' in label:
            helixstructure=_helix(cum,length,colors[1])
            cartoon.add_patch(helixstructure)
            x_label = cum + length/2 - 1
            if length > 2 and show_labels:
                cartoon.text((x_label), -0.3,label,name='Arial',fontweight='normal',size=9*fig_inches/14,multialignment='center')
        elif 'L' in label:
            loopstructure=_loop(cum,length,colors[2])
            cartoon.add_patch(loopstructure)

    # format of secondary cartoon
    cartoon.xaxis.set_ticks_position('none')
    cartoon.yaxis.set_ticks_position('none')
    cartoon.axis('off')
    
    # size
    cartoon.set_xlim(self.start_position-0.1,len(self.secondary)+self.start_position+0.2)
    cartoon.set_ylim(-2,2.5)
    
    # adjust proximity to heatmap
    box = cartoon.get_position()
    box.y0 = box.y0-bottom_space
    box.y1 = box.y1-bottom_space
    cartoon.set_position(box)
    
    return

def _sheet(starting_aa,length_aa,color='lightgreen'):
    dx=length_aa
    sheetstructure=patches.FancyArrow(starting_aa, 0.25, dx, 0, width=2, length_includes_head=True, 
            head_width=4, head_length=3, shape='full', overhang=0, head_starts_at_zero=False,ec='k',fc=color)
    return sheetstructure

def _helix(starting_aa,length_aa,color='lavender'):
    dx=length_aa #so i can overlap tip
    helixstructure = plt.Rectangle((starting_aa, -0.85), dx, 2.2, fc=color,ec='k')
    return helixstructure

def _loop(starting_aa,length_aa,color='k'):
    dx=length_aa
    loopstructure = plt.Rectangle((starting_aa, 0), dx, 0.5, fc=color)
    return loopstructure


# ### Grouped Heatmap

# In[7]:


#input variables wont have "_" and function variables will have "_"
def plot_heatmap_selection(self, selection=['E','Q','A','P','V','Y'], nancolor='green', **kwargs):
    '''
    Generate a heatmap plot enrichment scores of mutagenesis selection assays. 
    Only selected amino acids will be displayed

    Parameters
    ----------
    show : list of aa to show, default ['E','Q','A','P','V','Y']. 
    **kwargs
        colormap : object, default bluewhitered
        colorbar_scale : list, default [-1,0,1]
            give only 3 values (min,midpoint,max)
        neworder_aminoacids : str, default 'DEKHRGNQASTPCVYMILFW*'
        title : str, default 'Title'
        start_position : int, default 2
        outputfilepath : str, default ''
        outputfilename : str, default ''
        outputformat : str, default 'png'
        savefile : boolean, default False

    Returns
    ----------
    Exported png file to desired folder
    '''
    # load font parameters
    font_parameters()

    # update kwargs
    temp_kwargs = copy.deepcopy(default_kwargs)
    temp_kwargs.update(kwargs)

    # load labels
    temp_kwargs['color_sequencelabels'] = _labels(self.start_position)[0]
    temp_kwargs['number_sequencelabels'] = _labels(self.start_position)[1]  

    # Add group and pivot df
    df = _select_aa(self.dataframe_stopcodons,selection,values='Score_NaN')
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
                        cmap=cmap, edgecolors='k', linewidths=0.2,antialiased=True, color='darkgrey')

    # put the major ticks at the middle of each cell
    ax.set_xticks(np.arange(dataset.shape[1]) + 0.5, minor=False)
    ax.set_yticks(np.arange(dataset.shape[0]) + 0.5, minor=False)

    #position of axis labels
    ax.tick_params('x',direction='out', pad=-2.5)
    ax.tick_params('y',direction='out', pad=0.4)

    #second axis
    ax2 = ax.twiny()
    ax2.set_xticks(np.arange(dataset.shape[1]) + 0.5, minor=False)
    ax2.tick_params(direction='out', pad=4)

    # Set the limits of the new axis from the original axis limits
    ax2.set_xlim(ax.get_xlim())

    # want a more natural, table-like display
    ax.invert_yaxis()
    ax.xaxis.tick_top()

    #so labels of x and y do not show up and my labels show up instead
    ax.set_xticklabels(list(self.sequence), fontsize=6.5, fontname="Arial", color='k', minor=False)
    ax.set_yticklabels(list(df.T.columns), fontsize=6, fontname="Arial", color='k', minor=False)
    ax2.set_xticklabels(temp_kwargs['number_sequencelabels'], fontsize=10, fontname="Arial", color='k', minor=False)
    
    #align the labels of the y axis 
    for ylabel in ax.get_yticklabels():
        ylabel.set_horizontalalignment('center')

    #for coloring the residues that are 10,20...
    for xtick, color in zip(ax.get_xticklabels(), temp_kwargs['color_sequencelabels']):
        xtick.set_color(color)

    # for putting title on graph
    plt.title(temp_kwargs['title'], horizontalalignment='center', fontname="Arial", fontsize=12)
    
    # for color bar format
    cbar1.axis('off')
    cb = plt.colorbar(heatmap, fraction=1, pad=0, ax=[cbar1], aspect=5, ticks=[temp_kwargs['colorbar_scale'][0], np.mean(
        temp_kwargs['colorbar_scale']), temp_kwargs['colorbar_scale'][1]], orientation='vertical')
    cb.ax.set_yticklabels(cb.ax.get_yticklabels(), fontsize=8, fontname="Arial", color='k')
    cb.update_ticks()
    gs.update(hspace=0.1, wspace=0.1/len(dataset[0])*50)

    #remove ticks
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')
    ax2.yaxis.set_ticks_position('none')
    ax2.xaxis.set_ticks_position('none')

    # save file
    _savefile(fig,temp_kwargs)
    
    plt.show()
    return


# ### Subset Heatmap

# In[8]:


def plot_heatmapsubset(self, segment, ylabel_color='k', nancolor='green', **kwargs):
    '''
    Generate a heatmap plot enrichment scores of mutagenesis selection assays.

    Parameters
    ----------
    segment : list
        Segment is typed as [20,40] and includes both residues 20 and 40
    ylabel_color : str, default 'k'
        choose white if you don't want amino acid y axis label
    **kwargs
        colormap : object, default bluewhitered
        colorbar_scale : list, default [-1,0,1]
            give only 3 values (min,midpoint,max)
        neworder_aminoacids : str, default 'DEKHRGNQASTPCVYMILFW*'
        title : str, default 'Title'
        start_position : int, default 2
        outputfilepath : str, default ''
        outputfilename : str, default ''
        outputformat : str, default 'png'
        savefile : boolean, default False

    Returns
    ----------
    Exported png file to desired folder
    '''

    # load font parameters
    font_parameters()

    # update kwargs
    temp_kwargs = copy.deepcopy(default_kwargs)
    temp_kwargs.update(kwargs)

    # load labels
    temp_kwargs['color_sequencelabels'] = _labels(self.start_position)[0]
    temp_kwargs['number_sequencelabels'] = _labels(self.start_position)[1]

    # sort data in specified order by user
    dataset = _df_rearrange(self.dataframe_stopcodons, temp_kwargs['neworder_aminoacids'],values='Score_NaN').to_numpy()

    # select subset
    dataset = dataset[:, segment[0]-self.start_position:segment[1]-self.start_position+1]

    # the size can be changed
    figwidth = 2*len(dataset[0])/22
    figheight = 2
    fig = plt.figure(figsize=(figwidth, figheight))
    gs = gridspec.GridSpec(nrows=1, ncols=1)
    ax = plt.subplot(gs[0])  # needed to set autoscale off to avoid missalignment
    
    # Change color of values that are NaN
    cmap = temp_kwargs['colormap']
    cmap.set_bad(color=nancolor)  
    
    # main heatmap
    heatmap = ax.pcolormesh(dataset, vmin=temp_kwargs['colorbar_scale'][0], vmax=temp_kwargs['colorbar_scale'][1],
                        cmap=cmap, edgecolors='k', linewidths=0.2,antialiased=True, color='darkgrey')

    # put the major ticks at the middle of each cell
    ax.set_xticks(np.arange(dataset.shape[1]) + 0.5, minor=False,)
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
    ax.set_xticklabels(list(self.sequence)[segment[0]-self.start_position:segment[1]-self.start_position+1], fontsize=6.5, fontname="Arial", color='k', minor=False)
    ax.set_yticklabels(temp_kwargs['neworder_aminoacids'], fontsize=6,
                       fontname="Arial", color=ylabel_color, minor=False)

    ax2_label = (segment[1]-segment[0]+1)*['']
    ax2_label[0] = segment[0]
    ax2_label[-1] = segment[1]
    ax2.set_xticklabels(ax2_label, fontsize=7, fontname="Arial", color='k', minor=False)

    # align the labels of the y axis
    for ylabel in ax.get_yticklabels():
        ylabel.set_horizontalalignment('center')

    # remove ticks
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')
    ax2.yaxis.set_ticks_position('none')
    ax2.xaxis.set_ticks_position('none')

    # save file
    _savefile(fig,temp_kwargs)
    
    plt.show()
    return


# ## Mean Plots

# ### Bar graph Enrichment

# In[416]:


def plot_mean(self, show_cartoon=False, **kwargs):
    '''
    Plot in a bargraph the mean enrichment for each residue of the protein. Red for gain of function, blue for loss of function

    Parameters
    ----------
    **kwargs
        color : str, default 'k'
        scale : list, default [-1,1]
            give only 2 values (min,max)
        figsize : list, default (2.5,2)
        title : str, default 'Enrichment Scores'
        outputfilepath : str, default none
        outputfilename : str
        outputformat : str, default png
        savefile : bool, default False

    Returns
    ----------
    Exported png file to desired folder
    '''
    # update kwargs
    temp_kwargs = copy.deepcopy(default_kwargs)
    temp_kwargs.update(kwargs)
    temp_kwargs['figsize'] = kwargs.get('figsize', (3,2.5))
    temp_kwargs['yscale'] = kwargs.get('yscale', (-1,1))
    temp_kwargs['y_label'] = kwargs.get('y_label', r'$∆E^i_x$')

    # load parameters
    parameters_mean()
    
    # make pandas
    df = self.dataframe.groupby('Position',as_index=False).mean()
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
    ax.set_ylabel(temp_kwargs['y_label'], fontsize=10, fontname="Arial", color='k', labelpad=10, rotation=0)
    ax.set_xticks(np.arange(self.start_position, len(df)+self.start_position, 20))
    ax.set_xlabel('Residue', fontsize=10, fontname="Arial", color='k', labelpad=4)
    ax.set_xlim(self.start_position-0.1, len(df)+self.start_position-1+0.1)
    plt.title(temp_kwargs['title'], fontsize=12, fontname='Arial', color='k')
    
    # cartoon
    if show_cartoon:
        _generate_cartoon(self,gs,1,temp_kwargs['cartoon_colors'],
                            bottom_space=-0.78, show_labels=False)
    # Put text labels
    _inputtext(temp_kwargs['text_labels'])
    
    # save file
    _savefile(fig,temp_kwargs)
    plt.show()
    return


def color_data(row):
    if row['Score'] > 0:
        return 'red'
    else:
        return 'blue'


def parameters_mean():
    # normal font
    mpl.rcParams['font.family'] = 'sans-serif'
    mpl.rcParams['font.sans-serif'] = ['Arial']

    # math font
    mpl.rcParams['mathtext.fontset'] = 'custom'
    mpl.rcParams['mathtext.rm'] = 'Arial'
    mpl.rcParams['svg.fonttype'] = 'none'

    # Parameters for all graphs
    mpl.rcParams['xtick.labelsize'] = 9
    mpl.rcParams['ytick.labelsize'] = 9

    return


# ### Compare two proteins

# In[82]:


def plot_meandifferential(self, obj2, show_cartoon=False,**kwargs):
    '''
    Plot the mean difference between two experiments

    Parameters
    ----------
    **kwargs
        color : str, default 'k'
        scale : list, default [-1,1]
            give only 2 values (min,max)
        figsize : list, default (2.5,2)
        title : str, default 'Enrichment Scores'
        outputfilepath : str, default none
        outputfilename : str
        outputformat : str, default png
        savefile : bool, default False

    Returns
    ----------
    Exported png file to desired folder
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

    plt.show()
    return


# ### Bar graph Counts

# In[417]:


def plot_meancounts(self, positions, counts, show_cartoon=False, **kwargs):
    '''
    Plot in a bargraph the mean counts for each residue of the protein.
    Parameters
    ----------
    self
    positions : list, x coordinates
    counts : list, y coordinates
    show_cartoon : boolean, default false
    **kwargs
        color : str, default 'k'
        scale : list, default [-1,1]
            give only 2 values (min,max)
        figsize : list, default (2.5,2)
        title : str, default 'Enrichment Scores'
        outputfilepath : str, default none
        outputfilename : str
        outputformat : str, default png
        savefile : bool, default False

    Returns
    ----------
    Exported png file to desired folder
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
    plt.show()
    return

def _inputtext(text_entries):
    '''the user can input text as a variable by manually giving the coordinates'''
    if text_entries:
        for entry in text_entries:
            plt.text(entry[0],entry[1],entry[2])
    return


# ## Scatter

# In[168]:


def plot_scatter(self, obj2, mode='pointmutant', **kwargs):
    '''
    Generate a scatter plot between object and a second object of the same class.

    Parameters
    ----------
    obj2 : object from class Screen
        second dataset you want to compare to
    mode : str, default 'pointmutant'. Alternative mean of each position
    **kwargs
        title : str, default 'Title'
        start_position : int, default 2
        outputfilepath : str, default ''
        outputfilename : str, default ''
        outputformat : str, default 'png'
        savefile : boolean, default False
        x_label : str, default 'x_label'
        y_label : str, default 'y_label'
        tick_spacing : float, default 0.5
    Returns
    ----------
    Exported png file to desired folder
    '''

    # update kwargs
    temp_kwargs = copy.deepcopy(default_kwargs)
    temp_kwargs.update(kwargs)
    temp_kwargs['figsize'] = kwargs.get('figsize', (2,2))
    temp_kwargs['xscale'] = kwargs.get('xscale', (-2,2))
    temp_kwargs['yscale'] = kwargs.get('yscale', (-2,2))
    
    # Chose mode:
    if mode=='pointmutant':
        df = _process_bypointmutant(self,obj2)
    else: 
        df = _process_meanresidue(self,obj2)

    # create figure
    fig, ax = plt.subplots(figsize=temp_kwargs['figsize'])

    ###import parameters
    parameters()

    # Scatter data points
    plt.scatter(df['dataset_1'], df['dataset_2'], c='k', s=8, alpha=0.5, rasterized=True, label='_nolegend_')

    # Titles
    plt.title(temp_kwargs['title'], fontsize=12, fontname='Arial', color='k', pad=8)
    plt.ylabel(temp_kwargs['y_label'], fontsize=10, fontname="Arial", color='k', labelpad=0)
    plt.xlabel(temp_kwargs['x_label'], fontsize=10, fontname="Arial", color='k')


    ###correlation and R2
    slope, intercept, r_value, p_value, std_err = stats.linregress(df['dataset_1'], df['dataset_2'])
    R2 = str(round(r_value**2, 2))
    legend_label = "$R^2$ = {}".format(R2)
    # fit and graph line
    fit = np.polyfit(df['dataset_1'], df['dataset_2'], 1)
    plt.plot(np.unique(df['dataset_1']), np.poly1d(fit)(
        np.unique(df['dataset_1'])), color='r', linewidth=1, label = legend_label)
    plt.grid()

    # other graph parameters
    plt.xlim(temp_kwargs['xscale'])
    plt.ylim(temp_kwargs['yscale'])
    ax.xaxis.set_major_locator(ticker.MultipleLocator(temp_kwargs['tick_spacing']))  
    ax.yaxis.set_major_locator(ticker.MultipleLocator(temp_kwargs['tick_spacing']))  
    plt.gca().set_aspect('equal', adjustable='box')
    plt.draw()
    
    # Legend
    plt.legend(loc='upper left', handlelength=0, 
               handletextpad=0, frameon=False, fontsize=10)

    # save file
    _savefile(fig,temp_kwargs)
    plt.show()
    return

def _process_bypointmutant(self,obj):
    # truncate so both datasets have same length and delete stop codons
    minlength = min(len(self.dataset[0]), len(obj.dataset[0]))
    dataset_1 = self.dataset[0:20, 0:minlength].copy()
    dataset_2 = obj.dataset[0:20, 0:minlength].copy()

    # convert to dataframe and eliminate Nans
    df = pd.DataFrame()
    df['dataset_1'] = dataset_1.ravel()
    df['dataset_2'] = dataset_2.ravel()
    df.dropna(how='any', inplace=True)
    return df

def _process_meanresidue(self,obj):
    # truncate so both datasets have same length and delete stop codons
    dataset_1 = self.dataframe.groupby(['Position'],as_index=False).mean()
    dataset_2 = obj.dataframe.groupby(['Position'],as_index=False).mean()
    minlength = min(len(dataset_1), len(dataset_2))

    # convert to dataframe and eliminate Nans
    df = pd.DataFrame()
    df['dataset_1'] = list(dataset_1['Score'])[0:minlength]
    df['dataset_2'] = list(dataset_2['Score'])[0:minlength]
    df['Position'] = list(dataset_1['Position'])[0:minlength]
    df['d1 - d2'] = df['dataset_1'] - df['dataset_2']
    df.dropna(how='any', inplace=True)
    return df


# ## SNV

# ### Plot Histogram

# In[13]:


def plot_hist(self, population='All', bins=40, loc='upper left', **kwargs):
    '''
    Generate a histogram plot.

    Parameters
    ----------
    population : str, default 'All'. Other options are 'SNV' and 'nonSNV'
    **kwargs
        color : str, default 'k'
        scale : list, default [-4,2]
        figsize : default (2.5,2)
        title : str, default 'Title'
        start_position : int, default 2
        outputfilepath : str, default ''
        outputfilename : str, default ''
        outputformat : str, default 'png'
        savefile : boolean, default False

    Returns
    ----------
    Exported png file to desired folder

    '''
    # update kwargs
    temp_kwargs = copy.deepcopy(default_kwargs)
    temp_kwargs.update(kwargs)
    temp_kwargs['figsize'] = kwargs.get('figsize', (2,2))
    temp_kwargs['yscale'] = kwargs.get('yscale', (0,2))
    temp_kwargs['xscale'] = kwargs.get('xscale', (-2,2))

    # Select case input data
    df = self.dataframe['Score']
    if population == 'SNV':
        df = self.dataframe_SNV['Score']
    elif population == 'nonSNV':
        df = self.dataframe_nonSNV['Score']

    # create figure
    fig = plt.figure(figsize=temp_kwargs['figsize'])

    # Import parameters
    parameters()

    # plot figure
    plt.hist(df, density=True, bins=bins, color='k')

    # axes labels and title
    plt.xlabel(r'$∆E^i_x$' if temp_kwargs['x_label'] == 'x_label' else temp_kwargs['x_label'], fontsize=10, fontname="Arial", color='k', labelpad=0) 
    plt.ylabel('Probability density', fontsize=10, fontname="Arial", color='k', labelpad=3)
    plt.title(temp_kwargs['title'], fontsize=10, fontname='Arial', color='k')

    # axes limits. spacer will be 1 or the
    plt.xlim(temp_kwargs['xscale'])
    plt.xticks(np.arange(temp_kwargs['xscale'][0], temp_kwargs['xscale'][1]+temp_kwargs['tick_spacing'], temp_kwargs['tick_spacing']))
    plt.ylim(temp_kwargs['yscale'])
    plt.grid()

    # save file
    _savefile(fig,temp_kwargs)
    plt.show()
    return


# ### Internal SNV

# Bring the SNV code I made in cancer script

# In[14]:


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
    NonSNV = pd.concat([SNV, df], sort=False)[['Position','Variant', 'Score']]
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
    df = df[['Position','Variant', 'Score']].copy()

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
    aminoacids = list('FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG')

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
    df['Codons_'+namecolumn] = df.apply(lambda x: _aatocodons(x[namecolumn]), axis=1)

    return df


# ## Miniheatmap

# ### Mean substitution heatmap

# In[15]:


def plot_miniheatmap(self, offset=0, **kwargs):
    '''
    Generate a miniheatmap plot enrichment scores of mutagenesis selection assays.

    Parameters
    ----------
    self
    offset : int, if you want to study effects of a residue when is behind or in front of another residue
         offset of 1 means that you evaluate the effect of following residue n+1 on n
         offset of -1 means that you look at the previous residue (n-1 on n)
    **kwargs
        start_position : int, default is 1        
        colormap : object, default bluewhitered
        scale : list, default [-1,1]
            give only 3 values (min,midpoint,max)
        title : str, default 'Enrichment Scores'
        outputfilepath : str, default none
        outputformat : str, default png
        savefile : bool, default False

    Returns
    ----------
    Exported png file to desired folder
    '''

    # load font parameters
    font_parameters()

    # update kwargs
    temp_kwargs = copy.deepcopy(default_kwargs)
    temp_kwargs.update(kwargs)

    # load labels
    temp_kwargs['color_sequencelabels'] = _labels(self.start_position)[0]
    temp_kwargs['number_sequencelabels'] = _labels(self.start_position)[1]
    
    # do offset if appropriate 
    dataframe_stopcodons= _transform_dataset_offset(self,offset)
    
    # calculate condensed heatmap
    dataset = condense_heatmap(dataframe_stopcodons, temp_kwargs['neworder_aminoacids'])

    _plot_miniheatmap(dataset,temp_kwargs)
    
    return


def condense_heatmap(df, new_order):
    '''
    Converts the np.array with stored enrichment scores into the condensed heatmap
    '''
    # Convert dataset to df
    df = df.copy()
    df.drop(['Position'], axis=1, inplace=True)

    # Group by sequence and aminoacid, and then pivot table
    df_grouped = df.groupby(['Sequence', 'Aminoacid']).mean()
    df_pivoted = df_grouped.pivot_table(values='Score',
                                        index='Aminoacid',  columns='Sequence')
    df_pivoted.reset_index(drop=False, inplace=True)

    # Sort in y axis desired order
    df_pivoted['Aminoacid'] = pd.Categorical(df_pivoted['Aminoacid'], new_order)
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
        trimmedsequence = sequence[start_position-1+offset:len(dataset[0])+start_position-1+offset]
    else:
        sequence = 'X'*(np.absolute(offset))+sequence
        trimmedsequence = sequence[start_position-1:len(dataset[0])+start_position-1]

    return trimmedsequence


def _transform_dataset_offset(self, offset,stopcodons=True):
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
    #Generate new offset sequence
    df['Sequence'] = np.ravel([[aa]*aa_number for aa in offset_sequence])
    
    # Drop rows with X
    df.drop(df.index[df['Sequence'] == 'X'], inplace = True)

    return df


# ### Neighbor residues

# In[16]:


def plot_neighboreffect(self, offset=1, **kwargs):
   '''
   Generate a miniheatmap plot telling you the effect of having a residue in front or behind

   Parameters
   ----------
   self
   offset : int, if you want to study effects of a residue when is behind or in front of another residue
        offset of 1 means that you evaluate the effect of following residue n+1 on n
        offset of -1 means that you look at the previous residue (n-1 on n)
   **kwargs

   Returns
   ----------
   Exported png file to desired folder
   '''
   # load font parameters
   font_parameters()

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
   
   plt.show()
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
   mean = condense_heatmap(self.dataframe, aalist)
   
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

# In[20]:


def plot_correlation(self, **kwargs):
    '''
    Generate a miniheatmap plot enrichment scores of mutagenesis selection assays.

    Parameters
    ----------
    **kwargs
        start_position : int, default is 1        
        colormap : object, default bluewhitered
        scale : list, default [-1,0,1]
            give only 3 values (min,midpoint,max)
        title : str, default 'Enrichment Scores'
        outputfilepath : str, default none
        outputformat : str, default png
        savefile : bool, default False

    Returns
    ----------
    Exported png file to desired folder
    '''

    # load font parameters
    font_parameters()

    # update kwargs
    temp_kwargs = copy.deepcopy(default_kwargs)
    temp_kwargs.update(kwargs)

    # load labels
    temp_kwargs['color_sequencelabels'] = _labels(self.start_position)[0]
    temp_kwargs['number_sequencelabels'] = _labels(self.start_position)[1]

    # calculate correlation heatmap
    dataset = _calculate_correlation(self.dataframe_stopcodons, temp_kwargs['neworder_aminoacids'])

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
    cb.ax.set_yticklabels(cb.ax.get_yticklabels(), fontsize=7, fontname="Arial", color='k')
    cb.update_ticks()
    plt.text(len(dataset.columns)+1.2*coeff, len(dataset.columns)/2.5, 'R', horizontalalignment='center', fontsize=7, fontname="Arial", color='k')

    # for putting title on graph
    plt.title(temp_kwargs['title'], horizontalalignment='center',
              fontname="Arial", fontsize=10, pad=10)

    # save file
    _savefile(fig,temp_kwargs)
    
    plt.show()
    return


def _calculate_correlation(df, order_aminoacids):

    dataset = df.copy()
    dataset = dataset.pivot_table(values='Score', index='Position',  columns='Aminoacid')
    dataset = dataset.corr()
    dataset = dataset.reindex(index=order_aminoacids)[order_aminoacids]

    return dataset

def _calculate_correlation_byresidue(df):

    dataset = df.copy()
    dataset = dataset.pivot_table(values='Score', index='Position',  columns='Aminoacid')
    dataset = dataset.T.corr()

    return dataset


# ### Mean correlation

# In[21]:


def plot_meancorrelation(self, **kwargs):
    '''
    Genereates a bar plot of the mean correlation of each position with the rest'

    Parameters
    -----------
    **kwargs

    Returns
    ----------
    Exported png file to desired folder

    '''
    # Load parameters
    parameters()

    # Update kwargs
    temp_kwargs = copy.deepcopy(default_kwargs)
    temp_kwargs.update(kwargs)
    temp_kwargs['figsize'] = kwargs.get('figsize', (3,2))
    temp_kwargs['yscale'] = kwargs.get('yscale', (0,1))
    
    # Get data
    if '*' in temp_kwargs['neworder_aminoacids']: temp_kwargs['neworder_aminoacids'].remove('*')
    df = _calculate_correlation(self.dataframe, temp_kwargs['neworder_aminoacids']).mean()**2

    # Make figure
    fig, ax = plt.subplots(figsize=temp_kwargs['figsize'])
    ticks = np.arange(0, len(df))  # label locations
    width = 0.5
    labels = temp_kwargs['neworder_aminoacids']
    # Plot figure
    ax.bar(ticks, df, width, color='blue', ec='k',)

    # graph parameters
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels, fontsize=9, fontname="Arial", color='k', minor=False, rotation=0)
    ax.set_ylabel(r'$R^2$', fontsize=10, fontname="Arial", color='k', labelpad=12, rotation=0)
    ax.set_ylim(temp_kwargs['yscale'])
    plt.title(temp_kwargs['title'], horizontalalignment='center',
              fontname="Arial", fontsize=10, pad=5)

    # save file
    _savefile(fig,temp_kwargs)
    plt.show()
    return


# ### Group correlation

# In[22]:


def plot_topcorrelators(self, r2, groups=['DEHKR','QN','CASTG','ILMV','WYF'], output=False, **kwargs):
    '''Use logomaker to plot the most frequent residues'''
    # Update kwargs
    temp_kwargs = copy.deepcopy(default_kwargs)
    temp_kwargs.update(kwargs)
    
    # Apply parameters
    parameters()
    
    # Calculate
    if '*' in temp_kwargs['neworder_aminoacids']: temp_kwargs['neworder_aminoacids'].remove('*')
    df = _calculate_substitution_correlations(self,temp_kwargs['neworder_aminoacids'],groups='')
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

    # save file
    _savefile(fig,temp_kwargs)
    
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
    for combination in replacement_combinations: # Iterate over a combination
        temp_list = []
        for group,aa_selected in zip(groups,combination): # Iterate over a group of the combination
            for aa_nonselected in group: # Find correlation values from correlation plot
                if aa_nonselected == aa_selected:
                    temp_list.append(1)
                else:
                    temp_list.append(_find_correlation(aa_selected, aa_nonselected,corr_values))
        df[combination] = temp_list # Store in df
    return _polishdf(df)

def _polishdf(df):
    df_mean = df.copy()
    df_mean = df.mean().to_frame()
    df_mean.reset_index(drop=False,inplace=True)
    df_mean.rename(columns={0:'R2'},inplace=True)
    df_mean['Combinations'] = list(df_mean['index'].apply(lambda x: ''.join(x)))
    df_mean.drop(columns = ['index'],inplace=True)
    return df_mean

def _find_correlation(aa1, aa2, corr_values):
    return float(corr_values[aa1].loc[corr_values['Aminoacid']==aa2])


# ## PCA

# In[23]:


def plot_pca(self, mode='aminoacid',dimensions=[0, 1], **kwargs):
    '''
    Genereates a plot of PCA1 vs PCA2

    Parameters
    -----------
    dimensions : list, default [0,1]
        specify which two PCA dimensions to plot. By default PCA1 vs PCA2.
        Also available PCA3 (2) and PCA4 (3)
    **kwargs

    Returns
    ----------
    Exported png file to desired folder

    '''

    # load parameters
    parameters()

    # update kwargs
    temp_kwargs = copy.deepcopy(default_kwargs)
    temp_kwargs.update(kwargs)
    temp_kwargs['figsize'] = kwargs.get('figsize', (2,2))

    # calculate correlation heatmap. Choose mode
    dataset = self.dataframe.copy()
    if mode is 'aminoacid':
        dataset = _calculate_correlation(dataset, temp_kwargs['neworder_aminoacids'][:20])
        if '*' in temp_kwargs['neworder_aminoacids']: temp_kwargs['neworder_aminoacids'].remove('*')
        textlabels = temp_kwargs['neworder_aminoacids']
    elif mode is 'secondary':
        dataset = _calculate_correlation_bysecondary(dataset,self.secondary_dup)
        textlabels = list(dataset.columns)
    elif mode is 'residue':
        dataset = _calculate_correlation_byresidue(dataset)
        textlabels = list(dataset.columns)
        # plot using plot_clusters
    dimensionstoplot, variance = _calculate_clusters(dataset, dimensions)

    ###x and y
    x = dimensionstoplot.iloc[:, 0]
    y = dimensionstoplot.iloc[:, 1]

    fig, ax = plt.subplots(figsize=temp_kwargs['figsize'])
    ax.scatter(x, y, s=4, c='k')  # ,cmap='tab10')

    # labels
    plt.xlabel('PCA ' + str(dimensions[0]+1) + ': ' + str(int(variance[dimensions[0]]*100))+'%',
               fontsize=10, labelpad=5, fontweight='normal')
    plt.ylabel('PCA ' + str(dimensions[1]+1) + ': ' + str(int(variance[dimensions[1]]*100))+'%',
               fontsize=10, labelpad=-2, fontweight='normal')

    # label of data points
    texts = _auto_text(x, y, textlabels)
    if mode is not 'residue': adjust_text(texts, autoalign='xy')

    plt.title(temp_kwargs['title'], horizontalalignment='center',
              fontname="Arial", fontsize=10, pad=5)

    # save file
    _savefile(fig,temp_kwargs)
    
    plt.show()
    return

def _auto_text(x, y, textlabels,):
    '''auto anotates text labels'''
    texts = [plt.annotate(textlabels[i],  # this is the text
                      (x[i], y[i]),  # this is the point to label
                      textcoords="offset points",  # how to position the text
                      xytext=(2, 2),  # distance from text to points (x,y)
                      fontsize=8,
                      ha='center')  # horizontal alignment can be left, right or center
         for i in range(len(textlabels))]
    return texts

def _calculate_clusters(dataset, dimensions):
    '''input the dataframe that needs to be correlated, the dimensions, and will calculate PCA descomposition. '''

    # call pca model
    pca = PCA(n_components=6, random_state=554)  # 121

    # fit model to df. use aux function correlation_aminoacids
    model = pca.fit(dataset)

    # create df with PCA data
    df_aa = pd.DataFrame((model.components_).T, columns=['PCA1', 'PCA2', 'PCA3', 'PCA4', 'PCA5', 'PCA6'])

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
    df = df.groupby(['Secondary','Aminoacid'], as_index=False).mean()

    return df

def _calculate_correlation_bysecondary(df, secondary):
    dataset = _grouby_secondary(df, secondary)
    dataset = dataset.pivot_table(values='Score', index='Secondary',  columns='Aminoacid')
    dataset = dataset.T.corr()

    return dataset


# ## Secondary Structure

# In[24]:


def plot_secondary(self, **kwargs):
    '''
    Genereates a plot of data sorted by secondary elements (alpha helixes and beta sheets).

    Parameters
    -----------
    **kwargs

    Returns
    ----------
    Exported png file to desired folder

    '''
    # Load parameters
    parameters()

    # Update kwargs
    temp_kwargs = copy.deepcopy(default_kwargs)
    temp_kwargs.update(kwargs)
    temp_kwargs['figsize'] = kwargs.get('figsize', (2.5,2))
    temp_kwargs['yscale'] = kwargs.get('yscale', (-2,1))
    
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
    ax.set_xticklabels(labels, fontsize=9, fontname="Arial", color='k', minor=False, rotation=0)
    ax.set_ylabel(r'$∆E^i_x$', fontsize=10, fontname="Arial", color='k', labelpad=12, rotation=0)
    ax.set_ylim(temp_kwargs['yscale'])
    plt.title(temp_kwargs['title'], horizontalalignment='center',
              fontname="Arial", fontsize=10, pad=5)
    # plt.tight_layout()

    # save file
    _savefile(fig,temp_kwargs)
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

# In[25]:


def plot_roc(self, df_class=None, **kwargs):
    '''
    Generates ROC AUC plot

    Parameters
    -----------
    self
    df_class: Pandas dataframe
        A dataframe that contains a column of variants labeled 'Variant' with a column labeled 'Class'
        containing the true class of that mutation

    Returns
    --------
    Exported png file to desired folder
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
    temp_kwargs['figsize'] = kwargs.get('figsize', (2.5,2.5))

    ###import parameters
    parameters()

    # create figure
    fig, ax = plt.subplots(figsize=temp_kwargs['figsize'])
    lw = 2
    plt.plot(fpr, tpr, color='k', lw=lw, label='AUC = %0.2f' % auc)
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')

    # Graph limits
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    tick_spacing = 0.2
    ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(tick_spacing))  # Plt ticks
    ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(tick_spacing))  # Plt ticks

    # Axis labels
    plt.title(temp_kwargs['title'], fontsize=12, fontname='Arial', color='k')
    plt.ylabel('True Positive Rate', fontsize=12, fontname="Arial", color='k', labelpad=0)
    plt.xlabel('False Positive Rate', fontsize=12, fontname="Arial", color='k')

    # Legend
    plt.legend(loc='lower right', handlelength=0, handletextpad=0, frameon=False)

    # save file
    _savefile(fig,temp_kwargs)
    plt.show()
    return


def _rocauc(df):
    '''
    Calculate roc rates and auc.

    The input is a dataframe that contains [Variants,Class,Score]
    '''
    fpr, tpr, thresholds = metrics.roc_curve(df['Class'], df['Score'], drop_intermediate=True)
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

def _concattrueposneg(df_tp, df_tn, subset = 'Variant',keep='first'):
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
    Exported png file to desired folder
    '''
    # Concatenate tp and tn datasets
    df_true = pd.concat([df_tp, df_tn],sort=False)

    # Will keep a variant as true positive if found in both datasets (because could be a mistake in gnomAD)
    df_true.drop_duplicates(subset = subset, keep = keep, inplace=True)
    
    return df_true


# ## Cumulative

# In[26]:


def plot_cumulative(self, modes=None,cutoffs=None, hows=None,listdf=None,
                    labels=None,colors_1=None,colors_2=None,**kwargs):
    '''
    Generates a cumulative plot of mutations 

    Parameters
    -----------
    self
    listdf: List of Pandas dataframe containing other elements that you may want to compare to the screen data

    Returns
    --------
    Exported png file to desired folder
    '''

    # update kwargs
    temp_kwargs = copy.deepcopy(default_kwargs)
    temp_kwargs.update(kwargs)
    temp_kwargs['figsize'] = kwargs.get('figsize', (3.5,2.5))

    # import parameters
    parameters()
    
    # create figure
    fig, ax = plt.subplots(figsize=temp_kwargs['figsize'])
    lw = 2
    
    # Get data filtered
    for mode,cutoff,how,color in zip(modes,cutoffs,hows,colors_1):
        df = _filter(self,mode,cutoff,how)
        df_cumulative = _cumsum(df)
        plt.plot(df_cumulative['Position'], df_cumulative['Cumulative Norm'], color=color, lw=lw)

    # Graph limits
    plt.xlim(self.dataframe['Position'].min(), self.dataframe['Position'].max()+20)
    plt.ylim(0, 1.1)

    ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(temp_kwargs['tick_spacing']))  # Plt ticks

    # Axis labels
    plt.title(temp_kwargs['title'], fontsize=12, fontname='Arial', color='k')
    plt.ylabel('Cumulative', fontsize=12, fontname="Arial", color='k', labelpad=5)
    plt.xlabel('Position', fontsize=12, fontname="Arial", color='k', labelpad=0)
    
    # Plot other datasets
    if labels:
        for df,color in zip(listdf,colors_2):
            temp = _cumsum(df)
            plt.plot(temp['Position'], temp['Cumulative Norm'], color=color, lw=lw)
        # Legend
        plt.legend(labels=labels,loc='lower right', handlelength=0.5, handletextpad=0.5, frameon=False, fontsize=10)
    
    # x=y line
    plt.plot([0, df_cumulative['Position'].max()], [0, 1], color='silver', lw=lw, linestyle='--')

    # save file
    _savefile(fig,temp_kwargs)
    plt.show()
    return

def _filter(self,mode,cutoff,how='right'):
    '''Choose either GoF or LoF mutations'''
    
    # Select all, SNV, nonSNV
    if mode == 'all':
        df = self.dataframe
    elif mode == 'SNV':
        df = self.dataframe_SNV
    elif mode == 'nonSNV':
        df = self.dataframe_nonSNV        
    
    # Filter data
    if how == 'right':
        df_filtered = df.loc[df['Score']>cutoff].copy()
    elif how == 'left':
        df_filtered = df.loc[df['Score']<cutoff].copy()
    
    # Add column counts
    df_filtered['Counts'] = 1
    return df_filtered

def _cumsum(df_raw):
    ''' Add column with cumulative sum to df'''
    df = df_raw.copy()
    if 'Variant' in df.columns:
        df.drop_duplicates(subset = 'Variant', keep = 'first', inplace=True)
        df['Position'] = df['Variant'].str.extract('(\d+)').astype(str).astype(int)
        df['Counts'] = 1
    df.sort_values(by=['Position'],inplace=True)
    df = df.groupby('Position',as_index=False).sum()[['Position','Counts']].copy()
    df['Cumulative Norm'] = np.cumsum(df['Counts'])/df['Counts'].sum()
    return df

def plot_box(self, df_input, xcolumn, method = 'absolute mean',**kwargs):
    '''
    Genereates a boxplot. Data needs to be binned prior before using this function. 

    Parameters
    -----------
    method : str, default absolute mean
        others : mean, absolute individual, individual 
    **kwargs

    Returns
    ----------
    Exported png file to desired folder

    '''
    # Load parameters
    parameters()

    # Update kwargs
    temp_kwargs = copy.deepcopy(default_kwargs)
    temp_kwargs.update(kwargs)
    temp_kwargs['figsize'] = kwargs.get('figsize', (2.5,2))
    temp_kwargs['yscale'] = kwargs.get('yscale', (0,2))
    temp_kwargs['xscale'] = kwargs.get('xscale', (0,3))
    
    # Process data
    data = _merge_dataset(self.dataframe,df_input,method)
    
    # Make figure
    fig, ax = plt.subplots(figsize=temp_kwargs['figsize'])
    # Plot data
    ax = sns.boxplot(x=xcolumn, y='Score', data = data,color='white',fliersize=2)
    
    # color edges
    plt.setp(ax.artists, edgecolor = 'k', facecolor='w')
    plt.setp(ax.lines, color='k')
    
    # graph parameters
    plt.title(temp_kwargs['title'], fontsize=10, fontname='Arial', color='k', pad=8)
    plt.ylabel(temp_kwargs['y_label'], fontsize=10, fontname="Arial", color='k', labelpad=0)
    plt.xlabel(temp_kwargs['x_label'], fontsize=10, fontname="Arial", color='k')
    
    # axes limits
    plt.xlim(temp_kwargs['xscale'])
    plt.ylim(temp_kwargs['yscale'])
    plt.grid()
    
    # save file
    _savefile(fig,temp_kwargs)
    plt.show()
    
    return

def _merge_dataset(df,df_input,method):
    if method == 'mean':
        df = df[['Position','Score']].groupby(by='Position').mean()
        df.reset_index(drop=False,inplace=True)
        df_merged = pd.merge(df,df_input,how='left',on=['Position'])
    elif method == 'absolute mean':
        df = df[['Position','Score']].abs().groupby(by='Position').mean()
        df.reset_index(drop=False,inplace=True)
        df_merged = pd.merge(df,df_input,how='left',on=['Position'])
    elif method == 'absolute individual':
        df_merged = pd.merge(df,df_input,how='left',on=['Position','Aminoacid'])
        df_merged['Score'] = df_merged['Score'].abs() 
    else:
        df_merged = pd.merge(df,df_input,how='left',on=['Position','Aminoacid'])
    df_merged.dropna(how='any',inplace=True)
    return df_merged


# ## Map into Pymol

# In[35]:


def plot_pymol(self, pdb, residues=None, **kwargs):
    '''
    Color pymol structure residues. User can specify the residues to color, or can use the mutagenesis data.
    Activating mutations will be colored red and loss of function blue. Neutral mutations in green

    Parameters
    ----------
    pdb : str,
        User should specify the PDB chain in the following format 4G0N_A
    residues : list , default None
        If user decides to pass custom arguments, use the following format
        residues = ['1,2,3,4-10','12-15,23,24,35','48,49,50,52-60'] which are [blue,red,green]
    **kwargs
         start_position : int, default is 2
         gof : int, default is 1
             cutoff for determining gain of function mutations based on mutagenesis data
         lof : int, default is -1
             cutoff for determining loss of function mutations based on mutagenesis data
    Returns
    ----------
    Opens pymol and plots desired pdb structure
    '''
    # update kwargs
    temp_kwargs = copy.deepcopy(default_kwargs)
    temp_kwargs.update(kwargs)

    # Calculate residues only if they are not given by the user
    if residues is None:
        # remove *
        df = self.dataframe.copy()
        #df = convert_to_df(self.dataset,self.sequence,self.aminoacids,self.start_position)
        df = df[df['Aminoacid'] != '*']
        residues = pymol_fitness(df, temp_kwargs['gof'], temp_kwargs['lof'])

    # Start Pymol
    if not pymol._process_is_running():
        pymol.start()

    # Fetch structure
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
    light_parameters()

    pymol.deselect()
    # pymol.quit()
    return

# Convert fitness scores into pymol residues


def pymol_fitness(df, gof, lof):
    '''You input the dataframe. Removes stop codons. 
    Returns the positions that are going to be colored blue,red and white'''

    # adjust loss of funtion values so they don't move the mean too much
    df_grouped = df.groupby(['Position'], as_index=True).mean()
    df_grouped.reset_index(drop=False, inplace=True)

    # Color of mutations
    blue_mutations = df_grouped[df_grouped['Score'] < lof]
    red_mutations = df_grouped[df_grouped['Score'] > gof]
    white_mutations = df_grouped[df_grouped['Score'].between(lof, gof, inclusive=True)]

    # Pymol Format
    blue_pymol = array_to_pymol(blue_mutations['Position'])
    red_pymol = array_to_pymol(red_mutations['Position'])
    white_pymol = array_to_pymol(white_mutations['Position'])

    residues = [blue_pymol, red_pymol, white_pymol]

    for i, residue in enumerate(residues):
        if residue == '':
            residues[i] = '0'

    return residues


def array_to_pymol(array):
    '''Input an array with positions of aminoacids, return it in pymol format'''
    pymol = ''
    for aminoacid in array:
        pymol += str(aminoacid)+'+'

    # delete last '+'
    pymol = pymol[:-1]
    return pymol


def light_parameters():
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


# ## Internal Functions

# In[36]:


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
    df['Position'] = np.ravel([[i]*len(aminoacids) for i in range(start_position, len(dataset[0])+start_position)])
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
    trimmedsequence = sequence[start_position-1:len(dataset[0])+start_position-1]

    return trimmedsequence


def _transform_secondary(dataset, secondary, start_position,aminoacids):
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
    trimmedsecondary = secondary_list[start_position-1:len(dataset[0])+start_position-1]

    # Multiply each element by number of aminoacids. not use stop codon
    aminoacids = list(np.copy(aminoacids))
    if '*' in aminoacids: aminoacids.remove('*')
    secondary_dup = [x for item in trimmedsecondary for x in itertools.repeat(item, len(aminoacids))]

    return trimmedsecondary,secondary_dup


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
    df['Sequence'] = np.ravel([[i]*len(aminoacids) for i in sequence[:len(dataset[0])]])
    df['Score'] = np.copy(dataset.T).ravel()
    return df


def _df_rearrange(df,new_order,values = 'Score'):
    '''
    convert a df into a numpy array for mutagenesis data. 
    Allows the option of keeping NaN scores

    Returns copy
    '''
    dfcopy = df.copy()
    df_pivoted = dfcopy.pivot_table(values=values, index='Aminoacid',
                                    columns=['Position'],dropna=False)
    df_reindexed = df_pivoted.reindex(index=list(new_order))

    return df_reindexed


def _common(a, b):
    '''
    return common elements of two lists
    '''
    c = [value for value in a if value in b]
    return c

def _transpose(df,values = 'Score'):
    '''
    convert a df into a numpy array for mutagenesis data

    Returns copy
    '''
    df = df.pivot_table(values=values, index='Aminoacid',columns=['Position']).T
    return  df

def _select_aa(df,selection,values = 'Score'):
    '''returns copy'''
    df = _transpose(df.copy(),values) 
    
    df =  df[selection].T

    return df

def _savefile(fig,temp_kwargs):
    '''Save file function'''
    if temp_kwargs['savefile'] is True:
        filename = temp_kwargs['outputfilepath'] +             temp_kwargs['outputfilename']+"."+temp_kwargs['outputformat']
        fig.savefig(filename, format=temp_kwargs['outputformat'],
                    bbox_inches='tight', dpi=temp_kwargs['dpi'], transparent=True)    
    return


# # Kwargs

# In[397]:


default_kwargs = {'colormap': generatecolormap(),
                  'colorbar_scale': [-1, 1],
                  'color': 'k',
                  'title': 'Title',
                  'x_label': 'x_label',
                  'y_label': 'y_label',
                  'tick_spacing' : 1,
                  'outputfilepath': '',
                  'outputfilename': '',
                  'outputformat': 'png',
                  'dpi' : 600,
                  'aminoacids' : list('ACDEFGHIKLMNPQRSTVWY*'),
                  'neworder_aminoacids': list('DEKHRGNQASTPCVYMILFW*'),
                  'savefile': False,
                  'gof': 1,
                  'lof': -1,
                  'cartoon_colors' : ['lightgreen','lavender','k'],
                  'text_labels' : None,
                  
                 }


# # Define Class

# In[418]:


class Screen:
    '''
    Class to define screen objects.

    Parameters
    -----------
    dataset : array
        2D matrix containing the enrichment score of a selection assay.
        Nrows:21, Ncolumns: variable
        Make sure stop codon is in the last row
    sequence : str
        protein sequence in 1 letter code format
    **kwargs : dict
        aminoacids : str, default ACDEFGHIKLMNPQRSTVWY*
            order of amino acids in which the input dataset is
        start_position : int, default 2
            number of first residue you want to count from the sequence. 
            The last residue will be calculated based on the length of the input dataset 
        secondary : list, default None
            the format is the name of the secondary structure multiplied by the residue length of that motif
            example : [['β1']*(8),['L1']*(7),['α1']*(9)]
        roc_df: Pandas dataframe
            A dataframe that contains a column of variants labeled 'Variant' with a column labeled 'Class'
            containing the true class of that mutation
        fillna : int, default 0
            How to replace NaN values
    '''

    def __init__(self, dataset, sequence, **kwargs):
        self.dataset = dataset
        self.dataset_nonstop = dataset[:,0:20]
        self.aminoacids = list(kwargs.get('aminoacids', 'ACDEFGHIKLMNPQRSTVWY*'))
        self.start_position = int(kwargs.get('start_position', 2))
        self.sequence_raw = sequence
        self.sequence = _transform_sequence(self.dataset, sequence, self.start_position)
        self.dataframe_stopcodons, self.dataframe = _transform_dataset(
            dataset, self.sequence, self.aminoacids, self.start_position, kwargs.get('fillna', 0))
        self.dataframe_SNV = _select_SNV(self.dataframe)
        self.dataframe_nonSNV = _select_nonSNV(self.dataframe)
        
        # Optional parameters
        self.roc_df = kwargs.get('roc_df', None)
        self.secondary = kwargs.get('secondary', None)
        if self.secondary is not None:
            self.secondary, self.secondary_dup = _transform_secondary (self.dataset, self.secondary, self.start_position, self.aminoacids)

    # Associated functions
    kernel = plot_kernel
    heatmap = plot_heatmap
    heatmap_selection = plot_heatmap_selection
    subset = plot_heatmapsubset
    mean = plot_mean
    meancounts = plot_meancounts
    differential = plot_meandifferential
    scatter = plot_scatter
    histogram = plot_hist
    miniheatmap = plot_miniheatmap
    neighbor = plot_neighboreffect
    neighboreffect = plot_neighboreffect
    correlation = plot_correlation
    meancorrelation = plot_meancorrelation
    pca = plot_pca
    secondary_mean = plot_secondary
    roc = plot_roc
    boxplot = plot_box
    pymol = plot_pymol

