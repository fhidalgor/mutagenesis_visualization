#!/usr/bin/env python
# coding: utf-8

# # Import Modules

# In[1]:


# Regular libraries
import numpy as np
import pandas as pd
import itertools
from collections import defaultdict
from Bio.Seq import Seq
from pathlib import Path
from typing import Union


# # Internal Functions

# ## Parse dataset

# In[2]:


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
    '''DEPRECATED
    Save file function'''
    if temp_kwargs['savefile'] is True:
        filename = temp_kwargs['outputfilepath'] +             temp_kwargs['outputfilename']+"."+temp_kwargs['outputformat']
        fig.savefig(filename, format=temp_kwargs['outputformat'],
                    bbox_inches='tight', dpi=temp_kwargs['dpi'], transparent=True)
    return

def _save_work(fig, output_file, temp_kwargs):
    '''Save file function using pathlib'''
    if output_file:
        fig.savefig(Path(output_file), format=Path(output_file).suffix.strip('.'),
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


# ## SNV internal

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
    NonSNV = pd.concat([SNV, df], sort=False)[['Position', 'Variant', 'Score','Score_NaN']]
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
    df = df[['Position', 'Variant', 'Score','Score_NaN']].copy()

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


# ## Scatter Internal

# In[ ]:


def _process_bypointmutant(self, obj):
    # truncate so both datasets have same length and delete stop codons
    minlength = min(len(self.dataframe), len(obj.dataframe))
    df = pd.DataFrame()
    df['dataset_1'] = list(self.dataframe['Score_NaN'])[:minlength]
    df['dataset_2'] = list(obj.dataframe['Score_NaN'])[:minlength]
    df['Variant'] = list(self.dataframe['Variant'])[:minlength]
    
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


# In[ ]:


def _color_data(row, color_gof, color_lof):
    if row['Score'] > 0:
        return color_gof
    else:
        return color_lof


# # To manipulate reads

# In[1]:


def _translate_codons(df):
    '''Translate the index of the df from codons to AA'''
    list_aa = [str(Seq(codon).translate()) for codon in list(df.index)]
    return list_aa


def _is_DNA(df):
    '''Check if the index of the dataframe are the DNA codons'''
    aminoacids = 'DEFHIKLMNPQRSVWY*'
    for aa in aminoacids:
        if aa in ''.join(list(df.index)):
            return False
        return True

