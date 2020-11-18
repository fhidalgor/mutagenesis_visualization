#!/usr/bin/env python
# coding: utf-8

# # Test of code_utils
# 

# In[1]:


try:
    import import_notebook
except ModuleNotFoundError:
    pass

import pandas as pd
import numpy as np
import code_utils


# ## Internal functions

# ### Parse dataset

# In[ ]:


def test_common():
    # lists
    list_a = list('abcde')
    list_b = list('efghi')
    list_c = list('jklmn')
    list_d = list('mklj')

    # error message
    error_message = '_common not returning the common elements of two lists'

    # assert
    assert code_utils._common(list_a, list_b) == ['e'], error_message
    assert code_utils._common(list_a, list_c) == [], error_message
    assert code_utils._common(list_c, list_d) == list('jklm'), error_message


# ### SNV internal

# In[2]:


def test_aminoacids_snv():
    '''testing full capabilities of function'''
    
    # Create dict with codons
    codontable = code_utils._dict_codontoaa()
    
    # Create input arguments
    pairs = [['F', 'L'], ['I', 'M'], ['T', 'A'], ['S', 'R'],
         ['F', 'P'], ['I', 'G'], ['T', 'L'], ['S', 'H'], ['A','A']]
    
    # Expected answer
    expected_answer = [True]*4 + [False]*5
    
    # Calculate answer and assert
    for pair, exp_ans in zip(pairs, expected_answer):
        calculated_answer = code_utils._aminoacids_snv(
            pair[0], pair[1], codontable, same_aa_SNV=False)
        assert exp_ans == calculated_answer, 'Error in determining SNV'
    
    # Now change the same_aa_SNV parameter
    expected_answer = code_utils._aminoacids_snv('A', 'A', codontable, same_aa_SNV=True)
    
    assert  expected_answer == True, 'Error in determining SNV when the two amino acids are the same'


# In[19]:


def test_codons_snv():
    '''testing full capabilities of function'''
    
    # Create input arguments
    pairs = [['AAA', 'AAT'], ['ACA', 'ACT'], ['ATA', 'TTA'], ['CCC', 'CCT'],
         ['AAA', 'ACC'], ['CAA', 'CCC'], ['ATA', 'TTT'], ['CCC', 'AAA'], ['AAA','AAA']]
    
    # Expected answer
    expected_answer = [True]*4 + [False]*5
    
    # Calculate answer and assert
    for pair, exp_ans in zip(pairs, expected_answer):
        calculated_answer = code_utils._codons_pointmutants(
            pair[0], pair[1], same_codon_SNV=False)
        assert exp_ans == calculated_answer, 'Error in determining SNV'
    
    # Now change the same_aa_SNV parameter
    expected_answer = code_utils._codons_pointmutants('CAA', 'CAA', same_codon_SNV=True)
    
    assert  expected_answer == True, 'Error in determining SNV when the two codons are the same'


# ### Scatter Internal

# In[ ]:


def test_process_bypointmutant():
    '''testing that output type is a dataframe'''
    # Create mock objects
    self = type('', (), {})()
    obj = type('', (), {})()

    # Create dataframes as attributes of the objects
    self.dataframe = pd.DataFrame(np.array([[1, 2], [4, 5], [7, 8]]),
                                  columns=['Score_NaN', 'Variant'])
    obj.dataframe = pd.DataFrame(np.array([[7, 8], [9, 0]]),
                                 columns=['Score_NaN', 'Variant'])

    # Call the function we are testing
    df = code_utils._process_bypointmutant(self, obj)

    # Assert
    assert len(df) == 2, 'truncation of longer dataset is not working properly'


# In[ ]:


def test_process_meanresidue():
    '''testing full capabilities of function'''
    # Create mock objects
    self = type('', (), {})()
    obj = type('', (), {})()

    # Create dataframes as attributes of the objects
    self.dataframe = pd.DataFrame(np.array([[1, 2], [1, 6], [2, 8], [2, 4]]),
                                  columns=['Position', 'Score'])
    obj.dataframe = pd.DataFrame(np.array([[7, 8], [9, 0], [1, 6]]),
                                 columns=['Position', 'Score'])
    expected_answer = pd.DataFrame(np.array([[4, 6, 1, -2], [6, 8, 2, -2]]),
                                   columns=['dataset_1', 'dataset_2', 'Position', 'd1 - d2'])
    # Call the function we are testing
    df = code_utils._process_meanresidue(self, obj)

    # Assert
    assert df.equals(expected_answer), 'error in _process_meanresidue'


# In[ ]:


def test_color_data():
    '''testing full capabilities of function'''
    df = pd.DataFrame()
    df['Score'] = [1, 2, 3, 0, -1, -2, -3]
    df['Expected_Answer'] = ['red']*3+['blue']*4
    df['Calculated_answer'] = [code_utils._color_data(
        df.loc[i], 'red', 'blue') for i in range(0, len(df['Score']))]
    assert (df['Expected_Answer'] == df['Calculated_answer']
            ).all(), 'error when assigning a color'


# ## To manipulate reads

# In[ ]:


def test_translate_codons():
    '''testing full capabilities of function'''
    list_codons = ['AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 'AGC', 'AGG', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT', 'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC', 'CGG', 'CGT', 'CTA', 'CTC', 'CTG',
                   'CTT', 'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT', 'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT', 'TAA', 'TAC', 'TAG', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT', 'TGA', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT']
    list_aminoacids = ['K', 'N', 'K', 'N', 'T', 'T', 'T', 'T', 'R', 'S', 'R', 'S', 'I', 'I', 'M', 'I', 'Q', 'H', 'Q', 'H', 'P', 'P', 'P', 'P', 'R', 'R', 'R', 'R', 'L', 'L',
                       'L', 'L', 'E', 'D', 'E', 'D', 'A', 'A', 'A', 'A', 'G', 'G', 'G', 'G', 'V', 'V', 'V', 'V', '*', 'Y', '*', 'Y', 'S', 'S', 'S', 'S', '*', 'C', 'W', 'C', 'L', 'F', 'L', 'F']
    df = pd.DataFrame(index=list_codons)
    translation = code_utils._translate_codons(df)
    assert (translation ==
            list_aminoacids), 'error when translating the codons of the dataframe index'


# In[34]:


def test_is_DNA():
    '''testing full capabilities of function'''
    df = pd.DataFrame(index=['A', 'C', 'T', 'G', 'P', 'L'])
    df2 = pd.DataFrame(index=['ATA', 'CAT', 'TGG', 'TGT'])
    assert (code_utils._is_DNA(
        df) == False), 'error determining if the index of the dataframe contains DNA'
    assert (code_utils._is_DNA(df2) ==
            True), 'error determining if the index of the dataframe contains DNA'

