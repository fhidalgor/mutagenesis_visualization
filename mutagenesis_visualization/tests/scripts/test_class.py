#!/usr/bin/env python
# coding: utf-8

# # Test of code_class
# 

# In[2]:


import pandas as pd
import numpy as np

try:
    import mutagenesis_visualization.main.scripts.code_class as code_class
except ModuleNotFoundError:
    import import_notebook
    import os
    directory = os.getcwd()
    new_directory = directory.replace('tests','main')
    os.chdir(new_directory)
    import code_class as code_class
    os.chdir(directory)


# ## Counts

# In[12]:


def test_Counts():
    
    # fake dataframe
    df = pd.DataFrame(
        np.random.rand(21,10)*100, index=list('ACDEFGHIKLMNPQRSTVWY*')
    )
    
    # Define aux function
    def _test_Class_output(df):
        error = False
        try:
            code_class.Counts(df)
        except:
            error = True
        return error
    
    assert _test_Class_output(df)==False, 'Error when generating class Counts'


# ## Screen

# In[18]:


def test_Screen():
    
    # fake dataframe
    df = pd.DataFrame(
        np.random.rand(21,10)
    )
    sequence = 'MTEYKRVVVLL'
    # Define aux function
    def _test_Screen_output(df, sequence):
        error = False
        try:
            code_class.Screen(df, sequence)
        except:
            error = True
        return error
    
    assert _test_Screen_output(df,sequence)==False, 'Error when generating class Screen'

