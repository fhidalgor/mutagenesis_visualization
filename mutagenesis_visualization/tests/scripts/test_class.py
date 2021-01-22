#!/usr/bin/env python
# coding: utf-8

# # Test of code_class
# 

# In[3]:


import pandas as pd
import numpy as np

try:
    import mutagenesis_visualization.main.scripts.code_class as code_class
except ModuleNotFoundError:
    import import_notebook
    import os
    directory = os.getcwd()
    new_directory = directory.replace('tests', 'main')
    os.chdir(new_directory)
    import code_class as code_class
    os.chdir(directory)


# ## Counts

# In[20]:


def test_Counts():

    # fake dataframe
    df = pd.DataFrame(
        np.random.rand(21, 10) * 100, index=list('ACDEFGHIKLMNPQRSTVWY*')
    )
    aminoacids = list('ACDEFGHIKLMNPQRSTVWY*')

    # Define aux function
    def _test_Counts_output(parameters):
        error = False
        try:
            code_class.Counts(**parameters)
        except:
            error = True
        return error

#    assert _test_Class_output(df) == False, 'Error when generating class Counts'
    list_params = [{'df': df},
                   {'df': df, 'start_position': 1},
                   {'df': df, 'aminoacids': aminoacids}]    
        # Assert
    for parameters in list_params:
        assert _test_Counts_output(
            parameters
        ) == False, "class Counts failed with {} parameters".format(parameters)


# ## Screen

# In[6]:


def test_Screen():

    # fake dataframe
    df = pd.DataFrame(np.random.rand(21, 10))
    sequence = 'MTEYKRVVVLL'
    secondary = ['β1'] * len(sequence)

    # Define aux function
    def _test_Screen_output(parameters):
        error = False
        try:
            code_class.Screen(**parameters)
        except:
            error = True
        return error

    list_params = [{'dataset': df, 'sequence': sequence},
                   {'dataset': df, 'sequence': sequence}]

    # Assert
    for parameters in list_params:
        assert _test_Screen_output(
            parameters
        ) == False, "class Screen failed with {} parameters".format(parameters)


# In[7]:


test_Screen()


# In[ ]:




