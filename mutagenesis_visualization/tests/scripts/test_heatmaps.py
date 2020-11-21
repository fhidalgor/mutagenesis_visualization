#!/usr/bin/env python
# coding: utf-8

# # Test of code_heatmaps
# 

# ## Import modules

# In[1]:


# kwargs and parameters
try:
    import import_notebook
except ModuleNotFoundError:
    pass

import pandas as pd
try:
    from mutagenesis_visualization.main.scripts.code_heatmaps import (
        _hierarchical_sort, _helix, _labels, _sheet, _loop)
except ModuleNotFoundError:
    import import_notebook
    import os
    directory = os.getcwd()
    new_directory = directory.replace('tests','main')
    os.chdir(new_directory)
    from code_heatmaps import (
        _hierarchical_sort, _helix, _labels, _sheet, _loop)
    os.chdir(directory)


# ## Test main functions

# ## Test aux functions

# In[2]:


def test_hierarchical_sort():
    df = pd.DataFrame([[1,7,6,2],[0,0,0,0],[10,10,10,10],[1,1,1,1]])
    result = _hierarchical_sort(df.T)
    assert (result == [2,0,1,3]).all(), 'columns are not properly sorted out'
    


# In[3]:


def test_helix():
    '''testing function produces matplotlib object'''
    assert (str(type(_helix(0,5))) == "<class 'matplotlib.patches.Rectangle'>"), "function _helix failed"


# In[4]:


def test_labels():
    """testing function produces tuple"""
    assert (str(type(_labels(1)))) == "<class 'tuple'>", "function _labels failed"


# In[52]:


def test_sheet():
    """testing function prouduces matplotlib object"""
    assert (str(type(_sheet(1,5)))) == "<class 'matplotlib.patches.FancyArrow'>", "function _sheet failed"


# In[6]:


def test_loop():
    '''testing function produces matplotlib object'''
    assert (str(type(_loop(1,5)))) == "<class 'matplotlib.patches.Rectangle'>", "function _loop failed"

