#!/usr/bin/env python
# coding: utf-8

# ## Import Modules

# In[1]:


import pandas as pd
import numpy as np

try:
    from mutagenesis_visualization.main.scripts.code_3D import (
        _centroid, _color_3D_scatter
    )
except ModuleNotFoundError:
    import import_notebook
    import os
    directory = os.getcwd()
    new_directory = directory.replace('tests', 'main')
    os.chdir(new_directory)
    from code_3D import (
        _centroid, _color_3D_scatter
    )
    os.chdir(directory)


# ## Test 3D Scatter functions
# 

# In[5]:


def test_centroid():
    """tests the mean calculated in each of the three dimensions is correct"""
    df = pd.DataFrame({"x": [1, 2, 9], "y": [4, 5, 9], "z": [5, 8, 8]})
    assert _centroid(df) == (4.0, 6.0, 7.0), "function _centroid failed"


# In[ ]:


# what about scores between -1 and 1?
# check correct input for Aminoacid
# should there be 5 columns?


# In[2]:


def test_color_3D_scatter():
    """check that output is expected: returns new dataframe with Color column added"""
    df_1 = pd.DataFrame({
        "Position": [1, 2, 9], "Aminoacid": [1, 2, 3], "Score": [5, -8, 8]
    })
    df_solution = pd.DataFrame({
        "Position": [1, 2, 9], "Aminoacid": [1, 2, 3], "Score": [5, -8, 8],
        "Color": ["red", "blue", "red"]
    })
    df_calculated = _color_3D_scatter(df_1, 'mean', 1, -1)

    assert df_solution.equals(df_calculated) == True

