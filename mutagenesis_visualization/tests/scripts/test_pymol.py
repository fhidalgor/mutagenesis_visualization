#!/usr/bin/env python
# coding: utf-8

# # Import Modules

# In[1]:


import pandas as pd
import numpy as np

try:
    from mutagenesis_visualization.main.scripts.code_pymol import (
        _light_parameters, _array_to_pymol,
    )
except ModuleNotFoundError:
    import import_notebook
    import os
    directory = os.getcwd()
    new_directory = directory.replace('tests', 'main')
    os.chdir(new_directory)
    from code_3D import (
        _light_parameters, _array_to_pymol,
    )
    os.chdir(directory)

try:
    from ipymol import viewer as pymol
except ModuleNotFoundError:
    pass


# # Test map into Pymol functions
# 

# In[ ]:


def test_array_to_pymol():
    """test that output is same as expected"""
    assert _array_to_pymol([1, 2,
                            3]) == '1+2+3', "function _array_to_pymol failed"


# In[ ]:


def _test_light_parameters():
    """check that there is no error when function runs, checks output is function's output"""
    def _test_light_parameters_output():
        """checks if function runs or gives error"""
        try:
            _light_parameters()
        except:
            error = True
            return error

    assert _test_light_parameters_output() == _light_parameters(
    ), "function _light_parameters failed"

