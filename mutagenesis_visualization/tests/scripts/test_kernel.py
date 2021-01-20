#!/usr/bin/env python
# coding: utf-8

# # Test of code_kernel
# 

# ## Import modules

# In[1]:


import pandas as pd
import numpy as np

try:
    from mutagenesis_visualization.main.scripts.code_create_objects import (
        hras_RBD
    )
except ModuleNotFoundError:
    import import_notebook
    import os
    directory = os.getcwd()
    new_directory = directory.replace('tests', 'main')
    os.chdir(new_directory)
    from code_create_objects import (hras_RBD)
    os.chdir(directory)


# ## Test main functions

# In[2]:


def test_plot_kernel():
    # Get object
    obj_test = hras_RBD()

    # Define aux function
    def _test_plot_kernel_output(parameters):
        error = False
        try:
            obj_test.kernel(
                **parameters
            )  # pass dictionary as arguments of method
        except:
            error = True
        return error

    # Define dictionary of parameters
    # Each dict in the list will be a different set of parameters
    list_params = [{'show': False}, {'cumulative': True, 'show':
                                     False},]

    # Assert
    for parameters in list_params:
        assert _test_plot_kernel_output(
            parameters
        ) == False, "plot_kernel failed with {} parameters".format(parameters)


# In[10]:


def test_plot_histogram():
    # Get object
    obj_test = hras_RBD()

    # Define aux function
    def _test_plot_histogram_output(parameters):
        error = False
        try:
            obj_test.histogram(
                **parameters
            )  # pass dictionary as arguments of method
        except:
            error = True
        return error

    # Define dictionary of parameters
    # Each dict in the list will be a different set of parameters
    list_params = [{'show': False}, {'population': 'SNV', 'show':
                                     False},]

    # Assert
    for parameters in list_params:
        assert _test_plot_histogram_output(
            parameters
        ) == False, "plot_histogram failed with {} parameters".format(parameters)

