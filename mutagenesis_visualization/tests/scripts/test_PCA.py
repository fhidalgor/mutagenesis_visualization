#!/usr/bin/env python
# coding: utf-8

# # Test of code_PCA

# In[1]:


import pandas as pd
import numpy as np

try:
    from mutagenesis_visualization.main.scripts.code_create_objects import (
        hras_RBD,
        aph_obj,
    )
except ModuleNotFoundError:
    import import_notebook
    import os
    directory = os.getcwd()
    new_directory = directory.replace('tests', 'main')
    os.chdir(new_directory)

    import code_create_objects
    from code_create_objects import (
         hras_RBD,
         aph_obj,
     )
    os.chdir(directory)


# # Test of plot_correlation

# In[5]:


import traceback
def test_plot_correlation():
    
    # Create dictionary with objects
    dict_obj = {
        'obj_test_1': hras_RBD(),
        #'obj_test_2': aph_obj(),
    }

    # Define aux function
    def _test_plot_correlation(obj_test, parameters):
        error = False
        try:
            obj_test.correlation(
                **parameters
            )  # pass dictionary as arguments of method
        except Exception as e:
            tb = traceback.format_exc()
            error = True
        return error

    # Define dictionary of parameters
    # Each dict in the list will be a different set of parameters
    list_params = [
        {
            'show': False,
        },
        {
            'figsize': (3, 2.5),
            'y_label': r'$∆E^i_x$',
            'show': False,
        },
    ]

    # Assert
    for obj_label, obj_test in dict_obj.items(): # Loop over the dictionary
        for parameters in list_params: # Loop over the parameters
            assert _test_plot_correlation( # Assert that that set of parameters works on that object
                obj_test,
                parameters,
            ) == False, "plot_correlation failed with {} object and {} parameters".format(
                obj_label,
                parameters,
            )


# # Test of plot_individualcorrelation

# In[4]:


import traceback
def test_plot_individual_correlation():
    
    # Create dictionary with objects
    dict_obj = {
        'obj_test_1': hras_RBD(),
        'obj_test_2': aph_obj(),
    }

    # Define aux function
    def _test_plot_individual_correlation(obj_test, parameters):
        error = False
        try:
            obj_test.individual_correlation(
                **parameters
            )  # pass dictionary as arguments of method
        except Exception as e:
            tb = traceback.format_exc()
            print(e)
            print(tb)
            error = True
        return error

    # Define dictionary of parameters
    # Each dict in the list will be a different set of parameters
    list_params = [
        {
            'show': False,
        },
        {
            'figsize': (3, 2.5),
            'y_label': r'$∆E^i_x$',
            'show': False,
        },
    ]

    # Assert
    for obj_label, obj_test in dict_obj.items(): # Loop over the dictionary
        for parameters in list_params: # Loop over the parameters
            assert _test_plot_individual_correlation( # Assert that that set of parameters works on that object
                obj_test,
                parameters,
            ) == False, "plot_individual_correlation failed with {} object and {} parameters".format(
                obj_label,
                parameters,
            )


# # Test of plot_group_correlation

# In[10]:


import traceback
def test_plot_group_correlation():
    
    # Create dictionary with objects
    dict_obj = {
        'obj_test_1': hras_RBD(),
        #'obj_test_2': aph_obj(),
    }

    # Define aux function
    def _test_plot_group_correlation(obj_test, parameters):
        error = False
        try:
            obj_test.group_correlation(r2=0.75,
                groups = ['DE', 'HKR', 'QN', 'CST', 'AG', 'ILMV', 'WYF', 'P'],
                **parameters
            )  # pass dictionary as arguments of method
        except Exception as e:
            tb = traceback.format_exc()
            error = True
        return error

    # Define dictionary of parameters
    # Each dict in the list will be a different set of parameters
    list_params = [
        {
            'show': False,
        },
        {
            'figsize': (3, 2.5),
            'y_label': r'$∆E^i_x$',
            'show': False,
        },
    ]

    # Assert
    for obj_label, obj_test in dict_obj.items(): # Loop over the dictionary
        for parameters in list_params: # Loop over the parameters
            assert _test_plot_group_correlation( # Assert that that set of parameters works on that object
                obj_test,
                parameters,
            ) == False, "plot_group_correlation failed with {} object and {} parameters".format(
                obj_label,
                parameters,
            )


# # Test of plot_pca

# In[8]:


import traceback
def test_plot_pca():
    
    # Create dictionary with objects
    dict_obj = {
        'obj_test_1': hras_RBD(),
        'obj_test_2': aph_obj(),
    }

    # Define aux function
    def _test_plot_pca(obj_test, parameters):
        error = False
        try:
            obj_test.pca(
                **parameters
            )  # pass dictionary as arguments of method
        except Exception as e:
            tb = traceback.format_exc()
            print(e)
            print(tb)
            error = True
        return error

    # Define dictionary of parameters
    # Each dict in the list will be a different set of parameters
    list_params = [
        {
            'show': False,
        },
        {
            'figsize': (3, 2.5),
            'y_label': r'$∆E^i_x$',
            'show': False,
        },
    ]

    # Assert
    for obj_label, obj_test in dict_obj.items(): # Loop over the dictionary
        for parameters in list_params: # Loop over the parameters
            assert _test_plot_pca( # Assert that that set of parameters works on that object
                obj_test,
                parameters,
            ) == False, "plot_pca failed with {} object and {} parameters".format(
                obj_label,
                parameters,
            )

