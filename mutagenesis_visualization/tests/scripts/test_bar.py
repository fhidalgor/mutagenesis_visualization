#!/usr/bin/env python
# coding: utf-8

# # Test of code_bar

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


# # Test bar graph with mean enrichment

# In[5]:


def test_plot_mean():
    
    # Create dictionary with objects
    dict_obj = {
        'obj_test_1': hras_RBD(),
        'obj_test_2': aph_obj(),
    }

    # Define aux function
    def _test_plot_mean_output(obj_test, parameters):
        error = False
        try:
            obj_test.mean(
                **parameters
            )  # pass dictionary as arguments of method
        except:
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
            assert _test_plot_mean_output( # Assert that that set of parameters works on that object
                obj_test,
                parameters,
            ) == False, "plot_mean failed with {} object and {} parameters".format(
                obj_label,
                parameters,
            )


# In[6]:


test_plot_mean()


# In[9]:


hras_RBD().mean()


# In[8]:


help(hras_RBD().mean)


# ## Test Mean Differential Graph
# 

# In[38]:


def test_plot_meandifferential():
    
    # Create dictionary with objects
    dict_obj = {
        'obj_test_1': hras_RBD(),
        'obj_test_2': aph_obj(),
    }

    # Define aux function
    def _test_plot_meandifferential_output(obj_test, parameters):
        error = False
        try:
            obj_test.differential(obj_test,
                **parameters
            )  # pass dictionary as arguments of method
        except:
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
            assert _test_plot_meandifferential_output( # Assert that that set of parameters works on that object
                obj_test,
                parameters,
            ) == False, "plot_meandifferential failed with {} object and {} parameters".format(
                obj_label,
                parameters,
            )


# In[39]:


test_plot_meandifferential()


# ## Test Plot Mean Counts Graph

# In[8]:


import traceback
def test_plot_meancounts():
    
    # Create dictionary with objects
    dict_obj = {
        'obj_test_1': hras_RBD(),
        'obj_test_2': aph_obj(),
    }

    # Define aux function
    def _test_plot_meancounts_output(obj_test, parameters):
        error = False
        try:
            obj_test.meancounts(
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
            assert _test_plot_meancounts_output( # Assert that that set of parameters works on that object
                obj_test,
                parameters,
            ) == False, "plot_meancounts failed with {} object and {} parameters".format(
                obj_label,
                parameters,
            )


# In[9]:


test_plot_meancounts()


# ## Test Plot Positional Graph

# In[16]:


#does not work when the position is 1 --should make a note of this somewhere
def test_plot_position():
    
    # Create dictionary with objects
    dict_obj = {
        'obj_test_1': hras_RBD(),
        'obj_test_2': aph_obj(),
    }

    # Define aux function
    def _test_plot_position_output(obj_test, parameters):
        error = False
        try:
            obj_test.position(position = 15,
                **parameters
            )  # pass dictionary as arguments of method
        except:
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
            assert _test_plot_position_output( # Assert that that set of parameters works on that object
                obj_test,
                parameters,
            ) == False, "plot_position failed with {} object and {} parameters".format(
                obj_label,
                parameters,
            )


# In[17]:


test_plot_position()


# In[4]:


import traceback
def test_plot_library_representation():
    
    # Create dictionary with objects
    dict_obj = {
        'obj_test_1': hras_RBD(),
        'obj_test_2': aph_obj(),
    }

    # Define aux function
    def _test_plot_library_representation(obj_test, parameters):
        error = False
        try:
            obj_test.library_representation(
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
            assert _test_plot_library_representation( # Assert that that set of parameters works on that object
                obj_test,
                parameters,
            ) == False, "plot_library failed with {} object and {} parameters".format(
                obj_label,
                parameters,
            )


# In[5]:


test_plot_library_representation()

