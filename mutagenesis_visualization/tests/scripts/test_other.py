#!/usr/bin/env python
# coding: utf-8

# # Test for code_other

# In[1]:


import pandas as pd
import numpy as np
import traceback

try:
    from mutagenesis_visualization.main.scripts.code_create_objects import (
        hras_RBD,
        aph_obj,
    )
    from mutagenesis_visualization.main.scripts.code_other import plot_box

except ModuleNotFoundError:
    import import_notebook
    import os
    directory = os.getcwd()
    new_directory = directory.replace('tests', 'main')
    os.chdir(new_directory)

    from code_other import plot_box
    from code_create_objects import (
         hras_RBD,
         aph_obj,
     )
    os.chdir(directory)


# # Test for plot_rank

# In[2]:


def test_plot_rank():
    
    # Create dictionary with objects
    dict_obj = {
        'obj_test_1': hras_RBD(),
        'obj_test_2': aph_obj(),
    }

    # Define aux function
    def _test_plot_rank(obj_test, parameters):
        error = False
        try:
            obj_test.rank(
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
            'mode':'mean',
            'figsize': (3, 2.5),
            'y_label': r'$∆E^i_x$',
            'show': False,
        },
    ]

    # Assert
    for obj_label, obj_test in dict_obj.items(): # Loop over the dictionary
        for parameters in list_params: # Loop over the parameters
            assert _test_plot_rank( # Assert that that set of parameters works on that object
                obj_test,
                parameters,
            ) == False, "plot_rank failed with {} object and {} parameters".format(
                obj_label,
                parameters,
            )


# # Test for plot_miniheatmap

# In[3]:


def test_plot_miniheatmap():
    
    # Create dictionary with objects
    dict_obj = {
        'obj_test_1': hras_RBD(),
#        'obj_test_2': aph_obj(),
    }

    # Define aux function
    def _test_plot_miniheatmap(obj_test, parameters):
        error = False
        try:
            obj_test.miniheatmap(offset=0,
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
            'neworder_aminoacids':list('ACDEFGHIKLMNPQRSTVWY*'),
            'figsize': (3, 2.5),
            'y_label': r'$∆E^i_x$',
            'show': False,
        },
    ]

    # Assert
    for obj_label, obj_test in dict_obj.items(): # Loop over the dictionary
        for parameters in list_params: # Loop over the parameters
            assert _test_plot_miniheatmap( # Assert that that set of parameters works on that object
                obj_test,
                parameters,
            ) == False, "plot_miniheatmap failed with {} object and {} parameters".format(
                obj_label,
                parameters,
            )


# # Test for plot_neighboreffect

# In[5]:


def test_plot_neighboreffect():
    
    # Create dictionary with objects
    dict_obj = {
        'obj_test_1': hras_RBD(),
        'obj_test_2': aph_obj(),
    }

    # Define aux function
    def _test_plot_neighboreffect(obj_test, parameters):
        error = False
        try:
            obj_test.neighboreffect(
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
            assert _test_plot_neighboreffect( # Assert that that set of parameters works on that object
                obj_test,
                parameters,
            ) == False, "plot_neighboreffect failed with {} object and {} parameters".format(
                obj_label,
                parameters,
            )


# # Test for plot_secondary

# In[6]:


def test_plot_secondary():
    
    # Create dictionary with objects
    dict_obj = {
        'obj_test_1': hras_RBD(),
        'obj_test_2': aph_obj(),
    }

    # Define aux function
    def _test_plot_secondary(obj_test, parameters):
        error = False
        try:
            obj_test.secondary_mean(
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
            assert _test_plot_secondary( # Assert that that set of parameters works on that object
                obj_test,
                parameters,
            ) == False, "plot_secondary failed with {} object and {} parameters".format(
                obj_label,
                parameters,
            )


# # Test for plot_roc

# In[7]:


def test_plot_roc():
    
    # fake data for obj_test_1 (change for other objects)
    df_freq= pd.DataFrame()
    df_freq['Variant'] = hras_RBD().dataframe['Variant']
    df_freq['Class'] = np.random.randint(2, size=len(hras_RBD().dataframe))

    # Create dictionary with objects
    dict_obj = {
        'obj_test_1': hras_RBD(),
#        'obj_test_2': aph_obj(),
    }

    # Define aux function
    def _test_plot_roc(obj_test, parameters):
        error = False
        try:
            obj_test.roc(df_freq,
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
            assert _test_plot_roc( # Assert that that set of parameters works on that object
                obj_test,
                parameters,
            ) == False, "plot_roc failed with {} object and {} parameters".format(
                obj_label,
                parameters,
            )


# # Test for plot_cumulative

# In[8]:


def test_plot_cumulative():
    
    # Create dictionary with objects
    dict_obj = {
        'obj_test_1': hras_RBD(),
        'obj_test_2': aph_obj(),
    }

    # Define aux function
    def _test_plot_cumulative(obj_test, parameters):
        error = False
        try:
            obj_test.cumulative(
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
            'mode':'nonSNV',
            'show': False,
        },        {
            'mode':'SNV',
            'figsize': (3, 2.5),
            'y_label': r'$∆E^i_x$',
            'show': False,
        },
    ]

    # Assert
    for obj_label, obj_test in dict_obj.items(): # Loop over the dictionary
        for parameters in list_params: # Loop over the parameters
            assert _test_plot_cumulative( # Assert that that set of parameters works on that object
                obj_test,
                parameters,
            ) == False, "plot_cumulative failed with {} object and {} parameters".format(
                obj_label,
                parameters,
            )


# ## test plot_box

# In[12]:


def test_plot_box():
    # At some point, arguments cannot be vectors, only labels of dataframe
    # Define aux function
    def _test_plot_box():
        error = False
        try:
            plot_box([0,1,2,3],[4,3,2,1])
        except Exception as e:
            tb = traceback.format_exc()
            print(e)
            print(tb)
            error = True
        return error
    
    assert _test_plot_box() == False

