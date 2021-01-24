#!/usr/bin/env python
# coding: utf-8

# # Test of code_plotly

# In[3]:


import pandas as pd
import numpy as np
import os
import traceback
try:
    from mutagenesis_visualization.main.scripts.code_create_objects import (
        hras_RBD
    )
    location = os.path.dirname(os.path.realpath(__file__))
    pdb_path = os.path.join(location, '../../data', '5p21.pdb')

except ModuleNotFoundError:
    import import_notebook
    directory = os.getcwd()
    new_directory = directory.replace('tests', 'main')
    os.chdir(new_directory)
    pdb_path = '../../data/5p21.pdb'
    from code_create_objects import (hras_RBD)
    os.chdir(directory)


# In[4]:


'''class Student:
    def __init__(self, name, gpa, grad_date,units):
        self.name = name
        self.GPA = gpa
        self.grad_date = grad_date
        self.units = units
    
    def total_GPA(self):
        return self.GPA * self.units
    
Sage = Student('Sage',4.0, 'May 2021', 12)
Kate = Student('Kate', 3.8, 'Dec 2021', 14)'''


# ## Test heatmap

# In[5]:


def test_plot_heatmap_plotly():
    # Get object
    obj_test = hras_RBD()

    # Define aux function
    def _test_plot_heatmap_plotly_output(parameters):
        error = False
        try:
            obj_test.heatmap_plotly(
                **parameters
            )  # pass dictionary as arguments of method
        except:
            error = True
        return error

    # Define dictionary of parameters
    # Each dict in the list will be a different set of parameters
    list_params = [
        {'show': False},
        {
            'show': False,
            'title': 'Changing this parameter for fun',
            'return_plot_object': True,
        },
    ]

    # Assert
    for parameters in list_params:
        assert _test_plot_heatmap_plotly_output(
            parameters
        ) == False, "heatmap_plotly failed with {} parameters".format(
            parameters
        )


# ## Test scatter

# In[4]:


def test_plot_scatter_plotly():  # change
    # Get object
    obj_test = hras_RBD()

    # Define aux function
    def _test_plot_scatter_plotly_output(parameters):  # change
        error = False
        try:
            obj_test.scatter_plotly(
                obj_test, **parameters
            )  # pass dictionary as arguments of method
        except:
            error = True
        return error

    # Define dictionary of parameters
    # Each dict in the list will be a different set of parameters
    list_params = [
        {'show': False},
        {
            'show': False,
            'title': 'Changing this parameter for fun',
            'return_plot_object': True,
        },
        {'mode': 'mean', 'show': False, 'title': 'hello world'},
        {'mode': 'pointmutant', 'show': False, 'title': 'go bears'},
    ]

    # Assert
    for parameters in list_params:
        assert _test_plot_scatter_plotly_output(
            parameters
        ) == False, "scatter_plotly failed with {} parameters".format(
            parameters
        )


# ## Test Rank

# In[5]:


def test_plot_rank_plotly():
    # Get object
    obj_test = hras_RBD()

    # Define aux function
    def _test_plot_rank_plotly_output(parameters):
        error = False
        try:
            obj_test.rank_plotly(
                **parameters
            )  # pass dictionary as arguments of method
        except:
            error = True
        return error

    # Define dictionary of parameters
    # Each dict in the list will be a different set of parameters
    list_params = [
        {'show': False},
        {
            'show': False,
            'title': 'Changing this parameter for fun',
            'return_plot_object': True,
        },
        {'mode': 'mean', 'show': False, 'title': 'hello world'},
        {'mode': 'pointmutant', 'show': False, 'title': 'go bears'},
    ]

    # Assert
    for parameters in list_params:
        assert _test_plot_rank_plotly_output(
            parameters
        ) == False, "rank_plotly failed with {} parameters".format(parameters)


# ## Test Histogram 

# In[7]:


def test_plot_histogram_plotly():
    # Get object
    obj_test = hras_RBD()

    # Define aux function
    def _test_plot_histogram_plotly_output(parameters):
        error = False
        try:
            obj_test.histogram_plotly(
                **parameters
            )  # pass dictionary as arguments of method
        except:
            error = True
        return error

    # Define dictionary of parameters
    # Each dict in the list will be a different set of parameters
    list_params = [
        {'show': False},
        {
            'show': False,
            'title': 'Changing this parameter for fun',
            'return_plot_object': True,
        },
        {'mode': 'mean', 'show': False, 'title': 'hello world'},
        {'mode': 'pointmutant', 'show': False, 'title': 'go bears'},
    ]

    # Assert
    for parameters in list_params:
        assert _test_plot_histogram_plotly_output(
            parameters
        ) == False, "histogram_plotly failed with {} parameters".format(
            parameters
        )


# ## Test Mean

# In[9]:


def test_plot_mean_plotly():
    # Get object
    obj_test = hras_RBD()

    # Define aux function
    def _test_plot_mean_plotly_output(parameters):
        error = False
        try:
            obj_test.mean_plotly(
                **parameters
            )  # pass dictionary as arguments of method
        except:
            error = True
        return error

    # Define dictionary of parameters
    # Each dict in the list will be a different set of parameters
    list_params = [
        {'show': False},
        {
            'show': False,
            'title': 'Changing this parameter for fun',
            'return_plot_object': True,
        },
        {'mode': 'mean', 'show': False, 'title': 'hello world'},
        {'mode': 'A', 'show': False, 'title': 'go bears'},
    ]

    # Assert
    for parameters in list_params:
        assert _test_plot_mean_plotly_output(
            parameters
        ) == False, "mean_plotly failed with {} parameters".format(parameters)


# ## Test 3D Scatter

# In[11]:


def test_plot_scatter_3D_plotly():
    # Get object
    obj_test = hras_RBD()

    # Define aux function
    def _test_plot_scatter_3D_plotly_output(parameters):
        error = False
        try:
            obj_test.scatter_3D_plotly(
                **parameters
            )  # pass dictionary as arguments of method
        except:
            error = True
        return error

    # Define dictionary of parameters
    # Each dict in the list will be a different set of parameters
    list_params = [
        {
            'mode': 'A',
            'pdb_path': pdb_path,
            'title': 'Scatter 3D',
            'squared': True,
            'x_label': 'x',
            'y_label': 'y',
            'z_label': 'z',
            'show': False,
            'return_plot_object': True,
        },
        {
            'mode': 'mean',
            'pdb_path': pdb_path,
            'title': 'Scatter 3D',
            'squared': False,
            'x_label': 'x',
            'y_label': 'y',
            'z_label': 'z',
            'show': False,
            'return_plot_object': True,
        },
    ]

    # Assert
    for parameters in list_params:
        assert _test_plot_scatter_3D_plotly_output(
            parameters
        ) == False, "scatter_3D_plotly failed with {} parameters".format(
            parameters
        )


# ## test plot_scatter_3D_pdbprop_plotly

# In[15]:


def test_plot_scatter_3D_pdbprop_plotly():
    # Get object
    obj_test = hras_RBD()

    # Define aux function
    def _test_plot_scatter_3D_pdbprop_plotly_output(parameters):
        error = False
        try:
            obj_test.scatter_3D_pdbprop_plotly(
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
            'mode':'A',
            'plot': ['Distance', 'SASA', 'log B-factor'],
            'pdb_path': pdb_path,
            'show': False,
            'return_plot_object': True,
        },
        {
            'plot': ['Distance', 'SASA', 'log B-factor'],
            'pdb_path': pdb_path,
            'show': False,
            'return_plot_object': True,
        },
    ]

    # Assert
    for parameters in list_params:
        assert _test_plot_scatter_3D_pdbprop_plotly_output(
            parameters
        ) == False, "scatter_3D_pdbprop_plotly failed with {} parameters".format(
            parameters
        )

