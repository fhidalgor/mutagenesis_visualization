#!/usr/bin/env python
# coding: utf-8

# # Test of code_plotly

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


# In[12]:


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

# In[23]:


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
    list_params = [{'show': False}, {
        'show': False,
        'title': 'Changing this parameter for fun',
    },]

    
    # Assert
    for parameters in list_params:
        assert _test_plot_heatmap_plotly_output(
            parameters
        ) == False, "heatmap_plotly failed with {} parameters".format(
            parameters
        )


# ## Test scatter

# In[ ]:


def test_plot_heatmap_plotly(): # change
    # Get object
    obj_test = hras_RBD()

    # Define aux function
    def _test_plot_heatmap_plotly_output(parameters): # change
        error = False
        try:
            obj_test.heatmap_plotly( # change
                **parameters
            )  # pass dictionary as arguments of method
        except:
            error = True
        return error

    # Define dictionary of parameters
    # Each dict in the list will be a different set of parameters
    list_params = [{'show': False}, { # change
        'show': False,
        'title': 'Changing this parameter for fun',
    },]

    
    # Assert
    for parameters in list_params:
        assert _test_plot_heatmap_plotly_output( # change
            parameters
        ) == False, "heatmap_plotly failed with {} parameters".format( # change
            parameters
        )

