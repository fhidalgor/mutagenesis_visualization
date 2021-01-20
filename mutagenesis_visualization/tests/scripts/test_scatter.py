#!/usr/bin/env python
# coding: utf-8

# # Tests for plot_scatter

# In[2]:


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


# In[7]:


import traceback
def test_plot_scatter():
    
    # Create dictionary with objects
    dict_obj = {
        'obj_test_1': hras_RBD(),
        'obj_test_2': aph_obj(),
    }

    # Define aux function
    def _test_plot_scatter(obj_test, parameters):
        error = False
        try:
            obj_test.scatter(obj_test,
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
            assert _test_plot_scatter( # Assert that that set of parameters works on that object
                obj_test,
                parameters,
            ) == False, "plot_scatter failed with {} object and {} parameters".format(
                obj_label,
                parameters,
            )


# In[8]:


test_plot_scatter()

