#!/usr/bin/env python
# coding: utf-8

# # Test of code_create_objects
# 

# In[8]:


import pandas as pd
import numpy as np

try:
    import mutagenesis_visualization.main.scripts.code_create_objects as code_create_objects
    import mutagenesis_visualization.main.scripts.code_class as code_class

except ModuleNotFoundError:
    import import_notebook
    import os
    directory = os.getcwd()
    new_directory = directory.replace('tests','main')
    os.chdir(new_directory)
    import code_create_objects as code_create_objects
    import code_class as code_class
    os.chdir(directory)


# In[12]:


def test_create_objects():
    
    assert isinstance(code_create_objects.hras_RBD(), code_class.Screen)
    assert isinstance(code_create_objects.bla_obj(), code_class.Screen)
    assert isinstance(code_create_objects.sumo_obj(), code_class.Screen)
    assert isinstance(code_create_objects.mapk1_obj(), code_class.Screen)
    assert isinstance(code_create_objects.ube2i_obj(), code_class.Screen)
    assert isinstance(code_create_objects.tat_obj(), code_class.Screen)
    assert isinstance(code_create_objects.rev_obj(), code_class.Screen)
    assert isinstance(code_create_objects.asynuclein_obj(), code_class.Screen)
    assert isinstance(code_create_objects.aph_obj(), code_class.Screen)
    assert isinstance(code_create_objects.b11L5F_obj(), code_class.Screen)

    

