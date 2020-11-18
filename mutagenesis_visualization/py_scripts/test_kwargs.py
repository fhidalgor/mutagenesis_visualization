#!/usr/bin/env python
# coding: utf-8

# # Import modules 

# In[46]:


try:
    import import_notebook
except ModuleNotFoundError:
    pass

import matplotlib
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import rcParams
import code_kwargs


# # Kwargs

# In[7]:


#test that kwargs outputs a default dictionary
def test_kwargs(): 
    dictionary = type(code_kwargs.kwargs())
    assert dictionary == dict, 'error the kwarg type is not a dictionary'


# In[8]:


test_kwargs()


# In[51]:


#test that _generatecolormap creates a matplotlib color map
def test_generatecolormap(): 
    colormap = type(code_kwargs._generatecolormap())
    assert colormap == matplotlib.colors.LinearSegmentedColormap, 'error the colormap type is not a matplotlib color map'


# In[50]:


test_generatecolormap()


# In[54]:


#!pytest ../py_scripts/test_kwargs.py

