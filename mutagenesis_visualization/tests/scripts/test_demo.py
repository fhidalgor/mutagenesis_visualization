#!/usr/bin/env python
# coding: utf-8

# In[1]:


try:
    from mutagenesis_visualization.main.scripts.code_demo import (
        demo, demo_datasets, demo_pdbs, demo_fasta
    )
except ModuleNotFoundError:
    import import_notebook
    import os
    directory = os.getcwd()
    new_directory = directory.replace('tests', 'main')
    os.chdir(new_directory)
    from code_demo import (demo, demo_datasets, demo_pdbs, demo_fasta)
    os.chdir(directory)


# In[28]:


def test_demo():
    """
    This function will test that demo is capable of generating the
    types of figures ('heatmap', 'miniheatmap', 'mean', 'kernel', 'pca',
    'position', 'secondary_mean', 'correlation', 'individual_correlation') that demo()
    is supposed to. Will raise an error if at least one of the plots does not work.
    """
    def _test_output_demo(argument):
        """
        Aux function for test_demo.
        Will try to run a demo function, will return True if there is an error.
        """
        error = False
        try:
            demo(argument, show=False)
        except:
            error = True
        return error

    arguments = [
        'heatmap', 'miniheatmap', 'mean', 'kernel', 'pca', 'position',
        'secondary_mean', 'correlation', 'individual_correlation'
    ]
    solutions = [_test_output_demo(argument) for argument in arguments]

    failed = [arg for arg, sol in zip(arguments, solutions) if sol is True]

    assert any(solutions) == False, 'the following failed: {}'.format(
        str(failed)
    )


# In[ ]:


def test_demo_pdbs():
    """test that function returns dictionary"""
    assert (
        str(type(demo_pdbs()))
    ) == "<class 'dict'>", "function demo_pdbs failed"


# In[ ]:


def test_demo_datasets():
    """test that function returns dictionary"""
    assert (
        str(type(demo_datasets()))
    ) == "<class 'dict'>", "function demo_datasets failed"


# In[ ]:


def test_demo_fasta():
    """test that function returns dictionary"""
    assert (
        str(type(demo_fasta()))
    ) == "<class 'dict'>", "function demo_fasta failed"
